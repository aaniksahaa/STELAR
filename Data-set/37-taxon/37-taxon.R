library(ggplot2)
library(stringr)


currentDirectory = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currentDirectory)

data <- read.csv(file="37-taxon-fn-rates.csv", header=TRUE, sep=",")
print(data)

IlsLevel = c()
Genetrees = c()
Basepair = c()
MethodName = c()
Mean = c()
Std = c()
Order = c()


methods = c("astral","stelar","mpest","SuperTriplets")
#ilsLevels = c("scale2d","noscale","scale2u")
ilsLevels = c("noscale")
#ilsLevels2 = c("0.5X","1X","2X")
#basepairs= c("250b","500b","1000b","1500b","true")
basepairs = c("500b")
genes = c("25g", "50g", "100g", "200g", "400g", "800g")
#genes = c("200g")
order = c(1,2,3,4,5,6)
o = 1
for(gene in genes) {
for (method in methods) {
  l2 = 1
  for(ilsLevel  in ilsLevels){
      for(basepair in basepairs) {
          mydata = data[data$basepair == basepair & data$methodName == method & data$genetrees == gene & data$ilsLevel == ilsLevel,]
          IlsLevel = c(IlsLevel, ilsLevels2[l2])
          Genetrees = c(Genetrees, substr(gene, 1, nchar(gene)-1))
          Basepair = c(Basepair, basepair)
            if ( grepl("SuperTriplets", method, fixed = TRUE) == FALSE) {
            MethodName = c(MethodName, toupper(method))
          } 
          else {
            MethodName = c(MethodName, method)
          }
          Mean = c(Mean, mean(mydata$RobinsonFouldsDistance))
          Order = c(Order, order[o])
          Std = c(Std, sd(mydata$RobinsonFouldsDistance)/sqrt(length(mydata$RobinsonFouldsDistance)))
      }
    l2 = l2 + 1
  }
  
}
  o = o + 1
}
print(Mean)
#Mean[16] = 0.05	
#Std[16] = 0.008344859
#print(Mean)  
mydata = data.frame("ilsLevel" = IlsLevel, "genetrees" = Genetrees, "basepair" = Basepair, "methodName" = MethodName, "mean" = Mean, "std" = Std, "order" = Order)

print(mydata)
#nrow(data1)
p =  ggplot(mydata, aes(fill=mydata$methodName, y=mydata$mean, x=reorder(mydata$genetrees, mydata$order)))
p = p + geom_bar(position="dodge", stat="identity",colour="black", width = .8)
p = p + geom_errorbar(aes(ymin=mydata$mean-mydata$std, ymax=mydata$mean+mydata$std),
                      width=0.3,                    # Width of the error bars
                      position=position_dodge(.8))

p =  p + theme(axis.text= element_text(size = 12, face = "bold"))
p =  p + theme(axis.text.y = element_text(size = 12, face = "bold"))
#p = p + scale_y_continuous(limits=c(0,0.20), breaks=seq(0,.20,.04))
p = p +ylab("Avg. Robinson-Foulds distance")
p = p + xlab("number of genes") 
#p = p + xlab("sequence length")                    
#p = p + xlab("ILS level")
#p  = p + facet_grid(cols = vars(data1$basepair))
p =  p + theme(legend.position="right", legend.title=element_blank(), legend.text = element_text(size=8))
p =  p + theme(axis.text= element_text(size = 12, face = "bold"))
p =  p + theme(axis.text.y = element_text(size = 12, face = "bold"))

p

