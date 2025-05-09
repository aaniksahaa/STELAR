library(ggplot2)
currentDirectory = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currentDirectory)

data <- read.csv(file="11-RFRateFile.csv", header=TRUE, sep=",")
print(data)
genes = c(5,15,25,50,100)
modelconditions = c("true", "estimated")
methods = c("ASTRAL","STELAR","MPEST","SuperTriplets")


Mean = c()
Sd = c()
Method = c()
Gene = c()
Modelcondition = c()

for (gene in genes) {
  for( modelcondition in modelconditions) {
    for(method in methods) {
        mydata = data[data$method == method & data$modelCondition == modelcondition & data$gene == gene,]

        Mean = c(Mean, mean(mydata$RF_rate))
        Sd = c(Sd, sd(mydata$fn/sqrt(length(mydata$RF_rate))))

        Gene = c(Gene, gene)
        if (grepl("SuperTriplets", method, fixed = TRUE) == FALSE) {
          Method = c(Method, toupper(method))
        }
        else {
          Method = c(Method, method)
        }
        Modelcondition = c(Modelcondition, modelcondition)
    }
  }
}
print(Mean)
data = data.frame(
  "mean" = Mean,
  "sd" = Sd,
  "gene" = Gene,
  "modelcondition" = Modelcondition,
  "method" = Method
)
print(data)

p =  ggplot(data, aes(fill=data$method, y=data$mean, x=reorder(data$gene, data$gene)))
p = p + geom_bar(position="dodge", stat="identity",colour="black")
p =  p + theme(axis.text= element_text(size = 12, face = "bold"))
p =  p + theme(axis.text.y = element_text(size = 12, face = "bold"))
p = p + geom_errorbar(aes(ymin=data$mean-data$sd, ymax=data$mean+data$sd),
                      width=0.3,                    # Width of the error bars
                      position=position_dodge(.9))
p = p + ylab("Avg. Robinson-Foulds distance")
p = p + xlab("number of genes")
p  = p + facet_grid(cols = vars(data$modelcondition)) + theme(strip.text.x = element_text(size = 12, face="bold"))
p =  p + theme(legend.position="bottom", legend.title=element_blank())
#p = p + theme(axis.title.x=element_blank(),
# axis.text.x=element_blank(),
#axis.ticks.x=element_blank())
p

