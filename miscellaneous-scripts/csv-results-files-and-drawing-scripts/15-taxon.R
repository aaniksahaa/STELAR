library(ggplot2)
library(stringr)
library(sqldf)
library(dplyr)

currentDirectory = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currentDirectory)
data <- read.csv(file="15-taxa-fn-rate.csv", header=TRUE, sep=",")
print(data)
Method = c()
Mean = c()
Sd = c()
Labels = c()
Order = c()

count = 0
labels = c("100gene-100bp","1000gene-100bp","100gene-1000bp","1000gene-1000bp","100gene-true","1000gene-true")
methods = c("Astral","stelar","mpest","SuperTriplets")

for(label  in labels){
  count = 1
  #print(temp1)
  for(method in methods) {
    mydata = data[data$modelCondition == label & data$method == method,]
    if (grepl("SuperTriplets", method, fixed = TRUE) == FALSE) {
      Method = c(Method, toupper(method))
    } 
    else {
      Method = c(Method, method)
    }
    Labels = c(Labels,label)
    Mean = c(Mean, mean(mydata$RobinsonFouldsDistance))
    Sd = c(Sd, sd(mydata$RobinsonFouldsDistance/sqrt(length(mydata$RobinsonFouldsDistance))))
    Order = c(Order,count)
    count = count + 1
  }
}

data1 = data.frame(
  "method" = Method,
  "meanfn" = Mean,
  "sdfn" = Sd,
  "lab" = Labels,
  "count" = Order
)

data1$lab <- factor(data1$lab, levels = labels)
print(data1)

p =  ggplot(data1, aes(fill=data1$method, y=data1$meanfn,  x=data1$method))

p = p + geom_bar(position="dodge", stat="identity",colour="black")
p =  p + theme(axis.text= element_text(size = 12, face = "bold"))
p =  p + theme(axis.text.y = element_text(size = 12, face = "bold"))
p = p + geom_errorbar(aes(ymin=data1$meanfn-data1$sdfn, ymax=data1$meanfn+data1$sdfn),
                      width=0.3,                    # Width of the error bars
                      position=position_dodge(.8))+ylab("average FN rate")
p  = p + facet_wrap(~ data1$lab,scale='free_x')
p =  p + theme(legend.position="right",legend.title=element_blank())+ theme(strip.text.x = element_text(size = 12, face="bold"))
p = p +ylab("Avg. Robinson-Foulds distance")
p = p + theme(axis.title.x=element_blank(), axis.text.x=element_blank() , axis.ticks.x=element_blank())
p

