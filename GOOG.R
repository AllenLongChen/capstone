library(plyr)
setwd("~/Documents/Capstone/")
df<-read.csv("GOOG130814.csv",header=TRUE)

df[,"symbol1"]<-sapply(df[,"symbol"],function(x){substr(toString(x),1,5)})
df<-subset(df,symbol1!="GOOG7")

toDate<-function(x){
  x<-toString(x)
  as.Date(paste(substr(x,1,4),substr(x,5,6),substr(x,7,8),sep="-"),origin = "1970-01-01")
}

df[,"date"]<-sapply(df[,"date"],function(x){toDate(x)})
df[,"exdate"]<-sapply(df[,"exdate"],function(x){toDate(x)})
df[,"Texp"]<-(df[,"exdate"]-df[,"date"])/365

df[,"strike_price"]<-df[,"strike_price"]/1000
df<-df[,c("cp_flag","strike_price","Texp","best_bid","best_offer")]

df<-df[order(df[,"Texp"],df[,"strike_price"]),]
rownames(df)<-NULL
print(dim(df))

save(df,file="GOOG.RData")