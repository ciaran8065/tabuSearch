#setwd("C:/Users/Ciaran2/Documents/R")
source("R6tabu.R")
as.matrix(UScitiesD)

d<-as.matrix(UScitiesD)
row.names(d)<-c(1:10)
colnames(d)<-c(1:10)

evaluate4<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
    pathLength<-pathLength+(d[i,which(v==i)]) 
  }
  return(-1*pathLength)
}

tsptabu<-tabuObj$new(size=10,iters=10,objFunc=evaluate4,listSize=6,nRestarts=10,repeatAll=1)
tspres<-tsptabu$tabuSearchTSP(tsptabu,d)
tsptabu$TSPSumm(tspres,T)

#install.packages("TSP")
#library(TSP)


############################### UScitiesD
tspScores<-numeric(100)
R6tScores<-numeric(100)

for(i in 1:100){
  tsp<-TSP(as.matrix(UScitiesD))
  tspScores[i]<-tour_length(solve_TSP(tsp))
  tspres<-tsptabu$tabuSearchTSP(tsptabu,d)
  R6tScores[i]<-(-1*tspres$eUtilityKeep[which.max(tspres$eUtilityKeep)])
}
plot(tspScores)
plot(R6tScores)

hist(tspScores)
hist(R6tScores)

summary(tspScores)
summary(R6tScores)

########################## eurodist
tspScores<-numeric(100)
R6tScores<-numeric(100)
d<-as.matrix(eurodist)
row.names(d)<-c(1:21)
colnames(d)<-c(1:21)
tsptabu<-tabuObj$new(size=21,iters=10,objFunc=evaluate4,listSize=10,nRestarts=10,repeatAll=1)
for(i in 1:100){
  tsp<-TSP(as.matrix(eurodist))
  tspScores[i]<-tour_length(solve_TSP(tsp))
  tspres<-tsptabu$tabuSearchTSP(tsptabu,d)
  R6tScores[i]<-(-1*tspres$eUtilityKeep[which.max(tspres$eUtilityKeep)])
}
plot(tspScores)
plot(R6tScores)

hist(tspScores)
hist(R6tScores)

summary(tspScores)
summary(R6tScores)