#setwd("C:/Users/Ciaran2/Documents/R")
source("R6tabu.R")


#TESTING BASE TABU

evaluate<-function(config){#c is the configuration
  target<-876 #maximum amount to spend, maybe put this as a function input
  vCoef<-c(150,100,200,60,40,360,20,80,250,160,120,170,190,65,85,20,5,40,85,35) #coefficients of x1..x4
  value<-sum(config*vCoef)
  if(sum(config)==0 || value>target){
    return(0)
  }
  return(value)
}
vals<-numeric(100)
for(i in 1:100){
  tabu<-tabuObj$new(size=20,iters=10,objFunc=evaluate,listSize=4,nRestarts=4)
  res<-tabu$tabuSearch(tabu)
  #tabu$summ(res,T)
  vals[i]<-res$eUtilityKeep[which.max(res$eUtilityKeep)]
}
vals

#TESTING PROX TABU

evaluate2<-function(config){#c is the configuration
  limit<-400 #maximum amount to spend, maybe put this as a function input
  vCoef<-c(350,200,210,150,160,195) #coefficients of x1..x4
  value<-sum(config*vCoef)
  if(sum(config)==0)return(-1*.Machine$integer.max)#smallest integer number possible
  differ<-abs(value-limit) #proximity of the current value to the optimal value
  return(-1*differ)
}

ptabu<-tabuObj$new(size=6,iters=20,objFunc=evaluate2,listSize=4,nRestarts=4)
pres<-ptabu$tabuSearchProx(ptabu)
ptabu$ProxSumm(pres,T)

#TESTING KNAPSACK TABU

vals<-c(50,40,30,60,100,150,120,70,40)
ws<-c(50,30,60,50,50,50,50,200,10)

evaluate3<-function(conf,weights,values,limit){
  Wm<-limit
  vals<-values
  ws<-weights
  cVal<-sum(conf*vals) #config value
  cW<-sum(conf*ws) #config weight
  lim<-Wm+mean(ws)
  res<-numeric(3)
  res[1]<-ifelse(cW>lim,0,cVal) #if the weight is over the limit return 0, else return the value
  res[2]<-cW
  res[3]<-ifelse(cW>Wm && cW<lim,1,0) #if the weight is greater than Wm AND less than the limit return 1 else return 0
  return(res)
}

ktabu<-tabuObj$new(size=9,iters=50,objFunc=evaluate3,listSize=4,nRestarts=10)
kres<-ktabu$tabuSearchKnap(ktabu,weights=ws,values=vals,limit=400)
ktabu$KnapSumm(kres,evaluate3)

#TESTING TSP TABU

Xs<-sample(1:100,21,replace=T) 
Ys<-sample(1:100,21,replace=T) 
Ps<-cbind(Xs,Ys) 
d<-as.matrix(eurodist)
row.names(d)<-c(1:21)
colnames(d)<-c(1:21)

evaluate4<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
    pathLength<-pathLength+(d[i,which(v==i)]) 
  }
  return(-1*pathLength)
}

tsptabu<-tabuObj$new(size=21,iters=50,objFunc=evaluate4,listSize=15,nRestarts=30,repeatAll=1)
tspres<-tsptabu$tabuSearchTSP(tsptabu,d)
tsptabu$TSPSumm(tspres,T)


#getting the order of the path around europe
order<-tspres$configKeep[which.max(tspres$eUtilityKeep),]
nms<-row.names(as.matrix(eurodist))

getPath<-function(size,resM){
  visited<-numeric(size)
  current<-1
  for(i in 1:size){
    visited[i]<-resM[current] 
    current<-resM[current]
  }
  cycle<-c(1,visited)
  return(rbind(cycle))
}
