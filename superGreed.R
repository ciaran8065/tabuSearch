
tabuTSP<-function(size = 10, iters = 100, objFunc = NULL, config = NULL,
                  neigh = size^2, listSize = 9, nRestarts = 10, repeatAll = 1,
                  verbose = FALSE,dist=NULL)
{
  if (size < 2) { 
    stop("error: config too short!")
  }
  if(is.null(dist)){
    stop("error: must provide a distance matrix")
  }
  if (iters < 2) { 
    stop("error: not enough iterations!")
  }
  if (listSize >= size) { 
    stop("error: listSize too big!")
  }
  if (neigh > size^2) { 
    stop("error: too many neighbours!")
  }
  if (is.null(objFunc)) { 
    stop("A evaluation function must be provided. See the objFunc parameter.")
  }
  if (is.null(config)) {
    #config<-getGreedy(size,dist)
    values<-numeric(size)
    confs<-matrix(0,size,size)
    for(i in 1:size){
      vec<-getGreedy(size,dist,i)
      values[i]<-objFunc(vec,dist)
      confs[i,]<-vec
    }
    #values[which.max(values)]
    config<-confs[which.max(values),]
  }
  else if (size != length(config)) {
    stop("Length of the starting configuration != size")
  }
  if (repeatAll < 1) {
    stop("error: repeatAll must be > 0")
  }
  iter <- 1 
  configKeep <- matrix(0, repeatAll * iters * (nRestarts + 3), size) 
  eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3)) 
  for (j in 1:repeatAll) { 
    if (j > 1) { 
      config<-rbind(getRandom(size))
    }
    tabuList <- matrix(0, 1, size^2) 
    listOrder <- matrix(0, 1, listSize^2)
    eUtility <- objFunc(config,dist) 
    aspiration <- eUtility 
    
    
    preliminarySearch <- function() {
      configKeep[iter, ] <- config 
      eUtilityKeep[iter] <- eUtility 
      iter <- iter + 1 
      for (i in 2:iters) { 
        neighboursEUtility <- matrix(0, 1, size^2) 
        configTemp <- t(matrix(config, size, neigh)) 
        Neighbours <- c(1:(size^2))
        #print(configTemp)
        if(size>=60){
          configTemp<-findConfP(config)
        }else{
          configTemp<-findConf(config)
        }
        #print(configTemp)
        neighboursEUtility[Neighbours] <- apply(configTemp, 1,objFunc,d=dist)
        
        maxNontaboo <- max(neighboursEUtility[tabuList == 0]) 
        if(length(neighboursEUtility[tabuList==1])==0){ 
          maxTaboo<-(-2147483647)
        }else{
          maxTaboo <- max(neighboursEUtility[tabuList == 1]) 
        }
        move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration, 
                       
                       ifelse(length(which(neighboursEUtility == maxTaboo)) == 1, 
                              which(neighboursEUtility == maxTaboo), 
                              sample(which(neighboursEUtility == maxTaboo), 1)),
                       
                       ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1,
                              which(neighboursEUtility == maxNontaboo & tabuList==0), 
                              sample(which(neighboursEUtility == maxNontaboo & tabuList==0), 1)) 
                       
        ) 
        
        
        if (eUtility >= neighboursEUtility[move]) { 
          tabuList[move] <- 1 
          if (sum(tabuList) > listSize^2) { 
            tabuList[listOrder[1]] <- 0 
            listOrder[1:listSize] <- c(listOrder[2:listSize], 0) 
          }
          listOrder[min(which(listOrder == 0))] <- move
        }
        else if (neighboursEUtility[move] > aspiration)
          aspiration <- neighboursEUtility[move]
        eUtility <- neighboursEUtility[move] 
        
        config<-configTemp[move,]
        configKeep[iter, ] <- config 
        eUtilityKeep[iter] <- eUtility 
        iter <- iter + 1 
      }
      result = list(aspiration = aspiration, configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter)
      return(result)
    }
    
    
    if (verbose) 
      cat("Preliminary search stage...\n")
    
    result <- preliminarySearch() 
    aspiration <- result$aspiration 
    configKeep <- result$configKeep 
    eUtilityKeep <- result$eUtilityKeep 
    iter <- result$iter 
    temp_asp <- -2147483647
    restarts <- 0 
    
    
    while (temp_asp < aspiration & restarts < nRestarts) { 
      if (verbose)
        cat("Intensification stage...\n")
      eUtility <- max(eUtilityKeep[which(eUtilityKeep!=0)]) 
      temp_asp <- aspiration 
      config <- configKeep[max(which(eUtilityKeep == eUtility)), ] 
      result <- preliminarySearch() 
      aspiration <- result$aspiration 
      configKeep <- result$configKeep
      eUtilityKeep <- result$eUtilityKeep
      iter <- result$iter
      restarts <- restarts + 1 
    }
    if (verbose)
      cat("Diversification stage...\n")
    
    
    config<-rbind(getRandom(size))
    
    eUtility <- objFunc(config,dist) 
    
    frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0)) 
    
    tempTabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize)) 
    
    tabuList<-tl(tempTabuList) 
    
    listOrder <- sample(which(tabuList==1), listSize^2) 
    
    result <- preliminarySearch() 
    iter <- result$iter
    configKeep <- result$configKeep
    eUtilityKeep <- result$eUtilityKeep
    
  }
  
  endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, neigh = neigh, listSize = listSize, repeatAll = repeatAll)
  class(endResult) = "tabu"
  return(endResult)
  
}


evaluate<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
    pathLength<-pathLength+(d[i,which(v==i)]) 
  }
  return(-1*pathLength)
}

findConf<-function(v){ #non-parallel neighbour finding
  m<-matrix(v,1,length(v))
  
  for(i in 1:length(v)){
    for(j in 1:length(v)){
      m<-rbind(m,swap(v,i,j))
    }
  }
  return(m[2:nrow(m),])
}
findConfP<-function(v){ #parallel neighbour finding
  if(!"foreach" %in% rownames(installed.packages())){
    install.packages("foreach")
  }
  if("doParallel" %in% rownames(installed.packages())){
    install.packages("doParallel")
  }
  nc<-detectCores()-1
  cl<-makeCluster(nc)
  library(doParallel)
  library(foreach)
  m<-matrix(v,1,length(v))
  avec<-c(1:length(v))
  bvec<-c(1:length(v))
  mt<-
    foreach(b=bvec, .combine='rbind')%:%
    foreach(a=avec, .combine='rbind', .export='Pswap') %dopar%{
      Pswap(m,a,b)
    }
  stopCluster(cl)
  return(mt)
}


swap<-function(v,X1,X2){ 
  originalV<-v
  if(X1==X2){
    return(v)
  }
  X1p<-v[X1] 
  if(X1p==X2){ 
    return(v)
  }
  pX1<-which(v==X1) 
  pX2<-which(v==X2) 
  v[pX1]<-X1p 
  v[X1]<-X2 
  v[pX2]<-X1 
  bool<-T 
  for(i in 1:length(v)){
    if(v[i]==i){
      bool<-F
    }
  }
  
  if(bool==F){ 
    return(originalV)
  }
  return(v)
}

tl<-function(b){ 
  v<-numeric(length(b)^2)
  track<-1
  for(i in 1:length(b)){
    if(b[i]==1){
      v[track:(track+length(b)-1)]<-1
      track<-(track+length(b))
    }else{
      track<-(track+length(b))
    }
  }
  return(v)
}

generate<-function(size){
  v<-numeric(size) 
  t<-c(1:size) 
  for(i in 1:size){
    exc<-c(i,v,which(v==i)) 
    
    temp<-t[-exc] 
    
    if(i==size && (temp==size || length(temp)==0)){ 
      v<-generate(size) 
    }else{
      if(length(temp)==1){
        v[i]<-temp
      }else{
        v[i]<-sample(t[-exc],1) 
      }
    }
    
  }
  return(v)
}

getRandom<-function(size){
  loop<-T
  v<-numeric(size)
  while(loop){
    v<-generate(size)
    if(valid(v)==T){
      return(v)
    }
  }
}

valid<-function(b){ 
  if(length(b)<1){
    return(F)
  }else if(length(b)==1){
    return(T)
  }else{
    track<-1
    max<-length(b)
    bool<-T 
    visited<-numeric(length(b))
    current<-1
    for(i in 1:max){
      visited[i]<-b[current] 
      current<-b[current]
      
    }
    
  }
  
  return(!any(duplicated(visited)))
}
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

getGreedy<-function(size,d,start){
  diag(d)<-diag(d)+.Machine$integer.max
  v<-numeric(size)
  track<-start
  exc<-numeric(size)
  exc[1]<-start
  for(i in 1:size){
    exc<-c(track,exc)
    if(i==size){
      v[track]<-start
    }else{
      if(length(d[track,-exc])==1){ #it was being forced to a vector when there was length 1 and it was causing issues with "names()"
        vect2<-which(d[track,]==d[track,-exc]) #PROBLEM this can produce more than 1 number when there is 2 numbers exactly the same in the same row
        pos<-setdiff(vect2,exc) #This takes the differene between the 2 vectors i.e removes the excluded towns from vect2 (if an excluded town has exactly the same distance it can be included)
        v[track]<-pos
      }else{
        v[track]<-as.integer(names(which.min(d[track,-exc])))
      }
      track<-v[track]
    }
  }
  return(v)
}

summ<-function (object, verbose = FALSE, ...)
{
  tabuObject <- object
  nVars <- rowSums(tabuObject$configKeep)
  nSelect <- colSums(tabuObject$configKeep)
  uniqueConfig <- dim(unique(tabuObject$configKeep))[1]
  output <- paste("Tabu Settings", "\n", "  Type                                       = ",
                  tabuObject$type, "\n", "  No of algorithm repeats                    = ",
                  tabuObject$repeatAll, "\n", "  No of iterations at each prelim search     = ",
                  tabuObject$iters, "\n", "  Total no of iterations                     = ",
                  length(nVars), "\n", "  No of unique best configurations           = ",
                  uniqueConfig, "\n", "  Tabu list size                             = ",
                  tabuObject$listSize, "\n", "  Configuration length                       = ",
                  length(nSelect), "\n", "  No of neighbours visited at each iteration = ",
                  tabuObject$neigh, "\n", sep = "")
  maxObjFunction <- max(tabuObject$eUtilityKeep)
  optimumNVars <- nVars[which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))]
  optimumIteration <- which(tabuObject$eUtilityKeep == max(tabuObject$eUtilityKeep))
  if(length(optimumIteration)==1){
    optimumConfig <- matrix(tabuObject$configKeep[optimumIteration, ],1,length(nSelect))
  }else{
    optimumConfig<-tabuObject$configKeep[optimumIteration, ]
  }
  optionPart <- paste("Results:", "\n", "  Shortest path found              = ",
                      -1*maxObjFunction, "\n", "  Occurs # of times                = ",
                      length(optimumNVars), "\n", "  Optimum number of variables      = ",
                      unique(optimumNVars), "\n", sep = "")
  cat(output)
  cat(optionPart)
  if (verbose) {
    cat(paste("Optimum configuration found:", "\n"))
    print(unique(optimumConfig))
    cat(paste("Optimum path found:","\n"))
    print(getPath(length(nSelect),tabuObject$configKeep[which.max(tabuObject$eUtilityKeep),]))
  }
}


Xs<-sample(1:100,30,replace=T) 
Ys<-sample(1:100,30,replace=T) 
Ps<-cbind(Xs,Ys) 
d<-as.matrix(dist(Ps))

m<-as.matrix(eurodist)
rownames(m)<-c(1:21)
colnames(m)<-c(1:21)


#time<-proc.time()
res<-tabuTSP(size=21,iters=50,objFunc=evaluate,listSize=15,nRestarts=10,repeatAll=1,dist=m)
summ(res, verbose=T) #Worked for 125 towns, took approx 90 minutes
#proc.time()-time

plot(Ps,col="black",pch=16)
adjXs<-c(Xs,Xs[1])
adjYs<-c(Ys,Ys[1])
lines(adjXs,adjYs)
resM<-res$configKeep[which.max(res$eUtilityKeep),]

newPs<-function(Ps,resM){
  visited<-numeric(length(Ps[,1]))
  current<-1
  for(i in 1:length(Ps[,1])){
    visited[i]<-resM[current]
    current<-resM[current]
  }
  cycle<-c(1,visited)
  Xs<-Ps[cycle,1]
  Ys<-Ps[cycle,2]
  plot(Ps,col="black",pch=16)
  lines(Xs,Ys)
}
newPs(Ps,resM)