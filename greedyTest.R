############################ NEW GREEDY

getGreedy<-function(size,d){
  diag(d)<-diag(d)+.Machine$integer.max
  v<-numeric(size)
  track<-1
  exc<-numeric(size)
  exc[1]<-1
  for(i in 1:size){
    exc<-c(track,exc)
    if(i==size){
      v[track]<-1
    }else{
      if(length(d[track,-exc])==1){ #it was being forced to a vector when there was length 1 and it was causing issues with "names()"
        #cat("track: ",track,"\n")
        #cat("d[track,-exc]",d[track,-exc],"\n")
        #cat("exc: ",exc,"\n")
        #cat("length: ",length(d[track,-exc])==1,"\n")
        #cat("which: ",which(d[track,]==d[track,-exc]),"\n")
        #v[track]<-which(d[track,]==d[track,-exc])
        vect2<-which(d[track,]==d[track,-exc]) #PROBLEM this can produce more than 1 number when there is 2 numbers exactly the same in the same row
        pos<-setdiff(vect2,exc) #This takes the differene between the 2 vectors i.e removes the excluded towns from vect2 (if an excluded town has exactly the same distance it can be included)
        v[track]<-pos
      }else{
        v[track]<-as.integer(names(which.min(d[track,-exc])))
      }
      track<-v[track]
      #cat("V:",v,"\n")
    }
  }
  return(v)
}
evaluate<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
      pathLength<-pathLength+(d[i,which(v==i)]) 
    }
  return(-1*pathLength)
}

############################ OLD VERSION

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

findConfig<-function(size,d){
  initConfigs<-matrix(0,size,size)
  values<-numeric(size)
  for(i in 1:20){
    tConfig<-rbind(getRandom(size))
    initConfigs[i,]<-tConfig
    values[i]<-evaluate(tConfig,d)
  }
  return(initConfigs[which.max(values),])
  
}


########################### Testing
evaluate<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
    pathLength<-pathLength+(d[i,which(v==i)]) 
  }
  return(-1*pathLength)
}

greedyRes<-numeric(100)
oldRes<-numeric(100)
size<-20
for(i in 1:100){
  Xs<-sample(1:100,20)
  Ys<-sample(1:100,20)
  Ps<-cbind(Xs,Ys)
  d<-as.matrix(dist(Ps))
  
  oldRes[i]<-(-1*evaluate(findConfig(size,d),d))
  greedyRes[i]<-(-1*evaluate(getGreedy(size,d),d))
}

summary(greedyRes) #greedy algorithm results for initial configuration
summary(oldRes) #old algorithm results for initial configuration
