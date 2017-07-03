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
      if(length(d[track,-exc])==1){ #it was being forced to a vector when there was length 1 and it was causing issues
        v[track]<-which(d[track,]==d[track,-exc])
      }else{
        v[track]<-as.integer(names(which.min(d[track,-exc])))
      }
      track<-v[track]
      #cat("V:",v,"\n")
    }
  }
  return(v)
}
valid<-function(b){ #checks validity of path, i.e inf it contains loops return false
  if(length(b)<1){
    return(F)
  }else if(length(b)==1){
    return(T)
  }else{
    track<-1
    max<-length(b)
    bool<-T #result to return, will be set to false if loop detected
    visited<-numeric(length(b))#holds the ones already visited, b+1 because we want to come back to 1
    current<-1
    for(i in 1:max){
      visited[i]<-b[current] #the pointer to the next town
      current<-b[current]
      
    }
    
  }
  #print(visited)
  return(!any(duplicated(visited))) # ! used because I want it to return true is there is no duplicates
}
Xs<-sample(1:100,10)
Ys<-sample(1:100,10)
Ps<-cbind(Xs,Ys)
d<-as.matrix(dist(Ps))

vec<-getGreedy(10,d)
valid(vec)