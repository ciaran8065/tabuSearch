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
valids<-numeric(100)
for(i in 1:1000){
Xs<-sample(1:100,20)
Ys<-sample(1:100,20)
Ps<-cbind(Xs,Ys)
d<-as.matrix(dist(Ps))

vec<-getGreedy(20,d,10)
valids[i]<-valid(vec)
}
valids

evaluate<-function(v,d){
  pathLength<-0
  for(i in 1:length(v)){
    pathLength<-pathLength+(d[i,which(v==i)]) 
  }
  return(-1*pathLength)
}


Xs<-sample(1:100,100)
Ys<-sample(1:100,100)
Ps<-cbind(Xs,Ys)
d<-as.matrix(dist(Ps))

values<-numeric(100)
confs<-matrix(0,100,100)
for(i in 1:100){
  vec<-getGreedy(100,d,i)
  values[i]<-evaluate(vec,d)
  confs[i,]<-vec
}
values[which.max(values)]
confs[which.max(values),]