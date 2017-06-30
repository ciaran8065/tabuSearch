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
diag(d)<-diag(d)+.Machine$integer.max #set the main diagonal to the highest integer possible so min will ignore them
d
size<-100
v<-numeric(size)

#do a walk of the towns, i.e start at 1 and find it's closest go there next and repeat
#for the last one we must go back to 1 regardless of it's preference
track<-1
exc<-numeric(size)
exc[1]<-1
for(i in 1:size){
  exc<-c(track,exc)
  if(i==size){
    v[track]<-1
  }else{
  #cat("exc",exc,"\n")
  #cat("d[track, -exc]:",d[track,-exc],"\n")
  #cat("valid options:",names(which.min(d[track,-exc])),"\n")
  if(length(d[track,-exc])==1){
    v[track]<-which(d[track,]==d[track,-exc])
  }else{
    v[track]<-as.integer(names(which.min(d[track,-exc])))
  }
  track<-v[track]
  #cat("V:",v,"\n")
  }
}
v
evaluate(v,d) #calculating the length
valid(v) #ensuring the path is a valid cycle