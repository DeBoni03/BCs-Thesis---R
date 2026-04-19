ind.scores=function(c){
  max.1=c[[1]][4]
  min.1=c[[1]][4]
  max.2=c[[1]][5]
  min.2=c[[1]][5]#Initial values for min and max for CH and Dunn
  for (element in c){
    max.1=max(max.1,element[4])
    min.1=min(min.1,element[4])
    max.2=max(max.2,element[5])
    min.2=min(min.2,element[5])
  } #Final values for min and max for CH and Dunn
  score=sapply(c,function(el) {
    k1=c(2,8)#CVIs already in [0,1] where 1 means optimal conditions
    k2=c(3,6,9)#CVIs already in [0,1] where 0 means optimal conditions
    #k3=c(4:5,7) -> CVIs not yet in [0,1] where bigger values mean better conditions
    el[k2]=sapply(el[k2],function(ind.2){
      ind.2=1-ind.2
      return(ind.2)
      })
    el[4]=(el[4]-min.1)/(max.1-min.1)
    el[5]=(el[5]-min.2)/(max.2-min.2)
    el[7]=(el[7]+1)/2
    #Now all CVIs are in [0,1] where 1 means optimal conditions
    w=c(.15,.05,.2,.05,.05,.35,.1,.05)#Weights
    s=t(el[2:9])%*%w
    return(s)
    })
  return(score)
}
