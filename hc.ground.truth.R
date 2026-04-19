hc.ground.truth=function(data,k,ground_truth){
  crit=function(clust,k,s1,s2,d){ #Function to return CVIs for each clustering
    t=table(clust,s1)
    f=function(tab) {
      r=tryCatch({
        fisher.test(tab)
      }, error = function(e) {
        message("Error. Simulating p-value.")
        return(fisher.test(tab, simulate.p.value = TRUE))
      })
      return(r)
    }
    p=f(t)$p.value
    cs=chisq.test(t)$statistic
    chin=cs/(nrow(s2)*(k-1))
    cvi=cluster.stats(d,clustering=clust,alt.clustering = s1)
    cri=c(k,cvi$avg.silwidth,length(which(silhouette(clust,d)[,3]<0))/nrow(s2),cvi$ch,cvi$dunn2,cvi$vi/log(nrow(s2)),cvi$corrected.rand,chin,p) #Validation indexes
    names(cri)=c("K","Avg. Sil.","Neg.Sil.","CH","Dunn*","VI","ARI","Chi","P-value")
    return(cri)
  }
  library(cluster)
  library(clusterCrit)
  library(ggfortify)
  source("ind.scores.R")
  s1=data[,1:9]
  s1=as.integer(s1[ground_truth][[1]])#Ground truth
  s2=scale(data[,10:19])#Scaled data
  di=dist(s2)
  pc=prcomp(s2,rank.=2)#Principal components for Biplot
  cvis=list(a=NA,b=NA,c=NA,d=NA,e=NA,f=NA,g=NA,h=NA,i=NA,j=NA)#List of CVIs
  hc=list(a=NA,b=NA,c=NA,d=NA,e=NA,f=NA,g=NA,h=NA,i=NA,j=NA)
  hc$a=hclust(di,method = "single") #Testing clusterings for each combination of distance and linkage.
  cvis$a=crit(cutree(hc$a,k),k,s1,s2,di)
  hc$b=hclust(di,method = "complete")
  cvis$b=crit(cutree(hc$b,k),k,s1,s2,di)
  hc$c=hclust(di,method = "average")
  cvis$c=crit(cutree(hc$c,k),k,s1,s2,di)
  hc$d=hclust(di,method = "ward.D")
  cvis$d=crit(cutree(hc$d,k),k,s1,s2,di)
  hc$e=hclust(di,method = "ward.D2")
  cvis$e=crit(cutree(hc$e,k),k,s1,s2,di)
  di=dist(s2,method="manhattan")
  hc$f=hclust(di,method = "single")
  cvis$f=crit(cutree(hc$f,k),k,s1,s2,di)
  hc$g=hclust(di,method = "complete")
  cvis$g=crit(cutree(hc$g,k),k,s1,s2,di)
  hc$h=hclust(di,method = "average")
  cvis$h=crit(cutree(hc$h,k),k,s1,s2,di)
  hc$i=hclust(di,method = "ward.D")
  cvis$i=crit(cutree(hc$i,k),k,s1,s2,di)
  hc$j=hclust(di,method = "ward.D2")
  cvis$j=crit(cutree(hc$j,k),k,s1,s2,di)
  scores=ind.scores(cvis) #CVIs scores
  i=which(scores==max(scores))[1]
  hc=hc[[i]]
  print(plot(hc))
  rect.hclust(hc,k,border="red")#Dendrogram
  hc=cutree(hc,k)
  print(autoplot(pc,col=hc))
  print(autoplot(pc,col=as.numeric(s1)))#Biplots
  cri=cvis[[i]]
  l=list("Hierarchical Clustering"=hc,"CVI"=cri,"Which_One"=i,"Scores"=scores)
  return(l)
}
