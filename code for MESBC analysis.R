library(pheatmap)
library(scran)
library(igraph)
library(cluster)
library(mclust)
library(clue)
library(ggplot2)

##### Function of generating negative binomial distribution simulation data
gendata_NBP=function(mu,mu0,size,rownum,colnum){
  K=length(rownum)
  N=sum(rownum)
  P=sum(colnum)
  n=rownum
  p=colnum
  
  noise_data=matrix(rnbinom(N*P, mu = mu0, size = size),nrow=N,ncol=P) 
  
  for(i in 1:K){
    
    x=matrix(rnbinom(n[i]*p[i], mu = mu[i], size = size),nrow=n[i],ncol=p[i])
    
    if(i==1) noise_data[1:n[i],1:p[i]]=x
    if(i>1) noise_data[(sum(n[1:(i-1)])+1):sum(n[1:i]),(sum(p[1:(i-1)])+1):sum(p[1:i])]=x
  }
  
  return(noise_data)
  
}


##### MESBC function
MESBC <- function(data, K){
  m = nrow(data)
  p = ncol(data)
  n = m+p
  w = as.matrix(data+abs(min(0,range(data)[1])))
  d1 = apply(w,1,sum)
  d2 = apply(w,2,sum)
  D1 = diag(1/d1^0.5)
  D2 = diag(1/d2^0.5)
  w.standard = D1%*%w%*%D2
  mysvd = svd(w.standard)
  
  clust <- matrix(nrow=n,ncol=length(K))
  colnames(clust)=K
  for(k in K){
    if (k <= p & k <= m){
      u = mysvd$u[,1:k]; v = mysvd$v[,1:k]; lambda = mysvd$d[1:k]
    }
    if (k > p | k > m)
    { u = mysvd$u ; v = mysvd$v }
    u = as.matrix(u); v = as.matrix(v)
    Y = rbind(D1%*%u,D2%*%v)%*%diag(lambda)
    
    clus.out = kmeans(Y,centers=k,iter.max=100,nstart=100)
    clust[,as.character(k)] = clus.out$cluster
  }
  
  if(length(K)==1) clust=clust[,1]
  return(list(clus=clust, Y=rbind(D1%*%mysvd$u,D2%*%mysvd$v)%*%diag(mysvd$d)))
}


##### Function of calculating precision, sensitive, F1 and Accuracy
metric_micro=function(allclus,label,m,n){
  label2=levels(label)
  ktrue=length(label2)
  k=length(unique(allclus))
  ratio=matrix(nrow=k,ncol=ktrue,dimnames = list(paste0("clus",1:k),label2))
  for(j in 1:k)
    for(t in 1:ktrue)
      ratio[j,t]=length(intersect(which(allclus==j),which(label==label2[t])))/length(label==label2[t])
  if(k<=ktrue) matchout=solve_LSAP(ratio, maximum =T)
  if(k>ktrue){
    matchout=numeric()
    tmp=solve_LSAP(t(ratio), maximum =T)
    for(i in 1:ktrue)
      matchout[tmp[i]]=i
    d=k-ktrue
    res=(1:k)[-tmp]
    for(i in 1:d)
      matchout[res[i]]=ktrue+i
  }
  truemat=clusmat=matrix(0,nrow=m,ncol=n-m)
  for(i in 1:ktrue)
    truemat[which(as.numeric(label)[1:m]==i),which(as.numeric(label)[(m+1):n]==i)]=i
  for(i in 1:k)
    clusmat[which(allclus[1:m]==i),which(allclus[(m+1):n]==i)]=matchout[i]
  accurate=length(which(as.numeric(truemat)==as.numeric(clusmat)))/(m*(n-m))
  
  tp=tn=fp=fn=numeric()
  for(i in 1:ktrue){
    clus1=true1=rep(0,m*(n-m))
    clus1[which(as.numeric(clusmat)==i)]=i
    true1[which(as.numeric(truemat)==i)]=i
    tab=table(true1,clus1)
    if(sum(as.numeric(clusmat)==i)==0)  tab=cbind(tab,rep(0,2))
    tp[i]=tab[2,2]
    tn[i]=tab[1,1]
    fp[i]=tab[1,2]
    fn[i]=tab[2,1]
  }
  tp=sum(tp)
  tn=sum(tn)
  fp=sum(fp)
  fn=sum(fn)

  out=list(accurate=accurate,precision=tp/(tp+fp),sensitive=tp/(tp+fn),F1=2*(tp/(tp+fp)*tp/(tp+fn))/(tp/(tp+fp)+tp/(tp+fn)))
  out
}



##### a simple example
load("MESBC sample.Rdata")
# true block
pheatmap::pheatmap(simudata,show_colnames =F,show_rownames = F,legend =F,
                   cluster_rows = F,cluster_cols = F,
                   annotation_colors =F,annotation_legend = F)


# disturbed data
data = simudata[sample(1:nrow(simudata),nrow(simudata)),sample(1:ncol(simudata),ncol(simudata))]
pheatmap::pheatmap(data,show_colnames =F,show_rownames = F,legend =F,
                   cluster_rows = F,cluster_cols = F,
                   annotation_colors =F,annotation_legend = F)


# MESBC analysis
MESBC.result=MESBC(data, 5)
m=nrow(data)
p=ncol(data)
rowlabel=MESBC.result$clus[1:m]
collabel=MESBC.result$clus[(m+1):(m+p)]
MESBC.data=data[order(rowlabel),order(collabel)]
pheatmap::pheatmap(MESBC.data,show_colnames =F,show_rownames = F,legend =F,
                   cluster_rows = F,cluster_cols = F,
                   annotation_colors =F,annotation_legend = F)



##### simulation analysis
rbicnum=c(20,50,50,30,100)
cbicnum=c(20,20,20,40,50)
mu=rep(1,5)
size=1
N=100
ngene=sum(rbicnum)
ncond=sum(cbicnum)
rowend=colend=0
for(i in 1:length(rbicnum)){
  rowend[i+1]=sum(rbicnum[1:i])
  colend[i+1]=sum(cbicnum[1:i])
}


truebic2=list()
for(i in 1:(length(rowend)-1)){
  truebic2$bicr[[paste0("bic",i)]]=(rowend[i]+1):rowend[i+1]
  truebic2$bicc[[paste0("bic",i)]]=ngene+(colend[i]+1):colend[i+1]
}
K=2:10
metric=c("modularity","Silhouette coefficient","ARI",
         "recovery","relevance","precision","sensitive","F1","Accuracy")
metric.mat=matrix(nrow=N,ncol=length(metric),dimnames = list(1:N,metric))
ALLdata=ALLout=list()
ALLkbest=numeric()
true.mod=numeric()
mu0=0.1

for(iteration in 1:N){
  print(paste0("iteration",iteration))
  ALLdata[[iteration]]=ALLout[[iteration]]=list()
  ALLdata[[iteration]]$datao=data=gendata_NBP(mu,mu0,size,rbicnum,cbicnum)
  m=nrow(data)
  p=ncol(data)
  n=m+p
  
  ALLdata[[iteration]]$indr=indr=sample(1:m,m,replace = F)
  ALLdata[[iteration]]$indc=indc=sample(1:p,p,replace = F)
  truebic=truebic2
  for(i in 1:(length(rowend)-1)){
    for(j in 1:length(truebic2[["bicr"]][[paste0("bic",i)]]))
      truebic[["bicr"]][[paste0("bic",i)]][j]=which(indr==truebic2[["bicr"]][[paste0("bic",i)]][j])
    
    for(j in 1:length(truebic2[["bicc"]][[paste0("bic",i)]]))
      truebic[["bicc"]][[paste0("bic",i)]][j]=which((indc+ngene)==truebic2[["bicc"]][[paste0("bic",i)]][j])+ngene
    
  }
  data2=data[indr,indc]
  
  bcout=MESBC(data2,K)
  g <- buildSNNGraph(t(bcout[["Y"]]))
  distY=dist(bcout[["Y"]])
  
  grouptmp=numeric()
  for(i in 1:(length(rowend)-1))
    grouptmp[c(truebic[["bicr"]][[paste0("bic",i)]],truebic[["bicc"]][[paste0("bic",i)]])]=i
  true.mod[iteration]=modularity(g,grouptmp)
  mod=numeric()
  clus=list()
  
  # biclus
  for(k in K){
    clus$MESBC$row[[as.character(k)]]=bcout[["clus"]][1:m,as.character(k)]
    clus$MESBC$col[[as.character(k)]]=bcout[["clus"]][(m+1):n,as.character(k)]
    mod[k]=modularity(g,c(clus$MESBC$row[[as.character(k)]],clus$MESBC$col[[as.character(k)]]))
  }
  
  ALLout[[iteration]]=list(clus=clus,mod=mod)
  
  # analysis
  kbest=ALLkbest[iteration]=as.character(which(mod==max(mod,na.rm = T))[1])
  metric.mat[iteration,"modularity"]=max(mod,na.rm = T)
  metric.mat[iteration,"Silhouette coefficient"]=summary(silhouette(c(clus$MESBC$row[[kbest]],clus$MESBC$col[[kbest]]),distY))[["avg.width"]]
  metric.mat[iteration,"ARI"]=adjustedRandIndex(c(clus$MESBC$row[[kbest]],clus$MESBC$col[[kbest]]),factor(grouptmp))
  
  
  # relevance and recovery
  ktrue=length(truebic$bicr)
  clustmp=c(clus$MESBC$row[[kbest]],clus$MESBC$col[[kbest]])
  jactmp=matrix(nrow=ktrue,ncol=as.numeric(kbest),dimnames=list(paste0("bic",1:ktrue),paste0("clus",1:as.numeric(kbest))))
  for(j in 1:kbest){
    clusind=which(clustmp==j)
    for(t in 1:ktrue){
      bicind=c(truebic$bicr[[paste0("bic",t)]],truebic$bicc[[paste0("bic",t)]])
      jactmp[t,j]=length(intersect(clusind,bicind))/length(union(clusind,bicind))
    }
  }
  metric.mat[iteration,"recovery"]=mean(apply(jactmp,1,max))
  metric.mat[iteration,"relevance"]=mean(apply(jactmp,2,max))
  
  
  # precision, sensitive, F1 and Accuracy
  out=metric_micro(clustmp,as.factor(grouptmp),m,n)
  metric.mat[iteration,"precision"]=out$precision
  metric.mat[iteration,"sensitive"]=out$sensitive
  metric.mat[iteration,"F1"]=out$F1
  metric.mat[iteration,"Accuracy"]=out$accurate
  
}

df=data.frame(metrics=factor(rep(metric),levels = metric),
              score=rep(NA,length(metric)))

for(i in metric){
  mvalue=metric.mat[1,i]
  for(iteration in 2:N)
    mvalue=mvalue+metric.mat[iteration,i]
  mvalue=mvalue/N
  df[which(df$metrics==i),"score"]=mvalue
}

ggplot(df,aes(x=metrics,y=score))+
  geom_bar(stat="identity",position=position_dodge(0.63),width=0.6)+
  ylab("")+ylim(c(0,1))+
  scale_x_discrete(label = c("Modularity","Silhouette\ncoefficient","ARI",
                             "Relevance","Recovery","Precision","Recall",
                             "F1","Accuracy"))+
  theme_bw()+theme(legend.position = "top",
                   axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18),
                   axis.title.x = element_text(size = 19),axis.title.y = element_text(size = 19),
                   plot.title = element_text(size=21),
                   legend.text=element_text(size=28),legend.title=element_text(size=29),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())


