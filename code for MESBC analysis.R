library(pheatmap)

## MESBC function
MESBC <- function(data, K){
  m = nrow(data)
  p = ncol(data)
  n = m+p
  w = data+abs(min(0,range(data)[1]))
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
  return(clust)
}

## true block
pheatmap::pheatmap(simudata,show_colnames =F,show_rownames = F,legend =F,
                   cluster_rows = F,cluster_cols = F,
                   annotation_colors =F,annotation_legend = F)

## disturbed data
data = simudata[sample(1:nrow(simudata),nrow(simudata)),sample(1:ncol(simudata),ncol(simudata))]
pheatmap::pheatmap(data,show_colnames =F,show_rownames = F,legend =F,
                   cluster_rows = F,cluster_cols = F,
                   annotation_colors =F,annotation_legend = F)

## MESBC analysis
MESBC.result=MESBC(data, 5)
m=nrow(data)
p=ncol(data)
rowlabel=MESBC.result[1:m]
collabel=MESBC.result[(m+1):(m+p)]
MESBC.data=data[order(rowlabel),order(collabel)]
pheatmap::pheatmap(MESBC.data,show_colnames =F,show_rownames = F,legend =F,
                   cluster_rows = F,cluster_cols = F,
                   annotation_colors =F,annotation_legend = F)
