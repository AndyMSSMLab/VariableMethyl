enrich<-function(k,K,n,N){
  c((k/K)/(n/N), phyper(k-1, K, N - K, n, lower.tail = F))
}
