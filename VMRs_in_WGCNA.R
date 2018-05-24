run<-function(prefix, F_sd_thres, M_sd_thres){
  msign = paste0(prefix,"_M_Sign_sd")
  fsign = paste0(prefix,"_F_Sign_sd")
  msd = paste0(prefix,"_M_sd")
  fsd = paste0(prefix,"_F_sd")
  sdProbes = x[x$NULL_Sign_sd==1 & ((x[,msign] == 1 & x[,msd]>=M_sd_thres) | (x[,fsign] ==1 & x[,fsd]>=F_sd_thres)),1]
  beta = read.table(paste0(prefix,"_norm.txt"), as.is = T, header = T, check.names=F, sep = "\t")
  rownames(beta) = beta[,1]
  beta2 = beta[sdProbes,grep("Beta",colnames(beta),value=T)]
  
  library(WGCNA)
  library(cluster)

  datExpr = t(beta2)
  ADJ1=abs(cor(datExpr,use="p"))^6
  dissTOM=TOMdist(ADJ1)
  hierTOM = hclust(as.dist(dissTOM),method="average");
  minClusterSize = 15
  colorh1 = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,deepSplit=2, pamRespectsDendro = FALSE, minClusterSize=minClusterSize))
  datME=moduleEigengenes(datExpr,colorh1)$eigengenes
  datKME=signedKME(datExpr, datME, outputColumnName="MM.")
  rownames(datKME) = colnames(datExpr)

  g = apply(datKME, 2,function(t){t= unique(unlist(strsplit(unlist(x[names(t[abs(t)> 0.7]),c("ProbeID")]),","))); t[!is.na(t)] })
  sign = x[x$NULL_Sign_sd==1 & (x[,msign] == 1 | x[,fsign] ==1),]
  a = findDMR(sign)$dmrStart
  b = findDMR(sign)$dmrEnd
  coord = unlist(sapply(1:length(a), function(t){ rep(paste0(unique(sign[a[t]:b[t],2]),":",min(sign[a[t]:b[t],3]),"-",max(sign[a[t]:b[t],4])),length(a[t]:b[t])) }))
  names(coord) = rownames(sign)
  k = sapply(g, function(t){ unique(coord[t]) })
  k = k[sapply(k,length)>10]
  max.len <- max(sapply(k, length))
  corrected.list <- lapply(k, function(x) {c(x, rep(NA, max.len - length(x)))})
  mat <- do.call(cbind, corrected.list)
  write.table(mat, paste0(prefix,"_WGCNA_vmr.txt"), row.names = F, sep = "\t", quote = F)
}
