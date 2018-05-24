findDMR<-function(sign=NULL, chrCol="Chr_CpG", startCol = "Start_CpG", dmrSep = 1000){

  #Calculate pairwise distance of consecutive probes
  s = sign[,startCol]
  s2 = s[-1] - s[1:(length(s)-1)]

  #Check chromosome of consecutive probes
  c = sign[,chrCol]
  c2 = c[-1] == c[1:(length(c)-1)]

  #Consider seperate DMRs when two DMRs are locaated >dmrSep apart and are on different chromosomes
  dmrEnd = c(which(s2>dmrSep | c2 == F), nrow(sign))
  dmrStart = c(1,dmrEnd[1:(length(dmrEnd)-1)]+1)
  return(list(dmrStart = dmrStart,dmrEnd=dmrEnd))
}
