library(SpiecEasi)
library(huge)
library(pulsar)
library(igraph)

CalSE <- function(abu, t){
  ###############################################################
  #parameter estimate
  pargs <- list(rep.num=50, ncores=1)
  abu.mb <- spiec.easi(abu, method='mb', lambda.min.ratio=1e-5, nlambda=500,
                       sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs)
  abu.mb.refit <- getRefit(abu.mb)
  write.csv(as.matrix(abu.mb.refit), paste0("SpiecEasi_MB\\", t, "_",  as.character(k), ".csv"))
  
  #AUC PR
  mb.auc <- huge::huge.roc(abu.mb$est$path, graph, verbose=FALSE)
  mb.pr <- stars.pr(getOptMerge(abu.mb), graph, verbose=FALSE)
  # stars selected final network under: se.est$refit$stars
  #mb_auc = c(mb_auc, mb.auc$AUC)
  #mb_pr = c(mb_pr, mb.pr$AUC)
  abu.mb$refit
  print(paste(as.character(mb.auc$AUC), as.character(mb.pr$AUC), sep = "    "))
  
  ################################################################
  abu.gl <- spiec.easi(abu, method='glasso', lambda.min.ratio=1e-3, nlambda=30,
                       sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs)
  abu.gl.refit <- getRefit(abu.mb)
  write.csv(as.matrix(abu.mb.refit), paste0("SpiecEasi_GL\\", t, "_",  as.character(k), ".csv"))
  
  #AUC PR
  gl.auc <- huge::huge.roc(abu.gl$est$path, graph, verbose=FALSE)
  gl.pr <- stars.pr(getOptMerge(abu.gl), graph, verbose=FALSE)
  #gl_auc = c(gl_auc, gl.auc$AUC)
  #gl_pr = c(gl_pr, gl.pr$AUC)
  abu.gl$refit
  print(paste(as.character(gl.auc$AUC), as.character(gl.pr$AUC), sep = "    "))
  
  #################################################################
  save(abu.mb, file = paste("SpiecEasi_MB\\", t, "_", as.character(i), ".Rdata", sep=""))
  save(abu.gl, file = paste("SpiecEasi_GL\\", t, "_", as.character(i), ".Rdata", sep=""))

  return(c(mb.auc$AUC, mb.pr$AUC, gl.auc$AUC, gl.pr$AUC))
}

for(i in c("R")){
    s = c("S20")#"S30","S50", "S60", "S10", "S20", 
  for(j in s){
    setwd(paste("D:\\Fangsa\\causal_compare\\Data\\TS_", i, "_", j, sep = ""))
    for(k in seq(100)){
      print(k)
      index = c(); dt = c(); mb_auc = c(); mb_pr = c()
      gl_auc = c(); gl_pr = c()

  if(!dir.exists("SpiecEasi_MB")){
    dir.create("SpiecEasi_MB")
    dir.create("SpiecEasi_GL")
  }
  print(paste("#####################", as.character(k), "#######################"))
  #load interacton
  inter <- read.csv(paste("beta_" ,as.character(k),".csv", sep = ""),header = T, row.names = 1)
  diag(inter) <- 0
  graph <- as.matrix(apply(inter!=0, 1, as.integer))
  
  #load sequence data
  #abu <- read.table(paste(as.character(i), "\\counts.txt", sep = ""),
  #                  header = F, row.names = 1)
  abu_aa <- read.table(paste("AA_", as.character(k), ".txt", sep = ""),
                    header = T, row.names = 1)
  abu_aa <- as.matrix(t(abu_aa))
  res <- CalSE(abu_aa, "AA")
  index <- c(index, k); dt = c(dt, "AA"); mb_auc = c(mb_auc, res[1])
  mb_pr = c(mb_pr, res[2]); gl_auc = c(gl_auc, res[3]); gl_pr = c(gl_pr, res[4])
  
  abu_ra <- read.table(paste("RA_", as.character(k), ".txt", sep = ""),
                       header = T, row.names = 1)
  abu_ra <- as.matrix(t(abu_ra))
  res <- CalSE(abu_ra, "RA")
  index <- c(index, k); dt = c(dt, "RA"); mb_auc = c(mb_auc, res[1])
  mb_pr = c(mb_pr, res[2]); gl_auc = c(gl_auc, res[3]); gl_pr = c(gl_pr, res[4])
}
auc <- as.data.frame(cbind(index, dt, mb_auc, mb_pr, gl_auc, gl_pr))
names(auc) <- c("Data", "Data_type", "MB_AUC", "MB_pr", "GL_AUC", "GL_pr")
write.csv(auc, "SE_AUC.csv", row.names = F)
  }
}
