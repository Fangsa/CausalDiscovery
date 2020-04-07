library(Matrix)

prediction <- function(am_path, p_path = NULL){
  am <- read.csv(am_path)
  if(!is.null(p_path)){
    p <- read.csv(p_path)
    return(list(am, p))
  }
  return(m)
}

inter <- function(path){
  #im <- read.table(path)
  im <- read.csv(path, row.names = 1)
  return(im)
}

TFPN <- function(estGraph, trueGraph){
  tp = 0; tn = 0; fp = 0; fn = 0
  for(i in seq(dim(trueGraph)[1])){
    for(j in seq(dim(trueGraph)[1])){
      if(i != j){
        if(trueGraph[i, j]){
          if(estGraph[i, j]){
            tp  = tp + 1
          }else{
            fn = fn + 1
          }
        }
        else{
          if(estGraph[i, j]){
            fp = fp + 1
          }else{
            tn = tn + 1
          }
        }
      }
    }
  }
  return(list(tp = tp, tn = tn, fp = fp, fn = fn))
}

StructHamDist <- function(estGraph, trueGraph){
  m1 <- as.data.frame(estGraph); diag(m1)<-0
  m2 <- as.data.frame(trueGraph); diag(m2)<-0
  m1[m1 != 0] <- 1
  m2[m2 != 0] <- 1
  
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ds > 0] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ds < 0] <- m2[ds < 0]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}

evaluation <- function(estGraph, trueGraph){
  #Compute Structural Hamming Distance (SHD)
  StructHammingDist <- StructHamDist(estGraph, trueGraph)
  
  #structural intervention distance (SID)
  #StructInterventionDist=Inf
  #tryCatch({StructInterventionDist <-structIntervDist(trueGraph, estGraph, output = F)$sid},
  #                                   error = function(e){
  #                                     StructInterventionDist <- Inf
  #                                   })
  
  res <- TFPN(estGraph, trueGraph)
  tp <- res$tp; tn <- res$tn; fp <- res$fp; fn <- res$fn
  #Specificity
  Specificity<-tn/(tn + fp)
  
  #Sensitivity,recall
  Sensitivity<-tp/(tp + fn)
  
  #Precision
  Precision<-tp/(tp + fp)
  
  #F1
  f1<-2*tp/(2*tp + fp + fn)
  
  #accuracy
  accuracy <- (tp+tn)/(tp+tn+fp+fn)
  return(list(shd = StructHammingDist , spe = Specificity, sen = Sensitivity,
              pre = Precision, F1 = f1, acc = accuracy))
}

BuildEstGraph <- function(k, dir, dt = "AA"){
        #dir <- substr(dir, 3, 100)
  tryCatch({
  if(dir %in% c("ANM", "CDS", "IGCI")){#"eLSA", "TE", "FCI"
    am <- read.csv(paste(dir, "//",dt,"_adjacent_matrix_",as.character(k), ".csv", sep = ""), row.names = 1)
    p <- read.csv(paste(dir, "//",dt,"_p_",as.character(k), ".csv", sep = ""), row.names = 1)
    estGraph <- am*(p<=0.05) 
  }else if(dir %in% c("SparCC")){
    am <- read.table(paste(dir, "//", as.character(k), "_.txt", sep = ""), row.names = 1, sep = "\t")#"RA_am_",sparcc
    am[is.na(am)] <- 0
    p <-read.table(paste(dir, "//", as.character(k), "_pvals.two_sided.txt", sep = ""), row.names = 1)#"p_",
    estGraph <- am*(p<=0.01)
  }else if(dir %in% c("MDSINE")){
      am <- read.csv(paste(dir, "//",dt,"_adjacent_matrix_", as.character(k), ".csv", sep = ""), row.names = 1)
      p <- read.csv(paste(dir, "//",dt,"_p_", as.character(k), ".csv", sep = ""), row.names = 1)
      estGraph <- am*(p>=10)
  }else if(dir %in% c("CCDr", "GES", "GIES", "LiNGAM")){
    am <- read.csv(paste(dir, "//",dt,"_adjacent_matrix_", as.character(k), ".csv", sep = ""), row.names = 1)
    estGraph <- am
  }else if(dir %in% c("GS", "iamb", "mmhc", "mmpc", "PC")){
    am <- read.csv(paste(dir, "//",dt,"_adjacent_matrix_", as.character(k), ".csv", sep = ""))
    estGraph <- am[-1]
  }else if(dir %in% c("beem")){
    am <- read.csv(paste(dir, "//",dt,"_adjacent_matrix_",as.character(k), ".csv", sep = ""))
    p <- read.csv(paste(dir, "//",dt,"_p_",as.character(k), ".csv", sep = ""))
    estGraph <- p[-1]>1
  }else if(dir %in% c("SpiecEasi_GL", "SpiecEasi_MB")){
    am <- read.csv(paste(dir, "//",dt,"_", as.character(k), ".csv", sep = ""))
    estGraph <- am[-1]
  }
  return(estGraph)
  },error = function(e){estGraph = NULL})
}

for(i in c("R", "SF")){#
  for(j in c("S40", "S50", "S60")){#"S10", "S20", "S30", 
    setwd(paste("/Users/fangsa/compare_causal/causal/Data/TS_", i, "_", j, sep = ""))
    #dirs <- list.dirs(full.names = F, recursive = F)
    dirs <- c("SparCC", "SpiecEasi_GL", "SpiecEasi_MB", 
              "GS", "PC", "iamb", "mmpc", "CCDr", "GES", "GIES", "mmhc",
              "beem", "MDSINE", "ANM", "LiNGAM", "CDS", "IGCI" )
    shd<-c(); spe<-c();method<-c(); acc<-c()
    sen<-c(); pre<-c(); F1<-c();dts<-c(); ks <- c()
    for(dir in dirs){
      for(dt in c("AA", "RA")){#
        for(k in seq(100)){
          print(paste(i,j,dir,dt,k,sep = "   "))
          trueGraph <- inter(paste("beta_", as.character(k), ".csv", sep = ""))
          diag(trueGraph) <- 0
          if(dir == "SparCC" & dt == "AA"){
            next
          }else if(dir == "beem" & dt == "AA"){
            next
          }else{
            estGraph <- BuildEstGraph(k,dir, dt)
            if(is.null(estGraph)){
              next
            }
            diag(estGraph) <- 0
            res <- evaluation(estGraph, trueGraph)
            method<-c(method, dir);dts<-c(dts, dt); ks <- c(ks, k)
            shd<-c(shd, res$shd); spe<-c(spe, res$spe); acc<-c(acc, res$acc)
            sen<-c(sen, res$sen); pre<-c(pre, res$pre); F1<-c(F1, res$F1)
          }
        }
      }
    }
    table = as.data.frame(list(method=method,dt = dts, ks = ks, StructHammingDist=shd, Specificity=spe,Sensitivity=sen,
                               Precision=pre, F1 = F1, Accuracy = acc))
    write.csv(table, paste("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/", i,"_", j, "_evaluation.csv", sep = ""))
  }
}