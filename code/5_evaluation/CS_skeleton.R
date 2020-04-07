library(SID)
library(pcalg)
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
  #CS
  im <- read.table(path)
  #TS
  #im <- read.csv(path, row.names = 1)
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
  #return(shd)
}

StructHamDist_penalty <- function(estGraph, trueGraph, lambda=0.5){
  shd <- StructHamDist(estGraph, trueGraph)
  diag(estGraph) <- 0; diag(trueGraph) <- 0
  shdp <- shd + lambda*abs(sum(estGraph !=0) - sum(trueGraph !=0))
  print(shd)
  print(lambda*abs(sum(estGraph !=0) - sum(trueGraph !=0)))
  print(sum(trueGraph !=0))
  return(shdp)
}


evaluation <- function(estGraph, trueGraph){
  #Compute Structural Hamming Distance (SHD)
  StructHammingDist <- StructHamDist(estGraph, trueGraph)
  StructHammingDist_penalty <- StructHamDist_penalty(estGraph, trueGraph)
  
  
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
  print(tp+tn+fp+fn)
  return(list(shd = StructHammingDist , spe = Specificity, sen = Sensitivity,
              pre = Precision, F1 = f1, acc = accuracy, shdp = StructHammingDist_penalty))
}

BuildEstGraph <- function(k, dir, dt = "AA"){
  tryCatch({
    if(dir %in% c("ANM", "BF", "CDS", "RECI", "eLSA", "TE")){#, 
      am <- read.csv(paste(dir, "\\",dt,"_adjacent_matrix_",as.character(k), ".csv", sep = ""))
      p <- read.csv(paste(dir, "\\",dt,"_p_",as.character(k), ".csv", sep = ""))
      estGraph <- am[-1]*(p[-1]<=0.05) 
    }else if(dir %in% c("SparCC")){
      am <- read.table(paste(dir, "\\","RA_am_", as.character(k), ".txt", sep = ""), row.names = 1)
      p <-read.table(paste(dir, "\\","p_", as.character(k), ".txt", sep = ""), row.names = 1)
      estGraph <- am*(p<=0.05)
    }else if(dir %in% c("MDSINE")){
      am <- read.csv(paste(dir, "\\",dt,"_adjacent_matrix_", as.character(k), ".csv", sep = ""))
      p <- read.csv(paste(dir, "\\",dt,"_p_", as.character(k), ".csv", sep = ""))
      estGraph <- am[-1]*(p[-1]>=50)
    }else if(dir %in% c("CCDr", "GES", "GIES", "LiNGAM")){
      am <- read.csv(paste(dir, "\\",dt,"_adjacent_matrix_", as.character(k), ".csv", sep = ""), row.names = 1)
      estGraph <- am
    }else if(dir %in% c("GS", "iamb", "mmhc", "mmpc", "PC")){
      am <- read.csv(paste(dir, "\\",dt,"_adjacent_matrix_", as.character(k), ".csv", sep = ""))
      estGraph <- am[-1]
    }else if(dir %in% c("beem")){
      am <- read.csv(paste(dir, "\\RA_adjacent_matrix_",as.character(k), ".csv", sep = ""))
      p <- read.csv(paste(dir, "\\RA_p_",as.character(k), ".csv", sep = ""))
      estGraph <- am[-1]*(p[-1]>0.75) 
    }
    return(estGraph)
  },error = function(e){return("None")}
  )
}


simulate_model = c(); sample_type = c(); method = c();
method_class = c(); count_type = c()
shd<-c(); spe<-c(); acc<-c(); group <- c()
sen<-c(); pre<-c(); F1<-c(); shdp <- c()

for(i in c("gLV_CS")){
  for(j in c("R", "SF")){
    setwd(paste0("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\", i, "_", j))
    dirs = c("GS", "iamb", "mmpc", "PC", "mmhc", "CCDr", "GES", "GIES", "LiNGAM",
             "beem")#"TE", , "MDSINE" 
    for(dir in dirs){
      for(dt in c("RA", "AA")){
        for(k in seq(100)){
          #trueGraph <- inter(paste0("beta_", as.character(k), ".csv"))
          trueGraph <- inter(paste0("inter_", as.character(k), ".txt"))
          diag(trueGraph) <- 0
          trueGraph <- (trueGraph!=0)+t(trueGraph!=0)
          if(dir == "beem" & dt == "AA"){
            next
          }else{
            simulate_model = c(simulate_model, i); 
            sample_type = c(sample_type, j);
            method = c(method, dir);
            count_type = c(count_type, dt)
            estGraph <- BuildEstGraph(k, dir, dt)
            if(estGraph == "None"){
              print(paste(i,j,dir,dt,k,sep = "   "))
              shd<-c(shd, NA); spe<-c(spe, NA); acc<-c(acc, NA)
              sen<-c(sen, NA); pre<-c(pre, NA); F1<-c(F1, NA); group<-c(group, k)
              shdp <- c(shdp, NA)
              next
            }
            diag(estGraph) <- 0
            estGraph <- (estGraph!=0)+t(estGraph!=0)
            res <- evaluation(estGraph, trueGraph)
            
            shd<-c(shd, res$shd); spe<-c(spe, res$spe); acc<-c(acc, res$acc)
            sen<-c(sen, res$sen); pre<-c(pre, res$pre); F1<-c(F1, res$F1); group<-c(group, k)
            shdp <- c(shdp, res$shdp)
          }
        }
      }
    }
  }
}

table = as.data.frame(list(simulate_model = simulate_model, sample_type = sample_type,method=method,
                           count_type = count_type, Group=group, StructHammingDist=shd, Specificity=spe,Sensitivity=sen,
                           Precision=pre, F1 = F1, Accuracy = acc, shdp = shdp))
write.csv(table,"C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\evaluation_CS_skeleton.csv")