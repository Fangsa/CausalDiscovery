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

BuildEstGraph <- function(k, dir, dt = "AA"){
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
      #am <- read.csv(paste(dir, "//",dt,"_adjacent_matrix_",as.character(k), ".csv", sep = ""))
      #p <- read.csv(paste(dir, "//",dt,"_p_",as.character(k), ".csv", sep = ""))
      #estGraph <- p[-1]>1
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
    for(dt in c("AA", "RA")){#
      for(k in seq(100)){
        dirs <- c("SparCC", "SpiecEasi_GL", "SpiecEasi_MB", 
                  "GS", "PC", "iamb", "mmpc", "CCDr", "GES", "GIES", "mmhc",
                  "beem", "MDSINE", "ANM", "LiNGAM", "CDS", "IGCI" )
        result <- matrix(0, length(dirs), length(dirs)+1)
        row.names(result) <- dirs; names(result) <- c("True", dirs)
        for(s in seq(length(dirs)){
          dir1 <- dirs[s]
          print(paste(i , j, dir1, dt, k, sep = "   "))
          trueGraph <- inter(paste("beta_", as.character(k), ".csv", sep = ""))
          if(dir1 == "SparCC" & dt == "AA"){
            next
          }else{
            estGraph1 <- BuildEstGraph(k, dir1, dt)
            if(is.null(estGraph1)){
              next
            }
            res <- StructHamDist(estGraph, trueGraph)
            result[s, 1] <- res
          }
          if(s == length(dirs)){
            next
          }
          
          for(t in seq(s+1, length(dirs))){
            dir2 <- dirs[t]
            if(dir2 == "SparCC" & dt == "AA"){
              next
            }else{
              estGraph2 <- BuildEstGraph(k, dir2, dt)
              if(is.null(estGraph2)){
                next
              }
              res <- StructHamDist(estGraph1, estGraph2)
              result[s, t+1] <- res
            }
          }
        }
        write.csv(result, paste("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/", i,"_", j, "_", dt, 
                                "_", as.character(k), ".csv", sep = ""))
      }
    }
  }
}
