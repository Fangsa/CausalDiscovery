library(RTransferEntropy)
library(future)
plan(multiprocess)


#rule <-function(am, p){
#  m <- am*(p <= 0.05) > 0
#  for(i in seq(dim(m)[2])){
#    X <- i
#    index <- which(m[,i]>0, arr.ind = T)
#    if(index >= 2){
#      index_combn <- combn(index, 2)
#      for(j in seq(dim(index_combn)[2])){
#        Y <- index_combn[1, j]; Z <- index_combn[2, j]
#        XtoY <- m[Y, X]; XtoZ <- m[Z, X]
#        if(m[Y, Z] < XtoY & m[Y, Z] < XtoZ){
#          m[Y, Z] = 0
#        }else if(m[Z, Y] < XtoY & m[Z, Y] < XtoZ){
#          m[Z, Y] = 0
#        }
#      }
#    }
#  }
#  return(m)
#}

mkdir <- function(fold){
  if(!dir.exists(fold)){
    dir.create(fold)
  }
}

te <- function(dat){
  am <- as.data.frame(matrix(0, 10, 10)); p <- as.data.frame(matrix(1, 10, 10))
  row.names(am) = paste("sp", seq(0, 9), sep = "");names(am) = paste("sp", seq(0, 9), sep = "")
  row.names(p) = paste("sp", seq(0, 9), sep = "");names(p) = paste("sp", seq(0, 9), sep = "")
  
  for(i in seq(dim(dat)[1]-1)){
    for(j in seq(i+1, dim(dat)[1])){
      print(i)
      print(j)
      x <- dat[i,]; y <- dat[j,]
      #effective transfer entropy (Marschinski and Kantz 2002)
      shannon_te <- transfer_entropy(x, y, nboot = 1000, shuffles = 100)
      res <- as.data.frame(shannon_te$coef)
      am[i,j] <- res["Y->X", "ete"]; p[i,j] <- res["Y->X", "p-value"]
      am[j,i] <- res["X->Y", "ete"]; p[j,i] <- res["X->Y", "p-value"]
    }
    #m <- rule(am, p)
  }
  return(list(am = am, p = p))
}

for(i in c("beem")){#, , "soi""gLV"
  for(j in c("data")){
    setwd(paste("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\", i,"_", j, sep = ""))
    for(k in seq(1, 50)){
      aa <- read.table(paste("AA_", as.character(k), ".txt", sep = ""), row.names = 1)
      ra <- read.table(paste("RA_", as.character(k), ".txt", sep = ""), row.names = 1)
      
      mkdir("TE")
      res_aa <- te(dat = aa)
      write.csv(res_aa$am, paste("TE\\AA_adjacent_matrix_", as.character(k), ".csv", sep = ""), row.names = T, col.names = T, quote = F)
      write.csv(res_aa$p, paste("TE\\AA_p_", as.character(k), ".csv", sep = ""), row.names = T, col.names = T, quote = F)
      res_ra <- te(dat = ra)
      write.csv(res_ra$am, paste("TE\\RA_adjacent_matrix_", as.character(k), ".csv", sep = ""), row.names = T, col.names = T, quote = F)
      write.csv(res_ra$p, paste("TE\\RA_p_", as.character(k), ".csv", sep = ""), row.names = T, col.names = T, quote = F)
    }
  }
}

meta <- read.table("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\beem_data1\\TE\\metadata_1.txt", header = T)
for(i in c("beem")){#, , "soi""gLV"
  for(j in c("data1")){
    setwd(paste("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\", i,"_", j, sep = ""))
    for(k in seq(1, 50)){
      aa <- read.table(paste("AA_", as.character(k), ".txt", sep = ""), row.names = 1)
      ra <- read.table(paste("RA_", as.character(k), ".txt", sep = ""), row.names = 1)
      
      mkdir("TE")
      for(s in unique(meta$subjectID)){
        res_aa <- te(dat = aa[, meta$subjectID==s])
        write.csv(res_aa$am, paste("TE\\AA_adjacent_matrix_", as.character(k), "_", as.character(s),".csv", sep = ""), row.names = T, col.names = T, quote = F)
        write.csv(res_aa$p, paste("TE\\AA_p_", as.character(k),  "_", as.character(s),".csv", sep = ""), row.names = T, col.names = T, quote = F)
        res_ra <- te(dat = ra[, meta$subjectID==s])
        write.csv(res_ra$am, paste("TE\\RA_adjacent_matrix_", as.character(k), "_", as.character(s),".csv", sep = ""), row.names = T, col.names = T, quote = F)
        write.csv(res_ra$p, paste("TE\\RA_p_", as.character(k), "_", as.character(s),".csv", sep = ""), row.names = T, col.names = T, quote = F)
      }
    }
  }
}

