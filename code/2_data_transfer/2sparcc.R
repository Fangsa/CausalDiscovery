mkdir <- function(fold){
  if(!dir.exists(fold)){
    dir.create(fold)
  }
}


for(i in c("R")){#"R", 
  for(j in c("S10","S20")){#, "S30","S40","S50","S60"
    setwd(paste("D:\\Fangsa\\causal_compare\\Data\\TS_", i, "_", j, sep = ""))
    for(k in seq(100)){
      RA <- read.table(paste("RA_", as.character(k), ".txt", sep = ""), row.names = 1)
      RA <- RA*10000
      mkdir("SparCC")
      #mkdir("eLSA")
      RA <- rbind(as.character(seq(dim(RA)[2])), RA)
      row.names(RA)[1] <- "#OTU"
      write.table(RA, paste("SparCC\\RA_", as.character(k), ".txt", sep = ""), sep = "\t", row.names = T, col.names = F, quote = F)
    }
  }
}