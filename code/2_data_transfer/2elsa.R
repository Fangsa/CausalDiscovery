mkdir <- function(fold){
  if(!dir.exists(fold)){
    dir.create(fold)
  }
}


for(i in c("beem")){
  for(j in c("data1")){
    setwd(paste("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\", i,"_", j, sep = ""))
    mkdir("eLSA")
    for(k in seq(100)){
      RA <- read.table(paste("RA_", as.character(k), ".txt", sep = ""), row.names = 1)
      #RA <- RA*10000
      RA <- rbind(as.character(seq(dim(RA)[2])), RA)
      row.names(RA)[1] <- "#OTU"
      write.table(RA, paste("eLSA\\RA_", as.character(k), ".txt", sep = ""), sep = "\t", row.names = T, col.names = F, quote = F)
      
      AA <- read.table(paste("AA_", as.character(k), ".txt", sep = ""), row.names = 1)
      #RA <- RA*10000
      AA <- rbind(as.character(seq(dim(RA)[2])), AA)
      row.names(AA)[1] <- "#OTU"
      write.table(AA, paste("eLSA\\AA_", as.character(k), ".txt", sep = ""), sep = "\t", row.names = T, col.names = F, quote = F)
      
    }
  }
}