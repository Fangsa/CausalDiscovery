mkdir <- function(fold){
  if(!dir.exists(fold)){
    dir.create(fold)
  }
}

for(s in c("R")){
  for(t in c("S10", "S20")){#, "S30","S50", "S60", "S40", "S10", "S20"
    setwd(paste0("D:\\Fangsa\\causal_compare\\Data\\TS_", s, "_", t))
    mkdir("MDSINE")
    for(i in seq(100)){
      aa <- read.table(paste0("AA_", as.character(i), ".txt"))
      aa <- cbind(row.names(aa), aa)
      names(aa) <- c("#OTU ID", seq(dim(aa)[2]-1))
      write.table(aa, paste0("MDSINE\\AA_", as.character(i), ".txt"), row.names = F, quote = F, sep = "\t")
      
      ra <- read.table(paste0("RA_", as.character(i), ".txt"))
      ra <- cbind(row.names(ra), ra)
      names(ra) <- c("#OTU ID", seq(dim(ra)[2]-1))
      write.table(ra, paste0("MDSINE\\RA_", as.character(i), ".txt"), row.names = F, quote = F, sep = "\t")
      
      biomass <- read.table(paste0("biomass_", as.character(i), ".txt"), row.names = 1)
      names(biomass) <- c("biomass")
      write.table(biomass, paste0("MDSINE\\biomass_", as.character(i), ".txt"), row.names = F, quote = F)
    }
  }
}

