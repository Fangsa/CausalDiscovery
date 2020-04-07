setwd("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e")
for(i in c("R_100", "R_200","R_400","R_800","SF_200","SF_800")){
  for(j in seq(100)){
    count <- read.table(paste0("gLV_CS_",i,"\\Count_", as.character(j), ".txt"))
    row.names(count) <- paste0("sp", seq(10))
    names(count) <- seq(dim(count)[2])
    write.table(count, paste0("gLV_CS_", i, "\\AA_", as.character(j), ".txt"), sep = "\t", quote = F)
    RA <- as.data.frame(t(apply(count, 1, function(x)(x/colSums(count)))))
    write.table(RA, paste0("gLV_CS_", i, "\\RA_", as.character(j), ".txt"), sep = "\t", quote = F)
  }
}