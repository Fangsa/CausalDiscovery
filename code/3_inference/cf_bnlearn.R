library(bnlearn)

adjacent_matrix <- function(y_arcs){
  am <- matrix(0, nrow = 10, ncol = 10)
  row.names(am) <- paste("sp", seq(10), sep = "")
  colnames(am) <- paste("sp", seq(10), sep = "")
  if(dim(y_arcs)[1] > 0){
      for(i in seq(dim(y_arcs)[1])){
    #am[as.character(y_arcs$to[i]), as.character(y_arcs$from[i])] <- 1
    am[paste0("sp", as.character(as.integer(substr(y_arcs$to[1], 3, 5)))),
       paste0("sp", as.character(as.integer(substr(y_arcs$from[1], 3, 5))))] <- 1
  }
  }
  return(am)
}

f <- function(x, method = "PC", k, ab_type = "RA"){

  if(!dir.exists(method)){
    dir.create(method)
  }

  x <- t(x)
  x <- as.data.frame(x)
  if(method == "PC"){
    y <- pc.stable(x)
  }else if(method == "GS"){
    y <- gs(x)
  }else if(method == "iamb"){
    y <- iamb(x)
  }else if(method == "mmpc"){
    y <- mmpc(x)
  }else if(method == "mmhc"){
    y <- mmhc(x)
  }
  y_arcs <- as.data.frame(y$arcs)
  write.csv(y_arcs, paste(method,"\\", ab_type,"_edge_list_",as.character(k),".csv", sep = ""), row.names = F, quote = F)

  am <- adjacent_matrix(y_arcs)
  write.csv(am, paste(method,"\\", ab_type,"_adjacent_matrix_",as.character(k),".csv", sep = ""), quote = F)
}

#for(i in c("R", "SF")){#"hubbell", "soi"
#  if(i == "R"){
#    s = c("100", "200", "400", "800")
#  }else{
#    s = c("200", "800")
#  }
#  #s = c("TS_dense",  "TS")#"CS",
#  for(j in s){
#    setwd(paste("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\gLV_CS_", i, "_", j, sep = ""))
for(i in c("R")){#"hubbell", "soi"
  if(i == "R"){
    #s = c("100", "200", "400", "800")
    s = c("S20")#"S30","S50", "S60", "S10", "S20", 
  }
  for(j in s){
    setwd(paste("D:\\Fangsa\\causal_compare\\Data\\TS_", i, "_", j, sep = ""))
    for(k in seq(100)){
    print(k)
    aa <- read.table(paste("AA_", as.character(k), ".txt", sep = ""))
    ra <- read.table(paste("RA_", as.character(k), ".txt", sep = ""))

    
    #PC
    f(aa, method = "PC", k, ab_type = "AA")
    f(ra, method = "PC", k, ab_type = "RA")

    #GS
    f(aa, method = "GS", k, ab_type = "AA")
    f(ra, method = "GS", k, ab_type = "RA")

    #iamb
    f(aa, method = "iamb", k, ab_type = "AA")
    f(ra, method = "iamb", k, ab_type = "RA")

    #mmpc
    f(aa, method = "mmpc", k, ab_type = "AA")
    f(ra, method = "mmpc", k, ab_type = "RA")

    #mmhc
    f(aa, method = "mmhc", k, ab_type = "AA")
    f(ra, method = "mmhc", k, ab_type = "RA")
  }
  }
}
