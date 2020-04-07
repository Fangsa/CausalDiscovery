library(pcalg)
library(sparsebn)
library(CAM)
#library(CompareCausalNetworks)

library(MASS)
library(methods)


mkdir <- function(fold){
  if(!dir.exists(fold)){
    dir.create(fold)
  }
}

#ccdr
#' @sparsebn is an R package for learning sparse Bayesian networks and other graphical models
#' from high-dimensional data via sparse regularization.
#' @estimate.dag for directed acyclic graphs (Bayesian networks).
#' @estimate.precision for undirected graphs (Markov random fields).
#' @estimate.covariance for covariance matrices.
ccdr <- function(dat){
  dat <- sparsebnData(t(dat), type = "continuous")
  dags <- estimate.dag(dat)
  lambda <- select.parameter(dags,dat)
  if(!is.infinite(lambda)){
    am <- get.adjacency.matrix(dags[[lambda]])
  }else{
    am <- matrix(0, 10, 10)
  }
  am <- as(am, "matrix")
  return(am)
}

#GES
ges <- function(dat){
  score <- new("GaussL0penObsScore", data = t(dat))
  result <- pcalg::ges(score)
  gesmat <- as(result$essgraph, "matrix")
  gesmat[gesmat] <- 1
  return(gesmat)
}

#GIES
gies <- function(dat){
  score <- new("GaussL0penObsScore", data = t(dat))
  result <- pcalg::gies(score)
  gesmat <- as(result$essgraph, "matrix")
  gesmat[gesmat] <- 1
  return(gesmat)
}

#LiNGAM
lingam_new <- function(dat){
  estDAG <- try(lingam(t(dat)), silent = T)
  if(typeof(estDAG) == "list"){
    return(estDAG$Bpruned)
  }else if(typeof(estDAG) == "character"){
    return(matrix(0, nrow(dat), nrow(dat)))
  }
}

#CAM
#cam <- function(dat){
#  estDAG <- try(CAM(t(dat), scoreName = "SEMGAM", numCores = 4, output = FALSE,
#                variableSel = TRUE, variableSelMethod = selGamBoost, pruning = TRUE,
#                pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = str(0.001))), silent = T)
#  if(typeof(estDAG) == "character"){
#    return(matrix, nrow(dat), nrow(dat))
#  }else if(){
#  }
#}

#FCI
fci_new <- function(dat, numCores = 5){
  cov.mat <- cov(t(dat))
  suffStat <- list(C = cov2cor(cov.mat), n = 10^9)
  res <- fci(suffStat, indepTest=gaussCItest, alpha = 0.9999, p = dim(cov.mat)[1], doPdsep = FALSE, numCores = numCores)
  return(list(am = res@amat, p.mat = res@pMax))
}

#IDA
ida <- function(dat){

}

decorate <- function(m){
  row.names(m) <- paste("sp", seq(0, 9), sep = "")
  colnames(m) <- paste("sp", seq(0, 9), sep = "")
  return(m)
}

for(i in c("R")){#"hubbell", "soi"
  if(i == "R"){
    #s = c("100", "200", "400", "800")
    s = c("S20")#"S30", , "S20", "S40", "S50", "S60"
  }
  for(j in s){
    setwd(paste("D:\\Fangsa\\causal_compare\\Data\\TS_", i, "_", j, sep = ""))
    for(k in seq(100)){
      print(k)

      aa <- read.table(paste("AA_", as.character(k), ".txt", sep = ""))
      ra <- read.table(paste("RA_", as.character(k), ".txt", sep = ""))

      #Causal Additive Model

      #CCDr
      mkdir("CCDr")
      am_aa <- ccdr(aa)
      am_ra <- ccdr(ra)
      
      am_aa <- decorate(am_aa)
      am_ra <- decorate(am_ra)
      write.csv(am_aa, paste("CCDr\\AA_adjacent_matrix_",as.character(k),".csv", sep = ""), quote = F)
      write.csv(am_ra, paste("CCDr\\RA_adjacent_matrix_",as.character(k),".csv", sep = ""), quote = F)

      #GES
      if(sum(am_aa) == 0){
        print(paste(i, j, as.character(k), sep = "\t"))
      }else{
        am_aa <- ges(aa)
        
      }
      mkdir("GES")
    
      am_ra <- ges(ra)
      
      am_aa <- decorate(am_aa)
      am_ra <- decorate(am_ra)
      write.csv(am_aa, paste("GES\\AA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
      write.csv(am_ra, paste("GES\\RA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)

      #GIES
      mkdir("GIES")
      if(sum(am_aa) == 0){
        print(paste(i, j, as.character(k), sep = "\t"))
      }else{
        am_aa <- gies(aa)
      }
      
      am_ra <- gies(ra)
      
      am_aa <- decorate(am_aa)
      am_ra <- decorate(am_ra)
      write.csv(am_aa, paste("GiES\\AA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
      write.csv(am_ra, paste("GiES\\RA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)

      #LiNGAM
      mkdir("LiNGAM")
      am_aa <- lingam_new(aa)
      am_ra <- lingam_new(ra)
      
      am_aa <- decorate(am_aa)
      am_ra <- decorate(am_ra)
      write.csv(am_aa, paste("LiNGAM\\AA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
      write.csv(am_ra, paste("LiNGAM\\RA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)

      #FCI
      mkdir("FCI")
      am_aa <- fci_new(aa)
      am_ra <- fci_new(ra)
      
      am_aa$am <- decorate(am_aa$am)
      am_ra$am <- decorate(am_ra$am)
      am_aa$p <- decorate(am_aa$p)
      am_ra$p <- decorate(am_ra$p)
      write.csv(am_aa$am, paste("FCI\\AA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
      write.csv(am_aa$p.mat, paste("FCI\\AA_p_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
      write.csv(am_ra$am, paste("FCI\\RA_adjacent_matrix_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
      write.csv(am_ra$p.mat, paste("FCI\\RA_p_",as.character(k),".csv", sep = ""), row.names = T, quote = F)
    }
  }
}

