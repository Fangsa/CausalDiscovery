library(seqtime)
library(knitr)
library(ggplot2)
library(reshape2)

#' @title Generate a cross-section dataset
#'
#' @description Generate an abundance dataset,
#' without environmental perturbance, by using generalized Lotka-Volterra, soi and hubbell
#'
#' @param N Number of OTUs(species)
#' @param S Number of samples
#' @param A Interaction matrix (generated from generateA.R)
#' @param count Total number of individuals in dataset
#' @param mode Mode for generateAbundances; default value samples counts from Poisson distribution with lambda count/N
#' @return The abundance dataset


generateDataSet_new = function(samples, matrix, method = "ricker", count = 10000, mode = 4){
  N = length(matrix[,1])
  dataset = matrix(nrow = N, ncol = samples)
  if (method == "ricker"){
    for (i in 1:samples){
      y = generateAbundances(N, count=count, mode=mode)
      y = y*rbinom(length(y), 1, 0.7)
      
      series = ricker(N,A=A,y=(y/sum(y)),K=rep(0.15,N), sigma=-1,tend=500)
      dataset[,i] = series[,301]
    }
  }
  else if (method == "soi") {
    for (i in 1:samples){
      y = generateAbundances(N, count=count, mode=mode)
      y = y*rbinom(length(y), 1, 0.7)
      
      print(paste("soi", as.character(i)))
      series = soi(N, I=10000, A=A, m.vector=y, tend=500)
      dataset[,i] = series[,301]
    }
  }
  else if (method == "hubbell") {
    for (i in 1:samples){
      y = generateAbundances(N, count=count, mode=mode)
      y = y*rbinom(length(y), 1, 0.7)
      
      series = simHubbell(N=N, M=N,I=10000, y=(y/sum(y)), m=0.1, d=1, tskip=500, tend=1000)
      dataset[,i] = series[,301]
    }
  }
  else if(method == "glv"){
    for (i in 1:samples){
      #y = generateAbundances(N, count=count, mode=mode)
      y <- abs(rnorm(N, 0, 10))
      y = y*rbinom(length(y), 1, 0.7)
      
      series = glv(N, matrix, y=y)
      #print(i)
      dataset[,i] = series[,1001]
    }
  }
  dataset[dataset < 0] = 0.00000001 # glv produces very small negative values, these are set to a pseudocount value
  return(dataset)
}

setwd("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e")
N = 10
S = 800
for(i in seq(100)){
  timestart <- Sys.time()
  th <- TRUE
  while(th){
    tryCatch({
      A = generateA(N, "random", c =0.2)
      #optimal. increase the amount of negative edges
      A <- modifyA(A, perc=70, strength="uniform", mode="negpercent")
      write.table(A, paste("gLV_CS_R_800\\inter_", as.character(i),".txt",sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      pdf(file=paste("gLV_CS_R_800\\interaction_", as.character(i),".pdf", sep = ""))
      plotA(A, header="Klemm-Eguiluz interaction matrix")
      dev.off()
      
      #gLV
      print("start glv")
      print(as.character(i))
      
      dataset = generateDataSet_new(S, A, method = 'glv')
      write.table(dataset, paste("gLV_CS_R_800\\Count_", as.character(i), ".txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
      dataset = normalize(dataset)
      dataset = as.matrix(dataset)
      dataset = melt(dataset)
      colnames(dataset) = c("Species", "Sample", "Abundance")
      #pdf(file=paste("gLV_CS//", as.character(i),".pdf", sep = ""))
      ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")
      ggsave(paste("gLV_CS_R_800\\", as.character(i),".pdf", sep = ""))
      #dev.off()
      th <- FALSE
    }, error = function(e){th <- TRUE})
  }
  
  #Ricker
  #dataset = generateDataSet_new(S, A, method = "ricker")
  
  #soi
  #print("start soi")
  #print(as.character(i))
  
  #dataset = generateDataSet_new(samples = S, matrix=A, method = "soi")
  #write.table(dataset, paste("soi_CS\\Count_", as.character(i), ".txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
  
  #dataset <- as.data.frame(dataset)
  #row.names(dataset) <- paste("sp", seq(10), sep = "")
  #colnames(dataset) <- seq(dim(dataset)[2])
  #dataset <- as.matrix(dataset)
  
  #dataset = normalize(dataset)
  #dataset = melt(dataset)
  #colnames(dataset) = c("Species", "Sample", "Abundance")
  #pdf(file=paste("soi_CS//", as.character(i),".pdf", sep = ""))
  #ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")
  #ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")
  #ggsave(paste("soi_CS\\", as.character(i),".pdf", sep = ""))
  #dev.off()
  
  #hubbell
  #print("start hubbell")
  #print(as.character(i))
  
  #dataset = generateDataSet_new(samples = S, matrix=A, method = "hubbell")
  #write.table(dataset, paste("hubbell_CS\\Count_", as.character(i), ".txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")
  
  #dataset <- as.data.frame(dataset)
  #row.names(dataset) <- paste("sp", seq(10), sep = "")
  #colnames(dataset) <- seq(dim(dataset)[2])
  #dataset <- as.matrix(dataset)
  
  #dataset = normalize(dataset)
  #dataset = melt(dataset)
  #colnames(dataset) = c("Species", "Sample", "Abundance")
  #pdf(file=paste("hubbell//", as.character(i),".pdf", sep = ""))
  #ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")
  #ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")
  #ggsave(paste("hubbell_CS\\", as.character(i),".pdf", sep = ""))
  #dev.off()
  
  timeend <- Sys.time()
  print("############################################")
  print(timeend-timestart)
  print("############################################")
}
