library(beem)
library(seqtime)
suppressMessages(library(deSolve))  ## for simulation

setwd("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\")

for(iter in seq(100)){
  ## GLV function
  glv <- function(t, x, params){
    with(as.list(params, c(x)), {
      dx <- alpha * x + x * (beta %*% x) 
      list(dx)
    })
  }
  ## integration function
  integ <- function(t, init.x, model, params){
    ode(init.x, t, model, params)
  }
  
  
  indicator = 1
  while(indicator){
    #set.seed()
    p <- 10    ## number of species
    ntp <- 20  ## number of time points
    n <- 40    ## number of subjects
    
    th = TRUE
    while(th){
      tryCatch({
        ################
        alpha <- abs(rnorm(p, 0.2, 0.3)) ## growth parameters
        A = generateA(p, "random", c =0.2) ##random network
        #A = generateA(p, "klemm", c =0.2) ##scale-free network
        #optimal. increase the amount of negative edges
        A <- modifyA(A, perc=70, strength="uniform", mode="negpercent")
        beta <- A 
        write.csv(as.data.frame(beta), paste0("beem_data_R_S40\\beta_", as.character(iter), ".csv"))
        write.csv(as.data.frame(alpha), paste0("beem_data_R_S40\\alpha_", as.character(iter), ".csv"))
        params <- list(alpha=alpha, beta=beta)
        
        ## simulate absolute abundances
        dat <- foreach(i=1:n, .combine=rbind) %do%{
          init.x <- abs(rnorm(p, 0, 10))
          names(init.x) <- sprintf("sp%03d", 1:p)
          out <- integ(0:ntp, init.x, glv, params)
          out[-1,]}
        th = FALSE}, error = function(e){th = TRUE})
    }
    if(sum(dat<0) ==0){
      ### Negative abundance - not a realistic simulation
      #stop("Negative abundances observed. Please try with another set of parameters")
      indicator <- 0
    }
    ## relative abundances
    absAbun <- dat[,-1]
    relAbun <- prop.table(dat[,-1],1)
  }
  
  
  ## abundance trajectroy
  pdf(paste0("beem_data_R_S40\\Abundance_trajectroy_", as.character(iter), ".pdf"))
  par(mfrow = c(2,1))
  matplot(dat[,-1], ylab='Absolute abundance', cex=0.4)
  matplot(relAbun, ylab='Relative abundance', cex=0.4)
  dev.off()
  
  ## metadata
  metadata <- data.frame(sampleID=1:(n*ntp), isIncluded=1, subjectID=rep(1:n, each=ntp), measurementID=1:ntp, 
                         perturbID=0, exptblock=1, intv=c(rep(1,p), rep(-1, n*ntp-p)) ## These are optional
  )
  write.table(metadata, paste0("beem_data_R_S40\\metadata_", as.character(iter), ".txt"), quote = F, sep = "\t")
  ## counts
  counts <- as.data.frame(t(relAbun))
  write.table(as.data.frame(t(absAbun)), paste0("beem_data_R_S40\\AA_", as.character(iter), ".txt"), quote = F, sep = "\t")
  write.table(counts, paste0("beem_data_R_S40\\RA_", as.character(iter), ".txt"), quote = F, sep = "\t")
  
  
  res <- EM(dat=counts, meta=metadata, verbose = TRUE, ncpu = 8)
  
  
  biomassBEEM <- biomassFromEM(res)
  biomassTrue <- rowSums(dat[,-1])
  write.table(biomassTrue, paste0("beem_data_R_S40\\biomass_", as.character(iter), ".txt"))
  #scaling <- median(biomassTrue)/median(biomassBEEM) ## scaling constant
  #pdf(paste0("beem_output\\biomass_", as.character(iter), ".pdf"))
  #par(mfrow = c(2,1))
  #matplot(cbind(biomassBEEM*scaling, biomassTrue), log='y', xlab='Index', ylab='Biomass',cex=0.4)
  #plot(cbind(biomassBEEM*scaling, biomassTrue), xlab='BEEM',ylab='True')
  #dev.off()
  
  paramBEEM <- paramFromEM(res, counts, metadata) ## compute parameters (with the sparse mode)
  write.csv(paramBEEM, paste0("beem_output_R_S40\\paramBEEM_", as.character(iter), ".csv"))
  #pdf(paste0("beem_output1\\paramBEEM_", as.character(iter), ".pdf"))
  #par(mfrow = c(1,2))
  #plot(paramBEEM[paramBEEM$parameter_type=='growth_rate', 4], alpha, xlab='BEEM', ylab='True', main='Growth')
  #plot(paramBEEM[paramBEEM$parameter_type!='growth_rate', 4], c(beta), xlab='BEEM', ylab='True', main='Interaction')
  #dev.off()
  
  
  paramBEEM.nonSparse <- paramFromEM(res, counts, metadata, sparse=FALSE) ## compute parameter with non-sparse mode
  write.csv(paramBEEM, paste0("beem_output_R_S40\\paramBEEM.nonsparse_", as.character(iter), ".csv"))
  #pdf(paste0("beem_output1\\paramBEEM.nonsparse_", as.character(iter), ".pdf"))
  #par(mfrow = c(1,2))
  #plot(paramBEEM.nonSparse[paramBEEM$parameter_type=='growth_rate', 4], alpha, xlab='BEEM', ylab='True', main='Growth')
  #plot(paramBEEM.nonSparse[paramBEEM$parameter_type!='growth_rate', 4], c(beta), xlab='BEEM', ylab='True', main='Interaction')
  #dev.off()
}


