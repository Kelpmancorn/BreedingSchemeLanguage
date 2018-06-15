#'Create a founder population
#'
#'@param sEnv the environment that BSL functions operate in. Default is "simEnv" so use that to avoid specifying when calling functions
#'@param nInd population size
#'
#'@return modifies the list sims in environment sEnv by creating a founder population
#'
#'@export
initializePopulation <- function(sEnv=NULL, nInd=100,nPops=12,subsize=24){
  initializePopulation.func <- function(data, nInd,nPops,subsize){
    seed <- round(stats::runif(1, 0, 1e9))
    md <- data$mapData
    nQTL<-length(md$effectID)
    md$effects <- matrix(stats::rnorm(nQTL), nQTL) #resimulate the QTL effect just in case all gvalue is negative
    #geno <- data$founderHaps * 2 - 1  ####making diploid???????
    #data$founderHaps <- NULL
    ##to randome mate withini subpopulaitons
    geno2 <- NULL
    pedigree2 <- NULL
    for (i in 1:nPops){
      #geno <- geno[sample(nrow(geno), nrow(geno), replace=T),]  ##why do we sample???
      geno <- data$founderHaps * 2 - 1
      #geno <- geno[((i-1)*subsize+1):(subsize*i),][sample(subsize, subsize, replace=T),] #extract subpopulation from the whole 
      geno <- geno[((i-1)*subsize+1):(subsize*i),] # do not sample from subpops
      geno <- randomMate(popSize=nInd, geno=geno, pos=md$map$Pos)
      #pedigree <- cbind(-geno$pedigree, 0) # For founders, parents will be negative
      pedigree <- cbind(-(geno$pedigree+(i-1)*subsize/2), 0) #increment the pedigree id with subsize
      geno2 <- rbind(geno2,geno$progenies) 
      pedigree2 <- rbind(pedigree2,pedigree)
    }
    geno <- geno2
    pedigree <- pedigree2
    colnames(pedigree) <- 1:3
    data$founderHaps <- NULL
    # Genetic effects. This works even if locCov is scalar
    gValue <- calcGenotypicValue(geno=geno, mapData=md)
    covGval <- stats::var(gValue)
    if (any(is.na(covGval)) | sum(diag(covGval)) == 0){
      covGval <- diag(ncol(gValue))
    } else if (Matrix::rankMatrix(covGval) < ncol(gValue)){
      covGval <- diag(diag(covGval))
    }
    coef <- solve(chol(covGval)) %*% chol(data$varParms$locCov)
    md$effects <- md$effects %*% coef
    gValue <- gValue %*% coef
    # Year and location effects: create matrices with zero columns until phenotyped
    #locEffects <- matrix(0, nrow=nInd, ncol=0)
    #locEffectsI <- matrix(0, nrow=nInd, ncol=0)
    #yearEffects <- matrix(0, nrow=nInd, ncol=0)
    #yearEffectsI <- matrix(0, nrow=nInd, ncol=0)
    locEffects <- matrix(0, nrow=nInd*nPops, ncol=0)    # The sum of all subpopulaitons
    locEffectsI <- matrix(0, nrow=nInd*nPops, ncol=0)
    yearEffects <- matrix(0, nrow=nInd*nPops, ncol=0)
    yearEffectsI <- matrix(0, nrow=nInd*nPops, ncol=0)
    
    #genoRec <- data.frame(GID=1:nInd, pedigree=pedigree, popID=0, basePopID=0, hasGeno=FALSE)
    genoRec <- data.frame(GID=1:(nInd*nPops), pedigree=pedigree, popID=0, basePopID=0, hasGeno=FALSE) # GID should not be just nInd
    data$mapData <- md
    #data$nFounders <- nInd
    data$nFounders <- nInd*nPops #nFounders shoulb be sum of founders from all subpopulations
    data$geno <- geno; data$genoRec <- genoRec; data$gValue <- gValue
    data$locEffects <- locEffects; data$locEffectsI <- locEffectsI
    data$yearEffects <- yearEffects; data$yearEffectsI <- yearEffectsI
    return(data)
  }
  
  if(is.null(sEnv)){
    if(exists("simEnv", .GlobalEnv)){
      sEnv <- get("simEnv", .GlobalEnv)
    } else{
      stop("No simulation environment was passed")
    }
  } 
  parent.env(sEnv) <- environment()
  with(sEnv, {
    if(nCore > 1){
      snowfall::sfInit(parallel=T, cpus=nCore)
      sims <- snowfall::sfLapply(sims, initializePopulation.func, nInd=nInd,nPops,subsize)
      snowfall::sfStop()
    }else{
      sims <- lapply(sims, initializePopulation.func, nInd=nInd,nPops,subsize)
    }
  })
}
