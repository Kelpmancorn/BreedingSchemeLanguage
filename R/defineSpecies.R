#'Define and create species
#'
#'@param loadData if null create a new species (default), else the file name of previously created species like "fileName" (Do not put "RData" extension).
#'@param importFounderHap if null create new founder haplotypes (default), else the file name of externally generated founder haplotypes in hapmap format.
#'@param saveDataFileName the name to save the species data with double-quotation, like "species1_1". (default: "previousData")
#'@param nSim the number of simulation trials
#'@param nCore simulation processed in parallel over this number of CPUs (Check computer capacity before setting above 1.)
#'@param nChr the number of chromosomes
#'@param lengthChr the length of each chromosome (cM. all chromosomes are the same length)
#'@param effPopSize the effective population size in the base population
#'@param nMarkers the number of markers, which is used especially for genomic selection
#'@param nQTL the number of QTLs controlling the target trait
#'@param propDomi the probability of dominant QTL among the all QTL
#'@param nEpiLoci the expected number of epistatic loci for each effect
#'@param domModel the dominance model: "HetHom" means homozygotes have equal effect but opposite to that of heterozygotes, "Partial": zero means ancestral dominant over derived, one means derived dominant over ancestral, any value in between means partial dominance. At the moment, functionality for the "Partial" option is only available when using ImportFounderHap.
#'
#'@return An environment that contains a list sims with each object of the list being one replicate to initiate a simulation
#'
#'@export
#defineSpecies <- function(loadData=NULL, importFounderHap=NULL, saveDataFileName="previousData", nSim=1, nCore=1, nChr=7, lengthChr=150, effPopSize=100, nMarkers=1000, nQTL=50, propDomi=0, nEpiLoci=0, domModel="HetHom"){
defineSpecies <- function(loadData=NULL, importFounderHap=NULL, saveDataFileName="previousData", nSim=1, nCore=1, nPops=12, nPopsSamples=rep(10,12), nChr=31, lengthChr=1800, effPopSize=10, nMarkers=10000, nQTL=100, propDomi=0, nEpiLoci=0, domModel="HetHom"){
  defineSpecies.func <- function(simNum, nChr, lengthChr, effPopSize, nMarkers, nQTL, propDomi, nEpiLoci, founderHaps=NULL, domModel){
    seed <- round(stats::runif(1, 0, 1e9))
    ##why the nloci is like this???????
    nLoci <- nMarkers + nQTL * (nEpiLoci + 1) * 2
    if (is.null(founderHaps)){
      minMAF <- 0.01
      ##piecesPerM is number of fragments in one Morgen, which is different from GENOME software setting
      piecesPerM <- 10000
      nPiecesPerChr <- lengthChr / 100 * piecesPerM
      recBTpieces <- 1 / piecesPerM
      ##why 2* effpopsize???????????
      coalSim <- getCoalescentSim(nPops=nPops, nPopsSamples=nPopsSamples,effPopSize=2 * effPopSize, nMrkOrMut=nLoci, nChr=nChr, nPiecesPerChr=nPiecesPerChr, recBTpieces=recBTpieces, minMAF=minMAF, seed=seed)
      
      #if (nPops == 1){
        markers <- coalSim$markers
        map <- coalSim$map
        mapData <- makeMap(map=map, nLoci=nLoci, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, interactionMean=nEpiLoci)
      #  final <- list(POP1mapData=mapData, POP1founderHaps=markers)
      #}else{
        ###need to get extract subpopulations, need to figure out where to put????????
      #  final <- NULL
      #  for (i in 1:nPops){
      #    Popname <- getSubPop(pop=i, popsize=nPopsSamples[i],nChr=nChr, nPiecesPerChr=nPiecesPerChr,  recBTpieces=recBTpieces, nMrkOrMut=nLoci, minMAF=minMAF, tree=0)
      #    markers <- Popname$markers
      #    map <- Popname$map
      #    mapData <- makeMap(map=map, nLoci=nLoci, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, interactionMean=nEpiLoci)
      #    mapname <- paste("POP",i,"mapData",sep="")
      #    markername <- paste("POP",i,"founderHaps",sep="")
      #    final[[mapname]] <- mapData
      #    final[[markername]] <- markers
       # }
    }else{
      markers <- founderHaps$markers
      map <- founderHaps$map
      if (nrow(map) < nLoci) stop("Not enough loci in imported founder haplotypes for both markers and QTL")
      markers[is.na(markers)] <- 1 # No missing data
      mapData <- makeMap(map=map, nLoci=nLoci, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, interactionMean=nEpiLoci, qtlInfo=founderHaps$qtlInfo)
    }
    mapData$domModel <- domModel
    #return (final)
    return(list(mapData=mapData, founderHaps=markers))
  }#END defineSpecies.func
  
  if(is.null(loadData)){
    if (is.null(importFounderHap)){
    sims <- lapply(1:nSim, defineSpecies.func, nChr=nChr, lengthChr=lengthChr, effPopSize=effPopSize, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, nEpiLoci=nEpiLoci, domModel=domModel)
    }else{ # importFounderHap not NULL
      foundHap <- utils::read.table(file=paste(importFounderHap, ".hmp", sep=""))
      foundHap <- phasedHapMap2mat(foundHap)
      sims <- lapply(1:nSim, defineSpecies.func, nChr=nChr, lengthChr=lengthChr, effPopSize=effPopSize, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, nEpiLoci=nEpiLoci, founderHaps=foundHap, domModel=domModel)
    }
    save(sims, nSim, nCore, file=paste(saveDataFileName, ".RData", sep=""))
  }else{ # loadData not NULL
    load(paste(loadData, ".RData", sep=""))
    # Backward compatibility for versions that did not have domModel
    if (is.null(sims[[1]]$mapData$domModel)) sims <- lapply(sims, function(bsl){
      bsl$mapData$domModel <- "HetHom"
      return(bsl)
    })
  }
  # list of objects to remove before returning the environment
  toRemove <- c(setdiff(ls(), c("sims", "nSim", "nCore")), "toRemove")
  rm(list=toRemove)
  defineVariances(environment())
  return(environment())
}
