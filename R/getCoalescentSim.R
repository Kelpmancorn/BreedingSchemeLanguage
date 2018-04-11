#'getCoalescentSim
#'
#'@useDynLib BreedingSchemeLanguage coalescentSim
#'
#'@param nPopsSamples -pop[2] (GENOME)
#'@param effPopSize -N (GENOME)
#'@param nChr -c (GENOME)
#'@param nPiecesPerChr -pieces (GENOME)
#'@param recBTpieces -rec (GENOME)
#'@param nMrkOrMut parameter to calculate -s (GENOME)
#'@param minMAF -maf (GENOME)
#'@param seed -seed (GENOME)
#'@param tree -tree (GENOME)
#'
#getCoalescentSim <- function(nPopsSamples=NULL, effPopSize=100, nChr=7, nPiecesPerChr=15000, recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, seed=as.integer(Sys.time()), tree=0){

##number of population nPops is added to simulate subpopulations. len parameter is 10000 by default so it is not included in

getCoalescentSim <- function(nPops=12, nPopsSamples=rep(12,12), effPopSize=100, nChr=31, nPiecesPerChr=180,  recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, seed=as.integer(Sys.time()), tree=0){ ##effective population size should be larger than the sample size

##try to set up some trial parameters 
#gsub(",",";",c("nPops=12, nPopsSamples=rep(100,12), effPopSize=10, nChr=31, nPiecesPerChr=18, len=100, recBTpieces=0.0001, nMrkOrMut=10, minMAF=0.01, seed=as.integer(Sys.time()), tree=0"))
  
systemCall <- rep(NA, 21)
#systemCall <- rep(NA, 10)
#nPopsSamples <- effPopSize
systemCall[1] <- nPops     # pop1
systemCall[c(2:13)] <- nPopsSamples     #pop2
systemCall[14] <- effPopSize     # N
systemCall[15] <- nChr     #c
systemCall[16] <- nPiecesPerChr     # pieces
#systemCall[17] <- len     # len
systemCall[17] <- recBTpieces     # rec
systemCall[18] <- minMAF     #maf
systemCall[19] <- seed     # seed
systemCall[20] <- tree     # tree
#systemCall[21] <- round(8 * nMrkOrMut / nChr)     # s  ##why?????????????
systemCall[21] <- -1     # s  


doGenome <- function(call){
  .Call('coalescentSim',
        par1 = call[1],
        par2 = call[c(2:13)],
        par3 = as.character(call[14]),
        par4 = call[15],
        par5 = call[16],
        #par6 = call[17],
        par6 = as.character(call[17]),
        par7 = call[18],
        par8 = call[19],
        par9 = call[20],
        par10 = call[21])
}

if(file.exists("genomeOUT.txt")){
  file.remove("genomeOUT.txt")
}
sink("genomeOUT.txt")
doGenome(systemCall)
sink()
#}

##try out simulations
#getCoalescentSim(nPops=12, nPopsSamples=rep(100,12), effPopSize=10000, nChr=31, nPiecesPerChr=180,  recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, seed=1521156195, tree=0)##the computation time is a lot compared to running in linux

###build functions to return subpop###
#getSubPop <- function(pop=1, popsize=100,nChr=31, nPiecesPerChr=180,  recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, tree=0){
genOut <- readLines("genomeOUT.txt", n=-1)
##We don't want to remove the GENOME output
#file.remove("genomeOUT.txt")
strtPos <- grep("SNP genetic position", genOut)  ##get the linenumbers of genetic positions
map <- NULL

##from genetic to physical distance????????????????????
for (chr in 1:nChr){
  posVec <- round(as.numeric(strsplit(genOut[strtPos[chr]+1], split=" ")[[1]]) * 100 * (nPiecesPerChr - 1) * recBTpieces, 6)
  if (length(posVec > 0)) map <- rbind(map, cbind(chr, posVec))
}
colnames(map) <- c("Chr", "Pos")
strtInd <- NULL
for (pop in 1:length(nPopsSamples)){
  strtInd <- c(strtInd, grep(paste("POP", pop, ":", sep=""), genOut))
}
#strtInd <- c(grep(paste("POP", pop, ":", sep=""), genOut))
#genOut <- substring(genOut[strtInd], 7) ##get the line number of individuals, if populaiton ID >=10, it is not 7 any more, this will result in NA in the dataset when caculating freq
trim.leading <- function (x)  sub("^\\s+", "", x)
genOut <- trim.leading(substring(genOut[strtInd], 7))
nSNP <- nchar(genOut[1])
#markers <- t(array(as.numeric(unlist(strsplit(genOut, split=""))), c(nSNP, sum(popsize))))
markers <- t(array(as.numeric(unlist(strsplit(genOut, split=""))), c(nSNP, sum(nPopsSamples))))
freq <- colMeans(markers)   
##remove maf < minMAF
maf <- abs((freq > 0.5) - freq)
markers <- markers[, maf >= minMAF]
map <- map[maf >= minMAF,]
maf <- maf[maf >= minMAF]


##number of markers is not the limiting factor of speed, we just sample that numbers of markers from existed
if (length(nMrkOrMut) == 1){
  if (ncol(markers) < nMrkOrMut){
    print("Warning! Settings such that fewer markers simulated than needed")
    print(paste("There were", ncol(markers), "markers"))
  } else{
    keepMrk <- sort(sample(ncol(markers), nMrkOrMut))
    markers <- markers[, keepMrk]
    map <- map[keepMrk, ]
    maf <- maf[keepMrk]
  }
}
nMrk <- ncol(markers)
# Ensure that no two markers have _exactly_ the same position
for (chr in 1:nChr){
  mrkThisChr <- map[,"Chr"] == chr
  uniqPos <- unique(map[mrkThisChr,"Pos"])
  for (pos in uniqPos){
    mrkAtPos <- which(mrkThisChr & map[,"Pos"] == pos)
    if (length(mrkAtPos) > 1) map[mrkAtPos,"Pos"] <- map[mrkAtPos,"Pos"] + sort(round(stats::runif(length(mrkAtPos), 0, 100 * recBTpieces), 6))
  }
}
map <- as.data.frame(map)
return(list(markers=markers, map=map))
}

#Pop1 <- getSubPop(pop=1, nChr=31, nPiecesPerChr=180,  recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, tree=0)

#linux command
#./genome-linux-64bit  -pop 12 100 100 100 100 100 100 100 100 100 100 100 100 -N 100 -c 31 -pieces 180 -len 10000 -s -1 -rec 0.0001 -seed 1521156195 -maf 0.01

