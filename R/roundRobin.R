#'
#'@param n the number of pairs
#'
roundRobin <-function(n=20){
rounds <- list()
teams <- 1:n
r <- n-1
for( i in 1:r){
  round <- 
    data.frame(
      round = i,
      team1 = teams[1:(n/2)], 
      team2 = rev(teams)[1:(n/2)])
  rounds[[i]] <- round
  teams <- c( teams[1],  last(teams), head(teams[-1],-1) ) 
}

rr <- dplyr::bind_rows(rounds)
rr <- as.matrix(rr)[,2:3]
colnames(rr) <- NULL
return(rr)
}