#' Zero-inflated (semicontinuous) Data New Test Simulation
#'
#' This function highlights a new methodology for analysis of zero-inflated (semicontinous) data.
#' The proposed method separates and recombines the discrete and continuous parts
#' of the data in a new test specially made for zero-inflated data with higher expected power.
#' The function takes user input to generate datasets and outputs a graph of the
#' t statistic produced by the new test to hightlight and confirm the distribution of the new test.
#'
#' @param prob.x Probability of 0s in first data set
#' @param prob.y Probability of 0s in second data set
#' @param mu.x Normal Distribution mean for first data set
#' @param mu.y Normal Distribution mean for second data set
#' @param sd.x Normal Distribution standard deviation for first data set
#' @param sd.y Normal Distribution standard deviation for second data set
#' @param n1 Total number of data points for first data set
#' @param n2 Total number of data points for second data set
#'
#' @author Nicholas Sun, \email{ns2874@@columbia.edu}
#'
#' @export


########################################################################################
########################################################################################
SCtest_distribution <- function(prob.x, prob.y, mu.x, mu.y, sd.x, sd.y, n1, n2, simu.num)
{
  pval.vec <- matrix(NA, nrow=simu.num, ncol=1)
  for(i in 1:simu.num)
  {
  ### generate the data
  zeroLoc1 <- rbinom(n1, size=1, prob=prob.x)
  xnorm <- rnorm(n1, mean=mu.x, sd=sd.x)
  x <- zeroLoc1*xnorm
  zeroLoc2 <- rbinom(n2, size=1, prob=prob.y)
  ynorm <- rnorm(n2, mean=mu.y, sd=sd.y)
  y <- zeroLoc2*ynorm

  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  mu <- n1*n2/2
  sigma <- n1*n2*(n+1)/12

  ### permutate the sample
  dat <- c(x, y)
  ### the proposed test: combine the discrete and continuous parts
  n.x <- c(sum(x==0), sum(x!=0))
  n.y <- c(sum(y==0), sum(y!=0))
  temp.propose.part1 <- chisq.test(rbind(n.x, n.y))$p.value
  temp.propose.part2 <- t.test(x[x!=0], y[y!=0])$p.value
  temp.propose.test <- -2*log(temp.propose.part1)-2*log(temp.propose.part2)
  pval.vec[i,1] <- temp.propose.test
  }
  plot(density(pval.vec))
}





