#' Zero-inflated (semicontinuous) Simulation Data New Test
#'
#' This function highlights a new methodology for analysis of zero-inflated (semicontinous) data.
#' The proposed method separates and recombines the discrete and continuous parts
#' of the data in a new test specially made for zero-inflated data with higher expected power.
#' The function takes user input to generate datasets and outputs the p-value
#' in a two-sample difference of means test.
#'
#' @param prob.x Probability of 0s in first data set
#' @param prob.y Probability of 0s in second data set
#' @param mu.x Normal Distribution mean for first data set
#' @param mu.y Normal Distribution mean for second data set
#' @param sd.x Normal Distribution standard deviation for first data set
#' @param sd.y Normal Distribution standard deviation for second data set
#' @param n1 Total number of data points for first data set
#' @param n2 Total number of data points for second data set
#' @param B Number of permutations (should be over 10000)
#' @param alpha Alpha level. Defaults to 0.05
#' @return A 1x5 matrix of the power calculations for given input parameters
#'
#' @author Nicholas Sun, \email{ns2874@@columbia.edu}
#'
#' @export


########################################################################################
########################################################################################
SCtest_pvalue_generate <- function(prob.x, prob.y, mu.x, mu.y, sd.x, sd.y, n1, n2, B, alpha=0.05)
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
  perm.propose <- rep(NA, B)
  for(b in 1:B)
  {
    dat.perm <- dat[sample(1:n, n, replace=F)]
    x.perm <- dat.perm[1:n1]
    y.perm <- dat.perm[(n1+1):n]
    perm.t[b] <- t.test(x.perm, y.perm)$statistic
    perm.wil[b] <- (wilcox.test(x.perm, y.perm)$statistic-mu)/sqrt(sigma)
    n.x.perm <- c(sum(x.perm==0), sum(x.perm!=0))
    n.y.perm <- c(sum(y.perm==0), sum(y.perm!=0))
    perm.propose.part1 <- chisq.test(rbind(n.x.perm, n.y.perm))$p.value
    perm.propose.part2 <- t.test(x.perm[x.perm!=0], y.perm[y.perm!=0])$p.value
    perm.propose[b] <- -2*log(perm.propose.part1)-2*log(perm.propose.part2)

    pvalue <- mean(perm.propose > temp.propose.test)
  }
  ### calculate the average p-value
  return(pvalue)
}




