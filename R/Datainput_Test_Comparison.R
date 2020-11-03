#' Zero-inflated Data Test Comparison
#'
#' This function compares the power of multiple parametric and non-parametric tests
#' (ordinary t-test, wilcoxon test and permutated t-test, wilcoxon test)
#' against a new type of method of semicontinuous data analysis that combines
#' discrete and continuous data in a unique fashion for enhanced analysis.
#' Simulation type function utilized for proof of concept which can also
#' be applied to real data.
#'
#'
#' @param datax First set of data, should take form of numeric vector or dataframe
#' @param datay Second set of data, should take form of numeric vector or dataframe
#' @param B Number of permutations (should be over 10000)
#' @param alpha Alpha level. Defaults to 0.05
#'
#' @return A 1x5 matrix of the p-value (t-test, wilcoxon, permutated t-test, permutated wilcoxon, new test)
#'
#' @author Nicholas Sun, \email{ns2874@@columbia.edu}
#'
#' @export


########################################################################################
###### Compare the test for  semicontinuous data #######################################
### Compared Tests: ordinary t-test, wilcoxon test and permutated t-test, wilcoxon test
########################################################################################
SC_alltest_pvalue <- function(datax, datay, B, alpha=0.05)
{
    i=1
    x <- datax
    y <- datay

    n1 <- length(x)
    n2 <- length(y)
    n <- n1 + n2
    mu <- n1*n2/2
    sigma <- n1*n2*(n+1)/12

    ### Test 1: two sample t-test
    pval.vec[i,1] <- t.test(x, y)$p.value

    ### Test 2: wilcoxon test
    pval.vec[i,2] <- wilcox.test(x, y)$p.value

    ### permutate the sample
    dat <- c(x, y)
    temp.t.test <- t.test(x, y)$statistic
    temp.wil.test <- (wilcox.test(x, y)$statistic-mu)/sqrt(sigma)
    ### the proposed test: combine the discrete and continuous parts
    n.x <- c(sum(x==0), sum(x!=0))
    n.y <- c(sum(y==0), sum(y!=0))
    temp.propose.part1 <- chisq.test(rbind(n.x, n.y))$p.value
    temp.propose.part2 <- t.test(x[x!=0], y[y!=0])$p.value
    temp.propose.test <- -2*log(temp.propose.part1)-2*log(temp.propose.part2)
    perm.t <- rep(NA, B)
    perm.wil <- rep(NA, B)
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
    }
    ### Test 3: permutation t-test
    pval.vec[i,3] <- mean(perm.t < -abs(temp.t.test) | perm.t > abs(temp.t.test))
    ### Test 4: permutation wilcoxon test
    pval.vec[i,4] <-  mean(perm.wil < -abs(temp.wil.test) | perm.wil > abs(temp.wil.test))
    ### Test 5: the proposed test
    pval.vec[i,5] <- mean(perm.propose > temp.propose.test)

  return(pval.vec)
}


