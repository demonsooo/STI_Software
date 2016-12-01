#-----------------------------------------------------------------------
#     VERSION HISTORY OF PerformanceEvaluation_AR.R
#-----------------------------------------------------------------------
​
# Version: v1.0.0
# Date: 04/27/2016
# Author: Yijia Liu
# Comment: Initial version. 

# Version: v1.0.2
# Date: 08/09/2016
# Author: Yijia Liu
# Comment: Testing different combinations of sensitivity and specificity

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------

#' Learning Performance Evaluation of AR(1) model
#'
#' An internal function to Evaluate the performance of statistical methods for
#' betas estimated by EM algorithm and GLM function.
#'   Performance measures including:
#'     Bias
#'     Mean
#'     Standard Deviation
#'     Percentage Bias
#'     Standard Bias
#'     Mean Square Error
#'
#' @param beta.est A matrix of betas estimated from different datasets
#' @param sim the number of simulation performed
#' @param cvg the option to find 95\% CI coverage
#' @return A dataframe of the performance of statistical methods for given set of betas
#' @keywords Evaluation
PerformanceEvaluation_AR <- function(beta.est, sim = simN, cvg = FALSE) {
  Mean <- apply(beta.est,2,mean)
  Bias <- Mean - true
  ESE <- apply(beta.est,2,sd)
  Pct.Bias <- Bias / true * 100
  Std.Bias <- Bias / ESE * 100
  MSE <- Bias ^2 + ESE ^2

  # Combine the results into a dataframe
  if(cvg == FALSE) {
    anlys <- as.data.frame(rbind(Mean, Bias, ESE, Pct.Bias, Std.Bias, MSE))
  } else {
    # Calculate beta coverage
    beta.cvg <- apply(t(sapply(1:sim, function(j) {
      (CI.beta[2 * j - 1, ] <= true) & (CI.beta[2 * j, ] >= true)})),2,sum) / sim
    # Calculate empirical beta coverage
    beta.cvg.emp <- apply(t(sapply(1:sim, function(j) {
      (CI.beta.emp[2 * j - 1, ] <= true) & (CI.beta.emp[2 * j, ] >= true)})),2,sum) / sim
    anlys <- as.data.frame(rbind(Mean, Bias, ESE, Pct.Bias, Std.Bias, MSE, beta.cvg, beta.cvg.emp))
  }
  anlys
}

#' Bootstrapping AR(1) dataset
#'
#' An internal function to compute the confidence interval of resampled data using bootstrap methods
#'
#' @param bootN The number of resampling to be performed
#' @param data the dataset to be evaluated
#' @param sim the number of simulation performed
#' @return The confidence interval of resampled data
#' \item{CI.beta}{ 95\% Confidence Interval of beta }
#' \item{CI.beta.emp}{ 95\% Empirical Confidence Interval of beta }
#' \item{CI.prev}{ 95\% Confidence Interval of prevalence }
#' \item{CI.prev.emp}{ 95\% Empirical Confidence Interval of prevalence }
#' 
#' @keywords Bootstrap
BootstrapCI_AR <- function(bootN = 100, data, sim = 1) {
  seed.pool2 <- runif(bootN, 2321, 40245)
  if (sim == 1) {
    CI.beta <<- NULL
    CI.beta.emp <<- NULL
    CI.prev <<- NULL
    CI.prev.emp <<- NULL
  }
  beta.boot <<- NULL
  prev.boot <<- NULL
  for (i in 1:bootN) {
    set.seed(seed.pool2[i])
    indices <- sample(c(1:dim(data)[1]), dim(data)[1], replace = TRUE, prob = NULL)
    data.boot <- data[indices, ]
    beta.grp <- NULL
    for (grp in unique(data.boot$grpIdx)) {
      grpdata <- data.boot[data.boot$grpIdx == grp, ]
      if(n.cov <= 0) {
        x.cov <- NULL
      } else {
        x.cov <- grpdata[ , (ncol(grpdata)-n.cov+1):ncol(grpdata)]
      }
      beta <- CalculateBetaEM_AR(data = grpdata, x.cov = x.cov, sens = sens.grp[grp],
                  spec = spec.grp[grp])
      beta.grp <- rbind(beta.grp, beta)
    }
    beta <- apply(beta.grp, 2, mean)
    prev <- NULL
    for (y in c(0,1)) {
      for (y1 in c(0,1)) {
        prev.y1 <- 1/(1+exp(-beta[1] - beta[2]*y1))
        ey <- prev.y1 * sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) / 
              (prev.y1*sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) + 
              (1-prev.y1) * (1-spec.grp[grp])^y * spec.grp[grp]^(1-y))
        prev <- c(prev, ey)
      }
    }
    beta.boot <<- rbind(beta.boot, beta)
    prev.boot <<- rbind(prev.boot, prev)
  }
  Mean.beta <- apply(beta.boot, 2, mean)
  sd.beta <- apply(beta.boot, 2, sd)
  Mean.prev <- apply(prev.boot, 2, mean)
  sd.prev <- apply(prev.boot, 2, sd)
  CI.beta <<- rbind(CI.beta, rbind(Mean.beta - 1.96 * sd.beta, Mean.beta + 1.96 * sd.beta))
  CI.beta.emp <<- rbind(CI.beta.emp, rbind(apply(beta.boot, 2, quantile, 0.025), 
    apply(beta.boot, 2, quantile, 0.975)))
  CI.prev <<- rbind(CI.prev, rbind(Mean.prev - 1.96 * sd.prev, Mean.prev + 1.96 * sd.prev))
  CI.prev.emp <<- rbind(CI.prev.emp, rbind(apply(prev.boot, 2, quantile, 0.025), 
    apply(prev.boot, 2, quantile, 0.975)))
}
