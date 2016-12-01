#-----------------------------------------------------------------------
#     VERSION HISTORY OF PerformanceEvaluation_CS.R
#-----------------------------------------------------------------------

# Version: v1.0.0
# Date: 04/27/2016
# Author: Yijia Liu
# Comment: Initial version. 

# Version: v1.0.2
# Date: 08/09/2016
# Author: Yijia Liu
# Comment: Testing different combinations of sens and spec

# Version: v2.0.1
# Date: 10/05/2016
# Author: Yijia Liu
# Comment: Empirical and Normal CI for both estimated beta and prevalence

# Version: v2.0.2
# Date: 10/09/2016
# Author: Yijia Liu
# Comment: Improved generalization of calculations 
#          Automate CI construction

# Version: v2.0.3
# Date: 10/23/2016
# Author: Yijia Liu
# Comment: zero covariate case added, small calculation error fixed

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------
#' Learning Performance Evaluation of CS model
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
PerformanceEvaluation_CS <- function(beta.est, sim = simN, cvg = FALSE) {
  # Calculate Mean
  Mean <- apply(beta.est,2,mean)
  # Calculate Bias
  Bias <- Mean - true
  # Calculate Standard Deviation
  ESE <- apply(beta.est,2,sd)
  # Calculate Percentage Bias
  Pct.Bias <- Bias / true * 100
  # Calculate Standard Bias
  Std.Bias <- Bias / ESE * 100
  # Calculate Mean Square Error
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
  cbind(anlys)
}


#' Bootstrapping CS dataset
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
BootstrapCI_CS <- function(bootN = 100, data, sim = 1) {
  # Set the seed pool for resampling
  seed.pool2 <- runif(bootN, 2322, 40246)
  # Initialize bootstrap confidence intervals if it is the first simulation
  if (sim == 1) {
    CI.beta <<- NULL
    # 95% Empirical Confidence Interval of beta
    CI.beta.emp <<- NULL
    # 95% Confidence Interval of prevalence
    CI.prev <<- NULL
    # 95% Empirical Confidence Interval of prevalence
    CI.prev.emp <<- NULL
  }
  # Initialize Bootstraped beta
  beta.boot <<- NULL
  # Initialize Bootstraped prevalence
  prev.boot <<- NULL
  # Resamples and estimates beta from resampled data for bootN times
  for (i in 1:bootN) {
    set.seed(seed.pool2[i])
    # Choose resampling indicies
    indices <- sample(c(1:dim(data)[1]), dim(data)[1], replace = TRUE, prob = NULL)
    data.boot <- data[indices, ]
    # Estimate bootstrapped beta
    beta.grp <- NULL
    prev.grp <- NULL
    for (grp in unique(data.boot$grpIdx)) {
      # Estimate beta
      grpdata <- data.boot[data.boot$grpIdx == grp, ]
      if(n.cov == 0) {
        x.cov <- NULL
      } else {
        x.cov <- grpdata[ , (ncol(grpdata)-n.cov+1):ncol(grpdata)]
      }
      beta <- CalculateBetaEM_CS(data = grpdata, x.cov = x.cov, sens = sens.grp[grp],
        spec = spec.grp[grp], sens.x = sens.x.grp[grp], spec.x = spec.x.grp[grp])
      beta.grp <- rbind(beta.grp, beta)
       # Estimate bootstrapped prevalence for four cases
      prev <- NULL
      for (y in c(0,1)) {
        for (x in c(0,1)) {
          prev.x <- 1/(1+exp(-beta[1] - beta[2]*x))
          ey <- prev.x * sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) / 
                (prev.x*sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) + 
                (1-prev.x) * (1-spec.grp[grp])^y * spec.grp[grp]^(1-y))
          prev <- c(prev, ey)
        }
      }
      prev.grp <- rbind(prev.grp, prev)
    }
    beta <- apply(beta.grp, 2, mean)
    prev <- apply(prev.grp, 2, mean)
   
    beta.boot <<- rbind(beta.boot, beta)
    prev.boot <<- rbind(prev.boot, prev)
  }
  Mean.beta <- apply(beta.boot, 2, mean)
  sd.beta <- apply(beta.boot, 2, sd)
  Mean.prev <- apply(prev.boot, 2, mean)
  sd.prev <- apply(prev.boot, 2, sd)
  # 95% Confidence Interval of beta
  CI.beta <<- rbind(CI.beta, rbind(Mean.beta - 1.96 * sd.beta, Mean.beta + 1.96 * sd.beta))
  # 95% Empirical Confidence Interval of beta
  CI.beta.emp <<- rbind(CI.beta.emp, rbind(apply(beta.boot, 2, quantile, 0.025), 
    apply(beta.boot, 2, quantile, 0.975)))
  # 95% Confidence Interval of prevalence
  CI.prev <<- rbind(CI.prev, rbind(Mean.prev - 1.96 * sd.prev, Mean.prev + 1.96 * sd.prev))
  # 95% Empirical Confidence Interval of prevalence
  CI.prev.emp <<- rbind(CI.prev.emp, rbind(apply(prev.boot, 2, quantile, 0.025), 
    apply(prev.boot, 2, quantile, 0.975)))
}

