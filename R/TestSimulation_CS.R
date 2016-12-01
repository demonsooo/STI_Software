#-----------------------------------------------------------------------
#     VERSION HISTORY OF TestSimulation_n_cov.R
#-----------------------------------------------------------------------

# Version: v1.0.0
# Date: 04/19/2016
# Author: Yijia Liu
# Comment: Initial version. 

# Version: v1.0.1
# Date: 04/26/2016
# Author: Yijia Liu
# Comment: Categorical variable / dummy matrix added

# Version: v1.0.2
# Date: 05/06/2016
# Author: Yijia Liu
# Comment: Allowing specify "n.binom", "n.cont" and "n.cat" on data 
#          simulation

# Version: v2.0.1
# Date: 10/04/2016
# Author: Yijia Liu
# Comment: Finalized for package

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
#' Simulate CS dataset
#' 
#' an internal function to simulate binary status and observations with misclassification
#'
#' @param n Number of subjects in sample
#' @param p probability of $X_1$ to be 1
#' @param x.cov matrix of correctly measured covariates
#' @param sub.sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param sub.spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param sub.sens.x Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
#' @param sub.spec.x Specificity.x = Prob(X1.obs = 0 | X1 = 0)
#' @param k Seed for simulation
#'
#' @import data.table
#' @return A dataframe including subject ID, binary status and observations
#'     with misclassification
simulate_CS <- function(n = sampleN, p = 0.5, x.cov = x.cov, sub.sens, sub.spec,
  sub.sens.x, sub.spec.x, k = 2836) {
  # Set the seed for simulation
  set.seed(k)

  subID <- x.cov$subID
  # number of covariates
  n.cov <- dim(x.cov)[2] - 2
  # Simulate misclassified covariates
  #   x1: misclassified variable
  x1 <- rbinom(n=n, size=1, prob=p)

  # Convert x.cov to covariate matrix used in calculation
  if(n.cov == 0) {
    X <- cbind(rep(1,n), x1)
  } else {
    xcov <- x.cov[ , c(3:ncol(x.cov))]
    X <- cbind(rep(1,n), x1, ConvertCov(xcov))
  }

  # Define the parameters of the model
  b <- rep(1, dim(X)[2])
  b[1] <- 0
  e <- rnorm(n, 0, 0.1)

  # Simulate infection status
  #   prob.y = P(Y = 1 | X1, X2)
  prob.y <- 1 / (1 + exp(- X %*% b - e))
  y <- rbinom(n = n, size = 1, prob = prob.y)

  # Simulate misclassified X1
  #   prob.x1.obs = P(X1.obs = 1 | X1)
  prob.x1.obs <- sub.sens.x * (x1 == 1) + (1 - sub.spec.x) * (x1 == 0)
  x1.obs <- rbinom(n = n, size = 1, prob = prob.x1.obs)

  # Simulate misclassified test result
  #   prob.y.obs = P(Y.obs = 1 | Y)
  prob.y.obs <- sub.sens * (y == 1) + (1 - sub.spec) * (y == 0)
  y.obs <- rbinom(n = n, size = 1, prob = prob.y.obs)

  ## Return data frame
  data <- as.data.frame(cbind(subID, y, y.obs, x1, x1.obs))
  data <- merge(data, x.cov, by = c("subID"), all.x = TRUE)
  data
}

#' create dummy matrix
#' Create dummy matrix for categorical variables
#'
#' @param x.cov matrix of correctly measured covariates
#' @param cat.index index of categorical variable(s)
#'
#' @return A dummy matrix for categorical variables
create_dummy_mat <- function(x.cov = x.cov, cat.index = is.cat) {
  # initialize dummy variable matrix
  x.dummy <- NULL
  for(i in cat.index) {
    X <- as.factor(x.cov[ ,i])
    x.temp <- model.matrix(~X)
    x.temp <- x.temp[ ,-1]
    if(class(x.temp) == "numeric") {
      x.dummy <- c(x.dummy, x.temp)
    } else {
      x.dummy <- cbind(x.dummy, x.temp)
    }  
  }
  x.dummy
}

#' simulate covariates for CS model
#'
#' an internal function to simulate correctly measured covariates
#'
#' @param n  Number of subjects in sample
#' @param n.binom number of correctly measured binary covariates
#' @param n.cat number of correctly measured categorical covariates
#' @param n.cont number of correctly measured continuous covariates
#' @param n.lev number of levels for categorical covariate
#' @param grpIdx group index for each subject
#' @param subID ID for each subject
#' @param p probability of binary covariate to be 1
#' @param p.cat probability of different levels of categorical covariate
#' @param k Seed for simulation
#'
#' @return A dataframe including correctly measured covariates
#' \item{subID}{ID for each subject}
#' \item{grpIdx}{group index for each subject}
#' \item{$X_n$}{correctly measured covariates}
#' @keywords covaraites
sim_covariates_CS <- function(n = sampleN, n.binom = 1, n.cat = 1, n.cont = 1, n.lev = 3,
  grpIdx, subID, p = 0.6, p.cat = c(0.1,0.2,0.7), k = 2812) {
  # Set the seed for simulation
  set.seed(k)
  # initialize the matrix of covariates
  x.cov <- data.frame(row.names = (1:n))
  n.cov <- n.binom + n.cont + n.cat
  if(n.cov > 0) {
    for (i in 1: n.cov) {
      # simulate binary covariates
      if (i <= n.binom) {
        x.temp <- matrix(rbinom(n=n, size=1, prob=p))
      # simulate continuous covariates
      } else if (i <= n.binom + n.cont) {  
        x.temp <- matrix(rnorm(n),n,1)
      # simulate categorical covariates   
      } else {
        x.temp <- matrix(sample(LETTERS[1:n.lev], n, replace=TRUE, prob=p.cat))
      }
      colnames(x.temp) <- paste("x", i+1, sep="")
      x.cov <- cbind(x.cov, x.temp)
    }
  }
  x.cov <- cbind(subID, grpIdx, x.cov)
  x.cov
}
