#-----------------------------------------------------------------------
#     VERSION HISTORY OF TestSimulation_AR.R
#-----------------------------------------------------------------------

# Version: v1.0.0
# Date: 08/23/2016
# Author: Yijia Liu
# Comment: Seperate calculation for groups using different test tool

# Version: v2.0.0
# Date: 10/04/2016
# Author: Yijia Liu
# Comment: Finalized for package

# Version: v2.0.1
# Date: 10/23/2016
# Author: Yijia Liu
# Comment: zero covariate case added, small calculation error fixed

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------
#' simulate AR dataset
#'
#' an internal function to simulates binary status and observations
#' with misclassification
#'
#' @param n Number of subjects in sample
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param max.obs Maximum number of obervation for each subject
#' @param x.cov matrix of correctly measured covariates
#' @param n.vv number of varying (time dependent) variable
#' @param subID ID for each subject
#' @param k Seed for simulation
#'
#' @import data.table
#' @return A dataframe including subject ID, binary status and observations
#'     with misclassification
simulate_AR <- function(n = sampleN, sens = 0.7, spec = 0.9, max.obs = 5, x.cov = x.cov,
  n.vv, subID = subID, k = 2836) {
  # Set the seed for simulation
  set.seed(k)
  # Define the parameters of the model
  #   beta for intercept
  b0 <- 0
  #   beta for y(t-1)
  b1 <- 1
  #   normally distributed error term
  e <- rnorm(n, 0, 0.1)

  if(dim(x.cov)[2] > 3) {
    # Convert x.cov to covariate matrix used in calculation
    xcov <- x.cov[ , c(3:(ncol(x.cov)-1))]
    X <- ConvertCov(xcov)
    #   beta for X's
    bx <- rep(1, dim(X)[2])
  }

  # Initiates y(with first column zero) and probs
  y <- matrix(0, n, max.obs + 1)
  prob <- matrix(0, n, max.obs)

  # Simmulate misclassified Y
  for (i in 1:max.obs) {
    # Calculate prob of real yt respect to y(t-1)
    if(dim(x.cov)[2] > 3) {
      Xt <- X[((i-1)*n+1):(i*n), ]
      prob[ , i] <- 1/(1 + exp(-b0 - b1 * y[ ,i] - Xt %*% bx - e))
    } else {
      prob[ , i] <- 1/(1 + exp(-b0 - b1 * y[ ,i]))
    }
    # Simulate yt
    y[ , i + 1] <- rbinom(n = n,size = 1,prob = prob[ , i])
  }

  # Get rid of the first zero column
  y <- y[ , c(2:(max.obs + 1))]

  # Generate random number of observations
  missingN <- sample(3 : max.obs, n, replace = TRUE)
  for (i in 1:n) {
    if(missingN[i] != max.obs) {
      y[i, c(missingN[i]:max.obs)] <- NA
    }
  }

  # Unroll y into large vector
  y.all <- c(sapply(c(1 : max.obs),function(m) {temp <- y[ ,m]} ))
  # Calculate the probablity of observations
  p.obs <- sens*(y.all==1) + (1-spec)*(y.all==0)
  # Simulate imperfect observations
  y.all.obs <- rbinom(n=length(p.obs),size=1,prob=p.obs)

  # yt1: status at (t - 1)
  yt1 <- y.all[1:((max.obs - 1)*n)]
  # yt1.obs: observed status value at (t - 1)
  yt1.obs <- y.all.obs[1:((max.obs - 1)*n)]
  # yt1: status at t
  yt <- y.all[(n+1):(max.obs*n)]
  # yt1.obs: observed status value at t
  yt.obs <- y.all.obs[(n+1):(max.obs*n)]
  # generate t
  t <- rep(1:(max.obs - 1),each=n)
  if(n.vv > 0) {
    # get time-varying variable at t = 0
    xv.0 <- x.cov[1:n, 3:(n.vv+2)]
    data <- as.data.frame(cbind(subID, yt, yt.obs, yt1, yt1.obs, y[ ,1], yt1.obs[1:n], t, xv.0))
  } else {
    data <- as.data.frame(cbind(subID, yt, yt.obs, yt1, yt1.obs, y[ ,1], yt1.obs[1:n], t))
  }
  colnames(data)[c(6,7)] <- c("y0","y0.obs")
  data <- merge(data, x.cov, by = c("subID", "t"), all.x = TRUE)
  # remove rows with NA's 
  data <- na.omit(data)
}

#' simulate covariates for AR(1) model
#'
#' Simulates correctly measured covariates
#'
#' @param n  Number of subjects in sample 
#' @param n.binom number of correctly measured binary covariates
#' @param n.cat number of correctly measured categorical covariates
#' @param n.cont number of correctly measured continuous covariates
#' @param n.vv number of varying (time dependent) variable
#' @param n.lev number of levels for categorical covariate
#' @param max.obs Maximum number of obervation for each subject
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
#' \item{t}{time}
#' @keywords covaraites
sim_covariates_AR <- function(n = n, n.binom = 1, n.cat = 1, n.cont = 1, n.vv = 1, n.lev = 3,
  max.obs = 5, grpIdx, subID, p = 0.6, p.cat = c(0.1,0.2,0.7), k = 2812) {
  # Set the seed for simulation
  set.seed(k)
  # initialize the matrix of covariates
  x.cov <- data.frame(row.names = (1:(n*max.obs)))
  n.cov <- n.vv + n.binom + n.cont + n.cat
  if(n.cov > 0) {
    for (i in 1: n.cov) {
      # simulate time-dependent continuous covariates
      if (i <= n.vv) {
        x.vv <- matrix(rnorm(n * max.obs), n, max.obs)
        x.temp <- matrix(c(sapply(c(1 : max.obs),function(m) {temp <- x.vv[ ,m]})))
      # simulate fixed binary covariates
      } else if (i <= n.vv + n.binom) {
        x.temp <- matrix(rbinom(n=n, size=1, prob=p))
      # simulate fixed continuous covariates
      } else if (i <= n.vv + n.binom + n.cont) {  
        x.temp <- matrix(rnorm(n),n,1)
      # simulate fixed categorical covariates   
      } else {
        x.temp <- matrix(sample(LETTERS[1:n.lev], n, replace=TRUE, prob=p.cat))
      }
      colnames(x.temp) <- paste("x", i, sep="")
      x.cov <- cbind(x.cov, x.temp)
    }
  }
  t <- rep(0:(max.obs - 1),each=n)
  x.cov <- cbind(subID, grpIdx, x.cov, t)
}

#' parse covariates
#'
#' an internal function to parses the order of covariates
#'
#' @param x.cov original covariates matrix
#' @param vv.index index of time-dependent variables
#'
#' @return A matrix of covariates in the order of [subID, time-dependent variable, 
#'   fixed varlable, t]
parse_cov <- function(x.cov, vv.index) {
  xv <- x.cov[vv.index == 1]
  xf <- x.cov[vv.index == 0]
  sub <- xf[ , 1:2]
  x.cov <- cbind(sub, xv, xf[ ,3:ncol(xf)])
  x.cov
}
