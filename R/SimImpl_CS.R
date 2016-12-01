#-----------------------------------------------------------------------
#     VERSION HISTORY OF SimImpl
#-----------------------------------------------------------------------

# Version: v1.0.0
# Date: 10/04/2016
# Author: Yijia Liu
# Comment: Simulation Implementation of Cross Sectional model

# Version: v1.0.1
# Date: 10/07/2016
# Author: Yijia Liu
# Comment: Group Simulation and Distinction added
#          Iteration number and convergence critiria parametrized
#          The option to save estimated beta added
#          Empirical and Normal CI for both estimated beta and prevalence added
#          CI coverage included in the performance summary

# Version: v1.0.2
# Date: 10/09/2016
# Author: Yijia Liu
# Comment: Prediction for specific patient based on given dataset added
#          (the option to construct CI included)
#          Simulation of dataset example added
#          Improved generalization of calculations 

# Version: v1.0.3
# Date: 10/23/2016
# Author: Yijia Liu
# Comment: zero covariate case added, small calculation error fixed

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------

#' Simulation Implementation of CS model
#'
#' Implement simulation of Cross-sectional model
#'
#' @param simN Number of dataset simulation
#' @param sampleN Number of subjects in sample
#' @param n.binom number of binomial covariates
#' @param n.cat number of categorial covariates
#' @param n.cont number of continuous covariates
#' @param n.lev number of levels for categorical covariate
#' @param n.grp number of testing technology
#' @param saveEst the option to save intermediate estimation
#' @param sensitivity Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param specificity Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param sensitivity.x Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
#' @param specificity.x Specificity.x = Prob(X1.obs = 0 | X1 = 0)
#' @param bootN The number of resampling to be performed
#' @param seed the root seed of simulation
#' @return writes simulation results into csv file
#' \item{b_AR}{evaluation of GLM estimates for correctly measured $X_1$}
#' \item{b.obs_AR}{evaluation of GLM estimates for misclassified $X_1$}
#' \item{b.EM_AR}{evaluation of EM estimates for misclassified $X_1$}
#' @examples
#' SimImpl_CS()
#' @keywords simulation
#' @keywords implementation
#' @export
SimImpl_CS <- function(simN = 5, sampleN = 500, n.binom = 1, n.cat = 1, n.cont = 1, n.lev = 3,
  n.grp = 3, saveEst = TRUE, sensitivity = c(0.8,0.85,0.9), specificity = c(0.9,0.85,0.8), 
  sensitivity.x = c(0.8,0.85,0.9), specificity.x = c(0.9,0.85,0.8), seed = 616, bootN = 2) {
  # Make sure the dimension of sensitivity/specificity is consistent with n.grp
  # Stop simulation and print error message if not
  if(length(sensitivity) != n.grp | length(sensitivity.x) != n.grp |
    length(specificity) != n.grp | length(specificity.x) != n.grp ) {
    stop("the dimension of sensitivity/specificity does not conform with number of group")
  }

  options(warn = -1)

  set.seed(seed)  # Set the root seed
  seed.pool <- runif(simN+1,2321,40245)  # Set the seed pool for data generation
  n.cov <<- n.binom + n.cat + n.cont  # number of covariates
  n <<- sampleN
  # Sensitivity = Prob(Y = 1 | Y.obs = 1)
  sens.grp <<- sensitivity
  # Specificity = Prob(Y.obs = 0 | Y = 0)
  spec.grp <<- specificity
  # Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
  sens.x.grp <<- sensitivity.x
  # Specificity.x = Prob(X1.obs = 0 | X1 = 0)
  spec.x.grp <<- specificity.x

  # Initialize betas'
  b <- NULL
  b.obs <- NULL
  b.em <- NULL

  # Perform simN times of data generation & evaluation
  for(i in 1:simN) {
    data <- NULL
    for (grp in 1:n.grp) {
      # generate subject ID
      subID <- sample((sampleN*grp+1):(sampleN*(grp+1)), sampleN)
      # Simulates correctly measured covariates
      x.cov <- sim_covariates_CS(
        # Number of subjects in sample
        n = sampleN,
        # Number of binomial covariates
        n.binom = n.binom,
        # Number of categorial covariates
        n.cat = n.cat,
        # Number of continuous covariates
        n.cont = n.cont,
        # n.lev: number of levels for categorical covariate
        n.lev = n.lev,
        # group index for each subject
        grpIdx = grp,
        # subject ID
        subID = subID, 
        # Set the seed for simulation
        k = seed.pool[i + 1] + grp
      )   # End of simulation

      grpdata <- simulate_CS(
        # Number of subjects in sample
        n = sampleN,
        # matrix of correctly measured covariates
        x.cov = x.cov,
        # sensitivity of each subject in the data
        sub.sens = sens.grp[grp],
        # specificity of each subject in the data
        sub.spec = spec.grp[grp],
        # sensitivity.x of each subject in the data
        sub.sens.x = sens.x.grp[grp],
        # specificity.x of each subject in the data
        sub.spec.x = spec.x.grp[grp],
        # Set the seed for simulation
        k = seed.pool[i] + grp
      )
      data <- rbind(data, grpdata)
    }
    # shuffle the order of tuples in dataset
    data <- data[sample(nrow(data)),]
    
    # Estimate beta's from dataset
    for (grp in unique(data$grpIdx)) {
      grpbeta <- data.frame(cbind(calc_grp_CS(grpdata = data[data$grpIdx == grp, ],
        grp = grp), grp))
      b <- rbind(b, grpbeta[1, ])
      b.obs <- rbind(b.obs, grpbeta[2, ])
      b.em <- rbind(b.em, grpbeta[3, ])
    }
    
    #  Bootstrap Confidence Interval of beta and prevalence for resampled data
    if (bootN != 0) {
      BootstrapCI_CS(bootN, data, i)
    }
  }

  if(saveEst == TRUE) {
    beta.EM <<- b.em
    beta.obs <<- b.obs
    beta.true <<- b
  }  
  ###### Simulation Evaluation #######
  true <<- rep(1,ncol(b))
  true[1] <<- 0

  # Save the simulation results
  b_CS <<- round(PerformanceEvaluation_CS(b), 4)
  b.obs_CS <<- round(PerformanceEvaluation_CS(b.obs), 4)
  b.EM_CS <<- round(PerformanceEvaluation_CS(beta.EM, sim = simN, cvg = (bootN > 0)), 4)
}

#' Predict STI probability
#'
#' Predict STI infection probability based on provided dataset and patient info
#'
#' @param data dataset for beta estimation
#' @param patInfo patient's diagonistic information
#' @param n.cov number of covariates
#' @param n.grp n.grp: number of testing technology
#' @param sensitivity Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param specificity Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param sensitivity.x Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
#' @param specificity.x Specificity.x = Prob(X1.obs = 0 | X1 = 0)
#' @param getBeta the option to save beta estimated from given dataset
#' @param CI the option to get 95\% confidence interval
#' @param bootN the number of BootStrapping used to construct CI
#'
#' @return pecific patient's infection probability
#' \item{beta.est}{estimated beta from provided daaset}
#' \item{pat.CI.prev}{patient's 95\% CI of prevalence}
#' \item{pat.CI.prev.emp}{patient's 95\% empirical CI of prevalence}
#' @examples
#' data <- simExample_CS(sampleN=500)
#' patInfo <- data[5, ]
#' predict_CS(data,patInfo,3,3,c(0.8,0.85,0.9),c(0.9,0.85,0.8),c(0.8,0.85,0.9),c(0.9,0.85,0.8))
#' beta.est
#' @keywords prediction
#' @keywords prevalence
#' @export
predict_CS <- function (data, patInfo, n.cov, n.grp, sensitivity, specificity, sensitivity.x,
  specificity.x, getBeta = TRUE, CI = FALSE, bootN = 10) {
  if(length(sensitivity) != n.grp | length(sensitivity.x) != n.grp |
    length(specificity) != n.grp | length(specificity.x) != n.grp ) {
    stop("the dimension of sensitivity/specificity does not conform with number of group")
  }
  n.cov <<- n.cov
  # detect number of subjects in the given data
  n <<- dim(data)[1]
  # get patient's group
  grp <- patInfo$grp
  # get patient's observed infection status
  y <- patInfo$y.obs
  # get patient's observed x1
  x <- patInfo$x1.obs
  # Sensitivity = Prob(Y = 1 | Y.obs = 1)
  sens.grp <<- sensitivity
  # Specificity = Prob(Y.obs = 0 | Y = 0)
  spec.grp <<- specificity
  # Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
  sens.x.grp <<- sensitivity.x
  # Specificity.x = Prob(X1.obs = 0 | X1 = 0)
  spec.x.grp <<- specificity.x
  if(n.cov != 0) {
    # Extract patient's covariates
    x.cov <- patInfo[ , (ncol(patInfo)-n.cov+1):ncol(patInfo)]
    # Convert x.cov into matrix
    Xcov <- as.numeric(ConvertCov(x.cov))
  }
  
  # get EM estimation of beta for each group
  # Initialize estimation
  beta.EM <- NULL
  for (grp in unique(data$grpIdx)) {
    grpdata <- data[data$grpIdx == grp, ]
    if(n.cov == 0) {
      x.cov <- NULL
    } else {
      x.cov <- grpdata[ , (ncol(grpdata)-n.cov+1):ncol(grpdata)]
    }
    b <- CalculateBetaEM_CS(data = grpdata, x.cov = x.cov, sens = sens.grp[grp],
      spec = spec.grp[grp], sens.x = sens.x.grp[grp], spec.x = spec.x.grp[grp])
    beta.EM <- rbind(beta.EM, b)
  }
  # average beta from different groups
  beta <- apply(beta.EM, 2, mean)
  if(getBeta == TRUE) beta.est <<- beta
  # Calculate Empirical and Statistical 95% CI
  if(CI  == TRUE) {
    BootstrapCI_CS(bootN,data)
    if(y == 0) {
      if(x == 0) {
        case <- 1
      } else {
        case <- 2
      }
    } else {
      if(x == 0) {
        case <- 3
      } else {
        case <- 4
      }
    }
    if (n.cov != 0) {
      adj <- Xcov %*% beta[-c(1,2)]
      CI.prev <<- 1 / (1 + exp(log(1 / CI.prev[ , case] - 1) - adj))
      CI.prev.emp <<- 1 / (1 + exp(log(1 / CI.prev.emp[ , case] - 1) - adj))
    } else {
      CI.prev <<- CI.prev[ , case]
      CI.prev.emp <<- CI.prev.emp[ , case]
    }
  }
  # Calculate prevalence
  if(n.cov == 0) {
    prev.x <- 1 / (1 + exp(-beta[1] - beta[2] * x))
  } else {
    prev.x <- 1 / (1 + exp(-beta[1] - beta[2] * x - Xcov %*% beta[-c(1,2)]))
  }
  prev <- prev.x * sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) / 
          (prev.x * sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) + 
          (1-prev.x) * (1-spec.grp[grp])^y * spec.grp[grp]^(1-y))
  if(prev > 1) prev <- 1
  if(prev < 0) prev <- 0
  prev
}

#' Simulate an example of CS dataset
#'
#' Simulation an example for the users to prepare their data structure
#'
#' @param sampleN number of tuples in each group
#' @param sensitivity Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param specificity Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param sensitivity.x Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
#' @param specificity.x Specificity.x = Prob(X1.obs = 0 | X1 = 0)
#' @param n.grp number of testing technology
#' @param n.binom number of binomial covariates
#' @param n.cat number of categorial covariates
#' @param n.cont number of continuous covariates
#' @param seed the root seed of simulation
#'
#' @return an example of dataset
#' @examples
#' data <- simExample_CS()
#' head(data)
#' @keywords Example 
#' @keywords simulate 
#' @export
simExample_CS <- function(sampleN = 100, sensitivity = c(0.8,0.85,0.9),
  specificity = c(0.9,0.85,0.8), sensitivity.x = c(0.8,0.85,0.9),
  specificity.x = c(0.9,0.85,0.8), n.grp = 3, n.binom = 1, n.cat = 1,
  n.cont = 1, seed = 616) {
  set.seed(seed)  # Set the root seed
  # number of covariates
  n.cov <- n.binom + n.cat + n.cont
  n <- sampleN
  # Sensitivity = Prob(Y = 1 | Y.obs = 1)
  sens.grp <- sensitivity
  # Specificity = Prob(Y.obs = 0 | Y = 0)
  spec.grp <- specificity
  # Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
  sens.x.grp <- sensitivity.x
  # Specificity.x = Prob(X1.obs = 0 | X1 = 0)
  spec.x.grp <- specificity.x
  data <- NULL
  for (grp in 1:n.grp) {
    # generate subject ID
    subID <- sample((sampleN*grp+1):(sampleN*(grp+1)), sampleN)
    # Simulates covariates
    x.cov <- sim_covariates_CS(
      # Number of subjects in sample
      n = sampleN,
      # Number of binomial covariates
      n.binom = n.binom,
      # Number of categorial covariates
      n.cat = n.cat,
      # Number of continuous covariates
      n.cont = n.cont,
      # group index for each subject
      grpIdx = grp,
      # subject ID
      subID = subID, 
      # Set the seed for simulation
      k = seed + grp
    )   # End of simulation

    grpdata <- simulate_CS(
      # Number of subjects in sample
      n = sampleN,
      # matrix of correctly measured covariates
      x.cov = x.cov,
      # sensitivity of each subject in the data
      sub.sens = sens.grp[grp],
      # specificity of each subject in the data
      sub.spec = spec.grp[grp],
      # sensitivity.x of each subject in the data
      sub.sens.x = sens.x.grp[grp],
      # specificity.x of each subject in the data
      sub.spec.x = spec.x.grp[grp],
      # Set the seed for simulation
      k = seed + grp + 1
    )
    data <- rbind(data, grpdata)
  }
  # shuffle the order of tuples in dataset
  data <- data[sample(nrow(data)),]
  # Subtract informations provided by user
  if(n.cov == 0) {
    data.user <- data[, c(1,3,5,6)]
  } else {
    data.user <- cbind(data[, c(1,3,5,6)], data[ , (ncol(data)-n.cov+1):ncol(data)])
  }
  data.user
}


