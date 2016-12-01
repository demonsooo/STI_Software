#-----------------------------------------------------------------------
#     VERSION HISTORY OF SimImpl_AR
#-----------------------------------------------------------------------

# Version: v1.0.0
# Date: 09/21/2016
# Author: Yijia Liu
# Comment: Simulation Implementation of AR(1) model

# Version: v1.0.1
# Date: 10/12/2016
# Author: Yijia Liu
# Comment: Prediction for specific patient based on given dataset added
#          (the option to construct CI included)
#          Simulation of dataset example added
#          Improved generalization of calculations

# Version: v1.0.2
# Date: 10/23/2016
# Author: Yijia Liu
# Comment: zero covariate case added, small calculation error fixed

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------

#' Simulation Implementation
#' 
#' Implement simulation of AR(1) model
#'
#' @param simN Number of dataset simulation
#' @param sampleN Number of subjects in sample
#' @param max.obs maximum number of obervation for each subject
#' @param n.binom number of binomial covariates
#' @param n.cat number of categorial covariates
#' @param n.cont number of continuous covariates
#' @param n.v number of varying (time dependent) variable
#' @param n.lev number of levels for categorical covariate
#' @param n.grp number of testing technology
#' @param saveEst the option to save intermediate estimation
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param bootN The number of resampling to be performed
#' @param seed the root seed of simulation
#'
#' @return save simulation results to environment
#' \item{b_AR}{evaluation of GLM estimates for correctly measured Y}
#' \item{b.obs_AR}{evaluation of GLM estimates for misclassified Y}
#' \item{b.EM_AR}{evaluation of EM estimates for misclassified Y}
#' @examples
#' SimImpl_AR()
#' @keywords simulation
#' @keywords implementation
#' @export
SimImpl_AR <- function(simN = 3, sampleN = 500, max.obs = 5, n.binom = 0, n.cat = 0, n.cont = 0,
  n.v = 0, n.lev = 3, n.grp = 3, saveEst = TRUE, sens = c(0.7,0.8,0.9), spec = c(0.9,0.8,0.7),
  bootN = 10, seed = 616) {
  # Make sure the dimension of sensitivity/specificity is consistent with n.grp
  # Stop simulation and print error message if not
  if(length(sens) != n.grp | length(spec) != n.grp) {
    stop("the dimension of sensitivity/specificity does not conform with number of group")
  }

  options(warn = -1)

  set.seed(seed)  # Set the root seed
  seed.pool <- runif(simN+1,2321,40245)  # Set the seed pool for data generation
  subID <- sample(1:sampleN, sampleN)  # generate subject ID
  grpIdx <- sample(1:n.grp, sampleN, replace = TRUE)  # generate group index for each subject
  n.vv <<- n.v
  n.cov <<- n.binom + n.cat + n.cont + n.v
  # Sensitivity = Prob(Y = 1 | Y.obs = 1)
  sens.grp <<- sens
  # Specificity = Prob(Y.obs = 0 | Y = 0)
  spec.grp <<- spec

  # Initialize betas'
  beta <- NULL
  beta.obs <- NULL
  beta.EM <- NULL

  # Perform simN times of data generation & evaluation
  for(i in 1:simN) {
    # Simulates correctly measured covariates
    x.cov <- sim_covariates_AR(
    # Number of subjects in sample
    n = sampleN,
    # Number of binomial covariates
    n.binom = n.binom,
    # Number of categorial covariates
    n.cat = n.cat,
    # Number of continuous covariates
    n.cont = n.cont,
    # Number of varying variable
    n.vv = n.v,
    # n.lev: number of levels for categorical covariate
    n.lev = n.lev,
    # group index for each subject
    grpIdx = grpIdx,
    # subject ID
    subID = subID, 
    # Set the seed for simulation
    k = seed.pool[i + 1]
    )    # End of simulation

    if(n.cov > 0 && n.vv > 0) {
       # index of time-dependent variables, can be input manually for real data
      vv.index <- rep(0,ncol(x.cov))
      vv.index[3:(2+n.vv)] <- 1
      # Parse the order of covariates
      x.cov <- parse_cov(x.cov = x.cov, vv.index = vv.index)
    }
    # Set sens/spec for groups with different test tool
    # sub.sens: sensitivity of each subject in the dataset
    # sub.spec: specificity of each subject in the dataset
    if (n.grp == 1) {
      sub.sens <- sens
      sub.spec <- spec
    } else {
      idx <- as.factor(x.cov$grpIdx)
      grp.mat <- model.matrix(~idx-1)    
      sub.sens <- grp.mat %*% sens.grp 
      sub.spec <- grp.mat %*% spec.grp
    }

    # Simulate imperfect test result
    data <- simulate_AR(
      # Number of subjects in sample
      n = sampleN,
      # Sensitivity = Prob(Y = 1 | Y.obs = 1)
      sens = sub.sens,
      # Specificity = Prob(Y.obs = 0 | Y = 0)
      spec = sub.spec,
      # matrix of correctly measured covariates
      x.cov = x.cov,
      # Number of varying variable
      n.vv = n.v, 
      # ID for each subject
      subID = subID,
      # Set the seed for simulation
      k = seed.pool[i]
    )  # End of simulation
    
    for (grp in unique(grpIdx)) {
      grpbeta <- data.frame(cbind(calc_group_AR(grpdata = data[data$grpIdx == grp, ], grp = grp),
        grp))
      beta <- rbind(beta, grpbeta[1, ])
      beta.obs <- rbind(beta.obs, grpbeta[2, ])
      beta.EM <- rbind(beta.EM, grpbeta[3, ])
    }

    # Bootstrap Confidence Interval of beta and prevalence for resampled data
    if (bootN != 0) {
      BootstrapCI_AR(bootN, data, i)
    }
  }

  if(saveEst == TRUE) {
    beta.EM <<- beta.EM
    beta.obs <<- beta.obs
    beta.true <<- beta
  }

  ###### Simulation Evaluation #######
  true <<- rep(1,ncol(beta))
  true[1] <<- 0

  # Save the simulation results
  b_AR <<- round(PerformanceEvaluation_AR(beta), 4)
  b.obs_AR <<- round(PerformanceEvaluation_AR(beta.obs), 4)
  b.EM_AR <<- round(PerformanceEvaluation_AR(beta.EM, sim = simN, cvg = TRUE), 4)
}

#' Predict STI probability
#'
#' Predict STI infection probability based on provided dataset and patient info
#'
#' @param data dataset for beta estimation
#' @param patInfo patient's diagonistic information
#' @param n.cov number of covariates
#' @param n.grp n.grp: number of testing technology
#' @param n.v Number of time varying variable
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param getBeta the option to save beta estimated from given dataset
#' @param CI the option to get 95\% confidence interval
#' @param bootN the number of BootStrapping used to construct CI
#'
#' @return specific patient's infection probability
#' \item{beta.est}{estimated beta from provided daaset}
#' \item{pat.CI.prev}{patient's 95\% CI of prevalence}
#' \item{pat.CI.prev.emp}{patient's 95\% empirical CI of prevalence}
#' @examples
#' data <- simExample_AR()
#' patInfo <- data[5, ]
#' predict_AR(data,patInfo,4,3,1,c(0.8,0.85,0.9),c(0.9,0.85,0.8))
#' beta.est
#' @keywords prediction
#' @keywords prevalence
#' @export
predict_AR <- function (data, patInfo, n.cov, n.grp, n.v, sens, spec, getBeta = TRUE, CI = FALSE,
  bootN = 10) {
  # Make sure the dimension of sensitivity/specificity is consistent with n.grp
  # Stop simulation and print error message if not
  if(length(sens) != n.grp | length(spec) != n.grp) {
    stop("the dimension of sensitivity/specificity does not conform with number of group")
  }

  # detect group index in the given data
  grpIdx <<- unique(data$grpIdx)
  n.cov <<- n.cov
  # detect number of subjects in the given data
  n <<- dim(data)[1]
  # get patient's group
  grp <- patInfo$grp
  # get patient's observed infection status
  y <- patInfo$yt.obs
  # get patient's observed x1
  y1 <- patInfo$yt1.obs
  # Sensitivity = Prob(Y = 1 | Y.obs = 1)
  sens.grp <<- sens
  # Specificity = Prob(Y.obs = 0 | Y = 0)
  spec.grp <<- spec
  n.vv <<- n.v
  if(n.cov != 0) {
    # Extract patient's covariates
    x.cov <- patInfo[ , (ncol(patInfo)-n.cov+1):ncol(patInfo)]
    # Convert x.cov into matrix
    Xcov <- as.numeric(ConvertCov(x.cov))
  }

  # get EM estimation of beta for each group
  # Initialize estimation
  beta.EM <- NULL
  for (grp in unique(grpIdx)) {
    grpdata <- data[data$grpIdx == grp, ]
    if(n.cov == 0) {
      x.cov <- NULL
    } else {
      x.cov <- grpdata[ , (ncol(grpdata)-n.cov+1):ncol(grpdata)]
    }
    b <- CalculateBetaEM_AR(data = grpdata, x.cov = x.cov, sens = sens.grp[grp],
                spec = spec.grp[grp])
    beta.EM <- rbind(beta.EM, b)
  }
  # average beta from different groups
  beta <- apply(beta.EM, 2, mean)
  if(getBeta == TRUE) beta.est <<- beta
  # Calculate Empirical and Statistical 95% CI
  if(CI  == TRUE) {
    BootstrapCI_AR(bootN,data)
    # adjusted probability for this subject 
    if(y == 0) {
      if(y1 == 0) {
        case <- 1
      } else {
        case <- 2
      }
    } else {
      if(y1 == 0) {
        case <- 3
      } else {
        case <- 4
      }
    }
    if (n.cov != 0) {
      adj <- Xcov %*% beta[-c(1,2)]
      pat.CI.prev <<- 1 / (1 + exp(log(1 / CI.prev[ , case] - 1) - adj))
      pat.CI.prev.emp <<- 1 / (1 + exp(log(1 / CI.prev.emp[ , case] - 1) - adj))
    } else {
      pat.CI.prev <<- CI.prev[ , case]
      pat.CI.prev.emp <<- CI.prev.emp[ , case]
    }
  }
  # Calculate prevalence
  if(n.cov == 0) {
    prev.y1 <- 1 / (1 + exp(-beta[1] - beta[2] * y1))
  } else {
    prev.y1 <- 1 / (1 + exp(-beta[1] - beta[2] * y1 - Xcov %*% beta[-c(1,2)]))
  }
  prev <- prev.y1 * sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) / 
        (prev.y1*sens.grp[grp]^y * (1-sens.grp[grp])^(1-y) + 
        (1-prev.y1) * (1-spec.grp[grp])^y * spec.grp[grp]^(1-y))
  if(prev > 1) prev <- 1
  if(prev < 0) prev <- 0
  prev
}

#' Simulate an example of AR(1) dataset
#'
#' Simulation an example for the users to prepare their data structure
#'
#' @param sampleN number of tuples in each group
#' @param max.obs maximum number of obervation for each subject
#' @param n.binom number of binomial covariates
#' @param n.cat number of categorial covariates
#' @param n.cont number of continuous covariates
#' @param n.v number of varying variable
#' @param n.grp number of testing technology
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param seed the root seed of simulation
#'
#' @return an example of dataset
#' @examples
#' data <- simExample_AR()
#' head(data)
#' @keywords Example 
#' @keywords simulate 
#' @export
simExample_AR <- function(sampleN = 500, max.obs = 5, n.binom = 1, n.cat = 1, n.cont = 1,
  n.v = 1, n.grp = 3, sens = c(0.7,0.8,0.9), spec = c(0.9,0.8,0.7), seed = 616) {
  # Set the root seed
  set.seed(seed) 
  subID <- sample(1:sampleN, sampleN)  # generate subject ID
  grpIdx <- sample(1:n.grp, sampleN, replace = TRUE)  # generate group index for each subject
  n.cov <- n.binom + n.cat + n.cont + n.v  # number of covariates
  # Sensitivity = Prob(Y = 1 | Y.obs = 1)
  sens.grp <- sens
  # Specificity = Prob(Y.obs = 0 | Y = 0)
  spec.grp <- spec
  # Simulates correctly measured covariates
  x.cov <- sim_covariates_AR(
    # Number of subjects in sample
    n = sampleN,
    # Number of binomial covariates
    n.binom = n.binom,
    # Number of categorial covariates
    n.cat = n.cat,
    # Number of continuous covariates
    n.cont = n.cont,
    # Number of varying variable
    n.vv = n.v,
    # group index for each subject
    grpIdx = grpIdx,
    # subject ID
    subID = subID, 
    # Set the seed for simulation
    k = seed
   )   # End of simulation

  # Set sens/spec for groups with different test tool
  # sub.sens: sensitivity of each subject in the dataset
  # sub.spec: specificity of each subject in the dataset
  if (n.grp == 1) {
    sub.sens <- sens
    sub.spec <- spec
  } else {
    idx <- as.factor(x.cov$grpIdx)
    grp.mat <- model.matrix(~idx-1)    
    sub.sens <- grp.mat %*% sens.grp 
    sub.spec <- grp.mat %*% spec.grp
  }

  # Simulate imperfect test result
  data <- simulate_AR(
    # Number of subjects in sample
    n = sampleN,
    # Sensitivity = Prob(Y = 1 | Y.obs = 1)
    sens = sub.sens,
    # Specificity = Prob(Y.obs = 0 | Y = 0)
    spec = sub.spec,
    # matrix of correctly measured covariates
    x.cov = x.cov,
    # Number of varying variable
    n.vv = n.v, 
    # ID for each subject
    subID = subID,
    # Set the seed for simulation
    k = seed + 1
  )  # End of simulation
  # Subtract informations provided by user
  if(n.cov == 0) {
    data.user <- data[, c(1,2,4,6,8,9)]
  } else {
    data.user <- cbind(data[, c(1,2,4,6,8)], data[ , (ncol(data)-n.cov-n.v):ncol(data)])
  }
}
