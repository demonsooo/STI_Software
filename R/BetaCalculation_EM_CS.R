#-----------------------------------------------------------------------
#     VERSION HISTORY OF BetaCalculation_EM_CS.R
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

# Version: v2.0.0
# Date: 10/04/2016
# Author: Yijia Liu
# Comment: Finalized for package

# Version: v2.0.1
# Date: 10/09/2016
# Author: Yijia Liu
# Comment: Improved generalization of calculations

# Version: v2.0.2
# Date: 10/23/2016
# Author: Yijia Liu
# Comment: zero covariate case added, small calculation error fixed

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------

#' probability expectation matrix of CS dataset
#' 
#' an internal function to compute the probability expectation matrix of four cases used in E step
#' Define Four Cases:
#' C1: y = 1; x1 = 0
#' C2: y = 1; x1 = 1
#' C3: y = 0; x1 = 0
#' C4: y = 0; x1 = 1
#'
#' @param b the beta used to calculate negative log likelihood
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param sens.x Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
#' @param spec.x Specificity.x = Prob(X1.obs = 0 | X1 = 0)
#'
#' @return he expected probability of four cases for given beta
#' @keywords expectation
#' @keywords probability
ExpY_CS <- function(b, data, x.cov, sens, spec, sens.x, spec.x){
  # Eliminated because of unpublished paper
}

#' Negative Log-Likelihood
#'
#' An internal function to compute the expected value of negative log likelihood function
#' The E step in EM algorithm
#'
#' @param b the beta used to calculate negative log likelihood
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param EY the probability expectation matrix calculated in the ExpY function
#' @return expected value of negative log likelihood
#' @keywords LogLikelihood
NegativeLogLikelihood_CS <- function(b, data, EY, x.cov) {
  # Eliminated because of unpublished paper
}

#' Gradient of CS model
#'
#' An internal function to omputes the gradient of objective function(nLL) respect to beta
#'
#' @param b the beta used to calculate negative log likelihood
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param EY the probability expectation matrix calculated in the ExpY function
#' @return expected value of negative log likelihood
#' @keywords LogLikelihood
Gradient_CS <- function(b, data, EY, x.cov) {
  # Eliminated because of unpublished paper
}

#' convert covariates
#'
#' an internal function to convert x.cov to covariate matrix used in calculation
#'
#' @param x.cov matrix of correctly measured covariates
#'
#' @return processed covariate matrix used in calculation
ConvertCov <- function(x.cov) {
  # get type each covariate
  x.class <- sapply(x.cov, class)
  # detect the index of categorical variable
  is.cat <- which(x.class == 'factor') 
  # extract covariates that are not categorical variable
  not.cat <- which(x.class != 'factor')
  if (length(is.cat) == 0){ 
    X <- data.matrix(x.cov[ ,not.cat])
  } else {
    # Create dummy matrix for categorical variables
    x.dummy <- create_dummy_mat(x.cov = x.cov, cat.index = is.cat)
    if(class(x.dummy) == "numeric") {
      X <- t(data.matrix(c(x.cov[ ,not.cat], x.dummy)))
    } else {
      X <- data.matrix(cbind(x.cov[ ,not.cat], x.dummy))
    }
  }
  X
}

#' EM estimation
#'
#' Excecute EM algorithm to estimate beta
#'
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param sens.x Sensitivity.x = Prob(X1 = 1 | X1.obs = 1)
#' @param spec.x Specificity.x = Prob(X1.obs = 0 | X1 = 0)
#' @param MaxIt maximum number of iteration
#' @param crt criteria of convergence
#' @return optimal beta estimated by EM algorithm
#' @keywords EM
#' @keywords regression
CalculateBetaEM_CS <- function(data, x.cov, sens, spec, sens.x, spec.x,
  # Eliminated because of unpublished paper
}

#' separate group calculation
#' 
#' an internal function to extract data from different group and estimate beta's
#'
#' @param grpdata data from a specific group
#' @param grp group index
#' @param maxit maximum number of iteration
#' @param crt criteria of convergence
#'
#' @return optimal beta estimated by EM algorithm and GLM
calc_grp_CS <- function(grpdata, grp, maxit = 2000, crt = 10^(-8)) {
  y <- factor(c(grpdata$y))
  y.obs <- factor(c(grpdata$y.obs))
  x1 <- c(grpdata$x1)
  x1.obs <- c(grpdata$x1.obs)
  if(n.cov <= 0) {
    x.cov <- NULL
    # GLM estimates for correctly measured Y
    b <- c(glm(y~x1,family=binomial(link="logit"))$coefficients)
    # GLM estimates for misclassified Y
    b.obs <- c(glm(y.obs~x1.obs,family=binomial(link="logit"),
      x=T, y=T)$coefficients)
  } else {
    x.cov <- grpdata[ , (ncol(grpdata)-n.cov+1):ncol(grpdata)]
    # get type each covariate
    x.class <- sapply(x.cov, class)
    # detect the index of categorical variable
    is.cat <- which(x.class == 'factor')
    # extract covariates that are not categorical variable
    not.cat <- which(x.class != 'factor')
    ###### GLM #######
    rank <- if(length(is.cat)!=0) {factor(x.cov[ ,is.cat])
              } else {rep(0,length(y))}   ## Deal with the case with no categorical variables
    cov <- if(length(not.cat)!=0) {data.matrix(x.cov[ ,not.cat])
      } else {rep(0,length(y))} ## Deal with the case with only categorical variables
    # GLM estimates for correctly measured Y
    b <- c(glm(y~x1+cov+rank,family=binomial(link="logit"))$coefficients)
    # GLM estimates for misclassified Y
    b.obs <- c(glm(y.obs~x1.obs+cov+rank,family=binomial(link="logit"),
      x=T, y=T)$coefficients)
  }  
  ####### LL-EM #######
  b.EM <- CalculateBetaEM_CS(data = grpdata, x.cov = x.cov, sens = sens.grp[grp],
    spec = spec.grp[grp], sens.x = sens.x.grp[grp], spec.x = spec.x.grp[grp],
    MaxIt = maxit, crt = crt)
  # Retuen the results
  rbind(b,b.obs,b.EM)
}
