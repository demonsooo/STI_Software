#-----------------------------------------------------------------------
#     VERSION HISTORY OF BetaCalculation_EM_AR.R
#-----------------------------------------------------------------------
​
# Version: v1.0.0
# Date: 04/07/2016
# Author: Yijia Liu
# Comment: Initial version. 

# Version: v1.0.1
# Date: 04/13/2016
# Author: Yijia Liu
# Comment: Calculations for Y0 is corrected and some temporary variables
#          are added for acceleration

# Version: v2.0.1
# Date: 05/25/2016
# Author: Yijia Liu
# Comment: arbitrary number of fixed variables and one time-dependent variable

# Version: v2.0.1
# Date: 07/14/2016
# Author: Yijia Liu
# Comment: Complete EM algorithm of AR(1) model with arbitrary number of  
#          binary, categorical, and continuous covariates, change or not 
#          change within time

# Version: v3.0.0
# Date: 08/23/2016
# Author: Yijia Liu
# Comment: Seperate calculation for groups using different test tool

#-----------------------------------------------------------------------
#     BEGINNING OF THE CODE
#-----------------------------------------------------------------------
#' probability expectation matrix
#'
#' An internal function to compute the probability expectation matrix of four cases used in E step
#' Define Four Cases:
#' C1: y = 1; y(t-1) = 0
#' C2: y = 1; y(t-1) = 1
#' C3: y = 0; y(t-1) = 0
#' C4: y = 0; y(t-1) = 1
#'
#' @param b the beta used to calculate negative log likelihood
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#'
#' @return The probability expectation matrix of given beta
#' @keywords expectation
#' @keywords probability
ExpY_AR <- function(b, data, x.cov, sens = 0.7, spec = 0.9){
  # Eliminated because of unpublished paper
}

#' Negative Log-Likelihood of AR(1) model
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
NegativeLogLikelihood_AR <- function(b, data, x.cov, EY) {
  # Eliminated because of unpublished paper
}

#' Gradient of AR(1) model
#'
#' An internal function to omputes the gradient of objective function(nLL) respect to beta
#'
#' @param b the beta used to calculate negative log likelihood
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param EY the probability expectation matrix calculated in the ExpY function
#' @return expected value of negative log likelihood
#' @keywords LogLikelihood
Gradient_AR <- function(b, data, x.cov, EY) {
  # Eliminated because of unpublished paper
}

#' EM estimation
#'
#' Excecute EM algorithm to estimate beta
#'
#' @param data the dataset to be evaluated
#' @param x.cov matrix of correctly measured covariates
#' @param sens Sensitivity = Prob(Y = 1 | Y.obs = 1)
#' @param spec Specificity = Prob(Y.obs = 0 | Y = 0)
#' @param MaxIt maximum number of iteration
#' @param crt criteria of convergence
#' @return optimal beta estimated by EM algorithm
#' @keywords EM
#' @keywords regression
CalculateBetaEM_AR <- function(data, x.cov, sens = 0.7,spec = 0.9, MaxIt = 2000, crt=10^(-8)){
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
calc_group_AR <- function(grpdata, grp, maxit = 2000, crt = 10^(-8)) {
  yt <- factor(grpdata$yt)
  yt.obs <- factor(grpdata$yt.obs)
  yt1 <- grpdata$yt1
  yt1.obs <- grpdata$yt1.obs
  
  if(n.cov <= 0) {
    x.cov <- NULL
    # GLM estimates for correctly measured Y
    b <- c(glm(yt~yt1,family=binomial(link="logit"))$coefficients)
    # GLM estimates for misclassified Y
    b.obs <- c(glm(yt.obs~yt1.obs,family=binomial(link="logit"),
      x=T, y=T)$coefficients)
  } else {
    x.cov <- grpdata[ , (ncol(grpdata)-n.cov+1):ncol(grpdata)]
    # get type each covariate
    x.class <- sapply(x.cov, class)
    # detect the index of categorical variable
    is.cat <- which(x.class == 'factor')
    # extract covariates that are not categorical variable
    not.cat <- which(x.class != 'factor')
    rank <- if(length(is.cat)!=0) {factor(x.cov[ ,is.cat])
              } else {rep(0,length(yt))}   ## Deal with the case with no categorical variables
    cov <- if(length(not.cat)!=0) {data.matrix(x.cov[ ,not.cat])
              } else {rep(0,length(yt))} ## Deal with the case with only categorical variables
    ###### GLM #######
    # GLM estimates for correctly measured Y
    b <- c(glm(yt~yt1+cov+rank,family=binomial(link="logit"))$coefficients)
    # GLM estimates for misclassified Y
    b.obs <- glm(yt.obs~yt1.obs+cov+rank,family=binomial(link="logit"),x=T, y=T)$coefficients
  }
  ####### LL-EM #######
  b.EM <- CalculateBetaEM_AR(data = grpdata, x.cov = x.cov, sens = sens.grp[grp],
  						  spec = spec.grp[grp], MaxIt = maxit, crt = crt)
  # Retuen the results
  rbind(b,b.obs,b.EM)
}
