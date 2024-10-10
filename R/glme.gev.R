# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#----------------------------------------------------------------
# L-me under a preference (prior) function for stationary GEV ---
# ---------------------------------------------------------------


MS_pk = function(para=para, p=6,q=9){

  z=list()

  Bef <- function(x) { ((0.5+x)^(p-1)) * ((0.5-x)^(q-1)) }
  Be  <- integrate(Bef, lower=-0.5, upper=0.5)[1]$value

  if( abs(para[3]) < 0.5 ){
    z$pk.ms <- ((0.5+para[3])^(p-1))*((0.5-para[3])^(q-1))/ Be
  }else if ( abs(para[3]) >= 0.5 ) {
    z$pk.ms = 1e-20
  }

  # z$fmax= 3.11898237
  # z$maxit = -0.11538

  return(z)
}

# 주어진 매개변수에 대해 GEV 모델의 L-moments 계산
glme.sol <- function(para, slmgev=slmgev, eps=eps) {

  zz=rep(10000,3)

  if(para[2] <= 0) return(zz)
  if(abs(para[3]) > 1.0) return(zz)
  if(abs(para[3]) < 1e-5) para[3]= -(1e-5)

  emom=list()
  emom=lmomgev(vec2par(para,'gev'))
  if( is.na(emom$lambda[1]) ) return(zz)
  if( is.na(emom$lambda[2]) ) return(zz)
  if( is.na(emom$ratios[3]) ) return(zz)

  #  slmgev=lmoms(xdat)

  zz[1] = emom$lambdas[1]-slmgev$lambdas[1]
  zz[2] = emom$lambdas[2]-slmgev$lambdas[2]
  zz[3] = max( abs( emom$ratios[3]-slmgev$ratios[3] ) - eps, 0 )

  return(zz)
}


new_pf_norm =function(para=NULL, mu=NULL, std=NULL){

  Brone = 1 + dnorm(para[3], mean= mu, sd=std) #*2

  return(Brone)
}

pk.beta.stnary = function(para=NULL, lme.center=NULL,  p=NULL){

  pk.one = 1e-10
  ulim= 0.3
  aa= max(-1.0, lme.center[3]-ulim)
  bb= min(0.3, lme.center[3]+ulim)
  al=min(aa,bb)
  bl=max(aa,bb)

  if(lme.center[3] <= 0) {

    c1=10
    c2=5

    qlim= min( 0.0+abs(lme.center[3])*c1, c2 )
  }else{ qlim =0.0 }

  #  qlim=10

  p=p; q=p+qlim

  Bef <- function(x) { ((-al+x)^(p-1)) * ((bl-x)^(q-1)) }
  Be  <- integrate(Bef, lower=al, upper=bl)[1]$value

  if(lme.center[3] <= 0.1){
    if( (para[3] > al) & (para[3] < bl) ) {
      pk.one <- ((-al+para[3])^(p-1))*((bl-para[3])^(q-1))/ Be
    }
  }
  return(pk.one)
}


#' Calculate the likelihood for Generalized L-moments estimation of GEV distribution
#'
#' @description
#' This function calculates the likelihood (or more precisely, a penalized negative log-likelihood)
#' for the Generalized L-moments estimation of the Generalized Extreme Value (GEV) distribution.
#'
#' @param par A vector of GEV parameters (location, scale, shape).
#' @param xdat A numeric vector of data.
#' @param slmgev Sample L-moments of the data.
#' @param covinv Inverse of the covariance matrix of the sample L-moments.
#' @param lcovdet Log determinant of the covariance matrix.
#' @param mu Mean for the normal penalization (used when pen='norm').
#' @param std Standard deviation for the normal penalization (used when pen='norm').
#' @param lme L-moment estimates of the parameters.
#' @param pen Penalization method ('norm' or 'beta').
#'
#' @details
#' The function performs the following steps:
#' 1. Checks if the parameters are within valid ranges.
#' 2. Calculates the expected L-moments based on the current parameters.
#' 3. Computes the difference between expected and sample L-moments.
#' 4. Calculates the generalized L-moments distance.
#' 5. Applies a penalization term based on the specified method ('norm' or 'beta').
#' 6. Returns the sum of the L-moments distance and the penalization term.
#'
#' @return A numeric value representing the penalized negative log-likelihood.
#' A lower value indicates a better fit.
#'
#' @author [Jeong-Soo Park]
#' @export
glme.like = function(par=par, xdat=xdat, slmgev=slmgev, covinv=covinv,
                     lcovdet=lcovdet, mu=mu, std=std, lme=lme, pen=pen){

  if( par[2] <= 0) return(10^8)
  if( abs(par[3]) > 1) return(10^8)

  nsample=length(xdat)

  if( abs(par[3]) < 1e-5) par[3]= -(1e-5)

  emom= lmomgev( vec2par(par,type='gev') )

  if( is.na(emom$lambda[1]) ) return(10^8)
  if( is.na(emom$lambda[2]) ) return(10^8)
  if( is.na(emom$ratios[3]) ) return(10^8)

  # slmgev=lmoms(xdat)
  # cov=lmoms.cov(xdat, nmom=3)

  zvec=rep(NA,3)
  zvec= emom$lambdas[1:3] - slmgev$lambdas[1:3]

  if( any(is.na(zvec)) ) return(10^8)

  z= t(zvec) %*% covinv %*% zvec

  prob.norm =  z/2   + (3/2)*log( (2*pi) ) + lcovdet

  if(pen=='norm'){

    pk_beta = -log( new_pf_norm(para=par, mu= -0.55, std= 0.15) )

  }else if(pen=='beta'){

    pk_beta = -log( pk.beta.stnary(para= par, lme.center=lme, p=6) )
  }

  zz= prob.norm  + pk_beta

  return(zz)
}

#' Initialize parameters for Generalized L-moments estimation of GEV distribution
#'
#' @description
#' This function initializes parameters for the Generalized L-moments estimation
#' of the Generalized Extreme Value (GEV) distribution.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of initial parameter sets to generate.
#'
#' @details
#' The function generates `ntry` sets of initial parameters for the GEV distribution.
#' It uses L-moment estimates as a starting point and then generates additional
#' sets of parameters using random perturbations. This approach increases the
#' likelihood of finding a global optimum in the subsequent optimization process.
#'
#' The function ensures that:
#' 1. The scale parameter (sigma) is always positive.
#' 2. The shape parameter (xi) is constrained between -0.7 and 0.7.
#'
#' @return A matrix with `ntry` rows and 3 columns, where each row represents
#' a set of initial parameters (location, scale, shape) for the GEV distribution.
#'
#' @author [Jeong-Soo Park]
#' @export
init.glme <-function(xdat, ntry=ntry){

  init <-matrix(0, nrow=ntry, ncol=3)

  lmom_init = lmoms(xdat)
  lmom_est <- pargev(lmom_init)

  init[1,]=lmom_est$para

  sd1= max(abs(init[1,1])*0.1, 2)
  sd2= max(init[1,2]*2, 2)

  init[2:ntry,1] <- init[1,1]+rnorm(n=(ntry-1),mean=0, sd = sd1)
  init[2:ntry,2] <- runif(n=(ntry-1), min= init[1,2]*0.5, max= sd2)
  init[2:ntry,3] <- runif(n=(ntry-1), min= -.7, max=.5)

  for (i in 1:ntry){
    if(init[i,2] <= 0) init[i,2] = 1.0
    if(abs(init[i,3]) >= 0.7 ) init[i,3]= 0.69*sign(init[i,3])
  }

  return(init)
}


#' Generalized L-moments estimation for generalized extreme value distribution
#'
#' @description
#' This function estimates the Generalized L-moments of Generalized Extreme Value distribution using an updated algorithm.
#'
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for parameter estimation. Higher values increase the chance of finding a more accurate estimate by trying different initial conditions.
#' @param pen Penalization method used in estimation. Options are 'beta' and 'norm'.
#'
#' @details
#' The function uses an optimization approach to estimate the parameters of the Generalized Extreme Value distribution.
#' It implements a penalized likelihood method, where the penalization can be either 'beta' or 'norm'.
#' The function handles potential issues with covariance matrix calculation by using a bootstrap approach when necessary.
#'
#' @return The glme.gev function returns a list containing the following elements:
#' \itemize{
#'  \item glme - The estimated parameters of the Generalized Extreme Value distribution.
#'  \item lme - The L-moment estimates of the parameters.
#'  \item covinv - The inverse of the covariance matrix of the L-moments.
#'  \item lcovdet - The log determinant of the covariance matrix.
#'  \item nllh.pref - The negative log-likelihood of the preferred solution.
#'  \item pen - The penalization method used ('beta' or 'norm').
#' }
#'
#' @author [Jeong-Soo Park]
#' @export
glme.gev= function(xdat, ntry=10, pen='beta'){

  # pen='beta' or 'norm'

  z=list()
  k =list()

  # initial setting ------
  nsample=length(xdat)
  sinit=matrix(0, nrow=ntry, ncol=3)

  sinit <- init.glme(xdat, ntry=ntry)

  lmom_init = lmoms(xdat)
  lmom_est <- pargev(lmom_init)

  lme = lmom_est$para
  z$lme =  lmom_est$para

  precis=rep(NA, ntry)
  pk.ms=rep(NA, ntry)
  pk.ms.lme=rep(NA, ntry)
  isol=0
  sol=list()
  mindist=1000
  dist=rep(1000, ntry)

  covinv= matrix(NA, 3, 3)

  slmgev=lmoms(xdat)
  cov=lmoms.cov(xdat, nmom=3)

  covinv=solve(cov)
  detc = det(cov)

  #--------------------------------------------------
  if(detc <= 0){
    # cat("det <0, Bootstrap to calculate cov","\n")

    BB=200          # we need Bootstrap to calculate cov ---
    sam.lmom= matrix(NA,BB,3)

    for (ib in 1:BB){
      sam.lmom[ib,1:3]=lmoms(sample(xdat,size=nsample,replace=T), nmom=3)$lambdas
    }
    cov=cov(sam.lmom)
    covinv=solve(cov)
    detc=det(cov)
  }

  lcovdet=log(detc)
  z$covinv =covinv
  z$lcovdet =lcovdet

  #-------------------------------------------------------
  # estimating paras using nleqslv or optim
  tryCatch(
    for(i in 1:ntry){

      value=list()

      value <- try(
        optim(par=as.vector(sinit[i,1:3]), fn=glme.like,
              xdat=xdat, slmgev=slmgev, covinv=covinv,
              lcovdet=lcovdet, mu=mu, std=std, lme=lme, pen=pen)
      )

      if(is(value)[1]=="try-error"){
        k[[i]]$fvec <- 10^6
      }else{
        k[[i]] <- value
        k[[i]]$root = value$par
        k[[i]]$fvec = value$value
      }

      if( value$convergence != 0) {precis[i]=10^6
      }else{
        isol=isol+1
        precis[i] = k[[i]]$fvec
      }

    } #for
  ) #tryCatch

  if(isol==0) {
    cat("-- No solution was found in nleqslv or optim --","\n")
    z$glme = z$lme
    return
  }

  selc_num = which.min( precis )    #precis=k[[i]]$fvec

  x  <- k[[selc_num]]

  z$nllh.pref = k[[selc_num]]$fvec
  z$glme = x$root
  z$pen = pen

  return(z)
}




pargev.kfix= function (lmom, kfix= 0.1, checklmom = TRUE, ...)
{

  # modified from 'pargev' function in lmomco package

  para <- rep(NA, 3)
  names(para) <- c("xi", "alpha", "kappa")
  SMALL <- 1e-05
  EPS <- 1e-06
  MAXIT <- 20
  EU <- 0.57721566
  DL2 <- 0.69314718
  DL3 <- 1.0986123
  A0 <- 0.2837753
  A1 <- -1.21096399
  A2 <- -2.50728214
  A3 <- -1.13455566
  A4 <- -0.07138022
  B1 <- 2.06189696
  B2 <- 1.31912239
  B3 <- 0.25077104
  C1 <- 1.59921491
  C2 <- -0.48832213
  C3 <- 0.01573152
  D1 <- -0.64363929
  D2 <- 0.08985247
  if (length(lmom$L1) == 0) {
    lmom <- lmorph(lmom)
  }
  if (checklmom & !are.lmom.valid(lmom)) {
    warning("L-moments are invalid")
    return()
  }
  T3 <- lmom$TAU3
  if (T3 > 0) {
    Z <- 1 - T3
    G <- (-1 + Z * (C1 + Z * (C2 + Z * C3)))/(1 + Z * (D1 +
                                                         Z * D2))
    if (abs(G) < SMALL) {
      para[3] <- 0
      para[2] <- lmom$L2/DL2
      para[1] <- lmom$L1 - EU * para[2]
      return(list(type = "gev", para = para))
    }
  }
  else {
    G <- (A0 + T3 * (A1 + T3 * (A2 + T3 * (A3 + T3 * A4))))/(1 +
                                                               T3 * (B1 + T3 * (B2 + T3 * B3)))
    if (T3 >= -0.8) {
    }
    else {
      if (T3 <= -0.97)
        G <- 1 - log(1 + T3)/DL2
      T0 <- (T3 + 3) * 0.5
      CONVERGE <- FALSE
      for (it in seq(1, MAXIT)) {
        X2 <- 2^-G
        X3 <- 3^-G
        XX2 <- 1 - X2
        XX3 <- 1 - X3
        T <- XX3/XX2
        DERIV <- (XX2 * X3 * DL3 - XX3 * X2 * DL2)/(XX2 *
                                                      XX2)
        GOLD <- G
        G <- G - (T - T0)/DERIV
        if (abs(G - GOLD) <= EPS * G)
          CONVERGE <- TRUE
      }
      if (CONVERGE == FALSE) {
        warning("Noconvergence---results might be unreliable")
      }
    }
  }

  #  para[3] <- G
  para[3]=kfix
  G = kfix

  GAM <- exp(lgamma(1 + G))
  para[2] <- lmom$L2 * G/(GAM * (1 - 2^(-G)))
  para[1] <- lmom$L1 - para[2] * (1 - GAM)/G
  return(list(type = "gev", para = para, source = "pargev"))
}



