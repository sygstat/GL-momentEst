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

glme.like = function(par=par, xdat=xdat, slmgev=slmgev, covinv=covinv,
                     lcovdet=lcovdet, lamb=lamb,
                     pref=pref,p=p,q=q){

  if( par[2] <= 0) return(10^8)
  if( abs(par[3]) > 0.5) return(10^8)

  nsample=length(xdat)

  if( abs(par[3]) < 1e-5) par[3]= -(1e-5)

  emom= lmomgev( vec2par(par,type='gev') )

  if( is.na(emom$lambda[1]) ) return(10^8)
  if( is.na(emom$lambda[2]) ) return(10^8)
  if( is.na(emom$ratios[3]) ) return(10^8)

  # slmgev=lmoms(xdat)
  # cov=lmoms.cov(xdat, nmom=3)

  zvec=rep(NA,3)
  zvec= emom$lambdas[1:3]-slmgev$lambdas[1:3]
  #  zvec[3] =emom$ratios[3]-slmgev$ratios[3]

  if( any(is.na(zvec)) ) return(10^8)

  z= t(zvec) %*% covinv %*% zvec

  prob.norm =  z/2   + (3/2)*log( (2*pi) ) + lcovdet

  pk_beta = - log( MS_pk(par, p=p,q=q)$pk.ms )

  zz= prob.norm  + pk_beta* lamb

  return(zz)
}

init.glme <-function(xdat, ntry=ntry){

  init <-matrix(0, nrow=ntry, ncol=3)

  lmom_init = lmoms(xdat)
  lmom_est <- pargev(lmom_init)

  init[1,]=lmom_est$para

  sd1= max(abs(init[1,1])*0.1, 2)
  sd2= max(init[1,2]*2, 2)

  init[2:ntry,1] <- init[1,1]+rnorm(n=(ntry-1),mean=0, sd = sd1)
  init[2:ntry,2] <- runif(n=(ntry-1), min= init[1,2]*0.5, max= sd2)
  init[2:ntry,3] <- runif(n=(ntry-1), min= -0.5, max=0.5)

  for (i in 1:ntry){
    if(init[i,2] <= 0) init[i,2] = 1.0
    if(abs(init[i,3]) >= 0.5 ) init[i,3]= 0.48*sign(init[i,3])
  }

  return(init)
}


# GEV 모델의 매개변수를 추정하는 함수를 정의합니다.
#' Generalized L-moments estimation for generalized extreme value distribution
#'
#' @description
#' This function estimates the Generalized L-moments of Generalized Extreme Value distribution.
#' @param xdat A numeric vector of data to be fitted.
#' @param ntry Number of attempts for parameter estimation. Higher values increase the chance of finding a more accurate estimate by trying different initial conditions.
#' @param alg Choice of algorithm for estimation. Options include 'pmax' and 'like', each with distinct estimation methodologies.
#' @param pref Preference function used in estimation. Options 'mso' and 'msa' influence the prioritization in the estimation process.
#' @param ftol Tolerance threshold for the estimation function's precision. Determines the stopping point of the algorithm based on the change in function value between iterations.
#' @param eps Small positive value used in the algorithm to avoid numerical issues like division by zero or to initiate steps in iterative processes.
#'
#' @details
#' The equations for the L-moments for LME of the GEVD are
#' \deqn{ \underline{\bf \lambda} - \underline{\bf l} = \underline{\bf 0},}
#' where \eqn{ \underline{\bf \lambda} =(\lambda_1,\; \lambda_2,\; \lambda_3)^t } and \eqn{\underline{\bf l} =(l_1,\; l_2,\; l_3)^t}.
#' Next, we define the generalized L-moments distance (GLD) as;
#' \deqn{(\underline{\bf \lambda} -\underline{\bf l})^t V^{-1} (\underline{\bf \lambda} -\underline{\bf l}),}
#' where \eqn{V} is the variance-covariance matrix of the sample L-moments up to the third order.
#' The LME value is obtained by minimizing the following function:
#' \deqn{f(GLD(\mu,\sigma,\xi)) = \frac{1}{(2 \pi)^{3/2} |V| }exp \{-{1\over2}  (\underline{\bf \lambda} -\underline{\bf l})^t V^{-1} (\underline{\bf \lambda} -\underline{\bf l}) \}.}
#' This function is a PDF of the three-variate normal random vector.
#' As sample L-moments typically converge to a multivariate normal distribution as \eqn{n \rightarrow \infty} (Hosking 1990),
#' this function can be treated as an approximation to the likelihood function of \eqn{\underline{\bf \theta}=(\mu,\sigma,\xi)} given three sample L-moments, \eqn{\tilde L (\underline{\mathbf \theta} | \underline{\mathbf l})}.
#' Then, the GLME value of \eqn{\underline{\bf \theta}=(\mu,\sigma,\xi)} is obtained by minimizing the following function with respect to \eqn{\underline{\bf \theta}}:
#' \deqn{-ln ( f(GLD(\underline{\bf \theta})))\; - \alpha_n \; ln ( p(\xi)),}
#' where \eqn{\alpha_n} denotes a weight for the prior \eqn{p(\xi)} compared to the approximated likelihood of \eqn{GLD(\underline{\bf \theta})}.
#'
#' @return The glme.gev function returns a list containing the following elements:
#' \itemize{
#'  \item nsol - Number of solutions found. Indicates how many viable parameter sets were identified during the estimation process.
#'  \item glme.conF - The estimated parameters of the Generalized Extreme Value distribution using the function's primary method. Represents the core output of the function.
#'  \item glme.conT - A variant of the estimated parameters, adjusted for certain conditions (e.g., boundary constraints on parameters).
#'  \item glme.bc - Bias-corrected version of the estimated parameters, accounting for potential biases in the estimation process.
#'  \item lme.conF and lme.conT - Similar to glme.conF and glme.conT, but derived from a different aspect of the algorithm or a different estimation perspective.
#'  \item algorithm - The algorithm used for the estimation, reflecting the alg input parameter.
#'  \item pref.func - The preference function used, corresponding to the pref input parameter.
#'  \item nsample - The sample size of the data set used in the estimation process.
#'  \item estim.precis - The precision of the estimated parameters, indicating the accuracy of the estimation.
#'  \item pmax - Additional information related to the 'pmax' algorithm, if used.
#' }
#' @author Jeong-Soo Park
#' @export
glme.gev= function(xdat, ntry=15, alg='pmax', pref='msa',
                   ftol=1e-5, eps=0.03 ){

  # Two algorithms (alg) are available: 'pmax', 'like'
  # Preference functions (pref) are 'mso' and 'msa'

  z=list()
  k =list()

  if(alg=='like') lamb=eps
  if(pref=='mso'){
    p=6; q=9
  }else if(pref=='msa'){
    p=3; q=4.5
  }

  # 데이터와 초기 매개변수 설정
  nsample=length(xdat)
  init=matrix(0, nrow=ntry, ncol=3)

  init <- init.glme(xdat, ntry=ntry)

  lmom_init = lmoms(xdat)
  lmom_est <- pargev(lmom_init)

  precis=rep(NA, ntry)
  pk.ms=rep(NA, ntry)
  pk.ms.lme=rep(NA, ntry)
  isol=1
  sol=list()
  mindist=1000
  dist=rep(1000, ntry)

  covinv= matrix(NA, 3,3)
  slmgev=lmoms(xdat)
  cov=lmoms.cov(xdat, nmom=3)
  covinv=solve(cov)
  det=det(cov)
  if(det <= 0){
    # 비정상성인 공분산의 경우 처리
    covinv[,]=0
    for (i in 1:3){
      covinv[i,i]=1
    }
    lcovdet=0
  }else{
    lcovdet=log(det(cov))
  }

  # nleqslv 또는 optim 함수를 사용하여 매개변수 추정을 시도합니다.
  tryCatch(
    for(i in 1:ntry){

      value=list()

      if(alg=='pmax' ){

        value <- try(
          nleqslv(x=as.vector(init[i,1:3]),fn=glme.sol,
                  slmgev=slmgev, eps=eps)
        )

      }else if(alg=='like'){

        value <- try(
          optim(par=as.vector(init[i,1:3]), fn=glme.like,
                xdat=xdat, slmgev=slmgev, covinv=covinv,
                lcovdet=lcovdet, lamb=lamb,
                pref=pref, p=p, q=q)
        )
      }

      if(is(value)[1]=="try-error"){
        k[[i]]$fvec <- 10^6
      }else{
        k[[i]] <- value
        if(alg=='pmax' ){
          k[[i]]$root = value$x
        }else if(alg=='like'){
          k[[i]]$root = value$par
          k[[i]]$fvec = value$value
        }
      }

      precis[i]= mean( abs(k[[i]]$fvec) )

      if(alg=='pmax' ){
        if( abs(value$termcd) > 2 ) precis[i]=1000
      }else{
        if( value$convergence != 0) precis[i]=10^6
      }

      precis[is.na(precis[i])]=10^6

      fgtol=ftol
      if(alg=='like') fgtol= 10^6

      if(precis[i] < fgtol){

        pk.ms[i]= MS_pk( k[[i]]$root, p=p, q=q )$pk.ms

        if(isol == 1 ) {
          sol[[isol]] = k[[i]]
          sol[[isol]]$pk.ms = pk.ms[i]
        }

        if( isol >= 2 ) {
          for(isolm1 in 1:(isol-1) ) {
            dist[isolm1]= sum( abs(k[[i]]$root[1:3] - sol[[isolm1]]$root[1:3]) )
          }
          mindist=min(dist)

          if( mindist > 0.001 ) {
            sol[[isol]] = k[[i]]
            sol[[isol]]$pk.ms = pk.ms[i]
          }
        }

        if( mindist > 0.001 ) isol=isol+1

      }else if(precis[i] >= fgtol) {
        pk.ms[i] = -100.0
      }

    } #for
  ) #tryCatch

  nsol=isol-1
  if(nsol == 0) {
    cat("-- No solution was found in nleqslv or optim --","\n")
  }else{
    z$nsol= nsol
    #    z$sol = sol
  }

  if(alg=='pmax'){
    selc_num = which.max( pk.ms )
  }

  if(alg=='like'){
    selc_num = which.min( abs(precis) )    #precis=k[[i]]$fvec
  }

  x  <- k[[selc_num]]

  if(alg=='pmax' ) z$estim.precis <- precis[selc_num]
  if(alg=='like') z$nllh.pref = k[[selc_num]]$fvec

  glme <- x$root
  z$pmax = pk.ms[selc_num]

  z$glme.conF = glme
  z$glme.conT = glme
  z$glme.bc = glme

  if(glme[3] < 0){
    #beta1 = -0.64766
    beta1 = -0.50766
    #beta1 = -0.40766
    z$glme.bc[3] = glme[3] - beta1*glme[3]/sqrt(nsample)
    if( abs(z$glme.bc[3]) > 0.5 ){
      z$glme.bc[3] = sign(z$glme.bc[3])*0.499999
    }
    glme.bc= pargev.kfix(lmom_init, kfix=z$glme.bc[3])$para
    z$glme.bc = glme.bc
  }

  if( abs(glme[3]) > 0.5 ){
    glme[3] = sign(glme[3])*0.499999
    z$glme.conT= pargev.kfix(lmom_init, kfix=glme[3])$para
  }

  z$lme.conF =lmom_est$para
  z$lme.conT =lmom_est$para

  if(abs(lmom_est$para[3]) > 0.5){
    lmom_est$para[3] = sign(lmom_est$para[3])*0.499999
    z$lme.conT= pargev.kfix(lmom_init, kfix=lmom_est$para[3])$para
  }

  z$algorithm = alg
  z$pref.func = pref

  if( alg=='pmax' ) {
    z$eps.pmax= eps
  }else if(alg=='like'){
    z$weight.like= lamb
  }
  z$nsample = nsample

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



