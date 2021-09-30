#----- Select the work directory

source("functions.R")
library(survival)
library(MASS)
library(randtoolbox)
library(marqLevAlg)

# Linear mixed model
# Y = X\beta + Ub

# survival model
# h(t) = h_0(t) exp( Z\alpha + AP) # AP:associations partagées
#-- Arguments (e.g. progress data)
formFixed = pressure ~  obstime + I(obstime^2) + sex + ethnie + age_stand + sex*obstime + ethnie*obstime + age_stand*obstime
formRandom = ~ I(obstime) + I(obstime^2) 
formSlopeFixed = ~ I(2*obstime) + sex + ethnie + age_stand
indices_beta_slope <- c(2,3,7,8,9)
formSlopeRandom = ~ I(2*obstime)
formGroup = ~ subject
formSurv = Surv(time, event) ~ ethnie + sex + age_stand
survMod = "weibull"
param = "value"
timeVar= "obstime"
idVar = "subject"
data = progress.p
precision <- 0.01

# data management

## longitudinal part

list.long <- data.manag.long(formGroup, formFixed, formRandom, data)

list.init.long <- initial.long(formFixed, formRandom, idVar, data, precision,
                               ncol(list.long$X))

## survival part

list.surv <- data.manag.surv(formGroup, formSurv, data)

list.init.surv <- initial.surv(formSurv, list.surv$tmp, precision)


#--- shared current value case
lag = 0
data.id <- list.long$data_long[!duplicated(list.long$id), ]

list.data.time <- data.time(data.id, list.surv$Time, formFixed, formRandom)

list.data.slope.time <- data.time(data.id, list.surv$Time, formSlopeFixed, 
                                  formSlopeRandom)

list.GaussKronrod <- data.GaussKronrod(data.id, list.surv$Time, k = 15)

list.data.GK <- data.time(list.GaussKronrod$data.id2, c(t(list.GaussKronrod$st))
                          , formFixed, formRandom)

list.data.GK.slope <- data.time(list.GaussKronrod$data.id2, 
                                c(t(list.GaussKronrod$st)), formSlopeFixed, 
                                formSlopeRandom)


S <- 1000
nb.e.a <- 4
indices_beta_slope <- indices_beta_slope
Xtime <- list.data.time$Xtime
Utime <- list.data.time$Utime 
Xslope <- list.data.slope.time$Xtime
Uslope <- list.data.slope.time$Utime
Time <- list.surv$Time
Z <- list.surv$Z
Xs <- list.data.GK$Xtime
Us <- list.data.GK$Utime 
Xs.slope <- list.data.GK.slope$Xtime
Us.slope <- list.data.GK.slope$Utime
wk <- list.GaussKronrod$wk
st <- list.GaussKronrod$st
P <- list.GaussKronrod$P
X_base <- list.long$X   
U <- list.long$U
event <- list.surv$event
y <- list.long$y
offset <- list.long$offset
priorR.Sigma2 = diag(rep(1/precision, ncol(list.long$U)))
priorK.Sigma2 = ncol(list.long$U)
priorMean.beta = list.init.long$priorMean.beta
Ind <- list.long$I


log_density_contribution <- function(param,S, nb.e.a,  Xtime, Utime, indices_beta_slope, Xslope , Uslope ,
                        shape, Time, Z, Xs, Us, Xs.slope , Us.slope ,
                        wk, st_calc, P,X_base,U, event, y, Zq){
  
  C <- matrix(rep(0,9),nrow = nb.e.a-1, ncol = nb.e.a-1)
  C[1,] <- param[1:3]
  C[2,2:3] <- param[4:5]
  C[3,3] <- param[6]
  mu.log.sigma <- param[7]
  sigma2.log.sigma <- param[8]**2
  beta <- param[9:17]
  alpha <- param[18:21]
  alpha.sigma <- param[22]
  alpha.current <- param[23]
  alpha.slope <- param[24]
  shape <- param[25]
  
  
  b_al <- Zq[,1:(nb.e.a-1)]%*%C
  log.sigma_al <- matrix(rep(mu.log.sigma,S), nrow = S, byrow = T) + Zq[,nb.e.a]*sqrt(sigma2.log.sigma)
  
  CV <- (beta%*%Xtime)[1,1] + b_al%*%Utime
  slope <- (beta[indices_beta_slope]%*%Xslope)[1,1]+ b_al[,-1]%*%Uslope
  
  h <- shape*Time**(shape-1)*exp((alpha%*%Z)[1,1] + alpha.sigma*exp(log.sigma_al) + alpha.current*CV + alpha.slope*slope)
  
  current.GK <- matrix(rep(beta%*%t(Xs),S),nrow = S, byrow = T) +b_al%*%t(Us)
  slope.GK <- matrix(rep(beta[indices_beta_slope]%*%t(Xs.slope),S), nrow = S, byrow = TRUE) + b_al[,-1]%*%t(Us.slope)
  SurvLong <- exp(alpha.current*current.GK + alpha.slope*slope.GK)%*%(wk*shape*(st_calc**(shape-1)))

  etaBaseline <- (alpha%*%Z)[1,1] + alpha.sigma*exp(log.sigma_al)
  Surv <- exp((-exp(etaBaseline) * P * SurvLong))
  
  
  
  if(is.null(nrow(X_base))){
    CV <- (beta%*%X_base)[1,1] + b_al%*%U
    f_Y_b_sigma <- dnorm(x = y, mean = CV , sd = exp(log.sigma_al))
  }
  else{
    f_Y_b_sigma <- rep(1,S)
    for (k in 1:nrow(X_base)){
      CV <- (beta%*%X_base[k,])[1,1] + b_al%*%U[k,]
      f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y[k], mean = CV , sd = exp(log.sigma_al))
    }
  }
  
  log_dens_glob <- log(sum((h**event)*Surv*f_Y_b_sigma))-log(S)
  
  return(log_dens_glob )
}

log_vraisemblance <- function(param,S, nb.e.a,  Xtime, Utime, indices_beta_slope , Xslope , Uslope ,
                      Time, Z, Xs, Us, Xs.slope, Us.slope ,
                       wk, st_calc, P,X_base,U, event, y, Zq, offset, Ind, fonction){
  
  ll_glob <- 0
  for(i in 1:Ind){
    log_dens <- fonction(param,S = S, nb.e.a = nb.e.a, Xtime = Xtime[i,], Utime = Utime[i,], 
                            indices_beta_slope = indices_beta_slope, Xslope = Xslope[i,], 
                            Uslope = Uslope[i,],Time = Time[i], 
                            Z = Z[i,], Xs = Xs[(15*(i-1)+1):(15*i),], Us = Us[(15*(i-1)+1):(15*i),], 
                            Xs.slope = Xs.slope[(15*(i-1)+1):(15*i),], Us.slope = Us.slope[(15*(i-1)+1):(15*i),], 
                            wk = wk, st_calc = st_calc[i,], P = P[i],X_base = X_base[offset[i]:(offset[i+1]-1),],U = U[offset[i]:(offset[i+1]-1),], 
                            event = event[i], y = y[offset[i]:(offset[i+1]-1)], Zq = Zq)
    ll_glob <- ll_glob + log_dens
    
  }
  return(ll_glob)
  
}



Zq <- randtoolbox::sobol(S, nb.e.a, normal = TRUE, scrambling = 1)
matr_cov <- rWishart(1,priorK.Sigma2,priorR.Sigma2)[,,1]
C_chol <- chol(matr_cov)
mu <- rnorm(1, mean = log(list.init.long$sigma), sd = 1 )
tau <- rgamma(1,0.1,0.1)

progress.p.id <- progress.p[!duplicated(progress.p$ID),]
mod_weib  <- survreg(formSurv, data = progress.p.id, dist = "weibull")

binit <- c(C_chol[1,], C_chol[2,2:3], C_chol[3,3],0, sqrt(0.1), priorMean.beta, mod_weib$coefficients, 0,
           0,0,1/mod_weib$scale)

estimation_freq <- marqLevAlg(binit, fn = log_vraisemblance, minimize = FALSE, S = S, nb.e.a = nb.e.a,  Xtime = Xtime, Utime = Utime,
                         indices_beta_slope =indices_beta_slope, Xslope =Xslope, Uslope = Uslope,
                         Time = Time, Z = Z, Xs = Xs, Us = Us, Xs.slope = Xs.slope, Us.slope = Us.slope,
                         wk = wk, st_calc = st, P = P ,X_base = X_base,U = U, event = event, y = y , Zq = Zq,
                         offset = offset, Ind = list.long$I, fonction = log_density_contribution,
                    nproc = 10, clustertype = "SOCK", maxiter = 1000, print.info = TRUE, 
                    file = "Marq_Info_Slope.txt",blinding = TRUE, epsa = 1e-03,
                    epsb = 1e-03,
                    epsd = 1e-01)












