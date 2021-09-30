##### Packages
library(survival)
library(jagsUI)
library(MASS)
library(randtoolbox)

source("functions.RData")

#-- Arguments 
formFixed = pressure ~  obstime + I(obstime^2)  + sex + ethnie + age_stand + sex*obstime + ethnie*obstime + age_stand*obstime
formRandom = ~ obstime +  I(obstime^2)
formSlopeFixed = ~ I(2*obstime) + sex + ethnie + age_stand
indices_beta_slope <- c(2,3,7,8,9)
formSlopeRandom = ~ I(2*obstime)
formGroup = ~ subject
formSurv = Surv(time, event) ~ sex + ethnie + age_stand
survMod = "weibull"
param = "value"
timeVar= "obstime"
idVar = "subject"
data = progress.p
n.chains = 3
n.iter = 200000
n.burnin = 100000
n.thin = 20
n.adapt = 50000
precision = 0.01
parallel = TRUE
C = 1000


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


### list for jags

jags.data <- list(y = list.long$y,
                  X = list.long$X,
                  U = list.long$U,
                  ncX = ncol(list.long$X),
                  ncU = ncol(list.long$U),
                  I = list.long$I,
                  offset = list.long$offset,
                  priorMean.beta = list.init.long$priorMean.beta,
                  priorTau.beta = list.init.long$priorTau.beta,
                  priorMean.b = c(rep(0, ncol(list.long$U))),
                  priorMean.log.sigma = log(list.init.long$sigma),
                  priorR.Sigma2 = diag(rep(1/precision, ncol(list.long$U))),
                  priorK.Sigma2 = ncol(list.long$U),
                  indices_beta_slope = indices_beta_slope,
                  C = C,
                  precision = precision,
                  zeros = numeric(list.surv$nTime),
                  Time = list.surv$Time,
                  event = list.surv$event,
                  Z = list.surv$Z,
                  ncZ = ncol(list.surv$Z),
                  priorMean.alpha = list.init.surv$priorMean.alpha,
                  priorTau.alpha = list.init.surv$priorTau.alpha,
                  K = list.GaussKronrod$K,
                  P = list.GaussKronrod$P,
                  st = list.GaussKronrod$st,
                  wk = list.GaussKronrod$wk,
                  Xtime = list.data.time$Xtime,
                  Utime = list.data.time$Utime,
                  Xslope = list.data.slope.time$Xtime,
                  Uslope = list.data.slope.time$Utime,
                  ncXslope = ncol(list.data.slope.time$Xtime),
                  ncUslope = ncol(list.data.slope.time$Utime),
                  Xs = list.data.GK$Xtime,
                  Us = list.data.GK$Utime,
                  Xs.slope = list.data.GK.slope$Xtime,
                  Us.slope = list.data.GK.slope$Utime)


## parameters to save in the sampling step
parms_to_save <- c("alpha", "alpha.sigma", "beta", "covariance.b", "mu.log.sigma", "sigma2.log.sigma", "alpha.current", "alpha.slope", "shape")
## call JAGS using jagsUI

out_jags = jagsUI::jags(data = jags.data,
                        parameters.to.save = parms_to_save  ,
                        model.file = "VC.Slope.sigma.txt",
                        #inits = param.inits,
                        n.chains = n.chains,
                        parallel = parallel,
                        n.adapt = n.adapt,
                        n.iter = n.iter,
                        n.burnin = n.burnin,
                        n.thin = n.thin,
                        DIC = T)

save(out_jags, file = "VC.Slope.sigma.17092021.RData")