#----- Select the work directory
setwd("~/Progress/MC/Bases")


#----- Load data
load("data_listjm.RData")
setwd("~/Progress/MC/Modeles")

#----- create database
progress.long <- as.data.frame(cbind(ID = data_listjm$subject,
                                     pressure = data_listjm$y,
                                     obstime = data_listjm$time))

progress.id <- as.data.frame(cbind(ID = data_listjm$covariates[, 5],
                                   data_listjm$covariates[, -5]))

# survival data with fixed period for piecewise constant function
progress.surv <- as.data.frame(cbind(ID = data_listjm$survid,
                                     event = data_listjm$event,
                                     period = data_listjm$period,
                                     end = data_listjm$end,
                                     start = data_listjm$start,
                                     ethnie = data_listjm$ethnie,
                                     age = data_listjm$age,
                                     sex = data_listjm$sex,
                                     rand = data_listjm$rando))

# Get survival data without periods
library(dplyr)  

test <- progress.surv %>% group_by(ID) %>% slice(which.max(end))
progress.id$event <- test$event
progress.id$time <- test$end
# get an unique table
progress <- left_join(progress.long, progress.id, by = c("ID"))

# manage the type of variables

progress$Stroke <- factor(progress$event,
                          levels = c(0,1),
                          labels = c("No", "Yes"))

progress$subject <- as.integer(progress$ID)
progress$event <- as.integer(progress$event)

progress$ethnie <- factor(progress$ethnie,
                          levels = c(0,1),
                          labels = c("Asian", "Non_Asian"))

progress$sex <- factor(progress$sex,
                       levels = c(0,1),
                       labels = c("Female", "Male"))

progress$rand <- factor(progress$rand,
                        levels = c(0,1),
                        labels = c("Placebo","Treated"))


# delete work datasets
rm(progress.surv,progress.id,progress.long,test)

# traversal dataset (or survival dataset)
progress.id <- progress[!duplicated(progress$ID), ]

# keep only subjects receiving placebo
progress.p <- progress[which(progress$rand=="Placebo"), ]
progress.p.id <- progress.p[!duplicated(progress.p$ID), ]

# standardisation of quantitative covariates
progress.p$age_stand <- (progress.p$age-mean(progress.p.id$age))/sd(progress.p.id$age)
progress.p$visit1 <- ifelse(progress.p$obstime == 0, 1,0)

#---------------------------------------------------------------------------------------------------------------------------------------

library(survival)
library(jagsUI)

# Linear mixed model
# Y = X\beta + Ub

# survival model
# h(t) = h_0(t) exp( Z\alpha + AP) # AP:associations partagées
#-- Arguments (e.g. progress data)
formFixed = pressure ~ visit1 + sex + ethnie + age_stand + obstime + I(obstime^2)
formRandom = ~ obstime +  I(obstime^2)
formGroup = ~ subject
formSurv = Surv(time, event) ~ ethnie + sex + age_stand
survMod = "weibull"
param = "value"
timeVar= "obstime"
idVar = "subject"
data = progress.p
n.chains = 3
n.iter = 100000
n.burnin = 5000
n.thin = 15
n.adapt = 10000
precision = 0.01
parallel = TRUE
C = 1000


# data management

## longitudinal part

### longitudinal data
data_long <- data[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
y <- data_long[all.vars(formFixed)][, 1]
mfX <- model.frame(formFixed, data = data_long)
X <- model.matrix(formFixed, mfX)
mfU <- model.frame(formRandom, data = data_long)
U <- model.matrix(formRandom, mfU)
id <- as.integer(data_long[all.vars(formGroup)][,1])
offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
I <- length(unique(id))
if(!("id" %in% colnames(data_long)))
  data_long <- cbind(data_long, id = id)

### exemple of use lcmm packages to initiate parameter values of LMM
long_model <- lcmm::hlme(fixed = formFixed, 
                         random= formRandom,
                         subject = idVar, 
                         data=data)
priorMean.beta <- long_model$best[1:ncol(X)]
priorTau.beta <- diag(rep(precision,length(priorMean.beta)))
sigma <- long_model$best["stderr"]

### creation of a data list for jags
jags.data <- list(y = y,
                  X = X,
                  U = U,
                  ncX = ncol(X),
                  ncU = ncol(U),
                  #ncU2 = ncol(U) +1,
                  I = I,
                  offset = offset,
                  priorMean.beta = priorMean.beta,
                  priorTau.beta = priorTau.beta,
                  priorMean.b = c(rep(0, ncol(U))),
                  priorMean.log.sigma = log(sigma),
                  priorR.Sigma2 = diag(rep(1/precision, ncol(U))),
                  priorK.Sigma2 = ncol(U)
)


## survival part

### survival data
tmp <- data[c(all.vars(formGroup),all.vars(formSurv))]
tmp <- unique(tmp)
Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
event <- tmp[all.vars(formSurv)][, 2]   # vector of event indicator (delta)
nTime <- length(Time)                   # number of subject having Time
zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
# design matrice
mfZ <- model.frame(formSurv, data = tmp)
Z <- model.matrix(formSurv, mfZ)

### use survival::coxph function to initiate values
cat("> Initialisation of survival parameter values using 'survival' package. \n")
tmp_model <- survival::coxph(formSurv,
                             data = tmp,
                             x = TRUE)
priorMean.alpha <- c(0,tmp_model$coefficients)
priorTau.alpha <- diag(c(precision, precision*(1/diag(tmp_model$var))))

### Complete the jags data
jags.data <- c(jags.data,
               list(C = C,
                    precision = precision,
                    zeros = numeric(nTime),
                    Time = Time,
                    event = event,
                    Z = Z,
                    ncZ = ncol(Z),
                    priorMean.alpha = priorMean.alpha,
                    priorTau.alpha = priorTau.alpha,
                    priorTau.alphaA = 1/precision)
)

#--- shared current value case
lag = 0
data.id <- data_long[!duplicated(id), ]
if (!timeVar %in% names(data_long))
  stop("\n'timeVar' does not correspond to one of the columns in formulas")
if (param %in% c("value")) {
  data.id[[timeVar]] <- pmax(Time-lag, 0)
  mfX.id <- model.frame(formFixed, data = data.id)
  Xtime <- model.matrix(formFixed, mfX.id)
  mfU.id <- model.frame(formRandom, data = data.id)
  jags.data <- c(jags.data, list(Xtime = Xtime, Utime = Utime))
  
  #-- approxitmation of the intergral via the Gaussian quadrature (Gauss Kronrod rule)
  gaussKronrod <-
    function (k = 15) {
      sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
              0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
              -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
              0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
      wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
      wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
               0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
      if (k == 7)
        list(sk = sk[1:7], wk = wk7)
      else
        list(sk = sk, wk = wk15)
    }
  
  wk <- gaussKronrod()$wk
  sk <- gaussKronrod()$sk
  K <- length(sk)
  P <- Time/2
  st <- outer(P, sk + 1)
  id.GK <- rep(seq_along(Time), each = K)
  data.id2 <- data.id[id.GK, ]
  data.id2[[timeVar]] <- c(t(st))
  mfX <- model.frame(formFixed, data = data.id2)
  mfU <- model.frame(formRandom, data = data.id2)
  Xs <- model.matrix(formFixed, mfX)
  Us <- model.matrix(formRandom, mfU)
  
  jags.data <- c(jags.data, list(K = K, P = P, st = st, wk = wk, Xs = Xs, Us = Us))
}




## parameters to save in the sampling step
parms_to_save <- c("alpha", "alpha.beta","alpha.sigma", "beta", "covariance.b", "mu.log.sigma", "sigma2.log.sigma", "priorMean.b","alpha.assoc")
## call JAGS using jagsUI
out_jags = jagsUI::jags(data = jags.data,
                        parameters.to.save = parms_to_save  ,
                        model.file = "barrett.weibull.valeurcouranteJags.txt",
                        #inits = param.inits,
                        n.chains = n.chains,
                        parallel = parallel,
                        n.adapt = n.adapt,
                        n.iter = n.iter,
                        n.burnin = n.burnin,
                        n.thin = n.thin,
                        DIC = T)

setwd("~/Progress/MC/Resultats")

#resultats <- list( initialisation = param.inits, sortie = out_jags)

save(out_jags, file = "Res_MA4_cor_valcour.RData")

