################## Functions ##################

# Data management
##longitudinal part

data.manag.long <- function(formGroup, formFixed, formRandom, data){
  
  data_long <- data[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
  y <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  id <- as.integer(data_long[all.vars(formGroup)][,1])
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  I <- length(unique(id))
  if(!("id" %in% colnames(data_long))) #To have a column named "id"
    data_long <- cbind(data_long, id = id)
  
  list.long <- list("data_long"= data_long, "y" = y, "X" = X, "U" = U,
                    "id" = id, "offset"=offset, "I" = I)
  
  return(list.long)
}

initial.long <- function(formFixed, formRandom, idVar, data, precision, ncX){
  
  long_model <- lcmm::hlme(fixed = formFixed, 
                           random= formRandom,
                           subject = idVar, 
                           data=data)
  priorMean.beta <- long_model$best[1:ncX]
  priorTau.beta <- diag(rep(precision,length(priorMean.beta)))
  sigma <- long_model$best["stderr"]
  
  list.init.long <- list("long_model" = long_model, "priorMean.beta" = priorMean.beta,
                         "priorTau.beta" = priorTau.beta, "sigma" = sigma)
  
  return(list.init.long)
  
}


## Survival part
data.manag.surv <- function(formGroup, formSurv, data){
  
  tmp <- data[c(all.vars(formGroup),all.vars(formSurv))]
  tmp <- unique(tmp)
  Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  event <- tmp[all.vars(formSurv)][, 2]   # vector of event indicator (delta)
  nTime <- length(Time)                   # number of subject having Time
  zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
  # design matrice
  mfZ <- model.frame(formSurv, data = tmp)
  Z <- model.matrix(formSurv, mfZ)
  
  list.surv <- list("tmp" = tmp, "Time" = Time , "event"= event, "nTime" = nTime,
                    zeros = "zeros", "Z" = Z)
  
  return(list.surv)
  
}

initial.surv <- function(formSurv, data_surv, precision){
  
  cat("> Initialisation of survival parameter values using 'survival' package. \n")
  tmp_model <- survival::coxph(formSurv,
                               data = data_surv,
                               x = TRUE)
  
  priorMean.alpha <- c(0,tmp_model$coefficients)
  priorTau.alpha <- diag(c(precision, precision*(1/diag(tmp_model$var))))
  
  list.init.surv <- list("tmp_model" = tmp_model, 
                         "priorMean.alpha" = priorMean.alpha, 
                         "priorTau.alpha" = priorTau.alpha)
  
  return(list.init.surv)
}

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



data.time <- function(data.id, Time, formFixed, formRandom){
  
  if (!timeVar %in% names(data.id))
    stop("\n'timeVar' does not correspond to one of the columns in formulas")
  if (param %in% c("value")) {
    data.id[[timeVar]] <- Time
    mfX.id <- model.frame(formFixed, data = data.id)
    Xtime <- model.matrix(formFixed, mfX.id)
    mfU.id <- model.frame(formRandom, data = data.id)
    Utime <- model.matrix(formRandom, mfU.id)
    
    list.data.time <- list("Xtime" = Xtime, "Utime" = Utime)
  }
  return(list.data.time)
}


data.GaussKronrod <- function(data.id, Time, k = 15){
  
  wk <- gaussKronrod()$wk
  sk <- gaussKronrod()$sk
  K <- length(sk)
  P <- Time/2
  st <- outer(P, sk + 1)
  id.GK <- rep(seq_along(Time), each = K)
  data.id2 <- data.id[id.GK, ]
  
  list.GaussKronrod <- list(K = K, P = P, st = st, wk = wk, data.id2 = data.id2,
                            id.GK = id.GK)
  return(list.GaussKronrod)
 
}  









  
