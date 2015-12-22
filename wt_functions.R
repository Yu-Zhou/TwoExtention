# Weighted/composite quantile regression
# Goal: 0.2quartile + 0.8 median
library(rgenoud)
library(quantreg)
library(xtable)
library(faraway)
library(Matrix)
library(lpSolve)
library(parallel)

compo_taus = rbind(c(0.25,0.5),c(0.25,0.5),c(0.10,0.30),c(0.10,0.30))
#nrow(compo_taus)
compo_wts = rbind(c(0.2,0.8),c(0.8,0.2),c(0.2,0.8),c(0.8,0.2))
#nrow(compo_wts)

sum_names <- c("intercept","X1","X2","hatQ")
#eta_range <- 1:3
sq<-rep(0,3)
nvars<-length(sq)
Domains<-cbind(rep(-1,nvars),rep(1,nvars)) 
pop.size=3000
it.num=3
p_level=0
hard_limit=FALSE

Perf_for_eta_Weighted_ver <- function(eta, quantile_levels, setting_num, compo_taus, compo_wts){ 
  x1p <- runif(10000)
  x2p <- runif(10000)
  xp<-cbind(1,x1p,x2p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- y_function(x1=x1p, x2=x2p,a=g)
  m=mean(yg)
  quantile = c()
  for (i in 1:length(quantile_levels)){   
    quantile = c(quantile, quantile(yg,quantile_levels[i]))
  }
  compo = c()
  for (k in 1:nrow(compo_taus)){
    quant_compok <- myhogg(x=data.frame(runif(10000)*0.01), y=yg, taus = compo_taus[k,],
                     weights = compo_wts[k,],
                     R=cbind(c(1,-1,0,0),c(-1,1,0,0),c(0,0,1,-1)),
                     r=c(0,0,0,0),eps=1e-14)    
    compo = c(compo, quant_compok$coeff[1])
  }
  return( list(value=c(m,quantile,compo), 
               Crit = c('mean',quantile_levels,rep('compo',setting_num)))
  )
}
myhogg = function (x, y, taus , weights ,R = NULL, r = NULL, beta = 0.99995, eps = 1e-06) 
{
  n <- length(y)
  n2 <- nrow(R)
  m <- length(taus)
  p <- ncol(x) + m
  if (n != nrow(x)) 
    stop("x and y don't match n")
  if (m != length(weights)) 
    stop("taus and weights differ in length")
  if (any(taus < eps) || any(taus > 1 - eps)) 
    stop("taus outside (0,1)")
  W <- diag(weights)
  if (m == 1) 
    W <- weights
  x <- as.matrix(x)
  X <- cbind(kronecker(W, rep(1, n)), kronecker(weights, x))
  y <- kronecker(weights, y)
  rhs <- c(weights * (1 - taus) * n, sum(weights * (1 - taus)) * 
             apply(x, 2, sum))
  #  if (n2 != length(r)) 
  #    stop("R and r of incompatible dimension")
  #  if (ncol(R) != p) 
  #    stop("R and X of incompatible dimension")
  d <- rep(1, m * n)
  u <- rep(1, m * n)
  if (length(r)) {
    wn1 <- rep(0, 10 * m * n)
    wn1[1:(m * n)] <- 0.5
    wn2 <- rep(0, 6 * n2)
    wn2[1:n2] <- 1
    z <- .Fortran("rqfnc", as.integer(m * n), as.integer(n2), 
                  as.integer(p), a1 = as.double(t(as.matrix(X))), c1 = as.double(-y), 
                  a2 = as.double(t(as.matrix(R))), c2 = as.double(-r), 
                  rhs = as.double(rhs), d1 = double(m * n), d2 = double(n2), 
                  as.double(u), beta = as.double(beta), eps = as.double(eps), 
                  wn1 = as.double(wn1), wn2 = as.double(wn2), wp = double((p + 
                                                                             3) * p), it.count = integer(3), info = integer(1), 
                  PACKAGE = "quantreg")
  }
  else {
    wn <- rep(0, 10 * m * n)
    wn[1:(m * n)] <- 0.5
    z <- .Fortran("rqfnb", as.integer(m * n), as.integer(p), 
                  a = as.double(t(as.matrix(X))), c = as.double(-y), 
                  rhs = as.double(rhs), d = as.double(d), as.double(u), 
                  beta = as.double(beta), eps = as.double(eps), wn = as.double(wn), 
                  wp = double((p + 3) * p), it.count = integer(2), 
                  info = integer(1), PACKAGE = "quantreg")
  }
  if (z$info != 0) 
    warning(paste("Info = ", z$info, "in stepy: singular design: iterations ", 
                  z$it.count[1]))
  coefficients <- -z$wp[1:p]
  if (any(is.na(coefficients))) 
    stop("NA coefs:  infeasible problem?")
  list(coefficients = coefficients, nit = z$it.count, flag = z$info)
}
Mean_Est<-function(eta,x,y,a,prob){
  #mean estimator
  g<-as.numeric(I(x%*%eta>0))
  c<-a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  val <- mean(c*y*wts)
  val
}
Quant_Est <- function(eta,x,y,a,prob,tau){
  #quantile estimator
  g <- as.numeric( I(x %*% eta > 0))
  c <- a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  wts <- c*wts
  model <- rq(y ~ 1, weights=wts, tau=tau)
  coefficients(model)[1]
}
Qweighed_Est <- function(eta,x,y,a,prob,taus,taus_wts){
  g <- as.numeric( I(x %*% eta > 0))
  consis <- a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  wts <- consis*wts
  
  # generate pseudo outcome
  ys = y*wts
  xs = wts
  NewD = cbind(consis, ys,xs)
  NewD = NewD[which(NewD[,1]==1),]
  xs = data.frame(NewD[,3])
  modelhogg = myhogg(xs,y=NewD[,2], taus=taus,weights = taus_wts,
                     R = cbind(c(1,-1,0,0),c(0,0,1,-1),rep(0,4)), r=c(0,0,0,0))
  return(modelhogg$coef[3])
}


#### generative function, sd = 1 #####
y_function <- function(x1,x2,a){
  error <- rnorm(length(x1), sd=1)
  y = 1 + x1 + x2 + 
      3*a -2.5*a*x1 -2.5*a*x2 +
      (1+a+a*x1+a*x2)*error
  return(y)
}
qestimate<-function(x,y,a,prob,tau,quantile_levels,compo_taus,compo_wts){
  nvars<-length(sq)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Quant_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,tau=tau,print.level=0,max=TRUE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=sq,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0))==1) return(c(rep(NA,4+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf<-Perf_for_eta_Weighted_ver(eta, quantile_levels=quantile_levels,setting_num=4,compo_taus=compo_taus,compo_wts=compo_wts)  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names,perf$Crit)
  
  return(summary)
}
mestimate<-function(x,y,a,prob,tau,quantile_levels,compo_taus,compo_wts){
  nvars<-length(sq)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Mean_Est,nvars=nvars,x=x,y=y,a=a,
              prob=prob,print.level=0,max=TRUE,pop.size=pop.size,
              wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, 
              P6=50, P7=50, P8=50, P9=0,Domains=Domains,
              starting.values=sq,solution.tolerance=0.00001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0))==1) return(c(rep(NA,4+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf<-Perf_for_eta_Weighted_ver(eta, quantile_levels=quantile_levels,setting_num=4,compo_taus=compo_taus,compo_wts=compo_wts)  
  #perf <- Perf_for_eta_Weighted_ver(eta, quantile_levels=quantile_levels )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names,perf$Crit)
  return(summary)
}

Qweighted <- function(x,y,a,prob,taus,taus_wts,
                  p_level, quantile_levels, compo_taus,compo_wts){
  nvars<-length(sq) 
  Domains<-cbind(rep(-1,nvars),rep(1,nvars)) #?nvars
  
  est<-genoud(fn=Qweighed_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,taus=taus, taus_wts=taus_wts,
              print.level=p_level,max=TRUE,pop.size=pop.size,wait.generations=it.num,
              gradient.check=FALSE, BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              Domains=Domains,starting.values=sq,
              hard.generation.limit=hard_limit,solution.tolerance=0.00001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0))==1) return(c(rep(NA,4+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf<-Perf_for_eta_Weighted_ver(eta, quantile_levels=quantile_levels,setting_num=4,compo_taus=compo_taus,compo_wts=compo_wts)  
  #perf<-Perf_for_eta_Weighted_ver(eta, quantile_levels=quantile_levels )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names,perf$Crit)
  #names(summary)<-c(sum_names,"mean", paste0(quantile_levels,".Q"))
  return(summary)
  
}

########
sim_function <- function(sim_num=2, nobs ,quantile_levels,compo_taus, compo_wts, p_level=0, hard_limit=FALSE){
  mean_data <- Q1_data <- Q2_data <- Q3_data <- Q4_data<-sim_data_weighted<- NULL
  sim_data<-list() #store original data
  
  cat("COMPOESITE QUANTILE/ WEIGHTED SUM OF QUANTILE LOSS")
  #######################  data generation #######################
  for(i in 1:sim_num){
    ptm0 <- proc.time()
    ##
    x1 <- runif(nobs)
    x2 <- runif(nobs)
    x=cbind(1,x1,x2)
    tp <- exp(-0.5+x1+x2)/(1+exp(-0.5+x1+x2))  
    a<-rbinom(nobs,1,tp)                                                        
    y <- y_function(x1,x2,a)
    sim_data = c(sim_data,list(x=x, y=y, a=a))
    ############################## Propensity score model    ########################
    logit<-glm(a~x1+x2,family=binomial, epsilon=1e-14)
    ph.true<- as.vector(logit$fit)

    mean_summary <- mestimate(x,y,a,ph.true, quantile_levels=quantile_levels,
                              compo_taus=compo_taus,compo_wts=compo_wts)
    #summary_10 <- qestimate(x,y,a,prob = ph.true,tau=.10, quantile_levels=quantile_levels)
    summary_1 <-qestimate(x,y,a,prob = ph.true, tau=quantile_levels[1], quantile_levels=quantile_levels, compo_taus,compo_wts)
    summary_2 <-qestimate(x,y,a,prob = ph.true, tau=quantile_levels[2], quantile_levels=quantile_levels, compo_taus,compo_wts)
    summary_3 <-qestimate(x,y,a,prob = ph.true, tau=quantile_levels[3], quantile_levels=quantile_levels, compo_taus,compo_wts)
    summary_4 <-qestimate(x,y,a,prob = ph.true, tau=quantile_levels[4], quantile_levels=quantile_levels, compo_taus,compo_wts)
    
    Qweighted_result<- NULL
    for (k in 1:nrow(compo_taus)){
        tmpk <-Qweighted(x,y,a,prob=ph.true,taus=compo_taus[k,],taus_wts=compo_wts[k,],
                        p_level, quantile_levels=quantile_levels, compo_taus, compo_wts)
        Qweighted_result <- rbind(Qweighted_result,c(tmpk,setting=k))
    }
    
    mean_data <- rbind(mean_data, mean_summary)
    #Q10_data <- rbind(Q10_data,summary_10)
    Q1_data <- rbind(Q1_data,summary_1)
    Q2_data <- rbind(Q2_data,summary_2)
    Q3_data <- rbind(Q3_data,summary_3)
    Q4_data <- rbind(Q4_data,summary_4)
    sim_data_weighted <- rbind(sim_data_weighted, Qweighted_result)
    
    ptm1=proc.time() - ptm0
    jnk=as.numeric(ptm1[3])
    cat('It took ', jnk, "seconds,  #### Current working dir: ", i,'. \n')
  }
  
  return(list(mean_data, 
              Q1_data, Q2_data, Q3_data, Q4_data,
              sim_data_weighted,sim_data))
}




