# Two Covariates Model, expanded to Three Covariates Class
# trade off 3 to 2
# trade off
# expect it to increase bias, and decrease variance
# Goal: 0.25 and 0.50 quantiles
# sd :1

library(utils)
library(rgenoud)
library(quantreg)
library(faraway)
library(Matrix)
library(lpSolve)
library(parallel)

sum_names3 <- c("intercept","X1","X2","X3","hatQ")
sum_names2 <- c("intercept","X1","X2","hatQ")
sq3<-rep(0,4)
sq2<-rep(0,3)
pop.size=1000; it.num=2
p_level=1; hard_limit=FALSE



Perf_for_eta_2cov <- function(eta,quantile_levels,Cnobs){ 
  n=100000
  x1p <- runif(n)
  x2p <- runif(n)
  xp<-cbind(1,x1p,x2p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- y_function(x1=x1p, x2=x2p,a=g)
  m=mean(yg)
  quantile = c()
  for (i in 1:length(quantile_levels)){   quantile = c(quantile, quantile(yg,quantile_levels[i]))}
  
  n=3000
  x1p <- runif(n)
  x2p <- runif(n)
  xp<-cbind(1,x1p,x2p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- y_function(x1=x1p, x2=x2p,a=g)
  plot(x1p,x2p,col = g+1)
  sum=0
  COMB = data.frame(t(combn(x=1:n,2)))
  COMB$absd = abs(yg[COMB$X1] - yg[COMB$X2])
  S = sum(COMB$absd)
  
#   for (k in 1:dim(COMB)[2]){
#     i = COMB[1,k]; j = COMB[2,k]
#     sum = sum + abs(yg[i]-yg[j])
#   }
  val <- S/(n*(n-1)/2)

  return( list(value=c(m,quantile,val), Crit = c('mean',quantile_levels,'Exp.AbsoDiff')))
}

Perf_for_eta_3cov<- function(eta,quantile_levels){ 
  # eta has 4 parameters
  x1p <- runif(100000)
  x2p <- runif(100000)
  x3p <- runif(100000)
  xp<-cbind(1,x1p,x2p,x3p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- y_function(x1=x1p, x2=x2p,a=g)
  m=mean(yg)
  quantile = c()
  for (i in 1:length(quantile_levels)){   quantile = c(quantile, quantile(yg,quantile_levels[i]))}
  
  n=1000
  x1p <- runif(n)
  x2p <- runif(n)
  x3p <- runif(n)
  xp<-cbind(1,x1p,x2p,x3p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- y_function(x1=x1p, x2=x2p,a=g)
  sum=0
  COMB = C1000
  for (k in 1:dim(COMB)[2]){
    i = COMB[1,k]; j = COMB[2,k]
    sum = sum + abs(yg[i]-yg[j])
  }
  val <- sum/(n*(n-1)/2)
  return( list(value=c(m,quantile,val), Crit = c('mean',quantile_levels,'Exp.AbsoDiff')))
}

#### generative function #####
y_function <- function(x1,x2,a){
  error <- rnorm(length(x1), sd = 1)
  y = 1 + x1 + x2 + 3*a -2.5*a*x1 -2.5*a*x2+
    (1+a+a*x1+a*x2)*error
  return(y)
}

qestimate3cov<-function(x,y,a,prob,tau,quantile_levels){
  nvars<-length(sq3)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Quant_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,tau=tau,print.level=0,max=TRUE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=sq3,solution.tolerance=0.00001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0,0))==1) return(c(rep(NA,5+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf<-Perf_for_eta_3cov(eta, quantile_levels=quantile_levels )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names3,"mean", paste0(quantile_levels,".Q"))
  
  return(summary)
}
mestimate3cov<-function(x,y,a,prob,tau,quantile_levels){
  nvars<-length(sq3)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Mean_Est,nvars=nvars,x=x,y=y,a=a,
              prob=prob,print.level=0,max=TRUE,pop.size=pop.size,
              wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, 
              P6=50, P7=50, P8=50, P9=0,Domains=Domains,
              starting.values=sq3,solution.tolerance=0.00001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0,0))==1) return(c(rep(NA, 5+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf <- Perf_for_eta_3cov(eta, quantile_levels=quantile_levels )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names3,"mean", paste0(quantile_levels,".Q"))
  return(summary)
}
qestimate2cov<-function(x,y,a,prob,tau,quantile_levels,Cnobs){
  nvars<-length(sq2)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Quant_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,tau=tau,
              print.level=0,max=TRUE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              Domains=Domains,starting.values=sq2,solution.tolerance=0.00001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0))==1) return(c(rep(NA,5+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf<-Perf_for_eta_2cov(eta, quantile_levels=quantile_levels,Cnobs )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names2,"mean", paste0(quantile_levels,".Q"))
  
  return(summary)
}
mestimate2cov<-function(x,y,a,prob,tau,quantile_levels){
  nvars<-length(sq2)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Mean_Est,nvars=nvars,x=x,y=y,a=a,
              prob=prob,print.level=0,max=TRUE,pop.size=pop.size,
              wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, 
              P6=50, P7=50, P8=50, P9=0,Domains=Domains,
              starting.values=sq2,solution.tolerance=0.00001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0))==1) return(c(rep(NA, 5+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf <- Perf_for_eta_2cov(eta, quantile_levels=quantile_levels,Cnobs )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names2,"mean", paste0(quantile_levels,".Q"))
  return(summary)
}

AbsoDiff<-function(x,y,a,prob,quantile_levels, Cnobs){
  nvars<-length(sq2)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Abso_Est,nvars=nvars,x=x,y=y,a=a,prob=prob, Cnobs= Cnobs,
              print.level=p_level,max=FALSE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              Domains=Domains,starting.values=sq2,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  if (prod(eta==c(0,0,0))==1) return(c(rep(NA,5+length(quantile_levels)))) 
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf<-Perf_for_eta_2cov(eta, quantile_levels=quantile_levels,Cnobs )  
  summary<-c(eta,hatQ,perf$value)
  names(summary)<-c(sum_names2,"mean", paste0(quantile_levels,".Q"),"Abso")
  
  return(summary)
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

Abso_Est<-function(eta,x,y,a,prob, Cnobs){
  g<-as.numeric(I(x%*%eta>0))
  c<-a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  n=length(y)
  
  Cnobs = data.frame(t(Cnobs))
  Cnobs$wtsprod = wts[Cnobs$X1] * wts[Cnobs$X2]
  Cnobs$cprod = c[Cnobs$X1] * c[Cnobs$X2]
  Cnobs$absd = abs(y[Cnobs$X1] - y[Cnobs$X2])
  S = sum(Cnobs$cprod * Cnobs$wtsprod * Cnobs$absd)
    
#   for (k in 1:dim(Cnobs)[2]){
#     i = Cnobs[1,k]; j = Cnobs[2,k]
#     sum = sum + c[i]*c[j]*wts[i]*wts[j]*abs(y[i]-y[j])
#   }
  
  val <- S/(n*(n-1)/2)
  return(val)
}

########
sim_function <- function(sim_num=2, nobs ,taus, p_level=0, hard_limit=FALSE){
  cat('TRADE OFF EXPERIEMNT: TRUE 2COVS, PRETEND IT HAS 3 COVS')
  D_mean_truep_3cov <-D_mean_truep_2cov <-D_mean_fp_3cov<-D_mean_fp_2cov<-NULL
  D_Q1_truep_3cov <-D_Q1_truep_2cov <-D_Q1_fp_3cov <- D_Q1_fp_2cov <- NULL
  D_Q2_truep_3cov <-D_Q2_truep_2cov <-D_Q2_fp_3cov <-D_Q2_fp_2cov<- NULL
  D_abso<-NULL
  sim_data<-list() #store original data
  quantile_levels = taus
  #######################  data generation #######################
  for(i in 1:sim_num){
    ptm0 <- proc.time()
    ##
    x1 <- runif(nobs)
    x2 <- runif(nobs)
    x3 <- runif(nobs)
    x2cov=cbind(1,x1,x2)
    x = cbind(1,x1,x2,x3)
    
    tp <- exp(-0.5+x1+x2)/(1+exp(-0.5+x1+x2))  
    a<-rbinom(nobs,1,tp)                                                        
    y <- y_function(x1,x2,a)
    sim_data = c(sim_data,list( x = x, y=y, a=a))
    ############################## Propensity score model    ########################
    logit<-glm(a~x1+x2,family=binomial, epsilon=1e-14)
    ph.true<- as.vector(logit$fit)
    logit3cov<-glm(a~x1+x2+x3,family=binomial, epsilon=1e-14)
    ph.3cov<- as.vector(logit3cov$fit)
    
    mean_truep_3cov <- mestimate3cov(x,    y,a,ph.true,quantile_levels=quantile_levels)
    #mean_fp_3cov <-    mestimate3cov(x,    y,a,ph.3cov,quantile_levels=quantile_levels)
    mean_truep_2cov <- mestimate2cov(x2cov,y,a,ph.true,quantile_levels=quantile_levels)
    #mean_fp_3to2 <-    mestimate2cov(x2cov,y,a,ph.3cov,quantile_levels=quantile_levels)
    
    summary_1_truep_3cov <-qestimate3cov(x,     y, a,prob = ph.true, tau=taus[1], quantile_levels=quantile_levels)
#     summary_1_fp_3cov    <-qestimate3cov(x,     y, a,prob = ph.3cov, tau=taus[1], quantile_levels=quantile_levels)
    summary_1_truep_2cov <-qestimate2cov(x2cov, y, a,prob = ph.true, tau=taus[1], quantile_levels=quantile_levels)
#     summary_1_fp_3to2    <-qestimate2cov(x2cov, y, a,prob = ph.3cov, tau=taus[1], quantile_levels=quantile_levels)
    
    summary_2_truep_3cov <-qestimate3cov(x,     y, a,prob = ph.true, tau=taus[2], quantile_levels=quantile_levels)
#     summary_2_fp_3cov    <-qestimate3cov(x,     y, a,prob = ph.3cov, tau=taus[2], quantile_levels=quantile_levels)
    summary_2_truep_2cov <-qestimate2cov(x2cov, y, a,prob = ph.true, tau=taus[2], quantile_levels=quantile_levels)
#     summary_2_fp_3to2    <-qestimate2cov(x2cov, y, a,prob = ph.3cov, tau=taus[2],  quantile_levels=quantile_levels)
    
    D_mean_truep_3cov <- rbind(D_mean_truep_3cov, mean_truep_3cov)
    D_mean_truep_2cov <- rbind(D_mean_truep_2cov, mean_truep_2cov)
    #D_mean_fp_3cov <- rbind(mean_data_fp_3cov, mean_fp_3cov)
    #D_mean_fp_2cov <- rbind(D_mean_fp_2cov, mean_fp_3to2)    
    D_Q1_truep_3cov <-rbind(D_Q1_truep_3cov, summary_1_truep_3cov )
    D_Q1_truep_2cov <-rbind(D_Q1_truep_2cov, summary_1_truep_2cov)
    # D_Q1_fp_3cov <- rbind(D_Q1_fp_3cov,      summary_1_fp_3cov)
    # D_Q1_fp_2cov<-  rbind(D_Q1_fp_2cov,       summary_1_fp_3to2)
    D_Q2_truep_3cov <-rbind(D_Q2_truep_3cov, summary_2_truep_3cov )
    D_Q2_truep_2cov <-rbind(D_Q2_truep_2cov, summary_2_truep_2cov)
    # D_Q2_fp_3cov <- rbind(D_Q2_fp_3cov,      summary_2_fp_3cov)
    # D_Q2_fp_2cov <- rbind(D_Q2_fp_2cov,      summary_2_fp_3to2)
    
    ptm1=proc.time() - ptm0
    jnk=as.numeric(ptm1[3])
    cat('It took ', jnk, "seconds,  #### Current working dir: ", i,'. \n')
  }
  
  return(list(D_mean_truep_3cov = D_mean_truep_3cov,
              D_mean_truep_2cov = D_mean_truep_2cov, # D_mean_fp_3cov,D_mean_fp_2cov,
              D_Q1_truep_3cov = D_Q1_truep_3cov,
              D_Q1_truep_2cov = D_Q1_truep_2cov,   # D_Q1_fp_3cov,D_Q1_fp_2cov,
              D_Q2_truep_3cov = D_Q2_truep_3cov, 
              D_Q2_truep_2cov = D_Q2_truep_2cov,  # D_Q2_fp_3cov, D_Q2_fp_2cov,
              sim_data=sim_data))
}

sim_function_abso <- function(sim_num=2, nobs ,taus, p_level=0, hard_limit=FALSE){
  D_abso<-NULL
  sim_data<-list() #store original data
  quantile_levels = taus
  #######################  data generation #######################
  for(i in 1:sim_num){
    ptm0 <- proc.time()
    x1 <- runif(nobs)
    x2 <- runif(nobs)
    x2cov=cbind(1,x1,x2)
    Cnobs = combn(1:nobs, 2)
    tp <- exp(-0.5+x1+x2)/(1+exp(-0.5+x1+x2))  
    a<-rbinom(nobs,1,tp)                                                        
    y <- y_function(x1,x2,a)
    sim_data = c(sim_data,list( x = x2cov, y=y, a=a))
    ############################## Propensity score model    ########################
    logit<-glm(a~x1+x2,family=binomial, epsilon=1e-14)
    ph.true<- as.vector(logit$fit)
    
    summary_abso = AbsoDiff(x2cov, y,a,prob=ph.true, quantile_levels=c(0.10,0.30,0.25,0.5), Cnobs)
    D_abso = rbind(D_abso,summary_abso)
    
    ptm1=proc.time() - ptm0
    jnk=as.numeric(ptm1[3])
    cat('It took ', jnk, "seconds,  #### Current working dir: ", i,'. \n')
  }
  
  return(list(D_abso=D_abso,sim_data=sim_data))
}



#######
#######
# load("~/Documents/Academic/Research_DTR/Weighted_2COV/Data_ouput/Exp_2to3_sd05.RData")
# library(xtable)
# Read <- function(R,cores_num){
#   sim_data1<-sim_data2<-sim_data3<-sim_data4<-sim_data5<-sim_data6 <- NULL
#   ColNames3cov = colnames(R[[1]][[1]])
#   ColNames2cov = colnames(R[[1]][[2]])
#   for(i in 1:cores_num){
#     sim_data1<- rbind(sim_data1, R[[i]][[1]])
#     sim_data2<- rbind(sim_data2, R[[i]][[2]])
#     sim_data3<- rbind(sim_data3, R[[i]][[3]])
#     sim_data4<- rbind(sim_data4, R[[i]][[4]])
#     sim_data5<- rbind(sim_data5, R[[i]][[5]])
#     sim_data6<- rbind(sim_data6, R[[i]][[6]])
#   }
#   RR = (list(sim_data1,sim_data2,sim_data3,sim_data4,sim_data5,sim_data6))
# 
#   nobs=10000
#   x1 <- runif(nobs)
#   x2 <- runif(nobs)
#   xp = cbind(1,x1,x2)
#   
#   Cor_matrix = NULL
#   ao25 = rule(x1,x2, m=2.01)
#   for (j in 3:4){
#     summary_corr = c()
#     for (i in 1:dim(RR[[j]])[1]){
#       eta = RR[[j]][i,1:3]
#       tmpaction = as.numeric(I(xp%*%eta>0))
#       tmpcorr = cor(tmpaction,ao25,method="pearson")
#       summary_corr =c(summary_corr, tmpcorr)
#     }
#     Cor_matrix=cbind(Cor_matrix,summary_corr)
#   }
#   
#   ao50 = rule(x1,x2, m=2.69)
#   for (j in 5:6){
#     summary_corr = c()
#     for (i in 1:dim(RR[[j]])[1]){
#       eta = RR[[j]][i,1:3]
#       tmpaction = as.numeric(I(xp%*%eta>0))
#       tmpcorr = cor(tmpaction,ao50,method="pearson")
#       summary_corr =c(summary_corr, tmpcorr)
#     }
#     Cor_matrix=cbind(Cor_matrix,summary_corr)
#   }
#   
#   
# 
#   Part1 = rbind(meansd(RR[[1]]),meansd(RR[[3]]),meansd(RR[[5]]))
#   Part2 = rbind(meansd(RR[[2]]),meansd(RR[[4]]),meansd(RR[[6]]))
#   colnames(Part1) = ColNames3cov
#   colnames(Part2) = ColNames2cov
#   return (list(Table_3cov = Part1, Table_2cov = Part2, Cor_matrix = Cor_matrix))
# }
# 
# meansd <-function(S){
#   S=data.frame(S)
#   means = round(apply(S, 2,mean), digits=3)
#   sds = round(apply(S, 2,sd), digits=3)
#   paste0(means,"(",sds,")")
# }
# rule <- function(x1,x2,m){
#   tmp = 2+ m +(m-4.5)*(x1+x2)-(x1+x2)^2
#   return((tmp > 0)*1)
# }
# 
# 
# (T = Read(results,7))
# 
# library(xtable)
# xtable(T[[1]],caption="Summary of Estimated Rules with 3 Covariates")
# xtable(T[[2]],caption="Summary of Estimated Rules with 2 Covariates")
# xtable(summary(T[[3]]))
# 

# NEED   geom_contour
#result = sim_function_abso(sim_num=7,nobs=500,taus=c(0.1,0.3),p_level=0)
load("~/Documents/Academic/Research_DTR/Weighted_2COV/Data_ouput/Ustat_40.RData")
D_abso = NULL
for (i in 1:7){
  D_abso = rbind(D_abso, Results[[i]]$D_abso)  
}
dim(D_abso)
# load("~/Documents/Academic/Research_DTR/Weighted_2COV/Data_ouput/Ustat_8.RData")
# for (i in 1:7){
#   D_abso = rbind(D_abso, Results[[i]]$D_abso)  
# }
# dim(D_abso)
DR = t(D_abso[,1:3])

mypp <- function (n,DR) {
  x <- seq(0, 1, len=n)
  df <- expand.grid(x1=x, x2=x)
  d3 <- cbind(1,df)
  g = (as.matrix(d3)%*% DR >0)*1
  Rnum <- dim(DR)[2]
  z=(Rnum-apply(g,1,sum))/Rnum
  df$proportion <- z
  df
}
library(ggplot2)
p <- ggplot(mypp(200,DR), aes(x=x1,y=x2))
p + geom_tile(aes(fill=proportion)) + 
    labs(title = paste("Contour plot of proportion of estimated rules \n in agreement with the optimal rule, given covariates \n by 100 simulations, n=1000")) + 
  theme(plot.title = element_text(hjust = 0.5))

p <- ggplot(mypp(100,DR), aes(x=x1,y=x2,z = proportion))
p + stat_contour()

d = (data.frame(t(paste0(round(apply(D_abso,2,mean),2),"(",round(apply(D_abso,2,sd),2),")"))))
colnames(d) = colnames(D_abso)
library(xtable)
xtable(d)
