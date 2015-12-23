#Ustat
library(utils)
library(rgenoud)
sum_names2 <- c("intercept","X1","X2","hatQ")

ModelDec <- function(x1,x2,a){
  error <- rnorm(length(x1), sd = 1)
  y = (1 + 0.5*x1  +0.5*x2+ a*(2*x1-0.4))^2+
    (a*(-1.5)+2)*error
  return(y)
}
Perf_for_eta_2cov <- function(eta,quantile_levels,Cnobs){ 
  n= 10^6
  x1p <- runif(n)
  x2p <- runif(n)
  xp<-cbind(1,x1p,x2p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- ModelDec(x1=x1p, x2=x2p,a=g)
  m=mean(yg)
  qt = quantile(yg, quantile_levels)
  #for (i in 1:length(quantile_levels)){   quantile = c(quantile, quantile(yg,quantile_levels[i]))}
  n=3000
  x1p <- runif(n)
  x2p <- runif(n)
  xp<-cbind(1,x1p,x2p) 
  g<-as.numeric(I(xp%*%eta>0))
  yg <- ModelDec(x1=x1p, x2=x2p,a=g)
  plot(x1p,x2p,col = g+1)
  sum=0
  COMB = data.frame(t(combn(x=1:n,2)))
  COMB$absd = abs(yg[COMB$X1] - yg[COMB$X2])
  S = sum(COMB$absd)
  
  val <- S/(n*(n-1)/2)
  
  return( list(value=c(m,qt,val), Crit = c('mean',quantile_levels,'Exp.AbsoDiff')))
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



AbsoDiff<-function(x,y,a,prob,quantile_levels, Cnobs){
  nvars<- ncol(x)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  est<-genoud(fn=Abso_Est,nvars=nvars,x=x,y=y,a=a,prob=prob, Cnobs= Cnobs,
              print.level=p_level,max=FALSE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              Domains=Domains,starting.values=rep(0,nvars),solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
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
    a <- rbinom(nobs,1,tp)                                                        
    y <- ModelDec(x1,x2,a)
    plot3d(x1,x2,y, col =1+a)
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

R = sim_function_abso(sim_num=1,nobs=1000,taus=c(0.1,0.3,0.5),p_level=0,hard_limit=T)
