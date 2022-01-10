#' Bayesian and Classical Estimates of Topp-leone distribution
#' Bayesian and Classical Estimates of TLL distribution uniform prior distribution for scale parameter and gamma distribution for shape parameter under Linex loss function, Square error loss function and Entropy loss function.
#'
#' @param  x Numeric vector.
#'
#' @param k Linex loss function parameter \code{numeric}
#' @param k1 Entropy loss function \code{numeric}
#' @param p shape parameter of a gamma distribution \code{numeric}
#' @param n scale parameter of a gamma distrbution \code{numeric}
#' @param r number of iteration \code{numeric}
#' @param alpha1 shape \code{numeric}
#' @param theta1 rate \code{numeric}
#' @return Bayesian and ML estimates
#'
#' @examples
#'
#' sim_TLL(n=500, r =1000, alpha=1,theta=2, k=1,p=2,q=1,k1=1)
#'
#' @export
#'
#'
#TLL SIMULATION CODE


#function to do the simulation
sim_TLL<-function(n,r,alpha1,theta1, k,p,q,k1){
  LLS<- c()
  LLA<-c()
  MLS<- c()
  MLA<-c()
  ELS<-c()
  ELA<-c()
  MES<-c()
  MEA<-c()
  SIS<-c()
  MIS <- c()
  SIA <-c()
  MIA <-c()
  alphac<-c()
  thetac<-c()
  alpha<-theta<-c();
  mse <- function(estimate, truth) mean((estimate - truth)^2)
  bias <- function(estimate, truth) mean((estimate - truth))
  lind <- function(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
  {
    lin = u+0.5*((u11*s11)+(u22*s22))+(u1*r1*s11)+(u2*r1*s22)+0.5*((l30*u1*s11^2)+(l03*u2*s22^2))
    return(lin)
  }
  #pdf
  TLL <- function(x,theta,alpha){
    (((2*alpha*theta^2)/(theta+1))*(1+x)*((theta+1+theta*x)*exp(-2*theta*x))/(theta+1))*((1-(((theta+1+theta*x)*exp(-theta*x))/(theta+1))^2)^(alpha-1))
  }

  ##Log likelihood
  pdf<- function(x,par){
    alpha<- par[1]; theta<- par[2];
    LL <-  -sum(log( TLL(x,alpha,theta)))
    LL
  }


  #quantile function
  quantile<-function(w,alpha,theta){
    (-(1 + 1/theta + (1/theta)*lamW::lambertWm1((-(theta + 1))*(sqrt (1-w^(1/alpha)))*exp(-(theta+1)))))
  }
  set.seed(1080)
  for (i in 1:r){
    w<-stats::runif(n = n,min = 0,max = 1)
    x<-quantile(w,alpha1, theta1)
    hat<-try(stats::optim(c(alpha1,theta1),pdf,x=x,control = list(maxit = 1000)),silent=F)
    alpha_mle<-hat$par[1]
    theta_mle<-hat$par[2]
    alpha[i]<- c(alpha,hat$par[1])
    theta[i]<-theta_mle<-c(theta,hat$par[2])
    #obtain standard error
    rho1 <- -1/alpha_mle
    rho2 <- ((p-1)/theta_mle)-q
    t<- (1+theta_mle)
    t5<- (1+theta_mle)^5
    t4<- (1+theta_mle)^4
    t3<- (1+theta_mle)^3
    t2<- (1+theta_mle)^2
    t1 <- (theta_mle+1+theta_mle*x)
    e <- (-theta_mle*x)
    d<-(1-((t1^2*e^2)/t2))
    L02<- -(2*n)/theta_mle +sum(1/t2)+sum(t*((((t1*2)/t3)-(2/t2))/(t1)))-sum(t*(((1/t1)-(t1/t2))/t1^2))+sum((alpha_mle-1)*((((8*t1*e^2)/t3)-((2*e^2)/t2)-((6*t1^2*e^2)/t4)))/d)-sum((alpha_mle-1)*((2*t1^2*e^2)/t3-(2^t1*e^2)/t2)^2/d^2)
    L20<-  -n/alpha_mle^2
    L12 <- sum((((8*t1*e^2)/t3)-((e^2*2)/t2)-((t1^2*e^2)/t4)*6)/d)-sum((((2*t1^2*e^2)/t3)-((2*t1*e^2)/t2))^2/d^2)
    L21 <- 0
    L30 <- (2*n)/alpha_mle^3
    L03 <- (4*n)/theta_mle^3-2*sum(1/t3)+sum(t*((6/t3)-(6*t1)/t4)/t)-2*sum((t*((2*t1)/t3-2/t2))/t1^2)+2*sum(((2*t1)/t3-(6*t1)/t4)/t1)+2*sum(t1*((1/t)-(t1/t2))/t1^3)-2*sum(((1/t)-(t1/t2))/t1^2)+sum((1/d)*((alpha_mle-1)*((12*e^2)/t3-(36*t1*e^2)/t1^2+(24*t1^2*e^2)/t5)))-sum((1/d^2)*(3*(alpha_mle-1)*(((8*t1*e^2)/t3-(2*e^2)/t2)-(6*t1^2*e^2)/t4)*((2*t1^2*e^2)/t3-(2*t1*e^2)/t2)))+sum((2*(alpha_mle-1)*(((2*t1*e^2)/t3-(2*t1*e^2)/t2)^3/d^3)))
    phi22 <- -1/L02
    phi11 <- -1/L20

    # Linex Loss Function
    u = exp(-k*alpha_mle)
    u1 <- -k*exp(-k*alpha_mle)
    u11 <- k^2*exp(-k*alpha_mle)
    u2 = u22 =0
    u= exp(-k*theta_mle)
    u2 <- -k*exp(-k*theta_mle)
    u22 <- k^2*exp(-k*theta_mle)
    u1 = u11 = 0
    LLS[i]<- -(1/k)*log(exp(-k*alpha_mle)+0.5*u11*phi11+rho1*u1*phi11+0.5*L30*u1*phi11^2+0.5*L12*u1*phi11*phi22)
    LLA[i] <- -(1/k)*log(exp(-k*theta_mle)+0.5*u22*phi22+rho2*u2*phi22+0.5*L03*u2*phi22^2+0.5*L21*u2*phi11*phi22)

    #General Entropy Loss Function
    ##For alpha
    u = 1/(alpha_mle^k1)
    u1 = -k1*alpha_mle^(-k1-1)
    u11 = (-k1^2-k1)*alpha_mle^(-k1-2)
    u2=u22=0
    ELS[i]=lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
    u = 1/(theta_mle^k1)
    u2 = -k1*theta_mle^(-k1-1)
    u22 = (-k1^2-k1)*theta_mle^(-k1-2)
    u1=u11=0
    ELA[i] =lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)


    ##SELF
    u = alpha_mle
    u1 = 1
    u11 = 0
    u2=u22=0
    SIS[i]=lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
    u = theta_mle
    u2 = 1
    u22 = 0
    u1=u11=0
    SIA[i] =lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
  }
  #Bais and MSE for MLE
  Bias_mlea<-bias(alpha,alpha1)
  Bias_mleb<-bias(theta,theta1)
  msea_mle<- mse(alpha,alpha1)
  mseb_mle<- mse(theta,theta1)

  #Bais and MSE for LLF
  msela<- mse(LLS,alpha1)
  mselb<-mse(LLA,theta1)
  Bias_lla<-bias(LLS,alpha1)
  Bias_llb<-bias(LLA,theta1)

  #Bais and MSE for GELF
  msega<- mse(ELS,alpha1)
  msegb<-mse(ELA,theta1)
  Bias_ga<-bias(ELS,alpha1)
  Bias_gb<-bias(ELA,theta1)

  #Bais and MSE for SELF
  msesa<- mse(SIS,alpha1)
  msesb<-mse(SIA,theta1)
  Bias_sa<-bias(SIS,alpha1)
  Bias_sb<-bias(SIA,theta1)

  ## Outputs

  matrix1 <- matrix(c(mean(alpha),mean(theta),mean(LLS), mean(LLA),mean(ELS), mean(ELA),mean(SIS), mean(SIA),msea_mle,mseb_mle,msela,mselb,msega,msegb, msesa, msesb ,Bias_mlea, Bias_mleb, Bias_lla, Bias_llb, Bias_ga, Bias_gb, Bias_sa,Bias_sb  ),3, byrow = T)
  dimnames(matrix1) <- list(c("mean", "MSE","Bias"), c("mle-theta", "mle-alpha", "LINEX:bayes-theta", "LINEX:bayes-alpha","GELF:bayes-theta", "GELF:bayes-alpha","SELF:bayes-theta", "SELF:bayes-alpha"))
  return(matrix1)
}



#' Simulate data from TLL distribution
#' @param n  number of sample \code{Numeric}
#'
#' @param alpha which is the shape parameter  \code{numeric}
#' @param theta which is the rate  rate\code{numeric}
#'
#'
#' @examples
#' x = Rtll(n=50, alpha=2,theta=1)
#'
#' @export
#'
# quantile function
Rtll<-function(n,alpha,theta){
  w<-stats::runif(n ,min = 0,max = 1)
return((-(1 + 1/theta + (1/theta)*lamW::lambertWm1((-(theta + 1))*(sqrt (1-w^(1/alpha)))*exp(-(theta+1))))))
}


#' Bayesian and Classical Estimates of Topp-leone distribution
#' Bayesian and Classical Estimates of TLL distribution uniform prior distribution for scale parameter and gamma distribution for shape parameter under Linex loss function, Square error loss function and Entropy loss function.
#'
#' @param  x Numeric vector.
#'
#' @param k Linex loss function parameter \code{numeric}
#' @param k1 Entropy loss function \code{numeric}
#' @param p shape parameter of a gamma distribution \code{numeric}
#'#'@param alpha1 shape \code{numeric}
#' @param theta1 rate \code{numeric}
#' @return Bayesian and ML estimates on real dataset
#'
#' @examples
#'x = Rtll(n=50, alpha=1,theta=1)
#'real_TLL(x, k=1,p=1,q=1,k1=1)
#' @export
#'
#function to do the simulation
real_TLL<-function(x, k,p,q,k1){
  n = length(x)
  lind <- function(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
  {
    lin = u+0.5*((u11*s11)+(u22*s22))+(u1*r1*s11)+(u2*r1*s22)+0.5*((l30*u1*s11^2)+(l03*u2*s22^2))
    return(lin)
  }
  #pdf
  TLL <- function(x,theta,alpha){
    (((2*alpha*theta^2)/(theta+1))*(1+x)*((theta+1+theta*x)*exp(-2*theta*x))/(theta+1))*((1-(((theta+1+theta*x)*exp(-theta*x))/(theta+1))^2)^(alpha-1))
  }

  ##Log likelihood
  pdf<- function(x,par){
    alpha<- par[1]; theta<- par[2];
    LL <-  -sum(log( TLL(x,alpha,theta)))
    LL
  }
  hat<-try(stats::optim(c(0.5,0.5),pdf,x=x,control = list(maxit = 1000)),silent=F)
  alpha_mle<-hat$par[1]
  theta_mle<-hat$par[2]
  #obtain standard error
  rho1 <- -1/alpha_mle
  rho2 <- ((p-1)/theta_mle)-q
  t<- (1+theta_mle)
  t5<- (1+theta_mle)^5
  t4<- (1+theta_mle)^4
  t3<- (1+theta_mle)^3
  t2<- (1+theta_mle)^2
  t1 <- (theta_mle+1+theta_mle*x)
  e <- (-theta_mle*x)
  d<-(1-((t1^2*e^2)/t2))
  L02<- -(2*n)/theta_mle +sum(1/t2)+sum(t*((((t1*2)/t3)-(2/t2))/(t1)))-sum(t*(((1/t1)-(t1/t2))/t1^2))+sum((alpha_mle-1)*((((8*t1*e^2)/t3)-((2*e^2)/t2)-((6*t1^2*e^2)/t4)))/d)-sum((alpha_mle-1)*((2*t1^2*e^2)/t3-(2^t1*e^2)/t2)^2/d^2)
  L20<-  -n/alpha_mle^2
  L12 <- sum((((8*t1*e^2)/t3)-((e^2*2)/t2)-((t1^2*e^2)/t4)*6)/d)-sum((((2*t1^2*e^2)/t3)-((2*t1*e^2)/t2))^2/d^2)
  L21 <- 0
  L30 <- (2*n)/alpha_mle^3
  L03 <- (4*n)/theta_mle^3-2*sum(1/t3)+sum(t*((6/t3)-(6*t1)/t4)/t)-2*sum((t*((2*t1)/t3-2/t2))/t1^2)+2*sum(((2*t1)/t3-(6*t1)/t4)/t1)+2*sum(t1*((1/t)-(t1/t2))/t1^3)-2*sum(((1/t)-(t1/t2))/t1^2)+sum((1/d)*((alpha_mle-1)*((12*e^2)/t3-(36*t1*e^2)/t1^2+(24*t1^2*e^2)/t5)))-sum((1/d^2)*(3*(alpha_mle-1)*(((8*t1*e^2)/t3-(2*e^2)/t2)-(6*t1^2*e^2)/t4)*((2*t1^2*e^2)/t3-(2*t1*e^2)/t2)))+sum((2*(alpha_mle-1)*(((2*t1*e^2)/t3-(2*t1*e^2)/t2)^3/d^3)))
  phi22 <- -1/L02
  phi11 <- -1/L20

  # Linex Loss Function
  u = exp(-k*alpha_mle)
  u1 <- -k*exp(-k*alpha_mle)
  u11 <- k^2*exp(-k*alpha_mle)
  u2 = u22 =0
  u= exp(-k*theta_mle)
  u2 <- -k*exp(-k*theta_mle)
  u22 <- k^2*exp(-k*theta_mle)
  u1 = u11 = 0
  LLS<- -(1/k)*log(exp(-k*alpha_mle)+0.5*u11*phi11+rho1*u1*phi11+0.5*L30*u1*phi11^2+0.5*L12*u1*phi11*phi22)
  LLA <- -(1/k)*log(exp(-k*theta_mle)+0.5*u22*phi22+rho2*u2*phi22+0.5*L03*u2*phi22^2+0.5*L21*u2*phi11*phi22)

  #General Entropy Loss Function
  ##For alpha
  u = 1/(alpha_mle^k1)
  u1 = -k1*alpha_mle^(-k1-1)
  u11 = (-k1^2-k1)*alpha_mle^(-k1-2)
  u2=u22=0
  ELS=lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
  u = 1/(theta_mle^k1)
  u2 = -k1*theta_mle^(-k1-1)
  u22 = (-k1^2-k1)*theta_mle^(-k1-2)
  u1=u11=0
  ELA =lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)


  ##SELF
  u = alpha_mle
  u1 = 1
  u11 = 0
  u2=u22=0
  SIS=lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
  u = theta_mle
  u2 = 1
  u22 = 0
  u1=u11=0
  SIA =lind(u,u1,u2,u11,u22,s11=phi11,s22=phi22,r1=rho1,r2=rho2,l03=L03,l30=L30)
  ## Outputs
  matrix1 <- matrix(c(mean(alpha_mle),mean(theta_mle),mean(LLS), mean(LLA),mean(ELS), mean(ELA),mean(SIS), mean(SIA)),1, byrow = T)
  dimnames(matrix1) <- list(c("mean"), c("mle-theta", "mle-alpha", "LINEX:bayes-theta", "LINEX:bayes-alpha","GELF:bayes-theta", "GELF:bayes-alpha","SELF:bayes-theta", "SELF:bayes-alpha"))
  return(matrix1)
}


