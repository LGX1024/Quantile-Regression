#######################   BCQR   #####################
#######################  By LGX  #####################
######################################################
  library(mvtnorm)
  library(lpSolve)
  library(splines)
  
######################################################
################## input data ########################
######################################################
  q = 3; p = 2; n = 400
  beta = c(2.5,1.2,0.5)
  a1 <- function(x){
    5.5+10*exp(2*x-1)
    #x^2
  }
  a2 <- function(x){
    2-x*(2-x)
    #x+1
  }
  z_var = diag(1.5,q); x_var = diag(1,p)
  set.seed(1)
  z.train <- rmvnorm(n,rep(0,q),z_var)
  x.train <- rmvnorm(n,rep(0,p),x_var)
  u = runif(n)
  esp.train = rnorm(n,0,1)
  a <- cbind(a1(u),a2(u))
  y.train <- z.train%*%beta+diag(x.train%*%t(a))+esp.train 
#######################################################
##################      function       ################
################## quantile regression ################
#######################################################
  QR.LGX<-function(Covariate, Response, Tau){
    X <- Covariate; Y <- Response; tau <- Tau
    n = nrow(Y); p = ncol(X); 
    one.vec = rep(1,n); I.mat = diag(1,n)
    
    #objective
    obj <- c(rep(0,2),rep(0,2*p),tau*one.vec,(1-tau)*one.vec)
    
    #const.mat
    const_mat1 <- cbind(one.vec,-one.vec,X,-X,I.mat,-I.mat)
    const_mat  <- rbind(const_mat1,diag(1,2*(1+p+n)))
    
    #const.dir
    const_dir <- c(rep("=",n),rep(">=",2*(1+p+n)))
    #const.rhs
    const_rhs <- c(c(Y),rep(0,2*(1+p+n)))
    
    solution <- lp("min",obj,const_mat,const_dir,const_rhs)$solution
    
    b0.hat <- solution[1]-solution[2]
    sol <- solution[-c(1:2)]
    b.hat  <- sol[1:p]-sol[(p+1):(2*p)]
    resid <- sol[(2*p+1):(2*p+n)]-sol[(2*p+n+1):(2*p+2*n)]
                    
    fitvalue <- b0.hat+X%*%b.hat
    result <- list(b0.hat=b0.hat,b.hat=b.hat,resid=resid,fitted.values=fitvalue,Tau=tau)
    return(result)
  }
########################################################
############ composite quantile regression #############
########################################################
  CQR.LGX<-function(Covariate, Response, M){
    X <- Covariate; Y <- Response; m <- M
    n = nrow(Y); p = ncol(X); 
    one.vec = rep(1,n); I.mat = diag(1,n*m)
    tau <- c(1:m)/(m+1)
    #objective
    
    obj <- c(rep(0,2*m),rep(0,2*p),tau%x%one.vec,(1-tau)%x%one.vec)
    
    #const.mat
    const_mat1 <- cbind(diag(1,m)%x%one.vec,-diag(1,m)%x%one.vec,
                        rep(1,m)%x%X,-rep(1,m)%x%X,I.mat,-I.mat)
    const_mat  <- rbind(const_mat1,diag(1,2*(m+p+n*m)))
    
    #const.dir
    const_dir <- c(rep("=",n*m),rep(">=",2*(m+p+n*m)))
    #const.rhs
    const_rhs <- c(rep(c(Y),m),rep(0,2*(m+p+n*m)))
    
    solution <- lp("min",obj,const_mat,const_dir,const_rhs)$solution
    
    b0.hat <- solution[1:m]-solution[(m+1):(2*m)]
    sol <- solution[-c(1:(2*m))]
    b.hat  <- sol[1:p]-sol[(p+1):(2*p)]
    #resid <- matrix(sol[(2*p+1):(2*p+n*m)]-sol[(2*p+n*m+1):(2*p+2*n*m)],byrow = F,nrow = n)
    
    #fitvalue <- b0.hat+X%*%b.hat
    result <- list(b0.hat=b0.hat,b.hat=b.hat,Tau=tau)
    return(result)
  }
##############################################################
############ composite quantile regression by B-splines ######
##############################################################
  BCQR<-function(Z.Covariate, X.Covariate, Response, M, U, kn){
    Z <- Z.Covariate; X <- X.Covariate; 
    Y <- Response; m <- M
    
    n = nrow(Y); q = ncol(Z); p = ncol(X)
    one.vec = rep(1,n); I.mat = diag(1,(n*m)); 
    
    # B-spline basic function 
    B = bs(U,df=kn)
    N = ncol(B)
    W = matrix(0,nrow = n,ncol = N*p)
    for (i in 1:n) {
      W[i,] = X[i,]%x%B[i,]
    }
    tau<-c(1:m)/(m+1)
    n.all = 2*(m+q+N*p+n*m)
    #  objective in
    obj_mat <- c(rep(0,2*m),rep(0,2*q),rep(0,2*N*p),tau%x%one.vec,(1-tau)%x%one.vec)
    # constraint condition
    const_mat1 <- cbind(diag(1,m)%x%one.vec,-diag(1,m)%x%one.vec,
                        rep(1,m)%x%Z,-rep(1,m)%x%Z,rep(1,m)%x%W,-rep(1,m)%x%W,
                        I.mat,-I.mat)
    const_mat  <- rbind(const_mat1,diag(1,n.all))
    
    const_dir <- c(rep("=",n*m),rep(">=",n.all))
    const_rhs <- c(rep(c(Y),m),rep(0,n.all))
    
    solution <- lp("min",obj_mat,const_mat,const_dir,const_rhs)$solution
    # b0 hat
    b0.hat <- solution[1:m]-solution[(m+1):(2*m)]
    # beta hat
    sol_beta <- solution[-c(1:(2*m))]
    b.hat  <- sol_beta[1:q]-sol_beta[(q+1):(2*q)]
    # alpha hat
    sol_alpha <- sol_beta[-c(1:(2*q))]
    alpha.hat <- sol_alpha[1:(N*p)]-sol_alpha[(N*p+1):(2*N*p)]
    alpha.hat.mat <- matrix(alpha.hat ,byrow = F,nrow = N)
    # a(u) hat 
    a.hat <- B%*%alpha.hat.mat
    # resid
    sol_resid <- sol_alpha[-c(1:(2*N*p))]
    resid <- matrix(sol_resid[1:(m*n)]-sol_resid[(m*n+1):(2*m*n)],
                    byrow = F,nrow = n)
    # loss
    loss <- matrix(c(rep(tau,each=n)  *sol_alpha[(2*N*p+1):(2*N*p+n*m)]+
                     rep(1-tau,each=n)*sol_alpha[(2*N*p+n*m+1):(2*N*p+2*m*n)]),
                   byrow = F,nrow = n)
    # fitted value
    fitvalue <- Z%*%b.hat + W%*%alpha.hat
    #fitvalue <- Z%*%b.hat + diag(X%*%t(a.hat))
    
    result <- list(b0.hat=b0.hat,b.hat=b.hat,a.hat=a.hat,
                   resid=resid,fitted.values=fitvalue,LOSS=loss,Tau=tau,Basic=B)
    
    return(result)
  }
  
  ###########################################
  time1<-proc.time()
  solution<-BCQR(Z.Covariate = z.train, X.Covariate = x.train,
                 Response = y.train,M=3,U=u,kn=10)
  time2<-proc.time()
  time2-time1
  bhat<-solution$b.hat
  bhat
  ahat<-solution$a.hat
  par(mfrow=c(2,2))
  plot(u,a1(u),col="red",ylim = c(0,32)); par(new=T);plot(u,ahat[,1],ylim = c(0,32))
  plot(u,a2(u),col="red",ylim = c(0, 2)); par(new=T);plot(u,ahat[,2],ylim = c(0, 2))
  sqrt(sum((a1(u)-ahat[,1])^2)/n);sqrt(sum((a2(u)-ahat[,2])^2)/n)
  
  library(cqrReg)
  B <- solution$Basic
  N <- ncol(B)
  W = matrix(0,nrow = n,ncol = N*p)
  for (i in 1:n) {
    W[i,] = x.train[i,]%x%B[i,]
  }
  tau<-solution$Tau
  kk<-cqr.mm(cbind(z.train,W),y.train,tau = tau)
  kk$b; kk$beta[1:q]
  solution$b0.hat; solution$b.hat
  alpha<-matrix(kk$beta[-c(1:q)],nrow = N)
  ahat_cqr<-B%*%alpha
  
  plot(u,a1(u),col="red",ylim = c(0,32)); par(new=T);plot(u,ahat_cqr[,1],ylim = c(0,32))
  plot(u,a2(u),col="red",ylim = c(0,2));par(new=T);plot(u,ahat_cqr[,2],ylim = c(0,2))
  sqrt(sum((a1(u)-ahat_cqr[,1])^2)/n);sqrt(sum((a2(u)-ahat_cqr[,2])^2)/n)
  
  TT<-cbind(u,ahat,a1(u),a2(u))
  tt<-TT[-which(TT[,1]<=0.1),]
  plot(tt[,1],tt[,4],col="red",ylim = c(0,32)); par(new=T);plot(tt[,1],tt[,2],ylim = c(0,32))
  plot(tt[,1],tt[,5],col="red",ylim = c(0, 2)); par(new=T);plot(tt[,1],tt[,3],ylim = c(0, 2))
  sqrt(sum((tt[,4]-tt[,2])^2)/n);sqrt(sum((tt[,5]-tt[,3])^2)/n)
  
  
  
  
#################################################
############### repeat 100 times ################
#################################################
  tt=0
  amse<-matrix(0,ncol = 2,nrow = 100)
  time1<-proc.time()
  repeat{
    tt=tt+1
    z <- rmvnorm(n,rep(0,q),z_var)
    x <- rmvnorm(n,rep(0,p),x_var)
    u = runif(n)
    esp = rnorm(n,0,1)
    a <- cbind(a1(u),a2(u))
    y <- z%*%beta+diag(x%*%t(a))+esp
    bcqr<-BCQR(z,x,y,9,u,10)
    amse[tt,1]<-sqrt(sum((a1(u)-bcqr$a.hat[,1])^2)/n)
    amse[tt,2]<-sqrt(sum((a2(u)-bcqr$a.hat[,2])^2)/n)
    rm(x,y,z,u,esp,a,bcqr)
    if(tt==100)break()
  }
  time2<-proc.time()
  time2-time1
  apply(amse, 2, mean)
   