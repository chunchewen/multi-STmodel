#################################
# Simulation 07-30-21.r         #
# Bivariate NB Model for Sim 1  #
# Correlated Beta and Bs        #
# Georgia Adj Matrix            #
#################################
# Uses approximate conditional prior from the Corollary
# Infile:   Download Georgia County Adj Matrix (Adj_Ga_Ordered.txt)
# Outfile:  Save MCMC Samples: Ex:Beta.txt, Beta2.txt, R1.txt, etc 
# July 30, 2021                 #
#################################
##########################
# Load Required Packages #
##########################
library(MASS)         # For glm.nb function
library(BayesLogit)   # For rpg function
library(mvtnorm)      # For rmnvorm function
library(splines2)     # For splines
library(spam)         # For sparse matrices and rmvnorm.canonical function
library(MCMCpack)     # For RIwish function
library(tidyverse)    # For plotting
library(plotly)       # Plotting
library(coda)         # Diagnostics

out<-"...\\Directory\\"  # Store MCMC Samples

################
# Spatial Data #
################
set.seed(71454)		

A<-matrix(scan("...\\Adj_Ga_Ordered.txt"),159,159)    # Import Georgia County Adj Matrix
m<-apply(A,1,sum)	                                    # No. neighbors
ncounty<-ncol(A)
Q<-as.spam(diag(m))-as.spam(A)+diag(0.0001,ncounty)   # Ridge by 0.0001 to invert and generate data (negligible impact on results)

#################
# Temporal Data	#
#################         
lday<-300                 # Number of time points (lday = "last day")
nis<-rep(lday,ncounty)    # Number of obs for each county
id<-rep(1:ncounty,nis)    # County id indicator for each day
n<-length(id)             # Total sample size (N in paper)
days<-rep(1:lday,ncounty) # Repeat 1:lday for each ncounty 
pop<-rep(sample(10000:1000000,ncounty,replace=T),nis)  # Generate county population sizes
lpop<-log(pop)            # Offset

#############################
# Fixed-effect Spline Coefs #
# for Overall Time Trend    #
#############################
knots<-seq(10,290,by=10)                              # 29 interior knots
Xtmp<-bSpline(1:lday,knots=knots,intercept=T)         # Xstar in paper
k<-ncol(Xtmp)                                         # Total dim of X*
X<-apply(Xtmp,2,rep,ncounty)

D1 <-diag(1,nrow=k)[-k,] + cbind(0,diag(-1,nrow=k-1)) # First order difference matrix
K  <- as.spam(t(D1)%*%D1) + diag(0.0001,k)            # Ridge by 0.0001 to invert and generate data
Sigmabeta<-matrix(c(.50,.25,.25,.20),2,2)             # Cov b.w beta1 and beta2: rho ~ 0.79
precbeta<-(solve(Sigmabeta)%x%K)                      # Precision
beta<-rmvnorm.prec(1, Q=precbeta)[1,]                 # rmvnorm.prec assumes mean 0

Beta<-matrix(beta,k,2)
beta1<-truebeta1<-Beta[,1]-mean(Beta[,1])             # Center for identifiability                       
beta2<-truebeta2<-Beta[,2]-mean(Beta[,2])

plot(1:lday,exp(Xtmp%*%truebeta1),type="l")
plot(1:lday,exp(Xtmp%*%truebeta2),type="l")

###########################
# County-Level Covariates #
###########################
w1<-rep(rnorm(ncounty,0,1),nis)
w2<-rep(rbinom(ncounty,1,.75),nis)
W<-cbind(w1,w2)
alpha1<-truealpha1<-c(0.25,-0.25)
alpha2<-truealpha2<-c(-0.25,.25)
p<-ncol(W)

################################
# County-Specific Spline Coefs #
################################
Sigmab<-matrix(c(.20,.05,.05,.20),2,2)               # Cov b/w b1 and b2 across space and time: rho=0.25
precb<-(solve(Sigmab)%x%K%x%Q)                       # Precision
b<-rmvnorm.prec(1, Q=precb)[1,]                      # rmvnorm.prec assumes mean 0 (this takes a few minutes)
# write(b,paste(outb,"true_bs.txt",sep=""))          # Write random effect to file to avoid re-generating
#b<-scan(paste(outb,"true_bs.txt",sep=""))
Bmat<-matrix(b,ncounty*k,2)
b1<-Bmat[,1]
b2<-Bmat[,2]
Bs1<-matrix(b1,ncounty,k)
Bs2<-matrix(b2,ncounty,k)

Crow<-diag(k) - tcrossprod(rep(1, k))/k                   # Row centering matrix -- like A matrix in SST
Ccol<-diag(ncounty) - tcrossprod(rep(1, ncounty))/ncounty # Column centering matrix 
B1<-trueB1<-Ccol%*%Bs1%*%Crow                             # Center by row and column 
B2<-trueB2<-Ccol%*%Bs2%*%Crow
cov(cbind(c(trueB1),c(trueB2)))
tmp1<-c(Xtmp%*%t(B1))
tmp2<-c(Xtmp%*%t(B2))

######################################
# Generate Responses Y1 and Y2 Using #
# Pillow and Scott (2012) format     #
######################################
eta1<-W%*%alpha1+X%*%beta1+tmp1+lpop-log(100000)   # Incident case rate per 100k 
eta2<-W%*%alpha2+X%*%beta2+tmp2+lpop-log(1000000)  # Death rate per Million
r1<-1                         # Dispersion 
r2<-1.25

psi1<-exp(eta1)/(1+exp(eta1)) # Prob of success y1 (this uses Pillow and Scott (2012) format)
q1<-1-psi1                    # Note: rnbinom models q=1-p 
y1<-rnbinom(n,r1,q1)          # Y1: Incident case counts

psi2<-exp(eta2)/(1+exp(eta2)) # Prob of success for y2
q2<-1-psi2                    # Note: rnbinom models q=1-p 
y2<-rnbinom(n,r2,q2)          # Y2: Death counts

incidence1<-y1/pop*100000
incidence2<-y2/pop*1000000

orr1<-tapply(incidence1,days,mean)  # Daily statewide mean incidence rates across counties (for pop. avg comparisons)
orr2<-tapply(incidence2,days,mean)  # Daily statewide mean death rates across counties 
plot(1:lday,orr1,main="Observed Population Mean Trend, Y1")
plot(1:lday,orr2,main="Observed Population Mean Trend, Y2")

############################
# Plot Histogram of Counts #
############################
tmp<-table(y2)/n*100
barplot(tmp, ylab="Percent",xlab="Count",col="slateblue4")

################################
# MLEs from Fixed Effect Model #
################################
nb1<-glm.nb(y1~ W+ X-1+offset(lpop-log(100000)))
nb2<-glm.nb(y2~ W+ X-1+offset(lpop-log(1000000)))

###########
# Priors  #
###########
alpha0<-rep(0,p)              # Prior mean of alpha1 and alpha2
T0<-diag(0.01,p)              # Prior precision of alpha1 and alpha2

##################
# Initial Values #
##################
alpha1<-coef(nb1)[1:p]
beta1<-coef(nb1)[(p+1):(k+p)]
alpha2<-coef(nb2)[1:p]
beta2<-coef(nb2)[(p+1):(k+p)]
taub1<-taub2<-1			          # Random effect precision
r1<-r2<-1
l1<-l2<-rep(0,n)              # latent counts for updating r (see Zhou et al)
B1<-B2<-matrix(0,ncounty,k)   # Spline coefs
Q<-diag(m)-A
K<-as.spam(t(D1)%*%D1)
ssq1<-ssq2<-ssqb1<-ssqb2<-1   # Var beta1 and beta1, B1 and B2
rho<-rhob<-0                  # Corr(beta1,beta2) and corr(b1,b2)

############
# Num Sims #
############
nsim<-10500 #10500  			        # Number of MCMC iterations -- MCMC stabilizes quickly, so shorter runs poss. for prelim results
thin<-1                   # Thinning interval
burn<-500	                # Burn-in
lastit<-(nsim-burn)/thin	# Last stored value

########
# MCMC #
########
tmp.time<-proc.time()

for (i in 1:nsim){
  # Update latent counts, l1 and l2, using Chinese restaurant table distribution (c.f. Zhou and Carin, Dadaneh)
  for(j in 1:n) l1[j]<-sum(rbinom(y1[j],1,round(r1/(r1+1:y1[j]-1),6))) # Rounding avoids numerical stability
  for(j in 1:n) l2[j]<-sum(rbinom(y2[j],1,round(r2/(r2+1:y2[j]-1),6))) 
  
  # Update r1 and r2 from conjugate gamma distribution given l and psi
  tmp1<-c(Xtmp%*%t(B1))
  eta1<-W%*%alpha1+X%*%beta1+tmp1+lpop-log(100000)
  psi1<-exp(eta1)/(1+exp(eta1))
  r1<-rgamma(1,.01+sum(l1),0.01-sum(log(1-psi1)))
  
  tmp2<-c(Xtmp%*%t(B2))
  eta2<-W%*%alpha2+X%*%beta2+tmp2+lpop-log(1000000)
  psi2<-exp(eta2)/(1+exp(eta2))
  r2<-rgamma(1,.01+sum(l2),0.01-sum(log(1-psi2)))
  
  # Update z1 and z2
  eta1<-W%*%alpha1+X%*%beta1+tmp1+lpop-log(100000)
  w1<-rpg(n,y1+r1,eta1)                   # Polya weights
  z1<-(y1-r1)/(2*w1)-lpop+log(100000)     # Latent response, z (incorporating offsets)
  
  eta2<-W%*%alpha2+X%*%beta2+tmp2+lpop-log(1000000)
  w2<-rpg(n,y2+r2,eta2)                   
  z2<-(y2-r2)/(2*w2)-lpop+log(1000000)       
  
  # Update alpha1 and alpha2
  v<-solve(T0+crossprod(sqrt(w1)*W))      # posterior variance
  m<-v%*%(T0%*%alpha0+t(w1*W)%*%(z1-X%*%beta1-tmp1))
  alpha1<-c(rmvnorm(1,m,v))
  
  v<-solve(T0+crossprod(sqrt(w2)*W))
  m<-v%*%(T0%*%alpha0+t(w2*W)%*%(z2-X%*%beta2-tmp2))
  alpha2<-c(rmvnorm(1,m,v))
  
  # Update beta1 and beta2 using bivariate RW(1)
  # Beta 1
  priorprec<-1/((1-rho^2)*ssq1)*K                 # Conditional prior prec
  priormean<-rho*sqrt(ssq1/ssq2)*beta2            # Cond prior mean
  prec<-priorprec+crossprod(X*sqrt(w1))           # Post. Precision
  m<-c(priorprec%*%priormean)+t(w1*X)%*%(z1-W%*%alpha1-tmp1) # Post mean without multiplying by prec 
  beta1<-rmvnorm.canonical(1, m, prec)[1,]        # Sample beta1
  beta1<-beta1-mean(beta1)                        # Center
  
  # Beta 2
  priorprec<-1/((1-rho^2)*ssq2)*K 
  priormean<-rho*sqrt(ssq2/ssq1)*beta1
  prec<-priorprec+crossprod(X*sqrt(w2))
  m<-c(priorprec%*%priormean)+t(w2*X)%*%(z2-W%*%alpha2-tmp2)  
  beta2<-rmvnorm.canonical(1, m, prec)[1,]
  beta2<-beta2-mean(beta2)                        # Center
  
  # Update Sigma_beta
  Betamat<-cbind(beta1,beta2)
  nu<-4 + k-1                                       # post df 
  S<-diag(2)+t(Betamat)%*%K%*%Betamat
  Sigma<-riwish(nu,S)                               # Sigma_beta
  ssq1<-Sigma[1,1]
  ssq2<-Sigma[2,2]
  rho<-Sigma[1,2]/sqrt(ssq1*ssq2)
  
  # Update bs using Bivariate RWCAR Prior (uses approximation in corollary)
  # b_11 and b_21
  priorprec<-2/((1-rhob^2)*ssqb1)*Q                           # Prior Prec -- Note: multiply by 2
  priormean<-B1[,2]/2+rhob*sqrt(ssqb1/ssqb2)*(B2[,1]-B2[,2]/2)# Prior Mean w/o premultiplying by var
  prec<-priorprec+diag(tapply(X[,1]^2*w1,id,sum))             # Post Precision     
  tmp1<-c(Xtmp[,-1]%*%t(B1[,-1]))
  m<-c(priorprec%*%priormean)+tapply(w1*X[,1]*(z1-W%*%alpha1-X%*%beta1-tmp1),id,sum)    # Post mean
  B1[,1]<-rmvnorm.canonical(ncounty,m,prec)[1,]
  
  priorprec<-2/((1-rhob^2)*ssqb2)*Q                           # Prior Prec -- Note: multiply by 2
  priormean<-B2[,2]/2+rhob*sqrt(ssqb2/ssqb1)*(B1[,1]-B1[,2]/2)# Prior Mean w/o premultipliying by var
  prec<-priorprec+diag(tapply(X[,1]^2*w2,id,sum))             # Post Precision
  tmp2<-c(Xtmp[,-1]%*%t(B2[,-1]))                             # Or, could write Z as ncounty x lday  matrix and condition on Z[,1]
  m<-c(priorprec%*%priormean)+tapply(w2*X[,1]*(z2-W%*%alpha2-X%*%beta2-tmp2),id,sum)    # Post mean
  B2[,1]<-rmvnorm.canonical(ncounty,m,prec)[1,]
  
  # b_12, ..., b_1{k-1} and b_21,...b_2{k-1}
  for (j in 2:(k-1)){
    priorprec<-2/((1-rhob^2)*ssqb1)*Q                         # Prior Prec -- Note: multiply by 2
    priormean<-(B1[,j-1]+B1[,j+1])/2+rhob*sqrt(ssqb1/ssqb2)*(B2[,j]-(B2[,j-1]+B2[,j+1])/2)
    prec<-priorprec+diag(tapply(X[,j]^2*w1,id,sum))             
    tmp1<-c(Xtmp[,-j]%*%t(B1[,-j]))
    m<-c(priorprec%*%priormean)+tapply(w1*X[,j]*(z1-W%*%alpha1-X%*%beta1-tmp1),id,sum)
    B1[,j]<-rmvnorm.canonical(ncounty,m,prec)[1,]
    
    priorprec<-2/((1-rhob^2)*ssqb2)*Q                         # Prior Prec -- Note: multiply by 2
    priormean<-(B2[,j-1]+B2[,j+1])/2+rhob*sqrt(ssqb2/ssqb1)*(B1[,j]-(B1[,j-1]+B1[,j+1])/2)
    prec<-priorprec+diag(tapply(X[,j]^2*w2,id,sum),ncounty)            
    tmp2<-c(Xtmp[,-j]%*%t(B2[,-j]))
    m<-c(priorprec%*%priormean)+tapply(w2*X[,j]*(z2-W%*%alpha2-X%*%beta2-tmp2),id,sum)
    B2[,j]<-rmvnorm.canonical(ncounty,m,prec)[1,]
  }
  
  # b_1k and b_2k
  priorprec<-1/((1-rhob^2)*ssqb1)*Q                           # Prior Prec 
  priormean<-B1[,k-1]+rhob*sqrt(ssqb1/ssqb2)*(B2[,k]-B2[,k-1])# Prior Mean
  prec<-priorprec+diag(tapply(X[,k]^2*w1,id,sum))             # Post Precision
  tmp1<-c(Xtmp[,-k]%*%t(B1[,-k]))                             # Or, could write Z as ncounty x lday  matrix and condition on Z[,1]
  m<-c(priorprec%*%priormean)+tapply(w1*X[,k]*(z1-W%*%alpha1-X%*%beta1-tmp1),id,sum)
  B1[,k]<-rmvnorm.canonical(ncounty,m,prec)[1,]
  
  priorprec<-1/((1-rhob^2)*ssqb2)*Q                           # Prior Prec 
  priormean<-B2[,k-1]+rhob*sqrt(ssqb2/ssqb1)*(B1[,k]-B1[,k-1])# Prior Mean
  prec<-priorprec+diag(tapply(X[,k]^2*w2,id,sum))             # Post Precision
  tmp2<-c(Xtmp[,-k]%*%t(B2[,-k]))                             # Or, could write Z as ncounty x lday  matrix and condition on Z[,1]
  m<-c(priorprec%*%priormean)+tapply(w2*X[,k]*(z2-W%*%alpha2-X%*%beta2-tmp2),id,sum)
  B2[,k]<-rmvnorm.canonical(ncounty,m,prec)[1,]
  
  # Center Bs
  B1<-Ccol%*%B1%*%Crow
  B2<-Ccol%*%B2%*%Crow
  
  # Update Sigmab
  Bmat<-cbind(c(B1),c(B2))            # ncounty*k x 2 
  nu<-4+(ncounty-1)*(k-1)             # Posterior df (prior is 4)
  S<-diag(2)+t(Bmat)%*%(K%x%Q)%*%Bmat # Posterior scale
  Sigmab<-riwish(nu,S)
  ssqb1<-Sigmab[1,1]
  ssqb2<-Sigmab[2,2]
  rhob<-Sigmab[1,2]/sqrt(ssqb1*ssqb2)
  
  # Likelihood for WAIC
  tmp1<-c(Xtmp%*%t(B1))
  tmp2<-c(Xtmp%*%t(B2))
  eta1<-W%*%alpha1+X%*%beta1+tmp1+lpop-log(100000)   # Incident case rate per 100k 
  eta2<-W%*%alpha2+X%*%beta2+tmp2+lpop-log(1000000)  # Death rate per Million
  
  psi1<-exp(eta1)/(1+exp(eta1)) # Prob of success y1 (this uses Pillow and Scott (2012) format)
  psi2<-exp(eta2)/(1+exp(eta2)) # Prob of success for y2
  q1<-1-psi1                    # Note: rnbinom models q=1-p 
  q2<-1-psi2                    # Note: rnbinom models q=1-p 
  
  # Joint likelihood (for WAIC)
  L<-dnbinom(y1,r1,q1)*dnbinom(y2,r2,q2)
  
  # Store  
  if (i> burn & i%%thin==0) {
    write(alpha1, paste(out,"Alpha1.txt",sep=""),ncol=p,append=T)
    write(alpha2, paste(out,"Alpha2.txt",sep=""),ncol=p,append=T) 
    write(beta1, paste(out,"Beta1.txt",sep=""),ncol=k,append=T)
    write(beta2, paste(out,"Beta2.txt",sep=""),ncol=k,append=T)      
    write(r1, paste(out,"R1.txt",sep=""),ncol=1,append=T)     
    write(r2, paste(out,"R2.txt",sep=""),ncol=1,append=T) 
    write(c(Sigma),paste(out,"Sigma.txt",sep=""),ncol=4,append=T)
    write(c(Sigmab),paste(out,"Sigmab.txt",sep=""),ncol=4,append=T) 
    write(c(t(B1)),paste(out,"B1.txt",sep=""),ncol=ncounty*k,append=T)
    write(c(t(B2)),paste(out,"B2.txt",sep=""),ncol=ncounty*k,append=T)
    write(L,paste(out,"L.txt",sep=""),ncol=n,append=T)
  }
  
  if (i%%10==0) {
    print(i)    
    print(alpha1)
    print(alpha2)
    print(r1)
    print(r2)
    print(Sigma)
    print(Sigmab)
    print(L[1:5])
  }
  
}

run.time<-proc.time()-tmp.time # Originally: 4.35 hours, July 30: 5.02 hours

# Results
Alpha1<-matrix(scan(paste(out,"Alpha1.txt",sep="")),ncol=p,byrow=T)
Alpha2<-matrix(scan(paste(out,"Alpha2.txt",sep="")),ncol=p,byrow=T)
Beta1<-matrix(scan(paste(out,"Beta1.txt",sep="")),ncol=k,byrow=T)
Beta2<-matrix(scan(paste(out,"Beta2.txt",sep="")),ncol=k,byrow=T)
R1<-matrix(scan(paste(out,"R1.txt",sep="")),ncol=1,byrow=T)
R2<-matrix(scan(paste(out,"R2.txt",sep="")),ncol=1,byrow=T)
Sigmab<-matrix(scan(paste(out,"Sigmab.txt",sep="")),nrow=lastit,byrow=T)
Sigma<-matrix(scan(paste(out,"Sigma.txt",sep="")),nrow=lastit,byrow=T)
Sigmab<-matrix(scan(paste(out,"Sigmab.txt",sep="")),nrow=lastit,byrow=T)
Bs1<-matrix(scan(paste(out,"B1.txt",sep="")),nrow=lastit,byrow=T)
Bs2<-matrix(scan(paste(out,"B2.txt",sep="")),nrow=lastit,byrow=T)
Ls<-matrix(scan(paste(out,"L.txt",sep="")),nrow=lastit,byrow=T)

malpha1<-colMeans(Alpha1)
qalpha1<-apply(Alpha1,2,quantile,prob=c(0.025,0.975))
malpha2<-colMeans(Alpha2)
qalpha2<-apply(Alpha2,2,quantile,prob=c(0.025,0.975))

mbeta1<-colMeans(Beta1)
qbeta1<-apply(Beta1,2,quantile,prob=c(0.025,0.975))
mbeta2<-colMeans(Beta2)
qbeta2<-apply(Beta2,2,quantile,prob=c(0.025,0.975))

mr1<-mean(R1)
qr1<-quantile(R1,prob=c(0.025,0.975))
mr2<-mean(R2)
qr2<-quantile(R2,prob=c(0.025,0.975))

msigma<-colMeans(Sigma)
qsigma<-apply(Sigma,2,quantile,prob=c(0.025,0.975))

msigmab<-colMeans(Sigmab)
qsigmab<-apply(Sigmab,2,quantile,prob=c(0.025,0.975))

mb1<-matrix(apply(Bs1,2,mean),ncounty,k,byrow=T)
mb2<-matrix(apply(Bs2,2,mean),ncounty,k,byrow=T)  

########
# WAIC #
########
lhat<-apply(Ls,2,mean)
lpd<-sum(log(lhat))             # -179907.5
pwaic<-sum(apply(log(Ls),2,var))# 1772.519
waic<- -2*(lpd-pwaic)           # 363360 vs 363398  for indep: net diff of 38 points improvement

#######################
# Diagnostics: Fig S1 #
#######################

#pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\FigS1.pdf")
par(mfrow=c(3,2))
plot(1:lastit,Beta1[,1],type="l",col="lightgreen", main=expression(beta[11]),
     xlab="Iteration",ylab=expression(beta[11]))
legend(cex=.7,"topleft",legend=c("ESS = 2000","p = 0.78"))
abline(h=mbeta1[1])
effectiveSize(Beta1[,1])
geweke.diag(Beta1[,1])
2*(1-pnorm(geweke.diag(Beta1[,1])$z))

plot(1:lastit,Beta2[,1],type="l",col="lightgreen",main=expression(beta[21]),
     xlab="Iteration",ylab=expression(beta[21]))
legend(cex=.7,"topleft",legend=c("ESS = 2000","p = 0.88"))
abline(h=mbeta2[1])
effectiveSize(Beta2[,1])
geweke.diag(Beta2[,1])
2*(pnorm(geweke.diag(Beta2[,1])$z))

plot(1:lastit,R1,type="l",col="lightgreen",main=expression(r[1]),
     xlab="Iteration",ylab=expression(r[1]))
legend(cex=.7,"topleft",legend=c("ESS = 1131","p = 0.68"))
abline(h=mr1)
effectiveSize(R1)
geweke.diag(R1)
2*pnorm(geweke.diag(R1)$z)

plot(1:lastit,R2,type="l",col="lightgreen",main=expression(r[2]),
     xlab="Iteration",ylab=expression(r[2]))
legend(cex=.7,"topleft",legend=c("ESS = 917","p = 0.96"))
abline(h=mr2)
effectiveSize(R2)
geweke.diag(R2)
2*(1-pnorm(geweke.diag(R2[,1])$z))

plot(1:lastit,Sigma[,4],type="l",col="lightgreen",main=expression(Sigma[beta[22]]),
     xlab="Iteration",ylab=expression(Sigma[beta[22]]))
legend(cex=.7,"topleft",legend=c("ESS = 2000","p = 0.10"))
abline(h=msigma[4])
effectiveSize(Sigma[,4])
2*(1-pnorm(geweke.diag(Sigma[,4])$z))

plot(1:lastit,Sigmab[,1],type="l",col="lightgreen",main=expression(Sigma[b[11]]),
     xlab="Iteration",ylab=expression(Sigma[b[11]]))
legend(cex=.7,"topleft",legend=c("ESS = 462","p = 0.86"))
abline(h=msigmab[1])
effectiveSize(Sigmab[,1])
2*(pnorm(geweke.diag(Sigmab[,1])$z))
dev.off()

###############################################
# Pop Average and County-Specific Predictions #
###############################################
Pred1<-Pred2<-array(0,dim=c(ncounty,lday,lastit))
for (j in 1:lastit){
  for (h in 1:ncounty){
    Bs21<-matrix(Bs1[j,],ncol=k,byrow=T)
    Bs22<-matrix(Bs2[j,],ncol=k,byrow=T)
    beta1<-Beta1[j,]
    beta2<-Beta2[j,]
    r1<-R1[j]
    r2<-R2[j]
    alpha1<-Alpha1[j,]
    alpha2<-Alpha2[j,]
    Pred1[h,,j]<-r1*c(exp(alpha1[1]*mean(W[id==h,1])+alpha1[2]*mean(W[id==h,2])+Xtmp%*%beta1+Xtmp%*%Bs21[h,]))     # County x Day Predictions (per capita mean per 100k); don't include alpha and W for covariate-adjusted predictions
    Pred2[h,,j]<-r2*c(exp(alpha2[1]*mean(W[id==h,1])+alpha2[2]*mean(W[id==h,2])+Xtmp%*%beta2+Xtmp%*%Bs22[h,]))     # County x Day Predictions (per capita mean per million)
  } 
  print(j)
}

pred1<-colMeans(t(colMeans(Pred1[,,1:lastit])))     # PA prediction (average across counties and iterations)
pred2<-colMeans(t(colMeans(Pred2[,,1:lastit])))     # PA prediction (average across counties and iterations)

countypred1<-apply(Pred1[,,1:lastit],c(1,2),mean)   # County-specific predictions, i.e., post mean for each county/day
countypred2<-apply(Pred2[,,1:lastit],c(1,2),mean)   # County-specific predictions, i.e., post mean for each county/day

#####################################
# Credible Intervals for Mean Trend #
#####################################
# County specific
countylower1<-apply(Pred1[,,1:lastit],c(1,2),quantile,prob=0.025)  # Lower CI for each county 
countyupper1<-apply(Pred1[,,1:lastit],c(1,2),quantile,prob=0.975)  # Average lower bounds across county 

countylower2<-apply(Pred2[,,1:lastit],c(1,2),quantile,prob=0.025)  # Lower CI for each county 
countyupper2<-apply(Pred2[,,1:lastit],c(1,2),quantile,prob=0.975)  # Average lower bounds across county

# Population average
palower1<-apply(t(colMeans(Pred1[,,1:lastit])),2,quantile,prob=0.025)
paupper1<-apply(t(colMeans(Pred1[,,1:lastit])),2,quantile,prob=0.975)

palower2<-apply(t(colMeans(Pred2[,,1:lastit])),2,quantile,prob=0.025)
paupper2<-apply(t(colMeans(Pred2[,,1:lastit])),2,quantile,prob=0.975)

############################
# Population Average Plots #
############################

plot1a_data <- tibble(
  day = 1:lday,
  orr = orr1,
  postmean = pred1,
  palower = palower1,
  paupper = paupper1
)

plot1b_data <- tibble(
  day = 1:lday,
  orr = orr2,
  postmean = pred2,
  palower = palower2,
  paupper = paupper2
)

###############################
# Fig 1a: Pop Avg Plot for y1 #
###############################
#pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig1a.pdf")
ggplot(plot1a_data, aes(x = day, y = orr ))+
  geom_point(aes(color = "Daily Average"))+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = palower, ymax = paupper, color="95% Credible Interval"),alpha=0.3, show.legend = F)+
  scale_color_manual(breaks = c("Daily Average", 
                                "Posterior Mean", 
                                "95% Credible Interval"
  ), 
  values = c("Daily Average"="darkorange",
             "Posterior Mean"="blue4",
             "95% Credible Interval" = "gray40"
  ))+
  ylab("Cases per 100K")+
  #ylim(0,6)+
  guides( color = guide_legend(title="Incidence Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 1, 1),
                                 shape = c(19,NA, NA )),
                               reverse = F))+
  
  theme(legend.key = element_rect(fill = "white"), legend.position = c(0.85, 0.85),legend.text=element_text(size=11))

#dev.off()

###############################
# Fig 1b: Pop Avg Plot for y2 #
###############################
#pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig1b.pdf")
ggplot(plot1b_data, aes(x = day, y = orr ))+
  geom_point(aes(color = "Daily Average"))+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = palower, ymax = paupper, color="95% Credible Interval"),alpha=0.3, show.legend = F)+
  scale_color_manual(breaks = c("Daily Average", 
                                "Posterior Mean", 
                                "95% Credible Interval"
  ), 
  values = c("Daily Average"="darkorange",
             "Posterior Mean"="blue4",
             "95% Credible Interval" = "gray40"
  ))+
  ylab("Deaths per Million")+
  scale_y_continuous(name="Date",breaks=c(0,5,10),limits=c(0,10))+
  # ylim(0,10)+
  guides( color = guide_legend(title="Death Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 1, 1),
                                 shape = c(19,NA, NA )),
                               reverse = F))+
  #ggtitle("Population Average Death Trend Across All Counties")+
  theme(legend.key = element_rect(fill = "white"), legend.position = c(0.85, 0.85),legend.text=element_text(size=11))

#dev.off()

######################################
# Figs 2a-2d: County Specific Plots  #
######################################
#############
# Y1: Cases #
#############
countyid<-sample(1:ncounty,1)  # Random select 
print(countyid) 

# County 70, 151
countyid<-151 # 70 

simcurve1<-exp(W[id==countyid,]%*%truealpha1+Xtmp%*%truebeta1+Xtmp%*%trueB1[countyid,]) # True, Simulated Curve for Y1

# County inc data
county_inc <- tibble(
  day = 1:lday,
  orr1 = incidence1[id==countyid],
  postmean = countypred1[countyid, ],
  countyl = countylower1[countyid, ],
  countyu = countyupper1[countyid, ],
  simcurve = simcurve1
)


# Incidence Plot

# pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig2c.pdf")
ggplot(county_inc, aes(x = day, y = orr1 ))+
  geom_point(aes(color = "Daily Incidence Rate"), shape = 16)+
  geom_line(aes(day, simcurve, color="Simulated Curve"),size=1,linetype=2)+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = countyl, ymax = countyu, color="95% Credible Interval"),alpha=0.3, show.legend = F,
              linetype=1)+
  scale_color_manual(breaks = c("Daily Incidence Rate",
                                "Simulated Curve",
                                "Posterior Mean", 
                                "95% Credible Interval"), 
                     values = c("Daily Incidence Rate"="darkorange",
                                "Simulated Curve"="darkgreen",
                                "Posterior Mean"="blue4",
                                "95% Credible Interval" = "gray40"
                     ))+
  ylab("Cases per 100K")+
  xlab("Day")+
  guides( color = guide_legend(title="Incidence Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 2,1, 1),
                                 shape = c(16,NA,NA, NA )),
                               reverse = F))+
  # ylim(0, 15)+
  # scale_y_continuous(name="Date",breaks=c(0,10,20))+
  # ggtitle("Cases")+
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 11),
        plot.title=element_text(hjust=.5),
        legend.position = c(0.85, 0.85))

# dev.off()

###############
# Y2: Deaths  #
###############
simcurve2<-exp(W[id==countyid,]%*%truealpha2+Xtmp%*%truebeta2+Xtmp%*%trueB2[countyid,]) # True, simulated Curve for y2

county_death <- tibble(
  day = 1:lday,
  orr1 = incidence2[id==countyid],
  postmean = countypred2[countyid, ],
  countyl = countylower2[countyid, ],
  countyu = countyupper2[countyid, ],
  simcurve = simcurve2
)

# Deaths 

# pdf(file="C:\\Brian\\Covid\\Drafts\\Fig\\Fig2d.pdf")
ggplot(county_death, aes(x = day, y = orr1 ))+
  geom_point(aes(color = "Daily Death Rate"), shape = 16)+
  geom_line(aes(day, simcurve, color="Simulated Curve"),size=1,linetype=2)+
  geom_line(aes(day, postmean, color="Posterior Mean"),size=1,linetype=1)+
  geom_ribbon(aes(ymin = countyl, ymax = countyu, color="95% Credible Interval"),alpha=0.3, show.legend = F,linetype=1)+
  scale_color_manual(breaks = c("Daily Death Rate", 
                                "Simulated Curve",
                                "Posterior Mean", 
                                "95% Credible Interval"), 
                     values = c("Daily Death Rate"="darkorange",
                                "Simulated Curve"="darkgreen",
                                "Posterior Mean"="blue4",
                                "95% Credible Interval" = "gray40"
                     ))+
  ylab("Deaths per Million")+
  xlab("Day")+
  guides( color = guide_legend(title="Death Rate Trends",
                               override.aes = list(
                                 linetype = c(0, 2,1, 1),
                                 shape = c(16,NA,NA, NA )),
                               reverse = F))+
  
  #scale_y_continuous(breaks=c(0,10,20,30))+
  #ylim(0, 20)+
  
  #ggtitle("Deaths")+
  theme(legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 11),
        plot.title=element_text(hjust=.5),
        legend.position = c(0.85, 0.85))

# dev.off()

# Save files and workspaces
#save.image("C:\\Brian\\Covid\\R\\Workspaces\\Simulation\\Simulation_07-30-21.RData")
#load("C:\\Brian\\Covi\\R\\Workspaces\\Simulation\\Simulation_07-30-21.RData")
