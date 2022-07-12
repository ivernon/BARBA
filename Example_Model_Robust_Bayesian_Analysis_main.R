##################################################################################################################################
####################               Toy Model for Computer Model Analysis of a Bayesian Analysis               ####################
##################################################################################################################################

### Created using R version 4.1.2 (2021-11-01) -- "Bird Hippie" ###

### Adjust this address to point to directory where "ToyModelMCMC_for_BABA_paper1.R" is ###
#setwd("/Users/ianvernon/Work/A Bayesian Analysis of a Bayesian Analysis/Paper 1/Toy Model Code for release/Toy Model Code Prerelease V2/")

### plotting directory ###
plot_dir <- "plots/"


############################################################################################################################################
### Set up MCMC for Toy Model with contaminated Likelihood, mixture of exponential and Gaussian. Prior has fixed mean but free variance. ###

### Define Prior ###
prior_dens <- function(theta,nu=10,mu_0=5,...){
	dgamma(x=theta,shape=mu_0^2/nu,rate=mu_0/nu)	
}

### Define log prior without proportionality constant ###
log_prior_dens <- function(theta,nu=10,mu_0=5,...) {
	(mu_0^2/nu - 1)*log(theta) - mu_0/nu * theta           # log density of equation (19)
}


### Load data used in paper for toy model (given in Appendix A, Table 2)
load("Simulated Data and MCMC Output/Sim_Data1.Rdata") 			# load data used for paper
dat                                                         # the data z_i, given in Sup Mat Appendix A Table 2


### Define log likelihood ###
LogLikeHood <- function(theta,x=dat,epsilon=0,...){
	if(epsilon>1) epsilon <- 1	else             # catch values of epsilon outside of allowed [0,1] range
	if(epsilon<0) epsilon <- 0
	n <- length(x)
	n*log(theta) + sum(  log( (1-epsilon)*exp(-theta*x) + epsilon*2/pi*exp(-theta^2*x^2/pi) )  )		# log of equation (18)
}


### Evaluate Likelihood over range of theta and epsilon values, to make Figure1_left ###
th_seq <- seq(0.5,3.5,len=100)
ep_seq <- seq(0,1,len=21)
llh_mat <- matrix(0,length(ep_seq),length(th_seq))
for(i in 1:length(th_seq)) {
	for(j in 1:length(ep_seq)) llh_mat[j,i] <- LogLikeHood(theta=th_seq[i],epsilon=ep_seq[j])
	}

### Generate Figure1_left_log_likelihoods.pdf ###	
pdf(paste(plot_dir,"Figure1_left_log_likelihoods.pdf",sep=""),width=6.1,height=6.1)
plot(th_seq,llh_mat[1,],ty="n",col=2,ylim=range(llh_mat),xlab=expression(theta),ylab="Log Likelihood")
for(j in 1:length(ep_seq)) lines(th_seq,llh_mat[j,],col=rainbow(length(ep_seq),st=0,en=5/6)[j],lwd=2)
ea_s <- c(1,7,13,20)
legend("bottomright",legend=c(expression(epsilon*" = "*0),expression(epsilon*" = "*0.3),
								expression(epsilon*" = "*0.6),expression(epsilon*" = "*1)), 
								col=rainbow(length(ep_seq),st=0,en=5/6)[ea_s],lty=1,lwd=2)
dev.off()


### Define simple MCMC Metropolis Hastings Algorithm: single step ###
mcmc_step <- function(theta_i,sig_prop=0.3,...){
	
	theta_prop <- abs( rnorm(1,mean=theta_i,sd=sig_prop) ) 			              # folded normal proposal
	
	accept_ratio <- (log_prior_dens(theta_prop,...) + LogLikeHood(theta=theta_prop,...) 
					- log_prior_dens(theta_i,...) - LogLikeHood(theta=theta_i,...) )  # Classic MH acceptance ratio formula

	u <- runif(1)
	if(accept_ratio > log(u)) return(theta_prop) else return(theta_i)         # accept or reject proposed point theta_prop
}

### Define simple MCMC Metropolis Hastings Algorithm: multiple steps ###
mcmc_alg <- function(theta_init=1,Nsteps=100,thin=1,...){
	theta_draws <- rep(0,Nsteps)                   # set up empty MCMC chain
	theta_draws[1] <- theta_init                   # initialise first element of chain
	for(i in 1:(Nsteps-1)) theta_draws[i+1] <- mcmc_step(theta_i=theta_draws[i],...)      # run many MH MCMC steps
	cat("Accept Prob = ",1-sum(diff(theta_draws)==0)/length(diff(theta_draws)),"\n" )			# output acceptance probability
	out <- theta_draws[seq(1,Nsteps,thin)]         # thinning step if needed
	return(list(out,c(mean(out[-1]),sd(out[-1])))) # return the thinned set of theta values (thinned MCMC chain) and summaries (which discard initial element of chain)
}


### run the MCMC algorithm at a single point to prepare for Figure1_right ###
set.seed(seed=3)
out <- mcmc_alg(theta_init=0.5,Nsteps=200000,thin=10,sig_prop=0.9,epsilon=0,nu=1^2)[[1]]			


### Prep for Figure1_right: function to plot posterior and prior pdfs and means in conjugate epsilon = 0 case ###
plot_conj <- function(xdat,nu=10,mu_0=5,inc_prior=0,add=TRUE,inc_mean=1){
  n <- length(xdat)
  curve(dgamma(x,shape=mu_0^2/nu + n,rate=mu_0/nu + sum(xdat)),col="blue",add=add,lwd=3)    # plot posterior pdf
  if(inc_mean==1) abline(v=(mu_0^2/nu + n)/ (mu_0/nu + sum(xdat)),lty=2,col="blue",lwd=3)   # vertical line at pasterior mean
  if(inc_prior==1) curve(dgamma(x,shape=mu_0^2/nu,rate=mu_0/nu),col="red",add=add,lwd=3)    # plot prior pdf
  if(inc_mean==1) lines(c(mu_0,mu_0),c(-0.025,0.5),lty=2,col="red",lwd=3)                   # vertical line at prior mean
}

### Generate Figure1_right_posterior_prior_densities: plot the results of the MCMC at single point and compare to conjugate results ###
pdf(paste(plot_dir,"/Figure1_right_posterior_prior_densities.pdf",sep=""),width=7.3,height=5.9)
hist(out,br=seq(0,10,0.1),freq=FALSE,xlab=expression(theta),main="",col="grey80",xlim=c(1.8,6.3),
				ylab="Posterior or Prior Probability Density")
plot_conj(xdat=dat,inc_prior=1,nu=1^2)
legend("topright",legend=c("Conjugate Prior","Conjugate Posterior","Thinned MCMC Draws","Posterior Mean","Prior Mean"),lwd=c(2,2,5,2,2),col=c("red","blue","grey80","blue","red"),cex=0.85,lty=c(1,1,1,2,2))
dev.off()

############################################################################################################################################




#######################################################################################################
####### Create 2d space filling Lattice design ########

### Create simple 2-dim lattice design. Simple and ensures space filling
### but undesirable projection properties (in certain directions). OK here, but
### would typically use a maximin LHC via e.g. library(SLHD) for higher dimensional cases.

x <- c(0,0)
g1 <- c(5,2)   # define generators of latice
g2 <- c(2,5)   # define generators of latice

s1 <- seq(-15,15,len=41)
s2 <- seq(-15,15,len=41)

# set up lattice points 
xmat <- matrix(0,nrow=length(s1)*length(s2),ncol=2)
for(i in 1:length(s1)){
	for(j in 1:length(s2)){
		xmat[(i-1)*length(s1)+j,] <- x+s1[i]*g1 + s2[j]*g2
	}
}

# allowed points within  suitable box
al <-  xmat[,1] >= -10 & xmat[,1] <= 10 &   xmat[,2] >= -10 & xmat[,2] <= 10   

# rescale to correct range #
x_ran <- matrix(c(0.3,sqrt(4), 0,1),nrow=2, byrow=TRUE)    # Ranges of input space X
x_mean <- apply(x_ran,1,mean)                              # Center of ranges
x_dif <- apply(x_ran,1,diff)                               # Width of ranges

xmat_sc <- t(t(xmat[al,])/10*1.1*x_dif/2 + x_mean)			 # the 1.1 is a 20 percent enlargement so we have points slightly outside boundary


#######################################################################################################


#######################################################################################################
### Run the MCMC algorithm at each point in the design ###

x1 <- xmat_sc                              # the final design
nps <- nrow(x1)                            # number of run points 
y1 <- matrix(0,nrow=nps,ncol=2)            # set up matrix of MCMC output: just posterior means and variances
yfull <- matrix(0,nrow=nps,ncol=20000)     # set up matrix of MCMC output: full chain output

### Run MCMC at each design point ###
full_evaluation <- 0     # select this to load up pre evaluated 35 MCMC runs
# full_evaluation <- 1   # select this to perform fresh full 35 MCMC runs (takes 1.2 mins or so)

if(full_evaluation){
  for(i in 1:nps) {
		out <- mcmc_alg(theta_init=0,Nsteps=200000,thin=10,sig_prop=0.9,epsilon=x1[i,2],nu=x1[i,1]^2)   # MCMC algorithm
		yfull[i,] <- out[[1]]      # store full MCMC thinned chain output
		y1[i,] <- out[[2]]         # store just posterior mean and sd MCMC output
		cat(i,"\n")
	}
#save(x1,y1,nps,file="Simulated Data and MCMC Output/Lattice_design_and_MCMC_out_20000.Rdata")
#save(yfull,x1,y1,nps,file="Simulated Data and MCMC Output/Lattice_design_and_MCMC_out_20000_fullout.Rdata")
}
if(!full_evaluation){
load("Simulated Data and MCMC Output/Lattice_design_and_MCMC_out_20000.Rdata")
load("Simulated Data and MCMC Output/Lattice_design_and_MCMC_out_20000_fullout.Rdata")
}

### Check MCMC outputs ###
cbind(x1,y1)      # This gives the results in Table 1 of the Supp.Mat.

### End Run the MCMC algorithm at each point in the design ###
#######################################################################################################


#######################################################################################################
### Plot MCMC output: mixing and posteriors

### Generate Figure2_right.pdf: Check mixing plots for runs in corners of 2D space ###
pdf(paste(plot_dir,"Figure2_right_mixing_plots.pdf",sep=""),width=6.1,height=6.1)
op <- par(mfrow=c(2,2),mar=c(4,4,0.4,0.4))
for(i in c(1,26,10,35)){		# Run number for runs in topleft, topright, bottomleft, bottomright of input space.
	plot(yfull[i,0+seq(1,2000,1)],ty="l",col="grey40",ylab=expression(theta),xlab="step",ylim=c(0.7,5.7))
}
dev.off()
par(op)  

### Generate Figure2_left.pdf: Plot all 35 MCMC posteriors using density estimate ###
pdf(paste(plot_dir,"Figure2_left_MCMCposteriors.pdf",sep=""),width=6.1,height=6.1)
plot(1:5,ty="n",xlim=c(1.1,5.6),ylim=c(0,1.8),xlab=expression(theta),ylab=expression("MCMC Posterior  "*pi(theta*"|"*z,nu,epsilon)))
k <- 1
for(i in order(x1[,1])){
	lines(density(yfull[i,seq(100,20000,1)]),col=rainbow(nrow(yfull))[k],lwd=2)
	k <- k+1
}
dev.off()

### End Plot MCMC output: mixing and posteriors
#######################################################################################################




#######################################################################################################
### Simple emulation construction for posterior mean and posterior SD

library(pdist)    # Package pdist version 1.2

### Define simple constant plus GP emulator ###
### Note: this is a simple emulator provided for demonstration/tutorial purposes. For higher dimensional 
### applications we recommend any of the many GP packages (e.g RobustGaSP, hmer, MOGP, GPfit etc.)

simple_emulator <- function(x1,                  # matrix of design points: locations in input space where MCMC runs were performed
                            y1,                  # vector of run outputs (e.g. posterior means or posterior sds) corresponding to x1
                            xp,                  # prediction points where the emulator will be evaluated
                            theta=1,             # emulator correlation length parameter
                            nugget=0.1,          # emualtor nugget parameter
                            xran,                # physical range of X space, used to convert to [-1,1] scale
                            just_var=1,          # set to 1 if just want emulator means and variances, set to 0 if want means and full covariance matrix across xp
                            em_mcmc_notf=0       # set to 0 to get predictions for f, set to 1 to get predictions for noisy MCMC output, e.g. for diagnostics
                            ){     
	x1 <- scale(x1,center=apply(xran,1,mean),scale=apply(xran,1,diff)/2)     # scale design inputs to standard [-1,1] scale
	xp <- scale(xp,center=apply(xran,1,mean),scale=apply(xran,1,diff)/2)     # scale new prediction point inputs to standard [-1,1] scale
	D <- y1                             # store design outputs as data vector D
	mu <- mean(y1)                      # estimate mean as mean of runs XXX JP: is it simpler just to standardise the runs? Then specify mean 0, var 1 for emulator?
	sigma <- sd(y1)                     # estimate sigma as sd of runs
	Cmat <- function(d) sigma^2 * exp(-d^2/theta^2) * (1-nugget)    # define squared exponential covariance function: see equation (10)
	covBD <- Cmat( as.matrix(pdist(xp,x1)) )      # Covariance matrix between new outputs B = f(xp) and design points D = f(x1): see equation (11)
	varD  <- Cmat( as.matrix(dist(x1)) ) + sigma^2 * diag(nugget,nrow=nrow(x1))    # Variance matrix for design points D = f(x1), restoring nugget to diagonal: see equation (10)
	varB  <- Cmat( as.matrix(dist(xp)) ) + em_mcmc_notf * sigma^2 * diag(nugget,nrow=nrow(xp))		# Var for prediction points B = f(xp). Note: em_mcmc_notf = 1 => variance includes mcmc noise for diagnostics, em_mcmc_notf = 0 just emulating f
	EB    <- mu													# set prior expectation for B					
	ED 	  <- mu                         # set prior expectation for D
	varD_inv <- chol2inv(chol(varD))    # invert varD using Choleski decompostion
	ED_B   <- EB   + covBD %*% varD_inv %*% (D - ED)      # Mean update: see equation (7)
	VarD_B <- varB - covBD %*% varD_inv %*% t(covBD)      # Covariance matrix update: see equation (8)
	if(just_var)  return(list(ED_B,diag(VarD_B))) else    # return emulator means and variances at prediction points xp
	if(!just_var) return(list(ED_B,VarD_B))               # or return emulator means and full covariance matric for prediction points xp
}


### set up grid of prediction points xp to evaluate emulator at, for plotting ###
xran <- matrix(c(0.3,sqrt(4), 0,1),nrow=2, byrow=TRUE)
sd_seq <- seq(0.3,sqrt(4),len=40)
ep_seq <- seq(0,1,len=40)
xp <- as.matrix(expand.grid(sd_seq,ep_seq))




################## Posterior Mean emulator plots ##################

### can estimate nugget by averaging sample mean error (sd/sqrt(20000)) over all 35 runs 
(mean(y1[,2]/sqrt(20000)) / sd(y1[,1]) )^2   # gives 0.0000178..., a guide to choosing nugget

### emulate posterior mean ###
em_out <- simple_emulator(x1=x1,y1=y1[,1],xp=xp,theta=0.7,nugget=0.0000178,xran=xran)		# employ the emulator for posterior mean
em_exp <- matrix(em_out[[1]],length(sd_seq),length(ep_seq))      # store emulator expectations for posterior mean

### plot emulator expecation of posterior mean ###
sc <- 0.9
levs <- seq(2.1,4.8,0.1)
pdf(paste(plot_dir,"Figure3_left.pdf",sep=""),width=7.2*sc,height=6.5*sc)
filled.contour(sd_seq,ep_seq,em_exp,color.palette=terrain.colors,
               xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
               main="",level=levs,
               plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col="blue") 		   # plot original specifcation point in blue
                 points(1.5,0.5,pch=3,cex=2,col=1,lwd=2.5)									               # plot case 1 as single point
                 arrows(x0=0.8,y0=0.003,y1=0.997,code=3,angle=90,col=1,len=0.1,lwd=2.5)    # plot case 2 as line
                 arrows(x0=0.5,x1=1.9,y0=0.72,code=3,angle=90,col=1,len=0.1,lwd=2.5)       # plot case 3 as line
                 xp1 <- seq(0.70001,1.29999,len=200)                                       
                 lines(xp1,0.4*sqrt(1-(xp1-1)^2/0.3^2),col=1,lwd=2.5)                      # plot case 4 as semi-ellipsoid
               }
)
dev.off()



################## Posterior SD emulator plots ##################

# estimate using numerical estimation of sd estimators SE for sample size 20000 from normal dist ###
sd(apply(matrix(rnorm(20000*200,mean=0,sd=1),nrow=200,ncol=20000),1,sd)) # gives about 0.004..., a guide to choosing nugget

### emulate posterior SD ###
em_out <- simple_emulator(x1=x1,y1=y1[,2],xp=xp,theta=0.7,nugget=0.004,xran=xran,em_mcmc_notf=0)			# employ the emulator for posterior SD
em_exp <- matrix(em_out[[1]],length(sd_seq),length(ep_seq))   # store emulator expectations for posterior SD

### plot emulator expecation of posterior SD ###
levs <- seq(0.22,0.62,0.02)
pdf(paste(plot_dir,"Figure3_right.pdf",sep=""),width=7.2*sc,height=6.5*sc)
filled.contour(sd_seq,ep_seq,em_exp,color.palette=terrain.colors,
               xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
               main="",
               level=levs,
               plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col="blue")}
)
dev.off()


### End Emulation construction for posterior mean and posterior SD
#######################################################################################################################################




#######################################################################################################################################
### Perform Quantile Emulation for Supplementary Material 						

load("Simulated Data and MCMC Output/Lattice_design_and_MCMC_out_20000_fullout.Rdata") 			# load full posterior data

### evaluate quantiles from the 35 MCMC runs, discarding initial element in MCMC chain ###
q05 <- apply(yfull[,-1],1,quantile,probs=c(0.05))
q25 <- apply(yfull[,-1],1,quantile,probs=c(0.25))
q75 <- apply(yfull[,-1],1,quantile,probs=c(0.75))
q95 <- apply(yfull[,-1],1,quantile,probs=c(0.95))
IQR <- q75 - q25
Var_over_Mean <- y1[,2]^2/y1[,1]

### emulate and plot 5th, 25th, 75th and 95th quantiles, IQR and ratio of posterior Var over Mean ###
plot_nam <- c("Q05","Q25","Q75","Q95","IQR","Var_over_Mean")
qx <- cbind(q05,q25,q75,q95,IQR,Var_over_Mean)     # construct matrix of outputs of interest

### emualte each output of interest ###
for(i in 1:6){
	em_out <- simple_emulator(x1=x1,y1=qx[,i],xp=xp,theta=0.7,nugget=0.00001,xran=xran)    # perform emulation
	em_exp <- matrix(em_out[[1]],length(sd_seq),length(ep_seq))                            # convert to matrix for plotting
	levs <- seq(min(em_exp),max(em_exp),len=20)
	pdf(paste(plot_dir,"SupMat_Figure3_",plot_nam[i],".pdf",sep=""),width=7.2*sc,height=6.5*sc)
	filled.contour(sd_seq,ep_seq,em_exp,color.palette=terrain.colors,
					xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
						main="",level=levs,plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col="blue")})
	dev.off()
}

### End Perform Quantile Emulation for Supplementary Material 						
#######################################################################################################################################




##############################################################################################################################
### Predictions plus uncertainty for several example specifications: Cases 1 to 4 in paper ###
### 4 cases, 1 point, two lines, one half ellipse ###

library(MASS)     # Package MASS version 7.3-54
library(xtable)   # Package xtable version 1.8-4

### MCMC BIT: slow ###
full_evaluation <- 0     # select this to load up pre evaluated MCMC evaluations to compare with cases 2 and 3
# full_evaluation <- 1   # select this to perform full 2x40 MCMC evaluations to compare with cases 2 and 3

if(full_evaluation){
  load("Simulated Data and MCMC Output/Sim_Data1.Rdata")
  set.seed(1)
  ymcmc_case1 <- matrix(0,nrow=length(ep_seq),ncol=2)
  for(i in 1:length(ep_seq)) {
    ymcmc_case1[i,] <- mcmc_alg(theta_init=0,Nsteps=800000,thin=10,sig_prop=0.9,epsilon=ep_seq[i],nu=0.8^2)[[2]]
    cat(i,"\n")
  }
  ymcmc_case2 <- matrix(0,nrow=length(sd_seq),ncol=2)
  for(i in 1:length(sd_seq)) {
    ymcmc_case2[i,] <- mcmc_alg(theta_init=0,Nsteps=800000,thin=10,sig_prop=0.9,epsilon=0.72,nu=sd_seq[i]^2)[[2]]
    cat(i,"\n")
  }
  # save(ymcmc_case1,ymcmc_case2,file="Simulated Data and MCMC Output/Case2_3_straightlines_MCMCdiag_reps80000.Rdata")
}
if(!full_evaluation) load("Simulated Data and MCMC Output/Case2_3_straightlines_MCMCdiag_reps80000.Rdata")


### details of output and emulation ###
set.seed(1)
add_diag <- 1				    # if we want to add the diagnostic points
for(out_mv in 1:2){			# do full loop over both outputs: posterior mean (out_mv=1) and posterior sd (out_mv=2)
  thet <- 0.7				    # correlation length
  nug <- c(0.0000178,0.004)    # nuggets for posterior mean and posterior SD outputs

  # Case a) single point, but wants derivatives (see extended emulator below for derivative calculations)
  nsim <- 10000                         # number of draws for Fig. 5
  xp <- cbind(c(1.5,1.5),c(0.5,0.5))    # single input point in Case 1 in paper (emulator wants > 1 input point so double up: fix?)
  em_out <- simple_emulator(x1=x1,y1=y1[,out_mv],xp=xp,theta=thet,nugget=nug[out_mv],xran=xran,just_var=0,em_mcmc_notf=0)
  yp_casea    <- em_out[[1]][1,1]                           # store emulator mean
  yp_sd_casea <- sqrt(diag(em_out[[2]]))[1]                 # store emulator sd
  em_draws <- rnorm(nsim,mean=yp_casea,sd=yp_sd_casea)      # simulate lots of draws of f(x) for Fig. 5.
  

  ### next two cases: single lines ###
  nsim <- 10000               # number of draws foor Fig. 5
  max1 <- matrix(0,nrow=nsim,ncol=3)    # create matrix for storing max values over region: M_j = max_i f^(j)(x^(i)_E) (see equation (26))
  min1 <- matrix(0,nrow=nsim,ncol=3)    # create matrix for storing min values over region: m_j = max_i f^(j)(x^(i)_E) (see equation (26))
  midpoint <- 1:3                       # vector to store midpoint of cases 2, 3 and 4 for Fig. 5.
  
  for(case in 1:2){     # note that the cases enumerated in paper = case + 1
    
    if(case==1){
      # Case b) vary epsilon in [0.1], nu fixed at 0.8 #
      xp <- cbind(0.8,seq(0,1,len=101))    # points to evaluate emulator at
      xg <- xp[,2]                         # points used for plotting
    } else if(case==2){
      # Case c) vary nu in [0.5,1.9], epsilon fixed at 0.72 #
      xp <- cbind(seq(0.5,1.9,len=101),0.72)     # points to evaluate emulator at
      xg <- xp[,1]                               # points used for plotting
    }
    
    # emulate and plot emulator with uncertainty
    em_out <- simple_emulator(x1=x1,y1=y1[,out_mv],xp=xp,theta=thet,nugget=nug[out_mv],xran=xran,just_var=0,em_mcmc_notf=0)
    ypred <- em_out[[1]][,1]                           # store emulator expecations (predictions)
    yp_covar <- em_out[[2]]                            # store emulator covariance matrix over set of xp points
    yp_sd <- sqrt(diag(em_out[[2]]))                   # store emulator sds at xp points
    # Plot panels for Figure 4 (2 each loop over output out_mv: emulating mean or sd)
    pdf(paste(plot_dir,"Figure4_Emul_out",out_mv,"_case",case+1,".pdf",sep=""),width=0.9*6.1,height=0.9*6.1)
    ex <- c(expression(f[1]*(x)*" "==" "*E(theta*" | "*z,nu,epsilon)),expression(f[2]*(x)*" "==" "*SD(theta*" | "*z,nu,epsilon)))[out_mv]  # y-axis lable 
    plot(xg,ypred,ty="l",col=4,ylim=range(ypred+2*yp_sd,ypred-2*yp_sd),lwd=2,
         xlab=c(expression(epsilon),expression(nu))[case],ylab=ex)
    # Add black diagnostics points evaluated using full MCMC above
    if(add_diag==1){
      if(case==1) points(ep_seq,ymcmc_case1[,out_mv],pch=16)	
      if(case==2) points(sd_seq,ymcmc_case2[,out_mv],pch=16)
    }
    lines(xg,ypred,col="blue",lwd=2)               # emulator expectation 
    lines(xg,ypred+2*yp_sd,col="red",lwd=2)        # emulator credible interval (upper)
    lines(xg,ypred-2*yp_sd,col="red",lwd=2)        # emulator credible interval (lower)
    leg_loc <-matrix(c("bottomleft","topright","bottomleft","bottomright"),nrow=2,byrow=TRUE)
    if(add_diag==0) legend(leg_loc[out_mv,case],legend=c("Emulator Mean","Emulator 95% Credible Interval"),lty=c(1,1),col=c("blue","red"),cex=0.8,lwd=2) 
    if(add_diag==1) legend(leg_loc[out_mv,case],legend=c("Emulator Mean","Emulator 95% Credible Interval","Left-Out Diagnostic Points"),lty=c(1,1,0),col=c("blue","red",1),cex=0.8,lwd=2,pch=c(-1,-1,16)) 
    
    dev.off()
    
    # simulate from many possible f(x) curves across set of x_E points using multivariate normal (see equations (24) and (25)) #
    gp_reals <- mvrnorm(n=nsim,mu=ypred,Sigma=yp_covar)
    # Calculate max and min from each simulated f(x) realisation, and store n_S realisations of max and min:
    max1[,case] <- apply(gp_reals,1,max)
    min1[,case] <- apply(gp_reals,1,min)
    # Store emulator expecation at midpoint of region for later plot
    midpoint[case] <- ypred[ceiling(nrow(xp)/2)]
  }
  
  # case 3: region is half ellipsoid around conjugate analysis (defined by equation (23)) #
  case <- 3             # note that the cases enumerated in paper = case + 1, so this is case 4 in paper
  len <- 47             # create fine grid of points over X of size (len x len)
  xp <- expand.grid(seq(xran[1,1],xran[1,2],len=len),seq(xran[2,1],xran[2,2],len=len))  # create the grid
  a <- 0.3				      # define ellipse major (a) and minor (b) axes lengths
  b <- 0.4
  cen <- c(1,0)					# centre of ellipse: previous conjugate analysis, defined by equation (23)
  elip <- ( (xp[,1] - cen[1])^2/a^2 + (xp[,2]-cen[2])^2/b^2  ) < 1  # logical: grid points that lie inside ellipse
  # plot ellipse just to check #
  plot(1:10,xlim=xran[1,],ylim=xran[2,],ty="n")
  points(xp,cex=0.8,col="grey70")                          # plot all grid points
  points(xp[elip,],col=2,pch=16,cex=0.8)      # just plot points in red that are inside ellipse 
  sum(elip)							                      # check number of points in ellipse
  xp <- xp[elip,]                             # redefine prediction points to be just those inside ellipse
  
  # emulate at points xp inside ellipse #
  em_out <- simple_emulator(x1=x1,y1=y1[,out_mv],xp=xp,theta=thet,nugget=nug[out_mv],xran=xran,just_var=0,em_mcmc_notf=0)
  ypred <- em_out[[1]][,1]     # store emulator expectations (predictions)
  yp_covar <- em_out[[2]]      # store emulator covariance matrix over set of xp points
  # simulate from many possible f(x) curves across set of x_E points using multivariate normal (see equations (24) and (25)) #
  gp_reals <- mvrnorm(n=nsim,mu=ypred,Sigma=yp_covar)
  # Calculate max and min from each simulated f(x) realisation, and store n_S realisations of max and min:
  max1[,case] <- apply(gp_reals,1,max)
  min1[,case] <- apply(gp_reals,1,min)
  # Store emulator expecation at midpoint of region which is actually conjugate point #
  xp <- rbind(c(1,0),c(1,0))			# conjugate coordinate
  midpoint[case] <- simple_emulator(x1=x1,y1=y1[,out_mv],xp=xp,theta=thet,nugget=nug[out_mv],xran=xran,just_var=0,em_mcmc_notf=0)[[1]][1,1]
  
  
  ### Now do plot of imprecise ranges for figure 5, plus emulator uncertainty on the max and min estimates for SupMat Figure 2###
  boxnum <- 1000
  loc <- c(1,1.5,2)
  pdf(paste(plot_dir,"Figure_5_All_cases_out",out_mv,".pdf",sep=""),width=6.1,height=6.1)    # Figure 5
  ex <- c(expression(f[1]*(x)*" "==" "*E(theta*" | "*z,nu,epsilon)),expression(f[2]*(x)*" "==" "*SD(theta*" | "*z,nu,epsilon)))[out_mv]
  plot(c(1,1),xlim=c(0.3,2.4),ylim=range(min1,max1),ty="n",ylab=ex,xaxt="n",xlab=" ")
  for(case in 1:3){                       # Actually Cases 2:4 in paper
    arrows(loc[case],mean(min1[,case]),loc[case],mean(max1[,case]),code=3,angle=90,len=0.1,col="blue",lwd=3)  # Blue error bars giving ranges
    boxplot(cbind(min1[1:boxnum,case],max1[1:boxnum,case]),add=TRUE,col="red",at=c(0.1,0.1)+loc[case],outline=TRUE,lwd=1.,
            pars=list(boxwex=c(0.07,0.07,0.07),staplewex=c(1,1,1),outcex=c(0.5,0.5,0.5),
                      outpch=c(1,1,1),whisklty=1),names=c("",""),xaxt="n")    # red boxplots giving uncertainty on max and min due to emulation
    points(loc[case],midpoint[case],pch=16,col="blue",cex=1.5)                # blue point giving midpoint result
  }       
  # case a boxplot: Case 1 in paper
  loca <- 0.5
  arrows(loca,yp_casea,loca,yp_casea*1.001,code=3,angle=90,len=0.1,col="blue",lwd=3)    # Blue error bars giving ranges: zero ranges here (Case 1 in paper)
  points(loca,yp_casea,pch=16,col="blue",cex=1.5)                                       # blue point giving single Case 1 result
  boxplot(cbind(em_draws[1:boxnum],em_draws[1:boxnum]),add=TRUE,col=2,at=c(0.1,10)+loca,outline=TRUE,lwd=1.,
          pars=list(boxwex=c(0.07,0.07),staplewex=c(1,1),outcex=c(0.5,0.5),outpch=c(1,1),whisklty=1),names=c("",""),xaxt="n") # red boxplot giving uncertainty on result due to emulation
  axis(side=1, at=c(loca,loc), labels=paste("Case",1:4), las=2, tck=-0.01)
  legend("bottomright",legend=c("Expected Range","Emul. Uncertainty","Midpoint Value"),
         lty=c(1,1,0),pch=c(1,22,16),pt.cex=c(0,1.3,1.1),col=c("blue",1,"blue"),lwd=c(3,1.2,0),cex=0.8,pt.bg="red")
  dev.off()    # end Figure 5
  
  
  ### Now plot joint pdf of max and min estimates for SupMat Figure 4 ###
  if(out_mv==2){            # just do SD output
    grey.cols <- function(n) grey(n:1/n)
    hsca <- 10 								# divide range into rough h value
    for(case in 1:2){         # Case 2 and 3 in paper
      pdf(paste(plot_dir,"SupMat_Figure4_Bivar_den_out",out_mv,"_case_",case+1,".pdf",sep=""),
          width=0.9*8.0,height=0.9*7.1)    # SupMat Figure 4
      den <- kde2d(min1[,case],max1[,case],n=65,h=c(diff(range(min1[,case]))/hsca,diff(range(max1[,case]))/hsca),
                   lims=c(range(min1[1:boxnum,case]),range(max1[1:boxnum,case])))   # kernel density smoothing of bivariate density
      filled.contour(den,color.palette=grey.cols,lev=c(0,seq(max(den[[3]])/50,max(den[[3]]),len=13)),xlab="Minimum",ylab="Maximum")
      dev.off()
    } 
  }  # end SupMat Figure 4
}

### End Predictions plus uncertainty for several example specifications: Cases 1 to 4 ###
###########################################################################################################




#######################################################################################################################################
### Perform large grid of MCMC runs for comparison with emulator results ###			

### Evaluate 40x40 grid. Note: takes awhile! ###
full_evaluation <- 0     # select this to load up pre evaluated 
# full_evaluation <- 1   # select this to perform full 40x40 grid of MCMC evaluations (takes awhile: 50 mins approx on single core)
load("Simulated Data and MCMC Output/Sim_Data1.Rdata")
if(full_evaluation){
  y2_mat <- array(0,c(length(sd_seq),length(ep_seq),2))
  for(j in 1:length(sd_seq)) {
    for(i in 1:length(ep_seq)) {
      y2_mat[j,i,] <- mcmc_alg(theta_init=0,Nsteps=200000,thin=10,sig_prop=0.9,epsilon=ep_seq[i],nu=sd_seq[j]^2)[[2]]
      cat(i,j,"\n")
    }
  }
  #save(sd_seq,ep_seq,y2_mat,dat,file="Simulated Data and MCMC Output/Case1B_in_eps_var_out_mean_var_40_40_20000.Rdata")
} 
if(!full_evaluation) load("Simulated Data and MCMC Output/Case1B_in_eps_var_out_mean_var_40_40_20000.Rdata")  

### Posterior Mean plot using large 40x40 grid of MCMC runs: no emulation involved ###
levs <- seq(2.1,4.8,0.1)
pdf(paste(plot_dir,"SupMat_Figure2_topleft.pdf",sep=""),width=7.2*sc,height=6.5*sc)
filled.contour(sd_seq,ep_seq,y2_mat[,,1],color.palette=terrain.colors,
               xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
               main=""
               ,level=levs
               ,plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col="blue") 		
               }
)
dev.off()

### Posterior SD plot using large 40x40 grid of MCMC runs: no emulation involved ###
levs <- seq(0.22,0.62,0.02)
pdf(paste(plot_dir,"SupMat_Figure2_topright.pdf",sep=""),width=7.2*sc,height=6.5*sc)
filled.contour(sd_seq,ep_seq,y2_mat[,,2],color.palette=terrain.colors,
               xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
               main=""
               ,level=levs
               ,plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col="blue") 		
               }
)
dev.off()

### End Perform large grid of MCMC runs for comparison with emulator results ###			
#######################################################################################################################################



#######################################################################################################################################
### Perform Posterior Mean Emulator Diagnostics using large 40x40 grid of MCMC runs done above ###			

### load large 40x40 grid of MCMC runs done above ###
load("Simulated Data and MCMC Output/Case1B_in_eps_var_out_mean_var_40_40_20000.Rdata")		
dim(y2_mat)     # array with grid of 40x40 points in, giving both posterior means and SDs.

### confirm predictions points on 40x40 grid ###
xp <- as.matrix(expand.grid(sd_seq,ep_seq))

### evaluate posterior mean emulator on 40x40 grid of diagnostic points ###
em_out <- simple_emulator(x1=x1,y1=y1[,1],xp=xp,theta=0.7,nugget=0.0000178,xran=xran,em_mcmc_notf=1)		# new and more accurate spec (nugget from mcmc sig^2/n etc)
em_exp <- matrix(em_out[[1]],length(sd_seq),length(ep_seq))     # store emulator expectation
em_var <- matrix(em_out[[2]],length(sd_seq),length(ep_seq))			# store emulator variance

### construct emualtor diagnostics ###
em_diag <- (em_exp - y2_mat[,,1])/sqrt(em_var)            # construct diagnostic, note mcmc estimation error now already included in emulator variance

### Posterior Mean Emulator Diagnostic Plot Sup Mat Figure 2 bottom left ###
pdf(paste(plot_dir,"SupMat_Figure2_bottomleft_Diagnostics_mean1.pdf",sep=""),width=7.2*sc,height=6.5*sc)
filled.contour(sd_seq,ep_seq,em_diag,color.palette=terrain.colors,
               xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
               lev=seq(-3.5,3.5,0.1),
               main="",plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col="blue")})
dev.off()

### End Perform Posterior Mean Emulator Diagnostics using large 40x40 grid of MCMC runs done above ###			
#######################################################################################################################################



#######################################################################################################################################
### Perform Posterior SD Emulator Diagnostics using large 40x40 grid of MCMC runs done above ###			

### load large 40x40 grid of MCMC runs done above ###
load("Simulated Data and MCMC Output/Case1B_in_eps_var_out_mean_var_40_40_20000.Rdata")		# these are calculated below
dim(y2_mat)     # array with grid of 40x40 points in, giving both posterior means and SDs.

### confirm predictions points on 40x40 grid ###
xp <- as.matrix(expand.grid(sd_seq,ep_seq))

### evaluate posterior SD emulator on 40x40 grid of diagnostic points ###
em_out <- simple_emulator(x1=x1,y1=y1[,2],xp=xp,theta=0.7,nugget=0.004,xran=xran,em_mcmc_notf=1)			# new and more accurate spec (nugget from var normal sim)
em_exp <- matrix(em_out[[1]],length(sd_seq),length(ep_seq))     # store emulator expectation
em_var <- matrix(em_out[[2]],length(sd_seq),length(ep_seq))			# store emulator variance

### construct emualtor diagnostics ###
em_diag  <- (em_exp - y2_mat[,,2])/sqrt(em_var)                 # construct diagnostic, note mcmc estimation error now already included in emulator variance

### Posterior SD Emulator Diagnostic Plot Sup Mat Figure 2 bottom right ###
pdf(paste(plot_dir,"SupMat_Figure2_bottomright_Diagnostics_SD1.pdf",sep=""),width=7.2*sc,height=6.5*sc)
filled.contour(sd_seq,ep_seq,em_diag,color.palette=terrain.colors,
               xlab=expression("Prior Standard Deviation parameter "*nu),ylab=expression("Contamination parameter "*epsilon),
               lev=seq(-3.5,3.5,0.1),
               main="",plot.axes={axis(1);axis(2);points(1,0.005,pch=16,cex=1.8,col=4)})
dev.off()

### End Perform Posterior SD Emulator Diagnostics using large 40x40 grid of MCMC runs done above ###			
#######################################################################################################################################
  


#######################################################################################################################################
### Extended Emulator to provide partial derivative information ###

simple_emulator_partial_derivatives <- function(x1,                    # matrix of design points: locations in input space where MCMC runs were performed
                                                y1,                    # vector of run outputs (e.g. posterior means or posterior sds) corresponding to x1
                                                xp,                    # prediction points where the emulator will be evaluated
                                                theta=1,               # emulator correlation length parameter
                                                nugget=0.1,            # emulator nugget parameter
                                                xran,                  # physical range of X space, used to convert to [-1,1] scale
                                                just_var=1,            # set to 1 if just want emulator means and variances, set to 0 if want means and full covariance matrix across xp
                                                em_mcmc_notf=0,        # set to 0 to get predictions for f, set to 1 to get predictions for noisy MCMC output, e.g. for diagnostics
                                                partial_deriv_index=1  # which partial derivative to evaluate, set to 1 for x_1, 2 for x_2 etc.
                                                ){
  x1 <- scale(x1,center=apply(xran,1,mean),scale=apply(xran,1,diff)/2)     # scale design inputs to standard [-1,1] scale
  xp <- scale(xp,center=apply(xran,1,mean),scale=apply(xran,1,diff)/2)     # scale new prediction point inputs to standard [-1,1] scale
  D <- y1                             # store design outputs as data vector D
  mu <- mean(y1)                      # estimate mean as mean of runs XXX JP: is it simpler just to standardise the runs? Then specify mean 0, var 1 for emulator?
  sigma <- sd(y1)                     # estimate sigma as sd of runs
  Cmat <- function(d) sigma^2 * exp(-d^2/theta^2) * (1-nugget)    # define squared exponential covariance function: see equation (10)
  covBD <- Cmat( as.matrix(pdist(xp,x1)) )      # Covariance matrix between new outputs B = f(xp) and design points D = f(x1): see equation (11)
  varD  <- Cmat( as.matrix(dist(x1)) ) + sigma^2 * diag(nugget,nrow=nrow(x1))    # Variance matrix for design points D = f(x1), restoring nugget to diagonal: see equation (10)
  varB  <- Cmat( as.matrix(dist(xp)) ) + em_mcmc_notf * sigma^2 * diag(nugget,nrow=nrow(xp))		# Var for prediction points B = f(xp). Note: em_mcmc_notf = 1 => variance includes mcmc noise for diagnostics, em_mcmc_notf = 0 just emulating f

  ### partial derivatives covariances ###
  mdist1d <- function(x,y,i) outer(x[,i],y[,i],"-")			# 1d dist in ith dimension, including sign
  cov_dfdx_D    <- -2/theta^2 * mdist1d(xp,x1,partial_deriv_index) * covBD   # Cov of df/dx_i at prediction points xp, with elements of D = f(x1) (see equation (12))
  cov_dfdx_fx   <- -2/theta^2 * mdist1d(xp,xp,partial_deriv_index) * varB    # Cov of df/dx_i at prediction points xp, with elements of B = f(xp) at prediction points (see equation (12))
  cov_dfdx_dfdx <- ( 2/theta^2 - 4/theta^4 * mdist1d(xp,xp,partial_deriv_index)^2 ) * varB   # Cov of df/dx_i at prediction points xp with itself (generalisation of equation (12))
  
  ### Construct Cov matrix for Bfull = (B, dfdx) with the run outputs D. Note B = f(xp) ###
  cov_Bfull_D <- rbind( covBD, cov_dfdx_D )
  
  ### Construct Var matrix for Bfull = (B, dfdx). Note B = f(xp) ###
  Var_Bfull <- cbind( rbind( varB , cov_dfdx_fx) , rbind( t(cov_dfdx_fx) , cov_dfdx_dfdx ) )
  
  ### set prior expectations of Bfull = (B, dfdx) ###
  E_Bfull <- c( rep(mu,nrow(xp)), rep(0,nrow(xp)) )		# set prior expectation for Bfull: zero for the partial derivatives dfdx			
  ED 	  <-   mu                                       # set prior expectation for D
  
  ### Evaluate Emulator expecation and variance ###
  varD_inv <- chol2inv(chol(varD))                                   # invert varD using Choleski decompostion
  ED_Bfull   <- E_Bfull + cov_Bfull_D %*% varD_inv %*% (D - ED)      # Mean update: see equation (7)
  VarD_Bfull <- Var_Bfull - cov_Bfull_D %*% varD_inv %*% t(cov_Bfull_D)      # Covariance matrix update: see equation (8)

  ### Undo scaling effect on derivative terms in emualtor expectations and variances ###
  dim_scaling <- apply(xran,1,diff)/2     # the scale factor reduction in x inputs when converting to [-1,1] range
  ED_Bfull[(nrow(xp)+1):length(ED_Bfull)]    <- 1/dim_scaling[partial_deriv_index] * ED_Bfull[(nrow(xp)+1):length(ED_Bfull)]      # undo scaling on last half of expectation of Bfull = (B, dfdx)
  VarD_Bfull[,(nrow(xp)+1):length(ED_Bfull)] <- 1/dim_scaling[partial_deriv_index] * VarD_Bfull[,(nrow(xp)+1):length(ED_Bfull)]   # undo scaling on right half of var of Bfull = (B, dfdx)
  VarD_Bfull[(nrow(xp)+1):length(ED_Bfull),] <- 1/dim_scaling[partial_deriv_index] * VarD_Bfull[(nrow(xp)+1):length(ED_Bfull),]   # undo scaling on bottom half of var of Bfull = (B, dfdx)
  
  ### return either expectations and variances or expectations and full covariance matrix ###
  if(just_var)  return(list(ED_Bfull,diag(VarD_Bfull))) else    # return emulator means and variances at prediction points xp
  if(!just_var) return(list(ED_Bfull,VarD_Bfull))               # or return emulator means and full covariance matric for prediction points xp
}

### End Extended Emulator to provide partial derivative information ###
#######################################################################################################################################


#######################################################################################################################################
### Apply Partial Derivative Emulator to Case 1 in main paper.

### Matrix to store Table 2 contents ###
Table2_contents <- matrix(0,nrow=2,ncol=6)

### Set up Case 1: single point but partial derivatives requested ###
thet <- 0.7				         # correlation length
nug <- c(0.00002,0.004)    # nuggets
xp <- cbind(c(1.5,1.5),c(0.5,0.5))    # single input point in Case 1 in paper (emulator wants > 1 input point so double up)

### Apply the partial derivative emulator, loop over the two partial derivatives ###
for(out_mv in 1:2){     # loop over f1(x) = posterior mean and f2(x) = posterior SD
  for(k in 1:2){        # loop over two partial derivatives in nu and epsilon directions
    em_out <- simple_emulator_partial_derivatives(x1=x1,y1=y1[,out_mv],xp=xp,partial_deriv_index=k,theta=thet,nugget=nug[out_mv],xran=xran,just_var=1,em_mcmc_notf=0)
    ind_shift <- 3*(out_mv-1)
    Table2_contents[1,1+ind_shift] <- em_out[[1]][1,1]       # store emulator expectation of f(x) (same for each k)
    Table2_contents[2,1+ind_shift] <- sqrt(em_out[[2]])[1]   # store emulator sd of f(x) (same for each k)
    
    Table2_contents[1,1+k+ind_shift] <- em_out[[1]][3]       # store emulator expecation of df(x)/dx_k
    Table2_contents[2,1+k+ind_shift] <- sqrt(em_out[[2]][3]) # store emulator sd of df(x)/dx_k
  }
}
rownames(Table2_contents) <- c("E[.]","SD[.]")
colnames(Table2_contents) <- c(" $f_1(x)$","$\\frac{\\partial f_1(x)}{\\partial \\nu}$","$\\frac{\\partial f_1(x)}{\\partial \\epsilon}$",
                    "$f_2(x)$","$\\frac{\\partial f_2(x)}{\\partial \\nu}$","$\\frac{\\partial f_2(x)}{\\partial \\epsilon}$")
Table2_contents      # The contents of Table 2 in main paper, giving partial derivatives for Case 1.

### Code to produce latex table ###
xtab <- xtable(Table2_contents, caption="\\footnotesize{Results of the local sensitivity analysis corresponding to specification case 1. The first row gives the expectation of all requested quantities of interest, namely the posterior mean $f_1(x)$, posterior SD $f_2(x)$ and partial derivatives of each. The second row gives the  corresponding SD of each of these estimates, which could be reduced using further MCMC runs.}",label="tab_case1")
align(xtab) <- "|r|rrr|rrr|"
digits(xtab) <- 3
print(xtab,file=paste(plot_dir,"Table2.tex",sep=""),sanitize.text.function = function(x) {x},floating=TRUE,latex.environments="center")

### End Apply Partial Derivative Emulator to Case 1 in main paper.
#######################################################################################################################################











