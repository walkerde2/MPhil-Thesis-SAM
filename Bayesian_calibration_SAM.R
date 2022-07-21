###########################################
# proxy to proxy calibration of IHG vs SAM
###########################################
require(astrochron)
install.packages("BayesianTools")
require(BayesianTools)
install.packages("truncdist")
require(truncdist)
install.packages("zoo")
require(zoo)

# set working directory
setwd("C://Dissertation R//Calibration")

#####################
# upload corrected IHG
IHG_corr <- read.table(file="C://Dissertation R//R_Data//myIHG_corr.txt", header=TRUE)
lag <- 20 # set a lag time (in years) for the IHGv
# apply lag
IHG_corr[,1] <- IHG_corr[,1] - lag
#####################
# upload reconstructed SAM data (Dätwyler et al., 2020)
SAM <- read.table(file="C://Dissertation R//Calibration//SAM Bayesian Calibration//SAM_last1k_Datwyler2020.txt", header = TRUE)
# estimate mean and error across SAM ensembles
SAM_mean <- apply(SAM[,6:1005], 1, mean)
SAM_sd <- apply(SAM[,6:1005], 1, sd)
# bind data
SAM <- cbind(SAM$Year, SAM_mean, SAM_sd)
#####################
# or use SAM from Dalaiden et al., 2022
SAM_dal <- read.table(file="C://Dissertation R//Calibration//SAM Bayesian Calibration//SAM_DA_Dalaiden.txt", header = TRUE)
SAM_dal[,2] <- SAM_dal[,3]
SAM_dal[,3] <- 1
####################
# upload observational SAM index (Marshall, 2003)
sam_index <- read.table(file="C://Dissertation R//Calibration//SAM Bayesian Calibration//SAM_index_Marshall.txt", header=TRUE)
#####################
# upload volcanic reconstruction from Sigl. et al., 2015
GVF <- read.table(file="C://Dissertation R//Calibration//SAM Bayesian Calibration//Sigl2015_GVF.txt", header=TRUE)
#####################

# interpolate and shift data
# let's consider a calibration period from 1950 to a selected cutoff age
# according to Datwyler et al. the RE values show skill from around AD 1400 onwards,
# then the reconstruction suffers from a reduced number of available predictor datasets
# focus on 1400 - 1950 (apply cut off of 1400)
# subset both IHG and SAM to 1400 onwards

cutoff <- 1400
SAM_calib <- SAM[SAM[,1] >= cutoff & SAM[,1] <= 1950 - lag,]
IHG_calib <- IHG_corr[IHG_corr[,1] >= cutoff & IHG_corr[,1] <= 1950 - lag,]
# reconstruction period 
IHG_rec <- IHG_corr[IHG_corr[,1] < cutoff,]

####################################################################
# Bayesian linear regression (with or without bandpassing)
##### bandpass the target
bpass <- 0.1 # 1/years
y <- bandpass(cbind(SAM_calib[,1],  SAM_calib[,2]), fhigh=bpass, genplot=FALSE, verbose = FALSE)[,2]
y_err <- bandpass(cbind(SAM_calib[,1],  SAM_calib[,3]), fhigh=bpass, genplot=FALSE, verbose = FALSE)[,2]
# x <- bandpass(cbind(IHG_calib[,1], IHG_calib[,2]), fhigh=bpass,genplot=FALSE,verbose = FALSE)[,2]
# x_err <- bandpass(cbind(IHG_calib[,1], IHG_calib[,3]), fhigh=bpass, genplot=FALSE, verbose = FALSE)[,2]
# bandpass the whole target record
# smooth the sAM record - remove all variability (noise) shorter than a decade
SAM_bp <- bandpass(cbind(SAM[,1],  SAM[,2]), fhigh=bpass, genplot=FALSE, verbose = FALSE)
SAM_err_bp <- bandpass(cbind(SAM[,1],  SAM[,3]), fhigh=bpass, genplot=FALSE, verbose = FALSE)
#
# independent variable
x <- IHG_calib[,2]
x_err <- IHG_calib[,3]

# all of this prepares x and y values for linear regression
# linear model
lm <- lm(y ~ x)
summary(lm) # information about the linear regression

# set up a MCMC method for linear regression
trueA <- lm$coefficients[2] # slope1
trueB <- lm$coefficients[1] # intercept
trueSd <- sqrt(deviance(lm)/df.residual(lm)) # residual std error

# likelihood function
likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  # linear model
  pred = (a*x) + b + rnorm(1, mean=0, sd = sd)
  # predicted values should be close to observed values
  likelihood1 = dnorm(pred, mean = y, sd = y_err, log=TRUE) 
  # variance of predicted values should be close to observed values
  V <- var(pred) - var(y)
  likelihood2 = dnorm(V, mean = 0, sd = 1, log = T) 
  # sum the two likelihood indeces
  sumll = sum(likelihood1 + likelihood2)
  return(sumll)   
}

###############
# DE Sampler #
##############
ll <- likelihood
startvalue = c(trueA, trueB, trueSd)
# set priors for slope, intercept and model error
low <- c(-2,-7.5,0) 
up <- c(0,0,5)
# setup
bayesianSetup <- createBayesianSetup(likelihood = ll, 
                                     lower = low, upper = up,
                                     best = startvalue, catchDuplicates = FALSE)
settings = list(iterations = 100000) # monte carlo to propagate
# run chains
chain <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DE",settings = settings)
# the following line extends extends the chain from above
# chain <- runMCMC(bayesianSetup = chain, sampler = "DEzs",settings = settings)
#############################
# summary stats and plot chain
marginalPlot(chain)
# summary(chain)
# plot(chain)
c = getSample(chain, start = 0)
# ################
# burnIn = 90000 # burn in
# c <- c[-(1:burnIn),] # remove burn in steps

# reconstruct propagating IHG data error
recon <- sapply(1:dim(c)[1], function(i) {
  X <- rnorm(length(x), mean=x, sd=x_err)
  result <- (c[i,1]*X) + c[i,2]
})

# estimate confidence levels and wmean
quant <- t(apply(recon,1,quantile, c(0.5, 0.05, 0.95, 0.32, 0.68)))

#############################################################
# plot Datwyler vs reconstructed SAM over calibration period (1400 - 1950)

{
  par(mfrow=c(1,1))
  plot(SAM_calib[,1], y, type='n', col=1, xlim=c(cutoff,1950), ylim=c(6,-8),
       xlab="Years AD", ylab="SAM index", main="Reconstructed SAM record vs Datwyler et al. (2020)
       (Over the Calibration Period)")
  lines(SAM_calib[,1], y, lwd=2)
  abline(h=seq(-10,10,1), lwd=0.5, lty=3)
  polygon(c(SAM_calib[,1], rev(SAM_calib[,1])), c(y + 2*y_err, rev(y - 2*y_err)),
          border=NA, col=adjustcolor("grey", alpha=0.2))
  polygon(c(IHG_calib[,1], rev(IHG_calib[,1])), c(quant[,2], rev(quant[,3])),
          border=NA, col=adjustcolor("red", alpha=0.2))
  lines(IHG_calib[,1], quant[,"50%"], col=2, lwd=2)
  legend("topright", c("Datwyler et al. (2020)","Reconstructed SAM"),
         lty=c(1,1,1,1), col=c(1,2), cex=0.8)
}
# plot the same but with Dalaiden and Marshall too (observed)

Marshall_bp <- bandpass(cbind(sam_index[,1],  sam_index[,2]), fhigh=bpass, genplot=FALSE, verbose = FALSE)
Dal_bp <- bandpass(cbind(SAM_dal[,1],  SAM_dal[,2]), fhigh=bpass, genplot=FALSE, verbose = FALSE)

{
  par(mfrow=c(1,1))
  plot(SAM_calib[,1], y, type='n', col=1, xlim=c(cutoff,2000), ylim=c(6,-8),
     xlab="Years AD", ylab="SAM index",
     main="SAM records over the Calibration Period")
  lines(SAM_calib[,1], y, lwd=2)
  abline(h=seq(-10,10,1), lwd=0.5, lty=3)
  polygon(c(SAM_calib[,1], rev(SAM_calib[,1])), c(y + 2*y_err, rev(y - 2*y_err)),
        border=NA, col=adjustcolor("grey", alpha=0.2))
  polygon(c(IHG_calib[,1], rev(IHG_calib[,1])), c(quant[,2], rev(quant[,3])),
        border=NA, col=adjustcolor("red", alpha=0.2))
  lines(IHG_calib[,1], quant[,"50%"], col=2, lwd=2)
  # show data assimilated SAM from Dalaiden et al., 2020
  #polygon(c(SAM_dal[,1], rev(SAM_dal[,1])), c(SAM_dal[,2]+2*SAM_dal[,3], rev(SAM_dal[,2]-2*SAM_dal[,3])),
  #        border=NA, col=adjustcolor("orange", alpha=0.2))
  #lines(SAM_dal[,1], SAM_dal[,2], col="mediumpurple3", lwd=1)
  lines(Dal_bp[,1], Dal_bp[,2], col="mediumpurple3", lwd=2)
  # show observational annual SAM index
  #lines(sam_index[,1], sam_index[,2], col="dodgerblue", lwd=1)
  lines(Marshall_bp[,1], Marshall_bp[,2], col="dodgerblue", lwd=2)
  legend("topright", c("Datwyler et al. (2020)","Dalaiden et al. (2022)",
                       "Reconstructed", "Observed (Marshall et al. 2003)"),
       lty=c(1,1,1,1), col=c(1,"mediumpurple3",2,"dodgerblue"), cex=0.8)
}
############################################################
# # reconstruct past 2kyr
# recon_past <- sapply(1:dim(c)[1], function(i) {
#   result <- (c[i,1]*IHG_rec[,2]) + c[i,2]
# })
# reconstruct propagating IHG data error
recon_past <- sapply(1:dim(c)[1], function(i) {
  X <- rnorm(length(IHG_rec[,2]), mean=IHG_rec[,2], sd=IHG_rec[,3])
  result <- (c[i,1]*X) + c[i,2]
})

# estimate confidence levels and wmean
quant_past <- t(apply(recon_past,1,quantile, c(0.5, 0.05, 0.95, 0.32, 0.68)))
# # find posterior weighted mean
# wmean_past <- sapply(1:nrow(recon_past), function(i){ 
#   densit <- density(recon_past[i,])
#   weighted.mean(densit$x,densit$y)
# })
#############################################################
 #reconstruction plot

{
  par(mfrow=c(1,1))
  plot(SAM_bp[,1], SAM_bp[,2], type='n', col='grey40', xlim=c(0,2000), ylim=c(7,-10), lwd=0.5,
     xlab="Years AD", ylab="SAM index",
     main="SAM Reconstruction using IHG over the past 2,000 years
     (with Common Era Volcanism)")
  polygon(c(SAM_bp[,1], rev(SAM_bp[,1])), c(SAM_bp[,2]+ 2*SAM_err_bp[,2], rev(SAM_bp[,2] - 2*SAM_err_bp[,2])),
        border=NA, col=adjustcolor("grey20", alpha=0.2))
  # polygon(c(SAM_bp[,1], rev(SAM_bp[,1])), c(SAM_bp[,2]+ SAM_err_bp[,2], rev(SAM_bp[,2] - SAM_err_bp[,2])),
  #         border=NA, col=adjustcolor("grey20", alpha=0.2))
  lines(SAM_bp[,1], SAM_bp[,2], lwd=2, col='grey40')
  rect(900,7,1399,-8, col = "white", border = NA) # hide SAM record older than 1400 AD
  polygon(c(IHG_calib[,1], rev(IHG_calib[,1])), c(quant[,2], rev(quant[,3])),
        border=NA, col=adjustcolor("red", alpha=0.2))
  # polygon(c(IHG_calib[,1], rev(IHG_calib[,1])), c(quant[,4], rev(quant[,5])),
  #         border=NA, col=adjustcolor("darkred", alpha=0.2))
  lines(IHG_calib[,1], quant[,"50%"], col=2, lwd=2)
  polygon(c(IHG_rec[,1], rev(IHG_rec[,1])), c(quant_past[,2], rev(quant_past[,3])),
        border=NA, col=adjustcolor("orange", alpha=0.2))
  # polygon(c(IHG_rec[,1], rev(IHG_rec[,1])), c(quant_past[,4], rev(quant_past[,5])),
  #         border=NA, col=adjustcolor("darkorange", alpha=0.2))
  lines(IHG_rec[,1], quant_past[,"50%"], col="orange", lwd=2)
  abline(h=seq(-10,10,1), lwd=0.25, lty=3)
  abline(v=seq(-200,2100,100), lwd=0.25, lty=3)
  abline(v=1400, lty=2, lwd=1.5)
  abline(h=0, lty=2, lwd=1.5)
  # show observational annual SAM index
  #lines(sam_index[,1], sam_index[,2], col="dodgerblue", lwd=0.5)
  lines(Marshall_bp[,1], Marshall_bp[,2], col="dodgerblue", lwd=2)
  legend("bottomleft", c("Datwyler et al. (2020)", "Reconstructed (calibration interval)", "Reconstructed", "Observed (Marshall et al. 2003)"),
       lty=c(1,1,1,1), col=c('grey40',2,"orange","dodgerblue"), cex=0.8)
}

###########
# writeup #
write.table(IHG_calib,file="mySAM_calinterval.txt", row.names=FALSE,
            col.names=c("yrsAD","median","stdev"))
write.table(IHG_rec,file="mySAM_reconstruction.txt", row.names=FALSE,
            col.names=c("yrsAD","median","stdev"))
# File with total reconstruction (2000 years)
fullSAM_recon <- read.table("fullSAM_reconstruction.txt", header=TRUE)

###########

# plot alongside volcanic forcing
par(new=TRUE)
plot(GVF, type='h', xlim=c(0,2000), col='mediumpurple3', xaxt='n', yaxt='n', ylab="", xlab="")
lines(smooth.spline(GVF, spar=0.25), lwd=2, col='mediumpurple3')

#############################################################
# estimate two-tailed Student t test on median SAM
win=100 # as in Ortega et al., 2015
z <- rollapplyr(c(quant_past[,"50%"],quant[,"50%"]), win, function(x) t.test(x, alternative = c("two.sided"))$p.value)
t <- rollapplyr(c(quant_past[,"50%"],quant[,"50%"]), win, function(x) t.test(x, alternative = c("two.sided"))$estimate)
tim <- seq(from=IHG_rec[1,1]+(win/2)-1,to=IHG_calib[length(IHG_calib[,1]),1]-(win/2))
z <- cbind(tim,z,t)
# convert positive and negative values to -1/+1
z[z[,2]>0.05,] <- NA
z <- na.omit(z)
z[z[,3]<0,3] <- -1
z[z[,3]>0,3] <- 1
z <- z[,-2]
par(new=TRUE)
plot(z,type="p",xlim=c(0,2000), ylim=c(1.25,-1.25), cex=0.5, pch=15,
     xaxt='n', yaxt='n', ylab="", xlab="")
#############################################################

#############################################################
# create box plots using error propagation
# thin 'recon_past' to 1,000 ensembles
thin <- 150
myrecon2 <- recon_past[, seq(1, ncol(recon_past), thin)] # extract only every nth output
myrecon1 <- recon[, seq(1, ncol(recon_past), thin)]
myrecon <- rbind(myrecon2, myrecon1)
myrecon <- cbind(c(IHG_rec[,1],IHG_calib[,1]), myrecon)
# calculate distributions for different climate periods
RWP <- myrecon[myrecon[,1] < 540, -1]
LALIA <- myrecon[myrecon[,1] >= 540 & myrecon[,1] <= 900, -1]
MCA <- myrecon[myrecon[,1] > 900 & myrecon[,1] <= 1250, -1]
LIA <- myrecon[myrecon[,1] > 1250 & myrecon[,1] <= 1850, -1]
# generate boxplot
boxplot(as.vector(RWP), as.vector(LALIA), as.vector(MCA), as.vector(LIA), sam_index[,2],
        ylim=c(10,-10), at=seq(1:5), ylab="SAM index", col=c("orange", "skyblue", "red", "dodgerblue", "grey"), 
        names = c("RWP", "LALIA", "MCA", "LIA", "post-1850"), outline=FALSE, width=c(rep(0.2,5)))
abline(h=0,lty=3)
#############################################################

#############################################################
# calculate rates of change over rolling windows
dt <- 100 # window length in years
mytid <- seq(1, nrow(myrecon)-dt,by=10)
# calculate rates of change for each ensemble memeber
mytrend <- sapply(mytid, function(x){
  est <- sapply(2:ncol(myrecon), function(j){
    mydifff <- lm(myrecon[x:(x+dt),1] ~ myrecon[x:(x+dt),j])$coefficients[2]
  })
})
# estimate mean, and Cl
mytrend2 <- t(apply(mytrend, 2,quantile, probs =c(0.05, 0.32, 0.5, 0.68, 0.95), na.rm=TRUE))
trend <- dt*lm(sam_index[,2] ~ sam_index[,1])$coefficients[2] # calculate modern trend
# plot results
plot(mytid+(dt/2), mytrend2[,3], type='l', ylim=c(-25,25), xlim=c(0,2000),
     xlab="years AD", ylab="RoC (SAM/dt)")
lines(mytid+(dt/2), mytrend2[,1], col="grey")
lines(mytid+(dt/2), mytrend2[,5], col="grey")
abline(h=trend, col=2, lty=2)
abline(h=0, col=1, lty=2)
abline(v = seq(2000, 0, -50), lty = 3, lwd= 0.5)
#############################################################

# PAGES2K Temperature Record
setwd("C://Dissertation R//Calibration//Pages 2K//recons//recons")
PAGES2K_recon <- read.delim("t_recon.txt")
bpass <- 0.1 # 1/years
PAGES2K_bp <- bandpass(cbind(PAGES2K_recon[,1],  PAGES2K_recon[,3]), fhigh=bpass, genplot=F, verbose=F)

# Change Column Names
names(PAGES2K_recon)[1]<-("Year AD")
names(PAGES2K_recon)[2]<-("Raw Instrumental Data")
names(PAGES2K_recon)[3]<-("50th percentile")
names(PAGES2K_recon)[4]<-("2.5th percentile")
names(PAGES2K_recon)[5]<-("97.5th percentile")

y <- PAGES2K_recon

# plot PAGES2K
par(new=TRUE)
plot(PAGES2K_recon[,1], type='n', col=2, xlim=c(0,2000), ylim=c(-0.8,0.6),
     xlab="Year", ylab="PAGES2K Temperature Anomaly", main="PAGES2K Ensemble Temperature Reconstruction")
lines(PAGES2K_recon[,1], PAGES2K_recon[,3], type='l',
      col=2, lwd=1, xlim=c(0,2000))
abline(h = seq(-0.8, 0.6, 0.1), lty = "dotted", lwd= 0.5)
abline(v=1400, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)

# Plot PAGES2K with SAM reconstruction and volcanism
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM_bp[,1], SAM_bp[,2], type='n', col=1, xlim=c(0,2000), ylim=c(10,-10), lwd=0.5,
            xlab="Years AD", ylab="SAM index")
polygon(c(SAM_bp[,1], rev(SAM_bp[,1])), c(SAM_bp[,2]+ 2*SAM_err_bp[,2], rev(SAM_bp[,2] - 2*SAM_err_bp[,2])),
              border=NA, col=adjustcolor("grey20", alpha=0.2))
lines(SAM_bp[,1], SAM_bp[,2])
rect(900,7,1399,-8, col = "white", border = NA) # hide SAM record older than 1400 AD
polygon(c(IHG_rec[,1], rev(IHG_rec[,1])), c(quant_past[,2], rev(quant_past[,3])),
                  border=NA, col=adjustcolor("orange", alpha=0.2))
#(c(IHG_rec[,1], rev(IHG_rec[,1])), c(quant_past[,4], rev(quant_past[,5])),
  #               border=NA, col=adjustcolor("darkorange", alpha=0.2))
lines(IHG_rec[,1], quant_past[,"50%"], col="darkorange") 
par(new=TRUE)
plot(PAGES2K_recon[,1], PAGES2K_recon[,3], type='l',
     col='dodgerblue2', lwd=1, xlim=c(0,2000), ylim=c(-0.8,0.8), xlab="", ylab="",
     main="PAGES2K vs SAM Reconstruction", xaxt='n', yaxt='n' )
axis(4, col='dodgerblue2', ylab='PAGES 2K', seq(-0.5, 0.5, 0.25))
mtext("PAGES2K Temp Anomaly", side=4, line=3, srt=45)
# abline(h = seq(-0.8, 0.8, 0.2), lty = "dotted", lwd= 0.5)
abline(v=1400, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)
legend("bottomleft", c("Reconstructed SAM", "Reconstructed (Calibration Interval)", 
                       "PAGES 2K Temp Record (2019)", "Volcanism, Sigl et al. (2015)"),
                 lty=c(1,1,1), col=c("orange", 1,"dodgerblue2", "mediumpurple3"), cex=0.7)
par(new=TRUE)
plot(GVF, type='h', xlim=c(0,2000), col='mediumpurple3', xaxt='n', yaxt='n', ylab="", xlab="")
lines(smooth.spline(GVF, spar=0.25), lwd=2, col='mediumpurple3')

# need to adjust the polygons showing climate periods 
panel.first = rect(c(1,-0.8), -1e6, c(541,0.4), 1e6, col=rgb(0.5,0.5,0.5,1/4), border='tomato3') #MWP
panel.first = rect(c(536,-0.8), -1e6, c(660,0.4), 1e6, col=rgb(0.5,0.5,0.5,1/4), border='steelblue2') #LALIA
panel.first = rect(c(900,-0.8), -1e6, c(1250,0.4), 1e6, col=rgb(0.5,0.5,0.5,1/4), border='tomato3') #MCA
panel.first = rect(c(1250,-0.8), -1e6, c(1850,0.4), 1e6, col=rgb(0.5,0.5,0.5,1/4), border='steelblue2') #LIA

#Last Millennium Reanalysis Temperature Record 
library(ncdf4)
setwd("C://Dissertation R//Calibration//Pages 2K")
global_temp <- "gmt_MCruns_ensemble_full_LMRv2.0.nc"
time <- ncvar_get(global_temp, "time")
gmt <- ncvar_get(global_temp, "Global-Mean Air Temperature at Surface")
gm_temp <- nc_open(global_temp)
print(gm_temp)

NH_temp <- "nhmt_MCruns_ensemble_full_LMRv2.0.nc"
nht <- nc_open(NH_temp)
print(nht)

SH_temp <- "shmt_MCruns_ensemble_full_LMRv2.0.nc"
sht <- nc_open(SH_temp)
print(sht)

#############################################################

# volcanism stacks using fullSAM_recon

# create cluster subsets of SAM record - 150 years either side of eruption year zer0
# 540 CE eruption - year zero: 539 
clusterA <- fullSAM_recon[fullSAM_recon[,1] > 438,]
cluster1 <- clusterA[clusterA[,1] < 640,]
cluster1[,2]=(cluster1[,2] - mean(cluster1[,2]))/sd(cluster1[,2])
cluster540 <- cbind(mytime, cluster1[,2])

# 682 CE unknown eruption - year zero: 681
clusterB <- fullSAM_recon[fullSAM_recon[,1] > 580,]
cluster2 <- clusterB[clusterB[,1] < 782,]
cluster2[,2]=(cluster2[,2] - mean(cluster2[,2]))/sd(cluster2[,2])
cluster682 <- cbind(mytime, cluster2[,2])

# 1258 CE Samalas - year zero: 1257
clusterC <- fullSAM_recon[fullSAM_recon[,1] > 1156,]
cluster3 <- clusterC[clusterC[,1] < 1358,]
cluster3[,2]=(cluster3[,2] - mean(cluster3[,2]))/sd(cluster3[,2])
cluster1258 <- cbind(mytime, cluster3[,2])

# 1458 CE unknown (Kuwae?) - year zero: 1457 
clusterD <- fullSAM_recon[fullSAM_recon[,1] > 1356,]
cluster4 <- clusterD[clusterD[,1] < 1558,]
cluster4[,2]=(cluster4[,2] - mean(cluster4[,2]))/sd(cluster4[,2])
cluster1458 <- cbind(mytime, cluster4[,2])

# 1815 CE Tambora - year zero: 1814 
clusterE <- fullSAM_recon[fullSAM_recon[,1] > 1713,]
cluster5 <- clusterE[clusterE[,1] < 1915,]
cluster5[,2]=(cluster5[,2] - mean(cluster5[,2]))/sd(cluster5[,2])
cluster1815 <- cbind(mytime, cluster5[,2])

#all clusters
mytime <- seq(-100,100)
all.clusters <- cbind(mytime, cluster1[,2], cluster2[,2], cluster3[,2], cluster4[,2], cluster5[,2])
volc.stack <- rowMeans(all.clusters[,2:6])
volc.sd <- apply(all.clusters[,2:6], 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-2,2), ylab = "SAM index", lwd=2,
     xlab="Year AD (relative to eruption year zero", 
     main="SAM Response to Volcanic Forcing (All Clusters)")
lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
# lines(mytime, cluster540[,2], col=2)
# lines(mytime, cluster682[,2], col=3)
# lines(mytime, cluster1258[,2], col=4)
# lines(mytime, cluster1458[,2], col=5)
# lines(mytime, cluster1815[,2], col=6)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)
axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))
# legend("bottomleft", c("Stacked Data", "540 CE", 
                       "682 CE", "1258 CE", "1458CE", "1815 CE"),
       lty=c(1,1,1), col=c(1,2,3,4,5,6), cex=0.7)

# 540 CE and 1258 CE stack
clusters_540_1258 <- cbind(cluster1[,2], cluster3[,2])
volc.stack <- rowMeans(clusters_540_1258)
volc.sd <- apply(clusters_540_1258, 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5), lwd=2,
     xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
     main="SAM response to volcanic forcing"))
lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)

#540 CE and 682 CE stack
clusters_540_682 <- cbind(cluster1[,2], cluster2[,2])
volc.stack <- rowMeans(clusters_540_682)
volc.sd <- apply(clusters_540_682, 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5), lwd=2,
     xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
     main="SAM response to volcanic forcing")
#lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
#lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)
lines(mytime, cluster540[,2], col=2)
lines(mytime, cluster682[,2], col=3)
legend("bottomleft", c("Stacked Data", "540 CE", 
                       "682 CE"),
       lty=c(1,1,1), col=c(1, 2, 3), cex=0.7)
axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))

# 540 CE, 1258 CE, and 1458 CE stack - counteracts the negative reponse found without 1458 included 
clusters_540_1258_1458 <- cbind(cluster1[,2], cluster3[,2], cluster4[,2])
volc.stack <- rowMeans(clusters_540_1258_1458)
volc.sd <- apply(clusters_540_1258_1458, 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5))
lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)

# 540 CE, 1258 CE, and 1815 CE stack - again doesn't show a phat lot
clusters_540_1258_1815 <- cbind(cluster1[,2], cluster3[,2], cluster5[,2])
volc.stack <- rowMeans(clusters_540_1258_1815)
volc.sd <- apply(clusters_540_1258_1815, 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5))
lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)

# 540 CE, 1258 CE, and 682 CE stack - supports the pattern found with 540 and 1258
clusters_540_1258_682 <- cbind(cluster1[,2], cluster3[,2], cluster2[,2])
volc.stack <- rowMeans(clusters_540_1258_682)
volc.sd <- apply(clusters_540_1258_682, 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5),lwd=2,
     xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
     main="SAM response to volcanic forcing")
#lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
#lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
lines(mytime, cluster540[,2], col=2)
lines(mytime, cluster682[,2], col=3)
lines(mytime, cluster1258[,2], col=4)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)
legend("bottomleft", c("Stacked Data", "540 CE", 
                       "682 CE", "1258 CE"),
       lty=c(1,1,1), col=c(1, 2, 3, 4), cex=0.7)
axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))

# what about just 1458 CE and 1815 CE? - completely opposite signal to the other clusters
# send SAM into a positive phase (whereas the others become negative)
clusters_1458_1815 <- cbind(cluster4[,2], cluster5[,2])
volc.stack <- rowMeans(clusters_1458_1815)
volc.sd <- apply(clusters_1458_1815, 1, sd)
plot(mytime, volc.stack, type='l', ylim=c(-3,3), lwd=2,
     xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
     main="SAM response to volcanic forcing")
#lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
#lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
lines(mytime, cluster1458[,2], col=5)
lines(mytime, cluster1815[,2], col=6)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)
legend("bottomleft", c("Stacked Data", "1458 CE", 
                       "1815 CE"),
       lty=c(1,1,1), col=c(1,5,6), cex=0.7)
axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))



# Calculate the median of all these subsets
# so I want the average (median) values of all 4 cluster subsets combines
# this will be a median value for each year

median(cluster1[,2])

# export to tables to calculate median in excel 
setwd("C://Dissertation R//Calibration//Volcanism")
write.table(cluster1,file="cluster_1.txt", row.names=FALSE,
            col.names=c("yrsAD","median","stdev"))
write.table(cluster2,file="cluster_2.txt", row.names=FALSE,
            col.names=c("yrsAD","median","stdev"))
write.table(cluster3,file="cluster_3.txt", row.names=FALSE,
            col.names=c("yrsAD","median","stdev"))
write.table(cluster4,file="cluster_4.txt", row.names=FALSE,
            col.names=c("yrsAD","median","stdev"))
write.table(cluster5,file="cluster_5.txt", row.names=FALSE,
            col.names=c("yrsAD", "median", "stdev"))

# Upload All Clusters - the median, mean and stdev of each year across each subset

# Import Stackdata

all_clusters <- stackdata_500_1107_1207_1407
min(all_clusters[,3])
max(all_clusters[,3])
mean(all_clusters[,3])
View(all_clusters)

# Plot clusters as a stack - using mean
# stackdata_clusters1235 refers to
# 4 of the volcanic clusters, specifically:
# (year zero) 500 AD, 1107 AD, 1207 AD, 1407 AD
par(new=TRUE)
plot(all_clusters[,1], all_clusters[,3], type='l', col=1, xlim=c(-80,80), ylim=c(-8,-2), lwd=0.5,
     xlab="Years AD (relative to year zero of eruption cluster)", ylab="SAM index (averaged over four time periods)",
     main="SAM index response to volcanic forcing")
abline(v=0, lty=2, lwd=1.5)
abline(h=-4.8, lty=2, lwd=1.5)
lines(all_clusters[,1], all_clusters[,7], lty="dotted", col="grey20")
lines(all_clusters[,1], all_clusters[,7], lty="dotted", col="grey20")

# Plots all clusters separately
# Import all separate clusters
sep_clusters <- sep_clusters
plot(all_clusters[,1], all_clusters[,3], type='l', col=1, xlim=c(-200,200), ylim=c(-10,2), lwd=3,
     xlab="Years AD (relative to year zero of eruption cluster)", ylab="SAM index",
     main="SAM index response to volcanic forcing")
abline(v=0, lty=2, lwd=1.5)
abline(h=-4.8, lty=2, lwd=1.5)
lines(all_clusters[,1], sep_clusters[,2], col=2, lwd=2)
lines(all_clusters[,1], sep_clusters[,5], col="dodgerblue2", lwd=2)
lines(all_clusters[,1], sep_clusters[,8], col="mediumpurple", lwd=2)
lines(all_clusters[,1], sep_clusters[,11], col="limegreen", lwd=2)
lines(all_clusters[,1], sep_clusters[,14], col="orange", lwd=2)
legend("topright", c("Mean across clusters", "500 AD year zero", 
                       "1107 AD year zero", "1207 AD year zero", "1407 AD year zero", "1707 AD year zero"),
       lty=c(1,1,1), col=c(1, 2,"dodgerblue2", "mediumpurple", "limegreen", "orange"), cex=1)

# # Plot clusters excluding 1707 AD
# # upload stack data using import dataset in the global environment
# all_clusters_no1707 <- stackdata_500_1107_1207_1407
# plot(all_clusters[,1], all_clusters_no1707$mean.all, type='l', col=1, xlim=c(-100,100), ylim=c(-10,2), lwd=3,
#      xlab="Years AD (relative to year zero of eruption cluster)", ylab="SAM index",
#      main="SAM index response to volcanic forcing")
# abline(v=0, lty=2, lwd=1.5)
# abline(h=-4.8, lty=2, lwd=1.5)
# lines(all_clusters[,1], sep_clusters[,2], col=2, lwd=2)
# lines(all_clusters[,1], sep_clusters[,5], col="dodgerblue2", lwd=2)
# lines(all_clusters[,1], sep_clusters[,8], col="orange", lwd=2)
# lines(all_clusters[,1], sep_clusters[,11], col="limegreen", lwd=2)
# legend("topleft", c("Mean across clusters", "500 AD year zero", 
#                      "1107 AD year zero", "1207 AD year zero", "1407 AD year zero"),
#        lty=c(1,1,1), col=c(1, 2,"dodgerblue2", "orange", "limegreen"), cex=0.8)
# axis(side=1, at=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100))
# # lines(all_clusters[,1], all_clusters[,3], col="grey20", lwd=2) #average including 1707 AD cluster

#Plot clusters excluding 1707 AD and 1407 AD
# upload stack data using import dataset in the global environment
clusters_500_1107_1207 <- stackdata_500_1107_1207
plot(all_clusters[,1], clusters_500_1107_1207$mean.all, type='l', col=1, xlim=c(-100,100), ylim=c(-10,2), lwd=3,
     xlab="Years AD (relative to year zero of eruption cluster)", ylab="SAM index",
     main="SAM index response to volcanic forcing")
abline(v=0, lty=2, lwd=1.5)
abline(h=-4.8, lty=2, lwd=1.5)
lines(all_clusters[,1], sep_clusters[,2], col=2, lwd=2)
lines(all_clusters[,1], sep_clusters[,5], col="dodgerblue2", lwd=2)
lines(all_clusters[,1], sep_clusters[,8], col="orange", lwd=2)
legend("topleft", c("Mean across clusters", "500 AD year zero", 
                    "1107 AD year zero", "1207 AD year zero"),
       lty=c(1,1,1), col=c(1, 2,"dodgerblue2", "orange"), cex=0.8)
axis(side=1, at=c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100))
# lines(all_clusters[,1], all_clusters[,3], col="grey20", lwd=2) #average including 1707 AD cluster
























#############################################################

ENSO Record 
setwd("C://Dissertation R//Calibration//ENSO")
library(ncdf4)
ENSO_PHYDA <- nc_open("PHYDA_AprMar.nc")
lon <- ncvar_get(ENSO_PHYDA, "lon")
lat <- ncvar_get(ENSO_PHYDA, "lat")
t <- ncvar_get(ENSO_PHYDA, "time")
head(lon)

ENSO.array <- ncvar_get(ENSO_PHYDA, "Nino 3.4 (reconstruction mean)")

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('double Nino_3.4_mn[tmon]',
                  list(dim1,dim2),
                  longname = 'Nino 3.4 (reconstruction mean)')

# Create an empty netcdf file...
outnc <- nc_create('ENSO_PHYDA.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, ENSO_PHYDA)
nc_close(outnc) # terminate access to the file




