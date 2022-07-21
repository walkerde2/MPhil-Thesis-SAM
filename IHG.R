# set working directory
setwd("C://Dissertation R//IHG")

install.packages("dplR", dependencies=TRUE, type="win.binary")
install.packages('astrochron')
install.packages('Hmisc')

require(astrochron)
require(Hmisc)
require(dplR)
# first let's upload our intcal20 data
NH <- read.table("C://Dissertation R//R_Data//IntCal20copy.txt", header = TRUE)
SH <- read.table("C://Dissertation R//R_Data//SHCal20copy.txt", header = TRUE)

# let's convert BP to AD
NH[,1] = 1950 - NH[,1]
names(NH)[1]<-("Year (AD)")
SH[,1] = 1950 - SH[,1]
names(SH)[1]<-("Year (AD)")

# let's subset our data
NH_2k = NH[NH[,1] > 0,]
SH_2k = SH[SH[,1] > 0,]

# let's plot the data with uncertainty
plot(NH_2k[,1], NH_2k[,4], type='l', col = "blue",
     xlim=c(2000, 0), ylim = c(-40, 30), xlab = "years AD", ylab = "14C (per mil)")
arrows(NH_2k[,1], NH_2k[,4] + NH_2k[,5], NH_2k[,1], NH_2k[,4] - NH_2k[,5], code = 0, col = "blue")
lines(SH_2k[,1], SH_2k[,4], col= "red")
arrows(SH_2k[,1], SH_2k[,4] + SH_2k[,5], SH_2k[,1], SH_2k[,4] - SH_2k[,5], code = 0, col = "red")
abline(v = seq(2000, 0, -100), lty = 3, lwd= 0.5)
legend("topright", c("NH", "SH"), lty = c(1,1), col = c("blue", "red"), cex= 0.5)

# let' calculate the IHG using Monte Carlo propagation of error
IHG = sapply(1:10000, function(x){
  NH_sim =rnorm(n = dim(NH_2k)[1], mean = NH_2k[,4], sd = NH_2k[,5])
  SH_sim =rnorm(n = dim(SH_2k)[1], mean = SH_2k[,4], sd = SH_2k[,5])
  IHG = SH_sim - NH_sim
})

IHG_quant = t(apply(IHG, 1, quantile, probs = c(0.05, 0.5, 0.95)))
# assuming the error is normally distributed after a few thousand permutations
IHG_sd <- apply(IHG, 1, sd)

#######################################################
# let's plot the IHG alongside selected climate indeces
# upload SAM data (Datwyler2020)
SAM <- read.table(file="C://Dissertation R//R_Data//SAM2020.txt", header = TRUE)
# bandpass data to filter high-frequency variability
SAM_bp <- bandpass(cbind(SAM[,1], SAM[,2]), 
                   fhigh=0.04, genplot=FALSE, verbose=FALSE)

######
# Plot NH and SH C14 record
lag <- 20 # set a lag time (in years) for the IHG
par(mar=c(5.1, 4.1, 4.1, 5.1)) # change margin figure to accommodate second axis label
plot(NH_2k[,1] - lag, NH_2k[,4], type='l', lwd=1,
     xlab = "Years AD", ylab = "C14 (per mil)", 
     xlim=c(0, 2000), ylim=c(-30,20), col='orangered3',
     main="Radiocarbon record over the Common Era")
lines(SH_2k[,1], SH_2k[,4], lwd=1, col='steelblue3')
#lines(NH_2k[,1], IHG_quant[,3], lty=3, lwd=0.5)
#abline(v = seq(2000, 0, -100), lty = 3, lwd= 0.5)
abline(v = 1850, lty = 3, lwd= 2.5, col=2) # PI boundary
abline(h=0, lty=3, lwd=2, col=1)
par(new=TRUE)
max_a <- pmax(IHG_quant[,1])
max_b <- pmax(IHG_quant[,3])
# add IHG alongside
lines(NH_2k[,1], IHG_quant[,2], type='l', lwd=1, col=1)
polygon(c(NH_2k$`Year (AD)`-lag, rev(NH_2k$`Year (AD)`-lag)), 
        c(max_a ,rev(max_b)), col = rgb(.3, 0.8, 1, 0.2), border=NA )
axis(side=1, at=c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000))
legend("topleft", c("NH", "SH", "IHG"), lty = c(1,1), 
       col = c("orangered3", "steelblue3", 1), cex= 0.8)

# Plot IHG (NB the IHG likely has a lagged response to some of these indeces due to the longer response time of the SOcean to atmopsheric forcing)
lag <- 20 # set a lag time (in years) for the IHG
par(mar=c(5.1, 4.1, 4.1, 5.1)) # change margin figure to accommodate second axis label
plot(NH_2k[,1] - lag, IHG_quant[,2], type='l', lwd=1,
     xlab = "Years AD", ylab = "IHG (per mil)", 
     ylim=c(-14,6), xlim=c(0, 2000),
     main="Interhemispheric Gradient (IntCal20 minus SHCal20)")
#lines(NH_2k[,1], IHG_quant[,1], lty=3, lwd=0.5)
#lines(NH_2k[,1], IHG_quant[,3], lty=3, lwd=0.5)
#abline(v = seq(2000, 0, -100), lty = 3, lwd= 0.5)
abline(v = 1850, lty = 3, lwd= 2.5, col=2) # PI boundary
abline(h=0, lty=3, lwd=2, col=1)
par(new=TRUE)
# max_a <- pmax(IHG_quant[,1])
# max_b <- pmax(IHG_quant[,3])
polygon(c(NH_2k$`Year (AD)`-lag, rev(NH_2k$`Year (AD)`-lag)), 
        c(max_a ,rev(max_b)), col = rgb(.3, 0.8, 1, 0.2), border=NA )
axis(side=1, at=c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000))


#### SAM ####
par(new=TRUE)
plot(SAM[,1], SAM[,2], type='l', xlab="", ylab="",
     col='lightblue', lwd=0.25, xlim=c(2000,0),ylim=c(3,-8), xaxt='n', yaxt='n')
lines(SAM_bp, col=4, lwd=2) # show bandpassed data
axis(4, col=4)
mtext("SAM index", side=4, line=2)
legend("topright", c("IH 14C gradient", "SAM"), col=c(1,4), 
       lty=c(1,1), cex=0.7)
par(mar=c(5.1, 4.1, 4.1, 2.1)) # reset margins to default values
#############
#### AMO ####
par(new=TRUE)
plot(AMO[,2], AMO[,5], type='l', 
     col='springgreen3', lwd=0.5, xlim=c(2000,0), ylim=c(21.6,22.8), xlab="", ylab="", 
     xaxt='n', yaxt='n')
lines(AMO_bp, col='springgreen4', lwd=2) # show bandpassed data
axis(4, col='springgreen3')
mtext("AMO index", side=4, line=2)
legend("topright", c("IH 14C gradient", "AMO"), col=c(1,'springgreen3'), 
       lty=c(1,1), cex=0.7)
par(mar=c(5.1, 4.1, 4.1, 2.1)) # reset margins to default values
#############
#### AMOC ####
par(new=TRUE)
plot(AMOC[,1], AMOC[,2], type='l',
     col=2, lwd=1, xlim=c(2000,0), ylim=c(33,25), xlab="", ylab="", 
     xaxt='n', yaxt='n')
lines(AMOC2[,1]+15, AMOC2[,2]-3, col=2) # shallower core (shifted by 15 years)
axis(4, col=2)
mtext("AMOC", side=4, line=2)
legend("topright", c("IH 14C gradient", "AMOC (Thornalley et al., 2018)"), col=c(1,2), 
       lty=c(1,1), cex=0.7)
par(mar=c(5.1, 4.1, 4.1, 2.1)) # reset margins to default values
#############

#######################################################
# let's plot the IHG spectra alongside selected climate indicies spectra
IHG_spectral <- read.table(file="C:\\Dissertation R\\R_Data\\IHG.txt",header=TRUE)

# Change Column Names
names(IHG_spectral)[1]<-("Year")
names(IHG_spectral)[3]<-("Median")
names(IHG_spectral)[2]<-("0.05 quant")
names(IHG_spectral)[4]<-("0.95 quant")

SAM_1400_1950 <- cbind(SAM$Year, SAM$Final_DJF_SAM_reconstruction)
SAM_1400_1950 <- SAM_1400_1950[SAM_1400_1950[,1] > 1449,]
SAM_1400_1950 <- SAM_1400_1950[SAM_1400_1950[,1] < 1951,]
names(SAM_1400_1950)[1]<-("Year AD")
names(SAM_1400_1950)[2]<-("SAM")

# run the redfit() function
x <- IHG_spectral[,3]
t <- IHG_spectral[,1]
IHGspec <- redfit(x, t, nsim = 1000, tType="time") # nsim = number of simulate AR1 spectra to compute 
x <- SAM_1400_1950[,2]
t <- SAM_1400_1950[,1]
SAM_1400_1950_spec <- redfit(x, t, nsim = 1000, tType="time") 
x <- AMO_spectral[,2]
t <- AMO_spectral[,1]
AMO_obsv_spec <- redfit(x, t, nsim = 1000, tType="time") 
x <- SAM_obsv[,2]
t <- SAM_obsv[,1]
SAM_obsv_spec <- redfit(x, t, nsim = 1000, tType="time") 

# plot redfit results
# (note that we plot log_period against power [i.e. 1/T] use the reciprocal of frequency)
plot(1/IHGspec[["freq"]], IHGspec[["gxxc"]], type = "n", ylab = "Normalized Spectral Power", 
     xlab = "log Period (yr)", main="Power Spectrum IHG vs SAM", log="x",
     xlim=c(2,100), ylim=c(0,25))

# SAM vs IHG spectral plot
lines(1/IHGspec[["freq"]], IHGspec[["gxxc"]], col = "red") 
lines((1/IHGspec[["freq"]]), IHGspec[["ci95"]], col = "grey20")
lines(1/SAM_1400_1950_spec[["freq"]], SAM_1400_1950_spec[["gxxc"]], col = "royalblue2") 
legend("topleft", c("IHG", "SAM (1400-1950)", "95% Cl"), lwd = 2,
       col = c("red", "royalblue2", "grey20"),
       bg = "white", cex=0.75)

# AMO observations vs IHG spectral plot
lines(1/AMO_obsv_spec[["freq"]], AMO_obsv_spec[["gxxc"]], col = "springgreen3") 
legend("topleft", c("IHG", "SAM (1400-1950)", "95% Cl"), lwd = 2,
       col = c("red", "springgreen3", "grey20"),
       bg = "white", cex=0.75)

# SAM observations vs IHG spectral plot
lines(1/SAM_obsv_spec[["freq"]], SAM_obsv_spec[["gxxc"]], col = "purple1") 
legend("topleft", c("IHG", "SAM (1400-1950)", "95% Cl"), lwd = 2,
       col = c("red", "purple1", "grey20"),
       bg = "white", cex=0.75)


# Compare results to data assimilation products from LM reanalysis
# open netcdf data
# name of the file
setwd("C://Dissertation R//R_Data")
ncname1 <- "posterior_climate_indices_MCruns_ensemble_full_LMRv2.1.nc"
dname1 <- "amo"  # name of the variable to extract
dname2 <- "sam"
# open the NetCDF file
library(ncdf4)
ncin1 <- nc_open(ncname1)
# obtain some basic information abou the file
print(ncin1)

time <- ncvar_get(ncin1, "time")/365 # time vector
sam <- ncvar_get(ncin1, dname2)
nc_close(ncin1)
# calculate ensemble mean
# You probably want to start with the "grand mean," which is an average over the MCrun dimension
sam <- apply(sam, c(1,3), mean)
sam_sd <- apply(sam, 2, sd)
# bandpass high-frequency variability
sam_bp <- bandpass(cbind(seq(1,2001),sam[1,]),
                   fhigh=0.04, genplot=FALSE, verbose=FALSE)
########
# Plot
par(mar=c(5.1, 4.1, 4.1, 5.1)) # change margin figure to accommodate second axis label
plot(NH_2k[,1] - lag, IHG_quant[,2], type='l', 
     xlab = "years AD", ylab = "SH minus NH (per mil)", 
     ylim=c(-12,6), xlim=c(2000,0))
# lines(NH_2k[,1], IHG_quant[,1], lty=3, lwd=0.5)
# lines(NH_2k[,1], IHG_quant[,3], lty=3, lwd=0.5)
abline(v = seq(2000, 0, -50), lty = 3, lwd= 0.5)
# abline(v = 1850, lty = 3, lwd= 2.5, col=2) # PI boundary
#### SAM assimilation ####
par(new=TRUE)
plot(sam[1,], type='l', 
     col="dodgerblue3", lwd=0.25, xlim=c(2000,0), ylim=c(4,-12), xlab="", ylab="", 
     xaxt='n', yaxt='n')
arrows(time, sam[1,]+2*sam_sd, time, sam[1,]-2*sam_sd, code=0, col="dodgerblue3", lwd=0.35)
lines(sam_bp, col="dodgerblue3", lwd=2)
axis(4, col="dodgerblue3")
mtext("SAM index", side=4, line=2)
legend("topright", c("IH 14C gradient", "SAM (data assimilation)"), col=c(1,"dodgerblue3"), 
       lty=c(1,1), cex=0.5)
par(mar=c(5.1, 4.1, 4.1, 2.1)) # reset margins to default values
#################################################
# save assimilated SAM data as table
write.table(cbind(time, sam[1,], sam_sd),
            file="~/Desktop/Ellie/SAM_DA.txt", row.names=FALSE,
            col.names=c("Years_AD", "median", "stdev"))
#################################################

# plot volcanic forcing
volc <- read.table(file="C://Dissertation R//R_Data//Sigl2015_GVF.txt", header=TRUE)
par(new=TRUE)
plot(volc, type='h', xlim=c(2000,0))
lines(smooth.spline(volc, spar=0.25), lwd=2, col=1)

#############
# remove Industrial Era (i.e. post 1850 AD) linear trend due to fossil fuel emissions in NH 
# using linear regression against CO2 data
install.packages('Hmisc')
require('Hmisc')
CO2 <- read.table("C://Dissertation R//R_Data//WDC_CO2_last1k.txt", header=T)
# interpolate CO2 data to res IHG over the historical period
C20_int <- as.data.frame(approxExtrap(CO2[,1], CO2[,2], NH_2k[1850:1950,1])) 
# linear model
mymodel <- lm(IHG_quant[1850:1950,2] ~ C20_int[,2])
summary(mymodel)
# 
model <- C20_int[,2]*mymodel$coefficients[2] + mymodel$coefficients[1]
remo <- IHG_quant[1850:1950,2] - model
offset <- remo[1] - IHG_quant[1850,2]
# calculate offset between raw and detrended data at 1850 AD
plot(NH_2k[1850:1950,1], remo - offset, type='l', col=2,
     ylim=c(-12,6), xlim=c(2000,0), xlab = "years AD", ylab = "SH minus NH (per mil)")
lines(NH_2k[,1], IHG_quant[,2])
# replace historical data with detrended values
IHG_corr <- c(IHG_quant[1:1849,2], remo - offset)

#################################################
# save IHG data as table
write.table(IHG_quant,
            file="C://Dissertation R//R_Data//myIHG_quant.txt", row.names=FALSE,
            col.names=c("5Cl", "median", "95Cl"))
write.table(cbind(NH_2k[,1], IHG_corr, IHG_sd),
            file="C://Dissertation R//R_Data//myIHG_corr.txt", row.names=FALSE,
            col.names=c("Years_AD", "median", "stdev"))
#################################################

SAM_low <- bandpass(cbind(SAM[,1], SAM[,2]), fhigh=0.04, genplot=FALSE, verbose=FALSE)
test <- approx(NH_2k[,1], IHG_corr, SAM[,1])

cor.test(test$y, SAM_low[,2], use = "complete.obs")
ccf(test$y, SAM_low[,2], na.action=na.omit, lag.max = 50)
