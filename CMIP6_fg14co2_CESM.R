library(ncdf4) # load the package (make sure the package is installed first)
library(fields)

##################################################
# Open and extract a variable from a netcdf file
##################################################
setwd("C:\\Dissertation R\\CESM2\\CESM2\\piControl\\fg14co2abio")
# name of the file
ncname1 <- "fg14co2abio_Oyr_CESM2_piControl_r1i1p1f1_gr_0101-1200.nc"
dname1 <- "fg14co2abio"  # name of the variable to extract
# open the NetCDF file
ncin1 <- nc_open(ncname1)
# obtain some basic information abou the file
print(ncin1)

setwd("C:\\Dissertation R\\CESM2\\CESM2\\piControl\\fgco2abio")
# name of the file
ncname2 <- "fgco2abio_Oyr_CESM2_piControl_r1i1p1f1_gr_0101-1200.nc"
dname2 <- "fgco2abio"  # name of the variable to extract
# open the NetCDF file
ncin2 <- nc_open(ncname2)
# obtain some basic information abou the file
print(ncin2)

#########################################################################
# now let's extract some variables (e.g. longitude, latitude, time, ect.)
lon <- ncvar_get(ncin1, "lon") # longitude
nlon <- dim(lon) # number of lon grid points 
lat <- ncvar_get(ncin1,"lat") # latitude
nlat <- dim(lat) # number of lat grid points 
time <- ncvar_get(ncin1, "time") # time vector
dim(time) # length time (remember that this is the number of months!)
# now extract the flux data 
f14co2 <- ncvar_get(ncin1, dname1)
fco2abio <- ncvar_get(ncin2, dname2)
# create a time vector in years (see file info 'print(ncin)')
year <- seq(101,1200) # this is the time span of the simulation
#
# if we are done extracting variables from the netCDF file we should 
# terminates access to the netCDF
nc_close(ncin1) 
nc_close(ncin2) 
#########################################################################
rs14 = 0.85*(10^12)  # normalizing the flux to the pre-anthropogenic late 19th century atmospheric value of C12/C14 (based on Orr et al. 2017)
# ^ same as start of results to Rodgers et al (2011) - creating the scaling factor (rs)
# NB f14co2 variable is Surface Downward Mass Flux of Carbon-14 as Abiotic 14CO2 [kgC m-2 s-1] 
# convert downward fluxes to upward fluxes (i.e. sea-air flux)
f14co2_array <- (-1*f14co2*83.259093974539) # convert Kg C to moles C m-2 s-1 and multiple by C12/C14 pre-Suess ratio
f14co2_array_sb <- f14co2_array[,,1:100]
f14co2_array_avg <- apply(f14co2_array_sb, c(1,2), mean, na.rm = TRUE) #avg for panoply

fco2abio_array <- (-1*fco2abio*83.259093974539)
fco2abio_array_sb <- fco2abio_array[,,1:100]
fco2abio_array_avg <- apply(fco2abio_array_sb, c(1,2), mean, na.rm = TRUE) #avg for panoply

# calculate disequlibrium flux
dis = (f14co2_array - fco2abio_array)*rs14 # sea-air flux of 14CO2 - sea-air flux of CO2
dis_sb <- dis[,,1:100]
dis_sb_avg <- apply(dis_sb, c(1,2), mean, na.rm = TRUE) #avg for panoply

# plot - to check first timestep
image.plot(z=f14co2_array[,,1], y=lat, x=lon,
           xlab="Lon", ylab="Lat", main="14CO2 flux (mol C m-2 s-1)")
image.plot(z=fco2abio_array[,,1], y=lat, x=lon, 
           xlab="Lon", ylab="Lat", main="CO2 flux (mol C m-2 s-1)")
image.plot(z=dis[,,1], y=lat, x=lon, 
           xlab="Lon", ylab="Lat", main="Disequilibirium Flux (mol C m-2 s-1)")

# calculate disequilibrium flux
# ∆14C = 0 ‰ (as per Rodgers et al., 2011)
# first calculate mean state (1100 years)
DQ_array <- apply(dis[,,1:1100], 1:2, mean, na.rm=TRUE)
# creates a 2D array with the same spatial coordinates as the original area but with time removed (as it is averaged over all timesteps)
# Zonal  integral  of  the  disequilibrium  flux for the last 1000 model yrs - summing all longitudinal points
DQ_zonal <- sapply(1:length(lat), function(i){
  DQ <- sum(DQ_array[,i], na.rm=TRUE)
})
plot(lat, DQ_zonal, type='l', ylab = "mol C m-1 s-1", main="Zonal Integral DQ (as per Rodgers fig 2d)")
# negative flux = net input of co2 into the ocean all the time

# total disequlibrium into the ocean - calculating mass balance
DG_array_100 <- dis[,,1:100]
DG_array_glob <- apply(DG_array_100, 1:2, mean, na.rm=TRUE) # mean over tiem - 2D object
test <- sum(DG_array_glob, na.rm = TRUE) *(1/83.259093974539) *3.154e+7 # Kg C m-2 yr-1 (atomic mass of CO2 is 44)
# 1 kg of CO2 is equivalent to 0.27 kg of C
test2 <- test*1e-12 *1e+6 # in Gt C km-2 yr-1
test3 <- test2/(3.61e+08 - 0.5415e+8) # Gt C yr-1 (Oceans area minus sea-ice covered regions (ca. 15%))

# calculate percentage of disequilibrium flux between the Northern and Southern Hemispheres
sum(DQ_zonal[1:90], na.rm=TRUE)/sum(DQ_zonal, na.rm=TRUE) # SH
sum(DQ_zonal[91:180], na.rm=TRUE)/sum(DQ_zonal, na.rm=TRUE) # NH
sh <- (test3/100)*68
nh <- (test3/100)*31 #GtMCE
# interhemispheric flux 
abs(sh-nh)
# interhemispheric exchange time of order 1.3 yr (Geller et al., 1997)
# interhemispheric difference in the atmospheric radiocarbon inventory
igh <- abs(sh-nh)*1.3 # in GtMCE
# preindustrial tropospheric carbon inventory of 480 GtC –or 240 Gt C/hemisphere
(igh/(240)) * 1000 # interhemispheric gradient in per mil
# interhemispheric difference due to air sea fluxes contributes to 2.0 per mill

# NB: The mean IHG for the last 2kyr based on data from IntCal20/SHCal20 is -4.5 per mil +/-2.2 (1 sigma)
(4.5 - (igh/(240) * 1000))/4.5 # fractional distance from mean PI IHG
# how far in fractional % from the mean value calculated from tree rings

# disequlibrium flux SH (<25S) and NH (>25N)
# separate plots because overlap doesn't work
mylat1 = 40 
mylat2 = 85 
myLat1 = which.min(abs(mylat1 - lat)) 
myLat2 = which.min(abs(mylat2 - lat)) 
# disequilibrium for NH
DQ_NH <- sapply(1:length(year), function(i){
  DQ <- mean(dis[,myLat1:myLat2, i], na.rm=TRUE)
})

par(mar=c(11.1, 8.1, 11.1, 8.1)) 
plot(year, DQ_NH, type='l', lwd=0.5, xlim=c(100,1200), col='skyblue2', 
      xlab='Year', ylab='NH DQ flux',
      ylim=c(-16000,-13000))
lines(smooth.spline(cbind(year, DQ_NH), spar=0.4), lwd=2, col='navy')
#plot(smooth.spline(cbind(year, DQ_NH), spar=0.4), lwd=2, col='skyblue2',
#     xlim=c(100,1200), main="Disequilibirum Flux (CESM2)", 
#     xlab='Year', ylab='NH DQ flux', ylim=c(-15200,-14200), type='l')
abline(h=-14800.17, col="navy", lty=3, lwd=2)
# par(new=TRUE)

mylat1 = -85 # can change these boundaries but do not swap values
mylat2 = -40 
myLat1 = which.min(abs(mylat1 - lat)) 
myLat2 = which.min(abs(mylat2 - lat)) 
# disequilibrium for SH
DQ_SH <- sapply(1:length(year), function(i){
  DQ <- mean(dis[,myLat1:myLat2, i], na.rm=TRUE)
})
par(new=TRUE)
par(mar=c(11.1, 8.1, 11.1, 8.1))
plot(year, DQ_SH, col='pink', type='l', xlab="Year", ylab="SH DQ flux", 
     lwd=0.5, xlim=c(100,1200), ylim=c(-43000,-33000))
lines(smooth.spline(cbind(year, DQ_SH), spar=0.4), lwd=2, col='darkred')
# plot(smooth.spline(cbind(year, DQ_SH)), type='l', xlab="", ylab="", 
#      spar=0.4, col='pink', lwd=2, xlim=c(100,1200), 
#      xaxt='n', yaxt='n', ylim=c(-38000,-34000))
#axis(4, col="pink")
#mtext("SH DQ flux", side=4, line=3)
abline(h=-35414.23, col="darkred", lty=3, lwd=2)
#legend("topleft", c("SQ NH (>40 N)","DQ SH (>40 S)"), col=c('skyblue2', 'pink'), 
#       lty=c(1,1), cex=0.7)

par(new=TRUE)
plot(year, DQ_SH, type='l', lwd=0.5, xlim=c(100,1200), col='pink', 
     xlab="", ylab="", xaxt='n', yaxt='n',
     ylim=c(-43000,-33000))
lines(smooth.spline(cbind(year, DQ_SH), spar=0.4), lwd=2, col='darkred')
axis(4, col="darkred")
mtext("SH DQ flux", side=4, line=3)

legend("topleft", c("DQ NH", "DQ SH"), col=c('navy','darkred'), 
       lty=c(1,1), cex=0.7)

# ∆DQ
dDQ <- DQ_SH - DQ_NH
par(mar=c(11.1, 8.1, 11.1, 8.1))
plot(year, dDQ, col='plum', type='l', xlab="Year", ylab="DQ flux (NH-SH)", 
     lwd=0.5, xlim=c(100,1200))
lines(smooth.spline(cbind(year, dDQ), spar=0.4), lwd=2, col='mediumpurple4')
abline(h=-20604.76, col="mediumpurple4", lty=3, lwd=2)


# read previously calculated AMOC at 26°N
# AMOC_CESM2 <- read.table(file="C:\\Dissertation R\\CESM2\\CESM2\\piControl\\AMOC_Oyr_CESM2_esm-piControl_r1i1p1f1_gn_0001-1200.txt", header=TRUE)
# AMO_CESM2 <- read.table(file="C:\\Dissertation R\\CESM2\\CESM2\\piControl\\AMO_Oyr_CESM2_esm-piControl_r1i1p1f1_gn_0001-1200.txt", header=TRUE)
SAM_CESM2 <- read.table(file="C:\\Dissertation R\\CESM2\\CESM2\\piControl\\mySAM.txt", header=TRUE)
 
# CESM SAM 
par(mar=c(11.1, 8.1, 11.1, 8.1))
par(new=TRUE)
plot(SAM_CESM2, type='l', col='wheat1', xlab="", ylab="", 
     xlim=c(0,1200), ylim=c(40,32), lwd=0.5, xaxt='n', yaxt='n')
lines(smooth.spline(SAM_CESM2, spar=0.4), lwd=2, col='orange')
abline(h=36.91004, col="orange", lty=3, lwd=2)
axis(4, col="orange")
mtext("SAM", side=4, line=3)

legend("topleft", c("DQ SH", "SAM"), col=c('darkred','orange'), 
       lty=c(1,1), cex=0.7)

# Plot SAM vs DQ SH smoothed
SAM_smooth <- (smooth.spline(SAM_CESM2, spar=0.4))
DQ_SH_smooth <- (smooth.spline(cbind(year, DQ_SH), spar=0.4))
par(mar=c(11.1, 8.1, 11.1, 8.1))
plot(DQ_SH_smooth, col='darkred', type='l', xlab="Year", ylab="SH DQ flux", 
     xlim=c(100,1200), lwd=2, ylim=c(-37750, -34500))
abline(h=-35414.23, col="darkred", lty=3, lwd=2)
par(new=TRUE)
plot(SAM_smooth, type='l', col='orange',xlab="", ylab="", 
     xlim=c(0,1200), ylim=c(38,35), lwd=2, xaxt='n', yaxt='n')
abline(h=36.91004, col="orange", lty=3, lwd=2)
axis(4, col="orange")
mtext("SAM", side=4, line=3)
par(new=TRUE)
legend("topleft", c("DQ SH", "SAM"), col=c(1,'orange'), 
       lty=c(1,1), cex=0.7)


plot(SAM_CESM2, type='l', col='wheat1', xlab="Year", ylab="SAM", 
     xlim=c(0,1200), ylim=c(40,32), lwd=0.5)
lines(smooth.spline(SAM_CESM2, spar=0.4), lwd=2, col='orange')
abline(h=36.91004, col="orange", lty=3, lwd=2)

# need to adjust presentation 
par(mar=c(5.1, 4.1, 4.1, 5.1)) #change margin to fit axis

plot(year, DQ_SH, col='grey90', type='l', lwd=0.5, xlim=c(0,1200), ylim=c(-39000,-33000),
     xlab="Year", ylab="SH DQ flux", main = "SH DQ vs SAM (CESM2)")
lines(smooth.spline(cbind(year, DQ_SH), spar=0.4), lwd=2, col=1)
par(new=TRUE)
plot(SAM_CESM2, type='l', col='lavender', xlab="", ylab="", xlim=c(0,1200), ylim=c(42,32), lwd=0.5, xaxt='n', yaxt='n')
lines(smooth.spline(SAM_CESM2, spar=0.4), lwd=2, col='slateblue3')
axis(4, col=1)
mtext("SH DQ flux", side=4, line=2)
par(new=TRUE)
legend("topright", c("DQ SH", "SAM"), col=c(1,'slateblue3'), 
       lty=c(1,1), cex=0.7)

# spatial correlation between SAM and DQ - regress AMO on dDQ
mycor <- apply(dis, 1:2, 
                    function(x) cor(x, SAM_CESM2[101:1200,2], use="pairwise.complete.obs"))

image.plot(z=mycor, x=lon, y=lat, zlim=c(-0.75,0.75),
           xlab="Lon", ylab="Lat", 
           main="Correlation of DQ flux against SAM index - CESM2") 


# export for plotting in Panoply - SAM and DQ correlation
myArray <- mycor # your 2D correlation/regression array
x <- lon # lon coordinates
y <- lat # lat coordinates
# Now let's define and compile variable dimensions (lon/lat) and variable (i.e. corr coefficients)
dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('r', 'corr coefficient', # this is the temperature variable
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'Correlation of SAM index and fg14co2 field')
# Create an empty netcdf file...
outnc <- nc_create('myOutput.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, myArray)
nc_close(outnc) # terminate access to the file

#########################################

# f14co2_array export for plotting in Panoply
setwd("C:\\Dissertation R\\CESM2\\CESM2\\piControl\\fg14co2abio")
f14co2_array_CESM_pan <- f14co2_array_avg
x <- lon
y <- lat

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('f14co2', 'fg14co2 flux',
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'fg14co2 flux')

# Create an empty netcdf file...
outnc <- nc_create('f14co2_flux_cesm.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, f14co2_array_CESM_pan)
nc_close(outnc) # terminate access to the file

##########################################

# fco2abio_array export for plotting in Panoply
setwd("C:\\Dissertation R\\CESM2\\CESM2\\piControl\\fgco2abio")
fco2abio_array_CESM_pan <- fco2abio_array_avg
x <- lon
y <- lat

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('fco2', 'fco2abio flux',
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'fco2abio flux')

# Create an empty netcdf file...
outnc <- nc_create('fco2abio_flux_cesm.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, fco2abio_array_CESM_pan)
nc_close(outnc) # terminate access to the file

############################################

# dis export for plotting in Panoply
setwd("C:\\Dissertation R\\CESM2\\CESM2\\piControl\\dis")
dis_CESM_pan <- dis_sb_avg
x <- lon
y <- lat

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('dis', 'disequilibrium flux',
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'disequilibrium flux CESM')

# Create an empty netcdf file...
outnc <- nc_create('dis_flux_cesm.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, dis_CESM_pan)
nc_close(outnc) # terminate access to the file

