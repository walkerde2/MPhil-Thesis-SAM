library(ncdf4) # load the package (make sure the package is installed first)
library(fields)

##################################################
# Open and extract a variable from a netcdf file
##################################################
setwd("C:\\Dissertation R\\MRI-ESM\\MRI-ESM\\piControl\\fg14co2abio")
# name of the file
ncname1 <- "fg14co2abio_Oyr_MRI-ESM2-0_piControl_r1i2p1f1_gr_1850-2100.nc"
dname1 <- "fg14co2abio"  # name of the variable to extract
# open the NetCDF file
ncin1 <- nc_open(ncname1)
# obtain some basic information abou the file
print(ncin1)
#
setwd("C:\\Dissertation R\\MRI-ESM\\MRI-ESM\\piControl\\fgco2abio")
# name of the file
ncname2 <- "fgco2abio_Oyr_MRI-ESM2-0_piControl_r1i2p1f1_gr_1850-2100.nc"
dname2 <- "fgco2abio"  # name of the variable to extract
# open the NetCDF file
ncin2 <- nc_open(ncname2)
# obtain some basic information about the file
print(ncin2)


#########################################################################
# now let's extract some variables (e.g. longitude, latitude, time, ect.)
lon <- ncvar_get(ncin1, "lon") # longitude
nlon <- dim(lon) # number of lon grid points 
lat <- ncvar_get(ncin1,"lat") # latitude
nlat <- dim(lat) # number of lat grid points 
time <- ncvar_get(ncin1, "time") # time vector
dim(time) # length time (remember that this is the number of months!)
tunits <- ncatt_get(ncin1,"time","units") # info on time units
# now extract the slp data 
f14co2 <- ncvar_get(ncin1, dname1)
fco2abio <- ncvar_get(ncin2, dname2)
# create a time vector in years (see file info 'print(ncin)')
year <- seq(1850,2100) # this is the time span of the simulation

# if we are done extracting variables from the netCDF file we should 
# terminates access to the netCDF
nc_close(ncin1) 
nc_close(ncin2) 
#########################################################################
rs14 = 0.85*(10^12)  # normalizing the flux to the pre-anthropogenic late 19th century atmospheric value of C12/C14 (based on Orr et al. 2017)
# it equals 1/(1.170*(10^-12))

# NB f14co2 variable is Surface Downward Mass Flux of Carbon-14 as Abiotic 14CO2 [kgC m-2 s-1] 
# convert downward fluxes to upward fluxes (i.e. sea-air flux) and
# convert KgC m-2 s-1 to molesC m-2 s-1 using the multiplicative inverse of the Molar mass of C (12.0107 g/mol)
f14co2_array <- (-1*f14co2*83.259093974539) 
f14co2_array_sb <- f14co2_array[,,1:100]
f14co2_array_avg <- apply(f14co2_array_sb, c(1,2), mean, na.rm = TRUE) #avg for panoply

fco2abio_array <- (-1*fco2abio*83.259093974539)
fco2abio_array_sb <- fco2abio_array[,,1:100]
fco2abio_array_avg <- apply(fco2abio_array_sb, c(1,2), mean, na.rm = TRUE) #avg for panoply

dis_sb <- dis[,,1:100]
dis_sb_avg <- apply(dis_sb, c(1,2), mean, na.rm = TRUE) #avg for panoply

# disequlibrium
dis = (f14co2_array - fco2abio_array)*rs14 # sea-air flux of 14CO2 - sea-air flux of CO2
# plot
image.plot(z=f14co2_array[,,1], y=lat, x=lon, 
           xlab="Lon", ylab="Lat", main="14CO2 flux (mol C m-2 s-1) MRI-ESM")
image.plot(z=fco2abio_array[,,1], y=lat, x=lon, 
           xlab="Lon", ylab="Lat", main="CO2 flux (mol C m-2 s-1) MRI-ESM")
image.plot(z=dis[,,1], y=lat, x=lon, 
           xlab="Lon", ylab="Lat", main="Disequilibrium Flux (mol C m-2 s-1) MRI-ESM")

# take 50 timestep average of each array

#####################################
#want to avg these ^ for panoply use#
#####################################

# calculate disequilibrium flux
# ∆14C = 0 ‰ (as per Rodgers et al., 2011)
# first calculate mean state (1000 years)
DQ_array <- apply(dis[,,1:251], 1:2, mean, na.rm=TRUE)

# Zonal  integral  of  the  disequilibrium  flux for the last 100 model yrs
DQ_zonal <- sapply(1:length(lat), function(i){
  DQ <- sum(DQ_array[,i], na.rm=TRUE)
})
plot(lat, DQ_zonal, type='l', ylab = "mol C m-1 s-1", main="Zonal Integral of DQ MRI-ESM")

# total disequlibrium into the ocean
DG_array_100 <- dis[,,1:251]
DG_array_glob <- apply(DG_array_100, 1:2, mean, na.rm=TRUE)
test <- sum(DG_array_glob, na.rm = TRUE) * (1/83.259093974539) *3.154e+7 # Kg C m-2 yr-1 (atomic mass of CO2 is 44)
# 1 kg of CO2 is equivalent to 0.27 kg of C
test2 <- test*1e-12 *1e+6 # in Gt C km-2 yr-1
test3 <- test2/(3.61e+08 - 0.5415e+8) # Gt C yr-1 (Oceans area minus sea-ice covered regions (ca. 15%))

# calculate percentage of disequilibrium flux between the Northern and Southern Hemispheres
sum(DQ_zonal[1:90], na.rm=TRUE)/sum(DQ_zonal, na.rm=TRUE) # SH
sum(DQ_zonal[91:180], na.rm=TRUE)/sum(DQ_zonal, na.rm=TRUE) # NH
sh <- (test3/100)*77
nh <- (test3/100)*22
# interhemispheric flux 
abs(sh-nh)
# interhemispheric exchange time of order 1.3 yr (Geller et al., 1997)
# interhemispheric difference in the atmospheric radiocarbon inventory
igh <- abs(sh-nh)*1.3 # in GtMCE
# preindustrial tropospheric carbon inventory of 480 GtC –or 240 Gt C/hemisphere
(igh/(240)) * 1000 # interhemispheric gradient in per mil

# NB: The mean IHG for the last 2k based on data from IntCal20 is -4.5 per mil +/-2.2 (1 sigma)
(4.5 - (igh/(240) * 1000))/4.5 # fractional distance from mean PI IHG
# disequlibrium flux SH (<25S) and NH (>25N)
# separate plots - overlap doesn't work
par(mar=c(5.1, 4.1, 4.1, 5.1)) #change margin to fit axis
mylat1 = 40 
mylat2 = 85 
myLat1 = which.min(abs(mylat1 - lat)) 
myLat2 = which.min(abs(mylat2 - lat)) 

DQ_NH <- sapply(1:length(year), function(i){
  DQ <- mean(dis[,myLat1:myLat2, i], na.rm=TRUE)
})
plot(year, DQ_NH, type='l', lwd=0.5, xlim=c(1850,2100), col='skyblue2', main="NH (>25N) Disequilibirum (MRI-ESM)", xlab='Year', ylab='NH DQ flux')
lines(smooth.spline(cbind(year, DQ_NH), spar=0.4), lwd=2, col='navy')

#
mylat1 = -85 
mylat2 = -40 
myLat1 = which.min(abs(mylat1 - lat)) 
myLat2 = which.min(abs(mylat2 - lat)) 

DQ_SH <- sapply(1:length(year), function(i){
  DQ <- mean(dis[,myLat1:myLat2, i], na.rm=TRUE)
})
# par(new=TRUE)
plot(year, DQ_SH, col='pink', xlab="Year", ylab="SH DQ flux", type='l', lwd=0.5,
     xlim=c(1850,2100), ylim=c(-75000,-60000), main="SH (>25S) Disequilibrium (MRI-ESM)")
lines(smooth.spline(cbind(year, DQ_SH), spar=0.4), lwd=2, col='darkred')
# axis(4, col=2)
# mtext("SH_DQ", side=4, line=2)

# ∆DQ
dDQ <- DQ_SH - DQ_NH

SAM_MRI <- read.table(file="C:\\Dissertation R\\MRI-ESM\\MRI-ESM\\piControl\\mySAM.txt", header=TRUE)
SAM_MRI[,1] <- year

par(mar=c(5.1, 4.1, 4.1, 5.1)) #change margin to fit axis
plot(year, DQ_SH, col='grey85', type='l', lwd=0.5, xlim=c(1850,2100), ylim=c(-75000,-60000), main="SH DQ vs SAM (MRI-ESM)",
     xlab="Year", ylab="SH DQ flux")
lines(smooth.spline(cbind(year, DQ_SH), spar=0.4), lwd=2, col=1)
par(new=TRUE)
plot(SAM_MRI, type='l', col='lavender', xlab="", ylab="", xlim=c(1850,2100), ylim=c(34,22), lwd=0.15, xaxt='n', yaxt='n')
lines(smooth.spline(SAM_MRI, spar=0.4), col='slateblue3', lwd=2)
axis(4, col='slateblue3')
mtext("SAM", side=4, line=2)
legend("topright", c("DQ SH", "SAM"), col=c(1,'slateblue3'), 
       lty=c(1,1), cex=0.7)


# regress SAM on dDQ
mycor <- apply(dis, 1:2, 
               function(x) cor(x, SAM_MRI[,2], use="pairwise.complete.obs"))

image.plot(z=mycor, x=lon, y=lat, zlim=c(-0.75,0.75),
           xlab="Lon", ylab="Lat", 
           main="Correlation of DQ flux against SAM index - MRI-ESM") 

# export for plotting in Panoply
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
setwd("C:\\Dissertation R\\MRI-ESM\\MRI-ESM\\piControl\\fg14co2abio")
f14co2_array_MRI_pan <- f14co2_array_avg
x <- lon
y <- lat

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('f14co2', 'fg14co2 flux',
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'fg14co2 flux')

# Create an empty netcdf file...
outnc <- nc_create('f14co2_flux_mri.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, f14co2_array_MRI_pan)
nc_close(outnc) # terminate access to the file

##########################################

# fco2abio_array export for plotting in Panoply
setwd("C:\\Dissertation R\\MRI-ESM\\MRI-ESM\\piControl\\fgco2abio")
fco2abio_array_MRI_pan <- fco2abio_array_avg
x <- lon
y <- lat

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('fco2', 'fco2abio flux',
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'fco2abio flux')


# Create an empty netcdf file...
outnc <- nc_create('fco2abio_flux_mri.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, fco2abio_array_MRI_pan)
nc_close(outnc) # terminate access to the file

############################################

# dis export for plotting in Panoply
setwd("C:\\Dissertation R\\MRI-ESM\\MRI-ESM\\piControl\\dis")
dis_MRI_pan <- dis_sb_avg
x <- lon
y <- lat

dim1 <- ncdim_def('lon', 'degrees_E', as.double(x)) # dimension 1 (longitude)
dim2 <- ncdim_def('lat', 'degrees_N', as.double(y)) # dimension 2 (latitude)
varz <- ncvar_def('dis', 'disequilibrium flux',
                  list(dim1,dim2), -9.9999e+36,
                  longname = 'disequilibrium flux MRI')

# Create an empty netcdf file...
outnc <- nc_create('dis_flux_mri.nc',
                   varz, force_v4 = FALSE)
# ...and finally add the pre-compiled data array with relative information
ncvar_put(outnc, varz, dis_MRI_pan)
nc_close(outnc) # terminate access to the file

