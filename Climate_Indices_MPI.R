require(ncdf4)
require(astrochron)
setwd("C:\\Dissertation R\\MPI-ESM1")

ncname1 <- "ts_Ayr_MPI-ESM1-2-LR_past2k_r1i1p1f1_gn"
ncname2 <- "psl_Ayr_MPI-ESM1-2-LR_past2k_r1i1p1f1_gn"
ncfname1 <- paste(ncname1, ".nc", sep = "")
ncfname2 <- paste(ncname2, ".nc", sep = "")
dname1 <- "ts"  # note: tmp means temperature (not temporary)
dname2 <- "psl"  

ncin1 <- nc_open(ncfname1)
ncin2 <- nc_open(ncfname2)
# # obtain some basic info
print(ncin1)
print(ncin2)

# get variables
lo <- ncvar_get(ncin1, "lon"); nlo <- dim(lo)
la <- ncvar_get(ncin1,"lat"); nla <- dim(la)
time <- ncvar_get(ncin1, "time")
ts_array <- ncvar_get(ncin1, dname1)
psl_array <- ncvar_get(ncin2, dname2)
dim(ts_array) # check dimension
dim(psl_array) # check dimension

# time vector
time <- seq(1:dim(ts_array)[3])

################
# calculate ENSO
# # your grid coordinates over which ENSO will be estimated
lon1 = -120 # West
lon2 = -170 # West
lat1 = -5
lat2 = 5

# get index of nearest lon/lat values
lon1a = which.min(abs(lon1 - lo))
lon2a = which.min(abs(lon2 - lo)) 
lat1a = which.min(abs(lat1 - la)) 
lat2a = which.min(abs(lat2 - la)) 

# sapply
ENSO <- sapply(1:length(time), function(i){ 
  mean(ts_array[lon1a:lon2a,lat1a:lat2a,i], na.rm=T)
}); ENSO <- cbind(time, ENSO)

# bandpass ENSO
ENSO_bp <- bandpass(ENSO, fhigh=0.04, verbose = F)
par(mfrow=c(1,1))

# save assimilated ENSO data as table
# to convert to Celsius in excel
write.table(cbind(time, ENSO[,2]),
            file="MPIESM1_ENSO.txt", row.names=FALSE,
            col.names=c("Years_AD", "ts"))
write.table(cbind(ENSO_bp),
            file="MPIESM1_ENSObp.txt", row.names=FALSE,
            col.names=c("Years_AD", "ts"))

# load Celsius ENSO files
# import dataset from base txt
ENSO_C <- MPIESM1_ENSO
ENSO_bp_C <- MPIESM1_ENSObp

# plot ENSO (degrees Celsius)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(ENSO_C, type='l', lwd=1,
     xlab = "Years AD", ylab = "ENSO (°C)", 
     xlim=c(0, 1850), ylim=c(28.6, 25.6), col='grey',
     main="MPI-ESM1 ENSO Record")
abline(h=27.5721, lty=3, lwd=1, col=1)
abline(h=28.0721, lty=3, lwd=1, col=1)
abline(h=27.0721, lty=3, lwd=1, col=1)
lines(ENSO_bp_C, lwd=2, col='royalblue3')

# add volcanic years
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}

# plot ENSO bandpass only
par(new=F)
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(ENSO_bp_C[,1], ENSO_bp_C[,2], lwd=2, type='l',
     xlab = "Years AD", ylab = "ENSO (°C)", 
     xlim=c(0, 1850), ylim=c(28,27), col='skyblue2', 
     main="MPI-ESM1 ENSO Record")
abline(h=27.607, lty=3, lwd=1, col=1)
abline(h=28.107, lty=3, lwd=1, col=1)
abline(h=27.107, lty=3, lwd=1, col=1)
abline(h=27.207, lty=3, lwd=1, col=1)
axis(side=1, at=c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000))
# add volcanic years
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}

################
## calculate SAM
# convert Pa to hPa
psl_corr_array <- psl_array/100
#
# Gong and Wang (1999) define the SAM index to be the difference in monthly 
# and zonal-mean sea-level pressure (psl) between 40Â°S and 65Â°S
mylat1 = -40 
mylat2 = -65 
myLat1 = which.min(abs(mylat1 - la)) 
myLat2 = which.min(abs(mylat2 - la)) 
#
psl_array_40S <- sapply(1:(dim(psl_corr_array)[3]), function(i){
    result <- mean(psl_corr_array[,myLat1,i], na.rm=TRUE)
  })
psl_array_65S <- sapply(1:(dim(psl_corr_array)[3]), function(i){
    result <- mean(psl_corr_array[,myLat2,i], na.rm=TRUE)
  })
SAM <- cbind(time, psl_array_40S - psl_array_65S)

# bandpass SAM
SAM_bp <- bandpass(SAM, fhigh=0.04, verbose = F)
par(mfrow=c(1,1))

# save assimilated ENSO data as table
write.table(cbind(time, ENSO[,2]),
            file="MPIESM1_SAM.txt", row.names=FALSE,
            col.names=c("Years_AD", "psl"))
write.table(cbind(time, SAM[,2]),
            file="MPIESM1_SAM_bp.txt", row.names=FALSE,
            col.names=c("Years_AD", "psl"))

# plot SAM and bandpass
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM, type='l', lwd=1,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(0, 1850), ylim=c(34,20), col='grey',
     main="MPI-ESM1 SAM Record")
abline(h=27.47224, lty=3, lwd=1, col=1)
abline(h=30.06394, lty=3, lwd=1, col=1)
abline(h=24.88053, lty=3, lwd=1, col=1)
lines(SAM_bp, lwd=2, col='orange')

# add volcanic years
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}

# plot SAM bandpass only
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM_bp, type='l', lwd=2,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(0, 1850), ylim=c(30,25), col='orange',
     main="MPI-ESM1 SAM Record")
abline(h=27.47224, lty=3, lwd=1, col=1)
axis(side=1, at=c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000))
# add volcanic years
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}


###################

# TOGETHER plot MPI ENSO and SAM with volcanism 
par(new=F)
par(mar=c(5.1, 4, 4, 5.1)) # change margin figure to accommodate second axis label
plot(ENSO_bp_C[,1], ENSO_bp_C[,2], type='l', lwd=2,
     xlab = "Years AD", ylab = "ENSO", 
     xlim=c(0, 1850), ylim=c(28,26.5), col='skyblue2',
     main="MPI-ESM1 ENSO vs SAM Record")
abline(h=27.607, lty=3, lwd=1, col='skyblue2')
par(new=TRUE)
plot(SAM_bp[,1],SAM_bp[,2], type='l', xlab="", ylab="",
     col='orange', lwd=2, xaxt='n', yaxt='n', ylim=c(32, 25))
abline(h=27.47224, lty=3, lwd=1, col='orange')
axis(4, col='orange')
mtext("SAM", side=4, line=3)
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}
legend("topright", c("ENSO", "SAM"), col=c('skyblue2','orange'), 
       lty=c(1,1), cex=0.7)

###################

# zoom into the SAM record in volcanic years to explore patterns
# 25 year average shows little change
# need to see if closer years to the eruption show more response

# 540 CE
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM, type='l', lwd=1,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(439,639), ylim=c(34,20), col='orange',
     main="MPI-ESM1 SAM Record 540 CE eruption")
abline(v=(540), lty=3, lwd=2)
abline(h=27.47224, lty=3, lwd=1, col=1)
abline(h=30.06394, lty=3, lwd=1, col=1)
abline(h=24.88053, lty=3, lwd=1, col=1)

# 682 CE
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM, type='l', lwd=1,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(581,781), ylim=c(34,20), col='orange',
     main="MPI-ESM1 SAM Record 682 CE eruption")
abline(v=(682), lty=3, lwd=2)
abline(h=27.47224, lty=3, lwd=1, col=1)
abline(h=30.06394, lty=3, lwd=1, col=1)
abline(h=24.88053, lty=3, lwd=1, col=1)

# 1258 CE
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM, type='l', lwd=1,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(1157,1357), ylim=c(34,20), col='orange',
     main="MPI-ESM1 SAM Record 1258 CE eruption")
abline(v=(1258), lty=3, lwd=2)
abline(h=27.47224, lty=3, lwd=1, col=1)
abline(h=30.06394, lty=3, lwd=1, col=1)
abline(h=24.88053, lty=3, lwd=1, col=1)

# 1458 CE
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM, type='l', lwd=1,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(1357,1557), ylim=c(34,20), col='orange',
     main="MPI-ESM1 SAM Record 1458 CE eruption")
abline(v=(1458), lty=3, lwd=2)
abline(h=27.47224, lty=3, lwd=1, col=1)
abline(h=30.06394, lty=3, lwd=1, col=1)
abline(h=24.88053, lty=3, lwd=1, col=1)

# 1815 CE
par(mar=c(5.1, 4.1, 4.1, 5.1))
plot(SAM, type='l', lwd=1,
     xlab = "Years AD", ylab = "SAM", 
     xlim=c(1714,1850), ylim=c(34,20), col='orange',
     main="MPI-ESM1 SAM Record 1815 CE eruption")
abline(v=(1815), lty=3, lwd=2)
abline(h=27.47224, lty=3, lwd=1, col=1)
abline(h=30.06394, lty=3, lwd=1, col=1)
abline(h=24.88053, lty=3, lwd=1, col=1)

#####################################################
# Compare to proxy response
# Plot MPI ENSO vs proxy ENSO with standard deviation 
# proxy ENSO has come from the LMR data extraction script
{
  par(mfrow=c(1,1))
  par(mar=c(5.1, 4.1, 4.1, 5.1)) # change margin figure to accommodate second axis label
  nino34_bp2000 <- nino34_bp[nino34_bp[,1] < 1851,] # shorten to finish at 1850 
  plot(nino34_bp2000[,1], nino34_bp2000[,2], type='n', col='grey40', xlim=c(0,1850), ylim=c(2,-1.5), lwd=0.5,
       xlab="Years AD", ylab="Nino 3.4",
       main="Comparing ENSO records over the past 2,000 years
     (with Common Era Volcanism)")
  par(new=TRUE)
  lines(nino34_bp2000[,1], nino34_bp2000[,2], lwd=2, type='l', col='navy')
  axis(side=1, at=c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000))
  }
# ENSO_bp_C ends at 1850 so needs shifting 
{
  par(new=TRUE)
  plot(ENSO_bp_C[,1], ENSO_bp_C[,2], type='l', lwd=2,
     xlab = "", ylab = "", xaxt='n', yaxt='n', 
     ylim=c(28,26.5), xlim=c(0,1850), col=2)
  lines(ENSO_bp_C, col=2, lwd=2)
# abline(h=27.607, lty=3, lwd=1, col=2)
#  abline(h=27.207, lty=3, lwd=1, col=2)
#  abline(h=28.107, lty=3, lwd=1, col=2)
  axis(4, col=2)
  mtext("MPI-ESM1 ENSO", side=4, line=3) 
  legend("topleft", c("Proxy ENSO (Nino 3.4)",
        "Simulated ENSO (MPI-ESM1)"),
        lty=c(1,1), col=c("navy", 2), cex=0.8)
}

# Add volcanic years
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}

########################
# Plot MPI SAM vs reconstructed SAM with standard deviation 
SAMcomp <- mySAM_calib[mySAM_calib[,1] < 1851,] # shorten to finish at 1850 
{
    plot(mySAMrec[,1], mySAMrec[,2], lwd=2, type='n', col='grey20', 
       xlab="Year AD", ylab="SAM reconstruction",  
       ylim=c(7, -11), xlim=c(0,1851), 
       main="Comparing SAM records over the past 2,000 years (with Common Era Volcanism)")
  lines(SAMcomp[,1], SAMcomp[,2], col='orange', lwd=2)
  lines(mySAMrec[,1], mySAMrec[,2], col="orange", lwd=2)
  axis(side=1, at=c(0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000))
}

# Add volcanic years
{
  abline(v=(540), lty=3, lwd=0.5)
  abline(v=(682), lty=3, lwd=0.5)
  abline(v=(1258), lty=3, lwd=0.5)
  abline(v=(1458), lty=3, lwd=0.5)
  abline(v=(1815), lty=3, lwd=0.5)
}


{
  par(new=TRUE)
  plot(SAM_bp[,1], SAM_bp[,2], type='l', lwd=1,
       xlab = "", ylab = "", xaxt='n', yaxt='n', 
       col='skyblue2', ylim=c(30,22))
  lines(SAM_bp[,1], SAM_bp[,2], col='mediumpurple3', lwd=2)
  axis(4, col='mediumpurple3')
  mtext("MPI-ESM1 SAM", side=4, line=3) 
  legend("topright", c("Proxy SAM (IHG reconstruction)",
                        "Simulated SAM (MPI-ESM1)"),
         lty=c(1,1), col=c("orange", 'mediumpurple3'), cex=0.8)
}

#############################################################################
############# Stacking Model Data of SAM response to volcanism ##############
#############################################################################

# create cluster subsets of SAM record from MPI-ESM1 - 100 years either side of eruption date
mytime <- seq(-100,100)
# 540 CE eruption - year zero: 539 
clusterA <- SAM_bp[SAM_bp[,1] > 438,]
cluster1 <- clusterA[clusterA[,1] < 640,]
cluster1[,2]=(cluster1[,2] - mean(cluster1[,2]))/sd(cluster1[,2])
cluster540 <- cbind(mytime, cluster1[,2])

# 682 CE unknown eruption - year zero: 681
clusterB <- SAM_bp[SAM_bp[,1] > 580,]
cluster2 <- clusterB[clusterB[,1] < 782,]
cluster2[,2]=(cluster2[,2] - mean(cluster2[,2]))/sd(cluster2[,2])
cluster682 <- cbind(mytime, cluster2[,2])

# 1258 CE Samalas - year zero: 1257
clusterC <- SAM_bp[SAM_bp[,1] > 1156,]
cluster3 <- clusterC[clusterC[,1] < 1358,]
cluster3[,2]=(cluster3[,2] - mean(cluster3[,2]))/sd(cluster3[,2])
cluster1258 <- cbind(mytime, cluster3[,2])

# 1458 CE unknown (Kuwae?) - year zero: 1457 
clusterD <- SAM_bp[SAM_bp[,1] > 1356,]
cluster4 <- clusterD[clusterD[,1] < 1558,]
cluster4[,2]=(cluster4[,2] - mean(cluster4[,2]))/sd(cluster4[,2])
cluster1458 <- cbind(mytime, cluster4[,2])

# 1815 CE Tambora - year zero: 1814 
clusterE <- SAM_bp[SAM_bp[,1] > 1713,]
cluster5 <- clusterE[clusterE[,1] < 1915,]
cluster5[,2]=(cluster5[,2] - mean(cluster5[,2]))/sd(cluster5[,2])
cluster1815 <- cbind(mytime, cluster5[,2])

# all eruptions in one cluster
all.clusters <- cbind(mytime, cluster1[,2], cluster2[,2], 
                      cluster3[,2], cluster4[,2], cluster5[,2])
volc.stack <- rowMeans(all.clusters[,2:6])
volc.sd <- apply(all.clusters[,2:6], 1, sd)
{
plot(mytime, volc.stack, type='l', ylim=c(-2,2.5), ylab = "SAM index", lwd=2,
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
#legend("bottomleft", c("Stacked Data", "540 CE", 
#"682 CE", "1258 CE", "1458CE", "1815 CE"),
#lty=c(1,1,1), col=c(1,2,3,4,5,6), cex=0.7)
}

#540 CE and 682 CE stack
clusters_540_682 <- cbind(cluster1[,2], cluster2[,2])
volc.stack <- rowMeans(clusters_540_682)
volc.sd <- apply(clusters_540_682, 1, sd)
{
  plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5), lwd=2,
     xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
     main="SAM response to volcanic forcing")
#lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
#lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
abline(v=0, lty=2, lwd=1.5)
abline(h=0, lty=2, lwd=1.5)
lines(mytime, cluster540[,2], col=2)
lines(mytime, cluster682[,2], col=3)
legend("topleft", c("Stacked Data", "540 CE", 
                       "682 CE"),
       lty=c(1,1,1), col=c(1, 2, 3), cex=0.7)
axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))
}


#540 CE, 682 CE and 1258 CE stack
clusters_540_682_1258 <- cbind(cluster1[,2], cluster2[,2], cluster3[,2])
volc.stack <- rowMeans(clusters_540_682_1258)
volc.sd <- apply(clusters_540_682_1258, 1, sd)
{
  plot(mytime, volc.stack, type='l', ylim=c(-2,2.5), lwd=2,
       xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
       main="SAM response to volcanic forcing")
  #lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
  #lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
  abline(v=0, lty=2, lwd=1.5)
  abline(h=0, lty=2, lwd=1.5)
  lines(mytime, cluster540[,2], col=2)
  lines(mytime, cluster682[,2], col=3)
  lines(mytime, cluster1258[,2], col=4)
  legend("topleft", c("Stacked Data", "540 CE", 
                      "682 CE", "1258 CE"),
         lty=c(1,1,1), col=c(1, 2, 3, 4), cex=0.7)
  axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))
}

#1258 CE, 1458 CE and 1815  CE stack
clusters_1258_1458_1815 <- cbind(cluster3[,2], cluster4[,2], cluster5[,2])
volc.stack <- rowMeans(clusters_1258_1458_1815)
volc.sd <- apply(clusters_1258_1458_1815, 1, sd)
{
  plot(mytime, volc.stack, type='l', ylim=c(-2.5,2.5), lwd=2,
       xlab="Year AD (relative to eruption year zero)", ylab="SAM index",
       main="SAM response to volcanic forcing")
  #lines(mytime, volc.stack+volc.sd, type='l', lty=3, col=2)
  #lines(mytime, volc.stack-volc.sd, type='l', lty=3, col=2)
  abline(v=0, lty=2, lwd=1.5)
  abline(h=0, lty=2, lwd=1.5)
  lines(mytime, cluster1258[,2], col=4)
  lines(mytime, cluster1458[,2], col=5)
  lines(mytime, cluster1815[,2], col=6)
    legend("topleft", c("Stacked Data", "1258 CE", 
                      "1458 CE", "1915 CE"),
         lty=c(1,1,1), col=c(1, 4, 5, 6), cex=0.7)
  axis(side=1, at=c(-100, -75, -50, -25, 0, 25, 50, 75, 100))
}
