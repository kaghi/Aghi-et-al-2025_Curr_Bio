#Required package
install.packages("R.matlab")
install.packages("spatstat")
install.packages("spatialEco") 
install.packages("plyr")
install.packages("MonteCarlo")
install.packages("ggplot2")
install.packages("WriteXLS")
install.packages("entropy")
install.packages("LaplacesDemon")
install.packages("CDFt")
install.packages("lemon")
install.packages("imager")

library(R.matlab)
library(spatstat)
library(spatialEco)
library(plyr)
library(MonteCarlo)
library(ggplot2)
library(WriteXLS)
library(entropy)
library(LaplacesDemon)
library(CDFt)
library(lemon)
library(imager)
#Set working directory. Make sure you use "/" and not "\" for the address.
setwd("C:/Users/Lucifer/Documents/Zach QuaSOR data/All Coords for Krisha June 2020/All Coords for Krisha/All WT BrpCac LF_Mini_HF QuaSORvsSTORM QuaSOR Coords ONLY")

#Load in file
# read in our data
quasor1b<- readMat("170810 P3LA3M4 WT_QuaSOR_Coords_ONLY.mat_R_ready_stability_comparison.mat")

# check out data structure
str(quasor1b)

#The data should be part of a particular structure. Extract that data. #

IbdataEvok<- quasor1b$QuasorIbDataEvokedIb
IbdataSpont<- quasor1b$QuasorIbDataSpontIb

IsdataEvok <- quasor1b$QuasorIbDataEvokedIs
IsdataSpont <- quasor1b$QuasorIbDataSpontIs

borderq <- quasor1b$mask
borderqline <- which(borderq==1, arr.ind=TRUE)

borderqline_x <- c(borderqline[1:nrow(borderqline),1])
borderqline_y <- c(borderqline[1:nrow(borderqline),2])

lines(borderqline)

TrialNo<- quasor1b$trial
#This section then lets us use Ripley's L function on our data. 


######################################################################
#Compare Evoked vs. Spont data for both Ib and Is

xvech_Ib_evok <- IbdataEvok[1:nrow(IbdataEvok),1]
yvech_Ib_evok <- IbdataEvok[1:nrow(IbdataEvok),2]
time_points_evok <- IbdataEvok[1:nrow(IbdataEvok),3]

plotdata_Ib_evok <-data.frame(xvech_Ib_evok,yvech_Ib_evok)

xvech_Ib_spont <- IbdataSpont[1:nrow(IbdataSpont),1]
yvech_Ib_spont <- IbdataSpont[1:nrow(IbdataSpont),2]
plotdata_Ib_spont <-data.frame(xvech_Ib_spont,yvech_Ib_spont)

#Plot the coordinates of Ib evoked vs spontaneous

lines(quasor1b$maskpoly[])
plot(xvech_Ib_evok, yvech_Ib_evok, col = "Red")

par(pty='s')
plot(xvech_Ib_spont, yvech_Ib_spont, col = "Purple", pch = 1, cex = 0.5, asp = 1)
par(new = TRUE)
plot(xvech_Ib_evok, yvech_Ib_evok, col = "Green", pch = 1, cex = 0.5, asp = 1)
par(new = TRUE)
plot(borderqline, pch = 20, cex = 0.2, asp = 1)

dev.off()
#MPP Data
w_e <- owin(c(0,max(IbdataEvok)), c(0,max(IbdataEvok)))

w_s <- owin(c(0,max(IbdataSpont)), c(0,max(IbdataSpont)))

ppdata_Ib_evoked <- ppp(xvech_Ib_evok,yvech_Ib_evok,w_e)
ppdatalf_Ib_spont <- ppp(xvech_Ib_spont,yvech_Ib_spont,w_s)

#PCf computation 
g12evok <- pcf(ppdata_Ib_evoked, method="d", spar=0.7)
g12spont <- pcf(ppdatalf_Ib_spont, method="d", spar=0.7)

intpoint_evoked <- nndensity(ppdata_Ib_evoked)
lambda_evoked <- mean(intpoint_evoked$v)
Oringevoked <- eval.fv(lambda_evoked * g12evok)
plot(Oringevoked, ylab = "O-ring statistic O(r) ")


intpoint_spont <- nndensity(ppdatalf_Ib_spont)
lambda_spont <- mean(intpoint_spont$v)
Oringspont <- eval.fv(lambda_spont * g12spont)
plot(Oringspont, ylab = "O-ring statistic O(r) ")

#We do this because the first point is read as infinity so we remove it. 
revoked <- Oringevoked$iso[-1]
rspont <- Oringspont$iso[-1]
rvals <- Oringevoked$r[-1]

#Ks stat
evokedcdf <- ecdf(revoked)
plot(evokedcdf, col ="red", main = "Empirical CDFs of 0.2 Hz Stimulation vs. Spontaneous Coordinates ", ylab = "Cumulative probability", xlim=c(0,1) )
spontcdf <- ecdf(rspont)
lines(spontcdf, col = "blue")

#Add legend
legend("right", legend = c("0.2 Hz Stimulation", "Spontaneous"), col = c("red", "blue"), pch = c(16))

#Set your z to ensure that you are encapsulating the entire distribution!
z = seq(0, 1, by=0.01)
evokedvals <- evokedcdf(z)
spontvals <-spontcdf(z)

res <- CramerVonMisesTwoSamples(evokedvals,spontvals)
pvalue = 1/6*exp(-res)

file_name <- paste("WT_GC6_25C_MiniEvoked8_7_ME2_QuaSOR_Coords_ONLY.mat_R_ready_with_time_points_", "ProcessedStats.mat")
writeMat(con = "WT_GC6_25C_MiniEvoked8_7_ME2_QuaSOR_Coords_ONLY.mat_R_ready_with_time_points.mat", filename=file_name, Ib_evok_coord = plotdata_Ib_evok, Ib_spont_coord = plotdata_Ib_spont, Oringevok = Oringevoked, Oringspont = Oringspont, g_evok = g12evok, g_spont = g12spont, lambda_evok = lambda_evoked, lambda_spont = lambda_spont, intpoint_evok = intpoint_evoked, intpoint_spont = intpoint_spont, evokedvals = evokedvals, spontvals = spontvals, res = res, x_axis = z, p_value = pvalue, file = file_name)

##########################################################################
#########################################################################
#Compare two halves of data set to HF data

halfindlf<- nrow(IbdataEvok)/2
xvecfhalf <- IbdataEvok[1:halfindlf,1]
yvecfhalf <- IbdataEvok[1:halfindlf,2]
w <- owin(c(0,max(IbdataEvok)), c(0,max(IbdataEvok)))

plotdatafhalf <-data.frame(xvecfhalf,yvecfhalf)


xveclhalf <- IbdataEvok[halfindlf:nrow(IbdataEvok),1]
yveclhalf <- IbdataEvok[halfindlf:nrow(IbdataEvok),2]
plotdatalhalf <-data.frame(xveclhalf,yveclhalf)

plot(xvecfhalf,yvecfhalf, xlab = "X-coordinates", ylab = "Y-coordinates", col="Blue")
plot(xveclhalf, yveclhalf, xlab = "X-coordinates", ylab = "Y-coordinates", col="Red")
points(xveclhalf,yveclhalf)


ppdatafhalf <- ppp(xvecfhalf,yvecfhalf,w)
ppdatalhalf <- ppp(xveclhalf,yveclhalf,w)

g12fhalf <- pcf(ppdatafhalf, method="d", spar=0.7)
g12lhalf <- pcf(ppdatalhalf, method="d", spar=0.7)

intpointfhalf <- nndensity(ppdatafhalf)
lambdafhalf <- mean(intpointfhalf$v)
Oringfhalf <- eval.fv(lambdafhalf * g12fhalf)
plot(Oringfhalf, ylab = "O-ring statistic O(r) ")
points(Oringl50)

intpointlhalf <- nndensity(ppdatalhalf)
lambdalhalf <- mean(intpointlhalf$v)
Oringlhalf <- eval.fv(lambdalhalf * g12lhalf)
plot(Oringlhalf, ylab = "O-ring statistic")

#K-l divergence

#We do this because the first point is read as infinity so we remove it. 
fhalf <- Oringfhalf$iso[-1]
lhalf <- Oringlhalf$iso[-1]
rvalshalf <- Oringfhalf$r[-1]

# kldwithintrial <- KLD(fhalf,lhalf)
# plot(rvalshalf, kldwithintrial$mean.KLD)
# kldvaltrialval<- kldwithintrial$mean.sum.KLD
# 
# kldwithintrialdf <- as.data.frame(kldwithintrial)


#Ks stat
fhalfcdf <- ecdf(fhalf)
plot(fhalfcdf, col ="red", main = "Empirical CDFs of first half and last half of LF data against HF data ", ylab = "Cumulative probability", xlim=c(0,0.4) )
lhalfcdf <- ecdf(lhalf)
lines(lhalfcdf, col = "blue")
lines(hfcdf)

#Set your z to ensure that you are encapsulating the entire distribution!
z = seq(0, 1, by=0.01)
firsthalfvals <- fhalfcdf(z)
lasthalfvals <-lhalfcdf(z)

res <- CramerVonMisesTwoSamples(firsthalfvals,lasthalfvals)
pvalue = 1/6*exp(-res)

writeMat(con = "170810 P3LA3M4 WT_QuaSOR_Coords_ONLY.mat_R_ready_stability_comparison.mat", firsthalfvals = firsthalfvals, lasthalfvals = lasthalfvals, res = res, x_axis = z, p_value = pvalue)


#Add legend
legend("right", legend = c("First Half of LF Values", "Second Half of LF Values", "HF Values"), col = c("red", "blue", "black"), pch = c(16))
##########################################################################
#Examine trial to trial variation

xvecf50 <- IbdataLF[1:50,1]
yvecf50 <- IbdataLF[1:50,2]
w <- owin(c(0,max(IbdataLF)), c(0,max(IbdataLF)))

plotdataf50 <-data.frame(xvecf50,yvecf50)

lindst <- nrow(IbdataLF)-50
linden <- nrow(IbdataLF)
xvecl50 <- IbdataLF[lindst:linden,1]
yvecl50 <- IbdataLF[lindst:linden,2]
plotdatal50 <-data.frame(xvecl50,yvecl50)

plot(xvecf50,yvecf50, xlab = "X-coordinates", ylab = "Y-coordinates", col="Red")
plot(xvecl50, yvecl50, xlab = "X-coordinates", ylab = "Y-coordinates", col="Red")
points(xvecl50,yvecl50)


ppdataf50 <- ppp(xvecf50,yvecf50,w)
ppdatal50 <- ppp(xvecl50,yvecl50,w)

g12f50 <- pcf(ppdataf50, method="d", spar=0.7)
g12l50 <- pcf(ppdatal50, method="d", spar=0.7)

intpointf50 <- nndensity(ppdataf50)
lambdaf50 <- mean(intpointf50$v)
Oringf50 <- eval.fv(lambdaf50 * g12f50)
plot(Oringf50, ylab = "O-ring statistic O(r) ")
points(Oringl50)

intpointl50 <- nndensity(ppdatal50)
lambdal50 <- mean(intpointl50$v)
Oringl50 <- eval.fv(lambdal50 * g12l50)
plot(Oringl50, ylab = "O-ring statistic")

#K-l divergence

f50 <- Oringf50$iso[-1]
l50 <- Oringl50$iso[-1]
rvals50 <- Oringf50$r[-1]

kldwithintrial <- KLD(f50,l50)
plot(rvals50, kldwithintrial$mean.KLD)
kldvaltrialval<- kldwithintrial$mean.sum.KLD

kldwithintrialdf <- as.data.frame(kldwithintrial)


#Ks stat
f50cdf <- ecdf(f50)
plot(f50cdf, col ="red", main = "Empirical CDFs of first and last 100 values", ylab = "Cumulative probability")
l50cdf <- ecdf(l50)
lines(l50cdf)

#Set your z to ensure that you are encapsulating the entire distribution!
z = seq(0, 0.15, by=0.01)
f50vals <- f50cdf(z)
l50vals <-l50cdf(z)

res <- CramerVonMisesTwoSamples(f50vals,l50vals)
pvalue = 1/6*exp(-res)



###########Compare low vs. high frequency data
IbdataHF<- quasor1b$QuasorIbDataHighF
IbdataLF<- quasor1b$QuasorIbDataLowF

xvechf <- IbdataHF[1:nrow(IbdataHF),1]
yvechf <- IbdataHF[1:nrow(IbdataHF),2]
plotdatahf <-data.frame(xvechf,yvechf)

xveclf <- IbdataLF[1:nrow(IbdataLF),1]
yveclf <- IbdataLF[1:nrow(IbdataLF),2]
plotdatalf <-data.frame(xveclf,yveclf)

plot(xvechf, yvechf, col = "Red")
plot(xveclf, yveclf)
points(xvechf, yvechf)
#MPP Data
ppdatahf <- ppp(xvechf,yvechf,w)
ppdatalf <- ppp(xveclf,yveclf,w)

#PCf computation 
g12hf <- pcf(ppdatahf, method="d", spar=0.7)
g12lf <- pcf(ppdatalf, method="d", spar=0.7)

intpointhf <- nndensity(ppdatahf)
lambdahf <- mean(intpointhf$v)
Oringhf <- eval.fv(lambdahf * g12hf)
plot(Oringhf, ylab = "O-ring statistic O(r) ")

intpointlf <- nndensity(ppdatalf)
lambdalf <- mean(intpointlf$v)
Oringlf <- eval.fv(lambdalf * g12lf)
plot(Oringlf, ylab = "O-ring statistic O(r)")


#K-L divergence

hfdata <- Oringhf$iso[-1]
lfdata <- Oringlf$iso[-1]
rvalshflf <- Oringhf$r[-1]

plot(rvalshflf, hfdata, col = "red", main= "PDF of coordinate data using O-ring statistic", ylab = "O-ring stat (O(r))")
points(lfdata)

legend(95, 95, legend = c("High-Freq", "Low-Freq"), col = c("red", "black"))

kldhflf <- KLD(hfdata,lfdata)
plot(rvalshflf, kldhflf$mean.KLD, main = "Mean-sum K-L divergence between Low and High Frequency PDFs")

#Ks stat
hfcdf <- ecdf(hfdata)
plot(hfcdf, col ="red", main = "Empirical CDFs of High-frequency vs. Low-frequency data", ylab = "Cumulative probability")
lfcdf <- ecdf(lfdata)
lines(lfcdf)
lines(l50cdf)
lines()

z = seq(0, 1.5, by=0.01)
lfcdfvals <- lfcdf(z)
hfcdfvals <-hfcdf(z)

plot(lfcdfvals, col="red")                  
points(hfcdfvals)  

res <- CramerVonMisesTwoSamples(hfcdfvals, lfcdfvals)
pvalue = 1/6*exp(-res)



#####################


data <- data.frame("Cramer-von Mises 2-sample p-values"= c(3.14E-11, 6.62E-06,1.43E-06, 1.89E-06, 1.97E-10, 1.53E-06, 1.38E-08, 1.48E-04, 5.76E-07, 1.14E-08, 1.42E-04, 3.31E-08, 0.01, 0.06, 0.17, 0.1415, 0.118, 0.143, 0.114, 0.124, 0.106, 0.132, 0.16, 0.0983), "Comparison" = factor(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 2,2,2,2,2,2,2,2,2,2,2,2)), labels = "HF-LF, LF-LF")
hflfpvals <- c(3.14E-11, 6.62E-06,1.43E-06, 1.89E-06, 1.97E-10, 1.53E-06, 1.38E-08, 1.48E-04, 5.76E-07, 1.14E-08, 1.42E-04, 3.31E-08
)


lflfpvals <- c(0.01, 0.06, 0.17, 0.1415, 0.118, 0.143, 0.114, 0.124, 0.106, 0.132, 0.16, 0.0983
)

x = "Hf-Lf pvals"
qplot(hflfpvals, lflfpvals, geom = "boxplot")



hflfkldiv <- c(0.0602, 0.03508456, 0.066, 0.014, 0.019, 0.018, 0.04, 0.03, 0.08, 0.017, 0.0277, 0.015
)


lflfkldiv <- c(0.003416652, 0.05871603, 0.0085, 0.0097, 0.0389, 0.0359, 0.0122, 0.0258, 0.02, 0.053, 0.026
)

boxplot(hflfkldiv, lflfkldiv, names = c("HF-HF (n=12)", "HF-LF (n=12)"), ylab = "Mean-Sum K-L Divergence", main = "K-L Divergence Values")
#####################


xvec <- Ibdata[1:nrow(Ibdata),1]
yvec <- Ibdata[1:nrow(Ibdata),2]
zvecind <- Ibdata[1:2223,3]
w <- owin(c(0,300), c(0,300))

#make dataframe
plotdataf <-data.frame(xvec,yvec)
qplot(xvec,yvec, xlab = "X-coordinates", ylab = "Y-coordinates" )

plotdata <-matrix(c(xvec,yvec), nrow = nrow(Ibdata), ncol =2 )
qplot(mpg, data=plotdataf, geom="density", alpha=I(.5), 
      main="Distribution of Gas Milage", xlab="Miles Per Gallon", 
      ylab="Density")




ppdata <-ppp(xvec,yvec,w)
qplot(mpg, data=plotdata, geom="density", fill=gear, alpha=I(.5))
      
L <- Lest(ppdata)

LEnv <- plot(envelope(ppdata, Kest),main="Monte Carlo estimate of L-function") 


MC_result<-MonteCarlo(func=Lfunc, nrep=1000, )

Ltheo <- (L$theo)^2
Kcross <- (Ltheo)*pi
g12 <- pcf(ppdata, method="d", spar=0.7)


#theo isw the L(r) = r aka an estimate of a stationary Poisson process.
plot(g12)
O<- o.ring(ppdata)
intpoint <- nndensity(ppdata)
lambda2 <- mean(intpoint$v)
Oring <- eval.fv(lambda2 * g12)
plot(Oring, ylab = "O-ring statistic")

file_name <- paste("170804 P1LA4M4 WT SynapGC6f 0.2Hz WT_GC6f_02Hz_Mini_HF_13 All QuaSOR Coordinate Data Exportconverted.mat", "ProcessedStats.mat")
save(plotdataf, LEnv, Oring, g12, O, lambda2, lambda, intpoint, file = file_name)


#For any number of data points, we then determine if the statistic is reliable across multiple sets of data. 
my_vectorx <- vector("numeric")
my_vectory <- vector("numeric")

CleanTrialNo <- compact(TrialNo)

for (i in 1:200)
{ 
  mini[i] = laply((CleanTrialNo[i]),min)
  maxi[i] = laply((CleanTrialNo[i]),max) 
  my_vectorx[i]<- list(Ibdata[mini[i]:maxi[i],1])
  my_vectory[i]<- list(Ibdata[mini[i]:maxi[i],2])
   
}

plot(totaldata)

totaldata <- c(my_vectorx,my_vectory)

marks(ppdata$x)<- factor("xcoord")
K12 <- Kcross(ppdata, "x", "y")
g12 <- (pcf(ppdata, method="d", spar=0.7))

Oring <-eval.fv(lambda*g12)

plot(Oring)
Oring <- eval.fv(lambda2 * g12)
plot(Oring)
#Then calculated the K-S statistic
OringHF <- Oring
Oring

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}
cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
  # Copy a data.frame to clipboard
  write.table(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
}

write.excel(OringHFdf)
OringHFdf<- as.data.frame(OringHF)
copyf <- write.xls(OringHFdf)





