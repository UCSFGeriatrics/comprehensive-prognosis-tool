
library(pec)
#library(prodlim)
# library(lava)
library(riskRegression)
library(survival)
library(coxphw)
library(cmprsk)
library(timereg)
library(haven)
library(survey)
library(dplyr)

packageDescription("pec", fields = "Version")
# [1] "2020.11.17"

packageDescription("riskRegression", fields = "Version")
# [1] "2020.12.08"

# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)


load("clinicalhrsdeident_gdr_20210624.Rdata")

#Check type of variables
str(mydata)

# Check columns classes
sapply(mydata, class)

#Need to make heart3g, smoking3g, and stroke3g factor variables
mydata$heart3g <- factor(mydata$heart3g)
mydata$smoking3g <- factor(mydata$smoking3g)
mydata$stroke3g <- factor(mydata$stroke3g)
sapply(mydata, class)


###############################################################################################################################################
######### DEATH ##########

#############################################
# Fit model without weights

table(mydata$death)
summary(subset(mydata, death == 0)$time2death)
summary(subset(mydata, death == 1)$time2death)

fitdth <- coxph(Surv(time2death,death)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
                                       R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
                                       bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,data=mydata,x=TRUE,y=TRUE)
fitdth #Checked the estimates and se, and they are similar to the ones obtained in SAS PHREG

####### Calibration plot
# calfitdth10<-calPlot(list("Unweighted Cox regression Death until time=10"=fitdth),time=10, type="risk", data=mydata)
# calfitdth10

# Email from Dr. Thomas Gerds (author of package pec):
# pec::calPlot is outdated. please use
# x=riskRegression::Score(...,plots="cal") and then
# plotCalibration(x,cens.method="local"...) the handling of censored data in particular has improved for the calibration plots
# (not using pseudo values anymore)

xs=Score(list("Unweighted Cox Death t=10"=fitdth),
         formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")

# conf.int:	
# Either logical or a numeric value between 0 and 1. In right censored data,
# confidence intervals are based on Blanche et al (see references). 
# Setting FALSE prevents the computation confidence intervals.
# TRUE means compute 95 percent confidence intervals and corresponding p-values for AUC and Brier score.
# If set to 0.87, the level of significance is 13 percent. So, do not set it to 0.87

xs
calfitdth10<-plotCalibration(xs,cens.method="local") 

# cens.method	
# For right censored data only. How observed proportions are calculated. Either "jackknife" or "local":
# "jackknife" Compute a running mean of the jackknife pseudovalues across neighborhoods/groups of the predicted risks.
#             Here we rely on the assumption that censoring is independent of the event time and the covariates, see References.
# "local"     Compute the Kaplan-Meier estimator in absence of competing risks and the Aalen-Johansen estimator in
#             presence of competing risks locally like a running mean in neighborhoods of the predicted risks.
#             The widths of the neighborhoods are defined according to method.
names(calfitdth10)
str(calfitdth10)
calfitdth10$bandwidth # 0.04876844

# Extract data for plot
calfitdth10$plotFrames
plotfr<-as.data.frame(calfitdth10$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.10)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('plotfr', 'loessMod', 'smoothed'))


####### Calibration bar plot
# calfitdth10bars<-calPlot(list("Unweighted Cox regression Death until time=10"=fitdth),time=10, type="risk", bars=TRUE,hanging=FALSE, showFrequencies=TRUE, data=mydata)
# calfitdth10bars

calfitdth10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitdth10bars)
calfitdth10bars$plotFrames
str(calfitdth10bars$plotFrames)
plotfr<-as.data.frame(calfitdth10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100
plotfr

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.40)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')
plotfr

rm(list = c('xs','plotfr', 'loessMod', 'smoothed'))


#############################################
#Fit model with weights

#First attempt
R5WTRESP<-mydata$R5WTRESP
fitdthw <- coxph(Surv(time2death,death)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
                   R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
                   bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,data=mydata,
                 weight=R5WTRESP, x=TRUE,y=TRUE)
fitdthw

xs=Score(list("Weighted Cox Death t=10"=fitdthw),
         formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")
# Error in predictCox(object = object, newdata = newdata, times = times,  : 
#                       predictCox does not know how to handle Cox models fitted with weights

#Import data with SurvProb for death at 10 yrs
surv10_death <- read.csv('SurvPredDth10.csv')
cum10_death<-surv10_death
cum10_death$cif_death<-1-cum10_death$s_death #compute predicted event probabilities
str(cum10_death)
cum10_death<-cum10_death[,3]

####### Calibration plot
# calfitdthw10<-calPlot(list("Weighted Cox regression Death until time=10"=surv10_death),time=10, type="risk",
#                       formula=Surv(time2death,death)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
#                         R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
#                         bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,
#                       data=mydata)
# calfitdthw10

xs=Score(list("Weighted Cox Death t=10"=cum10_death),
                       formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")
xs
calfitdthw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
# calfitdthw10bars<-calPlot(list("Weighted Cox regression Death until time=10"=surv10_death),time=10, type="risk", bars=TRUE,hanging=FALSE,showFrequencies=TRUE,
#                      formula=Surv(time2death,death)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
#                        R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
#                        bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,
#                      data=mydata)
# calfitdthw10bars

calfitdthw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitdthw10bars)
calfitdthw10bars$plotFrames
str(calfitdthw10bars$plotFrames)
plotfr<-as.data.frame(calfitdthw10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100
plotfr

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')
plotfr

rm(list = c('xs','plotfr', 'loessMod', 'smoothed'))

###############################################################################################################################################
######### ADL with core and exit interviews ##########

cov <- subset(mydata, subgroup_adldep==1, select=c(age3_lin, age3_sp1, age3_sp2, gender, R5EATA, R5CANCRE, R5DIABE,
                                                   heart3g, R5HIBPE, R5MEALSA, R5MONEYA, h5lvalone, R5LUNGE, R5PUSHA,
                                                   R5WALKSA, bmi3_lin, bmi3_sp1, bmi3_sp2, smoking3g, stroke3g))
cov<-data.matrix(cov)
time <- mydata[which(mydata$subgroup_adldep==1), "time_adldepdth2"]
time<-time$time_adldepdth2

status <- subset(mydata, subgroup_adldep==1, select=status_adldepdth2)
status<-status$status_adldepdth2

#############################################
# Fit model without weights

# fitadl <- crr(time,status,cov, failcode=1, cencode=0)
# fitadl
# calfitADL10=calPlot(fitadl,time=10,data=cov)
# #Error in calPlot(fitadl, time = 10, data = cov) : Argument formula is missing and first model has no usable formula: 

#Another way using Formula wrapper for crr from cmprsk (FGR). Formula interface for Fine-Gray regression competing risk models
mydata2<- subset(mydata, subgroup_adldep==1)

fgradl <- FGR(Hist(time_adldepdth2,status_adldepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
                R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
                bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,data=mydata2,cause=1)
fgradl #Checked the estimates and se, and they are similar to the ones obtained in SAS PHREG
#predict.crr <- cmprsk:::predict.crr

####### Calibration plot
# calfitadl10<-calPlot(list("Unweighted Competing-risk ADL (core/exit) until time=10"=fgradl),time=10,data=mydata2)
# calfitadl10

xs=Score(list("Unweighted FG ADL t=10"=fgradl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
# se.fit: Logical or 0 or 1. If FALSE or 0 do not calculate standard errors.
xs
calfitadl10<-plotCalibration(xs,cens.method="local") 

# Extract data for plot
names(calfitadl10)
calfitadl10$plotFrames
plotfr<-as.data.frame(calfitadl10$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.10)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs','plotfr', 'loessMod', 'smoothed'))

####### Calibration bar plot
# calfitadl10bars<-calPlot(list("Unweighted Competing-risk ADL (core/exit) until time=10"=fgradl),time=10, bars=TRUE,hanging=FALSE, data=mydata2)
# calfitadl10bars

calfitadl10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitadl10bars)
calfitadl10bars$plotFrames
plotfr<-as.data.frame(calfitadl10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs','plotfr', 'loessMod', 'smoothed'))


#############################################
#Fit model with weights
R5WTRESP<-mydata2$R5WTRESP
fgradlw <- FGR(Hist(time_adldepdth2,status_adldepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
                 R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
                 bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,
               weight=R5WTRESP, data=mydata2,cause=1)
#It doesnt allow weights

#Import data with CumInc for ADL at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,2]
class(cif10_adl)

####### Calibration plot
# Note: Hist provides the suitable extensions for dealing with right censored and interval censored data from competing risks
# and other multi state models

# calfitadlw10<-calPlot(list("Weighted Competing-risk ADL (core/exit) until time=10"=cif10_adl),time=10, data=mydata2,cause=1,
#                      formula=Hist(time_adldepdth2,status_adldepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
#                        R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
#                        bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g)
# calfitadlw10

xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
# calfitadlw10bars<-calPlot(list("Weighted Competing-risk ADL (core/exit) until time=10"=cif10_adl),time=10, data=mydata2,cause=1, bars=TRUE,hanging=FALSE,showFrequencies=TRUE,
#                       formula=Hist(time_adldepdth2,status_adldepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
#                         R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
#                         bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g)
# calfitadlw10bars

calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 


# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.40)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs', 'mydata2', 'plotfr', 'loessMod', 'smoothed'))

###############################################################################################################################################
######### Walk with core and exit interviews ##########

cov <- subset(mydata, subgroup_walkdep==1, select=c(age3_lin, age3_sp1, age3_sp2, gender, R5EATA, R5CANCRE, R5DIABE,
                                                    heart3g, R5HIBPE, R5MEALSA, R5MONEYA, h5lvalone, R5LUNGE, R5PUSHA,
                                                    R5WALKSA, bmi3_lin, bmi3_sp1, bmi3_sp2, smoking3g, stroke3g))
cov<-data.matrix(cov)
time <- mydata[which(mydata$subgroup_walkdep==1), "time_walkdepdth2"]
time<-time$time_walkdepdth2

status <- subset(mydata, subgroup_walkdep==1, select=status_walkdepdth2)
status<-status$status_walkdepdth2

#Another way using Formula wrapper for crr from cmprsk (FGR). Formula interface for Fine-Gray regression competing risk models
mydata2<- subset(mydata, subgroup_walkdep==1)

#############################################
# Fit model without weights
fgrwalk <- FGR(Hist(time_walkdepdth2,status_walkdepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
                 R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
                 bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,data=mydata2,cause=1)
fgrwalk #Checked the estimates and se, and they are similar to the ones obtained in SAS PHREG
#predict.crr <- cmprsk:::predict.crr

####### Calibration plot
# calfitwalk10<-calPlot(list("Weighted FG Walk t=10"=fgrwalk),time=10,data=mydata2)
# calfitwalk10

xs=Score(list("Unweighted FG Walk t=10"=fgrwalk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
# se.fit: Logical or 0 or 1. If FALSE or 0 do not calculate standard errors.
xs
calfitwalk10<-plotCalibration(xs,cens.method="local") 

# Extract data for plot
names(calfitwalk10)
calfitwalk10$plotFrames
plotfr<-as.data.frame(calfitwalk10$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.10)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('plotfr', 'loessMod', 'smoothed'))

####### Calibration bar plot
# calfitwalk10bars<-calPlot(list("Weighted FG Walk t=10"=fgrwalk),time=10, bars=TRUE,hanging=FALSE, data=mydata2)
# calfitwalk10bars

calfitwalk10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalk10bars)
calfitwalk10bars$plotFrames
plotfr<-as.data.frame(calfitwalk10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs','plotfr', 'loessMod', 'smoothed'))


#############################################
#Fit model with weights
R5WTRESP<-mydata2$R5WTRESP
fgrwalkw <- FGR(Hist(time_walkdepdth2,status_walkdepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
                  R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
                  bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g,
                weight=R5WTRESP, data=mydata2,cause=1)
#It doesnt allow weights

#Import data with CumInc for Walk at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,2]
class(cif10_walk)

####### Calibration plot
# Note: Hist provides the suitable extensions for dealing with right censored and interval censored data from competing risks
# and other multi state models

# calfitwalkw10<-calPlot(list("Weighted Competing-risk Walk (core/exit) until time=10"=cif10_walk),time=10, data=mydata2,cause=1,
#                      formula=Hist(time_walkdepdth2,status_walkdepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
#                        R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
#                        bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g)
# calfitwalkw10

xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
# calfitwalkw10bars<-calPlot(list("Weighted Competing-risk Walk (core/exit) until time=10"=cif10_walk),time=10, data=mydata2,cause=1, bars=TRUE,hanging=FALSE,showFrequencies=TRUE,
#                       formula=Hist(time_walkdepdth2,status_walkdepdth2)~age3_lin+age3_sp1+age3_sp2+gender+R5EATA+R5CANCRE+R5DIABE+heart3g+
#                         R5HIBPE+R5MEALSA+R5MONEYA+h5lvalone+R5LUNGE+R5PUSHA+R5WALKSA+
#                         bmi3_lin+bmi3_sp1+bmi3_sp2+smoking3g+stroke3g)
# calfitwalkw10bars

calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100

# Make plot from extracted data
par(las=1)
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed <- predict(loessMod)
lines(smoothed, x=plotfr$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs', 'mydata2', 'plotfr', 'loessMod', 'smoothed'))

