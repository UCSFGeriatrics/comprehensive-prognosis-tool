R.version$version.string
# "R version 4.1.1 (2021-08-10)"

# install.packages("pec", lib="C:/R4.1.1.library")
# install.packages("survival", lib="C:/R4.1.1.library")
# install.packages("coxphw", lib="C:/R4.1.1.library")
# install.packages("haven", lib="C:/R4.1.1.library")
# install.packages("survey", lib="C:/R4.1.1.library")
# install.packages("dplyr", lib="C:/R4.1.1.library")
# install.packages("plotrix", lib="C:/R4.1.1.library")

library(prodlim, lib.loc="C:/R4.1.1.library")
library(pec, lib.loc="C:/R4.1.1.library")
library(riskRegression, lib.loc="C:/R4.1.1.library")
library(survival, lib.loc="C:/R4.1.1.library")
library(coxphw, lib.loc="C:/R4.1.1.library")
library(cmprsk, lib.loc="C:/R4.1.1.library")
library(timereg, lib.loc="C:/R4.1.1.library")
library(haven, lib.loc="C:/R4.1.1.library")
library(survey, lib.loc="C:/R4.1.1.library")
library(dplyr, lib.loc="C:/R4.1.1.library")
library(plotrix, lib.loc="C:/R4.1.1.library")


packageDescription("pec", fields = "Version")
# [1] "2021.10.11"

packageDescription("riskRegression", fields = "Version")
# [1] "2021.10.10"


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
################################################ DEATH 1-unavailable predictors ###############################################################

####################################
#Fit 16-predictor model with weights

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

rm(list = c('cum10_death','surv10_death', 'calfitdthw10', 'calfitdthw10bars', 'xs'))


####################################
#15-variable best model missing 1 variable (R5EATA): lowest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cum10_death <- read.csv('SurvPredDth10_1UnavailBest.csv')
str(cum10_death)
cum10_death<-cum10_death[,4]

xs=Score(list("Weighted Cox Death t=10"=cum10_death),
         formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")
xs
calfitdthw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitdthw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitdthw10bars)
calfitdthw10bars$plotFrames
str(calfitdthw10bars$plotFrames)
plotfr_1UnavailBest<-as.data.frame(calfitdthw10bars$plotFrames)
colnames(plotfr_1UnavailBest) <- c('Predicted', 'Observed')
plotfr_1UnavailBest$Predicted<-plotfr_1UnavailBest$Predicted*100
plotfr_1UnavailBest$Observed<-plotfr_1UnavailBest$Observed*100
plotfr_1UnavailBest

rm(list = c('cum10_death', 'calfitdthw10', 'calfitdthw10bars', 'xs'))


####################################
#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cum10_death <- read.csv('SurvPredDth10_1UnavailWorst.csv')
str(cum10_death)
cum10_death<-cum10_death[,4]

xs=Score(list("Weighted Cox Death t=10"=cum10_death),
         formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")
xs
calfitdthw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitdthw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitdthw10bars)
calfitdthw10bars$plotFrames
str(calfitdthw10bars$plotFrames)
plotfr_1UnavailWorst<-as.data.frame(calfitdthw10bars$plotFrames)
colnames(plotfr_1UnavailWorst) <- c('Predicted', 'Observed')
plotfr_1UnavailWorst$Predicted<-plotfr_1UnavailWorst$Predicted*100
plotfr_1UnavailWorst$Observed<-plotfr_1UnavailWorst$Observed*100
plotfr_1UnavailWorst

rm(list = c('cum10_death', 'calfitdthw10', 'calfitdthw10bars', 'xs'))


####################################
# Make plot from extracted data
par(las=1)

#Plot in png format
png('DeathBaseBestWorst_1unavail_R.png', width=4.5, height=3.25, units="in", res=800, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black",lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_1UnavailBest$Predicted, plotfr_1UnavailBest$Observed, col="blue", pch=20)
loessMod_1UnavailBest <- loess(Observed ~ Predicted, data=plotfr_1UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailBest <- predict(loessMod_1UnavailBest)
lines(smoothed_1UnavailBest, x=plotfr_1UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_1UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_1UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailWorst <- predict(loessMod_1UnavailWorst)
lines(smoothed_1UnavailWorst, x=plotfr_1UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 1 predictor unavailable", "Worst 1 predictor unavailable"), bty='n')
dev.off()


#Plot in pdf format
pdf('DeathBaseBestWorst_1unavail_R.pdf', width=4.88, height=3.25, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black", lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_1UnavailBest$Predicted, plotfr_1UnavailBest$Observed, col="blue", pch=20)
loessMod_1UnavailBest <- loess(Observed ~ Predicted, data=plotfr_1UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailBest <- predict(loessMod_1UnavailBest)
lines(smoothed_1UnavailBest, x=plotfr_1UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_1UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_1UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailWorst <- predict(loessMod_1UnavailWorst)
lines(smoothed_1UnavailWorst, x=plotfr_1UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 1 predictor unavailable", "Worst 1 predictor unavailable"), bty='n')
dev.off()

#Base
plotfr
#1UnavailBest
plotfr_1UnavailBest
#1UnavailWorst
plotfr_1UnavailWorst


rm(list = ls())


###############################################################################################################################################
################################################ DEATH 2-unavailable predictors ###############################################################

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

####################################
#Fit 16-predictor model with weights

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

rm(list = c('cum10_death','surv10_death', 'calfitdthw10', 'calfitdthw10bars', 'xs'))


####################################
#14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cum10_death <- read.csv('SurvPredDth10_2UnavailBest.csv')
str(cum10_death)
cum10_death<-cum10_death[,4]

xs=Score(list("Weighted Cox Death t=10"=cum10_death),
         formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")
xs
calfitdthw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitdthw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitdthw10bars)
calfitdthw10bars$plotFrames
str(calfitdthw10bars$plotFrames)
plotfr_2UnavailBest<-as.data.frame(calfitdthw10bars$plotFrames)
colnames(plotfr_2UnavailBest) <- c('Predicted', 'Observed')
plotfr_2UnavailBest$Predicted<-plotfr_2UnavailBest$Predicted*100
plotfr_2UnavailBest$Observed<-plotfr_2UnavailBest$Observed*100
plotfr_2UnavailBest

rm(list = c('cum10_death', 'calfitdthw10', 'calfitdthw10bars', 'xs'))


####################################
#14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cum10_death <- read.csv('SurvPredDth10_2UnavailWorst.csv')
str(cum10_death)
cum10_death<-cum10_death[,4]

xs=Score(list("Weighted Cox Death t=10"=cum10_death),
         formula=Surv(time2death,death)~1,data=mydata,conf.int=TRUE,times=10, plots="cal")
xs
calfitdthw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitdthw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitdthw10bars)
calfitdthw10bars$plotFrames
str(calfitdthw10bars$plotFrames)
plotfr_2UnavailWorst<-as.data.frame(calfitdthw10bars$plotFrames)
colnames(plotfr_2UnavailWorst) <- c('Predicted', 'Observed')
plotfr_2UnavailWorst$Predicted<-plotfr_2UnavailWorst$Predicted*100
plotfr_2UnavailWorst$Observed<-plotfr_2UnavailWorst$Observed*100
plotfr_2UnavailWorst

rm(list = c('cum10_death', 'calfitdthw10', 'calfitdthw10bars', 'xs'))


####################################
# Make plot from extracted data
par(las=1)

#Plot in png format
png('DeathBaseBestWorst_2Unavail_R.png', width=4.5, height=3.25, units="in", res=800, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black",lwd=.5)

#14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC
#points(plotfr_2UnavailBest$Predicted, plotfr_2UnavailBest$Observed, col="blue", pch=20)
loessMod_2UnavailBest <- loess(Observed ~ Predicted, data=plotfr_2UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailBest <- predict(loessMod_2UnavailBest)
lines(smoothed_2UnavailBest, x=plotfr_2UnavailBest$Predicted, col="blue", lty=3)

#14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC
loessMod_2UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_2UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailWorst <- predict(loessMod_2UnavailWorst)
lines(smoothed_2UnavailWorst, x=plotfr_2UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 2 predictors unavailable", "Worst 2 predictors unavailable"), bty='n')
dev.off()


#Plot in pdf format
pdf('DeathBaseBestWorst_2Unavail_R.pdf', width=4.88, height=3.25, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black", lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_2UnavailBest$Predicted, plotfr_2UnavailBest$Observed, col="blue", pch=20)
loessMod_2UnavailBest <- loess(Observed ~ Predicted, data=plotfr_2UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailBest <- predict(loessMod_2UnavailBest)
lines(smoothed_2UnavailBest, x=plotfr_2UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_2UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_2UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailWorst <- predict(loessMod_2UnavailWorst)
lines(smoothed_2UnavailWorst, x=plotfr_2UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 2 predictors unavailable", "Worst 2 predictors unavailable"), bty='n')
dev.off()

#Base
plotfr
#2UnavailBest
plotfr_2UnavailBest
#2UnavailWorst
plotfr_2UnavailWorst

rm(list = ls())


###############################################################################################################################################
################################ ADL with core and exit interviews 1 unavailable predictor ####################################################

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


####################################
#Fit 16-predictor model with weights

mydata2<- subset(mydata, subgroup_adldep==1)

#Import data with CumInc for ADL at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,2]
class(cif10_adl)

####### Calibration plot
xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100
plotfr

rm(list = c('cif10_adl', 'calfitadlw10', 'calfitadlw10bars', 'xs'))


####################################
#15-variable best model missing 1 variable (R5EATA): lowest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10_1UnavailBest.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,3]

xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr_1UnavailBest<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr_1UnavailBest) <- c('Predicted', 'Observed')
plotfr_1UnavailBest$Predicted<-plotfr_1UnavailBest$Predicted*100
plotfr_1UnavailBest$Observed<-plotfr_1UnavailBest$Observed*100

rm(list = c('cif10_adl', 'calfitadlw10', 'calfitadlw10bars', 'xs'))


####################################
#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10_1UnavailWorst.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,3]

xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr_1UnavailWorst<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr_1UnavailWorst) <- c('Predicted', 'Observed')
plotfr_1UnavailWorst$Predicted<-plotfr_1UnavailWorst$Predicted*100
plotfr_1UnavailWorst$Observed<-plotfr_1UnavailWorst$Observed*100

rm(list = c('cif10_adl', 'calfitadlw10', 'calfitadlw10bars', 'xs'))

####################################
# Make plot from extracted data
par(las=1)

#Plot in png format
png('ADLBaseBestWorst_1unavail_R.png', width=4.5, height=3.25, units="in", res=800, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black",lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_1UnavailBest$Predicted, plotfr_1UnavailBest$Observed, col="blue", pch=20)
loessMod_1UnavailBest <- loess(Observed ~ Predicted, data=plotfr_1UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailBest <- predict(loessMod_1UnavailBest)
lines(smoothed_1UnavailBest, x=plotfr_1UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_1UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_1UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailWorst <- predict(loessMod_1UnavailWorst)
lines(smoothed_1UnavailWorst, x=plotfr_1UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 1 predictor unavailable", "Worst 1 predictor unavailable"), bty='n')
dev.off()


#Plot in pdf format
pdf('ADLBaseBestWorst_1unavail_R.pdf', width=4.88, height=3.25, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black", lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_1UnavailBest$Predicted, plotfr_1UnavailBest$Observed, col="blue", pch=20)
loessMod_1UnavailBest <- loess(Observed ~ Predicted, data=plotfr_1UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailBest <- predict(loessMod_1UnavailBest)
lines(smoothed_1UnavailBest, x=plotfr_1UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_1UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_1UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailWorst <- predict(loessMod_1UnavailWorst)
lines(smoothed_1UnavailWorst, x=plotfr_1UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 1 predictor unavailable", "Worst 1 predictor unavailable"), bty='n')
dev.off()

#Base
plotfr
#1UnavailBest
plotfr_1UnavailBest
#1UnavailWorst
plotfr_1UnavailWorst

rm(list = c('plotfr_1UnavailBest', 'plotfr_1UnavailWorst', 'loessMod', 'loessMod_1UnavailBest', 'loessMod_1UnavailWorst','smoothed','smoothed_1UnavailBest', 'smoothed_1UnavailWorst'))


###############################################################################################################################################
################################ ADL with core and exit interviews 2 unavailable predictors ###################################################


####################################
#14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10_2UnavailBest.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,3]

xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr_2UnavailBest<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr_2UnavailBest) <- c('Predicted', 'Observed')
plotfr_2UnavailBest$Predicted<-plotfr_2UnavailBest$Predicted*100
plotfr_2UnavailBest$Observed<-plotfr_2UnavailBest$Observed*100

rm(list = c('cif10_adl', 'calfitadlw10', 'calfitadlw10bars', 'xs'))


####################################
#14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10_2UnavailWorst.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,3]

xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr_2UnavailWorst<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr_2UnavailWorst) <- c('Predicted', 'Observed')
plotfr_2UnavailWorst$Predicted<-plotfr_2UnavailWorst$Predicted*100
plotfr_2UnavailWorst$Observed<-plotfr_2UnavailWorst$Observed*100

rm(list = c('cif10_adl', 'calfitadlw10', 'calfitadlw10bars', 'xs'))

####################################
# Make plot from extracted data
par(las=1)

#Plot in png format
png('ADLBaseBestWorst_2Unavail_R.png', width=4.5, height=3.25, units="in", res=800, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black",lwd=.5)

#14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC
#points(plotfr_2UnavailBest$Predicted, plotfr_2UnavailBest$Observed, col="blue", pch=20)
loessMod_2UnavailBest <- loess(Observed ~ Predicted, data=plotfr_2UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailBest <- predict(loessMod_2UnavailBest)
lines(smoothed_2UnavailBest, x=plotfr_2UnavailBest$Predicted, col="blue", lty=3)

#14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC
loessMod_2UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_2UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailWorst <- predict(loessMod_2UnavailWorst)
lines(smoothed_2UnavailWorst, x=plotfr_2UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 2 predictors unavailable", "Worst 2 predictors unavailable"), bty='n')
dev.off()


#Plot in pdf format
pdf('ADLBaseBestWorst_2Unavail_R.pdf', width=4.88, height=3.25, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black", lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_2UnavailBest$Predicted, plotfr_2UnavailBest$Observed, col="blue", pch=20)
loessMod_2UnavailBest <- loess(Observed ~ Predicted, data=plotfr_2UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailBest <- predict(loessMod_2UnavailBest)
lines(smoothed_2UnavailBest, x=plotfr_2UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_2UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_2UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailWorst <- predict(loessMod_2UnavailWorst)
lines(smoothed_2UnavailWorst, x=plotfr_2UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 2 predictors unavailable", "Worst 2 predictors unavailable"), bty='n')
dev.off()

#Base
plotfr
#2UnavailBest
plotfr_2UnavailBest
#2UnavailWorst
plotfr_2UnavailWorst

rm(list = ls())

###############################################################################################################################################
################################ Walk with core and exit interviews 1 unavailable predictor ####################################################

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


####################################
#Fit 16-predictor model with weights

mydata2<- subset(mydata, subgroup_walkdep==1)

#Import data with CumInc for Walk at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,2]
class(cif10_walk)

####### Calibration plot
xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr) <- c('Predicted', 'Observed')
plotfr$Predicted<-plotfr$Predicted*100
plotfr$Observed<-plotfr$Observed*100
plotfr

rm(list = c('cif10_walk', 'calfitwalkw10', 'calfitwalkw10bars', 'xs'))


####################################
#15-variable best model missing 1 variable (R5EATA): lowest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10_1UnavailBest.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,3]

xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr_1UnavailBest<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr_1UnavailBest) <- c('Predicted', 'Observed')
plotfr_1UnavailBest$Predicted<-plotfr_1UnavailBest$Predicted*100
plotfr_1UnavailBest$Observed<-plotfr_1UnavailBest$Observed*100

rm(list = c('cif10_walk', 'calfitwalkw10', 'calfitwalkw10bars', 'xs'))


####################################
#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10_1UnavailWorst.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,3]

xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr_1UnavailWorst<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr_1UnavailWorst) <- c('Predicted', 'Observed')
plotfr_1UnavailWorst$Predicted<-plotfr_1UnavailWorst$Predicted*100
plotfr_1UnavailWorst$Observed<-plotfr_1UnavailWorst$Observed*100

rm(list = c('cif10_walk', 'calfitwalkw10', 'calfitwalkw10bars', 'xs'))

####################################
# Make plot from extracted data
par(las=1)

#Plot in png format
png('WalkBaseBestWorst_1unavail_R.png', width=4.5, height=3.25, units="in", res=800, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black",lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_1UnavailBest$Predicted, plotfr_1UnavailBest$Observed, col="blue", pch=20)
loessMod_1UnavailBest <- loess(Observed ~ Predicted, data=plotfr_1UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailBest <- predict(loessMod_1UnavailBest)
lines(smoothed_1UnavailBest, x=plotfr_1UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_1UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_1UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailWorst <- predict(loessMod_1UnavailWorst)
lines(smoothed_1UnavailWorst, x=plotfr_1UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 1 predictor unavailable", "Worst 1 predictor unavailable"), bty='n')
dev.off()


#Plot in pdf format
pdf('WalkBaseBestWorst_1unavail_R.pdf', width=4.88, height=3.25, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black", lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_1UnavailBest$Predicted, plotfr_1UnavailBest$Observed, col="blue", pch=20)
loessMod_1UnavailBest <- loess(Observed ~ Predicted, data=plotfr_1UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailBest <- predict(loessMod_1UnavailBest)
lines(smoothed_1UnavailBest, x=plotfr_1UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_1UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_1UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_1UnavailWorst <- predict(loessMod_1UnavailWorst)
lines(smoothed_1UnavailWorst, x=plotfr_1UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 1 predictor unavailable", "Worst 1 predictor unavailable"), bty='n')
dev.off()

#Base
plotfr
#1UnavailBest
plotfr_1UnavailBest
#1UnavailWorst
plotfr_1UnavailWorst

rm(list = c('plotfr_1UnavailBest', 'plotfr_1UnavailWorst', 'loessMod', 'loessMod_1UnavailBest', 'loessMod_1UnavailWorst','smoothed','smoothed_1UnavailBest', 'smoothed_1UnavailWorst'))


###############################################################################################################################################
################################ Walk with core and exit interviews 2 unavailable predictors ###################################################


####################################
#14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10_2UnavailBest.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,3]

xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr_2UnavailBest<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr_2UnavailBest) <- c('Predicted', 'Observed')
plotfr_2UnavailBest$Predicted<-plotfr_2UnavailBest$Predicted*100
plotfr_2UnavailBest$Observed<-plotfr_2UnavailBest$Observed*100

rm(list = c('cif10_walk', 'calfitwalkw10', 'calfitwalkw10bars', 'xs'))


####################################
#14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC

#Import data with predicted event probabilities for death at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10_2UnavailWorst.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,3]

xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

####### Calibration bar plot
calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr_2UnavailWorst<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr_2UnavailWorst) <- c('Predicted', 'Observed')
plotfr_2UnavailWorst$Predicted<-plotfr_2UnavailWorst$Predicted*100
plotfr_2UnavailWorst$Observed<-plotfr_2UnavailWorst$Observed*100

rm(list = c('cif10_walk', 'calfitwalkw10', 'calfitwalkw10bars', 'xs'))

####################################
# Make plot from extracted data
par(las=1)

#Plot in png format
png('WalkBaseBestWorst_2Unavail_R.png', width=4.5, height=3.25, units="in", res=800, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black",lwd=.5)

#14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC
#points(plotfr_2UnavailBest$Predicted, plotfr_2UnavailBest$Observed, col="blue", pch=20)
loessMod_2UnavailBest <- loess(Observed ~ Predicted, data=plotfr_2UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailBest <- predict(loessMod_2UnavailBest)
lines(smoothed_2UnavailBest, x=plotfr_2UnavailBest$Predicted, col="blue", lty=3)

#14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC
loessMod_2UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_2UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailWorst <- predict(loessMod_2UnavailWorst)
lines(smoothed_2UnavailWorst, x=plotfr_2UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 2 predictors unavailable", "Worst 2 predictors unavailable"), bty='n')
dev.off()


#Plot in pdf format
pdf('WalkBaseBestWorst_2Unavail_R.pdf', width=4.88, height=3.25, pointsize=6)
#16-predictor model
plot(plotfr$Predicted, plotfr$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch='')
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
lines(smoothed, x=plotfr$Predicted, col="black", lwd=.5)

#15-variable best model missing 1 variable (R5EATA): lowest average TCIC
#points(plotfr_2UnavailBest$Predicted, plotfr_2UnavailBest$Observed, col="blue", pch=20)
loessMod_2UnavailBest <- loess(Observed ~ Predicted, data=plotfr_2UnavailBest, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailBest <- predict(loessMod_2UnavailBest)
lines(smoothed_2UnavailBest, x=plotfr_2UnavailBest$Predicted, col="blue", lty=3)

#15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC
loessMod_2UnavailWorst <- loess(Observed ~ Predicted, data=plotfr_2UnavailWorst, span=0.50)
# as the span increases, the smoothing of the curve also increases
smoothed_2UnavailWorst <- predict(loessMod_2UnavailWorst)
lines(smoothed_2UnavailWorst, x=plotfr_2UnavailWorst$Predicted, col="red", lty=2)

legend('bottomright', lty=c(1,3,2), col=c("black", "blue", "red"), legend=c("Base", "Best 2 predictors unavailable", "Worst 2 predictors unavailable"), bty='n')
dev.off()

#Base
plotfr
#2UnavailBest
plotfr_2UnavailBest
#2UnavailWorst
plotfr_2UnavailWorst

rm(list = ls())

