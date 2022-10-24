# Grisell Workstation
install.packages("pec", lib="C:/Program Files/R/R-4.1.1/library")
install.packages("coxphw", lib="C:/Program Files/R/R-4.1.1/library")
install.packages("cmprsk", lib="C:/Program Files/R/R-4.1.1/library")
install.packages("timereg", lib="C:/Program Files/R/R-4.1.1/library")
install.packages("haven",lib="C:/Program Files/R/R-4.1.1/library")
install.packages("survey",lib="C:/Program Files/R/R-4.1.1/library")
install.packages("dplyr",lib="C:/Program Files/R/R-4.1.1/library")

# setwd(/path/Rfiles")

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
library(ggplot2)


packageDescription("pec", fields = "Version")
# [1] "2022.03.06"

packageDescription("riskRegression", fields = "Version")
# [1] "2022.03.22"

# R version 4.1.1 (2021-08-10) -- "Kick Things"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)



###############################################################################################################################################
######### DEATH ##########


#############################################
#Fit model with weights

#Import data with SurvProb for death at 10 yrs
cum10_death <- read.csv('SurvPredDth10_ageonly.csv')
str(cum10_death)
cum10_death<-cum10_death[,4]

####### Calibration plot
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
plotfr_death<-as.data.frame(calfitdthw10bars$plotFrames)
colnames(plotfr_death) <- c('Predicted', 'Observed')
plotfr_death$Predicted<-plotfr_death$Predicted*100
plotfr_death$Observed<-plotfr_death$Observed*100

loessMod <- loess(Observed ~ Predicted, data=plotfr_death, span=0.50)
plotfr_death$smoothed <- predict(loessMod)

plotfr_death

# Make plot from extracted data
par(las=1)
plot(plotfr_death$Predicted, plotfr_death$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')
#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
# as the span increases, the smoothing of the curve also increases
lines(plotfr_death$smoothed, x=plotfr_death$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')
plotfr_death

rm(list = c('xs','loessMod', 'calfitdthw10', 'calfitdthw10bars', 'cum10_death'))

###############################################################################################################################################
######### ADL with core and exit interviews ##########
mydata2<- subset(mydata, subgroup_adldep==1)


#############################################
#Fit model with weights

#Import data with CumInc for ADL at 10 yrs
cif10_adl <- read.csv('CumIncPredADL10_ageonly.csv')
str(cif10_adl)
cif10_adl<-cif10_adl[,3]
class(cif10_adl)

xs=Score(list("Weighted FG ADL t=10"=cif10_adl),
         formula=Hist(time_adldepdth2,status_adldepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitadlw10<-plotCalibration(xs,cens.method="local") 

calfitadlw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 


# Extract data for plot
names(calfitadlw10bars)
calfitadlw10bars$plotFrames
plotfr_adl<-as.data.frame(calfitadlw10bars$plotFrames)
colnames(plotfr_adl) <- c('Predicted', 'Observed')
plotfr_adl$Predicted<-plotfr_adl$Predicted*100
plotfr_adl$Observed<-plotfr_adl$Observed*100

#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
# as the span increases, the smoothing of the curve also increases
loessMod <- loess(Observed ~ Predicted, data=plotfr_adl, span=0.50)
plotfr_adl$smoothed <- predict(loessMod)

# Make plot from extracted data
par(las=1)
plot(plotfr_adl$Predicted, plotfr_adl$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')

lines(plotfr_adl$smoothed, x=plotfr_adl$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs','mydata2', 'loessMod', 'calfitadlw10', 'calfitadlw10bars', 'cif10_adl'))


###############################################################################################################################################
######### Walk with core and exit interviews ##########

mydata2<- subset(mydata, subgroup_walkdep==1)


#############################################
#Fit model with weights

#Import data with CumInc for Walk at 10 yrs
cif10_walk <- read.csv('CumIncPredWalk10_ageonly.csv')
str(cif10_walk)
cif10_walk<-cif10_walk[,3]
class(cif10_walk)

xs=Score(list("Weighted FG Walk t=10"=cif10_walk),
         formula=Hist(time_walkdepdth2,status_walkdepdth2)~1,data=mydata2,conf.int=TRUE,times=10, plots="cal", se.fit=1)
xs
calfitwalkw10<-plotCalibration(xs,cens.method="local") 

calfitwalkw10bars<-plotCalibration(xs,cens.method="local",bars=TRUE,hanging=FALSE, show.frequencies=TRUE) 

# Extract data for plot
names(calfitwalkw10bars)
calfitwalkw10bars$plotFrames
plotfr_walk<-as.data.frame(calfitwalkw10bars$plotFrames)
colnames(plotfr_walk) <- c('Predicted', 'Observed')
plotfr_walk$Predicted<-plotfr_walk$Predicted*100
plotfr_walk$Observed<-plotfr_walk$Observed*100

#lines(lowess(plotfr, f = .2, iter = 3, delta = 0.01 * diff(range(plotfr$Predicted))), col="blue")
# f: the smoother span. This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness
# loess, a newer formula based version of lowess (with different defaults!).
loessMod <- loess(Observed ~ Predicted, data=plotfr_walk, span=0.50)
# as the span increases, the smoothing of the curve also increases
plotfr_walk$smoothed <- predict(loessMod)


# Make plot from extracted data
par(las=1)
plot(plotfr_walk$Predicted, plotfr_walk$Observed, xlab='', ylab='',xaxt = "n",yaxt = "n", xlim=c(0,100), ylim=c(0,100), pch=20)
#lines(plotfr$Predicted, plotfr$Observed,lwd=0.6)
axis(1, at=seq(0, 100, by =25))
axis(2, at=seq(0, 100, by =25))
abline(a=0,b=1, lwd=0.6, col="gray")
title(xlab='Predicted risk ', ylab='Estimated actual risk')

lines(plotfr_walk$smoothed, x=plotfr_walk$Predicted, col="black")
legend('bottomright', lty=1, col="black", legend="Loess fit", bty='n')

rm(list = c('xs','mydata2', 'loessMod', 'calfitwalkw10', 'calfitwalkw10bars', 'cif10_walk'))


###############################################################################################################################################
######### Combine 3 plots in one with ggplot ##########

plotfr_adl$outcome<-'adl'
plotfr_adl$panel<-'A'

plotfr_walk$outcome<-'walk'
plotfr_walk$panel<-'B'

plotfr_death$outcome<-'death'
plotfr_death$panel<-'C'

plotfr<-rbind(plotfr_adl, plotfr_walk, plotfr_death)

str(plotfr)

# save(plotfr,file="CalPlotDataAgeOnly_gdr_20220510.Rdata")

load("CalPlotDataAgeOnly_gdr_20220510.Rdata")

# png('Fig1.CalibrationPlotsAgeOnly.png')

ggplot(data = plotfr, aes(x = Predicted, y = smoothed, color=outcome)) + 
  geom_abline(intercept=0, slope=1, lwd=0.3, col="gray") +
  geom_line(lwd=0.6) +
  geom_point(data = plotfr, aes(x = Predicted, y = Observed, color=outcome)) +
  theme_bw() +
  ## removed box around the facet title;decrease the margin between panels and text;
  theme(strip.background = element_blank())+
  theme(panel.spacing = unit(0.2, "lines"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~ panel, ncol=2,  strip.position = "top") +
  theme(strip.text = element_text(size = 14, face = "bold", hjust=0)) +
  labs(x="Predicted risk", y="Estimated actual risk") +  
  theme(axis.title.x = element_text(face="bold", vjust=-0.35, size=14)) +
  theme(axis.title.y = element_text(face="bold", vjust=0.35, size=14)) +
  theme(axis.text = element_text(size = 14)) + # Change Font Size of Axis Text
  scale_x_continuous(limits=c(0, 100), breaks=seq(0,100,by=25)) +
  scale_y_continuous(limits=c(0, 100), breaks=seq(0,100,by=25)) +
  scale_color_manual(name = "Outcome:",  
                   labels = c("ADL disability", 
                              "Walking disability", 
                              "Death"),  
                   values = c("adl"='cornflowerblue', "walk"='darkorange', "death"='darkgreen') ) +
  theme(legend.position = c(0.8,0.1), 
        legend.title = element_text(size = 14, face = "bold", vjust = 1), 
        legend.key.size = unit(1.2, "lines"),
        legend.text=element_text(size = 14) )   


# dev.off()

?ggsave
ggsave(filename="Fig1.CalibrationPlotsAgeOnly.png", width=10, height=10, units="in", device="png", dpi=300)
ggsave(filename="Fig1.CalibrationPlotsAgeOnly.pdf", width=10, height=10, units="in", device="pdf", dpi=300)
