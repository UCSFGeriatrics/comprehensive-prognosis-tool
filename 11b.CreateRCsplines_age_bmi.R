# install.packages("Hmisc",lib="C:/Program Files/R/R-4.1.1/library")
library(Hmisc)

packageDescription("Hmisc", fields = "Version")
# ver "4.5-0"
R.version$version.string
# R.version$version.string

setwd("path/Rfiles")

# Create splines for age

#ageint5 in original data
# n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
# 6646        0       33    0.997    77.85    6.544       70       71       73       77       81       86       89 
# lowest :  70  71  72  73  74, highest:  98  99 100 101 102

age3_lin<-seq(from=70, to=100, by=1)
age3_lin

spline_age<-rcspline.eval(age3_lin, knots=c(70, 75, 79, 89), inclx=TRUE)
colnames(spline_age)<-c('age3_lin', 'age3_sp1', 'age3_sp2')
mydata<-as.data.frame(spline_age)
write.csv(mydata, file='splinesAge_gdr_20220413.csv', row.names=FALSE)

# Create splines for BMI

# R5BMI in original data
# n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
# 6646        0      296        1    25.85     5.15    19.10    20.35    22.70    25.40    28.30    31.90    34.27 
# lowest : 12.9 13.3 13.8 14.0 14.1, highest: 52.1 53.1 54.5 58.2 75.5

bmi3_lin<-seq(from=14, to=50, by=0.1)
bmi3_lin

spline_bmi<-rcspline.eval(bmi3_lin, knots=c(19.1, 23.8, 27.1, 34.275), inclx=TRUE)
colnames(spline_bmi)<-c('bmi3_lin', 'bmi3_sp1', 'bmi3_sp2')
mydata<-as.data.frame(spline_bmi)
write.csv(mydata, file='splinesBMI_gdr_20220413.csv', row.names=FALSE)
