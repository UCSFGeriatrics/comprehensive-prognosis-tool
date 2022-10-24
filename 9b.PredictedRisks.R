R.version$version.string
# "R version 4.2.1 (2022-06-23 ucrt)"

# install.packages("Hmisc")
# install.packages("haven")
# install.packages('labelled')

library(Hmisc)
library(haven)
library(labelled)


packageDescription("Hmisc", fields = "Version")
# ver "4.7-1"
packageDescription("haven", fields = "Version")
# ver "2.5.1"
packageDescription("labelled", fields = "Version")
# ver "2.10.0"


#################################################### Prepare data #################################################### 

#Import data
# setwd("path/sasdata")
# mydata <- read_sas("clinicalhrsdeident_gdr_20210624.sas7bdat")
# 
# #Keep only variables of interest
# mydata<- subset(mydata, select = c(newid,status_adldepdth2, status_walkdepdth2, death, time_adldepdth2, time_walkdepdth2, time2death,
#                                    subgroup_adldep, subgroup_walkdep, R5WTRESP, RAEHSAMP, RAESTRAT, ageint5,
#                                    age3_lin, age3_sp1, age3_sp2, gender, R5EATA, R5CANCRE, R5DIABE, heart3g, R5HIBPE, R5MEALSA, R5MONEYA,
#                                    h5lvalone, R5LUNGE, R5PUSHA, R5WALKSA, R5BMI, bmi3_lin, bmi3_sp1, bmi3_sp2, smoking3g, stroke3g) )
# contents(mydata)
# 
# setwd("path/Rfiles")
# save(mydata,file="clinicalhrsdeident_gdr_20221008.Rdata")


#Load data
load("clinicalhrsdeident_gdr_20221008.Rdata")

#Check type of variables
str(mydata)

#Check columns classes
sapply(mydata, class)

#Create dummy variables for categorical variables

##Check if there are missing values: there are no missing values in predictors
describe(mydata)
contents(mydata)

mydata$heart1<-ifelse(mydata$heart3g==1, 1, 0)
mydata$heart2<-ifelse(mydata$heart3g==2, 1, 0)
table(mydata$heart3g,mydata$heart1)   
table(mydata$heart3g,mydata$heart2)   

mydata$smoking1<-ifelse(mydata$smoking3g==1, 1, 0)
mydata$smoking2<-ifelse(mydata$smoking3g==2, 1, 0)
table(mydata$smoking3g,mydata$smoking1)   
table(mydata$smoking3g,mydata$smoking2)   

mydata$stroke1<-ifelse(mydata$stroke3g==1, 1, 0)
mydata$stroke2<-ifelse(mydata$stroke3g==2, 1, 0)
table(mydata$stroke3g,mydata$stroke1)   
table(mydata$stroke3g,mydata$stroke2)   


#Create splines
##Note: the dataset of the manuscript already has the variables for the splines (age3_lin,age3_sp1,age3_sp2,bmi3_lin,bmi3_sp1,bmi3_sp2),
##but there are re-created below:

#Drop splines variables so they can be re-created
mydata2<- subset(mydata, select = -c(age3_lin,age3_sp1,age3_sp2,bmi3_lin,bmi3_sp1,bmi3_sp2) )
contents(mydata2)

#Create 3 restricted cubic splines for age using default knots in manuscript data

#Age
spline_age<-rcspline.eval(mydata$ageint5, knots=c(70, 75, 79, 89), inclx=TRUE) 
colnames(spline_age)<-c('age3_lin', 'age3_sp1', 'age3_sp2')
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots-1= 4-1=3splines (this includes linear time as the 1st column)

attributes(spline_age)
# $dim
# [1] 6646    3
# 
# $knots
# [1] 70 75 79 89

class(spline_age)

mydata2$age3_lin<-spline_age[,'age3_lin']
mydata2$age3_sp1<-spline_age[,'age3_sp1']
mydata2$age3_sp2<-spline_age[,'age3_sp2']

#Check that new variables are the same as old variables
summary(mydata2$age3_lin)
summary(mydata$age3_lin)

summary(mydata2$age3_sp1)
summary(mydata$age3_sp1)

summary(mydata2$age3_sp2)
summary(mydata$age3_sp2)

#Create 3 restricted cubic splines for BMI using default knots in manuscript data

#BMI
spline_bmi<-rcspline.eval(mydata$R5BMI, knots=c(19.1, 23.8, 27.1, 34.275), inclx=TRUE) 
colnames(spline_bmi)<-c('bmi3_lin', 'bmi3_sp1', 'bmi3_sp2')

attributes(spline_bmi)

# $knots
# [1] 19.100 23.800 27.100 34.275

mydata2$bmi3_lin<-spline_bmi[,'bmi3_lin']
mydata2$bmi3_sp1<-spline_bmi[,'bmi3_sp1']
mydata2$bmi3_sp2<-spline_bmi[,'bmi3_sp2']

#Check that new variables are the same as old variables
summary(mydata2$bmi3_lin)
summary(mydata$bmi3_lin)

summary(mydata2$bmi3_sp1)
summary(mydata$bmi3_sp1)

summary(mydata2$bmi3_sp2)
summary(mydata$bmi3_sp2)

rm(spline_age,spline_bmi)

#################################################### Get predicted Risks ADL outcome #################################################### 

#1)Import Baseline Survival functions
basesurv_adl<-read.csv('Weightbasesurv_ADL.csv')
#rows are the models, eg. final model, age-only, unavailable predictors models
#columns have baseline survival function from 1 to 19 

##Select 1st row (final model, no missing predictor) and columns for 5, 10, 14 years: So_ADL_t5, So_ADL_t10,  So_ADL_t14
basesurv_adl2<-subset(basesurv_adl, DELEVAR=="",select = c(modelnum,So_ADL_t5, So_ADL_t10,  So_ADL_t14) )


##Add baseline survival functions to mydata2
mydata2$So_ADL_t5<-basesurv_adl2$So_ADL_t5
mydata2$So_ADL_t10<-basesurv_adl2$So_ADL_t10
mydata2$So_ADL_t14<-basesurv_adl2$So_ADL_t14

contents(mydata2)
describe(mydata2)


#2)Import beta estimates
betas_adl<-read.csv('Weightbetasurv_ADL.csv')
#rows are the models, eg. final model, age-only, unavailable predictors models
#columns have beta estimates for each parameter in the model

#Select 1st row (final model, no missing predictor)
betas_adl2<-subset(betas_adl, DELEVAR=="" )
colnames(betas_adl2)

#Add betas estimates to mydata2
mydata2$age3_lin_beta<-betas_adl2$age3_lin_beta
mydata2$age3_sp1_beta<-betas_adl2$age3_sp1_beta
mydata2$age3_sp2_beta<-betas_adl2$age3_sp2_beta
mydata2$gender_beta<-betas_adl2$gender_beta
mydata2$R5EATA_beta<-betas_adl2$R5EATA_beta
mydata2$R5CANCRE_beta<-betas_adl2$R5CANCRE_beta
mydata2$R5DIABE_beta<-betas_adl2$R5DIABE_beta
mydata2$heart1_beta<-betas_adl2$heart1_beta
mydata2$heart2_beta<-betas_adl2$heart2_beta
mydata2$R5HIBPE_beta<-betas_adl2$R5HIBPE_beta
mydata2$R5MEALSA_beta<-betas_adl2$R5MEALSA_beta
mydata2$R5MONEYA_beta<-betas_adl2$R5MONEYA_beta
mydata2$h5lvalone_beta<-betas_adl2$h5lvalone_beta
mydata2$R5LUNGE_beta<-betas_adl2$R5LUNGE_beta
mydata2$R5PUSHA_beta<-betas_adl2$R5PUSHA_beta
mydata2$R5WALKSA_beta<-betas_adl2$R5WALKSA_beta
mydata2$bmi3_lin_beta<-betas_adl2$bmi3_lin_beta
mydata2$bmi3_sp1_beta<-betas_adl2$bmi3_sp1_beta
mydata2$bmi3_sp2_beta<-betas_adl2$bmi3_sp2_beta
mydata2$smoking1_beta<-betas_adl2$smoking1_beta
mydata2$smoking2_beta<-betas_adl2$smoking2_beta
mydata2$stroke1_beta<-betas_adl2$stroke1_beta
mydata2$stroke2_beta<-betas_adl2$stroke2_beta

contents(mydata2)
describe(mydata2)


#3)Compute xbeta
mydata2$xbeta<-mydata2$age3_lin*mydata2$age3_lin_beta+mydata2$age3_sp1*mydata2$age3_sp1_beta+mydata2$age3_sp2*mydata2$age3_sp2_beta+
               mydata2$gender*mydata2$gender_beta+mydata2$R5EATA*mydata2$R5EATA_beta+mydata2$R5CANCRE*mydata2$R5CANCRE_beta+
               mydata2$R5DIABE*mydata2$R5DIABE_beta+mydata2$heart1*mydata2$heart1_beta+mydata2$heart2*mydata2$heart2_beta+
               mydata2$R5HIBPE*mydata2$R5HIBPE_beta+mydata2$R5MEALSA*mydata2$R5MEALSA_beta+mydata2$R5MONEYA*mydata2$R5MONEYA_beta+
               mydata2$h5lvalone*mydata2$h5lvalone_beta+mydata2$R5LUNGE*mydata2$R5LUNGE_beta+mydata2$R5PUSHA*mydata2$R5PUSHA_beta+
               mydata2$R5WALKSA*mydata2$R5WALKSA_beta+mydata2$bmi3_lin*mydata2$bmi3_lin_beta+mydata2$bmi3_sp1*mydata2$bmi3_sp1_beta+
               mydata2$bmi3_sp2*mydata2$bmi3_sp2_beta+mydata2$smoking1*mydata2$smoking1_beta+mydata2$smoking2*mydata2$smoking2_beta+
               mydata2$stroke1*mydata2$stroke1_beta+mydata2$stroke2*mydata2$stroke2_beta

var_label(mydata2$xbeta) <- NULL #remove label

#4)Compute S_ADL_t5, S_ADL_t10, S_ADL_t14
mydata2$S_ADL_t5<-mydata2$So_ADL_t5^(exp(mydata2$xbeta))
mydata2$S_ADL_t10<-mydata2$So_ADL_t10^(exp(mydata2$xbeta))
mydata2$S_ADL_t14<-mydata2$So_ADL_t14^(exp(mydata2$xbeta))

#5)Compute Risk_ADL_t5, Risk_ADL_t10, Risk_ADL_t14
mydata2$Risk_ADL_t5<-1-mydata2$S_ADL_t5
mydata2$Risk_ADL_t10<-1-mydata2$S_ADL_t10
mydata2$Risk_ADL_t14<-1-mydata2$S_ADL_t14

contents(mydata2)
describe(mydata2)

#Check some Ids
check<-subset(mydata2, newid==16 | newid==30 | newid==847,select = c(newid,Risk_ADL_t5,Risk_ADL_t10,Risk_ADL_t14) )
check

rm(check,basesurv_adl,basesurv_adl2,betas_adl,betas_adl2)


#################################################### Get predicted Risks WALK outcome #################################################### 

#Remove beta and xbeta estimates for ADL outcome
mydata2<-subset(mydata2, select = -c(age3_lin_beta, age3_sp1_beta, age3_sp2_beta, gender_beta,   
                                     R5EATA_beta, R5CANCRE_beta, R5DIABE_beta, heart1_beta, heart2_beta, R5HIBPE_beta, R5MEALSA_beta, R5MONEYA_beta, 
                                     h5lvalone_beta, R5LUNGE_beta, R5PUSHA_beta, R5WALKSA_beta, bmi3_lin_beta, bmi3_sp1_beta, bmi3_sp2_beta, smoking1_beta, 
                                     smoking2_beta, stroke1_beta, stroke2_beta, xbeta) )


contents(mydata2)


#1)Import Baseline Survival functions
basesurv_walk<-read.csv('Weightbasesurv_WALK.csv')
#rows are the models, eg. final model, age-only, unavailable predictors models
#columns have baseline survival function from 1 to 19 

##Select 1st row (final model, no missing predictor) and columns for 5, 10, 14 years: So_WALK_t5, So_WALK_t10,  So_WALK_t14
basesurv_walk2<-subset(basesurv_walk, DELEVAR=="",select = c(modelnum,So_WALK_t5, So_WALK_t10,  So_WALK_t14) )


##Add baseline survival functions to mydata2
mydata2$So_WALK_t5<-basesurv_walk2$So_WALK_t5
mydata2$So_WALK_t10<-basesurv_walk2$So_WALK_t10
mydata2$So_WALK_t14<-basesurv_walk2$So_WALK_t14

contents(mydata2)
describe(mydata2)


#2)Import beta estimates
betas_walk<-read.csv('Weightbetasurv_WALK.csv')
#rows are the models, eg. final model, age-only, unavailable predictors models
#columns have beta estimates for each parameter in the model

#Select 1st row (final model, no missing predictor)
betas_walk2<-subset(betas_walk, DELEVAR=="" )
colnames(betas_walk2)

#Add betas estimates to mydata2
mydata2$age3_lin_beta<-betas_walk2$age3_lin_beta
mydata2$age3_sp1_beta<-betas_walk2$age3_sp1_beta
mydata2$age3_sp2_beta<-betas_walk2$age3_sp2_beta
mydata2$gender_beta<-betas_walk2$gender_beta
mydata2$R5EATA_beta<-betas_walk2$R5EATA_beta
mydata2$R5CANCRE_beta<-betas_walk2$R5CANCRE_beta
mydata2$R5DIABE_beta<-betas_walk2$R5DIABE_beta
mydata2$heart1_beta<-betas_walk2$heart1_beta
mydata2$heart2_beta<-betas_walk2$heart2_beta
mydata2$R5HIBPE_beta<-betas_walk2$R5HIBPE_beta
mydata2$R5MEALSA_beta<-betas_walk2$R5MEALSA_beta
mydata2$R5MONEYA_beta<-betas_walk2$R5MONEYA_beta
mydata2$h5lvalone_beta<-betas_walk2$h5lvalone_beta
mydata2$R5LUNGE_beta<-betas_walk2$R5LUNGE_beta
mydata2$R5PUSHA_beta<-betas_walk2$R5PUSHA_beta
mydata2$R5WALKSA_beta<-betas_walk2$R5WALKSA_beta
mydata2$bmi3_lin_beta<-betas_walk2$bmi3_lin_beta
mydata2$bmi3_sp1_beta<-betas_walk2$bmi3_sp1_beta
mydata2$bmi3_sp2_beta<-betas_walk2$bmi3_sp2_beta
mydata2$smoking1_beta<-betas_walk2$smoking1_beta
mydata2$smoking2_beta<-betas_walk2$smoking2_beta
mydata2$stroke1_beta<-betas_walk2$stroke1_beta
mydata2$stroke2_beta<-betas_walk2$stroke2_beta

contents(mydata2)
describe(mydata2)


#3)Compute xbeta
mydata2$xbeta<-mydata2$age3_lin*mydata2$age3_lin_beta+mydata2$age3_sp1*mydata2$age3_sp1_beta+mydata2$age3_sp2*mydata2$age3_sp2_beta+
  mydata2$gender*mydata2$gender_beta+mydata2$R5EATA*mydata2$R5EATA_beta+mydata2$R5CANCRE*mydata2$R5CANCRE_beta+
  mydata2$R5DIABE*mydata2$R5DIABE_beta+mydata2$heart1*mydata2$heart1_beta+mydata2$heart2*mydata2$heart2_beta+
  mydata2$R5HIBPE*mydata2$R5HIBPE_beta+mydata2$R5MEALSA*mydata2$R5MEALSA_beta+mydata2$R5MONEYA*mydata2$R5MONEYA_beta+
  mydata2$h5lvalone*mydata2$h5lvalone_beta+mydata2$R5LUNGE*mydata2$R5LUNGE_beta+mydata2$R5PUSHA*mydata2$R5PUSHA_beta+
  mydata2$R5WALKSA*mydata2$R5WALKSA_beta+mydata2$bmi3_lin*mydata2$bmi3_lin_beta+mydata2$bmi3_sp1*mydata2$bmi3_sp1_beta+
  mydata2$bmi3_sp2*mydata2$bmi3_sp2_beta+mydata2$smoking1*mydata2$smoking1_beta+mydata2$smoking2*mydata2$smoking2_beta+
  mydata2$stroke1*mydata2$stroke1_beta+mydata2$stroke2*mydata2$stroke2_beta

var_label(mydata2$xbeta) <- NULL #remove label

#4)Compute S_WALK_t5, S_WALK_t10, S_WALK_t14
mydata2$S_WALK_t5<-mydata2$So_WALK_t5^(exp(mydata2$xbeta))
mydata2$S_WALK_t10<-mydata2$So_WALK_t10^(exp(mydata2$xbeta))
mydata2$S_WALK_t14<-mydata2$So_WALK_t14^(exp(mydata2$xbeta))

#5)Compute Risk_WALK_t5, Risk_WALK_t10, Risk_WALK_t14
mydata2$Risk_WALK_t5<-1-mydata2$S_WALK_t5
mydata2$Risk_WALK_t10<-1-mydata2$S_WALK_t10
mydata2$Risk_WALK_t14<-1-mydata2$S_WALK_t14

contents(mydata2)
describe(mydata2)

#Check some Ids
check<-subset(mydata2, newid==16 | newid==30 | newid==847,select = c(newid,Risk_WALK_t5,Risk_WALK_t10,Risk_WALK_t14) )
check

rm(check,basesurv_walk,basesurv_walk2,betas_walk,betas_walk2)
contents(mydata2)


#################################################### Get predicted Risks DEATH outcome #################################################### 

#Remove beta and xbeta estimates for WALK outcome
mydata2<-subset(mydata2, select = -c(age3_lin_beta, age3_sp1_beta, age3_sp2_beta, gender_beta,   
                                     R5EATA_beta, R5CANCRE_beta, R5DIABE_beta, heart1_beta, heart2_beta, R5HIBPE_beta, R5MEALSA_beta, R5MONEYA_beta, 
                                     h5lvalone_beta, R5LUNGE_beta, R5PUSHA_beta, R5WALKSA_beta, bmi3_lin_beta, bmi3_sp1_beta, bmi3_sp2_beta, smoking1_beta, 
                                     smoking2_beta, stroke1_beta, stroke2_beta, xbeta) )



contents(mydata2)


#1)Import Baseline Survival functions
basesurv_death<-read.csv('Weightbasesurv_DEATH.csv')
#rows are the models, eg. final model, age-only, unavailable predictors models
#columns have baseline survival function from 1 to 19 

##Select 1st row (final model, no missing predictor) and columns for 5, 10, 14 years: So_DEATH_t5, So_DEATH_t10,  So_DEATH_t14
basesurv_death2<-subset(basesurv_death, DELEVAR=="",select = c(modelnum,So_DEATH_t5, So_DEATH_t10,  So_DEATH_t14) )


##Add baseline survival functions to mydata2
mydata2$So_DEATH_t5<-basesurv_death2$So_DEATH_t5
mydata2$So_DEATH_t10<-basesurv_death2$So_DEATH_t10
mydata2$So_DEATH_t14<-basesurv_death2$So_DEATH_t14

contents(mydata2)
describe(mydata2)


#2)Import beta estimates
betas_death<-read.csv('Weightbetasurv_DEATH.csv')
#rows are the models, eg. final model, age-only, unavailable predictors models
#columns have beta estimates for each parameter in the model

#Select 1st row (final model, no missing predictor)
betas_death2<-subset(betas_death, DELEVAR=="" )
colnames(betas_death2)

#Add betas estimates to mydata2
mydata2$age3_lin_beta<-betas_death2$age3_lin_beta
mydata2$age3_sp1_beta<-betas_death2$age3_sp1_beta
mydata2$age3_sp2_beta<-betas_death2$age3_sp2_beta
mydata2$gender_beta<-betas_death2$gender_beta
mydata2$R5EATA_beta<-betas_death2$R5EATA_beta
mydata2$R5CANCRE_beta<-betas_death2$R5CANCRE_beta
mydata2$R5DIABE_beta<-betas_death2$R5DIABE_beta
mydata2$heart1_beta<-betas_death2$heart1_beta
mydata2$heart2_beta<-betas_death2$heart2_beta
mydata2$R5HIBPE_beta<-betas_death2$R5HIBPE_beta
mydata2$R5MEALSA_beta<-betas_death2$R5MEALSA_beta
mydata2$R5MONEYA_beta<-betas_death2$R5MONEYA_beta
mydata2$h5lvalone_beta<-betas_death2$h5lvalone_beta
mydata2$R5LUNGE_beta<-betas_death2$R5LUNGE_beta
mydata2$R5PUSHA_beta<-betas_death2$R5PUSHA_beta
mydata2$R5WALKSA_beta<-betas_death2$R5WALKSA_beta
mydata2$bmi3_lin_beta<-betas_death2$bmi3_lin_beta
mydata2$bmi3_sp1_beta<-betas_death2$bmi3_sp1_beta
mydata2$bmi3_sp2_beta<-betas_death2$bmi3_sp2_beta
mydata2$smoking1_beta<-betas_death2$smoking1_beta
mydata2$smoking2_beta<-betas_death2$smoking2_beta
mydata2$stroke1_beta<-betas_death2$stroke1_beta
mydata2$stroke2_beta<-betas_death2$stroke2_beta

contents(mydata2)
describe(mydata2)


#3)Compute xbeta
mydata2$xbeta<-mydata2$age3_lin*mydata2$age3_lin_beta+mydata2$age3_sp1*mydata2$age3_sp1_beta+mydata2$age3_sp2*mydata2$age3_sp2_beta+
  mydata2$gender*mydata2$gender_beta+mydata2$R5EATA*mydata2$R5EATA_beta+mydata2$R5CANCRE*mydata2$R5CANCRE_beta+
  mydata2$R5DIABE*mydata2$R5DIABE_beta+mydata2$heart1*mydata2$heart1_beta+mydata2$heart2*mydata2$heart2_beta+
  mydata2$R5HIBPE*mydata2$R5HIBPE_beta+mydata2$R5MEALSA*mydata2$R5MEALSA_beta+mydata2$R5MONEYA*mydata2$R5MONEYA_beta+
  mydata2$h5lvalone*mydata2$h5lvalone_beta+mydata2$R5LUNGE*mydata2$R5LUNGE_beta+mydata2$R5PUSHA*mydata2$R5PUSHA_beta+
  mydata2$R5WALKSA*mydata2$R5WALKSA_beta+mydata2$bmi3_lin*mydata2$bmi3_lin_beta+mydata2$bmi3_sp1*mydata2$bmi3_sp1_beta+
  mydata2$bmi3_sp2*mydata2$bmi3_sp2_beta+mydata2$smoking1*mydata2$smoking1_beta+mydata2$smoking2*mydata2$smoking2_beta+
  mydata2$stroke1*mydata2$stroke1_beta+mydata2$stroke2*mydata2$stroke2_beta

var_label(mydata2$xbeta) <- NULL #remove label

#4)Compute S_DEATH_t5, S_DEATH_t10, S_DEATH_t14
mydata2$S_DEATH_t5<-mydata2$So_DEATH_t5^(exp(mydata2$xbeta))
mydata2$S_DEATH_t10<-mydata2$So_DEATH_t10^(exp(mydata2$xbeta))
mydata2$S_DEATH_t14<-mydata2$So_DEATH_t14^(exp(mydata2$xbeta))

#5)Compute Risk_DEATH_t5, Risk_DEATH_t10, Risk_DEATH_t14
mydata2$Risk_DEATH_t5<-1-mydata2$S_DEATH_t5
mydata2$Risk_DEATH_t10<-1-mydata2$S_DEATH_t10
mydata2$Risk_DEATH_t14<-1-mydata2$S_DEATH_t14

contents(mydata2)
describe(mydata2)

#Check some Ids
check<-subset(mydata2, newid==16 | newid==30 | newid==847,select = c(newid,Risk_DEATH_t5,Risk_DEATH_t10,Risk_DEATH_t14) )
check

rm(check,basesurv_death,basesurv_death2,betas_death,betas_death2)

mydata2<-subset(mydata2, select = -c(age3_lin_beta, age3_sp1_beta, age3_sp2_beta, gender_beta,   
                                     R5EATA_beta, R5CANCRE_beta, R5DIABE_beta, heart1_beta, heart2_beta, R5HIBPE_beta, R5MEALSA_beta, R5MONEYA_beta, 
                                     h5lvalone_beta, R5LUNGE_beta, R5PUSHA_beta, R5WALKSA_beta, bmi3_lin_beta, bmi3_sp1_beta, bmi3_sp2_beta, smoking1_beta, 
                                     smoking2_beta, stroke1_beta, stroke2_beta, xbeta) )

contents(mydata2)

#Check some Ids
check<-subset(mydata2, newid==16 | newid==30 | newid==847,select = c(newid,Risk_ADL_t5,Risk_ADL_t10,Risk_ADL_t14) )
check<-subset(mydata2, newid==16 | newid==30 | newid==847,select = c(newid,Risk_WALK_t5,Risk_WALK_t10,Risk_WALK_t14) )
check<-subset(mydata2, newid==16 | newid==30 | newid==847,select = c(newid,Risk_DEATH_t5,Risk_DEATH_t10,Risk_DEATH_t14) )

rm(check)

#Save predicted risks
mydata3<-subset(mydata2, select = c(newid,Risk_ADL_t5,Risk_ADL_t10,Risk_ADL_t14,Risk_WALK_t5,Risk_WALK_t10,Risk_WALK_t14,Risk_DEATH_t5,Risk_DEATH_t10,Risk_DEATH_t14) )
write.csv(mydata3, file='PredictedRisksResults.csv', row.names=FALSE)

