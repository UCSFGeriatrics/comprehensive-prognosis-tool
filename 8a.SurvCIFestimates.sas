***********************************************************************************************************************************************************************************;
*Purpose: Compute weighted predicted cumulative incidence for calibration plots                                                                                                    ;                                     
*Statistician: Grisell Diaz-Ramirez																																				   ;
*Started: 2019.10.07 																																							   ;
*Finished: 2021.07.13                                                                                                                                                              ;
*Date modified: 2021.10.15                                                                                                                                                         ;
***********************************************************************************************************************************************************************************;

options MERGENOBY=ERROR;

libname trk 'path';
libname harmo 'path';
libname rand 'path';
proc format cntlin=rand.sasfmts; run;
proc format cntlin=harmo.formats; run;
libname doi 'path';
libname dod 'path';
libname dob 'path';
libname fat 'path';
libname clinical 'path';


/******************************************************************************************************************************************************************/
*ADL: Compute weighted predicted cumulative incidence for Calibration plots;

data dataclinical_gdr_20210624; 
 set clinical.clinicalHRSdeident_gdr_20210624 (keep=newid RAEHSAMP RAESTRAT clustervar R5WTRESP ageint5 age_lin age_sp1 age_sp2 age_sp3 age3_lin age3_sp1 age3_sp2
                                gender R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g bmi_lin bmi_sp1 bmi_sp2 bmi_sp3 bmi3_lin bmi3_sp1 bmi3_sp2
                                smoking3g stroke3g volunteer R5PROXY status_adldepdth2 status_walkdepdth2 death time_adldepdth2 time_walkdepdth2 time2death
								normwgt_adl normwgt_walk normwgt_death R5BMI
                                subgroup_adldep subgroup_walkdep rename=(subgroup_adldep=subgroup_adl subgroup_walkdep=subgroup_walk)); 
 subgroup_death=1;
run;
/* 6646 observations and 72 variables */

*Confirmed that unweighted estimnates and SEs in R are very similar to those computed by PHREG;
proc phreg data = dataclinical_gdr_20210624; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
 where subgroup_adl=1;
run;

data covalladl (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_adl=1;
run;
/* 6001 observations and 21 variables*/

*So, since I only interested in the CumInc not the 95%CI I will omit the covsandwich(aggregate) suboption;
proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  baseline out=clinical.CumIncPredADL10(keep=newid time_adldepdth2 cif_adl) covariates=covalladl CIF=cif_adl timelist=10;
run;

data clinical.CumIncPredADL10;
 set clinical.CumIncPredADL10;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredADL10
            OUTFILE= "path\Rfiles\CumIncPredADL10.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


**********************************************************;
/*15-variable best model missing 1 variable (R5EATA): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g stroke3g*/

data covalladl (keep=newid age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_adl=1;
run;
/* 6001 observations and 20 variables*/


proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  baseline out=clinical.CumIncPredADL10_1UnavailBest(keep=newid time_adldepdth2 cif_adl) covariates=covalladl CIF=cif_adl timelist=10;
run;

data clinical.CumIncPredADL10_1UnavailBest;
 set clinical.CumIncPredADL10_1UnavailBest;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredADL10_1UnavailBest
            OUTFILE= "path\Rfiles\CumIncPredADL10_1UnavailBest.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalladl; run; quit;

**********************************************************;
/*15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi4g smoking3g stroke3g*/

data covalladl (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_adl=1;
run;
/* 6001 observations and 20 variables*/

proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  baseline out=clinical.CumIncPredADL10_1UnavailWorst(keep=newid time_adldepdth2 cif_adl) covariates=covalladl CIF=cif_adl timelist=10;
run;

data clinical.CumIncPredADL10_1UnavailWorst;
 set clinical.CumIncPredADL10_1UnavailWorst;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredADL10_1UnavailWorst
            OUTFILE= "path\Rfiles\CumIncPredADL10_1UnavailWorst.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalladl; run; quit;


**********************************************************;
/*14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g*/

data covalladl (keep=newid age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g);
 set dataclinical_gdr_20210624;
 where subgroup_adl=1;
run;
/* 6001 observations and 19 variables*/


proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  baseline out=clinical.CumIncPredADL10_2UnavailBest(keep=newid time_adldepdth2 cif_adl) covariates=covalladl CIF=cif_adl timelist=10;
run;

data clinical.CumIncPredADL10_2UnavailBest;
 set clinical.CumIncPredADL10_2UnavailBest;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredADL10_2UnavailBest
            OUTFILE= "path\Rfiles\CumIncPredADL10_2UnavailBest.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalladl; run; quit;


**********************************************************;
/*14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi4g smoking3g stroke3g*/

data covalladl (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_adl=1;
run;
/* 6001 observations and 19 variables*/

proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  baseline out=clinical.CumIncPredADL10_2UnavailWorst(keep=newid time_adldepdth2 cif_adl) covariates=covalladl CIF=cif_adl timelist=10;
run;

data clinical.CumIncPredADL10_2UnavailWorst;
 set clinical.CumIncPredADL10_2UnavailWorst;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredADL10_2UnavailWorst
            OUTFILE= "path\Rfiles\CumIncPredADL10_2UnavailWorst.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalladl; run; quit;


/******************************************************************************************************************************************************************/
*Walk: Compute weighted predicted cumulative incidence for Calibration plots and Tables;

*Confirmed that unweighted estimnates and SEs in R are very similar to those computed by PHREG;
proc phreg data = dataclinical_gdr_20210624; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
 where subgroup_walk=1;
run;

data covallwalk (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_walk=1;
run;
/* 6409 observations and 21 variables*/

*So, since I only interested in the CumInc not the 95%CI I will omit the covsandwich(aggregate) suboption;
proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  baseline out=clinical.CumIncPredwalk10(keep=newid time_walkdepdth2 cif_walk) covariates=covallwalk CIF=cif_walk timelist=10;
run;

data clinical.CumIncPredwalk10;
 set clinical.CumIncPredwalk10;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredwalk10
            OUTFILE= "path\Rfiles\CumIncPredwalk10.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


**********************************************************;
/*15-variable best model missing 1 variable (R5EATA): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g stroke3g*/

data covallwalk (keep=newid age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_walk=1;
run;
/* 6409 observations and 20 variables*/


proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  baseline out=clinical.CumIncPredWALK10_1UnavailBest(keep=newid time_walkdepdth2 cif_walk) covariates=covallwalk CIF=cif_walk timelist=10;
run;

data clinical.CumIncPredWALK10_1UnavailBest;
 set clinical.CumIncPredWALK10_1UnavailBest;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredWALK10_1UnavailBest
            OUTFILE= "path\Rfiles\CumIncPredWALK10_1UnavailBest.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covallwalk; run; quit;

**********************************************************;
/*15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi4g smoking3g stroke3g*/

data covallwalk (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_walk=1;
run;
/* 6409 observations and 20 variables*/

proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  baseline out=clinical.CumIncPredWALK10_1UnavailWorst(keep=newid time_walkdepdth2 cif_walk) covariates=covallwalk CIF=cif_walk timelist=10;
run;

data clinical.CumIncPredWALK10_1UnavailWorst;
 set clinical.CumIncPredWALK10_1UnavailWorst;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredWALK10_1UnavailWorst
            OUTFILE= "path\Rfiles\CumIncPredWALK10_1UnavailWorst.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covallwalk; run; quit;


**********************************************************;
/*14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g*/

data covallwalk (keep=newid age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g);
 set dataclinical_gdr_20210624;
 where subgroup_walk=1;
run;
/* 6409 observations and 19 variables*/


proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  baseline out=clinical.CumIncPredWALK10_2UnavailBest(keep=newid time_walkdepdth2 cif_walk) covariates=covallwalk CIF=cif_walk timelist=10;
run;

data clinical.CumIncPredWALK10_2UnavailBest;
 set clinical.CumIncPredWALK10_2UnavailBest;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredWALK10_2UnavailBest
            OUTFILE= "path\Rfiles\CumIncPredWALK10_2UnavailBest.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covallwalk; run; quit;

**********************************************************;
/*14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi4g smoking3g stroke3g*/

data covallwalk (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
 where subgroup_walk=1;
run;
/* 6409 observations and 19 variables*/

proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  baseline out=clinical.CumIncPredWALK10_2UnavailWorst(keep=newid time_walkdepdth2 cif_walk) covariates=covallwalk CIF=cif_walk timelist=10;
run;

data clinical.CumIncPredWALK10_2UnavailWorst;
 set clinical.CumIncPredWALK10_2UnavailWorst;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredWALK10_2UnavailWorst
            OUTFILE= "path\Rfiles\CumIncPredWALK10_2UnavailWorst.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covallwalk; run; quit;


/******************************************************************************************************************************************************************/
*Death: Compute weighted predicted survival probability for Calibration plots;

*Confirmed that unweighted estimnates and SEs in R are very similar to those computed by PHREG;
proc phreg data = dataclinical_gdr_20210624; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g / risklimits;
run;

data covalldeath (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
run;
/* 6646 observations and 21 variables*/

*So, since I only interested in the Survival probability not the 95%CI I will omit the covsandwich(aggregate) suboption;
proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /risklimits;
 baseline out=clinical.SurvPredDth10(keep=newid time2death s_death) covariates=covalldeath survival=s_death timelist=10;
run;

data clinical.SurvPredDth10;
 set clinical.SurvPredDth10;
 cif_death=1-s_death;
proc sort; by newid; run;


*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.SurvPredDth10
            OUTFILE= "path\Rfiles\SurvPredDth10.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


**********************************************************;
/*15-variable best model missing 1 variable (R5EATA): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g stroke3g*/

data covalldeath (keep=newid age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
run;
/* 6646 observations and 20 variables*/


proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g / risklimits;
  baseline out=clinical.SurvPredDth10_1UnavailBest(keep=newid time2death s_death) covariates=covalldeath survival=s_death timelist=10;
run;

data clinical.SurvPredDth10_1UnavailBest;
 set clinical.SurvPredDth10_1UnavailBest;
 cif_death=1-s_death;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.SurvPredDth10_1UnavailBest
            OUTFILE= "path\Rfiles\SurvPredDth10_1UnavailBest.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalldeath; run; quit;

**********************************************************;
/*15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi4g smoking3g stroke3g*/

data covalldeath (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
run;
/* 6646 observations and 20 variables*/

proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g / risklimits;
  baseline out=clinical.SurvPredDth10_1UnavailWorst(keep=newid time2death s_death) covariates=covalldeath survival=s_death timelist=10;
run;

data clinical.SurvPredDth10_1UnavailWorst;
 set clinical.SurvPredDth10_1UnavailWorst;
 cif_death=1-s_death;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.SurvPredDth10_1UnavailWorst
            OUTFILE= "path\Rfiles\SurvPredDth10_1UnavailWorst.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalldeath; run; quit;


**********************************************************;
/*14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g*/

data covalldeath (keep=newid age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g);
 set dataclinical_gdr_20210624;
run;
/* 6646 observations and 19 variables*/


proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g / risklimits;
  baseline out=clinical.SurvPredDth10_2UnavailBest(keep=newid time2death s_death) covariates=covalldeath survival=s_death timelist=10;
run;

data clinical.SurvPredDth10_2UnavailBest;
 set clinical.SurvPredDth10_2UnavailBest;
 cif_death=1-s_death;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.SurvPredDth10_2UnavailBest
            OUTFILE= "path\Rfiles\SurvPredDth10_2UnavailBest.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalldeath; run; quit;

**********************************************************;
/*14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi4g smoking3g stroke3g*/

data covalldeath (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g);
 set dataclinical_gdr_20210624;
run;
/* 6646 observations and 19 variables*/

proc phreg data = dataclinical_gdr_20210624 ; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g / risklimits;
  baseline out=clinical.SurvPredDth10_2UnavailWorst(keep=newid time2death s_death) covariates=covalldeath survival=s_death timelist=10;
run;

data clinical.SurvPredDth10_2UnavailWorst;
 set clinical.SurvPredDth10_2UnavailWorst;
 cif_death=1-s_death;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.SurvPredDth10_2UnavailWorst
            OUTFILE= "path\Rfiles\SurvPredDth10_2UnavailWorst.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

proc delete data=covalldeath; run; quit;


