/*********** Purpose: Compute "average for age" risk and compute Cumulative Incidence for calibration plots ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2022.04.20 ***********/
/*********** Date completed: 2022.05.20 ***********/

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
libname outdata2 "path";


/******************************************************************************************************************************************************************/
*Prepare data ;

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

*Create dummy variables for categorical variables and create categorical variable for age:70-74, 75-79, 80-84, 85+ ;
data dataclinical_gdr_20210624;
 set dataclinical_gdr_20210624;

 if heart3g=0 then do; heart1=0; heart2=0; end;
 else if heart3g=1 then do; heart1=1; heart2=0; end;
 else if heart3g=2 then do; heart1=0; heart2=1; end;

 if smoking3g=0 then do; smoking1=0; smoking2=0; end;
 else if smoking3g=1 then do; smoking1=1; smoking2=0; end;
 else if smoking3g=2 then do; smoking1=0; smoking2=1; end;

 if stroke3g=0 then do; stroke1=0; stroke2=0; end;
 else if stroke3g=1 then do; stroke1=1; stroke2=0; end;
 else if stroke3g=2 then do; stroke1=0; stroke2=1; end;

 if 70<=ageint5<=74 then age4g="70_74";
 else if 75<=ageint5<=79 then age4g="75_79";
 else if 80<=ageint5<=84 then age4g="80_84";
 else if ageint5>=85 then age4g="85+";

 label age4g="R age group at baseline.1.70-74, 2.75-79, 3.80-84, 4.85+";
run;
/* 6646 observations and 72+7=79 variables */

/*QC*/
proc freq data=dataclinical_gdr_20210624; tables age4g; run;
proc means data=dataclinical_gdr_20210624 n min max; var ageint5; class age4g; run;


/******************************************************************************************************************************************************************/

/*
 Fit a model with age as only predictor and add the betas and the baseline survival functions of this model to the files
 where we have the betas and baseline survival function for the final model and the unavailable predictor models.

 Helen would use these betas and baseline survival functions to compute the risk at 5, 10, and 14 years for ages 70-100 the same way she would do it for the final model and unavailable predictors
 models. This would be the "average for age" prediction.
 */

*Import data with age 70-100 and splines;
PROC IMPORT OUT= splinesAge 
            DATAFILE= "path\Rfiles\splinesAge_gdr_20220413.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;


*Covariate data with all Ids;
data covall (keep=newid age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart1 heart2 R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA
                       R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking1 smoking2 stroke1 stroke2);
 set dataclinical_gdr_20210624;
proc sort; by newid; run; /* 6646 observations and 24 variables*/

*Covariate data with all covariates=0 and only 1 Id;
data covall_base1obs;
 set covall;
 where newid=16;
 array base[*] age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart1 heart2 R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA
                       R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking1 smoking2 stroke1 stroke2;
 do i=1 to dim(base);
  base[i]=0;
 end; drop i;
run;
/* 1 observations and 24 variables */

proc means data=covall_base1obs n min max; run;

*Modify data with Wolbers approximation;
proc freq data=dataclinical_gdr_20210624; tables status_adldepdth2; run;
data dataclinical_gdr_20210624;
 set dataclinical_gdr_20210624;
 if status_adldepdth2=2 then do; status_adldepdth2=0; time_adldepdth2=19.1130249; end;
run;
proc freq data=dataclinical_gdr_20210624; tables status_adldepdth2; run;

*Modify data with Wolbers approximation;
proc freq data=dataclinical_gdr_20210624; tables status_walkdepdth2; run;
data dataclinical_gdr_20210624;
 set dataclinical_gdr_20210624;
 if status_walkdepdth2=2 then do; status_walkdepdth2=0; time_walkdepdth2=19.1130249; end;
run;
proc freq data=dataclinical_gdr_20210624; tables status_walkdepdth2; run;


*Get the weighted betas and baseline survival for 3 outcomes and model with age as only predictor;

%let NUMOUTCOMES=3; /*number of outcomes*/
%let ALLOUTCOME=status_adldepdth2 status_walkdepdth2 death;
%let ALLTIME=time_adldepdth2 time_walkdepdth2 time2death;
%let ALLLABEL= ADL WALK DEATH; /*labels for outcomes. I put them in capital letters so that it's recognized in proc sql with memname option*/
%let VARNAME=age3_lin age3_sp1 age3_sp2;

options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors/warnings to be written to the log*/
options nomlogic nomprint nomrecall nosymbolgen;
ods select none; /*to create output data sets through the ODS OUTPUT statement and suppress the display of all output*/

%macro betas_basesurv;

  %do j=1 %to &NUMOUTCOMES;
    %let OUTCOME=%scan(&ALLOUTCOME,&j); /*extract the jth outcome, jth time, jth label*/
    %let TIME=%scan(&ALLTIME,&j);
	%let LABEL=%scan(&ALLLABEL,&j);

	*1) Within SAS: Compute weighted betas (dataset called “WeightBetaSurvAvgAge_&label”):  ;
	proc phreg data = dataclinical_gdr_20210624 outest=WeightBetaSurvAvgAge_&label; 
	  weight R5WTRESP; /*use survey weight*/
	  model &time*&outcome(0) = &varname /risklimits;
	  where subgroup_&label=1;
	run;

	*2) Within SAS: Obtain So(t) (i.e. So_&label using above weighted betas and METHOD=PL;
	proc phreg data = dataclinical_gdr_20210624 inest=WeightBetaSurvAvgAge_&label ; 
	  model &time*&outcome(0) = &varname /risklimits maxit=0;
	 baseline out=WeightBaseSurvAvgAge_&label (keep=&time So_&label) covariates=covall_base1obs survival=So_&label timelist=1 to 19 by 1 / method=PL;
	run;
	/*WeightBaseSurvAvgAge_&label has 19 observations and 2 variables.*/

	*Keep variables of interest in WeightBetaSurvAvgAge_&label;
	data WeightBetaSurvAvgAge_&label;
	 set WeightBetaSurvAvgAge_&label (keep=&varname);
	run;

	*Rename variables in WeightBetaSurvAvgAge_&label;
	/* This code creates a macro variable &list with the list of variables in the form. */
	/* variable = variable_beta                                                         */
	/* This format could be used to add a suffix to all the variables.                  */
	proc sql noprint;
	   select cats(name,'=',name,'_beta')
	          into :list
	          separated by ' '
	          from dictionary.columns
	          where libname = 'WORK' and memname = "WEIGHTBETASURVAVGAGE_&label";
	quit;

	proc datasets library = work nolist;
	   modify WEIGHTBETASURVAVGAGE_&label;
	   rename &list;
	quit;

	/*Add DELEVAR	VARINMODEL	nummissvar	modelnum */
	data WeightBetaSurvAvgAge_&label;
	 set WeightBetaSurvAvgAge_&label;
	 DELEVAR="all except age";
     VARINMODEL="&varname";
	 nummissvar=15;
     modelnum=0;
	run;

	*Traspose Base survival data from long to wide format;
	data WeightBaseSurvAvgAge_&label  (keep=So_&label._t1-So_&label._t19 DELEVAR VARINMODEL nummissvar modelnum);
	 set WeightBaseSurvAvgAge_&label  end=last;
	 array surv[19] So_&label._t1-So_&label._t19;
	 retain So_&label._t1-So_&label._t19;
	 do i=1 to 19;
	  if &time=i then surv[i]=So_&label;
	 end; drop i;
	 if last;
     DELEVAR="all except age";
     VARINMODEL="&varname";
	 nummissvar=15;
     modelnum=0;
   run;


  %end; /*j loop*/



%mend betas_basesurv;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%betas_basesurv
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;



/******************************************************************************************************************************************************************/
*Compute manually the predicted risk at 5, 10, 14 years;

*ADL;
data adl_risk (keep=age3_lin risk_t5_adl risk_t10_adl risk_t14_adl);
 if ( _n_ eq 1 ) then set weightbasesurvavgage_adl;
 if ( _n_ eq 1 ) then set Weightbetasurvavgage_adl ;
 set splinesage ;
 XBETA=age3_lin*age3_lin_beta+age3_sp1*age3_sp1_beta+age3_sp2*age3_sp2_beta;
 risk_t5_adl=1-(So_ADL_t5**exp(XBETA));
 risk_t10_adl=1-(So_ADL_t10**exp(XBETA));
 risk_t14_adl=1-(So_ADL_t14**exp(XBETA));
run;
/* 31 observations and 4 variables*/

*Walk;
data walk_risk (keep=age3_lin risk_t5_walk risk_t10_walk risk_t14_walk);
 if ( _n_ eq 1 ) then set weightbasesurvavgage_walk;
 if ( _n_ eq 1 ) then set Weightbetasurvavgage_walk ;
 set splinesage ;
 XBETA=age3_lin*age3_lin_beta+age3_sp1*age3_sp1_beta+age3_sp2*age3_sp2_beta;
 risk_t5_walk=1-(So_WALK_t5**exp(XBETA));
 risk_t10_walk=1-(So_WALK_t10**exp(XBETA));
 risk_t14_walk=1-(So_WALK_t14**exp(XBETA));
run;
/* 31 observations and 4 variables*/

*Death;
data death_risk (keep=age3_lin risk_t5_death risk_t10_death risk_t14_death);
 if ( _n_ eq 1 ) then set weightbasesurvavgage_death;
 if ( _n_ eq 1 ) then set Weightbetasurvavgage_death ;
 set splinesage ;
 XBETA=age3_lin*age3_lin_beta+age3_sp1*age3_sp1_beta+age3_sp2*age3_sp2_beta;
 risk_t5_death=1-(So_DEATH_t5**exp(XBETA));
 risk_t10_death=1-(So_DEATH_t10**exp(XBETA));
 risk_t14_death=1-(So_DEATH_t14**exp(XBETA));
run;
/* 31 observations and 4 variables*/

*Merge 3 risk data;
data clinical.avg_risks_by_age70_100;
 merge adl_risk walk_risk death_risk;
 by age3_lin;
run;


*Export as csv file;
PROC EXPORT DATA= clinical.avg_risks_by_age70_100
            OUTFILE= "path\avg_risks_by_age70_100.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;


 *Append weighted betas and baseline survival function from models with age only to the weighted betas and baseline survival functions from unavailable predictors;

%macro finaldata;

%do j=1 %to &NUMOUTCOMES;

 %let OUTCOME=%scan(&ALLOUTCOME,&j); /*extract the jth outcome, jth time, jth label*/
 %let TIME=%scan(&ALLTIME,&j);
 %let LABEL=%scan(&ALLLABEL,&j);

 data clinical.Weightbetasurv_&label;
  set clinical.Weightbetasurvall_&label WeightBetaSurvAvgAge_&label;
 proc sort; by modelnum; run;

 data clinical.Weightbasesur_&label;
  set clinical.Weightbasesurvall_&label WeightBaseSurvAvgAge_&label ;
 proc sort; by modelnum; run;

 PROC EXPORT DATA= clinical.Weightbetasurv_&label
            OUTFILE= "path\Rfiles\Weightbetasurv_&label..csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;
 PROC EXPORT DATA= clinical.Weightbasesur_&label
            OUTFILE= "path\Rfiles\Weightbasesurv_&label..csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;

 %end; /*j loop*/

%mend finaldata;

%finaldata

;


/******************************************************************************************************************************************************************/
*Compute weighted predicted cumulative incidence for Calibration plots for model wih age only;

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


/*ADL*/
data covalladl (keep=newid age3_lin age3_sp1 age3_sp2);
 set dataclinical_gdr_20210624;
 where subgroup_adl=1;
run;
/* 6001 observations and 4 variables*/

*So, since I only interested in the CumInc not the 95%CI I will omit the covsandwich(aggregate) suboption;
proc phreg data = dataclinical_gdr_20210624 ; 
/*  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')*/
/*                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');*/
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 /eventcode=1 risklimits;
  where subgroup_adl=1;
  baseline out=clinical.CumIncPredADL10_ageonly(keep=newid time_adldepdth2 cif_adl) covariates=covalladl CIF=cif_adl timelist=10;
run;

data clinical.CumIncPredADL10_ageonly;
 set clinical.CumIncPredADL10_ageonly;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredADL10_ageonly
            OUTFILE= "path\Rfiles\CumIncPredADL10_ageonly.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


/*WALK*/
data covallwalk (keep=newid age3_lin age3_sp1 age3_sp2);
 set dataclinical_gdr_20210624;
 where subgroup_walk=1;
run;
/* 6409 observations and 4 variables*/

*So, since I only interested in the CumInc not the 95%CI I will omit the covsandwich(aggregate) suboption;
proc phreg data = dataclinical_gdr_20210624 ; 
/*  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')*/
/*                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');*/
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 /eventcode=1 risklimits;
  where subgroup_walk=1;
  baseline out=clinical.CumIncPredWALK10_ageonly(keep=newid time_walkdepdth2 cif_walk) covariates=covallwalk CIF=cif_walk timelist=10;
run;

data clinical.CumIncPredWALK10_ageonly;
 set clinical.CumIncPredWALK10_ageonly;
proc sort; by newid; run;

*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.CumIncPredWALK10_ageonly
            OUTFILE= "path\Rfiles\CumIncPredWALK10_ageonly.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


/*Death*/
data covalldeath (keep=newid age3_lin age3_sp1 age3_sp2);
 set dataclinical_gdr_20210624;
run;
/* 6446 observations and 4 variables*/

*So, since I only interested in the CumInc not the 95%CI I will omit the covsandwich(aggregate) suboption;
proc phreg data = dataclinical_gdr_20210624 ; 
/*  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')*/
/*                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');*/
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 / risklimits;
  baseline out=clinical.SurvPredDth10_ageonly(keep=newid time2death s_death) covariates=covalldeath survival=s_death timelist=10;
run;

data clinical.SurvPredDth10_ageonly;
 set clinical.SurvPredDth10_ageonly;
 cif_death=1-s_death;
proc sort; by newid; run;


*Export to R to do calibration plot;
PROC EXPORT DATA= clinical.SurvPredDth10_ageonly
            OUTFILE= "path\Rfiles\SurvPredDth10_ageonly.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


/******************************************************************************************************************************************************************/
*Calculate the predicted risk for each individual [using the full model, with their original covariate values] for 10-yr risk of each outcome
Then, categorize each participant as "lower than average", "about average" or "higher than average" using the +/-5% margin,
and calculate the percent that fall in each group, by age
;

*Get data with Predicted CIF for each outcome;
data cumincpredadl10;
 set clinical.cumincpredadl10;
proc sort; by newid; run;
data cumincpredwalk10;
 set clinical.cumincpredwalk10;
proc sort; by newid; run;
data survpreddth10;
 set clinical.survpreddth10;
proc sort; by newid; run;

*Merge 3 datasets;
data cuminc(drop=time_adldepdth2 time_walkdepdth2 time2death s_death);
 merge cumincpredadl10 cumincpredwalk10 survpreddth10;
 by newid;
proc sort; by newid; run;
/*6646 obs and 2+1+1=4 variables*/

*Get the age and weight of each Id;
data dataclinical_gdr_20210624; 
 set clinical.clinicalHRSdeident_gdr_20210624 (keep=newid age3_lin R5WTRESP);
run;
/* 6646 observations and 3 variables */

data cuminc;
 merge cuminc dataclinical_gdr_20210624;
 by newid;
proc sort; by age3_lin; run;
/*6646 obs and 4+2=6 variables*/

proc means data=cuminc n nmiss min max; var age3_lin; run; /*range: 70-102*/

*Get average by age prediction;
data avg_risks_by_age70_100;
 set clinical.avg_risks_by_age70_100 (keep=age3_lin risk_t10_adl risk_t10_walk risk_t10_death) ;
run;
/* 31 observations and 4 variables. */

*Merge avg_risks_by_age70_100 and cuminc;
data cuminc2;
 merge avg_risks_by_age70_100 cuminc;
 by age3_lin;
run;
/*6646 obs and 6+3=9 variables*/

proc contents data=cuminc2; run;
proc means data=cuminc2 n nmiss min max; run;


data cuminc3;
 set cuminc2;
 length risk_group5_adl risk_group5_walk risk_group5_death risk_group3_adl risk_group3_walk risk_group3_death $10;

 if cif_adl ne . and risk_t10_adl ne . then do;
  if cif_adl<risk_t10_adl-0.05 then risk_group5_adl="1.lower";
  else if risk_t10_adl-0.05<=cif_adl<=risk_t10_adl+0.05 then risk_group5_adl="2.similar";
  else if cif_adl>risk_t10_adl+0.05 then risk_group5_adl="3.higher";
 end;

 if cif_walk ne . and risk_t10_walk ne . then do;
  if cif_walk<risk_t10_walk-0.05 then risk_group5_walk="1.lower";
  else if risk_t10_walk-0.05<=cif_walk<=risk_t10_walk+0.05 then risk_group5_walk="2.similar";
  else if cif_walk>risk_t10_walk+0.05 then risk_group5_walk="3.higher";
 end;

 if cif_death ne . and risk_t10_death ne . then do;
  if cif_death<risk_t10_death-0.05 then risk_group5_death="1.lower";
  else if risk_t10_death-0.05<=cif_death<=risk_t10_death+0.05 then risk_group5_death="2.similar";
  else if cif_death>risk_t10_death+0.05 then risk_group5_death="3.higher";
 end;

 if cif_adl ne . and risk_t10_adl ne . then do;
  if cif_adl<risk_t10_adl-0.03 then risk_group3_adl="1.lower";
  else if risk_t10_adl-0.03<=cif_adl<=risk_t10_adl+0.03 then risk_group3_adl="2.similar";
  else if cif_adl>risk_t10_adl+0.03 then risk_group3_adl="3.higher";
 end;

 if cif_walk ne . and risk_t10_walk ne . then do;
  if cif_walk<risk_t10_walk-0.03 then risk_group3_walk="1.lower";
  else if risk_t10_walk-0.03<=cif_walk<=risk_t10_walk+0.03 then risk_group3_walk="2.similar";
  else if cif_walk>risk_t10_walk+0.03 then risk_group3_walk="3.higher";
 end;

 if cif_death ne . and risk_t10_death ne . then do;
  if cif_death<risk_t10_death-0.03 then risk_group3_death="1.lower";
  else if risk_t10_death-0.03<=cif_death<=risk_t10_death+0.03 then risk_group3_death="2.similar";
  else if cif_death>risk_t10_death+0.03 then risk_group3_death="3.higher";
 end;

run;
/*6646 obs and 9+6=15 variables*/

/*QC*/
proc means data=cuminc3 n nmiss min max; where age3_lin<=100; run;
proc freq data=cuminc3; tables risk_group5_adl risk_group5_walk risk_group5_death risk_group3_adl risk_group3_walk risk_group3_death; where age3_lin<=100; run;

proc means data=cuminc3 n nmiss min max; var cif_adl ; class age3_lin risk_group5_adl ; where age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var cif_walk ; class age3_lin risk_group5_walk ; where age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var cif_death ; class age3_lin risk_group5_death ; where age3_lin<=100; run;

proc means data=cuminc3 n nmiss min max; var cif_adl ; class age3_lin risk_group3_adl ; where age3_lin<=100; run;

/*Compute the 2x2 tables*/

*Death 5%;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group5_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_death ne . and risk_t10_death ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group5_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_death ne . and risk_t10_death ne .; 
run;


*ADL 5%;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group5_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_adl ne . and risk_t10_adl ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group5_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_adl ne . and risk_t10_adl ne .; 
run;


*Walk 5%;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group5_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_walk ne . and risk_t10_walk ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group5_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_walk ne . and risk_t10_walk ne .; 
run;


*Death 3%;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group3_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_death ne . and risk_t10_death ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group3_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_death ne . and risk_t10_death ne .; 
run;


*ADL 3%;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group3_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_adl ne . and risk_t10_adl ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group3_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_adl ne . and risk_t10_adl ne .; 
run;


*Walk 3%;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group3_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_walk ne . and risk_t10_walk ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group3_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_walk ne . and risk_t10_walk ne .; 
run;


*Combine some ages with low N in tables above;

*90 category;
proc means data=cuminc3 n nmiss min max; var cif_death ; class risk_group5_death; where age3_lin=90; run;
proc means data=cuminc3 n nmiss min max; var cif_adl  ; class risk_group5_adl; where age3_lin=90; run;
proc means data=cuminc3 n nmiss min max; var  cif_walk ; class risk_group5_walk; where age3_lin=90; run;

proc means data=cuminc3 n nmiss min max; var cif_death ; class risk_group5_death; where 89<=age3_lin<=91; run;
proc means data=cuminc3 n nmiss min max; var cif_adl  ; class risk_group5_adl; where 89<=age3_lin<=91; run;
proc means data=cuminc3 n nmiss min max; var  cif_walk ; class risk_group5_walk; where 89<=age3_lin<=91;  run;


proc means data=cuminc3 n nmiss min max; var cif_death ; class risk_group3_death; where age3_lin=90; run;
proc means data=cuminc3 n nmiss min max; var cif_adl  ; class risk_group3_adl; where age3_lin=90; run;
proc means data=cuminc3 n nmiss min max; var  cif_walk ; class risk_group3_walk; where age3_lin=90; run;

proc means data=cuminc3 n nmiss min max; var cif_death ; class risk_group3_death; where 89<=age3_lin<=91; run;
proc means data=cuminc3 n nmiss min max; var cif_adl  ; class risk_group3_adl; where 89<=age3_lin<=91; run;
proc means data=cuminc3 n nmiss min max; var  cif_walk ; class risk_group3_walk; where 89<=age3_lin<=91;  run;


*95 average category;
proc means data=cuminc3 n nmiss min max; var cif_death ; class risk_group5_death; where 92<=age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var cif_adl  ; class risk_group5_adl; where 92<=age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var  cif_walk ; class risk_group5_walk; where 92<=age3_lin<=100;  run;

proc means data=cuminc3 n nmiss min max; var cif_death ; class risk_group3_death; where 92<=age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var cif_adl  ; class risk_group3_adl; where 92<=age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var  cif_walk ; class risk_group3_walk; where 92<=age3_lin<=100;  run;


*Change data:
if 89<=age3_lin<=91 then age3_lin=90
else if 92<=age3_lin<=100 then age3_lin=96;


*Merge avg_risks_by_age70_100 and cuminc;
data cumincb;
 set cuminc;
 if 84<=age3_lin<=86 then age3_lin=85;
 else if 89<=age3_lin<=91 then age3_lin=90;
 else if 92<=age3_lin<=100 then age3_lin=96;
proc sort; by age3_lin; run;
/*6646 obs and 4+2=6 variables*/

data cuminc2b;
 merge avg_risks_by_age70_100 cumincb(in=A);
 by age3_lin;
 if A;
run;
/*6646 obs and 6+3=9 variables*/

proc contents data=cuminc2b; run;
proc means data=cuminc2b n nmiss min max; run;


data cuminc3b;
 set cuminc2b;
 length risk_group5_adl risk_group5_walk risk_group5_death risk_group3_adl risk_group3_walk risk_group3_death $10;

 if cif_adl ne . and risk_t10_adl ne . then do;
  if cif_adl<risk_t10_adl-0.05 then risk_group5_adl="1.lower";
  else if risk_t10_adl-0.05<=cif_adl<=risk_t10_adl+0.05 then risk_group5_adl="2.similar";
  else if cif_adl>risk_t10_adl+0.05 then risk_group5_adl="3.higher";
 end;

 if cif_walk ne . and risk_t10_walk ne . then do;
  if cif_walk<risk_t10_walk-0.05 then risk_group5_walk="1.lower";
  else if risk_t10_walk-0.05<=cif_walk<=risk_t10_walk+0.05 then risk_group5_walk="2.similar";
  else if cif_walk>risk_t10_walk+0.05 then risk_group5_walk="3.higher";
 end;

 if cif_death ne . and risk_t10_death ne . then do;
  if cif_death<risk_t10_death-0.05 then risk_group5_death="1.lower";
  else if risk_t10_death-0.05<=cif_death<=risk_t10_death+0.05 then risk_group5_death="2.similar";
  else if cif_death>risk_t10_death+0.05 then risk_group5_death="3.higher";
 end;

 if cif_adl ne . and risk_t10_adl ne . then do;
  if cif_adl<risk_t10_adl-0.03 then risk_group3_adl="1.lower";
  else if risk_t10_adl-0.03<=cif_adl<=risk_t10_adl+0.03 then risk_group3_adl="2.similar";
  else if cif_adl>risk_t10_adl+0.03 then risk_group3_adl="3.higher";
 end;

 if cif_walk ne . and risk_t10_walk ne . then do;
  if cif_walk<risk_t10_walk-0.03 then risk_group3_walk="1.lower";
  else if risk_t10_walk-0.03<=cif_walk<=risk_t10_walk+0.03 then risk_group3_walk="2.similar";
  else if cif_walk>risk_t10_walk+0.03 then risk_group3_walk="3.higher";
 end;

 if cif_death ne . and risk_t10_death ne . then do;
  if cif_death<risk_t10_death-0.03 then risk_group3_death="1.lower";
  else if risk_t10_death-0.03<=cif_death<=risk_t10_death+0.03 then risk_group3_death="2.similar";
  else if cif_death>risk_t10_death+0.03 then risk_group3_death="3.higher";
 end;

run;
/*6646 obs and 9+6=15 variables*/

/*QC*/
proc means data=cuminc3b n nmiss min max; where age3_lin<=100; run;
proc freq data=cuminc3b; tables risk_group5_adl risk_group5_walk risk_group5_death risk_group3_adl risk_group3_walk risk_group3_death; where age3_lin<=100; run;

proc means data=cuminc3b n nmiss min max; var cif_adl ; class age3_lin risk_group5_adl ; where age3_lin<=100; run;
proc means data=cuminc3b n nmiss min max; var cif_walk ; class age3_lin risk_group5_walk ; where age3_lin<=100; run;
proc means data=cuminc3b n nmiss min max; var cif_death ; class age3_lin risk_group5_death ; where age3_lin<=100; run;

proc means data=cuminc3b n nmiss min max; var cif_adl ; class age3_lin risk_group3_adl ; where age3_lin<=100; run;

/*Compute the 2x2 tables*/

*Death 5%;
proc freq data=cuminc3b;
 weight R5WTRESP;
 tables age3_lin*risk_group5_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_death ne . and risk_t10_death ne .; 
run;

proc freq data=cuminc3b;
 tables age3_lin*risk_group5_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_death ne . and risk_t10_death ne .; 
run;


*ADL 5%;
proc freq data=cuminc3b;
 weight R5WTRESP;
 tables age3_lin*risk_group5_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_adl ne . and risk_t10_adl ne .; 
run;

proc freq data=cuminc3b;
 tables age3_lin*risk_group5_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_adl ne . and risk_t10_adl ne .; 
run;


*Walk 5%;
proc freq data=cuminc3b;
 weight R5WTRESP;
 tables age3_lin*risk_group5_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_walk ne . and risk_t10_walk ne .; 
run;

proc freq data=cuminc3b;
 tables age3_lin*risk_group5_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_walk ne . and risk_t10_walk ne .; 
run;


*Death 3%;
proc freq data=cuminc3b;
 weight R5WTRESP;
 tables age3_lin*risk_group3_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_death ne . and risk_t10_death ne .; 
run;

proc freq data=cuminc3b;
 tables age3_lin*risk_group3_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_death ne . and risk_t10_death ne .; 
run;


*ADL 3%;
proc freq data=cuminc3b;
 weight R5WTRESP;
 tables age3_lin*risk_group3_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_adl ne . and risk_t10_adl ne .; 
run;

proc freq data=cuminc3b;
 tables age3_lin*risk_group3_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_adl ne . and risk_t10_adl ne .; 
run;


*Walk 3%;
proc freq data=cuminc3b;
 weight R5WTRESP;
 tables age3_lin*risk_group3_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_walk ne . and risk_t10_walk ne .; 
run;

proc freq data=cuminc3b;
 tables age3_lin*risk_group3_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 96) and cif_walk ne . and risk_t10_walk ne .; 
run;


/******************************************************************************************************************************************************************/
*Categorize each participant as "lower than average", "about average" or "higher than average" using +/- 2SD ;

/*
Steps:
1) Compute predicted risk/outcome for all Ids in the cohort using age-only model at 10 years
2) Compute weighted SD of predicted risk/outcome
3) Compute predicted risk/outcome for all Ids in the cohort using the full model, with their original covariate values at 10 years
4) Compute z-score=x-mu/SD
where x: value in 3
      mu: average risk for age (from age-only model)
      SD: computed in 2 
5) Classify Ids as:
 lower: z-score<-2SD
 similar: -2SD<=z-score<=2SD
 higher: z-score>2SD
*/

data clinicalhrs_gdr_20210624; 
 set clinical.clinicalhrs_gdr_20210624 (keep=hhidpn RAEHSAMP RAESTRAT clustervar R5WTRESP ageint5 age_lin age_sp1 age_sp2 age_sp3 age3_lin age3_sp1 age3_sp2
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

*Covariate data with all Ids;

*ADL;
data adl (keep=hhidpn age3_lin age3_sp1 age3_sp2 R5WTRESP RAEHSAMP RAESTRAT );
 set clinicalhrs_gdr_20210624;
 where subgroup_adl=1;
proc sort; by hhidpn; run; /* 6001 observations and 7 variables*/

*Walk;
data walk (keep=hhidpn age3_lin age3_sp1 age3_sp2 );
 set clinicalhrs_gdr_20210624;
 where subgroup_walk=1;
proc sort; by hhidpn; run; /* 6409 observations and 4 variables*/

*Death;
data death (keep=hhidpn age3_lin age3_sp1 age3_sp2  );
 set clinicalhrs_gdr_20210624;
proc sort; by hhidpn; run; /* 6646 observations and 4 variables*/

***********************************************************************************************;
*1) Compute predicted risk/outcome for all Ids in the cohort using age-only model at 10 years;

*Run macro above to get the weightbasesurvage and weightbetassurvage data below;
options source notes; 
ods select all; 

*ADL;
data adl_risk (keep=hhidpn risk_t10_adl R5WTRESP RAEHSAMP RAESTRAT);
 if ( _n_ eq 1 ) then set weightbasesurvavgage_adl;
 if ( _n_ eq 1 ) then set Weightbetasurvavgage_adl ;
 set adl ;
 XBETA=age3_lin*age3_lin_beta+age3_sp1*age3_sp1_beta+age3_sp2*age3_sp2_beta;
 risk_t10_adl=1-(So_ADL_t10**exp(XBETA));
proc sort; by hhidpn; run; /* 6001 observations and 5 variables*/

*Walk;
data walk_risk (keep=hhidpn risk_t10_walk );
 if ( _n_ eq 1 ) then set weightbasesurvavgage_walk;
 if ( _n_ eq 1 ) then set Weightbetasurvavgage_walk ;
 set walk ;
 XBETA=age3_lin*age3_lin_beta+age3_sp1*age3_sp1_beta+age3_sp2*age3_sp2_beta;
 risk_t10_walk=1-(So_WALK_t10**exp(XBETA));
proc sort; by hhidpn; run; /* 6409 observations and 2 variables*/

*Death;
data death_risk (keep=hhidpn risk_t10_death );
 if ( _n_ eq 1 ) then set weightbasesurvavgage_death;
 if ( _n_ eq 1 ) then set Weightbetasurvavgage_death ;
 set death ;
 XBETA=age3_lin*age3_lin_beta+age3_sp1*age3_sp1_beta+age3_sp2*age3_sp2_beta;
 risk_t10_death=1-(So_DEATH_t10**exp(XBETA));
proc sort; by hhidpn; run; /* 6646 observations and 2 variables*/

*Merge 3 risk data;
data clinical.risks_all_ageonlymodel;
 merge adl_risk walk_risk death_risk;
 by hhidpn;
proc sort; by hhidpn; run; /* 6646 observations and 5+2=7 variables*/


***********************************************************************************************;
*2) Compute weighted SD of predicted risk/outcome;

*Export to Stata to get the weighted SD;
PROC EXPORT DATA= clinical.risks_all_ageonlymodel
            OUTFILE= "path\statadata\risks_all_ageonlymodel.dta" 
            DBMS=STATA REPLACE;
RUN;


/*From Stata:
-------------------------------------
             |       Mean   Std. dev.
-------------+-----------------------
risk_t10_adl |   .4394924    .1126508
-------------------------------------

-------------------------------------
             |       Mean   Std. dev.
-------------+-----------------------
risk_t10_w~k |   .2882832     .066562
-------------------------------------

-------------------------------------
             |       Mean   Std. dev.
-------------+-----------------------
risk_t10_d~h |   .5219192    .1998384
-------------------------------------

*/

***********************************************************************************************;
*3) Get predicted risk/outcome for all Ids in the cohort using the full model, with their original covariate values at 10 years;

*Get data with Predicted CIF for each outcome;
data cumincpredadl10;
 set clinical.cumincpredadl10;
proc sort; by newid; run;
data cumincpredwalk10;
 set clinical.cumincpredwalk10;
proc sort; by newid; run;
data survpreddth10;
 set clinical.survpreddth10;
proc sort; by newid; run;

*Merge 3 datasets;
data cuminc(drop=time_adldepdth2 time_walkdepdth2 time2death s_death);
 merge cumincpredadl10 cumincpredwalk10 survpreddth10;
 by newid;
proc sort; by newid; run;
/*6646 obs and 2+1+1=4 variables*/

*Get the age and weight of each Id;
data dataclinical_gdr_20210624; 
 set clinical.clinicalHRSdeident_gdr_20210624 (keep=newid age3_lin R5WTRESP);
run;
/* 6646 observations and 3 variables */

data cuminc;
 merge cuminc dataclinical_gdr_20210624;
 by newid;
proc sort; by age3_lin; run;
/*6646 obs and 4+2=6 variables*/

proc means data=cuminc n nmiss min max; var age3_lin; run; /*range: 70-102*/

*Get average by age prediction;
data avg_risks_by_age70_100;
 set clinical.avg_risks_by_age70_100 (keep=age3_lin risk_t10_adl risk_t10_walk risk_t10_death) ;
run;
/* 31 observations and 4 variables. */

*Merge avg_risks_by_age70_100 and cuminc;
data cuminc2;
 merge avg_risks_by_age70_100 cuminc;
 by age3_lin;
run;
/*6646 obs and 6+3=9 variables*/

proc contents data=cuminc2; run;
proc means data=cuminc2 n nmiss min max; run;

***********************************************************************************************;
*4) Compute z-score=x-mu/SD
where x: value in 3
      mu: average risk for age (from age-only model)
      SD: computed in 2 ;

*5) Classify Ids as:
 lower: z-score<-2SD
 similar: -2SD<=z-score<=2SD
 higher: z-score>2SD;

/*From Stata:
-------------------------------------
             |       Mean   Std. dev.
-------------+-----------------------
risk_t10_adl |   .4394924    .1126508
-------------------------------------

-------------------------------------
             |       Mean   Std. dev.
-------------+-----------------------
risk_t10_w~k |   .2882832     .066562
-------------------------------------

-------------------------------------
             |       Mean   Std. dev.
-------------+-----------------------
risk_t10_d~h |   .5219192    .1998384

*/


data cuminc3;
 set cuminc2;
 length risk_group2SD_adl risk_group2SD_walk risk_group2SD_death $10;

 if cif_adl ne . and risk_t10_adl ne . then z_score_adl=(cif_adl-risk_t10_adl)/0.1126508;
 if cif_walk ne . and risk_t10_walk ne . then z_score_walk=(cif_walk-risk_t10_walk)/.066562;
 if cif_death ne . and risk_t10_death ne . then z_score_death=(cif_death-risk_t10_death)/.1998384;

 if z_score_adl ne . then do;
  if z_score_adl<-2*0.1126508 then risk_group2SD_adl="1.lower";
  else if -2*0.1126508<=z_score_adl<=2*0.1126508 then risk_group2SD_adl="2.similar";
  else if z_score_adl>2*0.1126508 then risk_group2SD_adl="3.higher";
 end;

 if z_score_walk ne . then do;
  if z_score_walk<-2*0.066562 then risk_group2SD_walk="1.lower";
  else if -2*0.066562<=z_score_walk<=2*0.066562 then risk_group2SD_walk="2.similar";
  else if z_score_walk>2*0.066562 then risk_group2SD_walk="3.higher";
 end;

 if z_score_death ne . then do;
  if z_score_death<-2*0.1998384 then risk_group2SD_death="1.lower";
  else if -2*0.1998384<=z_score_death<=2*0.1998384 then risk_group2SD_death="2.similar";
  else if z_score_death>2*0.1998384 then risk_group2SD_death="3.higher";
 end;

run;
/*6646 obs and 9+6=15 variables*/

/*QC*/
proc means data=cuminc3 n nmiss min max; where age3_lin<=100; run;
proc freq data=cuminc3; tables risk_group2SD_adl risk_group2SD_walk risk_group2SD_death; where age3_lin<=100; run;

proc means data=cuminc3 n nmiss min max; var cif_adl ; class age3_lin risk_group2SD_adl ; where age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var cif_walk ; class age3_lin risk_group2SD_walk ; where age3_lin<=100; run;
proc means data=cuminc3 n nmiss min max; var cif_death ; class age3_lin risk_group2SD_death ; where age3_lin<=100; run;


/*Compute the 2x2 tables*/

*Death 2SD;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group2SD_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_death ne . and risk_t10_death ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group2SD_death / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_death ne . and risk_t10_death ne .; 
run;


*ADL 2SD;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group2SD_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_adl ne . and risk_t10_adl ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group2SD_adl / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_adl ne . and risk_t10_adl ne .; 
run;


*Walk 2SD;
proc freq data=cuminc3;
 weight R5WTRESP;
 tables age3_lin*risk_group2SD_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_walk ne . and risk_t10_walk ne .; 
run;

proc freq data=cuminc3;
 tables age3_lin*risk_group2SD_walk / nocol nopercent ;
 where age3_lin in (70, 75, 80, 85, 90, 95, 100) and cif_walk ne . and risk_t10_walk ne .; 
run;


/******************************************************************************************************************************************************************************/
*Categorize each participant as "lower than average", "about average" or "higher than average" using weighted p25 and p75 from CIF obtained from final model: meeting 5.19.22 ;


*Get age and survey variables;
data dataclinical_gdr_20210624 (keep=newid RAEHSAMP RAESTRAT R5WTRESP age3_lin subgroup_adldep subgroup_walkdep subgroup_death
                                subgroup_adldep subgroup_walkdep rename=(subgroup_adldep=subgroup_adl subgroup_walkdep=subgroup_walk));
 set clinical.clinicalHRSdeident_gdr_20210624 ; 
 subgroup_death=1;
run;
/* 6646 observations and 8 variables */

*Get predicted risk/outcome for all Ids in the cohort using the full model, with their original covariate values at 10 years;
data cumincpredadl10;
 set clinical.cumincpredadl10;
proc sort; by newid; run;
data cumincpredwalk10;
 set clinical.cumincpredwalk10;
proc sort; by newid; run;
data survpreddth10;
 set clinical.survpreddth10;
proc sort; by newid; run;

*Merge 3 datasets;
data cuminc(drop=time_adldepdth2 time_walkdepdth2 time2death s_death);
 merge cumincpredadl10 cumincpredwalk10 survpreddth10;
 by newid;
proc sort; by newid; run;
/*6646 obs and 2+1+1=4 variables*/

*Merge CIF and age data;
data cuminc;
 merge cuminc dataclinical_gdr_20210624;
 by newid;
proc sort; by age3_lin; run;
/*6646 obs and 4+7=11 variables*/

proc means data=cuminc n nmiss min max; var age3_lin; run; /*range: 70-102*/

*Group age 90+;
data cuminc2;
 set cuminc;
 if age3_lin>100 then delete; /*3 ids deleted*/
 else if age3_lin>=90 then age3_lin=90;
proc sort; by age3_lin; run;
/*6646-3=6643 obs and 11 variables*/

proc means data=cuminc2 n nmiss min max; var age3_lin; run; /*range: 70-90*/

**************************************************************** ADL;
*Computed weighted p25 and p75 by age;
proc surveymeans data=cuminc2  percentile=(25 75) ;
 strata RAESTRAT;
 cluster RAEHSAMP;
 weight R5WTRESP;
 var cif_adl;
 domain age3_lin ;
 ods output domainquantiles=quant_adl;
run;

*Reshape from long to wide;
data quant_adl2;
 set quant_adl (keep=age3_lin Percentile Estimate);
 if Percentile=25 then number=1;
 else if Percentile=75 then number=2;
run;

data quant_adl3 (keep=age3_lin p25_adl p75_adl);
 set quant_adl2;
 by age3_lin;
 retain p25_adl p75_adl;

 array est(2) p25_adl p75_adl;

 if first.age3_lin then do;
  do i=1 to 2;
    est(i)=.;
  end;
 end;

 est(number)=Estimate;
 
 if last.age3_lin;

run;

*Get the N for each age-group;
proc means data=cuminc2  n ;
 var cif_adl;
 class age3_lin ;
 output out=n_adl;
run;
data n_adl (keep=age3_lin cif_adl rename=(cif_adl=n_adl));
 set n_adl;
 where _stat_="N" and _type_ ne 0;
run;
proc print data=n_adl; run;

*Merge percentiles and N;
data quant_adl4;
 merge quant_adl3 n_adl;
 by age3_lin;
run;


**************************************************************** Walk;
*Computed weighted p25 and p75 by age;
proc surveymeans data=cuminc2  percentile=(25 75) ;
 strata RAESTRAT;
 cluster RAEHSAMP;
 weight R5WTRESP;
 var cif_walk;
 domain age3_lin ;
 ods output domainquantiles=quant_walk;
run;

*Reshape from long to wide;
data quant_walk2;
 set quant_walk (keep=age3_lin Percentile Estimate);
 if Percentile=25 then number=1;
 else if Percentile=75 then number=2;
run;

data quant_walk3 (keep=age3_lin p25_walk p75_walk);
 set quant_walk2;
 by age3_lin;
 retain p25_walk p75_walk;

 array est(2) p25_walk p75_walk;

 if first.age3_lin then do;
  do i=1 to 2;
    est(i)=.;
  end;
 end;

 est(number)=Estimate;
 
 if last.age3_lin;

run;

*Get the N for each age-group;
proc means data=cuminc2  n ;
 var cif_walk;
 class age3_lin ;
 output out=n_walk;
run;
data n_walk (keep=age3_lin cif_walk rename=(cif_walk=n_walk));
 set n_walk;
 where _stat_="N" and _type_ ne 0;
run;
proc print data=n_walk; run;

*Merge percentiles and N;
data quant_walk4;
 merge quant_walk3 n_walk;
 by age3_lin;
run;

**************************************************************** Death;
*Computed weighted p25 and p75 by age;
proc surveymeans data=cuminc2  percentile=(25 75) ;
 strata RAESTRAT;
 cluster RAEHSAMP;
 weight R5WTRESP;
 var cif_death;
 domain age3_lin ;
 ods output domainquantiles=quant_death;
run;

*Reshape from long to wide;
data quant_death2;
 set quant_death (keep=age3_lin Percentile Estimate);
 if Percentile=25 then number=1;
 else if Percentile=75 then number=2;
run;

data quant_death3 (keep=age3_lin p25_death p75_death);
 set quant_death2;
 by age3_lin;
 retain p25_death p75_death;

 array est(2) p25_death p75_death;

 if first.age3_lin then do;
  do i=1 to 2;
    est(i)=.;
  end;
 end;

 est(number)=Estimate;
 
 if last.age3_lin;

run;

*Get the N for each age-group;
proc means data=cuminc2  n ;
 var cif_death;
 class age3_lin ;
 output out=n_death;
run;
data n_death (keep=age3_lin cif_death rename=(cif_death=n_death));
 set n_death;
 where _stat_="N" and _type_ ne 0;
run;
proc print data=n_death; run;

*Merge percentiles and N;
data quant_death4;
 merge quant_death3 n_death;
 by age3_lin;
run;

**************************************************************** All;
*Merge all 3 data;
data clinical.p25_p75_by_age70_90;
 merge quant_adl4 quant_walk4 quant_death4;
 by age3_lin;
run;


*Export as csv file;
PROC EXPORT DATA= clinical.p25_p75_by_age70_90
            OUTFILE= "path\p25_p75_by_age70_90.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;
