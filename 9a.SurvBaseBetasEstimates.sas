/*********** Purpose: Compute weighted betas and baseline functions for Cox and Competing-risk final and unavailable predictor models  ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2022.03.21 ***********/
/*********** Date completed: 2022.03.22 ***********/

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
*Prepare data for macro ;

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

*Create dummy variables for categorical variables;
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
run;
/* 6646 observations and 72+6=78 variables */


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

proc means data=covall_base1obs; run;

*Create data with only 1 variable containing predictors in final model and replace: ageint5,bmi4g,heart3g,stroke3g;
data final (keep=VARINMODEL2 rename=(VARINMODEL2=FINALMODEL));
 set outdata2.baTCIC_1UnavailPred (keep=VARINMODEL DELEVAR);
 where DELEVAR="";
 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "heart3g", "heart1 heart2");
 VARINMODEL2=tranwrd(VARINMODEL2, "smoking3g", "smoking1 smoking2");
 VARINMODEL2=tranwrd(VARINMODEL2, "stroke3g", "stroke1 stroke2");
run;

*Read data with final model (16 predictors) and 14 models (15 predictors) with 1 unavailable predictor and replace: ageint5,bmi4g,heart3g,stroke3g ;
data unavail1 (keep=VARINMODEL2 DELEVAR2 nummissvar rename=(VARINMODEL2=VARINMODEL DELEVAR2=DELEVAR)) ;
 if ( _n_ eq 1 ) then set final;
 set outdata2.baTCIC_1UnavailPred (keep=VARINMODEL DELEVAR);

 length DELEVAR2 $ 200;

 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "heart3g", "heart1 heart2");
 VARINMODEL2=tranwrd(VARINMODEL2, "smoking3g", "smoking1 smoking2");
 VARINMODEL2=tranwrd(VARINMODEL2, "stroke3g", "stroke1 stroke2");

 if DELEVAR="" then nummissvar=0;
 else nummissvar=1; 

 array x{23}  $ 32; /*there are 23 possible coefficients for 16 predictors in final model*/

 n=0;
 /*Fill in array x with all the coefficients in each reduced model */
 do i=1 to countw(VARINMODEL2,' ');
  temp=scan(VARINMODEL2,i,' ');
  if temp not in x then do; n+1; x{n}=temp; end;
 end;

 /*If the predictor from the final model is not found in array x then the predictor was deleted */
 do i=1 to countw(FINALMODEL,' ');
  temp=scan(FINALMODEL,i,' ');
  if temp not in x then DELEVAR2=catx(' ',DELEVAR2,temp);
 end;

run;
proc print data=unavail1; run;


*Read data with 91 models (14 predictors) with 2 unavailable predictor and replace: ageint5,bmi4g,heart3g,stroke3g ;
data unavail2 (keep=VARINMODEL2 DELEVAR2 nummissvar rename=(VARINMODEL2=VARINMODEL DELEVAR2=DELEVAR)) ;
 if ( _n_ eq 1 ) then set final;
 set outdata2.baTCIC_2UnavailPred (keep=VARINMODEL DELEVAR);

 length DELEVAR2 $ 200;

 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "heart3g", "heart1 heart2");
 VARINMODEL2=tranwrd(VARINMODEL2, "smoking3g", "smoking1 smoking2");
 VARINMODEL2=tranwrd(VARINMODEL2, "stroke3g", "stroke1 stroke2");

 nummissvar=2;

 array x{23}  $ 32; /*there are 23 possible coefficients for 16 predictors in final model*/

 n=0;
 /*Fill in array x with all the coefficients in each reduced model */
 do i=1 to countw(VARINMODEL2,' ');
  temp=scan(VARINMODEL2,i,' ');
  if temp not in x then do; n+1; x{n}=temp; end;
 end;

 /*If the predictor from the final model is not found in array x then the predictor was deleted */
 do i=1 to countw(FINALMODEL,' ');
  temp=scan(FINALMODEL,i,' ');
  if temp not in x then DELEVAR2=catx(' ',DELEVAR2,temp);
 end;

run;
proc print data=unavail2; run;

data unavail;
 set unavail1 unavail2;
 modelnum=_N_;
run;
/*15+91=106 models*/
proc print data=unavail; run;

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



/******************************************************************************************************************************************************************/
*Define macro variables;

proc sql noprint; select max(modelnum) format 3. into :modelnum from unavail; quit; /*create macro variable with total number of models*/
%put "&modelnum"; 

%let NUMOUTCOMES=3; /*number of outcomes*/
%let ALLOUTCOME=status_adldepdth2 status_walkdepdth2 death;
%let ALLTIME=time_adldepdth2 time_walkdepdth2 time2death;
%let ALLLABEL= ADL WALK DEATH; /*labels for outcomes. I put them in capital letters so that it's recognized in proc sql with memname option*/



/******************************************************************************************************************************************************************/
*Macro ;


/*Test macro: Cox model Death*/

*Compute weighted survival probability for all Ids via METHOD=PL; 

*1) Within SAS: Compute weighted betas (dataset called “WeightedBetasCoxDeath”):  ;
proc phreg data = dataclinical_gdr_20210624 outest=WeightBetaSurvDth; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = &varname /risklimits;
run;

*2) Within SAS: Obtain So(t) (i.e. So_death_weighted_pl using weighted betas above and METHOD=PL;
proc phreg data = dataclinical_gdr_20210624 inest=WeightBetaSurvDth ; 
  model time2death*death(0) = &varname /risklimits maxit=0;
 baseline out=WeightBaseSurvDth (keep=time2death So_death) covariates=covall_base1obs survival=So_death timelist=1 to 19 by 1 / method=PL;
run;
/*WeightBaseSurvDth has 19 observations and 2 variables.*/

*Keep variables of interest in clinical.BetasCoxDeath data;
data WeightBetaSurvDth;
 set WeightBetaSurvDth (keep=age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart1 heart2 R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA
                       R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking1 smoking2 stroke1 stroke2);
run;

*Rename variables in clinical.WEightedBetasCoxDeath;
/* This code creates a macro variable &list with the list of variables in the form. */
/* id = id_beta                                                                     */
/* This format could be used to add a suffix to all the variables.                  */
proc sql noprint;
   select cats(name,'=',name,'_beta')
          into :list
          separated by ' '
          from dictionary.columns
          where libname = 'WORK' and memname = 'WEIGHTBETASURVDTH';
quit;
%put "&list";

proc datasets library = work nolist;
   modify WeightBetaSurvDth;
   rename &list;
quit;

*Traspose data from long to wide format;
data WeightBaseSurvDth (keep=So_death_t1-So_death_t19);
 set WeightBaseSurvDth end=last;
 array surv[19] So_death_t1-So_death_t19;
 retain So_death_t1-So_death_t19;
 do i=1 to 19;
  if time2death=i then surv[i]=So_death;
 end; drop i;
 if last;
run;


data unavail;
 set unavail;
 where modelnum=1;
run;
%let NUMOUTCOMES=1; /*number of outcomes*/
%let ALLOUTCOME=death;
%let ALLTIME=time2death;
%let ALLLABEL= DEATH; /*labels for outcomes*/

/*Finish test*/


options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors/warnings to be written to the log*/
options nomlogic nomprint nomrecall nosymbolgen;
ods select none; /*to create output data sets through the ODS OUTPUT statement and suppress the display of all output*/

%macro betas_basesurv;
 %do i=1 %to &modelnum;
 /*For each model define VARNAME as the variables selected in the model*/
  data _null_;
   set unavail;
   where modelnum=&i;
   call symputx ('VARNAME' ,VARINMODEL);
  run;

  %do j=1 %to &NUMOUTCOMES;
    %let OUTCOME=%scan(&ALLOUTCOME,&j); /*extract the jth outcome, jth time, jth label*/
    %let TIME=%scan(&ALLTIME,&j);
	%let LABEL=%scan(&ALLLABEL,&j);

	*1) Within SAS: Compute weighted betas (dataset called “WeightBetaSurv_&label”):  ;
	proc phreg data = dataclinical_gdr_20210624 outest=WeightBetaSurv_&label; 
	  weight R5WTRESP; /*use survey weight*/
	  model &time*&outcome(0) = &varname /risklimits;
	  where subgroup_&label=1;
	run;

	*2) Within SAS: Obtain So(t) (i.e. So_&label using above weighted betas and METHOD=PL;
	proc phreg data = dataclinical_gdr_20210624 inest=WeightBetaSurv_&label ; 
	  model &time*&outcome(0) = &varname /risklimits maxit=0;
	 baseline out=WeightBaseSurv_&label (keep=&time So_&label) covariates=covall_base1obs survival=So_&label timelist=1 to 19 by 1 / method=PL;
	run;
	/*WeightBaseSurv_&label has 19 observations and 2 variables.*/

	*Keep variables of interest in WeightBetaSurv_&label;
	data WeightBetaSurv_&label;
	 set WeightBetaSurv_&label (keep=&varname);
	run;

	*Rename variables in WeightBetaSurv_&label;
	/* This code creates a macro variable &list with the list of variables in the form. */
	/* variable = variable_beta                                                         */
	/* This format could be used to add a suffix to all the variables.                  */
	proc sql noprint;
	   select cats(name,'=',name,'_beta')
	          into :list
	          separated by ' '
	          from dictionary.columns
	          where libname = 'WORK' and memname = "WEIGHTBETASURV_&label";
	quit;

	proc datasets library = work nolist;
	   modify WEIGHTBETASURV_&label;
	   rename &list;
	quit;

	/*Add model number so we can later merge it with data with model information: VARINMODEL DELEVAR*/
	data WeightBetaSurv_&label;
	 set WeightBetaSurv_&label;
     modelnum=&i;
	run;

	*Traspose Base survival data from long to wide format;
	data WeightBaseSurv_&label  (keep=So_&label._t1-So_&label._t19);
	 set WeightBaseSurv_&label  end=last;
	 array surv[19] So_&label._t1-So_&label._t19;
	 retain So_&label._t1-So_&label._t19;
	 do i=1 to 19;
	  if &time=i then surv[i]=So_&label;
	 end; drop i;
	 if last;
	run;

	/*Add model number so we can later merge it with data with model information: VARINMODEL DELEVAR*/
    data WeightBaseSurv_&label;
	 set WeightBaseSurv_&label;
     modelnum=&i;
	run;

	proc append base=WeightBetaSurvAll_&label data=WeightBetaSurv_&label force; run;
    proc delete data=WeightBetaSurv_&label; run; quit;

    proc append base=WeightBaseSurvAll_&label data=WeightBaseSurv_&label force; run;
    proc delete data=WeightBaseSurv_&label; run; quit;

  %end; /*j loop*/


%end; /*i loop*/

%mend betas_basesurv;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%betas_basesurv
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

/*
======MONITORING: 2022-03-22, 11:56======
======MONITORING: 2022-03-22, 11:58======
*/



/******************************************************************************************************************************************************************/
*Create final data ;

%macro finaldata;

%do j=1 %to &NUMOUTCOMES;

 %let OUTCOME=%scan(&ALLOUTCOME,&j); /*extract the jth outcome, jth time, jth label*/
 %let TIME=%scan(&ALLTIME,&j);
 %let LABEL=%scan(&ALLLABEL,&j);

 data clinical.Weightbetasurvall_&label;
  merge unavail Weightbetasurvall_&label;
  by modelnum;
 run;
 data clinical.Weightbasesurvall_&label;
  merge unavail Weightbasesurvall_&label;
  by modelnum;
 run;
 PROC EXPORT DATA= clinical.Weightbasesurvall_&label
            OUTFILE= "path\Rfiles\Weightbasesurvall_&label..csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;
 PROC EXPORT DATA= clinical.Weightbetasurvall_&label
            OUTFILE= "path\Rfiles\Weightbetasurvall_&label..csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;

 %end; /*j loop*/

%mend finaldata;

%finaldata


