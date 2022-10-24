/*********** Purpose: Compute IAUC of Cox and Competing-risk final and unvailable predictor models for full sample and subgroups ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2022.04.07 ***********/
/*********** Date completed: 2022.04.09 ***********/

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

proc contents data=clinical.clinicalHRSdeident_gdr_20210624; run;

data dataclinical_gdr_20210624; 
 set clinical.clinicalHRSdeident_gdr_20210624 (keep=newid RAEHSAMP RAESTRAT clustervar R5WTRESP race3 ageint5 age_lin age_sp1 age_sp2 age_sp3 age3_lin age3_sp1 age3_sp2
                                gender R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g bmi_lin bmi_sp1 bmi_sp2 bmi_sp3 bmi3_lin bmi3_sp1 bmi3_sp2
                                smoking3g stroke3g volunteer R5PROXY status_adldepdth2 status_walkdepdth2 death time_adldepdth2 time_walkdepdth2 time2death
								normwgt_adl normwgt_walk normwgt_death R5BMI
                                subgroup_adldep subgroup_walkdep rename=(subgroup_adldep=subgroup_adl subgroup_walkdep=subgroup_walk)); 
 subgroup_death=1;
run;
/* 6646 observations and 73 variables */

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
/* 6646 observations and 73+6=79 variables */

*Check N for subggroups;
proc means data=dataclinical_gdr_20210624 n; var ageint5 gender race3; run; /*race3 2 missing*/

*Death;
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5<80; run; /*4343/6646=65.35*/
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5>=80; run; /*2303/6646=34.65*/

proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5<75; run; /*2314/6646=34.82*/
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5>=75; run; /*4332/6646=65.18*/


proc freq data=dataclinical_gdr_20210624; tables gender race3; run; 
/*female=3841 (57.79%), male=2805 (42.21%), white=5397 (81.23%), black=722 (10.87%), 525 (7.90%), miss=2 */ 


*ADL;
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where subgroup_adl=1 ; run; 
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5<80 and subgroup_adl=1 ; run; /*4033/6001=67.21*/
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5>=80 and subgroup_adl=1; run; /*1968/6001=32.79*/

proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5<75 and subgroup_adl=1 ; run; /*2169/6001=36.14*/
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5>=75 and subgroup_adl=1; run; /*3832/6001=63.86*/


proc freq data=dataclinical_gdr_20210624; tables gender race3; where subgroup_adl=1 ; run; 
/*female=3436 (57.26%) , male=2565 (42.74%) , white=4955 (82.60%) , black=608 (10.14%) , 436 (7.27%) , miss=2 */ 


*Walk;
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where subgroup_walk=1 ; run; 
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5<80 and subgroup_walk=1 ; run; /*4233/6409=66.05*/
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5>=80 and subgroup_walk=1; run; /*2176/6409=33.95*/

proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5<75 and subgroup_walk=1 ; run; /*2263/6409=35.31*/
proc means data=dataclinical_gdr_20210624 n min max; var ageint5 ; where ageint5>=75 and subgroup_walk=1; run; /*4146/6409=64.69*/


proc freq data=dataclinical_gdr_20210624; tables gender race3; where subgroup_walk=1 ; run; 
/*female=3692 (57.61%)  , male=2717 (42.39%)  , white=5231 (81.65%)  , black=684 (10.68%)  , 492 (7.68%)  , miss=2 */ 


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
data unavail1 (keep=VARINMODEL2 DELEVAR2 nummissvar TCIC_avg model rename=(VARINMODEL2=VARINMODEL DELEVAR2=DELEVAR)) ;
 if ( _n_ eq 1 ) then set final;
 set outdata2.baTCIC_1UnavailPred (keep=VARINMODEL DELEVAR TCIC_avg);

 length DELEVAR2 model $ 200;

 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "heart3g", "heart1 heart2");
 VARINMODEL2=tranwrd(VARINMODEL2, "smoking3g", "smoking1 smoking2");
 VARINMODEL2=tranwrd(VARINMODEL2, "stroke3g", "stroke1 stroke2");

 if DELEVAR="" then do; nummissvar=0; model="final"; end;
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


*Keep model Best-case scenario for 1 missing predictor (minimum BIC);
data unavail1_nofinal; set unavail1; where model ne "final"; proc sort; by descending TCIC_avg ; run;
data unavail1_best (drop=TCIC_avg nummissvar);
 set unavail1_nofinal end=last;
 if last;
 model="best_1miss";
run; 
proc print data=unavail1_best; run;

*Keep model Worst-case scenario for 1 missing predictor (maximum BIC);
proc sort data=unavail1_nofinal; by TCIC_avg ; run;
data unavail1_worst (drop=TCIC_avg nummissvar);
 set unavail1_nofinal end=last;
 if last;
 model="worst_1miss";
run; 
proc print data=unavail1_worst; run;

*Stack data with final model, unavail1_best, unavail1_worst;
data final_model; set unavail1 (drop=TCIC_avg nummissvar); where model = "final"; run;

data unavail1_3models; 
 retain model DELEVAR VARINMODEL;
 set final_model unavail1_best unavail1_worst;
run;
proc print data=unavail1_3models; run;

proc delete data=unavail1_nofinal unavail1_best unavail1_worst final_model; run;quit;


*Read data with 91 models (14 predictors) with 2 unavailable predictor and replace: ageint5,bmi4g,heart3g,stroke3g ;
data unavail2 (keep=VARINMODEL2 DELEVAR2 nummissvar model TCIC_avg rename=(VARINMODEL2=VARINMODEL DELEVAR2=DELEVAR)) ;
 if ( _n_ eq 1 ) then set final;
 set outdata2.baTCIC_2UnavailPred (keep=VARINMODEL DELEVAR TCIC_avg);

 length DELEVAR2 model $ 200;

 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "heart3g", "heart1 heart2");
 VARINMODEL2=tranwrd(VARINMODEL2, "smoking3g", "smoking1 smoking2");
 VARINMODEL2=tranwrd(VARINMODEL2, "stroke3g", "stroke1 stroke2");

 nummissvar=2; model="";

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

*Keep model Best-case scenario for 2 missing predictor (minimum BIC);
proc sort data=unavail2; by descending TCIC_avg ; run;
data unavail2_best (drop=TCIC_avg nummissvar);
 set unavail2 end=last;
 if last;
 model="best_2miss";
run; 

*Keep model Worst-case scenario for 2 missing predictor (maximum BIC);
proc sort data=unavail2; by TCIC_avg ; run;
data unavail2_worst (drop=TCIC_avg nummissvar);
 set unavail2 end=last;
 if last;
 model="worst_2miss";
run; 

*Stack data unavail2_best, unavail2_worst;
data unavail2_2models; 
 retain model DELEVAR VARINMODEL;
 set unavail2_best unavail2_worst;
run;

proc delete data=unavail2_best unavail2_worst ; run;quit;


*Stack data with all 5 models;
data unavail;
 set unavail1_3models unavail2_2models;
 modelnum=_N_;
run;
/*3+2=5 models*/
proc print data=unavail; run;

proc delete data=unavail1_3models unavail2_2models unavail1 unavail2 final; run; quit;

*Add where conditions as columns to unavail data;
data unavail2;
 set unavail;
 length condition1-condition10 $100;
 condition1=""; /*full sample*/
 condition2="and ageint5<80"; 
 condition3="and ageint5>=80"; 
 condition4="and ageint5<75"; 
 condition5="and ageint5>=75"; 
 condition6="and gender=0"; /*female*/
 condition7="and gender=1"; /*male*/
 condition8="and race3=0" ;  /*0.white, 1.black, 2.hispanic/other*/
 condition9="and race3=1" ; 
 condition10="and race3=2" ; 
run;

proc print data=unavail2; run;
proc freq data=unavail2; tables condition1-condition10; run;

/******************************************************************************************************************************************************************/
*Define macro variables;

proc sql noprint; select max(modelnum) format 3. into :modelnum from unavail; quit; /*create macro variable with total number of models*/
%put "&modelnum"; 

%let NUMOUTCOMES=3; /*number of outcomes*/
%let ALLOUTCOME=status_adldepdth2 status_walkdepdth2 death;
%let ALLTIME=time_adldepdth2 time_walkdepdth2 time2death;
%let ALLLABEL= adl walk death; /*labels for outcomes*/
%let NUMCONDITIONS=10; /*number of subgroup analyses*/
%let ALLCONDLABEL= full less80 80plus less75 75plus female male white black hisp_other;


/******************************************************************************************************************************************************************/
*Macro ;

/*Test macro ADL: subgroup 80+, less75 has age3_sp2 beta=0*/

data _null_;
   set unavail2;
   where modelnum=1;
   call symputx ('VARNAME' ,VARINMODEL);
   call symputx ('cond1' ,condition1);
   call symputx ('cond2' ,condition2);
   call symputx ('cond3' ,condition3);
   call symputx ('cond4' ,condition4);
   call symputx ('cond5' ,condition5);
   call symputx ('cond6' ,condition6);
   call symputx ('cond7' ,condition7);
   call symputx ('cond8' ,condition8);
   call symputx ('cond9' ,condition9);
   call symputx ('cond10' ,condition10);
  run;

%let OUTCOME=%scan(&ALLOUTCOME,1); %put &outcome; /*extract the jth outcome, jth time, jth label*/
%let TIME=%scan(&ALLTIME,1);  %put &time;
%let LABEL=%scan(&ALLLABEL,1);  %put &label;

%let k=10;
%let condlab=%scan(&ALLCONDLABEL,&k);  %put &condlab;
%let condition=&&&cond&k; %put &condition; 
/*
&&&: Multiple ampersands can be used to allow the value of a macro variable to become another macro variable reference. 
The macro variable reference will be rescanned until the macro variable is resolved.
On the first scan:
[&&]: resolves to &
[&cond&k]: resolves to &cond1
On the second scan: &cond1 resolves to the actual condition1
*/

/*proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); --> we don't need this option since using only the weighted xb not the Standard errors to compute IAUC */

proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model &time*&outcome(0) = &varname /eventcode=1 risklimits;
  where subgroup_&label=1 &condition;
  ods output FailureSummary=CensUncens; 
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if &outcome=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    &outcome=0;
    &time=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007);
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc); 
  model &time*&outcome(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'Weighted'  pred=xb;
  ods output IAUC=iauc;
run;

*Merge data with IAUC and N;
data &label._&condlab; 
 merge iauc (keep=estimate rename=(estimate=iauc_&label._&condlab)) CensUncens(keep=total rename=(total=N_&label._&condlab) );
 modelnum=1;
run;
proc delete data=bsout iauc CensUncens; run; quit;

*if &k=1 then create outcome data;
/*data &label; set &label._&condlab; run; */
/*proc delete data=&label._&condlab; run; quit;*/

*else if &k ne 1 then merge;
data &label;
 merge &label &label._&condlab;
 by modelnum;
run;
proc delete data=&label._&condlab; run; quit;


*Check why subgroups 80+ and less75 have age3_sp2 beta=0. 4 defaults knots from R rcspline.eval function are: 70 75 79 89;

*less75;
proc means data=dataclinical_gdr_20210624; var age3_lin age3_sp1 age3_sp2; where ageint5<75; run;
*All Ids ageint5<75 have age3_sp2=0;

proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender / risklimits;
  where ageint5<75;
  output out=BSOUT xbeta=xb;
run;

*Change order of splines: same result as above;
proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_sp1 age3_sp2 age3_lin gender / risklimits;
  where ageint5<75;
  output out=BSOUT xbeta=xb;
run;

*75plus: no problem;
proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender / risklimits;
  where ageint5>=75;
  output out=BSOUT xbeta=xb;
run;

 *80+;
proc means data=dataclinical_gdr_20210624; var age3_lin age3_sp1 age3_sp2; where ageint5>=80; run;
proc means data=dataclinical_gdr_20210624; var age3_lin age3_sp1 age3_sp2; where ageint5>=80; class death; run;

*less80: no problem;
proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender / risklimits;
  where ageint5<80;
  output out=BSOUT xbeta=xb;
run;

/*
Cases:
1) age3_lin=80, sp1=2.76481994459834, sp2=0.342382271468144
2) age3_85=85, sp1=8.21218836565097, sp2=1.93240997229917
*/

data test;
 input newid age3_lin age3_sp1 age3_sp2 ;
 datalines;
 1 80 2.76481994459834 0.342382271468144
 2 85 8.21218836565097 1.93240997229917
 3 92 18 5.119113573
;
run;
proc print; run; 

*With this order, age3_sp2 is dropped because it's the last one;
proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2  / risklimits;
  where ageint5>=80;
  output out=BSOUT xbeta=xb;
  baseline out=SurvPred1 survival=s_death timelist=10 covariates=test;
run;

*Change order: age3_lin is dropped because it's the last one;
proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_sp1 age3_sp2 age3_lin   / risklimits;
  where ageint5>=80;
  output out=BSOUT xbeta=xb;
  baseline out=SurvPred2 survival=s_death timelist=10 covariates=test;
run;

*Change order: age3_sp1 is dropped because it's the last one;
proc phreg data = dataclinical_gdr_20210624 ; 
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) =  age3_lin age3_sp2 age3_sp1   / risklimits;
  where ageint5>=80;
  output out=BSOUT xbeta=xb;
  baseline out=SurvPred3 survival=s_death timelist=10 covariates=test;
run;

data survpred;
 set SurvPred1-SurvPred3;
run;

*You get exactly the same predictions;
proc print data=SurvPred; run;
proc means data=survpred; var s_death; class age3_lin; run;

*Correlation;
proc corr data=dataclinical_gdr_20210624; var age3_lin age3_sp1 age3_sp2; where ageint5>=80; run; /*super high correlation: >0.99*/
proc corr data=dataclinical_gdr_20210624; var age3_lin age3_sp1 age3_sp2; where ageint5<80; run; /*lower correlation compared with 80+ group: 0.70743-0.91586*/

proc delete data=test SurvPred1-SurvPred3 SurvPred bsout; run; quit;

/*Conclusion:
For Ids age>=80 one of the betas for the splines is dropped from the model, because for these Ids there is high correlation (>0.99) among the 3 splines,
so we have all the information that we need to fit the model with only two splines 
*/


options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors/warnings to be written to the log*/
options nomlogic nomprint nomrecall nosymbolgen;
ods select none; /*to create output data sets through the ODS OUTPUT statement and suppress the display of all output*/

%macro iauc_subgroups;
 %do i=1 %to &modelnum;
 /*For each model, define VARNAME as the variables selected in the model and cond1-cond10 for the subgroup analyses */
  data _null_;
   set unavail2;
   where modelnum=&i;
   call symputx ('VARNAME' ,VARINMODEL);
   call symputx ('cond1' ,condition1);
   call symputx ('cond2' ,condition2);
   call symputx ('cond3' ,condition3);
   call symputx ('cond4' ,condition4);
   call symputx ('cond5' ,condition5);
   call symputx ('cond6' ,condition6);
   call symputx ('cond7' ,condition7);
   call symputx ('cond8' ,condition8);
   call symputx ('cond9' ,condition9);
   call symputx ('cond10' ,condition10);
  run;

  %do j=1 %to &NUMOUTCOMES;
    %let OUTCOME=%scan(&ALLOUTCOME,&j); /*extract the jth outcome, jth time, jth label*/
    %let TIME=%scan(&ALLTIME,&j);
	%let LABEL=%scan(&ALLLABEL,&j);

    %do k=1 %to &NUMCONDITIONS; /*for each condition*/
	  %let condlab=%scan(&ALLCONDLABEL,&k); /*get the label of the condition*/
      %let condition=&&&cond&k; /*get the condition per se*/

	  %if &label ne death %then %do;
		proc phreg data = dataclinical_gdr_20210624 ; 
		  weight R5WTRESP; /*use survey weight*/
		  model &time*&outcome(0) = &varname /eventcode=1 risklimits;
		  where subgroup_&label=1 &condition;
		  ods output FailureSummary=CensUncens; 
		  output out=BSOUT xbeta=xb;
		run;
		data BSOUT; 
		  set BSOUT; 
		  if &outcome=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
		    &outcome=0;
		    &time=19.1130249;
		 end;
		run;
		*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007);
		proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc); 
		  model &time*&outcome(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
		  roc 'Weighted'  pred=xb;
		  ods output IAUC=iauc;
		run;
	 %end;

	 %else %if &label=death %then %do;
		proc phreg data = dataclinical_gdr_20210624 ; 
		  weight R5WTRESP; /*use survey weight*/
		  model &time*&outcome(0) = &varname / risklimits;
		  where subgroup_&label=1 &condition;
		  ods output CensoredSummary=CensUncens; 
		  output out=BSOUT xbeta=xb;
		run;
		*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007);
		proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc); 
		  model &time*&outcome(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
		  roc 'Weighted'  pred=xb;
		  ods output IAUC=iauc;
		run;
	 %end;

	*Merge data with IAUC and N;
	data &label._&condlab; 
	 merge iauc (keep=estimate rename=(estimate=iauc_&label._&condlab)) CensUncens(keep=total rename=(total=N_&label._&condlab) );
	 modelnum=&i;
	run;
	proc delete data=bsout iauc CensUncens; run; quit;

	%if &k=1 %then %do;
		data &label; set &label._&condlab; run; 
		proc delete data=&label._&condlab; run; quit;
    %end;
	%else %if &k ne 1 %then %do;
		data &label;
		 merge &label &label._&condlab;
		 by modelnum;
		run;
		proc delete data=&label._&condlab; run; quit;
	%end;

  %end; /*k loop: condition1-condition10*/

  proc append base=&label._all data=&label force; run;
  proc delete data=&label; run; quit;

%end; /*j loop: outcome1-outcome3*/

%end; /*i loop: model1-model5*/

%mend iauc_subgroups;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%iauc_subgroups
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

/*
======MONITORING: 2022-04-07, 16:21======
======MONITORING: 2022-04-07, 16:39======
*/



/******************************************************************************************************************************************************************/
*Create final data ;

%macro finaldata;

%do j=1 %to &NUMOUTCOMES;

 %let OUTCOME=%scan(&ALLOUTCOME,&j); /*extract the jth outcome, jth time, jth label*/
 %let TIME=%scan(&ALLTIME,&j);
 %let LABEL=%scan(&ALLLABEL,&j);

 data clinical.iauc_subgroups_&label;
  merge unavail &label._all;
  by modelnum;
 run;
 PROC EXPORT DATA= clinical.iauc_subgroups_&label
            OUTFILE= "path\iauc_subgroups_&label..csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
 RUN;

 %end; /*j loop*/

%mend finaldata;

%finaldata
