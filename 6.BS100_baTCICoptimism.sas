***********************************************************************************************************************************************************************************;
*Purpose: Compute optimism corrected stats                                                                                                                                         ;                                     
*Statistician: Grisell Diaz-Ramirez																																				   ;
*Started: 2020.02.05																																							   ;
*Finished: 2021.07.13                                                                                                                                                              ;
*Date modified: 2022.6.13                                                                                                                                                          ;
***********************************************************************************************************************************************************************************;


options MERGENOBY=ERROR;

libname harmo "path";
libname rand "path";
proc format cntlin=rand.sasfmts; run;
proc format cntlin=harmo.formats; run;

libname outdata "path";
libname outdata2 "path";
libname clinical "path";
libname savedata "path";
options nomlogic nomprint nomrecall nosymbolgen;
options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors to be written to the log*/
options noquotelenmax; /*supress warning: THE QUOTED STRING CURRENTLY BEING PROCESSED HAS BECOME MORE THAN 262 CHARACTERS LONG*/
ods select none; /*to create output data sets through the ODS OUTPUT statement and suppress the display of all output*/
options nothreads; /*to avoid Insufficient space error in Utility files*/
options errorabend; /*so the process stops as soon as an error is encountered*/


/***************************************************** Calculate Optimism Corrected C-statistic/IAUC for baTCIC models ******************************************************************/

/*
Steps:
We computed the bootstrap-based optimism-corrected C-statistic/IAUC as follows:
1) Obtain final model and corresponding C-statistic/IAUC on the case study data, namely C-statistic/IAUC-apparent
2) Obtain final models of each bootstrap sample and compute the C-statistic/IAUC of each bootstrap model, namely C-statistic/IAUC-boot
3) Compute the C-statistic/IAUC of each bootstrap model evaluated in the original case-study data, namely C-statistic/IAUC-original

In previous version, I used the variables selected in each bootstrap model to fit the model with these variables in the original data
This could be thought as bootstrap of the selection method.

A more correct way is:
"Freeze" the model obtained in 2 and obtain the C-statistic of each bootstrap model evaluated in the original data, this means:
3a) use the coefficient estimates from bootstrap model to obtain predictions in original data:
Model 1: PHREG option OUTEST=betas_boot using bootstrap data
Model 2: PHREG option INEST=betas_boot using original data and OUTPUT OUT=BSOUT xbeta=xb;

OR
Model 1: PHREG and STORE=item-store statement: requests that the fitted model be saved to an item store
Model 2: PLM restore=item-store created in 1)
         SCORE data=&DATAorig out=BSOUT predicted: score (Linear predictor) new observations based on the item store that was created in 1)

Both methods OUTEST/INTEST and STORE with PLM/SCORE give the same results

3b) use predictions in 3a) as covariate in model fitted in orginal data
Model 3: PHREG with MODEL statement and NOFIT option and ROC statement using PRED=predicted 
         This gives the same C-statistic/IAUC as having the MODEL statement with covariate predicted

3c) compute C-statitisc/IAUC of this model

4) Calculate the optimism in the fit of each bootstrap sample as: 
 C-statistic/IAUC-boot - C-statistic/IAUC-original
5) Average the optimism across 500 bootstrap sample, namely O
6) Compute the optimism-corrected C-statistic/IAUC of the case-study data as:
 C-statistic/IAUC-apparent - O
We then computed the location-shifted bootstrap confidence intervals of the optimism-corrected C-statistic/IAUC by subtracting the optimism estimate from the 2.5th and 97.5th percentiles of the C-statistic/IAUC-bootstrap distribution , , .

Note1: During model selection we used age as a linear function (ageint5), however, the baTCIC models for each bootstrap sample were fitted using
3 restricted cubic splines of age (age3_lin age3_sp1 age3_sp2).
Similarly for BMI, during model selection we used BMI as a categorical variable with 4 groups, however,
the baTCIC models for each bootstrap sample were fitted using 3 restricted cubic splines of BMI (bmi3_lin bmi3_sp1 bmi3_sp2).
The restricted cubic splines were estimated with 3 knots at default quantiles, we used the R function “rcspline.eval” in “Hmisc” package .

Note2: To avoid the computational constraint of fitting hundreds of competing-risk regression models.
The estimated C-statistic/IAUC-boot and C-statistic/IAUC-original were obtained from fitting Cox models instead of Competing-risk regression models
for ADL disability and mobility disability. To do this, we used a modified version of data where those who died were treated as being censored
at the longest possible time that any respondent was followed (i.e. 19 years) (Wolbers et al.).

*/

/*** Step 1 (doesn't change): Obtain final model and corresponding C-statistic/IAUC on the case study data, namely C-statistic-apparent ***/

*ADL;
proc phreg data=clinical.clinicalhrsdeident_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
        R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0') ;
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE
                           R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  id clustervar;
  where subgroup_adldep=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_adldepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_adldepdth2=0;
    time_adldepdth2=19.1130249;
 end;
run;
proc phreg data = BSOUT CONCORDANCE=HARRELL (SE) rocoptions(method=RECURSIVE iauc); 
  class gender  ;
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CompRiskC' pred=xb;
  ods output CONCORDANCE=concord_adl (keep=Estimate rename=(estimate=capp_adl))  IAUC=iaucr_adl (keep=Estimate rename=(estimate=iaucrapp_adl)) ;
run;
proc phreg data = BSOUT rocoptions(method=ipcw iauc); 
  class gender;
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CompRiskC'  pred=xb;
  ods output IAUC=iauc_adl (keep=Estimate rename=(estimate=iaucapp_adl));
run;
proc delete data=bsout; run; quit;

*Walk;
proc phreg data=clinical.clinicalhrsdeident_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
        R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0') ;
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE
                           R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  id clustervar;
  where subgroup_walkdep=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_walkdepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_walkdepdth2=0;
    time_walkdepdth2=19.1130249;
 end;
run;
proc phreg data = BSOUT CONCORDANCE=HARRELL (SE) rocoptions(method=RECURSIVE iauc); 
  class gender  ;
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CompRiskC' pred=xb;
  ods output CONCORDANCE=concord_walk (keep=Estimate rename=(estimate=capp_walk))  IAUC=iaucr_walk (keep=Estimate rename=(estimate=iaucrapp_walk)) ;
run;
proc phreg data = BSOUT rocoptions(method=ipcw iauc); 
  class gender;
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CompRiskC'  pred=xb;
  ods output IAUC=iauc_walk (keep=Estimate rename=(estimate=iaucapp_walk));
run;
proc delete data=bsout; run; quit;

*Death;
proc surveyphreg data=clinical.clinicalhrsdeident_gdr_20210624; 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
        R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0') ;
  cluster RAEHSAMP; /*survey cluster*/
  strata RAESTRAT; /*survey strata*/
  weight R5WTRESP; /*survey weight*/ 
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE
                           R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g / risklimits;
  output out=BSOUT xbeta=xb;
run;
proc phreg data = BSOUT CONCORDANCE=HARRELL (SE) rocoptions(method=RECURSIVE iauc); 
  class gender;
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted' pred=xb;
  ods output CONCORDANCE=concord_death (keep=Estimate rename=(estimate=capp_death))  IAUC=iaucr_death (keep=Estimate rename=(estimate=iaucrapp_death)) ;
run;
proc phreg data = BSOUT rocoptions(method=ipcw iauc); 
  class gender;
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc_death (keep=Estimate rename=(estimate=iaucapp_death));
run;
proc delete data=bsout; run; quit;

*Merge all 9 data;

options MERGENOBY=NOWARN;
data outdata.stats_app;
 merge concord_adl iaucr_adl iauc_adl concord_walk iaucr_walk iauc_walk concord_death iaucr_death iauc_death;
run;
*Calculate capp_avg and iaucrapp_avg;
data outdata.stats_app;
 set outdata.stats_app;
 capp_avg=mean(capp_adl,capp_walk,capp_death);
 iaucrapp_avg=mean(iaucrapp_adl,iaucrapp_walk,iaucrapp_death);
 iaucapp_avg=mean(iaucapp_adl,iaucapp_walk,iaucapp_death);
run;

proc delete data=concord_adl iaucr_adl iauc_adl concord_walk iaucr_walk iauc_walk concord_death iaucr_death iauc_death; run; quit;

options MERGENOBY=ERROR;

/*** Step 2(doesn't change): Obtain final models of each bootstrap sample and compute the C-statistic/IAUC of each bootstrap model, namely C-statistic/IAUC-boot ***/

*redefine VARINMODEL: 'ageint5' replace with 'age3_lin age3_sp1 age3_sp2', bmi4g replace with bmi3_lin bmi3_sp1 bmi3_sp2;
data baTCICbs; 
 set outdata.baTCIC_rep100;
 numVarsfin=countw(VARINMODEL, '');
 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 /*define VARCLASS variables*/
 VARCLASS=VARINMODEL;
 VARCLASS=tranwrd(VARCLASS, "ageint5", "");
 VARCLASS=tranwrd(VARCLASS, "bmi4g", "");
run;

/*QC*/
proc contents data=baTCICbs; run;
proc means data=baTCICbs n mean std stderr clm median p25 p75 maxdec=4; var numVarsfin; run;
proc print data=baTCICbs (obs=10); var replicate VARINMODEL VARINMODEL2 VARCLASS; run;

*Prepare BS data;
data bsdata;
 set savedata.clinical_bs500_gdr_20210707;
 where 0<replicate<=100;
proc sort; by newid; run;
/* 600100+640900+664600=1905600 observations and 5 variables */

proc contents data=bsdata; run;
proc freq data=bsdata; tables outcome replicate status; run;
proc contents data=clinical.clinicalhrsdeident_gdr_20210624; run;

*Merge with original dataset to get the covariates;
data bsdata2; 
  merge clinical.clinicalhrsdeident_gdr_20210624 bsdata (in=A);
  by newid;
  if A;
  if outcome='adl' then normwgt=normwgt_adl;
  else if outcome='walk' then normwgt=normwgt_walk;
  else if outcome='death' then normwgt=normwgt_death;
proc sort; by replicate newid; run;
/* 1905600 observations and 153+4+1=158 variables */

proc delete data=bsdata; run; quit;

%let DATA=bsdata3;
proc sql noprint; select max(replicate) format 3. into :B from bsdata2; quit; /*create macro variable with total number of bootstrap datasets*/
%put "&B"; 

%macro stats_bs;
 %do i=1 %to &B; /*do this for each bootstrap dataset up to the last bootstrap dataset*/
     proc sql noprint; select (&i-1)*(nobs/&B)+1 into :fobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the first id of lth bs sample*/
     proc sql noprint; select &i*(nobs/&B) into :lobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the last id of lth bs sample*/

	 *Select bs data;
     data bsdata3;
	 	set bsdata2 (FIRSTOBS=&fobs OBS=&lobs);
	 proc sort; by outcome; run;

	 *Define 'VARNAME' and 'VARCLASS' of baTCIC model for corresponding bs data;
	 data _null_;
	   set baTCICbs (keep=replicate VARINMODEL2 VARCLASS);
	   where replicate=&i;
	   call symputx ('VARNAME', VARINMODEL2);
	   call symputx ('VARCLASS', VARCLASS);
	 run;

     proc phreg data = &DATA ;
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &VARCLASS;
       weight normwgt; /*use normalized weight*/
       model time*status(0) = &VARNAME;
	   output out=BSOUT xbeta=xb;
      run;
      proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=RECURSIVE iauc); 
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &VARCLASS;
       model time*status(0) = &VARNAME / nofit;
	   roc 'CoxWeighted' pred=xb;
	   ods output CONCORDANCE=concord IAUC=iaucr;
      run;
	  proc phreg data = BSOUT rocoptions(method=IPCW iauc); 
	   by outcome; /*fit 3 regression model, one/outcome*/
	   class &VARCLASS;
	   model time*status(0) = &VARNAME / nofit;
	   roc 'CoxWeighted' pred=xb;
	   ods output IAUC=iauc;
	  run;

	  proc delete data=BSOUT bsdata3; run; quit;

	  data estimates;
	 	set concord(keep=outcome estimate rename=(estimate=c))
	    iaucr(keep=outcome estimate rename=(estimate=iaucr))
	    iauc(keep=outcome estimate rename=(estimate=iauc));
	  run;

      data estimates2 (keep=replicate cbs_adl cbs_walk cbs_death cbs_avg iaucrbs_adl iaucrbs_walk iaucrbs_death iaucrbs_avg iaucbs_adl iaucbs_walk iaucbs_death iaucbs_avg);
	    set estimates end=last;
	  	retain replicate cbs_adl cbs_walk cbs_death cbs_avg iaucrbs_adl iaucrbs_walk iaucrbs_death iaucrbs_avg iaucbs_adl iaucbs_walk iaucbs_death iaucbs_avg;
		replicate=&i;
		if outcome='adl' and c ne . then cbs_adl=c;
		else if outcome='adl' and iaucr ne . then iaucrbs_adl=iaucr;
		else if outcome='adl' and iauc ne . then iaucbs_adl=iauc;

		else if outcome='walk' and c ne . then cbs_walk=c;
		else if outcome='walk' and iaucr ne . then iaucrbs_walk=iaucr;
		else if outcome='walk' and iauc ne . then iaucbs_walk=iauc;

		else if outcome='death' and c ne . then cbs_death=c;
		else if outcome='death' and iaucr ne . then iaucrbs_death=iaucr;
		else if outcome='death' and iauc ne . then iaucbs_death=iauc;

		if last;
		if cbs_adl ne . and cbs_walk ne . and cbs_death ne . then cbs_avg=mean(cbs_adl, cbs_walk, cbs_death);
	    if iaucrbs_adl ne . and iaucrbs_walk ne . and iaucrbs_death ne . then iaucrbs_avg=mean(iaucrbs_adl, iaucrbs_walk, iaucrbs_death);
	    if iaucbs_adl ne . and iaucbs_walk ne . and iaucbs_death ne . then iaucbs_avg=mean(iaucbs_adl, iaucbs_walk, iaucbs_death);
     run;

     proc append base=stats_bs data=estimates2 force; run;
     proc delete data=concord iaucr iauc estimates estimates2; run; quit;
%end; /*i loop*/

%mend stats_bs;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%stats_bs
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
/*======MONITORING: 2021-07-13, 10:26======*/
/*======MONITORING: 2021-07-13, 10:42======*/


/*** Step 3(change as indicated above): Compute the C-statistic/IAUC of each bootstrap model evaluated in the original case-study data, namely C-statistic/IAUC-original ***/

*redefine VARINMODEL: 'ageint5' replace with 'age3_lin age3_sp1 age3_sp2', bmi4g replace with bmi3_lin bmi3_sp1 bmi3_sp2;
data baTCICbs; 
 set outdata.baTCIC_rep100;
 numVarsfin=countw(VARINMODEL, '');
 VARINMODEL2=VARINMODEL;
 VARINMODEL2=tranwrd(VARINMODEL2, "ageint5", "age3_lin age3_sp1 age3_sp2");
 VARINMODEL2=tranwrd(VARINMODEL2, "bmi4g", "bmi3_lin bmi3_sp1 bmi3_sp2");
 /*define VARCLASS variables*/
 VARCLASS=VARINMODEL;
 VARCLASS=tranwrd(VARCLASS, "ageint5", "");
 VARCLASS=tranwrd(VARCLASS, "bmi4g", "");
run;

/*QC*/
proc contents data=baTCICbs; run;
proc means data=baTCICbs n mean std stderr clm median p25 p75 maxdec=4; var numVarsfin; run;
proc print data=baTCICbs (obs=10); var replicate VARINMODEL VARINMODEL2 VARCLASS; run;

*Prepare BS data;
data bsdata;
 set savedata.clinical_bs500_gdr_20210707;
 where 0<replicate<=100;
proc sort; by newid; run;
/* 600100+640900+664600=1905600 observations and 5 variables */

proc contents data=bsdata; run;
proc freq data=bsdata; tables outcome replicate status; run;
proc contents data=clinical.clinicalhrsdeident_gdr_20210624; run;

*Merge with original dataset to get the covariates;
data bsdata2; 
  merge clinical.clinicalhrsdeident_gdr_20210624 bsdata (in=A);
  by newid;
  if A;
  if outcome='adl' then normwgt=normwgt_adl;
  else if outcome='walk' then normwgt=normwgt_walk;
  else if outcome='death' then normwgt=normwgt_death;
proc sort; by replicate newid; run;
/* 1905600 observations and 153+4+1=158 variables */

proc delete data=bsdata; run; quit;


*Prepare original data;
data finaldata;
 set clinical.clinicalhrsdeident_gdr_20210624;
 if status_adldepdth2=2 then do; status_adldepdth2=0; time_adldepdth2=19.1130249; end;
 if status_walkdepdth2=2 then do; status_walkdepdth2=0; time_walkdepdth2=19.1130249; end;
 rename subgroup_adldep=subgroup_adl;
 rename subgroup_walkdep=subgroup_walk;
 subgroup_death=1;
run;
/*6646 observations and 154 variables.*/

*Separate 3 outcomes in 3 datasets;
data adl (rename=(status_adldepdth2=status time_adldepdth2=time normwgt_adl=normwgt)) ;
 set finaldata(drop=status_walkdepdth2 death time_walkdepdth2 time2death normwgt_walk normwgt_death);
 length outcome $16;
 outcome="adl";
 where subgroup_adl=1;
run; /* 6001 observations and 149 variables.*/
data walk (rename=(status_walkdepdth2=status time_walkdepdth2=time normwgt_walk=normwgt)) ;
 set finaldata(drop=status_adldepdth2 death time_adldepdth2 time2death normwgt_adl normwgt_death);
 length outcome $16;
 outcome="walk";
 where subgroup_walk=1;
run; /*6409 observations and 149 variables.*/
data death (rename=(death=status time2death=time normwgt_death=normwgt)) ;
 set finaldata(drop=status_adldepdth2 status_walkdepdth2 time_adldepdth2 time_walkdepdth2 normwgt_adl normwgt_walk );
 length outcome $16;
 outcome="death";
 where subgroup_death=1;
run; /*6646 observations and 149 variables.*/

*Stack 3 datasets;
data all;
 set adl walk death;
proc sort; by outcome newid; run; 
/*sort by outcome because this is important to use in BY statement of the PHREG procedure*/

/*6001+6409+6646=19056 and 149 variables*/

proc delete data=finaldata adl walk death; run; quit;

%let DATAorig=all;
%let DATAboot=bsdata3;
proc sql noprint; select max(replicate) format 3. into :B from bsdata2; quit; /*create macro variable with total number of bootstrap datasets*/
%put "&B"; 


%macro stats_bs_ori;
 %do i=1 %to &B; /*do this for each bootstrap baTCIC model*/

     proc sql noprint; select (&i-1)*(nobs/&B)+1 into :fobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the first id of lth bs sample*/
     proc sql noprint; select &i*(nobs/&B) into :lobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the last id of lth bs sample*/

	 *Select bs data;
     data bsdata3;
	 	set bsdata2 (FIRSTOBS=&fobs OBS=&lobs);
	 proc sort; by outcome; run;

	 *Define 'VARNAME' and 'VARCLASS' of baBIC model for corresponding bs_ori data;
	 data _null_;
	   set baTCICbs (keep=replicate VARINMODEL2 VARCLASS);
	   where replicate=&i;
	   call symputx ('VARNAME', VARINMODEL2);
	   call symputx ('VARCLASS', VARCLASS);
	 run;

	 /*Get bootstrap model fitted in bootstrap data  */
     proc phreg data = &DATAboot;
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &VARCLASS;
       weight normwgt; /*use normalized weight*/
       model time*status(0) = &VARNAME;
	   store bootmodel; /*requests that the fitted model be saved to an item store  */
      run;

     /*Get linear predictions (predicted) using fitted model above in original data  */
	 proc plm restore=bootmodel;
	   score data=&DATAorig out=BSOUT predicted; /* score (Linear predictor) new observations based on the item store bootmodel that was created above*/
	 run;

	  /*Using linear predictions "predicted" above compute the C-statistic/IAUC-original*/
      proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=RECURSIVE iauc); 
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &VARCLASS;
       model time*status(0) = &VARNAME / nofit;
	   roc 'CoxWeighted' pred=predicted;
	   ods output CONCORDANCE=concord IAUC=iaucr;
      run;
	  proc phreg data = BSOUT rocoptions(method=IPCW iauc); 
	   by outcome; /*fit 3 regression model, one/outcome*/
	   class &VARCLASS;
	   model time*status(0) = &VARNAME / nofit;
	   roc 'CoxWeighted' pred=predicted;
	   ods output IAUC=iauc;
	  run;
	  proc delete data=BSOUT bsdata3; run; quit;

	  data estimates;
	 	set concord(keep=outcome estimate rename=(estimate=c))
	    iaucr(keep=outcome estimate rename=(estimate=iaucr))
	    iauc(keep=outcome estimate rename=(estimate=iauc));
	  run;

      data estimates2 (keep=replicate cbs_ori_adl cbs_ori_walk cbs_ori_death cbs_ori_avg iaucrbs_ori_adl iaucrbs_ori_walk iaucrbs_ori_death iaucrbs_ori_avg iaucbs_ori_adl iaucbs_ori_walk iaucbs_ori_death iaucbs_ori_avg);
	    set estimates end=last;
	  	retain replicate cbs_ori_adl cbs_ori_walk cbs_ori_death cbs_ori_avg iaucrbs_ori_adl iaucrbs_ori_walk iaucrbs_ori_death iaucrbs_ori_avg iaucbs_ori_adl iaucbs_ori_walk iaucbs_ori_death iaucbs_ori_avg;
		replicate=&i;
		if outcome='adl' and c ne . then cbs_ori_adl=c;
		else if outcome='adl' and iaucr ne . then iaucrbs_ori_adl=iaucr;
		else if outcome='adl' and iauc ne . then iaucbs_ori_adl=iauc;

		else if outcome='walk' and c ne . then cbs_ori_walk=c;
		else if outcome='walk' and iaucr ne . then iaucrbs_ori_walk=iaucr;
		else if outcome='walk' and iauc ne . then iaucbs_ori_walk=iauc;

		else if outcome='death' and c ne . then cbs_ori_death=c;
		else if outcome='death' and iaucr ne . then iaucrbs_ori_death=iaucr;
		else if outcome='death' and iauc ne . then iaucbs_ori_death=iauc;

		if last;
		if cbs_ori_adl ne . and cbs_ori_walk ne . and cbs_ori_death ne . then cbs_ori_avg=mean(cbs_ori_adl, cbs_ori_walk, cbs_ori_death);
	    if iaucrbs_ori_adl ne . and iaucrbs_ori_walk ne . and iaucrbs_ori_death ne . then iaucrbs_ori_avg=mean(iaucrbs_ori_adl, iaucrbs_ori_walk, iaucrbs_ori_death);
	    if iaucbs_ori_adl ne . and iaucbs_ori_walk ne . and iaucbs_ori_death ne . then iaucbs_ori_avg=mean(iaucbs_ori_adl, iaucbs_ori_walk, iaucbs_ori_death);
     run;

     proc append base=stats_bs_ori data=estimates2 force; run;
     proc delete data=concord iaucr iauc estimates estimates2; run; quit;
%end; /*i loop*/

%mend stats_bs_ori;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%stats_bs_ori 
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

/*
======MONITORING: 2022-06-13, 10:35======
======MONITORING: 2022-06-13, 11:16======
*/
options source notes; 
ods select all; 

data outdata2.stats_bs_ori; set stats_bs_ori; proc sort;  by replicate; run;

*From previous version;
data outdata2.stats_bs;
 set outdata.stats_optimis (keep=replicate cbs_adl cbs_walk cbs_death cbs_avg iaucrbs_adl iaucrbs_walk iaucrbs_death iaucrbs_avg iaucbs_adl iaucbs_walk iaucbs_death iaucbs_avg) ;
proc sort;  by replicate; run;

*Merge stats_bs with stats_bs_ori;
data outdata2.stats_optimis;
 merge outdata2.stats_bs outdata2.stats_bs_ori;
 by replicate;
run;
/* 100 observations and 25 variables*/


/*** Step 4: Calculate the optimism in the fit of each bootstrap sample as: abs(C-statistic/IAUC-boot - C-statistic/IAUC-original) ***/
data op_baTCIC (keep=replicate optimism_c_adl optimism_c_walk optimism_c_death optimism_c_avg optimism_iaucr_adl optimism_iaucr_walk optimism_iaucr_death optimism_iaucr_avg
                     optimism_iauc_adl optimism_iauc_walk optimism_iauc_death optimism_iauc_avg);
 set outdata2.stats_optimis;
 optimism_c_adl=abs(cbs_adl-cbs_ori_adl);
 optimism_c_walk=abs(cbs_walk-cbs_ori_walk);
 optimism_c_death=abs(cbs_death-cbs_ori_death);
 optimism_c_avg=abs(cbs_avg-cbs_ori_avg);

 optimism_iaucr_adl=abs(iaucrbs_adl-iaucrbs_ori_adl);
 optimism_iaucr_walk=abs(iaucrbs_walk-iaucrbs_ori_walk);
 optimism_iaucr_death=abs(iaucrbs_death-iaucrbs_ori_death);
 optimism_iaucr_avg=abs(iaucrbs_avg-iaucrbs_ori_avg);

 optimism_iauc_adl=abs(iaucbs_adl-iaucbs_ori_adl);
 optimism_iauc_walk=abs(iaucbs_walk-iaucbs_ori_walk);
 optimism_iauc_death=abs(iaucbs_death-iaucbs_ori_death);
 optimism_iauc_avg=abs(iaucbs_avg-iaucbs_ori_avg);
run;


*Create permanent data set with baTCIC optimism;
data outdata2.op_baTCIC; set op_baTCIC; run;

*QC;
proc print data=outdata2.stats_optimis (obs=3); var replicate cbs_adl cbs_walk cbs_death cbs_avg; run;
proc print data=outdata2.stats_optimis (obs=3); var replicate cbs_ori_adl cbs_ori_walk cbs_ori_death cbs_ori_avg; run;
proc print data=op_baTCIC (obs=3); run;
proc means data=op_baTCIC; run;


/*** Step 5: Average the optimism across 100 bootstrap sample, namely O ***/
proc means data=outdata2.op_baTCIC stackodsoutput n mean std stderr clm maxdec=4; 
 var optimism_c_adl optimism_c_walk optimism_c_death optimism_c_avg optimism_iaucr_adl optimism_iaucr_walk optimism_iaucr_death optimism_iaucr_avg
     optimism_iauc_adl optimism_iauc_walk optimism_iauc_death optimism_iauc_avg;
 ods output summary=outdata2.opAvg_baTCIC (keep=Variable Mean);
run;

/*Note:
Include the average optimism for each outcome.

To compute the average Optimism for the 3 outcomes:
1-) Obtain average C-stat-BS across outcomes for each BS (C-stat-BS-avg)
2-) Obtain average C-stat-original across outcomes for each BS (C-stat-original-avg)
3-) For each BS: Compute Absolute difference: C-stat-BS-avg - C-stat-Original-avg
4-) Compute Average (Absolute difference: C-stat-BS-avg - C-stat-original-avg) across 500 BS

The corrected C-stat final models = C-stat of original sample (without Wolbers approximation) – degree of optimism (using Wolbers approximation). 
*/


/*** Step 6: Compute the optimism-corrected C-statistic/IAUC of the case-study data as: C-statistic/IAUC-apparent - O ***/

data stats_correct (rename=(Mean=Optimism));
  if _N_=1 then set outdata.stats_app; 
  set outdata2.opAvg_baTCIC;
run;

data outdata2.stats_correct (keep=Variable Optimism c_correct iaucr_correct iauc_correct);
 set stats_correct;
/* retain c_correct_adl c_correct_walk c_correct_death c_correct_avg iaucr_correct_adl iaucr_correct_walk iaucr_correct_death iaucr_correct_avg;*/
 if Variable="optimism_c_adl" then c_correct=capp_adl-Optimism;
 else if Variable="optimism_c_walk" then c_correct=capp_walk-Optimism;
 else if Variable="optimism_c_death" then c_correct=capp_death-Optimism;
 else if Variable="optimism_c_avg" then c_correct=capp_avg-Optimism;
 else if Variable="optimism_iaucr_adl" then iaucr_correct=iaucrapp_adl-Optimism;
 else if Variable="optimism_iaucr_walk" then iaucr_correct=iaucrapp_walk-Optimism;
 else if Variable="optimism_iaucr_death" then iaucr_correct=iaucrapp_death-Optimism;
 else if Variable="optimism_iaucr_avg" then iaucr_correct=iaucrapp_avg-Optimism;
 else if Variable="optimism_iauc_adl" then iauc_correct=iaucapp_adl-Optimism;
 else if Variable="optimism_iauc_walk" then iauc_correct=iaucapp_walk-Optimism;
 else if Variable="optimism_iauc_death" then iauc_correct=iaucapp_death-Optimism;
 else if Variable="optimism_iauc_avg" then iauc_correct=iaucapp_avg-Optimism;
run;

proc print data=outdata2.stats_correct; run;

PROC EXPORT DATA= outdata2.stats_correct
            OUTFILE= "V:\Health and Retirement Study\Grisell\AlexSei\R01eprognosisClinicalPaper\sasdata\BackwardSelectionHRS\bootstrapping\baTCIC\20220613_NewVersion\stats_correct.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


/****************************************************************************************************************************************************************************************/
*Compute Location-shifted bootstrap confidence interval;

/*Note:
We computed the location-shifted bootstrap confidence intervals of the optimism-corrected C-statistic/IAUC by subtracting the optimism estimate from
the 2.5th and 97.5th percentiles of the C-statistic/IAUC-bootstrap distribution
*/


/*
Reference: as of 2.26.2021
Noma H, Shinozaki T, Iba K, Teramukai S, Furukawa TA. Confidence intervals of prediction accuracy measures for multivariable prediction models based
on the bootstrap-based optimism correction methods. arXiv preprint arXiv:2005.01457. https://arxiv.org/ftp/arxiv/papers/2005/2005.01457.pdf
*/

/*Method description:
Algorithm 1 (Location-shifted bootstrap confidence interval)
1. For a multivariable prediction model, let theta_hat_app be the apparent predictive measure for the derivation population and
   let theta_hat be the optimism-corrected predictive measure obtained from the Harrell’s bias correction, 0.632, or 0.632+ method.
2. In the computational processes of theta_hat, we can obtain a bootstrap estimate of the sampling distribution of theta_hat_app from the B bootstrap samples.
   Compute the bootstrap confidence interval of theta_app from the bootstrap distribution, (theta_hat_app_L, theta_hat_app_U); 
   for the 95% confidence interval, they are typically calculated by the 2.5th and 97.5th percentiles of the bootstrap distribution.
3. Calculate the bias estimate by optimism, delta_hat=theta_hat_app-theta_hat
4. Then, the location-shifted bootstrap confidence interval is computed as (theta_hat_app_L- theta_hat, theta_hat_app_U- theta_hat)
*/


/**********************************
2. In the computational processes of theta_hat, we can obtain a bootstrap estimate of the sampling distribution of hat_app from the B bootstrap samples.
   Compute the bootstrap confidence interval of theta_app from the bootstrap distribution, (theta_hat_app_L, theta_hat_app_U); 
   for the 95% confidence interval, they are typically calculated by the 2.5th and 97.5th percentiles of the bootstrap distribution.
*/

*Custom percentiles baTCIC;
proc stdize data=outdata2.stats_optimis PctlMtd=ORD_STAT outstat=baTCIC pctlpts=2.5, 97.5;
 var cbs_adl cbs_walk cbs_death iaucrbs_adl iaucrbs_walk iaucrbs_death iaucbs_adl iaucbs_walk iaucbs_death;
run;
 
data baTCIC;
 set baTCIC;
 where _type_ =: 'P';
run;

proc print data=baTCIC noobs; run;

proc transpose data=baTCIC out=wide ;
   by _type_;
   var cbs_adl cbs_walk cbs_death iaucrbs_adl iaucrbs_walk iaucrbs_death iaucbs_adl iaucbs_walk iaucbs_death;
run;

proc print data=wide noobs; run;

data wide (drop=_NAME_ _type_);
 set wide (rename=(col1=P2_5));
 length Variable $20;
 if _NAME_="cbs_adl" then Variable="optimism_c_adl";
 else if _NAME_="cbs_walk" then Variable="optimism_c_walk";
 else if _NAME_="cbs_death" then Variable="optimism_c_death";
 else if _NAME_="iaucrbs_adl" then Variable="optimism_iaucr_adl";
 else if _NAME_="iaucrbs_walk" then Variable="optimism_iaucr_walk";
 else if _NAME_="iaucrbs_death" then Variable="optimism_iaucr_death";
 else if _NAME_="iaucbs_adl" then Variable="optimism_iauc_adl";
 else if _NAME_="iaucbs_walk" then Variable="optimism_iauc_walk";
 else if _NAME_="iaucbs_death" then Variable="optimism_iauc_death";

 if _type_="P97_5000" and Variable="optimism_c_adl" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_c_walk" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_c_death" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_iaucr_adl" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_iaucr_walk" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_iaucr_death" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_iauc_adl" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_iauc_walk" then P97_5=P2_5;
 else if _type_="P97_5000" and Variable="optimism_iauc_death" then P97_5=P2_5;
run;

proc print data=wide noobs; run;

data wide2_5 (drop=p97_5);
 set wide ;
 where P97_5=.;
proc sort; by Variable; run;
data wide97_5 (drop=p2_5);
 set wide ;
 where P97_5 ne .;
proc sort; by Variable; run;
data baTCIC2;
 retain Variable P2_5 P97_5;
 merge wide2_5 wide97_5;
 by Variable;
run;
proc print data=baTCIC2; run;
proc contents data=baTCIC2; run;

proc delete data=wide wide2_5 wide97_5; run;

proc print data=outdata2.stats_correct; run;
proc contents data=outdata2.stats_correct; run;
proc sort data=outdata2.stats_correct out=stats_correct; by Variable; run;

data baTCIC3;
 merge stats_correct baTCIC2;
 by Variable;
 P2_5correct=P2_5-Optimism; /*Optimism: optimism estimate*/
 P97_5correct=P97_5-Optimism;
 if Variable in ("optimism_c_avg", "optimism_iaucr_avg", "optimism_iauc_avg") then delete;
run;
proc print data=baTCIC3; run;
proc print data=baTCIC; run;

*Save permanent data;
data outdata2.stats_correct_percentiles;
 set baTCIC3;
run;

PROC EXPORT DATA= outdata2.stats_correct_percentiles
            OUTFILE= "path\stats_correct_percentiles.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
