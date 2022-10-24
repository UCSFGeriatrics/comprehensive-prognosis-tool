***********************************************************************************************************************************************************************************;
*Purpose: Perform baTCIC backward elimination with time information and without using Competing-risk for BS datasets                                                               ;                                     
*Statistician: Grisell Diaz-Ramirez																																				   ;
*Started: 2020.02.05																																							   ;
*Finished: 2021.07.09                                                                                                                                                              ;
***********************************************************************************************************************************************************************************;

/*Check system options specified at SAS invocation*/
proc options option=work; run;
%put %sysfunc(pathname(work));
proc options option=utilloc; run;
proc options group=(memory performance); run;

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


*****************************************************PREPARE DATA FOR MACROs ********************************************************;
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


***************************************************** DEFINE ARGUMENTS FOR MACROS ********************************************************;

*Note: all the macro variables for outcomes need to be in the same order, e.g. : adl walk death; 
%let DATA=bsdata3;
%let NUMOUTCOMES=3; /*number of outcomes*/
%let ALLLABELS= adl walk death; /*labels for outcomes*/
%let NUMPRED=38; /*number of predictors WITHOUT including predictors that are forced in (e.g. ageint5, gender*/
%let NUMPREDPLUSONE=39; /*number of predictors defined in &NUMPRED (e.g.38) plus one*/

%let Cvars=C_adl C_walk C_death; /*name of variables for the C-stat*/
%let IAUCvars=iauc_adl iauc_walk iauc_death; /*name of variables for the IAUC stat*/
%let AICvars=AIC_adl AIC_walk AIC_death; /*name of variables for the AIC stat*/
%let TCICvars=TCIC_adl TCIC_walk TCIC_death; /*name of variables for the TCIC stat*/
%let TCICindvars=TCIC_adl TCIC_walk TCIC_death; /*name of variables for the TCIC in the individual outcome selection. This is to use with data: tcic_rep100_all*/
%let LOGLvars=LOGL_adl LOGL_walk LOGL_death; /*name of variables for the LOGL stat*/
%let UNCENSvars=UNCENS_adl UNCENCS_walk UNCENS_death; /*name of variables for the UNCENS stat*/
%let TCICfullvars=TCICfull_adl TCICfull_walk TCICfull_death; /*name of variables for the TCIC full stat*/
%let TCICbestvars=TCICbest_adl TCICbest_walk TCICbest_death; /*name of variables for the TCIC best stat*/

*Times that takes to collect each of the predictors in the model;
%let R5BATHASEC=5.45; %let R5BEDASEC=6.02;  ;%let R5DRESSASEC=7.21; %let R5EATASEC=6.12; %let R5TOILTASEC=7.23;
%let alcohol3gsec=10.76; %let R5ARTHRESEC=6.06; %let R5CANCRESEC=7.13; %let ageint5sec=3.36; %let R5DIABESEC=5.07;
%let R5FALLSEC=6.23; %let R5HEARINGSEC=4.44; %let heart3gsec=14.19; %let R5HIBPESEC=5.1;  %let R5MAPASEC=7.55;  %let R5MEDSASEC=5.46;  %let R5MONEYASEC=9.15;
%let R5MEALSASEC=6.07; %let R5PHONEASEC=6.15; %let R5SHOPASEC=6.15; %let R5URINASEC=8.33; %let H5LVALONESEC=1;
%let R5LUNGESEC=6.4; %let mstatsec=7.04; %let R5ARMSASEC=6.48; %let R5CHAIRASEC=7.08; %let otherclim3gsec=12.69;
%let R5DIMEASEC=5.55; %let R5LIFTASEC=8.1; %let R5PUSHASEC=7.33; %let R5SITASEC=7.08; %let R5STOOPASEC=6.13;
%let R5WALKSASEC=5.18; %let pain3gsec=7.35; %let bmi4gsec=4.25; %let gendersec=3.41; %let smoking3gsec=8.49;
%let stroke3gsec=6.95; %let volunteersec=9.27; %let R5PROXYSEC=1;
%let TOTALSEC=264.01; /*total time cost when all 40 variables are in the model*/
%let TOTALDF=51; /*DF when all 40 variables are in the model*/

proc sql noprint; select max(replicate) format 3. into :B from bsdata2; quit; /*create macro variable with total number of bootstrap datasets*/
%put "&B"; 


***************************************************** BACKWARD ELIMINATION 3 OUTCOMES **************************************************************;
/*Macro replicate, for each bootstrap dataset:
1-) Calls macro best_TCIC (calls macro DeleteOneVar)
2-) Creates one dataset with results from backward elimination
3-) Selects the 'best model' and stack results in one common dataset
*/

%macro replicate;
 %do l=1 %to &B; /*do this for each bootstrap dataset up to the last bootstrap dataset*/
     proc sql noprint; select (&l-1)*(nobs/&B)+1 into :fobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the first id of lth bs sample*/
     proc sql noprint; select &l*(nobs/&B) into :lobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the last id of lth bs sample*/

     data bsdata3;
	 	set bsdata2 (FIRSTOBS=&fobs OBS=&lobs);
	 proc sort; by outcome; run;

	 *sasfile TEMP.bsdata3 load;

	%PUT "&l" ======MORNITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
	%best_TCIC;
	%PUT "&l" ======MORNITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

	*Keep baTCIC model for lth replicate;
	data TCIC_&l; set TCIC (drop=DELELIST); replicate=&l; run;
	*sasfile TEMP.TCIC close;
    *sasfile TEMP.bsdata3 close; 
	proc delete data=bsdata3; run; quit;

	proc sort data=TCIC_&l; by descending TCIC_avg; run;
	data TCIC_&l;
	 set TCIC_&l point=nobs nobs=nobs;
	 output;
	 stop;
	run; 

    proc append base=TCICbs data=TCIC_&l force; run; /*TCICbs will have the information from all baBIC and onemin models for each bootstrap dataset*/
    *sasfile TEMP.TCICbs load;

	proc delete data=TCIC_&l; run; quit;

 %end; /*l loop*/

%mend replicate;


%macro best_TCIC;
	%let BASE=%sysfunc(compbl(ageint5 gender R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g smoking3g stroke3g volunteer R5PROXY)); /*all 40 initial variables*/
    %let CLASS=%sysfunc(compbl(gender R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g smoking3g stroke3g volunteer R5PROXY)); /*all 40-ageint5(continuous variable)=39 initial variables*/
    %let DELE=%sysfunc(compbl(R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g smoking3g stroke3g volunteer R5PROXY)); /*all 40 variables-2(ageint5 gender)=38*/
	%let BASESEC=264.01; 
	/*264.01 is the initial total time cost when all 38 variables are in the model. As we go from one model to a smaller model the basesec variable will be replaced for a new basesec based on the variables left in the model*/

	/*Run PHREG for the 1st time with all the predictors in the model*/

     proc phreg data = &DATA;
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &CLASS;
       weight normwgt; /*use normalized weight*/
       model time*status(0) = &BASE;
	   output out=BSOUT xbeta=xb;
	   ods output FITSTATISTICS=FITS1 CensoredSummary=CensUncens; /*CensoredSummary if Cox, FailureSummary if Competing-risks*/
      run;
      proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=RECURSIVE iauc); 
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &CLASS;
       model time*status(0) = &BASE / nofit;
	   roc 'NoCompRiskC' pred=xb;
	   ods output CONCORDANCE=concord IAUC=iauc;
      run;
	  proc delete data=BSOUT; run; quit;

	data estimates;
 	 set FITS1
         concord(keep=outcome estimate rename=(estimate=c))
         iauc(keep=outcome estimate rename=(estimate=iauc))
         CensUncens (keep=outcome event rename=(event=Uncens));
	run;

	data _null_;
	   set outdata.tcic_rep&B._all (keep=replicate &TCICindvars);
	   where replicate=&l;
	   call symputx (' best_adl'  ,TCIC_adl);
	   call symputx (' best_walk' ,TCIC_walk);
	   call symputx (' best_death' ,TCIC_death);
	run;

	data BASE_TCIC (keep= VARINMODEL DELEVAR DELELIST VARSEC &Cvars &IAUCvars &AICvars &TCICvars &TCICfullvars C_avg iauc_avg TCIC_avg); 
	  length VARINMODEL $ 600 DELELIST $600 DELEVAR $ 60;
	  set estimates end=last;
	  retain VARINMODEL DELEVAR VARSEC TOTALSEC DF &Cvars &IAUCvars &AICvars &TCICvars &LOGLvars &UNCENSvars &TCICfullvars &TCICbestvars;
	  format &Cvars &IAUCvars &AICvars &TCICvars &LOGLvars &TCICfullvars &TCICbestvars 10.4;
	  VARINMODEL="&BASE";
	  DELEVAR=" ";
	  DELELIST=" "; 
	  VARSEC=&BASESEC;
	  TOTALSEC=&TOTALSEC;
	  DF=&TOTALDF;

	  array label{&NUMOUTCOMES} &ALLLABELS;
	  array Cs{&NUMOUTCOMES} &Cvars;
	  array IAUCs{&NUMOUTCOMES} &IAUCvars;
	  array AICs{&NUMOUTCOMES} &AICvars;
	  array TCICs{&NUMOUTCOMES} &TCICvars;
	  array LOGLs{&NUMOUTCOMES} &LOGLvars;
	  array UNCENSs {&NUMOUTCOMES} &UNCENSvars;
	  array TCICfull{&NUMOUTCOMES} &TCICfullvars;
	  array TCICbest{&NUMOUTCOMES} &TCICbestvars;
	  array best{&NUMOUTCOMES} (&best_adl &best_walk &best_death);

	  do k=1 to &NUMOUTCOMES;
	   TCICbest{k}=best{k};
	   if outcome=vname(label{k}) then do;
	    if c ne . then Cs{k}=c;
	    else if iauc ne . then IAUCs{k}=iauc;
		else if Uncens ne . then UNCENSs{k}=Uncens;
		else if CRITERION='AIC' then AICs{k}=WITHCOVARIATES;
        else if CRITERION='-2 LOG L' then LOGLs{k}=WITHCOVARIATES;
		else if CRITERION='SBC' then do;
		 TCICfull{k}=WITHCOVARIATES; /*with time cost, the TCIC_full_model = BIC_full_model because time for k (fitted model, i.e. full model) = time for M (full model*/
		end;
		TCICs{k}=abs((LOGLs{k}+(DF*(VARSEC/TOTALSEC))*log(UNCENSS{k}))-TCICbest{k})/abs(TCICfull{k}-TCICbest{k});
	   end;
	  end; /*k loop*/
	 if last;
     C_avg=mean(of &Cvars);
	 iauc_avg=mean(of &IAUCvars);
     TCIC_avg=mean(of &TCICvars);
    run;

	proc delete data=FITS1 CensUncens concord iauc estimates; run; quit;

	data _null_;
	   set BASE_TCIC (keep= &TCICfullvars);
	   call symputx ('full_adl'  ,TCICfull_adl);
	   call symputx ('full_walk'  ,TCICfull_walk);
	   call symputx ('full_death'  ,TCICfull_death);
	run;

   proc append base=TCIC data=BASE_TCIC(drop=&TCICfullvars) force; run; /*TCIC is the final dataset with backward elimination steps from 38 variables to 2*/
   *sasfile TEMP.TCIC load;
   proc delete data=BASE_TCIC; run; quit;

	%do i=1 %to &NUMPRED;
	  %DeleteOneVar;
	     data _null_;
	       set AVGCTABLE2 (keep=VARINMODEL VARINCLASS DELELIST VARSEC);
	       call symputx ('BASE' ,VARINMODEL);
	       call symputx ('CLASS',VARINCLASS);
		   call symputx ('DELE' ,DELELIST);
		   call symputx ('BASESEC', VARSEC);
	     run;

        proc append base=TCIC data=AVGCTABLE2 (drop=VARINCLASS); run;
	    proc delete data=AVGCTABLE2; run; quit;
	%end; /*i loop*/
%mend best_TCIC;

%macro DeleteOneVar;
    %do j=1 %to %eval(&NUMPREDPLUSONE-&i);
        %let DELEVAR=%scan(&DELE,&j); /*select the jth word to delete. &DELE is defined in 'best_TCIC' macro*/
		%let VARNAME=%sysfunc(compbl(%sysfunc(tranwrd(&BASE,&DELEVAR,%str( ))))); /*select the final set of variables to run model by replacing the deleted variable with blank*/
		%let VARCLASS=%sysfunc(compbl(%sysfunc(tranwrd(&CLASS,&DELEVAR,%str( ))))); /*select the final set of categorical variables to run model by replacing the deleted variable with blank*/
		%let VARSEC=%sysevalf(&BASESEC-&&&DELEVAR.SEC); /*&&&DELEVAR.SEC: this is to tell sas that the variable delevar.sec keep changing as we delete a different variable each time*/

     proc phreg data = &DATA ;
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &VARCLASS;
       weight normwgt; /*use normalized weight*/
       model time*status(0) = &VARNAME;
	   output out=BSOUT xbeta=xb;
	   ods output FITSTATISTICS=FITS1 CensoredSummary=CensUncens; /*CensoredSummary if Cox, FailureSummary if Competing-risks*/
      run;
      proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=RECURSIVE iauc); 
	   by outcome; /*fit 3 regression model, one/outcome*/
       class &VARCLASS;
       model time*status(0) = &VARNAME / nofit;
	   roc 'NoCompRiskC' pred=xb;
	   ods output CONCORDANCE=concord IAUC=iauc;
      run;
	  proc delete data=BSOUT; run; quit;

		data estimates;
	 	 set FITS1
	         concord(keep=outcome estimate rename=(estimate=c))
	         iauc(keep=outcome estimate rename=(estimate=iauc))
	         CensUncens (keep=outcome event rename=(event=Uncens));
		run;

		data FITS2 (keep= VARINMODEL VARINCLASS DELEVAR DELELIST VARSEC &Cvars &IAUCvars &AICvars &TCICvars C_avg iauc_avg TCIC_avg); 
		  length VARINMODEL $ 600 VARINCLASS $600 DELELIST $600 DELEVAR $ 60;
		  set estimates end=last;
		  retain VARINMODEL DELEVAR VARSEC &Cvars &IAUCvars &AICvars &TCICvars &LOGLvars &UNCENSvars &TCICfullvars &TCICbestvars;
		  format &Cvars &IAUCvars &AICvars &TCICvars &LOGLvars &TCICfullvars &TCICbestvars 10.4;
		  VARINMODEL="&VARNAME"; /*variables in the model*/
		  VARINCLASS="&VARCLASS"; /*variables in class statement*/
		  DELEVAR="&DELEVAR"; /*deleted variable*/
		  DELELIST=compbl(tranwrd("&DELE","&DELEVAR",' ')); 
		  /*'DELELIST' contain the list of variables that we need to start with in subsequent run. That is, the 2nd time we run macro 'DeleteOneVar' we have 36 variables in macro variable 'DELE',
		  at this step we call it 'DELELIST' and it contains the 'DELE' list minos the variable deleted in the previous run. Later, 'DELELIST' is redefined as 'DELE'*/
	  	  VARSEC=&VARSEC; /*time cost of reduced model*/
		  TOTALSEC=&TOTALSEC;
		  DF=&TOTALDF;

		  array label{&NUMOUTCOMES} &ALLLABELS;
		  array Cs{&NUMOUTCOMES} &Cvars;
		  array IAUCs{&NUMOUTCOMES} &IAUCvars;
		  array AICs{&NUMOUTCOMES} &AICvars;
		  array TCICs{&NUMOUTCOMES} &TCICvars;
		  array LOGLs{&NUMOUTCOMES} &LOGLvars;
		  array UNCENSS {&NUMOUTCOMES} &UNCENSvars;
		  array TCICfull{&NUMOUTCOMES} &TCICfullvars;
		  array full{&NUMOUTCOMES} (&full_adl &full_walk &full_death);
		  array TCICbest{&NUMOUTCOMES} &TCICbestvars;
		  array best{&NUMOUTCOMES} (&best_adl &best_walk &best_death);

		  do k=1 to &NUMOUTCOMES;
		   TCICbest{k}=best{k};
		   TCICfull{k}=full{k};
		   if outcome=vname(label{k}) then do;
		    if c ne . then Cs{k}=c;
		    else if iauc ne . then IAUCs{k}=iauc;
			else if Uncens ne . then UNCENSS{k}=Uncens;
			else if CRITERION='AIC' then AICs{k}=WITHCOVARIATES;
	        else if CRITERION='-2 LOG L' then LOGLs{k}=WITHCOVARIATES;
		    TCICs{k}=abs((LOGLs{k}+(DF*(VARSEC/TOTALSEC))*log(UNCENSS{k}))-TCICbest{k})/abs(TCICfull{k}-TCICbest{k});
           end;
		  end; /*k loop*/
		 if last;
	     C_avg=mean(of &Cvars);
	     iauc_avg=mean(of &IAUCvars);
         TCIC_avg=mean(of &TCICvars);     
		run;

/*	    %if &j=1 %then %do;*/ /*in the 1st line of AVGCTABLE create AVGCTABLE dataset*/
	     proc append base=AVGCTABLE data=FITS2 force; run;
/*	     sasfile AVGCTABLE load;*/
/*	    %end;*/
/*	    %else %do;*/ /*for the rest of lines of AVGCTABLE: keep updating AVGCTABLE dataset in memory*/
/*	     proc append base=AVGCTABLE data=FITS2 force; run;*/
/*	    %end; */
	    proc delete data=FITS1 concord iauc CensUncens estimates FITS2; run; quit;
	  %end; /*j loop*/
	  *sasfile TEMP.AVGCTABLE close;

   proc sort data= AVGCTABLE; by descending TCIC_avg; run;

   data AVGCTABLE2;
    set AVGCTABLE point=nobs nobs=nobs;
    output;
    stop;
   run; /*choose the last observation that has minimum TCIC*/

   proc delete data= AVGCTABLE; run; quit;
%mend DeleteOneVar;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%replicate;
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

***************************************************** END BACKWARD 3 OUTCOMES ******************************************************************;

/*Create permanent dataset with the best (i.e. minimum TCIC) and 1 min combined model for each bs sample*/
data outdata2.baTCIC_rep&B; set TCICbs; run;
