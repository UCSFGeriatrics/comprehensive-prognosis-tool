***********************************************************************************************************************************************************************************;
*Purpose: Perform TCIC backward elimination with time information and without using Competing-risk for BS datasets by outcome: ADL                                                 ;                                     
*Statistician: Grisell Diaz-Ramirez																																				   ;
*Started: 2020.02.05																																							   ;
*Finished: 2021.07.07                                                                                                                                                              ;
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
libname clinical "path";
libname savedata "path";
options nomlogic nomprint nomrecall nosymbolgen;
options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors to be written to the log*/
options noquotelenmax; /*supress warning: THE QUOTED STRING CURRENTLY BEING PROCESSED HAS BECOME MORE THAN 262 CHARACTERS LONG*/
ods select none; /*to create output data sets through the ODS OUTPUT statement and suppress the display of all output*/
options nothreads; /*to avoid Insufficient space error in Utility files*/
options errorabend; /*so the process stops as soon as an error is encountered*/


*****************************************************PREPARE DATA FOR MACROs ********************************************************;

*Only use 100 BS samples and outcome='adl';
data bsdata;
 set savedata.clinical_bs500_gdr_20210707;
 where 0<replicate<=100 and outcome="adl";
proc sort; by newid; run;
/* 600100 observations and 5 variables */
proc contents data=bsdata; run;
proc freq data=bsdata; tables outcome replicate status; run;
proc contents data=clinical.clinicalhrsdeident_gdr_20210624; run;

*Merge with original dataset to get the covariates;
data bsdata2; 
  merge clinical.clinicalhrsdeident_gdr_20210624 bsdata (in=A);
  by newid;
  if A;
proc sort; by replicate newid; run;
/* 600100 observations and 153+4=157 variables */

proc delete data=bsdata; run; quit;

**************************************************** DEFINE ARGUMENTS FOR MACROS: outcome, best_tcic, DeleteOneVar, cstat ********************************************************;

%let DATA=bsdata3;
%let NORMWGT= normwgt_adl; /*normalized weight variables for each outcome*/
%let LABEL= adl ; /*labels for outcomes*/
%let NUMPRED=38; /*number of predictors WITHOUT including predictors that are forced in (e.g. ageint5, gender*/
%let NUMPREDPLUSONE=39; /*number of predictors defined in &NUMPRED (e.g.38) plus one*/
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

***************************************************** BACKWARD ELIMINATION INDIVIDUAL OUTCOME **************************************************************;
/*Macro outcome, for each bootstrap dataset and each outcome: 
1-) Calls macro best_tcic (calls macro DeleteOneVar)
2-) Creates one dataset with results from backward elimination
3-) Create a dataset with one line corresponding to the best model for each individual outcome and each replicate
*/

%macro outcome;
 %do i=1 %to &B; /*do this for each bootstrap dataset up to the last bootstrap dataset*/
  proc sql noprint; select (&i-1)*(nobs/&B)+1 into :fobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the first id of ith replicate*/
  proc sql noprint; select &i*(nobs/&B) into :lobs from dictionary.tables where libname='WORK' and memname='BSDATA2'; quit; /*create macro variable with the last id of ith replicate*/

  data bsdata3;
   set bsdata2 (FIRSTOBS=&fobs OBS=&lobs);
  run;

  *sasfile TEMP.bsdata3 load;

  %PUT "&i" ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
  
	%best_tcic;

    *sasfile TEMP.TCIC_&label close;

	/*Create new TCIC_&label dataset to merge with other outcomes*/
	data TCIC_&label (rename=(VARINMODEL=VARINMODEL_&label VARSEC=VARSEC_&label)); 
     set TCIC_&label (drop=DELELIST DELEVAR);
     replicate=&i;
    run;

	/*For each outcome and each replicate: select the best individual model*/
	proc sort data=TCIC_&label; by descending TCIC_&label; 
    data TCIC_&label;
     set TCIC_&label point=nobs nobs=nobs; 
	 output;
	 stop;
    run;


  *sasfile TEMP.bsdata3 close;
  proc delete data=bsdata3; run; quit;

  proc append base=TCICreplicate_&label data=TCIC_&label force; run; /*TCICreplicate will have the information from all best individual models for each outcome and for each bootstrap dataset*/
  proc delete data=TCIC_&label; run; quit;
  %PUT "&i" ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

 %end; /*i loop*/
%mend outcome;

%macro best_tcic;
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
	/*264.01 is the initial total time cost when all 40 variables are in the model. As we go from one model to a smaller model the basesec variable will be replaced for a new basesec based on the variables left in the model*/


    	proc phreg data = &DATA ; 
      	 class &CLASS;
         weight &NORMWGT; /*use normalized weight*/
      	 model time*status(0) = &BASE;
		 output out=BSOUT xbeta=xb;
	  	 ods output FITSTATISTICS=FITS1 CensoredSummary=CensUncens; /*CensoredSummary if Cox, FailureSummary if Competing-risks*/
        run;
        proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=RECURSIVE iauc); 
          class &CLASS;
          model time*status(0) = &BASE / nofit;
		  roc 'NoCompRiskC' pred=xb;
		  ods output CONCORDANCE=concord IAUC=iauc;
        run;
	   proc delete data=BSOUT; run; quit;

	*Need this dataset for selection with time cost:;
	data _null_; set CensUncens; call symputx ('Uncens' ,Event); run;

	data _null_; set concord; call symputx ('c' ,Estimate); run;
	data _null_; set iauc; call symputx ('iauc' ,Estimate); run;

	data base_&label (keep=VARINMODEL DELEVAR DELELIST VARSEC AIC_&label TCIC_&label C_&label iauc_&label); 
      length VARINMODEL $ 600 DELELIST $600 DELEVAR $ 60;
      set FITS1 end=last;
      retain VARINMODEL DELEVAR VARSEC C_&label iauc_&label AIC_&label TCIC_&label LOGL_&label;
      format AIC_&label TCIC_&label LOGL_&label 10.4; 
	  VARINMODEL="&BASE";
	  DELEVAR=" ";
	  DELELIST=" "; 
	  VARSEC=&BASESEC;
	  TOTALSEC=&TOTALSEC;
	  Uncens=&Uncens;
	  DF=&TOTALDF;
	  C_&label=&c;
	  iauc_&label=&iauc;
	  if CRITERION='AIC' then AIC_&label=WITHCOVARIATES;
	  if CRITERION='-2 LOG L' then LOGL_&label=WITHCOVARIATES;
      TCIC_&label=LOGL_&label+(DF*(VARSEC/TOTALSEC))*log(Uncens);  
      if last;
     run;
	
    proc delete data=FITS1 CensUncens concord iauc; run; quit;

   proc append base=TCIC_&label data=base_&label force; run; /*TCIC is the final dataset with backward elimination steps from 38 variables to 2*/
   *sasfile TEMP.TCIC_&label load;
   proc delete data=base_&label; run; quit;

	%do j=1 %to &NUMPRED;
	  %DeleteOneVar;
	     data _null_;
	       set CTABLE_&label (keep=VARINMODEL VARINCLASS DELELIST VARSEC);
	       call symputx ('BASE' ,VARINMODEL);
	       call symputx ('CLASS',VARINCLASS);
		   call symputx ('DELE' ,DELELIST);
		   call symputx ('BASESEC', VARSEC);
	     run;

      proc append base=TCIC_&label data=CTABLE_&label(drop=VARINCLASS) force; run;
	  proc delete data=CTABLE_&label; run; quit;
	%end;
%mend best_tcic;

%macro DeleteOneVar;
 	%do k=1 %to %eval(&NUMPREDPLUSONE-&j);
        %let DELEVAR=%scan(&DELE,&k); /*select the jth word to delete. &DELE is defined in 'best_tcic' macro*/
		%let VARNAME=%sysfunc(compbl(%sysfunc(tranwrd(&BASE,&DELEVAR,%str( ))))); /*select the final set of variables to run model by replacing the deleted variable with blank*/
		%let VARCLASS=%sysfunc(compbl(%sysfunc(tranwrd(&CLASS,&DELEVAR,%str( ))))); /*select the final set of categorical variables to run model by replacing the deleted variable with blank*/
		%let VARSEC=%sysevalf(&BASESEC-&&&DELEVAR.SEC); /*&&&DELEVAR.SEC: this is to tell sas that the variable delevar.sec keep changing as we delete a different variable each time*/

    	proc phreg data = &DATA ; 
      	 class &VARCLASS;
         weight &NORMWGT; /*use normalized weight*/
      	 model time*status(0) = &VARNAME;
		 output out=BSOUT xbeta=xb;
	  	 ods output FITSTATISTICS=FITS1 CensoredSummary=CensUncens; /*CensoredSummary if Cox, FailureSummary if Competing-risks*/
        run;
        proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=RECURSIVE iauc); 
          class &VARCLASS;
          model time*status(0) = &VARNAME / nofit;
		  roc 'CompRiskC' pred=xb;
		  ods output CONCORDANCE=concord IAUC=iauc;
        run;
	   proc delete data=BSOUT; run; quit;

	*Need this dataset for selection with time cost:;
	data _null_; set CensUncens; call symputx ('Uncens' ,Event); run;

	data _null_; set concord; call symputx ('c' ,Estimate); run;
	data _null_; set iauc; call symputx ('iauc' ,Estimate); run;

	data FITS2 (keep=VARINMODEL VARINCLASS DELEVAR DELELIST VARSEC  AIC_&label TCIC_&label C_&label iauc_&label); 
      length VARINMODEL $ 600 DELELIST $600 VARINCLASS $600 DELEVAR $ 60;
      set FITS1 end=last;
      retain VARINMODEL VARINCLASS DELEVAR VARSEC C_&label iauc_&label AIC_&label TCIC_&label LOGL_&label; 
      format AIC_&label TCIC_&label LOGL_&label 10.4; 
	  VARINMODEL="&VARNAME"; /*variables in the model*/
	  VARINCLASS="&VARCLASS"; /*variables in class statement*/
	  DELEVAR="&DELEVAR"; /*deleted variable*/
	  DELELIST=compbl(tranwrd("&DELE","&DELEVAR",' ')); 
	  /*'DELELIST' contain the list of variables that we need to start with in subsequent run. That is, the 2nd time we run macro 'DeleteOneVar' we have 35 variables in macro variable 'DELE',
	  at this step we call it 'DELELIST' and it contains the 'DELE' list minos the variable deleted in the previous run. Later, 'DELELIST' is redefined as 'DELE'*/
	  VARSEC=&VARSEC; /*time cost of reduced model*/
	  TOTALSEC=&TOTALSEC;
	  Uncens=&Uncens;
	  DF=&TOTALDF;
	  C_&label=&c;
	  iauc_&label=&iauc;
	  if CRITERION='AIC' then AIC_&label=WITHCOVARIATES;
	  if CRITERION='-2 LOG L' then LOGL_&label=WITHCOVARIATES;
      TCIC_&label=LOGL_&label+(DF*(VARSEC/TOTALSEC))*log(Uncens);  
      if last;
     run;

    /*%if &k=1 %then %do;*/ /*in the 1st line of CTABLE_&label create CTABLE_&label dataset*/
     proc append base=CTABLE_&label data=FITS2 force; run;
     /*sasfile CTABLE_&label load;*/
    /*%end;*/
   /* %else %do;*/ /*for the rest of lines of CTABLE_&label: keep updating CTABLE_&label dataset in memory*/
    /* proc append base=CTABLE_&label data=FITS2 force; run;*/
    /*%end; */
    proc delete data=FITS1 FITS2 CensUncens concord iauc; run; quit;


  %end;
  *sasfile TEMP.CTABLE_&label close;

  proc sort data=CTABLE_&label; by descending TCIC_&label; run;
  data CTABLE_&label;
   set CTABLE_&label point=nobs nobs=nobs;
   output;
   stop;
  run; /*choose the last observation with minimum TCIC*/
%mend DeleteOneVar;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%outcome;
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

data savedata.TCIC_rep&B._&label; set TCICreplicate_&label; proc sort; by replicate; run;


***************************************************** END BACKWARD ELIMINATION INDIVIDUAL OUTCOME ******************************************************************;

