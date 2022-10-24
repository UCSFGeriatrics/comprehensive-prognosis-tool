/*********** Purpose: Apply Unavailable predictor method to final baTCIC model:deleting any 2 predictors ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2020.05.19 ***********/
/*********** Date completed: 2021.09.14 ***********/

libname savedata "path";
libname outdata "path";
libname outdata2 "path";
libname harmo 'path';
libname rand 'path';
proc format cntlin=rand.sasfmts; run;
proc format cntlin=harmo.formats; run;

options nothreads; /*to avoid Insufficient space error in Utility files*/
options nomlogic nomprint nomrecall nosymbolgen;
options errorabend; /*so the process stops as soon as an error is encountered and you don't end up wasting additional system resources */
options noquotelenmax; /*supress warning: THE QUOTED STRING CURRENTLY BEING PROCESSED HAS BECOME MORE THAN 262 CHARACTERS LONG*/
ods select none; /*to create output data sets through the ODS OUTPUT statement and suppress the display of all output*/

/*Check system options specified at SAS invocation*/
proc options option=work; run;
%put %sysfunc(pathname(work));
proc options option=utilloc; run;
proc options group=(memory performance); run;

options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors to be written to the log*/

*****************************************************PREPARE DATA FOR MACROs ********************************************************;

data finaldata; set savedata.clinicalhrsdeident_gdr_20210624; run; /*6646 observations and 153 variables.*/
data finaldata;
 set finaldata;
 rename subgroup_adldep=subgroup_adl;
 rename subgroup_walkdep=subgroup_walk;
 subgroup_death=1;
run;
/*6646 observations and 153+1=154 variables.*/

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


*Stack 2 datasets for competing-risk;
data compdata;
 set adl walk;
proc sort; by outcome newid; run; 
/*sort by outcome because this is important to use in BY statement of the PHREG procedure*/

/*6001+6646=12410 observations and 149 variables.*/

proc delete data=finaldata adl walk; run; quit;

*Create data with all models that have 1 predictor less than final one-minute model;
data baTCIC (keep=nummod VARINMODEL VARSEC);
 set outdata2.batcic_1unavailpred;
 retain nummod;
 nummod+1;
 where delevar ne "";
run;


***************************************************** DEFINE ARGUMENTS FOR MACROS ********************************************************;

*Note: all the macro variables for outcomes need to be in the same order, e.g. : adl walk death; 
%let DATACOMP=compdata; /*data for Competing-risk regression*/
%let DATACOX=death;  /*data for Cox regression*/
%let NUMOUTCOMES=3; /*number of outcomes*/
%let ALLLABELS= adl walk death; /*labels for outcomes*/
%let NUMPRED=13; /*number of predictors WITHOUT including predictors that are forced in (e.g. ageint5, gender*/
%let NUMPREDPLUSONE=14; /*number of predictors defined in &NUMPRED (e.g.16) plus one*/

%let Cvars=C_adl C_walk C_death; /*name of variables for the C-stat*/
%let IAUCvars=iauc_adl iauc_walk iauc_death; /*name of variables for the IAUC stat*/
%let AICvars=AIC_adl AIC_walk AIC_death; /*name of variables for the AIC stat*/
%let TCICvars=TCIC_adl TCIC_walk TCIC_death; /*name of variables for the TCIC stat*/
%let TCICindvars=TCIC_adl TCIC_walk TCIC_death; /*name of variables for the TCIC in the individual outcome selection*/
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
proc sql noprint; select nobs format 2. into :nmodels from dictionary.tables where libname='WORK' and memname='BATCIC'; quit; 
%put "&nmodels"; 

***************************************************** BACKWARD ELIMINATION 3 OUTCOMES **************************************************************;

%macro models;
 %do l=1 %to &nmodels; /*do this for each of the models that have 1 predictor less than the final model up to the last model*/
  data _null_;
   set baTCIC;
   where nummod=&l;
   call symputx ('INIVARNAME' ,VARINMODEL);
   call symputx ('INIVARSEC' ,VARSEC);
  run;
 
  %let nummod=&l;
  %let INIVARCLASS=%sysfunc(compbl(%sysfunc(tranwrd(&INIVARNAME,ageint5,%str( ))))); /*select the final set of categorical variables to run model by replacing ageint5 with blank*/

  %PUT "&l" ======MORNITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
  %best_TCIC;
  %PUT "&l" ======MORNITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;

  data TCIC_&l; set TCIC; nummod=&l; run;
  proc delete data=TCIC; run; quit;

  %if &l=1 %then %do; /*in the 1st model create TCICall dataset*/
   proc append base=TCICall data=TCIC_&l force; run; /*BICsim will have the information from all best individual models for each outcome and for each simulated dataset*/
   sasfile WORK.TCICall load;
  %end;
  %else %do; /*for the rest of the models: keep updating TCICall dataset in memory*/
   proc append base=TCICall data=TCIC_&l; run;
  %end; 
  proc delete data=TCIC_&l; run; quit;

 %end; /*l loop*/

%mend models;

%macro best_TCIC;
 
	*Get the TCIC=SBC of full models with 40 predictors;
	%let BASEFULL=%sysfunc(compbl(ageint5 gender R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g smoking3g stroke3g volunteer R5PROXY)); /*all 40 initial variables*/
    %let CLASSFULL=%sysfunc(compbl(gender R5BATHA R5BEDA R5DRESSA R5EATA R5TOILTA alcohol3g R5ARTHRE R5CANCRE R5DIABE
 								R5FALL R5HEARING heart3g R5HIBPE R5MAPA R5MEALSA R5MEDSA R5MONEYA R5PHONEA R5SHOPA R5URINA H5LVALONE
 								R5LUNGE mstat R5ARMSA R5CHAIRA otherclim3g R5DIMEA R5LIFTA R5PUSHA R5SITA R5STOOPA
 								R5WALKSA pain3g bmi4g smoking3g stroke3g volunteer R5PROXY)); /*all 40-ageint5(continuous variable)=39 initial variables*/

	/*Competing-risk regression models*/
	proc phreg data = &DATACOMP ;
	  by outcome; /*fit 2 regression models, eg.: adl, walk, one/outcome*/
      class &CLASSFULL;
      weight normwgt; /*use normalized weight*/
      model time*status(0) = &BASEFULL /eventcode=1;
	  ods output FITSTATISTICS=FITS1;
    run;

	/*Cox regression model*/
     proc phreg data = &DATACOX ;
	   by outcome; /*fit 1 regression model*/
       class &CLASSFULL;
       weight normwgt; /*use normalized weight*/
       model time*status(0) = &BASEFULL;
	   ods output FITSTATISTICS=FITS1_cox ;
      run;

	 data fits;
 	  set FITS1 FITS1_cox;
     run;

	*Create variables with BIC/TCIC of full model;
	data _null_;
	   set fits;
	   if outcome='adl' and criterion='SBC' then call symputx ('full_adl'  ,WithCovariates);
	   else if outcome='walk' and criterion='SBC' then call symputx ('full_walk'  ,WithCovariates);
	   else if outcome='death' and criterion='SBC' then call symputx ('full_death'  ,WithCovariates);
	run;

	*Create variables with TCIC of best individual models;
	data _null_;
	   set outdata.TCIC_unionCompRiskNormWgt (keep= &TCICindvars);
	   call symputx (' best_adl'  ,TCIC_adl);
	   call symputx (' best_walk' ,TCIC_walk);
	   call symputx (' best_death' ,TCIC_death);
	run;

	proc delete data=FITS1 FITS1_cox fits; run; quit;

	/*Start with 15 predictors*/
	%let BASE=%sysfunc(compbl(&INIVARNAME)); /*all 15 initial variables*/
    %let CLASS=%sysfunc(compbl(&INIVARCLASS)); /*all 15-ageint5(continuous variable)=14 initial variables*/
    %let DELE=%sysfunc(compbl(%sysfunc(tranwrd(&CLASS,gender,%str( ))))); /*all 15 variables-2(ageint5 gender)=13 initial variables*/
	%let BASESEC=&INIVARSEC; /*&INIVARSEC has the initial total time cost for the 15-pred model. As we go from one model to a smaller model the BASESEC variable will be replaced for a new basesec based on the variables left in the model*/

	%DeleteOneVar;
     data _null_;
       set AVGCTABLE (keep=VARINMODEL VARINCLASS DELELIST VARSEC);
       call symputx ('BASE' ,VARINMODEL);
       call symputx ('CLASS',VARINCLASS);
	   call symputx ('DELE' ,DELELIST);
	   call symputx ('BASESEC', VARSEC);
     run;
     
	proc append base=TCIC data=AVGCTABLE (drop=VARINCLASS DELELIST) force; run;
    proc delete data=AVGCTABLE; run; quit;

%mend best_TCIC;

%macro DeleteOneVar;
    %do j=1 %to %eval(&NUMPREDPLUSONE-1);
        %let DELEVAR=%scan(&DELE,&j); /*select the jth word to delete. &DELE is defined in 'best_TCIC' macro*/
		%let VARNAME=%sysfunc(compbl(%sysfunc(tranwrd(&BASE,&DELEVAR,%str( ))))); /*select the final set of variables to run model by replacing the deleted variable with blank*/
		%let VARCLASS=%sysfunc(compbl(%sysfunc(tranwrd(&CLASS,&DELEVAR,%str( ))))); /*select the final set of categorical variables to run model by replacing the deleted variable with blank*/
		%let VARSEC=%sysevalf(&BASESEC-&&&DELEVAR.SEC); /*&&&DELEVAR.SEC: this is to tell sas that the variable delevar.sec keep changing as we delete a different variable each time*/

	/*Competing-risk regression models*/
     proc phreg data = &DATACOMP ;
	  by outcome; /*fit 2 regression models, eg.: adl, walk, one/outcome*/
      class &VARCLASS;
      weight normwgt; /*use normalized weight*/
      model time*status(0) = &VARNAME /eventcode=1;
	  ods output FITSTATISTICS=FITS1 FailureSummary=CensUncens; /*CensoredSummary if Cox, FailureSummary if Competing-risks*/
	  output out=BSOUT xbeta=xb;
    run;
	data BSOUT; 
 	 set BSOUT; 
 	 if status=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
  	  status=0;
  	  time=19.1130249;
     end;
    run;
	proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=ipcw iauc); 
	  by outcome; /*fit 2 regression models, one/outcome*/
      class &VARCLASS;
      model time*status(0) = &VARNAME / nofit;
	  roc 'Normwgt' pred=xb;
	  ods output CONCORDANCE=concord IAUC=iauc;
    run;
	proc delete data=BSOUT; run; quit; 

	/*Cox regression model*/
     proc phreg data = &DATACOX ;
	   by outcome; /*fit 1 regression model*/
       class &VARCLASS;
       weight normwgt; /*use normalized weight*/
       model time*status(0) = &VARNAME;
	   output out=BSOUT xbeta=xb;
	   ods output FITSTATISTICS=FITS1_cox CensoredSummary=CensUncens_cox; /*CensoredSummary if Cox, FailureSummary if Competing-risks*/
      run;
      proc phreg data = BSOUT CONCORDANCE=HARRELL rocoptions(method=ipcw iauc); 
	   by outcome; /*fit 1 regression model*/
       class &VARCLASS;
       model time*status(0) = &VARNAME / nofit;
	   roc 'Normwgt' pred=xb;
	   ods output CONCORDANCE=concord_cox IAUC=iauc_cox;
      run;
	  proc delete data=BSOUT; run; quit;

		data estimates;
	 	 set FITS1 FITS1_cox
	         concord(keep=outcome estimate rename=(estimate=c)) 	concord_cox(keep=outcome estimate rename=(estimate=c))
	         iauc(keep=outcome estimate rename=(estimate=iauc)) 	iauc_cox(keep=outcome estimate rename=(estimate=iauc))
	         CensUncens (keep=outcome event rename=(event=Uncens))	CensUncens_cox (keep=outcome event rename=(event=Uncens));
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

	    %if &j=1 %then %do; /*in the 1st line of AVGCTABLE create AVGCTABLE dataset*/
	     proc append base=AVGCTABLE data=FITS2 force; run;
	     sasfile WORK.AVGCTABLE load;
	    %end;
	    %else %do; /*for the rest of lines of AVGCTABLE: keep updating AVGCTABLE dataset in memory*/
	     proc append base=AVGCTABLE data=FITS2; run;
	    %end; 
	   proc delete data=FITS2 FITS1 FITS1_Cox CensUncens CensUncens_Cox concord concord_cox iauc iauc_cox estimates; run; quit;
	  %end; /*j loop*/
	sasfile WORK.AVGCTABLE close;

%mend DeleteOneVar;

%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;
%models;
%PUT ======MONITORING: %SYSFUNC(DATE(),YYMMDD10.), %LEFT(%SYSFUNC(TIME(),HHMM8.))======;


***************************************************** END BACKWARD 3 OUTCOMES ******************************************************************;

*Create permanent dataset;
sasfile WORK.TCICall close;
data outdata2.baTCIC_UnavailPred&NUMPREDPLUSONE.nCr2; set TCICall; run;
proc delete data=TCICall; run; quit;

*Leave only unique models. This is a combination: n!/(r!(n-r)!)=14!/(2!(14-2)!)=91 different models;

**1-) Create new VARINMODEL2 variable with sorted predictors within each string;
data baTCIC_2UnavailPred;
 set outdata2.baTCIC_UnavailPred&NUMPREDPLUSONE.nCr2;
 array predict[14] $32 _temporary_; /*create temporary array with 14 variables: predict1-predict14*/
 call missing(of predict[*]); /*Assigns missing values to the specified character or numeric variables.*/
 do i = 1 to dim(predict) until(p eq 0);
     call scan(VARINMODEL,i,p,l); /*CALL SCAN(<string>, count, position, length): Returns the position and length of the nth word from a character string.*/
     predict[i] = substrn(VARINMODEL,p,l); /*Returns a substring: divide the string into words, use position(p) and length (l) given by function call scan */
 end;
 call sortc(of predict[*]); /*Sorts the values of character arguments in an array.*/
 length VARINMODEL2 $600;
 VARINMODEL2 = catx(' ',of predict[*]); /*Removes leading and trailing blanks, inserts delimiters, and returns a concatenated character string.*/
 drop i p l;
proc sort; by VARINMODEL2 nummod; run;
/* 14 models *13 predictors=182 observations and 20 variables.*/

*QC;
proc print data=baTCIC_2UnavailPred; var nummod VARINMODEL2 VARINMODEL; run;
*14-pred_baTCIC model in baTCIC backward selection;
proc print data=baTCIC_2UnavailPred; var nummod VARINMODEL2 VARINMODEL DELEVAR;
where VARINMODEL="ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g "; run;


**2-) Leave only unique VARINMODEL2/models;
data outdata2.baTCIC_2UnavailPred (drop=VARINMODEL2);
 set baTCIC_2UnavailPred;
 by VARINMODEL2;
 if first.VARINMODEL2;
proc sort; by nummod; run;
/*91 observations and 19 variables*/

ods select all; /*to print results below*/
ods listing close; /*turn of the output window / "listing" output, so I don't get WARNING: Data too long for column*/
ods csv file="path\Results_baTCIC_2UnavailPred.csv";
proc print data=outdata2.baTCIC_2UnavailPred; run;
ods csv close;

proc means data=outdata2.baTCIC_2UnavailPred stackodsoutput n mean clm stderr min max maxdec=4;
 var C_adl C_walk C_death C_avg iauc_adl iauc_walk iauc_death iauc_avg;
 ods output summary=outdata2.baTCIC_2UnavailPredSummary;
run;

ods csv file="path\Results_baTCIC_2UnavailPredSummary.csv";
proc print data=outdata2.baTCIC_2UnavailPredSummary; run;
ods csv close;

*QC;
*All of them have 16-2=14 predictors as expected;
data test;
 set outdata2.baTCIC_2UnavailPred;
 delims = ' ';   
 numPred = countw(VARINMODEL, delims);
run;
proc means data=test; var numPred; run;





