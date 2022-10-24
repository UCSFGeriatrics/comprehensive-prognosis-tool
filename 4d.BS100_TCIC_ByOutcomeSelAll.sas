***********************************************************************************************************************************************************************************;
*Purpose: Create data with best models for individual selection by outcomes using bs dataset datasets and compute summary statistic                                                ;                                     
*Statistician: Grisell Diaz-Ramirez																																				   ;
*Started: 2020.02.05																																							   ;
*Finished: 2021.07.09                                                                                                                                                              ;
***********************************************************************************************************************************************************************************;


/*Check system options specified at SAS invocation*/
proc options option=work; run;
proc options option=utilloc; run;
proc options group=(memory performance); run;

libname outdata "path";

%let B=100;

/*Merge all 4 datasets with the best individual model within each outcome for each bs dataset */
data outdata.tcic_rep&B._all (drop=delims);
   merge outdata.TCIC_rep&B._adl outdata.TCIC_rep&B._walk outdata.TCIC_rep&B._death;
   by replicate;
   delims = ' ';
   numVarsfinadl=countw(VARINMODEL_adl, delims);
   numVarsfinwalk=countw(VARINMODEL_walk, delims);
   numVarsfindeath=countw(VARINMODEL_death, delims);
run; /* 100 observations and 22 variables.*/

*Summary stats of bs dataset models by outcome;

proc means data=outdata.tcic_rep&B._all stackodsoutput n mean std stderr clm median p25 p75 maxdec=4; 
 var numVarsfinadl numVarsfinwalk numVarsfindeath C_adl C_walk C_death iauc_adl iauc_walk iauc_death VARSEC_adl VARSEC_walk VARSEC_death ;
 ods output summary=outdata.tcic_rep&B._all_gralstats(rename=(Mean=Mean_bs_ByOutcome));
proc sort; by variable; run;

PROC EXPORT DATA= outdata.tcic_rep&B._all_gralstats
            OUTFILE="V:\Health and Retirement Study\Grisell\AlexSei\R01eprognosisClinicalPaper\sasdata\BackwardSelectionHRS\bootstrapping\Individual\tcic_rep&B._all_gralstats.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
