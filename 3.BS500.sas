/*********** Purpose: Create 500 bootstrap datasets ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2020.05.19 ***********/
/*********** Date completed: 2021.07.07 ***********/

options MERGENOBY=ERROR;
libname harmo "path";
libname rand "path";
proc format cntlin=rand.sasfmts; run;
proc format cntlin=harmo.formats; run;
libname clinical "path";
libname out "path"; 
options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors to be written to the log*/
options errorabend; /*so the process stops as soon as an error is encountered*/


/********************************************************************************************************************************************************************************/
/*** Create input data for each outcome ***/

proc contents data=clinical.clinicalhrsdeident_gdr_20210624; run;

*Death;
data input_death (keep=newid outcome RAESTRAT RAEHSAMP time2death death rename=(time2death=time death=status)); 
 set clinical.clinicalhrsdeident_gdr_20210624;
 length outcome $10;
 outcome='death';
proc sort; by RAESTRAT RAEHSAMP; run;
/*  6646 observations and 6 variables  */
proc freq data=input_death; tables status; run;

*ADL;
data input_adl (keep=newid outcome RAESTRAT RAEHSAMP time_adldepdth2 status_adldepdth2 rename=(time_adldepdth2=time status_adldepdth2=status));
 set clinical.clinicalhrsdeident_gdr_20210624 ;
 length outcome $10;
 outcome='adl';
 where subgroup_adldep=1;
 if status_adldepdth2=2 then do; status_adldepdth2=0; time_adldepdth2=19.1130249; end; *Apply Wolbers approximation to the Competing-risk setting where death is switched to censored and time-to-event=longest possible time-to-event;
proc sort; by RAESTRAT RAEHSAMP; run;
/*  6001 observations and 6 variables  */
proc freq data=input_adl; tables status; run;

*WALK;
data input_walk(keep=newid outcome RAESTRAT RAEHSAMP time_walkdepdth2 status_walkdepdth2 rename=(time_walkdepdth2=time status_walkdepdth2=status));
 set clinical.clinicalhrsdeident_gdr_20210624;
 length outcome $10;
 outcome='walk';
 where subgroup_walkdep=1;
 if status_walkdepdth2=2 then do; status_walkdepdth2=0; time_walkdepdth2=19.1130249; end; *Apply Wolbers approximation to the Competing-risk setting where death is switched to censored and time-to-event=longest possible time-to-event;
proc sort; by RAESTRAT RAEHSAMP; run;
/* 6409 observations and 6 variables */
proc freq data=input_walk; tables status; run;

/********************************************************************************************************************************************************************************/
/*** Create bootstrap samples using PROC SURVEYSELECT for each outcome***/

******************** Death;
*Create sample size input data set to provide the stratum sample sizes in the bootstrapping; 
proc freq data=input_death noprint;
  tables RAESTRAT*RAEHSAMP/out=nsize_death(rename=(count=_nsize_));
run;

*Run SURVEYSELECT to generate data with replicated Ids;
proc surveyselect data=input_death out=bsample_death method=urs sampsize=nsize_death reps=500 outhits seed=20210707;
 strata RAESTRAT RAEHSAMP;
run;
/*Total Sample size is : 6646*500=3,323,000 observations and 10 variables.  */

/*
Notes: 
sampsize=nsize: so that the bootstrap samples have the same number of Ids within each combination of RAESTRAT*RAEHSAMP
outhits:
	includes a distinct copy of each selected unit in the OUT= output data set when the same sampling unit is selected more than once.
	By default, the output data set contains a single copy of each unit selected, even when a unit is selected more than once, 
	and the variable NumberHits records the number of hits (selections) for each unit.
	If you specify the OUTHITS option, the output data set contains m copies of a sampling unit for which NumberHits is m;
	for example, the output data set contains three copies of a unit that is selected three times (NumberHits is 3).
*/


******************** ADL;
*Create sample size input data set to provide the stratum sample sizes in the bootstrapping; 
proc freq data=input_adl noprint;
  tables RAESTRAT*RAEHSAMP/out=nsize_adl(rename=(count=_nsize_));
run;

*Run SURVEYSELECT to generate data with replicated Ids;
proc surveyselect data=input_adl out=bsample_adl method=urs sampsize=nsize_adl reps=500 outhits seed=20210707;
 strata RAESTRAT RAEHSAMP;
run;
/*Total Sample size is : 6001*500=3,000,500 observations and 10 variables.  */


******************** WALK;
*Create sample size input data set to provide the stratum sample sizes in the bootstrapping; 
proc freq data=input_walk noprint;
  tables RAESTRAT*RAEHSAMP/out=nsize_walk(rename=(count=_nsize_));
run;

*Run SURVEYSELECT to generate data with replicated Ids;
proc surveyselect data=input_walk out=bsample_walk method=urs sampsize=nsize_walk reps=500 outhits seed=20210707;
 strata RAESTRAT RAEHSAMP;
run;
/*Total Sample size is : 6409*500=3,204,500 observations and 10 variables.  */


/********************************************************************************************************************************************************************************/
 *Create permanent dataset with BS samples of all outcomes;

data out.clinical_bs500_gdr_20210707;
 set bsample_adl (keep=replicate newid status time outcome) bsample_walk (keep=replicate newid status time outcome) bsample_death (keep=replicate newid status time outcome);
proc sort; by newid; run;
/*
NOTE: There were 3000500 observations read from the data set WORK.BSAMPLE_ADL.
NOTE: There were 3204500 observations read from the data set WORK.BSAMPLE_WALK.
NOTE: There were 3323000 observations read from the data set WORK.BSAMPLE_DEATH.
NOTE: The data set OUT.CLINICAL_BS500_GDR_20210707 has 9528000 observations and 5 variables.
*/


proc contents data=out.clinical_bs500_gdr_20210707; run;
