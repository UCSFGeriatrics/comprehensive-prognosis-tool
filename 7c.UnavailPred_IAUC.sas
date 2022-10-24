/*********** Purpose: Compute IAUC using IPCW method for best/worst case models missing 1 or 2 predictors   ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2021.12.07 ***********/
/*********** Date completed: 2021.12.07  ***********/

libname harmo "path";
libname rand "path";
proc format cntlin=rand.sasfmts; run;
proc format cntlin=harmo.formats; run;
libname clinical "path";
libname out "path"; 


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

/****************************************************************************************************************************************************************************************/
/*ADL*/

**********************************************************;
/*15-variable best model missing 1 variable (R5EATA): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g stroke3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_adldepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_adldepdth2=0;
    time_adldepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;


**********************************************************;
/*15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi4g smoking3g stroke3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_adldepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_adldepdth2=0;
    time_adldepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;


**********************************************************;
/*14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_adldepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_adldepdth2=0;
    time_adldepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;

**********************************************************;
/*14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi4g smoking3g stroke3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_adl=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_adldepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_adldepdth2=0;
    time_adldepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_adldepdth2*status_adldepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;

/****************************************************************************************************************************************************************************************/
/*Walking*/

**********************************************************;
/*15-variable best model missing 1 variable (R5EATA): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g stroke3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_walkdepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_walkdepdth2=0;
    time_walkdepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;


**********************************************************;
/*15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi4g smoking3g stroke3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_walkdepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_walkdepdth2=0;
    time_walkdepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;


**********************************************************;
/*14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_walkdepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_walkdepdth2=0;
    time_walkdepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;

**********************************************************;
/*14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi4g smoking3g stroke3g*/

proc phreg data = dataclinical_gdr_20210624 covsandwich(aggregate); 
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  weight R5WTRESP; /*use survey weight*/
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /eventcode=1 risklimits;
  where subgroup_walk=1;
  output out=BSOUT xbeta=xb;
run;
data BSOUT; 
  set BSOUT; 
  if status_walkdepdth2=2 then do; /*the status and time need to be changed here since CONCORDANCE option below doesnt understand status=2*/
    status_walkdepdth2=0;
    time_walkdepdth2=19.1130249;
 end;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time_walkdepdth2*status_walkdepdth2(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;

/****************************************************************************************************************************************************************************************/
/*Death*/

**********************************************************;
/*15-variable best model missing 1 variable (R5EATA): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g stroke3g*/

proc surveyphreg data = dataclinical_gdr_20210624 ;
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  cluster RAEHSAMP; /*survey cluster*/
  strata RAESTRAT; /*survey strata*/
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /risklimits;
  output out=BSOUT xbeta=xb;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;


**********************************************************;
/*15-variable worst model missing 1 variable (R5WALKSA): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA bmi4g smoking3g stroke3g*/

proc surveyphreg data = dataclinical_gdr_20210624 ;
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  cluster RAEHSAMP; /*survey cluster*/
  strata RAESTRAT; /*survey strata*/
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /risklimits;
  output out=BSOUT xbeta=xb;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;


**********************************************************;
/*14-variable best model missing 2 variables (R5EATA stroke3g): lowest average TCIC*/
/*ageint5 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5LUNGE R5PUSHA R5WALKSA bmi4g smoking3g*/

proc surveyphreg data = dataclinical_gdr_20210624 ;
  class gender (ref='0') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5LUNGE (ref='0.no') R5PUSHA (ref='0.No') R5WALKSA (ref='0.No') smoking3g (ref='0');
  cluster RAEHSAMP; /*survey cluster*/
  strata RAESTRAT; /*survey strata*/
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5LUNGE R5PUSHA R5WALKSA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g /risklimits;
  output out=BSOUT xbeta=xb;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;

**********************************************************;
/*14-variable worst model missing 2 variables (R5WALKSA R5LUNGE): highest average TCIC*/
/*ageint5 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE R5PUSHA bmi4g smoking3g stroke3g*/

proc surveyphreg data = dataclinical_gdr_20210624 ;
  class gender (ref='0') R5EATA (ref='0.No') R5CANCRE (ref='0.no') R5DIABE (ref='0.no') heart3g (ref='0') R5HIBPE (ref='0.no')
                    R5MEALSA (ref='0.No') R5MONEYA (ref='0.No') H5LVALONE (ref='0.no') R5PUSHA (ref='0.No') smoking3g (ref='0') stroke3g (ref='0');
  cluster RAEHSAMP; /*survey cluster*/
  strata RAESTRAT; /*survey strata*/
  weight R5WTRESP; /*use survey weight*/
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender R5EATA R5CANCRE R5DIABE heart3g R5HIBPE R5MEALSA R5MONEYA H5LVALONE
                                               R5PUSHA bmi3_lin bmi3_sp1 bmi3_sp2 smoking3g stroke3g /risklimits;
  output out=BSOUT xbeta=xb;
run;
*Compute IAUC using IPCW  (Inverse probability of censoring weighting, Uno et al. 2007): it takes longer to compute, but it is the only one that provides SEs;
proc phreg data = BSOUT rocoptions(method=ipcw(cl seed=20210709) iauc outauc=aucdata); 
  class gender;
  model time2death*death(0) = age3_lin age3_sp1 age3_sp2 gender  / nofit;
  roc 'CoxWeighted'  pred=xb;
  ods output IAUC=iauc;
run;
proc means data=aucdata n mean stderr clm std p25 median p75 min max maxdec=4; var _AUC_ _StdErr_ _UpperAUC_ _LowerAUC_; run;
proc means data=aucdata stackodsoutput n mean maxdec=4;
 var _StdErr_;
 ods output summary=se(keep=Mean rename=Mean=SE_mean);
run;
data iauc2;
 merge iauc(keep=estimate rename=estimate=IAUC) se; /*only 1 obs in each data*/
 lcl_95=IAUC - quantile('normal' , 0.975)*SE_mean;
 ucl_95=IAUC + quantile('normal' , 0.975)*SE_mean;
proc print noobs data=iauc2; var IAUC lcl_95 ucl_95; run;
proc delete data=BSOUT aucdata se iauc iauc2 ; run; quit;
