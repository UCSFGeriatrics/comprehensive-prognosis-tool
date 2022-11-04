# Title: A comprehensive prognosis tool for older adults: Predicting Death, ADL Disability, and walking Disability simultaneously.
## Lee AK, Diaz-Ramirez LG, Boscardin WJ, Smith AK, Lee SJ. J Am Geriatr Soc. 2022 July 6. https://doi.org/10.1111/jgs.17932

## Description of data and SAS and R codes for reproducing the analyses of this article

### File name: 1.TCICIndividual_CompRiskSvy
File format: .sas
Description: SAS code to perform TCIC backward elimination by Outcome using HRS original dataset with time information and considering Competing-risks. It uses Cox regression for Death, and Competing-risk regression for rest of outcomes.

### File name: 2.baTCIC_CompRiskSvy
File format: .sas
Description: SAS code to perform best average TCIC backward elimination using HRS original dataset. It uses Cox regression for Death, and Competing-risk regression for rest of outcomes.

### File name: 3.BS500
File format: .sas
Description: SAS code to create 500 bootstrap samples.

### File name: 4a.BS100_TCIC_ByOutcomeSeladl
File format: .sas
Description: SAS code to perform TCIC backward elimination for ADL outcome using bootstrap datasets.

### File name: 4b.BS100_TCIC_ByOutcomeSelwalk
File format: .sas
Description: SAS code to perform TCIC backward elimination for Walking outcome using bootstrap datasets.

### File name: 4c.BS100_TCIC_ByOutcomeSeldeath
File format: .sas
Description: SAS code to perform TCIC backward elimination for Mortality outcome using bootstrap datasets.

### File name: 4d.BS100_TCIC_ByOutcomeSelAll
File format: .sas
Description: SAS code to create data with best models for individual selection by outcomes using bootstrap datasets and compute summary statistics.

### File name: 5.BS100_baTCICsel
File format: .sas
Description: SAS code to perform best average TCIC backward elimination using bootstrap datasets.

### File name: 6.BS100_baTCICoptimism
File format: .sas
Description: SAS code to compute optimism corrected statistics.

### File name: 7a.UnavailPred14baTCIC
File format: .sas
Description: SAS code to apply one unavailable predictor method to final baTCIC model.

### File name: 7b.UnavailPred14nCr2baTCIC
File format: .sas
Description: SAS code to apply two unavailable predictor metho to final baTCIC model.

### File name: 7c.UnavailPred_IAUC
File format: .sas
Description: SAS code to compute IAUC using IPCW method for best/worst case models missing 1 or 2 predictors.

### File name: 8a.SurvCIFestimates
File format: .sas
Description: SAS code to compute weighted predicted cumulative incidence for calibration plots.

### File name: 8b.CalibrationPlotsFinalModel
File format: .R
Description: R script to do calibration plots for final model.

### File name: 8c.CalibrationPlotsUnavaliablePred
File format: .R
Description: R script to do calibration plots for unavailable predictors models.

### File name: 9a.SurvBaseBetasEstimates
File format: .sas
Description: SAS code to compute weighted betas and baseline functions for Cox and Competing-risk final and unavailable predictor models.

### File name: Weightbasesurv_ADL
File format: .csv
Description: Data for ADL outcome with baseline survival function from time 1 to 19 for final, age-only, and unavailable predictor models.

### File name: Weightbetasurv_ADL
File format: .csv
Description: Data for ADL outcome with beta estimates for final, age-only, and unavailable predictor models.

### File name: Weightbasesurv_WALK
File format: .csv
Description: Data for Walking outcome with baseline survival function from time 1 to 19 for final, age-only, and unavailable predictor models.

### File name: Weightbetasurv_WALK
File format: .csv
Description: Data for Walking outcome with beta estimates for final, age-only, and unavailable predictor models.

### File name: Weightbasesurv_DEATH
File format: .csv
Description: Data for Mortality outcome with baseline survival function from time 1 to 19 for final, age-only, and unavailable predictor models.

### File name: Weightbetasurv_DEATH
File format: .csv
Description: Data for Mortality outcome with beta estimates for final, age-only, and unavailable predictor models.

### File name: 9b.PredictedRisks
File format: .R
Description: R script to compute predicted risks at 5, 10, and 14 years for ADL, Walking and Mortality outcomes.

### File name: 10.Subgroups_IAUC
File format: .sas
Description: SAS code to compute IAUC of Cox and Competing-risk final and unavailable predictors models for full sample and subgroups.

### File name: 11a.AveragePersonRisks
File format: .sas
Description: SAS code to compute "average for age" risk and compute Cumulative Incidence for calibration plots.

### File name: 11b.CreateRCsplines_age_bmi
File format: .R
Description: Compute restricted cubic splines for age and BMI.

### File name: 12.CalibrationPlotsAgeOnlyModel
File format: .R
Description: R script to do calibration plots for model with age only as predictor.
