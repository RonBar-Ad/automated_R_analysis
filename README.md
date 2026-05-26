# automated_R_analysis
R scripts that automate parts of statistical analysis. Aim to automated model selection and fitting for high-dimensionality ordinal-response dataset.

# File description:

DataPreparation reads two CSV files, a wide-format Qualtrics output and a long-format in-VR questionnaire output, and merges them, accounting for missing data and recoding values for readability to prepare for analysis, then exports the result into a merged CSV.

AutoSignificanceAnalysis reads a long-pivoted CSV and finds significant predictors for each outcome. It is designed for ordinal regression (Likert scale responses) with a random effect of participant ID, assuming data duplication in the pivot process. The variable that was pivoted is called Environment, and is included as a fixed effect predictor in all regression models to control its effect.

AutoModelSelection reads a wide CSV and pools factor levels, excluding columns with factor level frequencies below a given threshold, then fits ordinalNet and glmnet cross-validated LASSO-penalty dimension reduction for each outcome on all predictors. Each model is verified by fitting a mixed-effect ordinal regression. The results are then exported to a word document with model details. It is designed for ordinal regression (Likert scale responses) with a random effect of participant ID, assuming data duplication in the pivot process. The variable that was pivoted is called Environment, and is included as a fixed effect predictor in all regression models to control its effect.
