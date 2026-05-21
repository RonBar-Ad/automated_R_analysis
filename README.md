# automated_R_analysis
R scripts that automate parts of statistical analysis. Aim to automated model selection and fitting for high-dimensionality ordinal-response dataset.

# File description:

Automated_Analysis reads a long-pivoted CSV and finds significant predictors for each outcome. It is designed for ordinal regression (Likert scale responses) with a random effect of participant ID, assuming data duplication in the pivot process. The variable that was pivoted is called Environment, and is included as a fixed effect predictor in all regression models to control its effect.

Automated_LASSO reads a long-pivoted CSV and fits ordinalNet and glmnet cross-validated LASSO-penalty dimension reduction for each outcome on all predictors, then verifies by fitting a mixed-effect ordinal regression. The results are then exported to a word document with model details. It is designed for ordinal regression (Likert scale responses) with a random effect of participant ID, assuming data duplication in the pivot process. The variable that was pivoted is called Environment, and is included as a fixed effect predictor in all regression models to control its effect.
