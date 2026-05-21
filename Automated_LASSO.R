# Automated Dimension Reduction Ordinal Regressions
##### Author: Ron Bar-Ad
##### Last Update: 01/05/2026
##### Description:
# This is a script that reads a long-pivoted CSV and fits ordinalNet and glmnet
# cross-validated LASSO-penalty dimension reduction for each outcome on all
# predictors, then verifies by fitting a mixed-effect ordinal regression.
# The results are then exported to a word document with model details.
# It is designed for ordinal regression (Likert scale responses) with a random
# effect of participant ID, assuming data duplication in the pivot process. The
# variable that was pivoted is called Environment, and is included as a fixed
# effect predictor in all regression models to control its effect.
##### Note on use:
# To adapt this script for a different long-pivoted dataset, search "Environment"
# and replace all with the name of the pivot variable. Then change the outcomes
# list, and the rest should work as-is.
##### Process:
# 1. Fit glmnet LASSO with and without outcomes as potential predictors,
#     verify with mixed effect ordinal regression clmm.
# 2. Fit ordinalNet LASSO with and without outcomes as potential predictors,
#     verify with mixed effect ordinal regression clmm.
# 3. Export as word document.



library(tidyverse) # for reading csv and data manipulation
library(broom) # for data manipulation and presentation
library(dplyr) # for data manipulation and presentation
library(flextable) # for output to word doc, makes tables
library(officer) # for output to word doc, makes word doc
library(glmnet) # for lasso
library(ordinalNet) # for ordinal net best model
library(ordinal) # for mixed effects ordinal regression clmm


# Set the random seed for reproduceable cross-validation results.
set.seed(20070830)

# For import and export
setwd("C:\\Data")

# Import merged data from Qualtrics and in-VR questionnaire
data.all <- read_csv("Data.csv")

# Remove unnecessary information:
##    Anything qualitative (___Other, Qual___);
##    Anything empty (ExpStudUnsure, LivedExpEarlyNon, LivedExpEarlyUnsure,
###      LivedExpLateNon, LivedExpLateUnsure).
data.quant <- subset(data.all, select=-c(ExpWorkOther, ExpWorkPlanAndOther,
                                         ExpStudUnsure, ExpStudOther, ExpStudUniOther,
                                         ExpStudCurrentOther, LivedExpEarlyNon, LivedExpEarlyUnsure,
                                         LivedExpLateNon, LivedExpLateUnsure, QualProfessional,
                                         QualResearch, QualGeneral))
# Set nominal variables as factors
data.quant <- data.quant %>% 
  mutate_at(vars(-c(3,5,6,20:39, 43:60)), as.factor)

# Pivot data to long format
##   This will allow the use of Environment as a predictor with 3 levels
##   and each perception (dependent variables) as a single outcome,
##   but will duplicate rows from Qualtrics responses.
##   Shouldn't affect models, since data proportions are  the same,
##   but may interfere in analysis that doesn't involve the in-VR questionnaires.
data.long <- pivot_longer(data=data.quant, cols=43:60, names_to=c("Environment", ".value"), names_sep="_")

# Set environment factor levels. 0 = no trees, 1 = 315 trees, 2 = 634 trees.
data.long$Environment[data.long$Environment == "NT"] <- "0"
data.long$Environment[data.long$Environment == "ST"] <- "1"
data.long$Environment[data.long$Environment == "MT"] <- "2"
data.long$Environment <- as.factor(data.long$Environment)

# Define our outcome variables, this will come in handy when running our Lasso
outcomes <- c("Greenness", "Beauty", "Density", "Safety", "BuildingHeight", "RoadWidth")

# Function getHighLowLinear:
##   Returns a vector where all numeric values are divided into quartiles
##   based on all possible values. e.g. a 2 on a 1-to-7 scale will be Low,
##   even if all other values are 1 (which would also return Low).
##   Arguments:
###     col: vector, a dataframe column of numeric data
###     minim: num, the lower bound of possible answers
###     maxim: num, the upper bound of possible answers
getHighLowLinear <- function(col, minim, maxim) {
  
  bounds <- c(minim, maxim)
  
  # Set low point at 1st quartile of possible answers
  low_point <- summary(bounds)[["1st Qu."]]
  # Set mid point at mean. With only min and max as data points,
  # mean and median are the same.
  mid_point <- summary(bounds)[["Mean"]]
  # Set high point at 3rd quartile of possible answers.
  high_point <- summary(bounds)[["3rd Qu."]]
  
  # Replace values. For each value in the vector,
  # if it fits a given case, it is replaced by the given string,
  # otherwise move on to next.
  # If none match (impossible in this case), returns "Failed".
  case_when(
    # Up to 1st quartile is low
    col <= low_point ~ "Low",
    # 1st quartile to mean is midlow
    col > low_point & col <= mid_point ~ "MidLow",
    # mean to 3rd quartile is midhigh
    col > mid_point & col <= high_point ~ "MidHigh",
    # above 3rd quartile is high
    col > high_point ~ "High",
    # if this happens, something catastrophic has gone wrong.
    TRUE ~ "Failed"
  )
}


# Function getHighLowOnCurve:
##   Returns a vector where all numeric values are divided into quartiles
##   based on present values. e.g. a 2 on a 1-to-7 scale will be high
##   if all other values are 1 (which would return Low).
##   Arguments:
###     col: vector, a dataframe column of numeric data
##   Note: because this function relies on available data, it needs to handle
##   different amounts of variance. getHighLowLinear always has four quartiles,
##   but columns with no variance cannot be divided that way.
getHighLowOnCurve <- function(col) {
  
  # Handling different levels of variance.
  
  # Less than 2 unique values, everything is Medium.
  if (length(unique(col)) < 2) {
    return (rep("Medium"), length(col))
    
  # Exactly 2 unique values, the lower one is Low and the higher one is High.
  # values assigned string based on position in sorted version of vector.
  } else if (length(unique(col)) == 2) {
    case_when(col == sort(unique(col))[1] ~ "Low", col == sort(unique(col))[2] ~ "High", TRUE ~ "Failed")
    
  # Exactly 3 unique values, the lowest is Low, the middle is Medium, the highest is High.
  # values assigned string based on position in sorted version of vector.
  } else if (length(unique(col)) == 3) {
    u <- sort(unique(col))
    case_when(col == u[1] ~ "Low", col == u[2] ~ "Medium", col == u[3] ~ "High")
    
  # More than 3 values, we set quartiles and switch per case.
  } else {
    
    # Set low point at 1st quartile of data
    low_point <- summary(col)[["1st Qu."]]
    # Set mid point at mean of data.
    mid_point <- summary(col)[["Mean"]]
    # Set high point at 3rd quartile of data
    high_point <- summary(col)[["3rd Qu."]]
    
    # Replace values. For each value in the vector,
    # if it fits a given case, it is replaced by the given string,
    # otherwise move on to next.
    # If none match (impossible in this case), returns "Failed".
    case_when(
      # Up to 1st quartile is low
      col <= low_point ~ "Low",
      # 1st quartile to mean is midlow
      col > low_point & col <= mid_point ~ "MidLow",
      # mean to 3rd quartile is midhigh
      col > mid_point & col <= high_point ~ "MidHigh",
      # above 3rd quartile is high
      col > high_point ~ "High",
      # if this happens, something catastrophic has gone wrong.
      TRUE ~ "Failed"
    )
  }
}


# Define our high-low datasets

data.preds <- data.long[,setdiff(colnames(data.long),c(outcomes, "ID"))]
data.out   <- data.long[,setdiff(colnames(data.long), colnames(data.preds))]

# hiloC is the data divided by quartiles of each column's values.
# Apply function getHighLowOnCurve, but only to numeric columns,
# then turn result into factor.
data.hiloC <- lapply(data.preds, function(x) if (is.numeric(x)) {as.factor(getHighLowOnCurve(x))} else {as.factor(x)})
# Make data frame for analysis.
data.hiloC <- as.data.frame(data.hiloC)
data.hiloC <- merge(data.hiloC, data.out)

# Define the bounds for our different variables
# Columns on a 1-7 scale
Sevenpoint <- colnames(data.preds[,23:38])
# Columns with a 1-5 scale
Fivepoint <- colnames(data.preds[,19:22])
# Continuous value columns, no upper bound.
Continuous <- c("Age", "ExpWorkPlanYears", "ExpWorkOtherYears")

# hiloL is the data divided by quartiles of possible values.
# Apply function getHighLowLinear, but only to numeric columns.
# Continuous columns with no limit (e.g. Age) are divided on a curve.
data.hiloL <- as.data.frame(Map(function(x, nm) {
  
  # Likert scale questions get min value 1 and max 7
  if (nm %in% Sevenpoint) {
    as.factor(getHighLowLinear(x, 1, 7))
  
  # Five-point scale questions get min value 1 and max 5
  } else if (nm %in% Fivepoint) {
    as.factor(getHighLowLinear(x,1,5))
  
  # Continuous vars get separated on a curve
  } else if (nm %in% Continuous) {
    as.factor(getHighLowOnCurve(as.numeric(x)))
  
  # Outcomes are kept numeric
  } else if (nm %in% outcomes) {
    x
    
  # Anything non-numeric is kept as-is.  
  } else {
    as.factor(x)
  }
}, data.preds, names(data.preds)))
data.hiloL <- merge(data.hiloL, data.out)

# The 3 datasets to use. Will automatically loop through all 3.
#datasets <- list ("data.long" = data.long)
datasets <- list(
  "data.long" = data.long,
  "data.hiloC" = data.hiloC,
  "data.hiloL" = data.hiloL
)


##############################################################################
# Utility Functions                                                          #
##############################################################################

# Function map_to_original:
##   Maps expanded names back to original names. The regsubsets function
##   returns predictors with factor levels appended to their names, e.g.
##   GenderFemale and GenderMale as distinct predictors. This function finds
##   what the variables were called beforehand.
##   Arguments:
###     exp_name: str, possibly-expanded column name;
###     original_names: vector, names of predictors passed in args.
map_to_original <- function(exp_name, original_names) {
  
  # Names that are already fine get returned as they are.
  if (exp_name %in% original_names) return(exp_name)
  
  # Loop to check if the variable name starts with any original name.
  for (orig in original_names) {
    # If the expanded name starts with an original name, return
    if (startsWith(exp_name, orig)) {
      return(orig)
    }
  }
  # If the variable couldn't be found, something has gone horribly wrong.
  warning("Expanded variable name",exp_name,"could not be mapped to original variable name.")
  return(NULL)
}

# Function calculate_model_snr:
##   Calculates Signal-to-Noise Ratio (SNR) for a given linear regression model,
##   since LASSO is found to be more effective than best-subset at low SNR and
##   best-subset more effective than LASSO at high SNR.
##   Args:
###     model: lm object, a given linear regression.
calculate_model_snr <- function(model) {
  
  # Calculate Variances
  var_signal <- var(fitted(model))
  var_noise <- var(residuals(model))
  
  # If noise variance is 0 (perfect fit), return infinite positive value.
  if (var_noise == 0 || is.na(var_noise)) {
    return(Inf) 
  }
  
  # Otherwise, return SNR.
  return(var_signal / var_noise)
}


remove_intercepts <- function(sum, preds) {
  # Define iterator.
  i <- 1
  
  # Iterate through rows of model details.
  # Each row name is a variable,
  # but some have factor levels appended to the name.
  # For each name, check if it starts with any predictor,
  # or if it already is a predictor,
  # and if not, remove it from the list.
  # While loop runs as long as there are Intercepts remaining.
  while (i <= nrow(sum)) {
    
    # Set removal condition.
    flag <- FALSE
    
    # Cycle through predictors, checking if row name
    # is in the list / starts with a predictor name (is a factor level)
    for (j in preds) {
      
      # startsWith will return TRUE for two equal strings
      # so it also works as if ([rowname] %in% predictors)
      if (startsWith(row.names(sum)[i], j)){
        
        # If it's a predictor, we signal not to remove it,
        # then exit the for loop.
        flag <- TRUE
        next
      }
    }
    # If it wasn't a variable in the predictors list, remove it.
    if (!flag) {
      
      # Remove the row. Treating summ as data.frame for cases where
      # only 1 column is non-intercept, which would remove the
      # rownames completely and flatten to vector if it were a matrix
      sum <- as.data.frame(sum[-c(i),])
      
      # The next row in the matrix is now at the index this one was in,
      # so we move our iterator back to keep it where it is.
      i <- i - 1
    }
    
    # Move iterator on to next row.
    i <- i + 1
  }
  return(sum)
}


##############################################################################
# Export                                                                     #
##############################################################################

format_p <- function(p, threshold = 0.01) {
  ifelse(p < threshold, 
         return(sprintf("< %.2f", threshold)),
         return(sprintf("%.3f", p))) 
}

# Function export_results:
##   Formats results of best-subsets and lassos as tables with flextable,
##   then exports to word document with officer.
##   Arguments:
###     results_list: list, contains results of all selected models;
###     output_file: str, name of file to export, default Automated_Report.docx,
export_results <- function(results_list, output_file = "Automated_Report_LASSO.docx") {
  
  # Initialising document
  doc <- read_docx()
  
  # Title
  title_fpar <- fpar("Automated Regression Analysis Report")
  doc <- body_add_fpar(doc, title_fpar, style="Normal")
  
  # Subtitle with date and time.
  subtitle_fpar <- fpar(paste("Generated on:", Sys.time()))
  doc <- body_add_fpar(doc, subtitle_fpar, style="Normal")
  doc <- body_add_break(doc)
  
  # Get datasets in list
  # Processing from datasets used in case they for some reason conflict with
  # datasets in global list.
  datasets_processed <- unique(sapply(names(results_list), function(x) strsplit(x, "_")[[1]][1]))
  
  # Loop through list of datasets, 
  for (ds_name in datasets_processed) {
    
    # Only see results for this dataset.
    ds_results <- results_list[names(results_list) %>% grep(ds_name, ., value = TRUE)]
    
    # Add header for dataset
    ds_header <- fpar(paste("Dataset:", ds_name))
    doc <- body_add_fpar(doc, ds_header, style = "heading 1")
    doc <- body_add_break(doc)
    
    # Loop through each model in dataset
    for (key in names(ds_results)) {
      # Split the model into parts
      res <- ds_results[[key]]
      outcome <- strsplit(key, "_")[[1]][2]
      mode <- strsplit(key, "_")[[1]][3]  # "WITH", "WITHOUT", "LASSO", "ON"
      
      # Add outcome as header
      outcome_header <- fpar(paste("Outcome:", outcome, "| Method:", mode))
      doc <- body_add_fpar(doc, outcome_header, style = "heading 2")
      
      if (mode == "ON") {
        ##### Format for ordinal net
        
        # If there's a note, write it in and skip to the next model
        if (!is.null(res$note)) {
          note_txt <- fpar(paste("Note:", res$note))
          doc <- body_add_fpar(doc, note_txt, style = "Normal")
          doc <- body_add_break(doc)
          next
        }
        
        # Add deviance and misclassification
        devPct_txt <- fpar(paste("Optimal deviance percentage:",
                                 round(res$devPct, 5)))
        doc <- body_add_fpar(doc, devPct_txt, style = "Normal")
        
        misclass_txt <- fpar(paste("Misclassification:", round(res$misclass, 4)))
        doc <- body_add_fpar(doc, misclass_txt, style = "Normal")
        
        # Add predictors list
        if (res$model_size > 0) {
          pred_txt <- fpar(paste("Selected Predictors (", res$model_size, "): ", 
                                 paste(res$selected_predictors, collapse = ", ")))
          doc <- body_add_fpar(doc, pred_txt, style = "Normal")
        } else {
          pred_txt <- fpar("Selected Predictors: None (Intercept only)")
          doc <- body_add_fpar(doc, pred_txt, style = "Normal")
        }
        
        # Add formula
        if (!is.null(res$formula)) {
          form_txt <- fpar(paste("Formula:", res$formula))
          doc <- body_add_fpar(doc, form_txt, style = "Normal")
        }
        
        # Add predictors to table
        if (nrow(res$significant_df) > 0) {
          res$significant_df$estimate <- round(res$significant_df$estimate,2)
          res$significant_df$p.value <- sapply(res$significant_df$p.value, format_p)
          significants <- res$significant_df[res$significant_df$p.value < 0.05,]
          if (nrow(significants) > 0){
            ft <- flextable(res$significant_df[res$significant_df$p.value < 0.05,])
            ft <- theme_vanilla(ft)
            ft <- set_caption(ft, caption = "Significant Predictors (p < 0.05) from Ordinal Refit")
            doc <- body_add_flextable(doc, ft)
          } else {
            note_txt <- fpar("No significant predictors found in ordinal refit.")
            doc <- body_add_fpar(doc, note_txt, style = "Normal")
          }
        } else {
          note_txt <- fpar("No significant predictors found in OLS refit (p >= 0.05).")
          doc <- body_add_fpar(doc, note_txt, style = "Normal")
        }
        
      } else if (mode == "LASSO") {
        ##### Format for LASSO
        
        # If there's a note, write it in and skip to the next model
        if (!is.null(res$note)) {
          note_txt <- fpar(paste("Note:", res$note))
          doc <- body_add_fpar(doc, note_txt, style = "Normal")
          doc <- body_add_break(doc)
          next
        }
        
        # Add lambda and cross-validation error
        lambda_txt <- fpar(paste("Optimal Lambda:", round(res$lambda, 5)))
        doc <- body_add_fpar(doc, lambda_txt, style = "Normal")
        
        cv_txt <- fpar(paste("Cross-Validated MSE:", round(res$cv_error, 4)))
        doc <- body_add_fpar(doc, cv_txt, style = "Normal")
        
        # Add predictors list
        if (res$model_size > 0) {
          pred_txt <- fpar(paste("Selected Predictors (", res$model_size, "): ", 
                                 paste(res$selected_predictors, collapse = ", ")))
          doc <- body_add_fpar(doc, pred_txt, style = "Normal")
        } else {
          pred_txt <- fpar("Selected Predictors: None (Intercept only)")
          doc <- body_add_fpar(doc, pred_txt, style = "Normal")
        }
        
        # Add formula
        if (!is.null(res$formula)) {
          form_txt <- fpar(paste("Formula:", res$formula))
          doc <- body_add_fpar(doc, form_txt, style = "Normal")
        }
        
        res$significant_df <- res$significant_df[!is.na(res$significant_df$p.value),]
        
        # Add predictors to table
        if (nrow(res$significant_df) > 0) {
          res$significant_df$estimate <- round(res$significant_df$estimate,2)
          res$significant_df$p.value <- sapply(res$significant_df$p.value, format_p)
          if (any(res$significant_df$p.value < 0.05)) {
            ft <- flextable(res$significant_df[res$significant_df$p.value < 0.05,])
            ft <- theme_vanilla(ft)
            ft <- set_caption(ft, caption = "Significant Predictors (p < 0.05) from Ordinal Refit")
            doc <- body_add_flextable(doc, ft)
          } else {
            note_txt <- fpar("No significant predictors found in ordinal refit.")
            doc <- body_add_fpar(doc, note_txt, style = "Normal")
          }
        } else {
          note_txt <- fpar("No significant predictors found in ordinal refit (p >= 0.05).")
          doc <- body_add_fpar(doc, note_txt, style = "Normal")
        }
        
      } else {
        ##### Best-Subset
        # Both with and without outcomes as predictors
        
        # Loop through the 3 metrics (Adj R-sq, Cp, BIC)
        for (metric in names(res$models)) {
          model_info <- res$models[[metric]]
          
          # Metric Header
          metric_header <- fpar(paste("Best Model by", metric))
          doc <- body_add_fpar(doc, metric_header, style = "heading 3")
          
          # Formula and Size
          formula_txt <- fpar(paste("Formula:", model_info$formula))
          doc <- body_add_fpar(doc, formula_txt, style = "Normal")
          
          size_txt <- fpar(paste("Number of Predictors:", model_info$size))
          doc <- body_add_fpar(doc, size_txt, style = "Normal")
          
          # Significant Predictors Table
          if (nrow(model_info$significant_df) > 0) {
            model_info$significant_df$estimate <- round(model_info$significant_df$estimate,2)
            model_info$significant_df$p.value <- sapply(model_info$significant_df$p.value, format_p)
            ft <- flextable(model_info$significant_df)
            ft <- theme_vanilla(ft)
            ft <- set_caption(ft, caption = "Significant Predictors (p < 0.05)")
            doc <- body_add_flextable(doc, ft)
          } else {
            note_txt <- fpar("No significant predictors found (p >= 0.05).")
            doc <- body_add_fpar(doc, note_txt, style = "Normal")
          }
          
          doc <- body_add_break(doc)
        }
      }
      
      # Add page break between outcomes
      doc <- body_add_break(doc)
    }
  }
  
  # 4. Save
  print(doc, target = output_file)
  cat("Unified report saved to:", output_file, "\n")
}

##############################################################################
# Analysis Functions                                                         #
##############################################################################

# Function get_lasso_results.
##   Runs a lasso to find the best-fit model, runs the model regression,
##   returns a list of significant predictors, effect sizes, p-values, and
##   model details, all to be turned into a table later. Includes error handling.
##   Arguments:
###     data: dataframe, the dataset being analysed;
###     outcome_var: str, the name of the variable being predicted;
###     alpha: num, the penalty applied to predictors,
####       set to 1 for lasso or 0 for ridge regression;
###     nfolds: num, number of folds for cross validation,
####       defaults to 10 as industry standard;
###     exclude_outcomes: bool, whether to exclude dependent variables as potential predictors,
####       defaults FALSE, but automation runs once with FALSE once with TRUE;
###     all_outcomes: str vector, list of names of other oucomes,
####       used to exclude/include, but also for error handling checks.
get_lasso_results <- function(data.lasso, outcome_var, alpha.lasso = 1, nfolds = 10,
                              exclude_outcomes = FALSE, all_outcomes = outcomes) {
  
  ########## 1. Transform dataframe to matrix.
  
  ##### 1.1. Outcome variable.
  
  # Check if outcome variable is included in dataframe.
  if (!(outcome_var %in% colnames(data.lasso))) {
    # If not, say so and exit function.
    warning(paste("Outcome", outcome_var, "not found."))
    return(NULL)
  }
  
  # Get outcome variable as numeric.
  # Remove any frame wrappers, value-key pairs, etc.
  y <- as.numeric(unlist(data.lasso[[outcome_var]]))
  
  # If number of non-NA unique values in outcome variable is 1 or 0, there won't
  # be enough variance for a LASSO, so return empty list. Page of word doc will
  # be empty.
  if (sum(!is.na(unique(y))) <= 1) {
    warning(paste("Outcome ", outcome_var, " has <= 1 non-NA unique values."))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, 
                note = paste("Outcome has only", n_unique, "unique value(s)")))
  }
  
  # If more than half of the outcome variable's values are missing (NA),
  # results would be meaningless, and LASSO won't run, so return empty list.
  # Page of word doc will be empty.
  if (sum(is.na(y)) >= length(y) * 0.5) {
    warning(paste("Outcome ", outcome_var, " is more than half NA values."))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, 
                note = "Outcome has >50% missing values"))
  }
  
  ##### 1.2. Predictors.
  
  # Get list of predictor variable names.
  # If we're excluding outcomes as predictors, we use the list of all outcomes
  # to create a list of column names excluding them. Otherwise, we use the
  # name of the outcome we're looking to predict and exclude only it.
  ifelse(exclude_outcomes, predictors <- setdiff(colnames(data.lasso), all_outcomes), 
         predictors <- setdiff(colnames(data.lasso), outcome_var))
  
  # Exclude from list any columns that:
  ##   have less than 2 non-NA values (not enough variance for LASSO)
  ##   have more NAs than actual values (meaningless, won't work)
  predictors <- predictors[sapply(data.lasso[predictors], function(x) {
    sum(!is.na(unique(x))) > 1 && sum(is.na(x)) < length(x) * 0.5
  })]
  
  # Remove ID, to be treated as random effect in mixed ordinal regression
  predictors <- predictors[predictors != "ID"]
  
  # If we've exculded all out predictors, return empty list, there's nothing
  # to run lasso on.
  if (length(predictors) == 0) {
    warning(paste("Outcome ", outcome_var, " has no predictors."))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, 
                note = "No valid predictors"))
  }
  
  # Get predictor values as vector.
  X <- data.lasso[, predictors, drop = FALSE]
  
  # Convert to numeric for lasso,
  # remove any frame wrappers, value-key pairs, etc.
  # Factors get turned into level codes (e.g. 1,2,3)
  X <- lapply(X, function(col) as.numeric(unlist(col)))
  
  # Make sure no NAs get in the function
  # Only select rows with no missing values
  complete_rows <- complete.cases(X, y)
  
  # If fewer than 10 rows without NA values, not enough predictors
  # Word doc will have blank page.
  if (sum(complete_rows) < 10) {
    warning(paste("Outcome ", outcome_var, " had too few complete cases after making numeric."))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, 
                note = "Not enough complete cases for cross-validation"))
  }
  
  # Select only complete rows
  # + make sure X is a matrix
  X <- as.matrix(as.data.frame(X)[complete_rows, ])
  y <- y[complete_rows]
  
  # If the predictor matrix isn't a matrix somehow, LASSO won't work.
  # Word doc will have a blank page.
  if (!is.numeric(X)) {
    warning(paste("Matrix conversion failed for", outcome_var, "after cleaning. Skipping."))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, 
                note = "Matrix conversion failed"))
  }
  
  
  ########## 2. Run LASSO.
  
  ##### 2.1. Set up penalties.
  # We want Environment to have no penalty so it's
  # included as a control more often.
  
  # Find where Environment is in predictors matrix.
  if ("Environment" %in% colnames(X)){
    index <- which(colnames(X) == "Environment")
  } else {
    warning(paste("Environment not included in predictors for", outcome_var))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, 
                note = "Environment not included in predictors."))
  }
  
  # Set penalty to alpha passed in args, likely 1 for LASSO
  penalties <- c(rep(alpha.lasso, length(colnames(X))))
  # Set penalty for Environment to 0.
  penalties[index] <- 0
  
  ##### 2.2. Define cross-validation folds.
  # We want our k-fold cross validation not to draw from the same
  # participant over and over, causing data leakage. Instead,
  # we make sure to get folds from grouped data per participant ID.
  
  # Get participant IDs.
  obs <- unique(data.lasso$ID)
  # Randomly shuffle IDs.
  obs_rand <- sample(obs)
  # Assign each ID into K groups (passed in args, default 5).
  f_assign <- rep(1:nfolds, length.out = length(obs))
  
  # Make a map to match IDs to their fold numbers
  fold_map <- data.frame(ID = obs_rand, FoldNumber = f_assign)
  # Map the fold numbers onto the row numbers by ID
  merged_folds <- merge(data.lasso[,"ID",drop=FALSE], fold_map,
                     by = "ID", sort=FALSE)
  # Turn into a vector for cv.glmnet.
  folds_v <- merged_folds$FoldNumber
  
  ##### 2.3. Fit the LASSO.
  
  # Try to fit a lasso:
  ##   cv.glmnet is a LASSO function with cross-validation.
  ##   Pass arguments:
  ###     X: matrix (int), all predictor values
  ###     y: vector (int), outcome variable values
  ###     alpha: num, from function args, penalty applied to predictors
  ###     nfolds: num, from function args, number of cross-validation folds
  cv_fit <- tryCatch({
    cv.glmnet(X, y, alpha = alpha.lasso, foldid = folds_v,
              penalty.factor = penalties)
    
    # If it fails, print error message and return empty list.
    # Word doc will have a blank page.
  }, error = function(e) {
    warning(paste("LASSO failed for", outcome_var, ":", e$message))
    return(NULL)
  })
  
  # If function threw no error but lasso still failed, return empty list.
  # Word doc will have a blank page.
  if (is.null(cv_fit)) {
    warning(paste("Outcome ", outcome_var, " glmnet returned null."))
    return(list(method = "LASSO", lambda = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = NA, note = "CV glmnet failed"))
  }
  
  
  
  ########## 4. Return list of LASSO best model details.
  
  ##### 4.1. Refit linear regression of best model.
  
  # Get best model performance.
  lambda_opt <- cv_fit$lambda.min
  # Get best model predictors.
  coef_opt <- coef(cv_fit, s = "lambda.min")
  
  # Get names of predictors:
  ##   For factors converted to numeric level codes,
  ##   variable names are the original column names.
  selected_vars <- rownames(coef_opt)[-1][coef_opt[-1, 1] != 0]
  
  # Only run linear regression if there are predictors to choose from.
  if (length(selected_vars) > 0) {
    
    # Ensure environment is included as fixed effect.
    if (!"Environment" %in% selected_vars) selected_vars <- c(selected_vars,
                                                              "Environment")
    
    # Make formula from outcome variable and the LASSO best model predictors
    formula_str <- paste("as.ordered(",outcome_var, ") ~",
                         paste(selected_vars, collapse = " + "),
                         " + (1 | ID)")

    
    # Try to fit linear model with above formula.
    final_lm <- tryCatch({
      #lm(as.formula(formula_str), data = data.lasso)
      clmm(as.formula(formula_str), data = data.lasso)

      # If it fails, print error and return empty list.
      # Word doc will have a blank page.
    }, error = function(e) {
      warning(paste("OLS refit failed:", e$message))
      return(NULL)
    })
    
    ##### 4.2. Return model details as list.
    
    # Only return details if the model didn't fail.
    if (!is.null(final_lm)) {
      
      if (!"Environment" %in% predictors) predictors <- c(predictors, "Environment")
      
      sig_coefs <- remove_intercepts(summary(final_lm)$coefficients, predictors)
      
      # Format results for table output.
      sig_coefs <- data.frame(
        term = row.names(sig_coefs),
        estimate = sig_coefs[,1],
        p.value = sig_coefs[,4],
        row.names = NULL
      )
      
      # Return a list with the details of the model:
      ##   method: clarify LASSO to juxtapose the best-subset models;
      ##   exclude_outcomes: clarify whether outcomes were excluded or not;
      ##   lambda: performance of best model
      ##   selected_predictors: which predictors were included in the best model
      ##   model_size: number of predictors in the best model
      ##   formula: formula of best model
      ##   significant_df: predictor effects and p-values
      ##   cv_error: cross-validation performance.
      return(list(method = "LASSO", exclude_outcomes = exclude_outcomes, 
                  lambda = lambda_opt, selected_predictors = selected_vars, 
                  model_size = length(selected_vars), formula = formula_str, 
                  significant_df = sig_coefs, cv_error = cv_fit$cvm[cv_fit$lambda == lambda_opt]))
      
      # If the linear model was empty, return empty list.
    } else {
      warning(paste("Outcome ", outcome_var, " OLS refit empty."))
      return(list(method = "LASSO", exclude_outcomes = exclude_outcomes,
                  lambda = lambda_opt, selected_predictors = selected_vars, 
                  model_size = length(selected_vars), formula = formula_str, 
                  significant_df = data.frame(), cv_error = cv_fit$cvm[cv_fit$lambda == lambda_opt],
                  note = "OLS refit failed"))
    }
    
    # If the list of predictors in the best model was empty, return empty list.
  } else {
    warning(paste("Outcome ", outcome_var, " LASSO selected no predictors."))
    return(list(method = "LASSO", exclude_outcomes = exclude_outcomes,
                lambda = lambda_opt, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), cv_error = cv_fit$cvm[cv_fit$lambda == lambda_opt],
                note = "LASSO selected no predictors"))
  }
}

# Function get_best_subset:
##   Runs a best-subset selection to find the best-fit model,
##   runs the model regression, and returns a list of significant predictors,
##   effect sizes, p-values, and model details, all to be turned
##   into a table later. Includes error handling.
##   Arguments:
###     data: dataframe, the dataset being analysed;
###     outcome_var: str, the name of the variable being predicted;
###     nvmax: num, the maximum size of subsets to consider, non-MIO best-subset
####       selection is considered incalculable at > 30 predictors, therefore
####       defaults to NULL, which causes error, value must be set in function call;
###     method: str, which type of selection to use, default is forward as exhaustive
####       is impossible with > 30 predictors;
###     exclude_outcomes: bool, whether to exclude outcomes from predictors;
###     all_outcomes, vector, contains column names of dependent variables.
get_best_subset <- function(data.bs, outcome_var, nvmax = NULL, 
                            method = "forward", alpha.bs = 0.5,
                            exclude_outcomes = FALSE, all_outcomes = outcomes) {
  
  ########## 1. Find best-subset model.
  
  # Check if outcome variable is included in dataframe.
  if (!(outcome_var %in% colnames(data.bs))) {
    # If not, say so and exit function.
    warning(paste("Outcome", outcome_var, "not found."))
    return(NULL)
  }
  
  # Get list of predictor variable names.
  # If we're excluding outcomes as predictors, we use the list of all outcomes
  # to create a list of column names excluding them. Otherwise, we use the
  # name of the outcome we're looking to predict and exclude only it.
  ifelse(exclude_outcomes, predictors <- setdiff(colnames(data.bs), all_outcomes), 
         predictors <- setdiff(colnames(data.bs), outcome_var))
  
  if (length(predictors) == 0) {
    warning(paste("No predictors available for", outcome_var, "after exclusions."))
    return(NULL)
  }
  
  # Create formula with all predictor names
  form <- as.formula(paste(outcome_var, "~", paste(predictors, collapse = " + ")))
  
  # Try to find best subset.
  ##   Pass arguments:
  ###     form: formula, containing every predictor;
  ###     data: dataframe, from args;
  ###     nvmax: num, from args;
  ###     method: str, from args;
  ###     really.big: bool, TRUE dismisses warnings about long calculation time.
  rs_fit <- tryCatch({
    regsubsets(form, data = data.bs, 
               nvmax = nvmax, 
               method = method,
               really.big = TRUE)
    
    # if failed, print error message and return empty list.
  }, error = function(e) {
    warning(paste("Failed to run regsubsets for", outcome_var, ":", e$message))
    return(NULL)
  })
  
  # If failed but no error, return empty list.
  if (is.null(rs_fit)) {
    warning("Outcome",outcome_var," regsubsets returned empty object.")
    return(NULL)
  } 
  
  
  
  ########## 2. Fit findings to linear regression.
  
  ##### 2.1. Set up to build formula.
  
  # The results of the best-subset selection
  summ <- summary(rs_fit)
  
  # Find best model by AdjR2, CP, and BIC
  best_indices <- c(
    "Adj R-sq" = which.max(summ$adjr2),
    "Cp"       = which.min(summ$cp),
    "BIC"      = which.min(summ$bic)
  )
  
  # Get the expanded variable names in the best-subset model,
  # these will need to be mapped back to original column names.
  expanded_names <- colnames(summ$which)
  
  ##### 2.2. Fit regressions for each metric's best model.
  
  # List to hold results.
  results_list <- list()
  
  # Fit a different regression for each metric's best subset of predictors.
  # Cycling through AdjR2, CP, and BIC.
  # Uses iterator "metric" that just holds a num (0,1,2)
  for (metric in names(best_indices)) {
    
    # Get which metric we're using
    idx <- best_indices[metric]
    
    # Get the model's result for that metric
    selected <- summ$which[idx, ]
    
    # Get the expanded variable names to later be mapped back.
    # This is now the list of predictors for the metric: the best subset.
    selected <- expanded_names[selected]
    
    # Remove intercept.
    if ("(Intercept)" %in% selected){
      selected <- selected[selected != "(Intercept)"]
    }
    
    # If the model returned no predictors, fit a null model.
    if (length(selected) == 0) {
      model_formula <- paste(outcome_var, "~ 1")
      fit <- lm(as.formula(model_formula), data = data.bs)
      
      # If the model did return predictors, get fit a model for them.
    } else {
      # Map factor level expanded variable names to original ones.
      selected <- sapply(selected, map_to_original, original_names = predictors)
      
      # Multiple levels of one factor will create duplicates, remove them.
      selected <- unique(selected)
      
      # Create formula from selected predictors.
      model_formula <- as.formula(paste(outcome_var, "~", paste(selected, collapse = " + ")))
      
      # Try to fit the linear model.
      fit <- tryCatch({
        # Using the formula of predictors from the best-subset selection.
        lm(model_formula, data = data.bs)
        
      # On failure, report warning and return null.
      }, error = function(e) {
        warning(paste("Failed to fit linear model:", model_formula, "\nError:", e$message))
        return(NULL)
      })

      # If model fails with no error, return null.
      if (is.null(fit)) {
        warning(paste("Model", model_formula, "returned null."))
        return(NULL)
      }
    }
    
    
    
    ########## 3. Return results in list to output as table.
    
    # Get only the predictors, effects, and p-values of significant predictors.
    sig_coefs <- tidy(summary(fit)) %>%
      filter(term != "(Intercept)", p.value < 0.05) %>%
      dplyr::select(term, estimate, p.value)
    
    # Store result in list, to iterate to next metric.
    results_list[[metric]] <- list(
      formula = model_formula,
      size = length(selected),
      predictors = selected,
      significant_df = sig_coefs
    )
  }
  
  # Return the results of the 3 metrics.
  return(list(
    outcome = outcome_var,
    exclude_outcomes = exclude_outcomes,
    dataset_name = deparse(substitute(data.bs)),
    models = results_list
  ))
}


# Function get_ordinal_net:
##   Runs an ordinal net to find the best-fit ordinal regression model,
##   runs that ordinal regression model, and returns a list of significant
##   predictors, effect sizes, p-values, and model details, all to be turned
##   into a table later. Includes error handling.
##   Arguments:
###     data: dataframe, the dataset being analysed;
###     outcome_var: str, the name of the variable being predicted;
###     exclude_outcomes: bool, whether to exclude outcomes from predictors;
###     nfolds.on: num, number of cross-validation folds;
###     nfolds.on.cv: num, number of cross-validation folds used to tune lambda;
###     all_outcomes: vector (str), column names of all dependent variables.
get_ordinal_net <- function(data.on, outcome_var, exclude_outcomes = FALSE,
                            folds.on = 5, alpha.on = 1,
                            all_outcomes = outcomes) {
  
  ########## 1. Transform dataframe to matrix.
  
  ##### 1.1. Outcome variable.
  
  # Check if outcome variable is included in dataframe.
  if (!(outcome_var %in% colnames(data.on))) {
    # If not, say so and exit function.
    warning(paste("Outcome", outcome_var, "not found."))
    return(NULL)
  }
  
  # Get outcome variable as vector.
  # Remove any frame wrappers, value-key pairs, etc.
  # Make sure it's ordinal.
  y <- as.ordered(unlist(data.on[[outcome_var]]))
  
  # If number of non-NA unique values in outcome variable is 1 or 0, there won't
  # be enough variance for a LASSO, so return empty list. Page of word doc will
  # be empty.
  if (sum(!is.na(unique(y))) <= 1) {
    warning(paste("Outcome ", outcome_var, " has <= 1 non-NA unique values."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = paste("Outcome has only", n_unique, "unique value(s)")))
  }
  
  # If more than half of the outcome variable's values are missing (NA),
  # results would be meaningless, and ordinalNet won't run, so return empty list.
  # Page of word doc will be empty.
  if (sum(is.na(y)) >= length(y) * 0.5) {
    warning(paste("Outcome ", outcome_var, " is more than half NA values."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA, 
                note = "Outcome has >50% missing values"))
  }
  
  ##### 1.2. Predictors.
  
  # Get list of predictor variable names.
  # If we're excluding outcomes as predictors, we use the list of all outcomes
  # to create a list of column names excluding them. Otherwise, we use the
  # name of the outcome we're looking to predict and exclude only it.
  ifelse(exclude_outcomes, predictors <- setdiff(colnames(data.on), all_outcomes), 
         predictors <- setdiff(colnames(data.on), outcome_var))
  
  # Remove ID, to be treated as random effect in mixed ordinal regression
  predictors <- predictors[predictors != "ID"]
  
  # Make sure environment is included if it somehow wasn't.
  if (!"Environment" %in% predictors) predictors <- c(predictors, "Environment")
  
  # Exclude from list any columns that:
  ##   have less than 2 non-NA values (not enough variance for LASSO)
  ##   have more NAs than actual values (meaningless, won't work)
  predictors <- predictors[sapply(data.on[predictors], function(x) {
    sum(!is.na(unique(x))) > 1 && sum(is.na(x)) < length(x) * 0.5
  })]
  
  # If we've exculded all out predictors, return empty list, there's nothing
  # to run lasso on.
  if (length(predictors) == 0) {
    warning(paste("Outcome ", outcome_var, " has no predictors."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = "No valid predictors"))
  }
  
  # Get predictor values as vector.
  X <- data.on[, predictors, drop = FALSE]
  # Convert to numeric for ordinalNet,
  # remove any frame wrappers, value-key pairs, etc.
  # Factors get turned into level codes (e.g. 1,2,3)
  X <- lapply(X, function(col) as.numeric(unlist(col)))
  # Make it a matrix for sure
  #X <- do.call(cbind, X)
  
  # Make sure no NAs get in the function
  # Only select rows with no missing values
  complete_rows <- complete.cases(X, y)
  
  # If fewer than 10 rows without NA values, not enough predictors
  # Word doc will have blank page.
  if (sum(complete_rows) < 10) {
    warning(paste("Outcome ", outcome_var, " had too few complete cases after making numeric."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = "Not enough complete cases for cross-validation"))
  }
  
  # Select only complete rows
  # + make sure X is a matrix
  X <- as.matrix(as.data.frame(X)[complete_rows, ])
  y <- y[complete_rows]
  
  # If the predictor matrix isn't a matrix somehow, LASSO won't work.
  # Word doc will have a blank page.
  if (!is.numeric(X)) {
    warning(paste("Matrix conversion failed for", outcome_var, "after cleaning. Skipping."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = "Matrix conversion failed"))
  }
  
  
  
  ########## 2. Run Ordinal Net

  
  ##### 2.1. Set up penalties.
  # We want Environment to have no penalty so it's
  # included as a control more often.
  
  # Find where Environment is in predictors matrix.
  if ("Environment" %in% colnames(X)){
    index <- which(colnames(X) == "Environment")
  } else {
    warning(paste("Environment not included in predictors for", outcome_var))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = "Environment not included in predictors."))
  }
  
  # Set penalty to alpha passed in args, likely 1 for LASSO
  penalties <- c(rep(alpha.on, length(colnames(X))))
  # Set penalty for Environment to 0.
  penalties[index] <- 0
  
  ##### 2.2. Define cross-validation folds.
  # We want our k-fold cross validation not to draw from the same
  # participant over and over, causing data leakage. Instead,
  # we make sure to get folds from grouped data per participant ID.
  
  # Get participant IDs.
  obs <- unique(data.on$ID)
  # Randomly shuffle IDs.
  obs_rand <- sample(obs)
  # Assign each ID into K groups (passed in args, default 5).
  f_assign <- rep(1:folds.on, length.out = length(obs))
  
  # Make empty list of dimensions NumRows x K.
  row_folds <- vector("list", folds.on)
  # For each group, add the row number corresponding to each of the IDs
  # within that group (e.g. ID 1 is on rows 1,2,3, so whichever group
  # it's in will include 1,2,3)
  for (i in 1:folds.on){
    row_folds[[i]] <- which(data.on$ID %in% obs_rand[f_assign == i])
  }
  
  ##### 2.3. Fit the ordinalNet.
  
  # Try to fit an ordinalNet:
  ##   ordinalNetCV is an ordinalNet function with cross-validation.
  ##   Pass arguments:
  ###     X: matrix (int), all predictor values
  ###     y: vector (int), outcome variable values
  ###     nFolds: num, number of cross-validation folds.
  ###     nFoldsCV: num, number of cross-validation folds used to tune lambda
  ###     maxiterOut: num, max number of outer-loop iterations. set to 500 to
  ####       ensure precise deviance assessment
  ###     alpha: num, penalty for predictors. Set to 1 for LASSO penalty.
  on_fit <- tryCatch({
    # ordinalNetCV(X, y, nFolds = folds.on, nFoldsCV = folds.on.cv,
    #              maxiterOut = 500, alpha = 1)
    ordinalNetTune(X, y, folds = row_folds, family = "cumulative",
                   link = "logit", penaltyFactors = penalties,
                   parallelTerms = TRUE, nonparallelTerms = FALSE)
    
    # If it fails, print error message and return empty list.
    # Word doc will have a blank page.
  }, error = function(e) {
    warning(paste("ordinal net failed for", outcome_var, ":", e$message))
  })
  
  
  
  # If function threw no error but ordinalNet still failed, return empty list.
  # Word doc will have a blank page.
  if (is.null(on_fit)) {
    warning(paste("Outcome ", outcome_var, " ordinalNetCV returned null."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = "ordinalNetCV failed"))
  }
  
  # If best model does not beat null model, return empty list.
  # Word doc will have blank page.
  if(which.max(on_fit$devPct) == 1 | max(on_fit$devPct) < 0){
    warning(paste("ordinal net for", outcome_var, "found no models beat null model."))
    return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
                devPct = NA, selected_predictors = character(0), 
                model_size = 0, formula = paste(outcome_var, "~ 1"), 
                significant_df = data.frame(), misclass = NA,
                note = "No ordinal model beat null model."))
  }
  
  ########## 4. Return list of ordinalNet best model details.
  
  ##### 4.1. Refit ordinal regression of best model.

  
  # find highest percentage of deviance explained
  #bestdev <- which(on_fit$devPct == max(on_fit$devPct), arr.ind=TRUE)[1]
  bestbic <- which.min(on_fit$fit$bic)
  coef_mat <- on_fit$fit$coefs[bestbic,]
  
  # remove intercepts
  coef_mat <- coef_mat[(on_fit$fit$nLev):length(coef_mat)]
  
  # get non-0 coefficients
  selected <- unique(predictors[coef_mat != 0])
  
  # Make sure Environment is controlled as a fixed effect.
  if (!"Environment" %in% selected) selected <- c(selected, "Environment")
  
  # Create formula, ensuring outcome variable treated as ordinal and
  # with the random effect of participant ID accounted for.
  # Ensures no overfitting due to data duplication in long format.
  form <- paste("as.ordered(", outcome_var, ") ~", paste(selected, collapse = " + "), "+ (1|ID)")
  
  # Try to fit ordinal regression.
  om <- tryCatch({
    clmm(as.formula(form), data.on, Hess = TRUE)
  }, error = function(e) {
    # If it fails, print error and return empty list.
    # Word doc will have a blank page.
    warning(paste("Ordinal regression refit failed:", e$message))
    return(NULL)
  }) 
  
  ##### 4.2. Return model details as list.
  
  # Get details of model.
  summ <- summary(om)$coefficients
  
  
  ##### 4.2.1. Remove intercepts
  
  if (!"Environment" %in% predictors) predictors <- c(predictors, "Environment")
  
  summ <- remove_intercepts(summ,predictors)
  
  # Format results for table output.
  summ <- data.frame(
    term = row.names(summ),
    estimate = summ[,1],
    p.value = summ[,4],
    row.names = NULL
  )
  
  # Return a list with the details of the model:
  ##   method: clarify ordinalNet, not linear LASSO or best-subset;
  ##   exclude_outcomes: clarify whether outcomes were excluded or not;
  ##   devPct: performance of best model by %age deviance explained
  ##   selected_predictors: which predictors were included in the best model
  ##   model_size: number of predictors in the best model
  ##   formula: formula of best model
  ##   significant_df: predictor effects and p-values
  ##   misclass: cross-validation by misclassification.
  return(list(method = "OrdinalNetCV", exclude_outcomes = exclude_outcomes, 
              devPct = which.max(on_fit$devPct), selected_predictors = selected, 
              model_size = length(selected), formula = form, significant_df = summ,
              misclass = on_fit$misclass[which.max(on_fit$devPct)]))
}




##############################################################################
# Automation                                                                 #
##############################################################################

# A list that will hold results to be exported.
all_results <- list()

# Configuration for best-subset selection
nvmax_limit <- 15  # Higher = slower to compute. >30 incalculable.
method_to_use <- "forward"  # "exhaustive" is impossible with >30 predictors.

# Get the maximum number of models generated.
# We will iterate through each outcome in each dataset twice,
##   once with other outcomes as predictors and once without. This is entirely
##   for progress tracking and is not used in calculations.
# There are 6 outcomes and 3 datasets, so 6*3*2 = 36.
total_combinations <- length(datasets) * length(outcomes) * 2
# A counter of generated models, again for progress tracking.
counter <- 0

########## 1. Main automation loop.
# For each dataset, iterates through outcomes, and for each one,
##   calls get_best_subsets twice: with and without outcomes as predictors,
##   and calls get_lasso_results twice: with and without outcomes as predictors.
# Stores results in list all_results.
for (ds_name in names(datasets)) {
  
  # Get the dataset from iterating through names.
  ds <- datasets[[ds_name]]
  cat("\nProcessing", ds_name, "...\n") # Progress tracking.
  
  # Nested loop, iterating through outcomes for each environment.
  # Calls get_best_subsets and LASSO twice per outcome, once with and once without
  ##   other outcomes as predictors.
  for (outcome in outcomes) {
    # Make sure outcome is in dataset.
    if (!(outcome %in% colnames(ds))) {
      # If not, skip this iteration and continue to next outcome.
      warning(paste("Outcome", outcome, "not found in", ds_name))
      next
    }
    
    # Progress update.
    counter <- counter + 1
    cat(sprintf("[%d/%d] %s | %s (with outcomes)...\n", counter, total_combinations, ds_name, outcome))
    
    
    ##### 1.1. Lasso with predictors as outcomes.
    
    # Calls get_lasso_results (see description above)
    # Args: current dataset, current outcome, penalty (alpha) at industry standard,
    ##   exclude_outcomes false for first run, and outcome list defined globally.
    lasso_w <- get_lasso_results(ds, outcome, alpha.lasso = 1,
                                 exclude_outcomes = FALSE,
                                 nfolds = 5,
                                 all_outcomes = outcomes)
    
    # Add results to export list with key of current dataset, current outcome,
    ##   method lasso, and outcomes included in predictors
    all_results[[paste(ds_name, outcome, "LASSO", "WITH", sep = "_")]] <- lasso_w

    
    ##### 1.2. Lasso without predictors as outcomes.
    
    # Calls get_lasso_results (see description above)
    # Args: current dataset, current outcome, penalty (alpha) at industry standard,
    ##   exclude_outcomes true for this run, and outcome list defined globally.
    lasso_wo <- get_lasso_results(ds, outcome, alpha.lasso = 1,
                                  exclude_outcomes = TRUE,
                                  nfolds = 5,
                                  all_outcomes = outcomes)
    
    # Add results to export list with key of current dataset, current outcome,
    ##   method lasso, and outcomes excluded from predictors
    all_results[[paste(ds_name, outcome, "LASSO", "WITHOUT", sep = "_")]] <- lasso_wo
    
    
    ##### 1.2. Ordinal net with predictors as outcomes.

    for(o in outcomes) {
      data.long[[o]] <- as.ordered(data.long[[o]])
    }

    # Alpha = 1 for LASSO.
    on_with <- get_ordinal_net(data.long, outcome, exclude_outcomes = FALSE,
                               folds.on = 5, all_outcomes = outcomes,
                               alpha.on = 1)

    all_results[[paste(ds_name, outcome, "ON", "WITH", sep = "_")]] <- on_with

    ##### 1.4. Ordinal net without predictors as outcomes.

    # Alpha = 1 for LASSO.
    on_wo <- get_ordinal_net(data.long, outcome, exclude_outcomes = TRUE,
                             folds.on = 5, alpha.on = 1, all_outcomes = outcomes)

    all_results[[paste(ds_name, outcome, "ON", "WITHOUT", sep = "_")]] <- on_wo
   }
}

# Progress update.
cat("\n========================================\n")
cat("Analysis complete!\n")
cat("Total results generated:", length(all_results), "\n")


# Run the unified function
export_results(all_results)

# To make the tables readable, in Word press Alt+F11
# Insert -> Module
# Copy paste the below, then F5 to run.
# Sub SetTableToFitToWindow()
#     Dim t As Table
#     For Each t In ActiveDocument.Tables
#         t.AutoFitBehavior wdAutoFitWindow
#     Next t
# End Sub
