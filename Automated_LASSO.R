# Automated Dimension Reduction Ordinal Regressions
##### Author: Ron Bar-Ad
##### Last Update: 01/05/2026
##### Description:
# This is a script that reads a CSV and fits ordinalNet and glmnet
# cross-validated LASSO-penalty dimension reduction for each outcome on all
# predictors, then verifies by fitting a mixed-effect ordinal regression.
# The results are then exported to a word document with model details.
# It is designed for ordinal regression (Likert scale responses) with a random
# effect of participant ID, assuming data duplication in the pivot process. The
# variable that was pivoted is called Environment, and is included as a fixed
# effect predictor in all regression models to control its effect.
# It is also designed to allow for multiple datasets, for example with and
# without factor level pooling.
##### Note on use:
# To adapt this script for a different long-pivoted dataset, search "Environment"
# and replace all with the name of the pivot variable. Then change the outcomes
# list, and the rest should work as-is.
# To add more than one dataset, simply read the data in, prepare it however it
# needs to be prepared, and add it to the list "datasets" (defined right before
# the utility functions).
##### Process:
# 1. Data preparation: remove factors with < 2 levels of meaningful frequency,
#     pool factor levels with neighbours where they fall under the frequency
#     threshold, calculate SUS usability score and TLX task load score.
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
setwd("C:\\Users\\ronba\\Nextcloud\\PhD\\Experiments_2-3\\")

# Import merged data from Qualtrics and in-VR questionnaire
data.all <- read_csv("C:\\Users\\ronba\\Nextcloud\\PhD\\Experiments_2-3\\Data\\!MasterData.csv")


# Define our outcome variables, this will come in handy when running our Lasso
outcomes <- c("Greenness", "Beauty", "Density", "Safety", "BuildingHeight", "RoadWidth")
#outcomes <- c("UseScore", "TLXScore")



# Define our continuous variables (do not mutate to factor)
not_factor <- c("Age", "ExpWorkPlanYears", "ExpWorkOtherYears")

# Define our wide-format outcome variables (so they aren't included in
# factor level merging).
for (col in colnames(data.all)) {
  for (out in outcomes) {
    if (endsWith(col, out)){
      not_factor <- c(not_factor, col)
    }
  }
}

# Set nominal variables as factors
# data.quant <- data.quant %>% 
#   mutate_at(vars(-c(3,5,6,20:39, 43:60)), as.factor)
data.quant <- data.all %>% mutate_at(vars(-any_of(not_factor)), as.factor)

# Stick ID in there too.
not_factor <- c(not_factor, "ID")


# FACTOR LEVEL MERGING.

# Set the minimum threshold for level frequency. Anything under this
# is merged to avoid model confusion. Agresti (2007) suggests a threshold of 5,
# Kuhn and Johnson (2019) suggest a 95:5 ratio. Since 5% of the sample in
# this dataset is lower than 5, 5 is the default threshold.
threshold <- 5
# Empty list to be filled with columns to remove entirely. Any factor with
# only 2 levels where one level has frequency < threshold is removed
# or it will cause spurious correlations with no meaning.
remove <- c()

# Iterating through our columns, testing one by one if they should
# be removed, pooled, or left as they are.
for (col in colnames(data.quant)) {
  
  # Leave the non-factors (and outcome vars) alone
  if (!col %in% not_factor) {
    
    # For readability.
    tab <- table(data.quant[[col]])
    
    # If there's only one factor level, it's a useless variable,
    # and can be removed.
    if (length(tab) < 2) remove <- c(remove, col)
    
    # If there are only 2 levels and one is below the threshold, it will
    # skew all analysis, but pooling would make it useless. It is therefore
    # removed.
    if (any(tab < threshold) & length(tab) < 3) {
      remove <- c(remove, col)
    } else {
      
      # If there are enough levels for pooling,
      # Test if there are any levels with frequencies below the threshold,
      # and pool them to their neighbours.
      
      # Loop count is set to the initial number of levels in the factor to
      # ensure the loop doesn't get stuck at the bottom and never end.
      count <- length(table(data.quant[[col]]))
      
      # Use a While loop so that number of levels is continually updated,
      # where a for loop sets its parameters before it begins. 
      while (any(table(data.quant[[col]]) < threshold & count > 0)) {
        
        # If we've pooled enough levels in past loops 
        # that there's only 2 levels left and one of them is *still* < 5,
        # this won't be a useful variable and we can remove it.
        if (length(table(data.quant[[col]])) < 3) {
          
          remove <- c(remove, col)
          # move on to next column
          next
        
        # If there are levels to pool, we get to pooling.
        } else {
          
          # Only bother looping through 1 at a time if there's more than 1
          # level that's < threshold.
          if (length(which(table(data.quant[[col]]) < threshold)) > 1) {
            
            # Start looping from last level.
            level <- length(table(data.quant[[col]]))
            
            # Using a while loop to iterate backwards, ensuring
            # that levels are all merged from one direction.
            while (level > 1) {
              
              # For readability
              tab <- table(data.quant[[col]])
              
              # If we're below the threshold, pool with smallest neighbour.
              if (tab[level] < threshold) {
                
                # Return the lowest neighbour's index number, unless we're
                # at the start or end, in which case we simply return the
                # index of the only neighbour.
                if (level == length(tab)) {
                  neighbour <- level-1
                } else if (level <= 1) {
                  neighbour <- level+1
                } else{
                  # Don't ask me why it's a nested if I tried a bunch of other
                  # stuff that didn't work for some reason.
                  # Smallest neighbour's index is returned.
                  if (tab[level-1] < tab[level+1]) neighbour <- level-1 else neighbour <- level+1
                }
                
                # Both level names pasted together to form new level.
                # E.g. if levels were 5 and 4, new level is 4&5.
                to_merge <- names(tab[c(neighbour,level)])
                new_level <- paste(to_merge, collapse="&")
                # Replace both levels with new level.
                levels(data.quant[[col]])[levels(data.quant[[col]]) %in% to_merge] <- new_level
              }
              # Move on to next levels in factor.
              level <- level -1
            }
          # If there's only 1 level < threshold, we just find it and replace it
          # with its smallest neighbour without all the loops and checks.
          } else {
            # For readability
            tab <- table(data.quant[[col]])
            # Set condition, will be true when index is collected.
            found <- FALSE
            # Iterator to go through the factor levels.
            level <- 1
            
            # Look I know this one really didn't need to be a while loop, but
            # I tried a bunch of other ways about this and none of them worked
            # for reasons I did not understand. It's a while loop because that's
            # the only way I got the index number to work.
            while (!found) {
              # If the level is below the threshold, that's our target.
              if (tab[level] < threshold){
                # Flip the condition, end the loop.
                found <- TRUE
              # If it's not our target, iterate onward.
              } else {
                level <- level +1
              }
            }
            
            # Find our smallest neighbour.
            # If it's the first level, its only neighbour is +1.
            if (level <= 1) {
              neighbour <- level + 1
            
            # If it's the last level, its only neighbour is -1.
            } else if (level >= length(tab)) {
              neighbour <- level -1
              
            # If it's somewhere in the middle, we find the smallest of +/- 1.
            } else {
              if (tab[level+1] < tab[level-1]) neighbour <- level+1 else neighbour <- level-1
            }
            # Get the names of the levels.
            to_merge <- names(tab[c(neighbour,level)])
            # Paste them together, so that levels "5" and "4" become level "4&5".
            new_level <- paste(to_merge, collapse="&")
            # Replace old levels with new level.
            levels(data.quant[[col]])[levels(data.quant[[col]]) %in% to_merge] <- new_level
          }
        }
        # Onwards through the table.
        count <- count - 1
      }
    }
  }
}

# We now have a data.quant with pooled factor levels and a list of columns to
# remove. We remove the columns.
data.quant <- data.quant[,!colnames(data.quant) %in% remove]

# Pivot data to long format
##   This will allow the use of Environment as a predictor with 3 levels
##   and each perception (dependent variables) as a single outcome,
##   but will duplicate rows from Qualtrics responses.
##   Shouldn't affect models, since data proportions are  the same,
##   but may interfere in analysis that doesn't involve the in-VR questionnaires.
to_pivot <- c()
for (col in 1:length(data.quant)){
  for (out in outcomes){
    if (endsWith(colnames(data.quant)[col],out)){
      to_pivot <- c(to_pivot, col)
    }
  }
}

data.long <- pivot_longer(data=data.quant, cols=all_of(to_pivot), names_to=c("Environment", ".value"), names_sep="_")

# Set environment factor levels. 0 = no trees, 1 = 315 trees, 2 = 634 trees.
data.long$Environment[data.long$Environment == "NT"] <- "0"
data.long$Environment[data.long$Environment == "ST"] <- "1"
data.long$Environment[data.long$Environment == "MT"] <- "2"
data.long$Environment <- as.factor(data.long$Environment)

##### Defining usability score from SUS.
# divide into positive and negative
pos <- c("UseFreq", "UseEasy", "UseInt", "UseQuick", "UseConf")
neg <- c("UseComplex", "UseTech", "UseInc", "UseCumb", "UseLearn")
# take only SUS columns
data.use <- data.all[, colnames(data.all) %in% c(pos, neg)]
# for positive columns subtract 1
data.use[pos] <- lapply(
  data.use[pos], function(x) ifelse(x>1,x-1,1))
# for negative columns subtract from 5
data.use[neg] <- lapply(
  data.use[neg], function(x) 5-x)
# score is the sum of all columns * 2.5, giving a range up to 100
data.use <- data.use %>% mutate(UseScore = (UseFreq + UseEasy + UseInt + UseQuick +
                                              UseConf + UseComplex + UseTech + UseInc +
                                              UseCumb + UseLearn) * 2.5)
# To prepare for long-pivot, we repeat each UseScore as many times as our pivot
# increased the number of rows.
use_score <- c()
reps <- nrow(data.long) / nrow(data.use)
for (i in 1:nrow(data.use)){
  use_score <- c(use_score, rep(data.use$UseScore[i], reps))
}
##### Defining TLX score from NASA TLX
data.tlx <- data.all[,startsWith(colnames(data.all), "Task")]
data.tlx <- data.tlx %>% mutate(TLXScore = (((((1+TaskMental)+(1+TaskPhysical)+(1+TaskRush)
                                               +(7-TaskAccom)+(1+TaskHard)+(1+TaskStress))/6)-1)/6)*99+1)
# The same pivot repetition for TLX score
tlx_score <- c()
reps <- nrow(data.long) / nrow(data.tlx)
for (i in 1:nrow(data.tlx)) {
  tlx_score <- c(tlx_score, rep(data.tlx$TLXScore[i], reps))
}

# And bring them in.
data.long <- data.long[,!colnames(data.long) %in% c(colnames(data.use), colnames(data.tlx))]
data.long$UseScore <- use_score
data.long$TLXScore <- tlx_score

# The datasets to use. Will automatically loop through all.
datasets <- list(
  "data.long" = data.long
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

# Function remove_intercepts:
##   Handles output object of clmm, removes intercepts from coefficients matrix.
##   Args:
###     sum: matrix, the coefficients matrix of a clmm object summary.
###     preds: vector (str), list of predictors to be kept in coefficients list,
####       to identify intercepts.
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

# Function format_p:
##   Formats p-values from long numerics or scientific notation to either
##   3-digit numbers or a string of "< [threshold]" with whatever the threshold
##   for p-values too small to care about it.
##   Args:
###     p: num, a p-value from some model.
###     threshold: num, the minimum threshold below which p-values are just
####       replaced with a string. Defaults to 0.01 to keep to 3 digits.
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
export_results <- function(results_list, output_file = "Automated_Report.docx") {
  
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
        bic_txt <- fpar(paste("BIC:",
                                 round(res$bic, 5)))
        doc <- body_add_fpar(doc, bic_txt, style = "Normal")
        
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
        # if (is.character(res$lambda)) {
        #   lambtxt <- strsplit(res$lambda,": ")[[1]][1]
        #   lambnum <- as.numeric(strsplit(res$lambda, ": ")[[1]][2])
        # }
        lambda_txt <- fpar(round(res$lambda,3))
        doc <- body_add_fpar(doc, lambda_txt, style = "Normal")
        
        cv_txt <- fpar(paste("Cross-validated error:", round(res$cv_error, 4)))
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
  
  # Get best model predictors.
  coef_opt <- coef(cv_fit, s = "lambda.min")
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  
  
  
  opt_vars <- rownames(coef_opt)[-1][coef_opt[-1, 1] != 0]
  ose_vars <- rownames(coef_1se)[-1][coef_1se[-1, 1] != 0]
  
  # Get names of predictors:
  ##   For factors converted to numeric level codes,
  ##   variable names are the original column names.
  opt_vars <- rownames(coef_opt)[-1][coef_opt[-1, 1] != 0]
  ose_vars <- rownames(coef_1se)[-1][coef_1se[-1, 1] != 0]
  
  if (length(opt_vars) > 5 & length(ose_vars[ose_vars != "Environment"]) > 0 & length(ose_vars) < length(opt_vars)){
    selected_vars <- ose_vars
    lambda_opt <- cv_fit$lambda.1se
  } else {
    selected_vars <- opt_vars
    lambda_opt <- cv_fit$lambda.min
  }
  
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
  ###     maxiterOut: num, max number of outer-loop iterations. set to 2000 to
  ####       ensure precise deviance assessment
  ###     alpha: num, penalty for predictors. Set to 1 for LASSO penalty.
  on_fit <- tryCatch({
    # ordinalNetCV(X, y, nFolds = folds.on, nFoldsCV = folds.on.cv,
    #              maxiterOut = 500, alpha = 1)
    ordinalNetTune(X, y, folds = row_folds, family = "cumulative",
                   link = "logit", penaltyFactors = penalties,
                   parallelTerms = TRUE, nonparallelTerms = FALSE,
                   maxiterOut = 2000)
    
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

  
  # find lowest BIC
  
  model_bic <- sort(on_fit$fit$bic)
  bestbic <- which(on_fit$fit$bic == model_bic[1])
  bic1se <- which(on_fit$fit$bic == model_bic[2])
  
  coef_opt <- on_fit$fit$coefs[bestbic,]
  coef_1se <- on_fit$fit$coefs[bic1se,]
  
  # remove intercepts
  coef_opt <- coef_opt[(on_fit$fit$nLev):length(coef_opt)]
  coef_1se <- coef_1se[(on_fit$fit$nLev):length(coef_1se)]
  
  # get non-0 coefficients
  coef_opt <- unique(predictors[coef_opt != 0])
  coef_1se <- unique(predictors[coef_1se != 0])
  
  if ((length(coef_opt) > 5 & length(coef_1se) > 0 & length(coef_1se) < length(coef_opt)) | length(coef_opt) < 1 & length(coef_1se) > 0) {
    selected <- coef_1se
    misclass <- on_fit$misclass[bic1se]
    bic <- model_bic[2]
  } else {
    selected <- coef_opt
    misclass <- on_fit$misclass[bestbic]
    bic <- model_bic[1]
  }
  
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
              bic = bic, selected_predictors = selected, 
              model_size = length(selected), formula = form, significant_df = summ,
              misclass = misclass))
}




##############################################################################
# Automation                                                                 #
##############################################################################

# A list that will hold results to be exported.
all_results <- list()

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
    
    
    #### 1.2. Ordinal net with predictors as outcomes.

    for(o in outcomes) {
      ds[[o]] <- as.ordered(ds[[o]])
    }

    # Alpha = 1 for LASSO.
    on_with <- get_ordinal_net(ds, outcome, exclude_outcomes = FALSE,
                               folds.on = 5, all_outcomes = outcomes,
                               alpha.on = 1)

    all_results[[paste(ds_name, outcome, "ON", "WITH", sep = "_")]] <- on_with

    ##### 1.4. Ordinal net without predictors as outcomes.

    # Alpha = 1 for LASSO.
    on_wo <- get_ordinal_net(ds, outcome, exclude_outcomes = TRUE,
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
