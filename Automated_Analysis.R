##### AUTO-ANALYSER
##### Author: Ron Bar-Ad
##### Last Update: 05/05/2026
##### Description:
# This is a script that reads a long-pivoted CSV and finds significant predictors
# for each outcome. It is designed for ordinal regression (Likert scale responses)
# with a random effect of participant ID, assuming data duplication in the
# pivot process. The variable that was pivoted is called Environment, and
# is included as a fixed effect predictor in all regression models to control
# its effect.
##### Note on use:
# To adapt this script for a different long-pivoted dataset, search "Environment"
# and replace all with the name of the pivot variable. Then change the outcomes
# list, and the rest should work as-is.
##### Process:
# 1. Fit all possible univariate regressions. Keep only significant predictors
#     for each outcome.
# 2. Randomly pair up significant predictors, those that lose significance are
#     removed.
# 3. Remaining predictors fitted to multivariate regressions.
# 4. Significant results reported.


##############################################################################
# Packages                                                                   #
##############################################################################

library(tidyverse) # for reading csv and data handling
library(ordinal) # for mixed effect ordinal regression clmm
library(flextable) # for presentation.

##############################################################################
# Define Data                                                                #
##############################################################################

# Read in csv
#data.use <- read_csv("C:\\Data\\Long.csv")

# Define outcome variables
outcomes <- c("Greenness", "Beauty", "Density", "Safety", "BuildingHeight", "RoadWidth")

# Define predictors as anything that isn't an outcome
#predictors <- setdiff(colnames(data.use), outcomes)
all_predictors <- colnames(data.use)

# Remove ID and Environment from predictors, they will be controlled effects.
all_predictors <- all_predictors[all_predictors != "ID"]
preds_no_env <- all_predictors[!all_predictors %in% c("ID", "Environment")]


##############################################################################
# Functions                                                                  #
##############################################################################

# Function remove_intercepts
##  Removes intercept values from clmm object after ordinal regression fitted,
##  so that only values for named predictors remain.
##  Args:
###   coefs: matrix, the coefficients matrix of a clmm summary object,
####     contains predictor names as row names, effect estimates and p-values.
###   preds: vector (str), list of predictors, which row names are compared to.
remove_intercepts <- function(coefs, preds = all_predictors) {
  # Define iterator.
  iterator <- 1
  
  # Iterate through rows of model details.
  # Each row name is a variable,
  # but some have factor levels appended to the name.
  # For each name, check if it starts with any predictor,
  # or if it already is a predictor,
  # and if not, remove it from the list.
  # While loop runs as long as there are Intercepts remaining.
  while (iterator <= nrow(coefs)) {
    
    # Set removal condition.
    flag <- FALSE
    
    # Cycle through predictors, checking if row name
    # is in the list / starts with a predictor name (is a factor level)
    for (j in preds) {
      
      # startsWith will return TRUE for two equal strings
      # so it also works as if ([rowname] %in% predictors)
      if (startsWith(row.names(coefs)[iterator], j)){
        
        # If it's a predictor, we signal not to remove it,
        # then exit the for loop.
        flag <- TRUE
        next
      }
    }
    # If it wasn't a variable in the predictors list, remove it.
    if (!flag) {
      
      # Remove the row. Treating coefs as data.frame for cases where
      # only 1 column is non-intercept, which would remove the
      # rownames completely and flatten to vector if it were a matrix
      coefs <- as.data.frame(coefs[-c(iterator),])
      
      # The next row in the matrix is now at the index this one was in,
      # so we move our iterator back to keep it where it is.
      iterator <- iterator - 1
    }
    
    # Move iterator on to next row.
    iterator <- iterator + 1
  }
  return(coefs)
}

# Function format_p
##  Rounds p-values within a threshold so that anything under 0.01 is replaced
##  with "< 0.01" and anything else is rounded to three significant figures,
##  used for handling small numbers in scientific notation.
##  Args:
###   p.val: num, a given p-value.
###   threshold: num, the minimum value to represent numerically.
format_p <- function(p.val, threshold = 0.01) {
  # Simple if-else statement: if value is under the threshold:
  ifelse(p.val < threshold, 
         # Return a string "<" whatever the threshold is.
         return(sprintf("< %.2f", threshold)),
         # Otherwise just round to the nearest 3 sig figs.
         return(sprintf("%.3f", p.val))) 
}

# Function map_to_original
##  Takes a list of coefficient names with their factor levels appended and
##  returns the original predictor name. E.g. "StartLoc_X-14" -> "StartLoc_X"
##  Args:
###   exp_name: str, coefficient name to be converted.
###   original_names: vector (str), list of predictors to compare to.
map_to_original <- function(exp_name, original_names = all_predictors) {
  
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

# Function pairwise
##  Iterates through a list of significant predictors and pairs them off
##  into 2-predictor ordinal regressions, removing any that are no longer
##  significant when paired up. Does this until it has included every
##  predictor in at least one pair, then returns remaining predictors.
##  Args:
###   out_var: str, the name of the outcome variable being predicted.
###   preds: vector (str), the names of predictors to be tested.
###   pr: vector (str), the names of all possible predictors (for calling map_to_original)
pairwise <- function(out_var, preds, pr = all_predictors) {
  
  # Start our loop at preds[2] so we can pair it up with preds[1]
  innerLoop <- 2
  
  # Continue looping until every predictor in the list has been paired up
  # at least once. E.g. if preds[1] stays on and preds[2] is removed, the
  # next predictor will move from preds[3] to preds[2], and innerLoop won't
  # increase. When innerLoop = the number of remaining predictors, run one last
  # time, then end loop.
  while (innerLoop <= length(preds)) {
    
    # Progress report expressed as %age of loops to run.
    cat(paste(out_var,":", round((innerLoop/length(preds))*100,3),"%\n"))
    
    # Create an ordinal regression formula using two predictors.
    form <- paste("as.ordered(",out_var,") ~",preds[innerLoop-1], "+", preds[innerLoop])
    
    # If Environment is one of the predictors, it doesn't need to be
    # forced in as a controlled fixed effect.
    if (! "Environment" %in% c(preds[innerLoop-1], preds[innerLoop])){
      # If not in, force in Environment. ID included as random effect.
      form <- paste(form, "+ Environment + (1 | ID)")
    } else {
      # If Environment is in, only include ID.
      form <- paste(form, "+ (1 | ID)")
    }
    
    # Fit ordinal regression with mixed effects.
    fit <- clmm(as.formula(form), data.use)
    
    # Make sure environment included in predictors list,
    # for map_to_original to work.
    if (! "Environment" %in% pr) pr <- c(pr, "Environment")
    
    # Get the coefficients from the model without the intercepts.
    summ <- remove_intercepts(summary(fit)$coefficients, pr)
    
    # If there's any NA values, model was overfitted. We skip over it
    # and handle it later.
    if (!any(is.na(summ[,4]))) {
      
      # If there are no significant results, we can remove both predictors
      if (!any(summ[,4] < 0.05)) {
        # List of predictors is now 2 shorter.
        preds <- preds[-c(innerLoop-1,innerLoop)]
        # We iterate backwards 1 step to stay where we are, the next pair
        # will be in the same spot.
        innerLoop <- innerLoop - 1
      
      # If there are significant results, but not all results are significant,
      # need to figure out which predictors to remove.
      } else if (!all(summ[,4] < 0.05)) {
        
        # Make a list of the non-factor-level names of non-significant predictors
        to_remove <- sapply(row.names(summ[!summ[,4] < 0.05,]), function(x) {
          map_to_original(x, pr)})
        # And remove duplicates and headers.
        to_remove <- unique(to_remove)
        
        # Make a list of non-factor-level names of significant predictors
        to_keep <- sapply(row.names(summ[summ[,4] < 0.05,]), function(x) {
          map_to_original(x,pr)})
        # And remove duplicates and headers.
        to_keep <- unique(to_keep)
        
        # We make sure not to remove any predictors that have a factor level
        # worth keeping. e.g. if StartLoc_X-14 is non significant, but
        # StartLoc_X-188 is significant, we don't want to remove StartLoc_X,
        # so filter list of removal predictors with any worth keeping.
        to_remove <- to_remove[!to_remove %in% to_keep]
        
        # Environment is a controlled predictor no matter what, if it needs
        # to be removed it will still be added in the next iteration, and the
        # loop will never end. Ignoring Environment, if there are still
        # predictors to remove, remove them and iterate a step backwards.
        if (length(to_remove[to_remove != "Environment"]) > 0) {
          
          # Filter predictors to only ones not in the removal list.
          preds <- preds[!preds %in% to_remove]
          
          # Iterate one step backwards, this will cause the loop to skip
          # predictors if to_remove has > 2 names, but this function will
          # be called again in a randomised order several times to account
          # for this.
          innerLoop <- innerLoop - 1
        }
      }
    }
    
    # Iterate forwards one step at the end of each loop, progressing through
    # the list of significant predictors.
    innerLoop <- innerLoop + 1
  }
  
  # Return the reduced list once all loops are done.
  return(preds)
}


fit_clmm <- function(var, preds, df){
  # Make sure environment is not included, to account for it as fixed effect.
  if ("Environment" %in% preds) preds <- preds[preds != "Environment"]
  # Build our mixed effect formula with environment as fixed effect.
  form <- paste("as.ordered(",var,") ~",paste(sig_preds, collapse=" + "), "+ Environment + (1 | ID)")
  # Fit the model to the formula
  fit <- clmm(as.formula(form), df)
  # Return the model coefficients
  return(summary(fit))
}
##############################################################################
# Automation                                                                 #
##############################################################################

##### 1. Univariate regressions.

# We pair every outcome with every predictor, one by one, adding significant
# results to a list for step 2.

# Progress tracking.
start.1 <- Sys.time()
cat("1. ----- Started automation. ------\n")

# Empty list to be filled with significant results.
data.pair <- data.frame(matrix(nrow=0,ncol=5))
colnames(data.pair) <- c("Outcome", "Predictor", "Effect", "P.value", "AIC")

# Outer loop: iterating through the outcomes.
for (i in 1:length(outcomes)) {
  
  # Progress update.
  cat(paste("Generating univariate models,", round((i/length(outcomes)*100),3), "% done.\n"))
  
  # For readability, instead of typing data.use[[outcomes[[i]]]] every time.
  out <- outcomes[i]
  
  # Make sure predictors do not include outcome variable.
  preds_no_env <- setdiff(colnames(data.use), c("ID", "Environment", outcomes))
  all_predictors <- setdiff(colnames(data.use), c("ID", outcomes))
  
  # Make sure our outcome variable is treated as ordinal.
  data.use[[out]] <- as.ordered(data.use[[out]])
  
  # Inner loop: iterating through the predictors per outcome.
  for (pred.j in 1:length(preds_no_env)) {
    
    # For readability.
    pred <- preds_no_env[pred.j]
    
    # Skip if the predictor being estimated is Environment (always fixed),
    # ID (always random effect) or somehow the predictor.
    if (pred %in% c("ID", "Environment", out)) next
    
    # Put together a formula with the outcome and predictor, accounting
    # for the fixed effect of Environment and the random effect of ID.
    form <- paste(out, "~", pred, "+ Environment + (1 | ID)")
    
    # Fit a mixed effects ordinal regression to the formula.
    model <- clmm(formula = form, data = data.use)
    
    # Check we actually have results, otherwise skip to next.
    if (!is.null(summary(model)$coefficients) & !any(is.na(summary(model)$coefficients[,4]))) {
      
      # Removing intercepts for just predictors.
      sigs <- remove_intercepts(summary(model)$coefficients, all_predictors)
      
      # If there's anything other than intercepts, we keep the predictor
      # names, effect estimates, and p.values.
      if (nrow(sigs) > 0){
        
        # Only keep significant findings.
        sigs <- sigs[sigs[,4] < 0.05,]
        
        # If there are any, save them as results.
        if (nrow(sigs[!startsWith(row.names(sigs), "Environment"),]) > 0) {
          for (row in row.names(sigs)){
            data.pair[nrow(data.pair)+1,] <- c(out, row, round(sigs[row,1], 2),
                                               format_p(sigs[row,4]), 
                                               model$info["AIC"])
          }
        }
      }
    }
  }
}

# Progress update.
end <- Sys.time()
cat(paste("1. ------ Univariate models complete. Time taken in seconds:", round(end-start.1,3), "-----\n"))


##### 2. Pairwise and multivariate regression.

# We'll loop through each outcome, for each of its significant predictors,
# and pair them up at random a few times, leaving only the predictors that
# remain significant when not univariate. Then, those will be fitted together
# for a final model for each outcome.

# Progress tracking.
start <- Sys.time()
cat("2. ----- Fitting pairwise regressions. -----\n")

# Fresh new data frame for the final results.
sig <- data.frame(matrix(nrow=0,ncol = 5))
# Give it the 4 column names we'll be reporting.
colnames(sig) <- colnames(data.pair)

all_predictors <- setdiff(colnames(data.use),"ID")
preds_no_env <- setdiff(colnames(data.use), c("ID", "Environment"))

# Outer loop: for each outcome, we fit as many pairwise models as there are
# variable pairs, repeating the process up to 5 times or until there are < 3
# predictors in the model.
for (var in outcomes) {
  # Progress update.
  cat(paste("Fitting pairwise models", round((which(outcomes == var)/length(outcomes))*100,3),
            "% done.\n"))
  
  # Get the predictors for this outcome as a list of original predictor names.
  sig_preds <- sapply(data.pair$predictor[data.pair$outcome == var], function(x){
    map_to_original(x,all_predictors)})
  # Get rid of duplicates and headers.
  sig_preds <- unique(sig_preds)
  
  # Loop through a maximum of five times. This may be computationally heavy,
  # each loop contains many smaller loops. E.g. if an outcome has 40 significant
  # univariate predictors, and 5 of them remain significant when paired up,
  # the first outer loop will have 40 inner loops, and the second, third, 
  # fourth, and fifth will have 5 inner loops, for a total of 60 inner loops
  # for just one outcome.
  loop <- 5
  
  # We also stop if we have fewer than 3 significant predictors, that's enough.
  while (length(sig_preds) > 3 & loop > 0) {
    # Randomise the list of predictors so you don't end up pairing off the same
    # two significant predictors over and over. This accounts for the potential
    # to slip through the cracks in pairwise (see comments above).
    sig_preds <- sample(sig_preds)
    # Run the two-predictor regressions.
    sig_preds <- pairwise(var, sig_preds, all_predictors)
    # Keep the loop moving.
    loop <- loop-1
  }
  
  # Progress update.
  cat(paste("Fitting multivariate model for", var,".\n"))
  
  # Make sure environment is not included, to account for it as fixed effect.
  if ("Environment" %in% sig_preds) sig_preds <- sig_preds[sig_preds != "Environment"]
  
  # Fit the final multivariate mixed effect ordinal regression for the outcome.
  # Get rid of the intercepts so we only have predictors to report.
  
  summ <- fit_clmm(var, sig_preds, data.use)
  aic <- summ$info["AIC"]
  summ <- remove_intercepts(summ$coefficients, all_predictors)
  
  tries <- length(sig_preds)
  while (any(is.na(summ[,4])) & tries > 0) {
    summ <- fit_clmm(var, sig_preds[-c(tries)], data.use)
    aic <- summ$info["AIC"]
    summ <- remove_intercepts(summ$coefficients, all_predictors)
    tries <- tries - 1
  }
  
  if (!any(is.na(summ[,4]))) {
  
    # Only keep significant predictors in the new model.
    summ <- summ[summ[,4] < 0.05,]
    
    # If there are none, there's no point adding anything to the new data frame,
    # and this last step gets skipped.
    if (nrow(summ) > 0) {
      
      # For each coefficient we'll report, add a new row to the data frame.
      for (co in row.names(summ)) {
        # We add a row with the outcome variable, the predictor, the effect, and
        # the p-value.
        sig[nrow(sig)+1,] <- c(var, co, round(summ[co,1],2),format_p(summ[co,4]),aic)
      }
    }
  }
}

# Final progress report.
end <- Sys.time()
cat(paste("2. ----- Pairwise and multivariate finished. Time taken in seconds:", round(end-start,3), "-----\n"))
cat(paste("All finished! Total time taken was: Time taken in seconds:", round(end-start.1,3), ".\n"))

# Display final results as flextable.
flextable(sig) %>% theme_vanilla() %>%
  set_caption("Significant predictors per outcome.")
