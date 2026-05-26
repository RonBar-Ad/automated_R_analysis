# BAR-AD VR EXPERIMENT data preparation
# written by Ron Bar-Ad in 2025.
# see online repository at https://drive.proton.me/urls/DMR6T0XPR4#Fl0xwnhCAgIm

# Data prep does the following:
# 1. merges data from in-VR questionnaires and qualtrics
# 2. standardises qualtrics response values
# 3. gives columns reasonable names
# 4. exports to !MasterData.csv

# Note: due to the questionnaire's logic design, questions that did not apply
# to a respondent were not shown, and would therefore be NA. All NAs can 
# therefore simply be replaced by 0 or equivalent "no" in categorical variables.


library(tidyverse) # for dataframe manipulation
library(mice) # for multivariate imputation (used for ConceptBuild)

# Read in data from files.
path.t <- "C:\\BarAdVR\\Data\\" # change path to file location.
data.n <- read_csv(paste(path.t, "Exp2_VRQAnswers.csv",sep=""))
data.i <- read_csv(paste(path.t, "Exp2_Qualtrics_Numerical.csv",sep=""))

# Naming columns for easy identification.
colnames(data.i) <- c("Gender", "Age", "ExpWorkYN", "ExpWorkOther",
                      "ExpWorkPlanAndOther", "ExpWorkPlanYears", "ExpWorkOtherYears",
                      "ExpStudYes", "ExpStudYesOther", "ExpStudNo", "ExpStudUnsure",
                      "ExpStudOther", "ExpStudUniYN", "ExpStudUniOther",
                      "ExpStudCurrentYN", "ExpStudCurrentOther", "ExpStudHobbyYN",
                      "LivedExpEarlyUrb", "LivedExpEarlySub", "LivedExpEarlyRur",
                      "LivedExpEarlyNon", "LivedExpEarlyUnsure", "LivedExpLateUrb",
                      "LivedExpLateSub", "LivedExpLateRur", "LivedExpLateNon",
                      "LivedExpLateUnsure", "LivedExpLife", "ConceptPerc",
                      "ConceptGreen", "ConceptBuild", "ConceptIVE", "UseFreq",
                      "UseComplex", "UseEasy", "UseTech", "UseInt", "UseInc", "UseQuick",
                      "UseCumb", "UseConf", "UseLearn", "TaskMental", "TaskPhysical",
                      "TaskRush", "TaskAccom", "TaskHard", "TaskStress",
                      "QualProfessional", "QualResearch", "QualGeneral")

# merge with qualtrics table for 1 master dataframe
data.i <- rowid_to_column(data.i, "ID")
data <- merge(data.i, data.n, by="ID")

# 0 = No work exp, 1 = other work exp, 2 = UrbPlan work exp, 3 = UrbPlan and other work exp
data$ExpWorkYN[data$ExpWorkYN == 2] <- 0
data$ExpWorkYN[data$ExpWorkYN == 5] <- 3
data$ExpWorkYN[data$ExpWorkYN == 1] <- 2
data$ExpWorkYN[data$ExpWorkYN == 6] <- 1

# Encoding the multiple choice so that 1 is yes and 0 is no
data["ExpWorkPlanYears"][is.na(data["ExpWorkPlanYears"])] <- 0
data["ExpWorkOtherYears"][is.na(data["ExpWorkOtherYears"])] <- 0
data["ExpStudYes"][is.na(data["ExpStudYes"])] <- 0
data["ExpStudYesOther"][is.na(data["ExpStudYesOther"])] <- 0
data["ExpStudNo"][is.na(data["ExpStudNo"])] <- 0
data["ExpStudUnsure"][is.na(data["ExpStudUnsure"])] <- 0

# 0 = no (including NA), 1 = college, 2 = uni
data$ExpStudUniYN[is.na(data$ExpStudUniYN)] <- 0
data$ExpStudUniYN[data$ExpStudUniYN == 4] <- 0
data$ExpStudUniYN[data$ExpStudUniYN == 1] <- 3
data$ExpStudUniYN[data$ExpStudUniYN == 2] <- 1
data$ExpStudUniYN[data$ExpStudUniYN == 3] <- 2


# 0 = no (i.e. na), 1 is college, 2 is undergrad, 3 is postgrad
data$ExpStudCurrentYN[is.na(data$ExpStudCurrentYN)] <- 0
data$ExpStudCurrentYN[data$ExpStudCurrentYN == 5] <- 0

data$ExpStudHobbyYN[data$ExpStudHobbyYN == 1] <- 0
data$ExpStudHobbyYN[data$ExpStudHobbyYN == 2] <- 1

# any missing answers for lived exp become 0
for(i in 21:31){
  data[,i][is.na(data[,i])] <- 0
}

# livedexplife: 1 = urb, 2= sub, 3 = rur, 4 = urb & sub
# concepts: 1 = unaware, 2 = awareness, 3 = understanding, 4 = strong understanding, 5 = unsure

# Usability for some reason encoded to skip 2 in qualtrics ¯\_(ツ)_/¯
for (i in 36:45){
  data[,i][data[,i] == 3] <- 2
  data[,i][data[,i] == 4] <- 3
  data[,i][data[,i] == 5] <- 4
  data[,i][data[,i] == 6] <- 5
  
  # imputed data from means for first two participants
  data[,i][is.na(data[,i])] <- round(mean(as.numeric(data[,i]), na.rm=TRUE),digits = 0)
}

# Task eval:
# 0 = strongly disagree, then 1-4, then 5 = strongly agree
for(i in 46:51){
  data[,i][data[,i] == 1] <- 0
  data[,i][data[,i] == 2] <- 1
  data[,i][data[,i] == 3] <- 2
  data[,i][data[,i] == 4] <- 3
  data[,i][data[,i] == 5] <- 4
  data[,i][data[,i] == 6] <- 5
  data[,i][data[,i] == 7] <- 6
  
  # imputed data from means for first two participants
  data[,i][is.na(data[,i])] <- round(mean(as.numeric(data[,i]), na.rm=TRUE),digits = 0)
}

# Re-encoding start location as 1 column, where 1 = less dense area, 2 = centre,
# and 3 = more dense area
data$StartLoc[data$StartLoc_X == -327] <- 1
data$StartLoc[data$StartLoc_X == -41] <- 2
data$StartLoc[data$StartLoc_X == 188] <- 3

# Removing the unnecessary start location variables.
data <- data[,! colnames(data) %in% c("StartLoc_X", "StartLoc_Z")]

# The concept variables are 1-4 ordinal, but have an "Unsure / Prefer not to say" option
# at level 5. Luckily for these variables, there's only 1 instance of this between
# all of them. We can consider this missing data and impute it based on patterns
# of other participants with similar answers to other questions.
# Excluding non-background variables from pattern.
data.tmp <- data[,2:23]
data.tmp$ConceptBuild <- as.numeric(data.tmp$ConceptBuild)
data.tmp$ConceptBuild[data.tmp$ConceptBuild == 5] <- NA
data.mice <- mice(data.tmp, m=1, method = "pmm", seed = 20070830)
data.tmp <- complete(data.mice)
# Whichever one is different in the imputed dataset is brought in to replace the
# non-ordinal value.
data$ConceptBuild[data$ConceptBuild == 5] <- data.tmp$ConceptBuild[data.tmp$ConceptBuild != data$ConceptBuild]

# Changing the order slightly for readability.
data <- data[,c(1:40, 59, 41:58)]


# And export to MasterData file.
write_csv(data, paste(path.t, "!MasterData.csv" , sep=""))

