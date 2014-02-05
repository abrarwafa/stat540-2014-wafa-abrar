prDat <- read.table("GSE4051_MINI.txt", row.names = 1)
str(prDat)


# How many rows are there? Hint: nrow(), dim().
# 39
nrow(prDat)
dim(prDat)

# How many columns or variables are there? Hint: ncol(), length(), dim().
# 6
ncol(prDat)
length(prDat)

# Inspect the first few observations or the last few or a random sample. 
# Hint: head(), tail(), x[i, j] combined with sample().
head(prDat)
# sample devStage gType crabHammer eggBomb poisonFang
# Sample_20     20      E16    wt     10.220   7.462      7.370
# Sample_21     21      E16    wt     10.020   6.890      7.177
# Sample_22     22      E16    wt      9.642   6.720      7.350
# Sample_23     23      E16    wt      9.652   6.529      7.040
# Sample_16     16      E16 NrlKO      8.583   6.470      7.494
# Sample_17     17      E16 NrlKO     10.140   7.065      7.005
tail(prDat)
# sample devStage gType crabHammer eggBomb poisonFang
# Sample_38     38  4_weeks    wt      9.767   6.608      7.329
# Sample_39     39  4_weeks    wt     10.200   7.003      7.320
# Sample_11     11  4_weeks NrlKO      9.677   7.204      6.981
# Sample_12     12  4_weeks NrlKO      9.129   7.165      7.350
# Sample_2       2  4_weeks NrlKO      9.744   7.107      7.075
# Sample_9       9  4_weeks NrlKO      9.822   6.558      7.043
prDat[sample(nrow(prDat), size = 6), ]

# What does row correspond to -- different genes or different mice?
# Different mice

# What are the variable names? Hint: names(), dimnames().
names(prDat)
#  "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"    "poisonFang"
dimnames(prDat)
# [[1]]
# [1] "Sample_20" "Sample_21" "Sample_22" "Sample_23" "Sample_16" "Sample_17"
# [7] "Sample_6"  "Sample_24" "Sample_25" "Sample_26" "Sample_27" "Sample_14"
# [13] "Sample_3"  "Sample_5"  "Sample_8"  "Sample_28" "Sample_29" "Sample_30"
# [19] "Sample_31" "Sample_1"  "Sample_10" "Sample_4"  "Sample_7"  "Sample_32"
# [25] "Sample_33" "Sample_34" "Sample_35" "Sample_13" "Sample_15" "Sample_18"
# [31] "Sample_19" "Sample_36" "Sample_37" "Sample_38" "Sample_39" "Sample_11"
# [37] "Sample_12" "Sample_2"  "Sample_9" 
# 
# [[2]]
# [1] "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"    "poisonFang"

# What "flavor" is each variable, i.e. numeric, character, factor? Hint: str().
str(prDat)
#'data.frame':  39 obs. of  6 variables:
# $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
# $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
# $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
# $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
# $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
# $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...

# For sample, do a sanity check that each integer between 1 and the number of rows in the dataset occurs exactly once. Hint: a:b, seq(), seq_len(), sort(), table(), ==, all(), all.equal(), identical().
identical(sort(prDat$sample), seq_len(nrow(prDat))) #TRUE

# For each factor variable, what are the levels? Hint: levels(), str().
levels(prDat$devStage)
# "4_weeks" "E16"     "P10"     "P2"      "P6" 
str(prDat$devStage)
# Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
levels(prDat$gType)
# "NrlKO" "wt"
str(prDat$gType)
# Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...

# How many observations do we have for each level of devStage? For gType? Hint: summary(), table().
summary(prDat$devStage)
# 4_weeks     E16     P10      P2      P6 
# 8       7       8       8       8 
table(prDat$gType)
# NrlKO    wt 
# 19    20 

# Perform a cross-tabulation of devStage and gType. Hint: table().
with(prDat, table(devStage, gType))
#or 
table(c(prDat["devStage"], prDat["gType"]))
#          gType
# devStage  NrlKO wt
# 4_weeks     4  4
# E16         3  4
# P10         4  4
# P2          4  4
# P6          4  4
# If you had to take a wild guess, what do you think the intended experimental design was? What actually happened in real life?
# 
# For each quantitative variable, what are the extremes? How about average or median?
quantile(prDat$sample)
# 0%  25%  50%  75% 100% 
# 1.0 10.5 20.0 29.5 39.0 
summary(prDat$crabHammer)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.214   8.938   9.611   9.428   9.830  10.340 
quantile(prDat$eggBomb)
# 0%    25%    50%    75%   100% 
# 6.1380 6.2780 6.7570 7.0945 8.1730 
summary(prDat$poisonFang)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.735   7.188   7.350   7.379   7.476   8.584 

# Create a new data.frame called weeDat only containing observations for 
# which expression of poisonFang is above 7.5.
weeDat <- subset(prDat, poisonFang>7.5)
weeDat

# For how many observations poisonFang > 7.5? 
# How do they break down by genotype and developmental stage?
str(weeDat) # 9 observations
summary(weeDat$devStage)
with(weeDat, table(devStage, gType))

# Print the observations with row names "Sample_16" and "Sample_38" to screen, 
# showing only the 3 gene expression variables.
prDat[c ("Sample_16", "Sample_38"), c("crabHammer", "eggBomb", "poisonFang") ]

# Which samples have expression of eggBomb less than the 0.10 quantile?
prDat[prDat$eggBomb < quantile(prDat$eggBomb, 0.1), ]

