Seminar 01
========================================================

This is an R Markdown document for Seminar 01 by Abrar Wafa


```r
prDat <- read.table("GSE4051_MINI.txt", row.names = 1)
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```


How many rows are there? Hint: nrow(), dim().

```r
nrow(prDat)
```

```
## [1] 39
```

```r
dim(prDat)
```

```
## [1] 39  6
```


How many columns or variables are there? Hint: ncol(), length(), dim().

```r
ncol(prDat)
```

```
## [1] 6
```

```r
length(prDat)
```

```
## [1] 6
```


Inspect the first few observations or the last few or a random sample. 

```r
head(prDat)
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_20     20      E16    wt     10.220   7.462      7.370
## Sample_21     21      E16    wt     10.020   6.890      7.177
## Sample_22     22      E16    wt      9.642   6.720      7.350
## Sample_23     23      E16    wt      9.652   6.529      7.040
## Sample_16     16      E16 NrlKO      8.583   6.470      7.494
## Sample_17     17      E16 NrlKO     10.140   7.065      7.005
```

```r
tail(prDat)
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_38     38  4_weeks    wt      9.767   6.608      7.329
## Sample_39     39  4_weeks    wt     10.200   7.003      7.320
## Sample_11     11  4_weeks NrlKO      9.677   7.204      6.981
## Sample_12     12  4_weeks NrlKO      9.129   7.165      7.350
## Sample_2       2  4_weeks NrlKO      9.744   7.107      7.075
## Sample_9       9  4_weeks NrlKO      9.822   6.558      7.043
```

```r
prDat[sample(nrow(prDat), size = 6), ]
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_25     25       P2    wt      9.038   6.170      7.449
## Sample_36     36  4_weeks    wt      9.960   7.866      6.993
## Sample_19     19      P10 NrlKO      9.771   7.081      7.586
## Sample_11     11  4_weeks NrlKO      9.677   7.204      6.981
## Sample_6       6      E16 NrlKO     10.340   7.017      6.735
## Sample_3       3       P2 NrlKO      9.414   6.166      7.200
```


What does row correspond to -- different genes or different mice?
Different mice

What are the variable names? Hint: names(), dimnames().

```r
names(prDat)
```

```
## [1] "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```

```r
dimnames(prDat)
```

```
## [[1]]
##  [1] "Sample_20" "Sample_21" "Sample_22" "Sample_23" "Sample_16"
##  [6] "Sample_17" "Sample_6"  "Sample_24" "Sample_25" "Sample_26"
## [11] "Sample_27" "Sample_14" "Sample_3"  "Sample_5"  "Sample_8" 
## [16] "Sample_28" "Sample_29" "Sample_30" "Sample_31" "Sample_1" 
## [21] "Sample_10" "Sample_4"  "Sample_7"  "Sample_32" "Sample_33"
## [26] "Sample_34" "Sample_35" "Sample_13" "Sample_15" "Sample_18"
## [31] "Sample_19" "Sample_36" "Sample_37" "Sample_38" "Sample_39"
## [36] "Sample_11" "Sample_12" "Sample_2"  "Sample_9" 
## 
## [[2]]
## [1] "sample"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```


What "flavor" is each variable, i.e. numeric, character, factor? 

```r
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```


For sample, do a sanity check that each integer between 1 and the number of rows in the dataset occurs exactly once.

```r
identical(sort(prDat$sample), seq_len(nrow(prDat)))
```

```
## [1] TRUE
```


For each factor variable, what are the levels? Hint: levels(), str().

```r
levels(prDat$devStage)
```

```
## [1] "4_weeks" "E16"     "P10"     "P2"      "P6"
```

```r
str(prDat$devStage)
```

```
##  Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
```

```r
levels(prDat$gType)
```

```
## [1] "NrlKO" "wt"
```

```r
str(prDat$gType)
```

```
##  Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
```



How many observations do we have for each level of devStage? For gType? Hint: summary(), table().

```r
summary(prDat$devStage)
```

```
## 4_weeks     E16     P10      P2      P6 
##       8       7       8       8       8
```

```r
table(prDat$gType)
```

```
## 
## NrlKO    wt 
##    19    20
```


Perform a cross-tabulation of devStage and gType. Hint: table().

```r
with(prDat, table(devStage, gType))
```

```
##          gType
## devStage  NrlKO wt
##   4_weeks     4  4
##   E16         3  4
##   P10         4  4
##   P2          4  4
##   P6          4  4
```

or

```r
table(c(prDat["devStage"], prDat["gType"]))
```

```
##          gType
## devStage  NrlKO wt
##   4_weeks     4  4
##   E16         3  4
##   P10         4  4
##   P2          4  4
##   P6          4  4
```


For each quantitative variable, what are the extremes? How about average or median?

```r
quantile(prDat$sample)
```

```
##   0%  25%  50%  75% 100% 
##  1.0 10.5 20.0 29.5 39.0
```

```r
summary(prDat$crabHammer)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    8.21    8.94    9.61    9.43    9.83   10.30
```

```r
quantile(prDat$eggBomb)
```

```
##    0%   25%   50%   75%  100% 
## 6.138 6.278 6.757 7.095 8.173
```

```r
summary(prDat$poisonFang)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    6.74    7.19    7.35    7.38    7.48    8.58
```


Create a new data.frame called weeDat only containing observations for 
which expression of poisonFang is above 7.5.

```r
weeDat <- subset(prDat, poisonFang > 7.5)
weeDat
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_24     24       P2    wt      8.869   6.587      7.508
## Sample_26     26       P2    wt      9.611   6.870      7.511
## Sample_27     27       P2    wt      8.613   6.800      7.843
## Sample_8       8       P2 NrlKO      9.116   6.264      8.016
## Sample_7       7       P6 NrlKO      8.803   6.188      7.754
## Sample_33     33      P10    wt      9.004   7.082      8.086
## Sample_34     34      P10    wt      8.519   6.757      8.584
## Sample_15     15      P10 NrlKO      9.746   7.226      7.786
## Sample_19     19      P10 NrlKO      9.771   7.081      7.586
```


For how many observations poisonFang > 7.5? 
How do they break down by genotype and developmental stage?

```r
str(weeDat)  # 9 observations
```

```
## 'data.frame':	9 obs. of  6 variables:
##  $ sample    : int  24 26 27 8 7 33 34 15 19
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 4 4 4 4 5 3 3 3 3
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 1 1 2 2 1 1
##  $ crabHammer: num  8.87 9.61 8.61 9.12 8.8 ...
##  $ eggBomb   : num  6.59 6.87 6.8 6.26 6.19 ...
##  $ poisonFang: num  7.51 7.51 7.84 8.02 7.75 ...
```

```r
summary(weeDat$devStage)
```

```
## 4_weeks     E16     P10      P2      P6 
##       0       0       4       4       1
```

```r
with(weeDat, table(devStage, gType))
```

```
##          gType
## devStage  NrlKO wt
##   4_weeks     0  0
##   E16         0  0
##   P10         2  2
##   P2          1  3
##   P6          1  0
```



Print the observations with row names "Sample_16" and "Sample_38" to screen, showing only the 3 gene expression variables.

```r
prDat[c("Sample_16", "Sample_38"), c("crabHammer", "eggBomb", "poisonFang")]
```

```
##           crabHammer eggBomb poisonFang
## Sample_16      8.583   6.470      7.494
## Sample_38      9.767   6.608      7.329
```


Which samples have expression of eggBomb less than the 0.10 quantile?

```r
prDat[prDat$eggBomb < quantile(prDat$eggBomb, 0.1), ]
```

```
##           sample devStage gType crabHammer eggBomb poisonFang
## Sample_25     25       P2    wt      9.038   6.170      7.449
## Sample_14     14       P2 NrlKO      9.572   6.138      7.250
## Sample_3       3       P2 NrlKO      9.414   6.166      7.200
## Sample_35     35      P10    wt      8.449   6.155      7.201
```

