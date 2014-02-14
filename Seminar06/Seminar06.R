library(limma)
library(lattice)
source("../Seminar05/Seminar05.R")

prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)


# Optional take-home exercise
set.seed(400)
m <- 1000 # genes
n <- 3 #sample size
g <- 2 #gType
x <- matrix(rnorm(m * n, mean=2, sd=5), nrow = m)
colnames(x)=paste0("sample",seq_len(ncol(x)))
xDat <- data.frame(x, group = rep(paste("g", 1:g, sep=""), each=m/g))

obsVars <- apply(x, 1, var)
summary(obsVars)
densityplot(~obsVars, n=200)

#Fit a linear model
wtDes <- subset(prDes, gType == "wt")
str(wtDes)
wtDat <- subset(prDat, select = prDes$gType == "wt")
str(wtDat, max.level = 0)

wtDesMat <- model.matrix(~devStage, wtDes)
wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)
topTable(wtEbFit)
colnames(coef(wtEbFit))
dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit))), number=nrow(wtDat))
str(dsHits)
stripplot(gExp ~ devStage | gene, prepareData((rownames(dsHits)[c(3,6,9)])), group = gType, jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE, subset=prDes$gType == "wt")

#Exercise
wetFit <- lm(formula = gExp ~ devStage, data = prepareData((rownames(dsHits)[3])), subset=gType =="wt")
summary(wetFit)


# Toptabel Boss
AdjHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit))), number=nrow(wtDat), p.value=1e-05)
str(AdjHits)
AdjHits[63, c("F","adj.P.Val","devStageP6")]

P2hit <- topTable(wtEbFit, coef = "devStageP2", number=nrow(wtDat), sort.by="none")
P10hit <- topTable(wtEbFit, coef = "devStageP10", number=nrow(wtDat), sort.by="none")

xyplot(P10hit$t~P2hit$t, aspect=1, scales=list(limit=(c(range(P2hit,P10hit)))), panel=function(...){panel.smoothScatter(...); panel.abline(0,1, col=2)})

densityplot(~P2hit$adj.P.Val + P10hit$adj.P.Val, auto.key=TRUE)

addmargins(table(P2hit$adj.P.Val<1e-03,P10hit$adj.P.Val<1e-03, dnn=c("P2","P10")))

P10BY <- topTable(wtEbFit, coef = "devStageP10", number=nrow(wtDat), sort.by="none", adjust.method="BY")
Pvalues <- data.frame(raw=P10hit$P.Value, BH=P10hit$adj.P.Val, BY=P10BY$adj.P.Val)
pairs(Pvalues)

# Perform inference for some contrasts/Take home exercise
colnames(wtDesMat) 
(cont.matrix <- makeContrasts(P2VsI = devStageP2 - Intercept, P6VsP2 = devStageP6 - devStageP2, P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - devStageP10, levels = wtDesMat))
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
wtEbFitCont <- eBayes(wtFitCont)
ContHit <- topTable(wtEbFitCont)

cutoff <- 1e-04
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)
(hits <- rownames(prDat)[which(wtResCont[, "P2VsI"] != 0 & wtResCont[, "P6VsP2"] !=  0 & wtResCont[, "P10VsP6"] ==  0 & wtResCont[, "fourweeksVsP10"] ==  0 )])
stripplot(gExp ~ devStage | gene, prepareData(hits[1:5]), group = gType, jitter.data = TRUE, auto.key = TRUE, type = c('p', 'a'), grid = TRUE, subset=prDes$gType == "wt")
