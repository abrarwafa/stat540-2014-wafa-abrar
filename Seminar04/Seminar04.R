library(lattice)
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)

# Two sample tests -- one gene
set.seed(987)
(theGene <- sample(1:nrow(prDat), 1))
pDat <- data.frame(prDes, gExp = unlist(prDat[theGene, ]))
str(pDat)

aggregate(gExp ~ gType, pDat, FUN = mean)
stripplot(gType ~ gExp, pDat)
t.test(gExp ~ gType, pDat)
ttRes <- t.test(gExp ~ gType, pDat)
str(ttRes)
ttRes$statistic
ttRes$p.value

wtRes <- wilcox.test(gExp ~ gType, pDat)
wtRes$statistic
wtRes$p.value

ktRes <- ks.test(pDat$gExp[pDat$gType == "wt"],pDat$gExp[pDat$gType == "NrlKO"])
ktRes$statistic
ktRes$p.value

# Can you pull test statistics and/or p-values from the different approaches into an common object, like a readable table?
# Yes
cbind(c(ttRes$p.value,wtRes$p.value, ktRes$p.value), c(ttRes$statistic,wtRes$statistic, ktRes$statistic))

# Are you getting the same message from the various approaches?
# Yes all show there is no difference in means


# Apply() for computing on rows and columns of matrices
kDat <- readRDS("../data/GSE4051_MINI.rds")
kMat <- as.matrix(kDat[c('crabHammer', 'eggBomb', 'poisonFang')])
str(kMat)
median(kMat[ , 'eggBomb'])
apply(kMat, 2, median)
apply(kMat, 2, quantile, probs = 0.5)
apply(kMat, 2, quantile, probs = c(0.25, 0.75))
apply(kMat, 1, min)
colnames(kMat)[apply(kMat, 1, which.min)]
all.equal(rowSums(kMat), apply(kMat, 1, sum))

# Computing on groups of observations with aggregate()
aggregate(eggBomb ~ devStage, kDat, FUN = mean)
aggregate(eggBomb ~ gType * devStage, kDat, FUN = mean)
aggregate(eggBomb ~ gType * devStage, kDat, FUN = range)


keepGenes <- c("1431708_a_at", "1424336_at", "1454696_at",
               "1416119_at", "1432141_x_at", "1429226_at" )
miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),gene = factor(rep(rownames(miniDat), each = ncol(miniDat)), levels = keepGenes))
miniDat <- suppressWarnings(data.frame(prDes, miniDat))
str(miniDat)
stripplot(gType ~ gExp | gene, miniDat,scales = list(x = list(relation = "free")),group = gType, auto.key = TRUE)
someDat <- droplevels(subset(miniDat, gene == keepGenes[1]))
t.test(gExp ~ gType, someDat)


# The plyr package
library(plyr)
d_ply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x), .print = TRUE)
ttRes <- dlply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x))
names(ttRes)
ttRes[["1454696_at"]]

ttRes <- ddply(miniDat, ~ gene, function(z) {
  zz <- t.test(gExp ~ gType, z)
  round(c(tStat = zz$statistic, pVal = zz$p.value), 4)
})
ttRes
