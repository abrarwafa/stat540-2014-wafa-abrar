library(lattice)
kDat <- read.table("../data/GSE4051_MINI.txt", header = TRUE, row.names = 1) 
str(kDat)
table(kDat$devStage)
table(kDat$gType)
with(kDat, table(devStage, gType))
xyplot(  eggBomb + poisonFang ~ crabHammer, kDat, outer = TRUE, grid = TRUE, groups = gType, auto.key = TRUE)
nDat <- with(kDat, data.frame(sample, devStage, gType, crabHammer,probeset = factor(rep(c("eggBomb", "poisonFang"), each = nrow(kDat))),geneExp = c(eggBomb, poisonFang)))
str(nDat)
xyplot(geneExp ~ crabHammer | probeset, nDat,grid = TRUE,groups = devStage, auto.key = TRUE)
stripplot(~ geneExp, nDat)
oDat <-with(kDat,data.frame(sample, devStage, gType,probeset = factor(rep(c("crabHammer", "eggBomb","poisonFang"), each = nrow(kDat))),geneExp = c(crabHammer, eggBomb, poisonFang)))
stripplot(~ geneExp, oDat)
stripplot( ~ geneExp | probeset, oDat, layout = c(nlevels(oDat$probeset), 1))
stripplot(geneExp ~ devStage | probeset, oDat, layout = c(nlevels(oDat$probeset), 1), groups = gType, auto.key = TRUE, grid = TRUE,type = 'b')
jBw <- 0.2
jn <- 400
densityplot(~ geneExp, oDat,groups = gType, auto.key = TRUE,bw = jBw, n = jn,main = paste("bw =", jBw, ", n =", jn))


prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat)
set.seed(1)
(yo <- sample(1:nrow(prDat), size = 50))
hDat <- prDat[yo, ]
prDes <- readRDS("../data/GSE4051_design.rds")
hDat <- as.matrix(t(hDat))
rownames(hDat) <- with(prDes,paste(devStage, gType, sidChar, sep="_"))
str(hDat)
heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8))
library(RColorBrewer)
display.brewer.all()

set.seed(924)
(yo <- sample(1:ncol(prDat), size = 2))
y <- prDat[[yo[1]]]
z <- prDat[[yo[2]]]
str(y)
xyplot(y ~ z, asp = 1)
smoothScatter(y ~ z, asp = 1)
install.packages("hexbin")

library(hexbin)
hexbinplot(y ~ z)

set.seed(3)
(yo <- sample(1:ncol(prDat), size = 4))
pairDat <- subset(prDat, select = yo)
pairs(pairDat, panel = function(...) smoothScatter(..., add=TRUE))


prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level=0)
set.seed(4)
(yo <- sample(1:ncol(prDat), size = 20))
hDat <- prDat[yo, ]
str(hDat)

