
library(SpotClean)
library(S4Vectors)
library(SummarizedExperiment)
# Load 10x Visium data
mbrain_raw <- read10xRaw("raw_feature_bc_matrix/")
mbrain_slide_info <- read10xSlide("spatial/tissue_positions_list.csv", 
                                  "spatial/tissue_lowres_image.png", 
                                  "spatial/scalefactors_json.json")

# Visualize raw data
mbrain_obj <- createSlide(count_mat = mbrain_raw, 
                          slide_info = mbrain_slide_info)
visualizeSlide(slide_obj = mbrain_obj)

#Visualize examplel gene
g <- 'Slc17a7'
visualizeHeatmap(mbrain_obj,g)

# Decontaminate raw data using spotclean
decont_obj <- spotclean(mbrain_obj)
g <- 'Slc17a7'
visualizeHeatmap(mbrain_obj,g)
visualizeHeatmap(decont_obj,g)


dim(mbrain_obj)
dim(decont_obj)
dim(mbrain_obj[, colnames(decont_obj)])

visualizeHeatmap(mbrain_obj[, colnames(decont_obj)],g)
visualizeHeatmap(decont_obj,g)

original <- mbrain_obj[g, colnames(decont_obj)]
new <- decont_obj[g,]
x <- original@assays@data$raw
y <- new@assays@data$decont
plot(x,y, main=g, pch=".")
cor.test(as.numeric(x),as.numeric(y))

############ STdeconvolve
library(STdeconvolve)
class(decont_obj)
cd <- decont_obj@assays@data$decont
cd <- floor(cd)

counts2 <- cleanCounts(cd, 
                       min.lib.size = 10^3.75,
                       min.reads = 2000)
## feature select for genes
corpus2 <- restrictCorpus(counts2, removeAbove=1.0, removeBelow = 0.05)
## choose optimal number of cell-types
ldas2 <- fitLDA(t(as.matrix(corpus2)), Ks = c(25))
## get best model results
optLDA2 <- optimalModel(models = ldas2, opt = "25")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results2 <- getBetaTheta(optLDA2, perc.filt = 0.05, betaScale = 1000)
deconProp2 <- results2$theta
deconGexp2 <- results2$beta
## visualize deconvolved cell-type proportions
pos.info <- mbrain_slide_info$slide
pos <- pos.info[, c('imagerow', 'imagecol')]
rownames(pos) <- pos.info$barcode
## fix position
pos <- pos[, c(2,1)]
pos[,2] <- -pos[,2] 
pos[,2] <- pos[,2] - min(pos[,2])
par(mfrow=c(1,1))
plot(pos)
ppos2 <- pos[rownames(deconProp2),]
colnames(ppos2) <- c('x', 'y')
par(mfrow=c(1,1))
plot(ppos2)

dim(deconProp2)
dim(ppos2)
vizAllTopics(deconProp2, ppos2, r=2.5, lwd=0)	 

## save file
save(deconProp2, deconGexp2, ppos2, file="results2.RData")
