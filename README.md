# Reprogramming of epidermal keratinocytes by PITX1 transforms the cutaneous cellular landscape and promotes wound healing 

### Overmiller et al. 2024, JCI Insight

# Introduction

This document includes the code used to analyze datasets and
generate figures related to single-cell RNA sequencing and Xenium _in situ_ analysis from 
Overmiller et al., 2024 (doi: XXX).

## Data

Single cell RNA-seq and Xenium _in situ_ data used in this R code available here -- 
GSE280088 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280088)

## Single cell RNA-seq

### Combined Healthy Skin and Buccal Mucosa

	library(dplyr)
	library(ggplot2)
	library(Matrix)
	library(patchwork)
	library(gdata)
	library(Seurat)
	library(SeuratObject)
	library(SeuratWrappers)
	library(BiocManager)
	library(metap)
	library(cowplot)
	library(sctransform)
	library(xlsx)
	library(glmGamPoi)
	library(clustree)
	library(biomaRt)
	library(monocle3)
	library(magrittr)
	library(slingshot)
	library(paletteer)
	library(nichenetr)
	library(multinichenetr)
	library(tidyverse)
	library(circlize)
	library(scales)
	library(kableExtra)
	library(knitr)
	library(SoupX)
	library(scDblFinder)
	library(viridis)
	library(CellChat)
	library(uwot)
	library(ComplexHeatmap)
	library(enrichR)
	
	#Set current working directory
	setwd("/data/overmilleram/scRNAseq/Skin & Oral/")
	################
	
	#Preprocess data
	
	#Use SoupX to remove ambient RNA
	toc.fcs1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_1/filtered_feature_bc_matrix/")
	tod.fcs1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_1/raw_feature_bc_matrix/")
	toc.fcs2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_2/filtered_feature_bc_matrix/")
	tod.fcs2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_2/raw_feature_bc_matrix/")
	toc.fcs3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_3/filtered_feature_bc_matrix/")
	tod.fcs3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_3/raw_feature_bc_matrix/")
	toc.fcs4 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_4/filtered_feature_bc_matrix/")
	tod.fcs4 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_4/raw_feature_bc_matrix/")
	toc.fcs5 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_5/filtered_feature_bc_matrix/")
	tod.fcs5 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Control_5/raw_feature_bc_matrix/")
	toc.mcs1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Control_1/filtered_feature_bc_matrix/")
	tod.mcs1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Control_1/raw_feature_bc_matrix/")
	toc.mcs2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Control_2/filtered_feature_bc_matrix/")
	tod.mcs2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Control_2/raw_feature_bc_matrix/")
	toc.mcs3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Control_3/filtered_feature_bc_matrix/")
	tod.mcs3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Control_3/raw_feature_bc_matrix/")
	toc.fds1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_1/filtered_feature_bc_matrix/")
	tod.fds1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_1/raw_feature_bc_matrix/")
	toc.fds2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_2/filtered_feature_bc_matrix/")
	tod.fds2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_2/raw_feature_bc_matrix/")
	toc.fds3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_3/filtered_feature_bc_matrix/")
	tod.fds3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_3/raw_feature_bc_matrix/")
	toc.fds4 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_4/filtered_feature_bc_matrix/")
	tod.fds4 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Female_Dox_4/raw_feature_bc_matrix/")
	toc.mds1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_1/filtered_feature_bc_matrix/")
	tod.mds1 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_1/raw_feature_bc_matrix/")
	toc.mds2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_2/filtered_feature_bc_matrix/")
	tod.mds2 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_2/raw_feature_bc_matrix/")
	toc.mds3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_3/filtered_feature_bc_matrix/")
	tod.mds3 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_3/raw_feature_bc_matrix/")
	toc.mds4 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_4/filtered_feature_bc_matrix/")
	tod.mds4 <- Seurat::Read10X("../../Sample Matrices/D0 Skin/8w_Male_Dox_4/raw_feature_bc_matrix/")
	toc.fbm1 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_1/filtered_feature_bc_matrix/")
	tod.fbm1 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_1/raw_feature_bc_matrix/")
	toc.fbm2 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_2/filtered_feature_bc_matrix/")
	tod.fbm2 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_2/raw_feature_bc_matrix/")
	toc.fbm3 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_3/filtered_feature_bc_matrix/")
	tod.fbm3 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_3/raw_feature_bc_matrix/")
	toc.fbm4 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_4/filtered_feature_bc_matrix/")
	tod.fbm4 <- Seurat::Read10X("../../Sample Matrices/Oral/Female_BM_Control_4/raw_feature_bc_matrix/")
	toc.mbm1 <- Seurat::Read10X("../../Sample Matrices/Oral/Male_BM_Control_1/filtered_feature_bc_matrix/")
	tod.mbm1 <- Seurat::Read10X("../../Sample Matrices/Oral/Male_BM_Control_1/raw_feature_bc_matrix/")
	toc.mbm2 <- Seurat::Read10X("../../Sample Matrices/Oral/Male_BM_Control_2/filtered_feature_bc_matrix/")
	tod.mbm2 <- Seurat::Read10X("../../Sample Matrices/Oral/Male_BM_Control_2/raw_feature_bc_matrix/")
	toc.mbm3 <- Seurat::Read10X("../../Sample Matrices/Oral/Male_BM_Control_3/filtered_feature_bc_matrix/")
	tod.mbm3 <- Seurat::Read10X("../../Sample Matrices/Oral/Male_BM_Control_3/raw_feature_bc_matrix/")
	
	sc.fcs1 <- SoupChannel(tod.fcs1,toc.fcs1, calcSoupProfile = T)
	sc.fcs2 <- SoupChannel(tod.fcs2,toc.fcs2, calcSoupProfile = T)
	sc.fcs3 <- SoupChannel(tod.fcs3,toc.fcs3, calcSoupProfile = T)
	sc.fcs4 <- SoupChannel(tod.fcs4,toc.fcs4, calcSoupProfile = T)
	sc.fcs5 <- SoupChannel(tod.fcs5,toc.fcs5, calcSoupProfile = T)
	sc.mcs1 <- SoupChannel(tod.mcs1,toc.mcs1, calcSoupProfile = T)
	sc.mcs2 <- SoupChannel(tod.mcs2,toc.mcs2, calcSoupProfile = T)
	sc.mcs3 <- SoupChannel(tod.mcs3,toc.mcs3, calcSoupProfile = T)
	sc.fds1 <- SoupChannel(tod.fds1,toc.fds1, calcSoupProfile = T)
	sc.fds2 <- SoupChannel(tod.fds2,toc.fds2, calcSoupProfile = T)
	sc.fds3 <- SoupChannel(tod.fds3,toc.fds3, calcSoupProfile = T)
	sc.fds4 <- SoupChannel(tod.fds4,toc.fds4, calcSoupProfile = T)
	sc.mds1 <- SoupChannel(tod.mds1,toc.mds1, calcSoupProfile = T)
	sc.mds2 <- SoupChannel(tod.mds2,toc.mds2, calcSoupProfile = T)
	sc.mds3 <- SoupChannel(tod.mds3,toc.mds3, calcSoupProfile = T)
	sc.mds4 <- SoupChannel(tod.mds4,toc.mds4, calcSoupProfile = T)
	sc.fbm1 <- SoupChannel(tod.fbm1,toc.fbm1, calcSoupProfile = T)
	sc.fbm2 <- SoupChannel(tod.fbm2,toc.fbm2, calcSoupProfile = T)
	sc.fbm3 <- SoupChannel(tod.fbm3,toc.fbm3, calcSoupProfile = T)
	sc.fbm4 <- SoupChannel(tod.fbm4,toc.fbm4, calcSoupProfile = T)
	sc.mbm1 <- SoupChannel(tod.mbm1,toc.mbm1, calcSoupProfile = T)
	sc.mbm2 <- SoupChannel(tod.mbm2,toc.mbm2, calcSoupProfile = T)
	sc.mbm3 <- SoupChannel(tod.mbm3,toc.mbm3, calcSoupProfile = T)
	
	seur.fcs1 <- CreateSeuratObject(sc.fcs1$toc)
	seur.fcs2 <- CreateSeuratObject(sc.fcs2$toc)
	seur.fcs3 <- CreateSeuratObject(sc.fcs3$toc)
	seur.fcs4 <- CreateSeuratObject(sc.fcs4$toc)
	seur.fcs5 <- CreateSeuratObject(sc.fcs5$toc)
	seur.mcs1 <- CreateSeuratObject(sc.mcs1$toc)
	seur.mcs2 <- CreateSeuratObject(sc.mcs2$toc)
	seur.mcs3 <- CreateSeuratObject(sc.mcs3$toc)
	seur.fds1 <- CreateSeuratObject(sc.fds1$toc)
	seur.fds2 <- CreateSeuratObject(sc.fds2$toc)
	seur.fds3 <- CreateSeuratObject(sc.fds3$toc)
	seur.fds4 <- CreateSeuratObject(sc.fds4$toc)
	seur.mds1 <- CreateSeuratObject(sc.mds1$toc)
	seur.mds2 <- CreateSeuratObject(sc.mds2$toc)
	seur.mds3 <- CreateSeuratObject(sc.mds3$toc)
	seur.mds4 <- CreateSeuratObject(sc.mds4$toc)
	seur.fbm1 <- CreateSeuratObject(sc.fbm1$toc)
	seur.fbm2 <- CreateSeuratObject(sc.fbm2$toc)
	seur.fbm3 <- CreateSeuratObject(sc.fbm3$toc)
	seur.fbm4 <- CreateSeuratObject(sc.fbm4$toc)
	seur.mbm1 <- CreateSeuratObject(sc.mbm1$toc)
	seur.mbm2 <- CreateSeuratObject(sc.mbm2$toc)
	seur.mbm3 <- CreateSeuratObject(sc.mbm3$toc)
	
	seur.list <- list(seur.fcs1=seur.fcs1, seur.fcs2=seur.fcs2, seur.fcs3=seur.fcs3, seur.fcs4=seur.fcs4, seur.fcs5=seur.fcs5,
	                  seur.mcs1=seur.mcs1, seur.mcs2=seur.mcs2, seur.mcs3=seur.mcs3,
	                  seur.fds1=seur.fds1, seur.fds2=seur.fds2, seur.fds3=seur.fds3, seur.fds4=seur.fds4, 
	                  seur.mds1=seur.mds1, seur.mds2=seur.mds2, seur.mds3=seur.mds3, seur.mds4=seur.mds4,
	                  seur.fbm1=seur.fbm1, seur.fbm2=seur.fbm2, seur.fbm3=seur.fbm3, seur.fbm4=seur.fbm4,
	                  seur.mbm1=seur.mbm1, seur.mbm2=seur.mbm2, seur.mbm3=seur.mbm3)
	
	seur.list <- lapply(X = seur.list, 
	                    FUN = function(x) {
	                      x = NormalizeData(x)
	                      x = FindVariableFeatures(x, 
	                                               selection.method = "vst",
	                                               nfeatures = 5000)
	                      x = ScaleData(x, verbose = TRUE)
	                      x = RunPCA(x, 
	                                 npcs = 50, verbose = TRUE)
	                      x = FindNeighbors(x, 
	                                        reduction = "pca", 
	                                        dims = 1:50)
	                      x = FindClusters(x,
	                                       algorithm = 3,
	                                       resolution = 1)
	                    })
	
	list2env(seur.list, 
	         .GlobalEnv)
	
	sc.fcs1 <- setClusters(sc.fcs1, seur.fcs1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fcs2 <- setClusters(sc.fcs2, seur.fcs2$seurat_clusters) %>% autoEstCont() %>% adjustCounts() 
	sc.fcs3 <- setClusters(sc.fcs3, seur.fcs3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fcs4 <- setClusters(sc.fcs4, seur.fcs4$seurat_clusters) %>% autoEstCont() %>% adjustCounts() 
	sc.fcs5 <- setClusters(sc.fcs5, seur.fcs5$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mcs1 <- setClusters(sc.mcs1, seur.mcs1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mcs2 <- setClusters(sc.mcs2, seur.mcs2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mcs3 <- setClusters(sc.mcs3, seur.mcs3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds1 <- setClusters(sc.fds1, seur.fds1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds2 <- setClusters(sc.fds2, seur.fds2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds3 <- setClusters(sc.fds3, seur.fds3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds4 <- setClusters(sc.fds4, seur.fds4$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mds1 <- setClusters(sc.mds1, seur.mds1$seurat_clusters) %>% autoEstCont(forceAccept = T) %>% adjustCounts() # high amounts of ambient RNA
	sc.mds2 <- setClusters(sc.mds2, seur.mds2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mds3 <- setClusters(sc.mds3, seur.mds3$seurat_clusters) %>% autoEstCont(forceAccept = T) %>% adjustCounts() # high amounts of ambient RNA
	sc.mds4 <- setClusters(sc.mds4, seur.mds4$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fbm1 <- setClusters(sc.fbm1, seur.fbm1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fbm2 <- setClusters(sc.fbm2, seur.fbm2$seurat_clusters) %>% autoEstCont() %>% adjustCounts() 
	sc.fbm3 <- setClusters(sc.fbm3, seur.fbm3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fbm4 <- setClusters(sc.fbm4, seur.fbm4$seurat_clusters) %>% autoEstCont() %>% adjustCounts() 
	sc.mbm1 <- setClusters(sc.mbm1, seur.mbm1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mbm2 <- setClusters(sc.mbm2, seur.mbm2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mbm3 <- setClusters(sc.mbm3, seur.mbm3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	
	fcs1 <- CreateSeuratObject(sc.fcs1, min.cells = 3, min.features = 200)
	fcs2 <- CreateSeuratObject(sc.fcs2, min.cells = 3, min.features = 200)
	fcs3 <- CreateSeuratObject(sc.fcs3, min.cells = 3, min.features = 200)
	fcs4 <- CreateSeuratObject(sc.fcs4, min.cells = 3, min.features = 200)
	fcs5 <- CreateSeuratObject(sc.fcs5, min.cells = 3, min.features = 200)
	mcs1 <- CreateSeuratObject(sc.mcs1, min.cells = 3, min.features = 200)
	mcs2 <- CreateSeuratObject(sc.mcs2, min.cells = 3, min.features = 200)
	mcs3 <- CreateSeuratObject(sc.mcs3, min.cells = 3, min.features = 200)
	fds1 <- CreateSeuratObject(sc.fds1, min.cells = 3, min.features = 200)
	fds2 <- CreateSeuratObject(sc.fds2, min.cells = 3, min.features = 200)
	fds3 <- CreateSeuratObject(sc.fds3, min.cells = 3, min.features = 200)
	fds4 <- CreateSeuratObject(sc.fds4, min.cells = 3, min.features = 200)
	mds1 <- CreateSeuratObject(sc.mds1, min.cells = 3, min.features = 200)
	mds2 <- CreateSeuratObject(sc.mds2, min.cells = 3, min.features = 200)
	mds3 <- CreateSeuratObject(sc.mds3, min.cells = 3, min.features = 200)
	mds4 <- CreateSeuratObject(sc.mds4, min.cells = 3, min.features = 200)
	fbm1 <- CreateSeuratObject(sc.fbm1, min.cells = 3, min.features = 200)
	fbm2 <- CreateSeuratObject(sc.fbm2, min.cells = 3, min.features = 200)
	fbm3 <- CreateSeuratObject(sc.fbm3, min.cells = 3, min.features = 200)
	fbm4 <- CreateSeuratObject(sc.fbm4, min.cells = 3, min.features = 200)
	mbm1 <- CreateSeuratObject(sc.mbm1, min.cells = 3, min.features = 200)
	mbm2 <- CreateSeuratObject(sc.mbm2, min.cells = 3, min.features = 200)
	mbm3 <- CreateSeuratObject(sc.mbm3, min.cells = 3, min.features = 200)
	
	#Use scDblFinder to determine # of doublets
	
	set.seed(42)
	
	fcs1sce <- scDblFinder(GetAssayData(fcs1, slot = "counts"))
	fcs2sce <- scDblFinder(GetAssayData(fcs2, slot = "counts"))
	fcs3sce <- scDblFinder(GetAssayData(fcs3, slot = "counts"))
	fcs4sce <- scDblFinder(GetAssayData(fcs4, slot = "counts"))
	fcs5sce <- scDblFinder(GetAssayData(fcs5, slot = "counts"))
	mcs1sce <- scDblFinder(GetAssayData(mcs1, slot = "counts"))
	mcs2sce <- scDblFinder(GetAssayData(mcs2, slot = "counts"))
	mcs3sce <- scDblFinder(GetAssayData(mcs3, slot = "counts"))
	fds1sce <- scDblFinder(GetAssayData(fds1, slot = "counts"))
	fds2sce <- scDblFinder(GetAssayData(fds2, slot = "counts"))
	fds3sce <- scDblFinder(GetAssayData(fds3, slot = "counts"))
	fds4sce <- scDblFinder(GetAssayData(fds4, slot = "counts"))
	mds1sce <- scDblFinder(GetAssayData(mds1, slot = "counts"))
	mds2sce <- scDblFinder(GetAssayData(mds2, slot = "counts"))
	mds3sce <- scDblFinder(GetAssayData(mds3, slot = "counts"))
	mds4sce <- scDblFinder(GetAssayData(mds4, slot = "counts"))
	fbm1sce <- scDblFinder(GetAssayData(fbm1, slot = "counts"))
	fbm2sce <- scDblFinder(GetAssayData(fbm2, slot = "counts"))
	fbm3sce <- scDblFinder(GetAssayData(fbm3, slot = "counts"))
	fbm4sce <- scDblFinder(GetAssayData(fbm4, slot = "counts"))
	mbm1sce <- scDblFinder(GetAssayData(mbm1, slot = "counts"))
	mbm2sce <- scDblFinder(GetAssayData(mbm2, slot = "counts"))
	mbm3sce <- scDblFinder(GetAssayData(mbm3, slot = "counts"))
	
	fcs1$scDblFinder.class <- fcs1sce$scDblFinder.class
	fcs2$scDblFinder.class <- fcs2sce$scDblFinder.class
	fcs3$scDblFinder.class <- fcs3sce$scDblFinder.class
	fcs4$scDblFinder.class <- fcs4sce$scDblFinder.class
	fcs5$scDblFinder.class <- fcs5sce$scDblFinder.class
	mcs1$scDblFinder.class <- mcs1sce$scDblFinder.class
	mcs2$scDblFinder.class <- mcs2sce$scDblFinder.class
	mcs3$scDblFinder.class <- mcs3sce$scDblFinder.class
	fds1$scDblFinder.class <- fds1sce$scDblFinder.class
	fds2$scDblFinder.class <- fds2sce$scDblFinder.class
	fds3$scDblFinder.class <- fds3sce$scDblFinder.class
	fds4$scDblFinder.class <- fds4sce$scDblFinder.class
	mds1$scDblFinder.class <- mds1sce$scDblFinder.class
	mds2$scDblFinder.class <- mds2sce$scDblFinder.class
	mds3$scDblFinder.class <- mds3sce$scDblFinder.class
	mds4$scDblFinder.class <- mds4sce$scDblFinder.class
	fbm1$scDblFinder.class <- fbm1sce$scDblFinder.class
	fbm2$scDblFinder.class <- fbm2sce$scDblFinder.class
	fbm3$scDblFinder.class <- fbm3sce$scDblFinder.class
	fbm4$scDblFinder.class <- fbm4sce$scDblFinder.class
	mbm1$scDblFinder.class <- mbm1sce$scDblFinder.class
	mbm2$scDblFinder.class <- mbm2sce$scDblFinder.class
	mbm3$scDblFinder.class <- mbm3sce$scDblFinder.class
	
	Idents(fcs1) <- "scDblFinder.class"
	Idents(fcs2) <- "scDblFinder.class"
	Idents(fcs3) <- "scDblFinder.class"
	Idents(fcs4) <- "scDblFinder.class"
	Idents(fcs5) <- "scDblFinder.class"
	Idents(mcs1) <- "scDblFinder.class"
	Idents(mcs2) <- "scDblFinder.class"
	Idents(mcs3) <- "scDblFinder.class"
	Idents(fds1) <- "scDblFinder.class"
	Idents(fds2) <- "scDblFinder.class"
	Idents(fds3) <- "scDblFinder.class"
	Idents(fds4) <- "scDblFinder.class"
	Idents(mds1) <- "scDblFinder.class"
	Idents(mds2) <- "scDblFinder.class"
	Idents(mds3) <- "scDblFinder.class"
	Idents(mds4) <- "scDblFinder.class"
	Idents(fbm1) <- "scDblFinder.class"
	Idents(fbm2) <- "scDblFinder.class"
	Idents(fbm3) <- "scDblFinder.class"
	Idents(fbm4) <- "scDblFinder.class"
	Idents(mbm1) <- "scDblFinder.class"
	Idents(mbm2) <- "scDblFinder.class"
	Idents(mbm3) <- "scDblFinder.class"
	
	fcs1 <- subset(fcs1, idents = "singlet")
	fcs2 <- subset(fcs2, idents = "singlet")
	fcs3 <- subset(fcs3, idents = "singlet")
	fcs4 <- subset(fcs4, idents = "singlet")
	fcs5 <- subset(fcs5, idents = "singlet")
	mcs1 <- subset(mcs1, idents = "singlet")
	mcs2 <- subset(mcs2, idents = "singlet")
	mcs3 <- subset(mcs3, idents = "singlet")
	fds1 <- subset(fds1, idents = "singlet")
	fds2 <- subset(fds2, idents = "singlet")
	fds3 <- subset(fds3, idents = "singlet")
	fds4 <- subset(fds4, idents = "singlet")
	mds1 <- subset(mds1, idents = "singlet")
	mds2 <- subset(mds2, idents = "singlet")
	mds3 <- subset(mds3, idents = "singlet")
	mds4 <- subset(mds4, idents = "singlet")
	fbm1 <- subset(fbm1, idents = "singlet")
	fbm2 <- subset(fbm2, idents = "singlet")
	fbm3 <- subset(fbm3, idents = "singlet")
	fbm4 <- subset(fbm4, idents = "singlet")
	mbm1 <- subset(mbm1, idents = "singlet")
	mbm2 <- subset(mbm2, idents = "singlet")
	mbm3 <- subset(mbm3, idents = "singlet")
	
	#File cleanup
	rm(fcs1sce, fcs2sce, fcs3sce, fcs4sce, fcs5sce, fcs6sce,
	   mcs1sce, mcs2sce, mcs3sce, mcs4sce, mcs5sce,
	   fds1sce, fds2sce, fds3sce, fds4sce, fds5sce, fds6sce,
	   mds1sce, mds2sce, mds3sce, mds4sce, mds5sce, mds6sce,
	   sc.fcs1, sc.fcs2, sc.fcs3, sc.fcs4, sc.fcs5, sc.fcs6,
	   sc.mcs1, sc.mcs2, sc.mcs3, sc.mcs4, sc.mcs5,
	   sc.fds1, sc.fds2, sc.fds3, sc.fds4, sc.fds5, sc.fds6,
	   sc.mds1, sc.mds2, sc.mds3, sc.mds4, sc.mds5, sc.mds6,
	   seur.fcs1, seur.fcs2, seur.fcs3, seur.fcs4, seur.fcs5, seur.fcs6,
	   seur.mcs1, seur.mcs2, seur.mcs3, seur.mcs4, seur.mcs5,
	   seur.fds1, seur.fds2, seur.fds3, seur.fds4, seur.fds5, seur.fds6,
	   seur.mds1, seur.mds2, seur.mds3, seur.mds4, seur.mds5, seur.mds6,
	   tod.fcs1, tod.fcs2, tod.fcs3, tod.fcs4, tod.fcs5, tod.fcs6,
	   tod.mcs1, tod.mcs2, tod.mcs3, tod.mcs4, tod.mcs5,
	   tod.fds1, tod.fds2, tod.fds3, tod.fds4, tod.fds5, tod.fds6,
	   tod.mds1, tod.mds2, tod.mds3, tod.mds4, tod.mds5, tod.mds6,
	   toc.fcs1, toc.fcs2, toc.fcs3, toc.fcs4, toc.fcs5, toc.fcs6,
	   toc.mcs1, toc.mcs2, toc.mcs3, toc.mcs4, toc.mcs5,
	   toc.fds1, toc.fds2, toc.fds3, toc.fds4, toc.fds5, toc.fds6,
	   toc.mds1, toc.mds2, toc.mds3, toc.mds4, toc.mds5, toc.mds6,
	   fcs1sce, fcs2sce, fcs3sce, fcs4sce, fcs5sce, fcs6sce,
	   mcs1sce, mcs2sce, mcs3sce, mcs4sce, mcs5sce,
	   fds1sce, fds2sce, fds3sce, fds4sce, fds5sce, fds6sce,
	   mds1sce, mds2sce, mds3sce, mds4sce, mds5sce, mds6sce,
	   fcbm1sce, fcbm2sce, fchp1sce, fchp2sce, mcbm1sce, mchp1sce,
	   sc.fcs1, sc.fcs2, sc.fcs3, sc.fcs4, sc.fcs5, sc.fcs6,
	   sc.mcs1, sc.mcs2, sc.mcs3, sc.mcs4, sc.mcs5,
	   sc.fds1, sc.fds2, sc.fds3, sc.fds4, sc.fds5, sc.fds6,
	   sc.mds1, sc.mds2, sc.mds3, sc.mds4, sc.mds5, sc.mds6,
	   sc.fcbm1, sc.fcbm2, sc.fchp1, sc.fchp2, sc.mcbm1, sc.mchp1,
	   seur.fcs1, seur.fcs2, seur.fcs3, seur.fcs4, seur.fcs5, seur.fcs6,
	   seur.mcs1, seur.mcs2, seur.mcs3, seur.mcs4, seur.mcs5,
	   seur.fds1, seur.fds2, seur.fds3, seur.fds4, seur.fds5, seur.fds6,
	   seur.mds1, seur.mds2, seur.mds3, seur.mds4, seur.mds5, seur.mds6,
	   seur.fcbm1, seur.fcbm2, seur.fchp1, seur.fchp2, seur.mcbm1, seur.mchp1,
	   tod.fcs1, tod.fcs2, tod.fcs3, tod.fcs4, tod.fcs5, tod.fcs6,
	   tod.mcs1, tod.mcs2, tod.mcs3, tod.mcs4, tod.mcs5,
	   tod.fds1, tod.fds2, tod.fds3, tod.fds4, tod.fds5, tod.fds6,
	   tod.mds1, tod.mds2, tod.mds3, tod.mds4, tod.mds5, tod.mds6,
	   tod.fcbm1, tod.fcbm2, tod.fchp1, tod.fchp2, tod.mcbm1, tod.mchp1,
	   toc.fcs1, toc.fcs2, toc.fcs3, toc.fcs4, toc.fcs5, toc.fcs6,
	   toc.mcs1, toc.mcs2, toc.mcs3, toc.mcs4, toc.mcs5,
	   toc.fds1, toc.fds2, toc.fds3, toc.fds4, toc.fds5, toc.fds6,
	   toc.mds1, toc.mds2, toc.mds3, toc.mds4, toc.mds5, toc.mds6, 
	   toc.fcbm1, toc.fcbm2, toc.fchp1, toc.fchp2, toc.mcbm1, toc.mchp1, 
	   fcbm1sce, fcbm2sce, fcbm3sce, fchp1sce, fchp2sce, fchp3sce, mcbm1sce, mcbm2sce, mchp1sce, mchp2sce,
	   fcbm1.desouped, fcbm2.desouped, fcbm3.desouped, fchp1.desouped, fchp2.desouped, fchp3.desouped, mcbm1.desouped, mcbm2.desouped, mchp1.desouped, mchp2.desouped,
	   sc.fcbm1, sc.fcbm2, sc.fcbm3, sc.fchp1, sc.fchp2, sc.fchp3, sc.mcbm1, sc.mcbm2, sc.mchp1, sc.mchp2,
	   seur.fcbm1, seur.fcbm2, seur.fcbm3, seur.fchp1, seur.fchp2, seur.fchp3, seur.mcbm1, seur.mcbm2, seur.mchp1, seur.mchp2,
	   tod.fcbm1, tod.fcbm2, tod.fcbm3, tod.fchp1, tod.fchp2, tod.fchp3, tod.mcbm1, tod.mcbm2, tod.mchp1, tod.mchp2,
	   toc.fcbm1, toc.fcbm2, toc.fcbm3, toc.fchp1, toc.fchp2, toc.fchp3, toc.mcbm1, toc.mcbm2, toc.mchp1, toc.mchp2,
	   fbm1sce, fbm2sce, fbm3sce, fbm4sce, fbm5sce, fbm6sce,
	   mbm1sce, mbm2sce, mbm3sce, mbm4sce, mbm5sce,
	   fds1sce, fds2sce, fds3sce, fds4sce, fds5sce, fds6sce,
	   mds1sce, mds2sce, mds3sce, mds4sce, mds5sce, mds6sce,
	   sc.fbm1, sc.fbm2, sc.fbm3, sc.fbm4, sc.fbm5, sc.fbm6,
	   sc.mbm1, sc.mbm2, sc.mbm3, sc.mbm4, sc.mbm5,
	   sc.fds1, sc.fds2, sc.fds3, sc.fds4, sc.fds5, sc.fds6,
	   sc.mds1, sc.mds2, sc.mds3, sc.mds4, sc.mds5, sc.mds6,
	   seur.fbm1, seur.fbm2, seur.fbm3, seur.fbm4, seur.fbm5, seur.fbm6,
	   seur.mbm1, seur.mbm2, seur.mbm3, seur.mbm4, seur.mbm5,
	   seur.fds1, seur.fds2, seur.fds3, seur.fds4, seur.fds5, seur.fds6,
	   seur.mds1, seur.mds2, seur.mds3, seur.mds4, seur.mds5, seur.mds6,
	   tod.fbm1, tod.fbm2, tod.fbm3, tod.fbm4, tod.fbm5, tod.fbm6,
	   tod.mbm1, tod.mbm2, tod.mbm3, tod.mbm4, tod.mbm5,
	   tod.fds1, tod.fds2, tod.fds3, tod.fds4, tod.fds5, tod.fds6,
	   tod.mds1, tod.mds2, tod.mds3, tod.mds4, tod.mds5, tod.mds6,
	   toc.fbm1, toc.fbm2, toc.fbm3, toc.fbm4, toc.fbm5, toc.fbm6,
	   toc.mbm1, toc.mbm2, toc.mbm3, toc.mbm4, toc.mbm5,
	   toc.fds1, toc.fds2, toc.fds3, toc.fds4, toc.fds5, toc.fds6,
	   toc.mds1, toc.mds2, toc.mds3, toc.mds4, toc.mds5, toc.mds6,
	   fbm1sce, fbm2sce, fbm3sce, fbm4sce, fbm5sce, fbm6sce,
	   mbm1sce, mbm2sce, mbm3sce, mbm4sce, mbm5sce,
	   fds1sce, fds2sce, fds3sce, fds4sce, fds5sce, fds6sce,
	   mds1sce, mds2sce, mds3sce, mds4sce, mds5sce, mds6sce,
	   fcbm1sce, fcbm2sce, fchp1sce, fchp2sce, mcbm1sce, mchp1sce,
	   sc.fbm1, sc.fbm2, sc.fbm3, sc.fbm4, sc.fbm5, sc.fbm6,
	   sc.mbm1, sc.mbm2, sc.mbm3, sc.mbm4, sc.mbm5,
	   sc.fds1, sc.fds2, sc.fds3, sc.fds4, sc.fds5, sc.fds6,
	   sc.mds1, sc.mds2, sc.mds3, sc.mds4, sc.mds5, sc.mds6,
	   sc.fcbm1, sc.fcbm2, sc.fchp1, sc.fchp2, sc.mcbm1, sc.mchp1,
	   seur.fbm1, seur.fbm2, seur.fbm3, seur.fbm4, seur.fbm5, seur.fbm6,
	   seur.mbm1, seur.mbm2, seur.mbm3, seur.mbm4, seur.mbm5,
	   seur.fds1, seur.fds2, seur.fds3, seur.fds4, seur.fds5, seur.fds6,
	   seur.mds1, seur.mds2, seur.mds3, seur.mds4, seur.mds5, seur.mds6,
	   seur.fcbm1, seur.fcbm2, seur.fchp1, seur.fchp2, seur.mcbm1, seur.mchp1,
	   tod.fbm1, tod.fbm2, tod.fbm3, tod.fbm4, tod.fbm5, tod.fbm6,
	   tod.mbm1, tod.mbm2, tod.mbm3, tod.mbm4, tod.mbm5,
	   tod.fds1, tod.fds2, tod.fds3, tod.fds4, tod.fds5, tod.fds6,
	   tod.mds1, tod.mds2, tod.mds3, tod.mds4, tod.mds5, tod.mds6,
	   tod.fcbm1, tod.fcbm2, tod.fchp1, tod.fchp2, tod.mcbm1, tod.mchp1,
	   toc.fbm1, toc.fbm2, toc.fbm3, toc.fbm4, toc.fbm5, toc.fbm6,
	   toc.mbm1, toc.mbm2, toc.mbm3, toc.mbm4, toc.mbm5,
	   toc.fds1, toc.fds2, toc.fds3, toc.fds4, toc.fds5, toc.fds6,
	   toc.mds1, toc.mds2, toc.mds3, toc.mds4, toc.mds5, toc.mds6, 
	   toc.fcbm1, toc.fcbm2, toc.fchp1, toc.fchp2, toc.mcbm1, toc.mchp1, 
	   fcbm1sce, fcbm2sce, fcbm3sce, fchp1sce, fchp2sce, fchp3sce, mcbm1sce, mcbm2sce, mchp1sce, mchp2sce,
	   fcbm1.desouped, fcbm2.desouped, fcbm3.desouped, fchp1.desouped, fchp2.desouped, fchp3.desouped, mcbm1.desouped, mcbm2.desouped, mchp1.desouped, mchp2.desouped,
	   sc.fcbm1, sc.fcbm2, sc.fcbm3, sc.fchp1, sc.fchp2, sc.fchp3, sc.mcbm1, sc.mcbm2, sc.mchp1, sc.mchp2,
	   seur.fcbm1, seur.fcbm2, seur.fcbm3, seur.fchp1, seur.fchp2, seur.fchp3, seur.mcbm1, seur.mcbm2, seur.mchp1, seur.mchp2,
	   tod.fcbm1, tod.fcbm2, tod.fcbm3, tod.fchp1, tod.fchp2, tod.fchp3, tod.mcbm1, tod.mcbm2, tod.mchp1, tod.mchp2,
	   toc.fcbm1, toc.fcbm2, toc.fcbm3, toc.fchp1, toc.fchp2, toc.fchp3, toc.mcbm1, toc.mcbm2, toc.mchp1, toc.mchp2)
	
	#Add condition identity
	fcs1$condition <- "control"
	fcs2$condition <- "control"
	fcs3$condition <- "control"
	fcs4$condition <- "control"
	fcs5$condition <- "control"
	mcs1$condition <- "control"
	mcs2$condition <- "control"
	mcs3$condition <- "control"
	fds1$condition <- "dox"
	fds2$condition <- "dox"
	fds3$condition <- "dox"
	fds4$condition <- "dox"
	mds1$condition <- "dox"
	mds2$condition <- "dox"
	mds3$condition <- "dox"
	mds4$condition <- "dox"
	fbm1$condition <- "control"
	fbm2$condition <- "control"
	fbm3$condition <- "control"
	fbm4$condition <- "control"
	mbm1$condition <- "control"
	mbm2$condition <- "control"
	mbm3$condition <- "control"
	
	#Add tissue ID
	fcs1$tissue <- "skin"
	fcs2$tissue <- "skin"
	fcs3$tissue <- "skin"
	fcs4$tissue <- "skin"
	fcs5$tissue <- "skin"
	mcs1$tissue <- "skin"
	mcs2$tissue <- "skin"
	mcs3$tissue <- "skin"
	fds1$tissue <- "skin"
	fds2$tissue <- "skin"
	fds3$tissue <- "skin"
	fds4$tissue <- "skin"
	mds1$tissue <- "skin"
	mds2$tissue <- "skin"
	mds3$tissue <- "skin"
	mds4$tissue <- "skin"
	fbm1$tissue <- "buccal_mucosa"
	fbm2$tissue <- "buccal_mucosa"
	fbm3$tissue <- "buccal_mucosa"
	fbm4$tissue <- "buccal_mucosa"
	mbm1$tissue <- "buccal_mucosa"
	mbm2$tissue <- "buccal_mucosa"
	mbm3$tissue <- "buccal_mucosa"
	
	#Add Pitx1 induction identity
	fcs1$induction <- "8_week_skin"
	fcs2$induction <- "8_week_skin"
	fcs3$induction <- "8_week_skin"
	fcs4$induction <- "8_week_skin"
	fcs5$induction <- "8_week_skin"
	mcs1$induction <- "8_week_skin"
	mcs2$induction <- "8_week_skin"
	mcs3$induction <- "8_week_skin"
	fds1$induction <- "8_week_skin"
	fds2$induction <- "8_week_skin"
	fds3$induction <- "8_week_skin"
	fds4$induction <- "8_week_skin"
	mds1$induction <- "8_week_skin"
	mds2$induction <- "8_week_skin"
	mds3$induction <- "8_week_skin"
	mds4$induction <- "8_week_skin"
	fbm1$induction <- "oral"
	fbm2$induction <- "oral"
	fbm3$induction <- "oral"
	fbm4$induction <- "oral"
	mbm1$induction <- "oral"
	mbm2$induction <- "oral"
	mbm3$induction <- "oral"
	
	#Add skin state
	fcs1$state <- "skin_control_8_week"
	fcs2$state <- "skin_control_8_week"
	fcs3$state <- "skin_control_8_week"
	fcs4$state <- "skin_control_8_week"
	fcs5$state <- "skin_control_8_week"
	mcs1$state <- "skin_control_8_week"
	mcs2$state <- "skin_control_8_week"
	mcs3$state <- "skin_control_8_week"
	fds1$state <- "skin_dox_8_week"
	fds2$state <- "skin_dox_8_week"
	fds3$state <- "skin_dox_8_week"
	fds4$state <- "skin_dox_8_week"
	mds1$state <- "skin_dox_8_week"
	mds2$state <- "skin_dox_8_week"
	mds3$state <- "skin_dox_8_week"
	mds4$state <- "skin_dox_8_week"
	fbm1$state <- "control_buccal_mucosa"
	fbm2$state <- "control_buccal_mucosa"
	fbm3$state <- "control_buccal_mucosa"
	fbm4$state <- "control_buccal_mucosa"
	mbm1$state <- "control_buccal_mucosa"
	mbm2$state <- "control_buccal_mucosa"
	mbm3$state <- "control_buccal_mucosa"
	
	#Add sex identity
	fcs1$sex <- "female"
	fcs2$sex <- "female"
	fcs3$sex <- "female"
	fcs4$sex <- "female"
	fcs5$sex <- "female"
	mcs1$sex <- "male"
	mcs2$sex <- "male"
	mcs3$sex <- "male"
	fds1$sex <- "female"
	fds2$sex <- "female"
	fds3$sex <- "female"
	fds4$sex <- "female"
	mds1$sex <- "male"
	mds2$sex <- "male"
	mds3$sex <- "male"
	mds4$sex <- "male"
	fbm1$sex <- "female"
	fbm2$sex <- "female"
	fbm3$sex <- "female"
	fbm4$sex <- "female"
	mbm1$sex <- "male"
	mbm2$sex <- "male"
	mbm3$sex <- "male"
	
	#Add batch #
	fcs1$batch <- "1"
	fcs2$batch <- "1"
	fcs3$batch <- "1"
	fcs4$batch <- "2"
	fcs5$batch <- "2"
	mcs1$batch <- "1"
	mcs2$batch <- "2"
	mcs3$batch <- "2"
	fds1$batch <- "1"
	fds2$batch <- "1"
	fds3$batch <- "2"
	fds4$batch <- "2"
	mds1$batch <- "1"
	mds2$batch <- "1"
	mds3$batch <- "2"
	mds4$batch <- "2"
	fbm1$batch <- "1"
	fbm2$batch <- "2"
	fbm3$batch <- "3"
	fbm4$batch <- "4"
	mbm1$batch <- "2"
	mbm2$batch <- "3"
	mbm3$batch <- "4"
	
	#Add ID#
	fcs1$sampleid <- "fcs1"
	fcs2$sampleid <- "fcs2"
	fcs3$sampleid <- "fcs3"
	fcs4$sampleid <- "fcs4"
	fcs5$sampleid <- "fcs5"
	mcs1$sampleid <- "mcs1"
	mcs2$sampleid <- "mcs2"
	mcs3$sampleid <- "mcs3"
	fds1$sampleid <- "fds1"
	fds2$sampleid <- "fds2"
	fds3$sampleid <- "fds3"
	fds4$sampleid <- "fds4"
	mds1$sampleid <- "mds1"
	mds2$sampleid <- "mds2"
	mds3$sampleid <- "mds3"
	mds4$sampleid <- "mds4"
	fbm1$sampleid <- "fbm1"
	fbm2$sampleid <- "fbm2"
	fbm3$sampleid <- "fbm3"
	fbm4$sampleid <- "fbm4"
	mbm1$sampleid <- "mbm1"
	mbm2$sampleid <- "mbm2"
	mbm3$sampleid <- "mbm3"
	
	#Add %mitochondrial gene to metadata, filter based on nFeature_RNA & percent.mt, normalize data, 
	fcs1[["percent.mt"]] <- PercentageFeatureSet(fcs1, pattern = "^mt-")
	fcs2[["percent.mt"]] <- PercentageFeatureSet(fcs2, pattern = "^mt-")
	fcs3[["percent.mt"]] <- PercentageFeatureSet(fcs3, pattern = "^mt-")
	fcs4[["percent.mt"]] <- PercentageFeatureSet(fcs4, pattern = "^mt-")
	fcs5[["percent.mt"]] <- PercentageFeatureSet(fcs5, pattern = "^mt-")
	mcs1[["percent.mt"]] <- PercentageFeatureSet(mcs1, pattern = "^mt-")
	mcs2[["percent.mt"]] <- PercentageFeatureSet(mcs2, pattern = "^mt-")
	mcs3[["percent.mt"]] <- PercentageFeatureSet(mcs3, pattern = "^mt-")
	fds1[["percent.mt"]] <- PercentageFeatureSet(fds1, pattern = "^mt-")
	fds2[["percent.mt"]] <- PercentageFeatureSet(fds2, pattern = "^mt-")
	fds3[["percent.mt"]] <- PercentageFeatureSet(fds3, pattern = "^mt-")
	fds4[["percent.mt"]] <- PercentageFeatureSet(fds4, pattern = "^mt-")
	mds1[["percent.mt"]] <- PercentageFeatureSet(mds1, pattern = "^mt-")
	mds2[["percent.mt"]] <- PercentageFeatureSet(mds2, pattern = "^mt-")
	mds3[["percent.mt"]] <- PercentageFeatureSet(mds3, pattern = "^mt-")
	mds4[["percent.mt"]] <- PercentageFeatureSet(mds4, pattern = "^mt-")
	fbm1[["percent.mt"]] <- PercentageFeatureSet(fbm1, pattern = "^mt-")
	fbm2[["percent.mt"]] <- PercentageFeatureSet(fbm2, pattern = "^mt-")
	fbm3[["percent.mt"]] <- PercentageFeatureSet(fbm3, pattern = "^mt-")
	fbm4[["percent.mt"]] <- PercentageFeatureSet(fbm4, pattern = "^mt-")
	mbm1[["percent.mt"]] <- PercentageFeatureSet(mbm1, pattern = "^mt-")
	mbm2[["percent.mt"]] <- PercentageFeatureSet(mbm2, pattern = "^mt-")
	mbm3[["percent.mt"]] <- PercentageFeatureSet(mbm3, pattern = "^mt-")
	
	fcs1[["percent.ribo"]] <- PercentageFeatureSet(fcs1, pattern = "Rp[sl]")
	fcs2[["percent.ribo"]] <- PercentageFeatureSet(fcs2, pattern = "Rp[sl]")
	fcs3[["percent.ribo"]] <- PercentageFeatureSet(fcs3, pattern = "Rp[sl]")
	fcs4[["percent.ribo"]] <- PercentageFeatureSet(fcs4, pattern = "Rp[sl]")
	fcs5[["percent.ribo"]] <- PercentageFeatureSet(fcs5, pattern = "Rp[sl]")
	mcs1[["percent.ribo"]] <- PercentageFeatureSet(mcs1, pattern = "Rp[sl]")
	mcs2[["percent.ribo"]] <- PercentageFeatureSet(mcs2, pattern = "Rp[sl]")
	mcs3[["percent.ribo"]] <- PercentageFeatureSet(mcs3, pattern = "Rp[sl]")
	fds1[["percent.ribo"]] <- PercentageFeatureSet(fds1, pattern = "Rp[sl]")
	fds2[["percent.ribo"]] <- PercentageFeatureSet(fds2, pattern = "Rp[sl]")
	fds3[["percent.ribo"]] <- PercentageFeatureSet(fds3, pattern = "Rp[sl]")
	fds4[["percent.ribo"]] <- PercentageFeatureSet(fds4, pattern = "Rp[sl]")
	mds1[["percent.ribo"]] <- PercentageFeatureSet(mds1, pattern = "Rp[sl]")
	mds2[["percent.ribo"]] <- PercentageFeatureSet(mds2, pattern = "Rp[sl]")
	mds3[["percent.ribo"]] <- PercentageFeatureSet(mds3, pattern = "Rp[sl]")
	mds4[["percent.ribo"]] <- PercentageFeatureSet(mds4, pattern = "Rp[sl]")
	fbm1[["percent.ribo"]] <- PercentageFeatureSet(fbm1, pattern = "Rp[sl]")
	fbm2[["percent.ribo"]] <- PercentageFeatureSet(fbm2, pattern = "Rp[sl]")
	fbm3[["percent.ribo"]] <- PercentageFeatureSet(fbm3, pattern = "Rp[sl]")
	fbm4[["percent.ribo"]] <- PercentageFeatureSet(fbm4, pattern = "Rp[sl]")
	mbm1[["percent.ribo"]] <- PercentageFeatureSet(mbm1, pattern = "Rp[sl]")
	mbm2[["percent.ribo"]] <- PercentageFeatureSet(mbm2, pattern = "Rp[sl]")
	mbm3[["percent.ribo"]] <- PercentageFeatureSet(mbm3, pattern = "Rp[sl]")
	
	fcs1 <- subset(fcs1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fcs2 <- subset(fcs2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5)
	fcs3 <- subset(fcs3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5)
	fcs4 <- subset(fcs4, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fcs5 <- subset(fcs5, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mcs1 <- subset(mcs1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mcs2 <- subset(mcs2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mcs3 <- subset(mcs3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds1 <- subset(fds1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds2 <- subset(fds2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds3 <- subset(fds3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds4 <- subset(fds4, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds1 <- subset(mds1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds2 <- subset(mds2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds3 <- subset(mds3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds4 <- subset(mds4, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fbm1 <- subset(fbm1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fbm2 <- subset(fbm2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5)
	fbm3 <- subset(fbm3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5)
	fbm4 <- subset(fbm4, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mbm1 <- subset(mbm1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mbm2 <- subset(mbm2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mbm3 <- subset(mbm3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	
	tissue.list <- list(fcs1 = fcs1, fcs2 = fcs2, fcs3 = fcs3, fcs4 = fcs4, fcs5 = fcs5,
	                    mcs1 = mcs1, mcs2 = mcs2, mcs3 = mcs3,
	                    fds1 = fds1, fds2 = fds2, fds3 = fds3, fds4 = fds4,
	                    mds1 = mds1, mds2 = mds2, mds3 = mds3, mds4 = mds4,
	                    fbm1 = fbm1, fbm2 = fbm2, fbm3 = fbm3, fbm4 = fbm4, 
	                    mbm1 = mbm1, mbm2 = mbm2, mbm3 = mbm3)
	
	saveRDS(tissue.list,
	        file = "healthytissue.list.rds",
	        compress = T)
	
	tissue.list <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.list.rds')
	
	#File cleanup
	rm(fcs1, fcs2, fcs3, fcs4, fcs5,
	   mcs1, mcs2, mcs3,
	   fds1, fds2, fds3, fds4,
	   mds1, mds2, mds3, mds4,
	   fbm1, fbm2, fbm3, fbm4,
	   mbm1, mbm2, mbm3)
	
	########
	
	tissue.list <- c(readRDS('/data/overmilleram/scRNAseq/Skin/healthyskinlist.rds'), 
	                 readRDS('/data/overmilleram/scRNAseq/Oral/healthyoral.list.rds'))
	
	#Merge Seurat objects to proceed with CellCycleScoring
	tissue.merge <- merge(x = tissue.list[[1]],
	                      y = tissue.list[-1],
	                      add.cell.ids = names(tissue.list),
	                      merge.data = F,
	                      project = 'skin_oral')
	
	tissue.merge <- NormalizeData(tissue.merge)
	
	tissue.merge <- JoinLayers(tissue.merge) # need to join layers as CellCycleScoring can't do separate layers (GetAssayData() doesn't work for multiple layers in v5 assay)
	
	m.s.genes <- cc.genes$s.genes %>% convert_human_to_mouse_symbols
	m.g2m.genes <- cc.genes$g2m.genes %>% convert_human_to_mouse_symbols
	
	tissue.merge <- CellCycleScoring(tissue.merge, 
	                                 s.features = m.s.genes,
	                                 g2m.features = m.g2m.genes, 
	                                 set.ident = T)
	
	tissue.merge[['CC.Diff']] <- tissue.merge[['S.Score']] - tissue.merge[['G2M.Score']]
	
	tissue.merge[['RNA']] <- split(tissue.merge[['RNA']], f = tissue.merge$sampleid)
	
	tissue.merge <- FindVariableFeatures(tissue.merge, 
	                                     selection.method = 'vst', 
	                                     nfeatures = 5000)
	
	tissue.merge <- SCTransform(tissue.merge,
	                            method = 'glmGamPoi',
	                            vars.to.regress = 'CC.Diff',
	                            vst.flavor = 'v2',
	                            ncells = length(colnames(tissue.merge))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.merge <- RunPCA(tissue.merge, npcs = 100)
	
	saveRDS(tissue.merge,
	        '/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sct.merge.rds',
	        compress = F)
	
	tissue.integrated <- IntegrateLayers(tissue.merge, 
	                                     method = HarmonyIntegration, 
	                                     orig.reduction = 'pca',
	                                     new.reduction = 'harmony',
	                                     assay = 'SCT',
	                                     theta = 4,
	                                     sigma = 0.1,# variable that increasing diversity of clusters as increased (i.e., higher theta = more, unique clusters)
	                                     max.iter.cluster = 100,
	                                     npcs = 100, 
	                                     verbose = T)
	
	tissue.integrated <- FindNeighbors(tissue.integrated, 
	                                   reduction = 'harmony', 
	                                   assay = 'SCT',
	                                   dims = 1:100)
	
	tissue.integrated <- RunUMAP(tissue.integrated, 
	                              reduction = 'harmony', 
	                              reduction.name = 'umap.harmony',
	                              assay = 'SCT',
	                              dims = 1:100)
	
	tissue.integrated <- FindClusters(tissue.integrated, 
	                                  method = 'igraph',
	                                  algorithm = 4,
	                                  resolution = 0.8)
	
	DimPlot(tissue.integrated,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 0.75,
	        reduction = 'umap.harmony',
	        label.size = 8,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Harmony2 UMAP.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	FeaturePlot(tissue.integrated,
	            features = 'Mlana',
	            cols = c('grey', 'red'),
	            min.cutoff = 'q10',
	            order = T,
	            raster = F,
	            pt.size = 1)
	
	rm(tissue.merge)
	
	tissue.integrated <- JoinLayers(tissue.integrated, assay = 'RNA')
	
	saveRDS(tissue.integrated,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony2.integrated.rds",
	        compress = FALSE)
	
	# find cluster markers
	
	tissue.integrated.markers <- FindAllMarkers(tissue.integrated,
	                                            only.pos = T, 
	                                            test.use = "MAST",
	                                            latent.vars = "sex",
	                                            min.pct = 0.50,
	                                            logfc.threshold = 2.00,
	                                            return.thresh = 0.05,
	                                            assay = "RNA",
	                                            densify = T) 
	
	tissue.integrated.cellcounts <- table(tissue.integrated@meta.data$seurat_clusters,
	                                      tissue.integrated@meta.data$state)
	
	tissue.integrated.cellcounts2 <- table(tissue.integrated@meta.data$seurat_clusters,
	                                       tissue.integrated@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.integrated.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           row.names = F,
	           append = FALSE)
	write.xlsx(tissue.integrated.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.integrated.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	########
	
	p1 <- VlnPlot(tissue.integrated, 
	        features = c('Krt5', 'Sbsn', 'Krt79', 'Krt4', 'Tgm3', 'Muc19', 'Bpifb2', 'Col1a1', 'Dcn', 'Lum', 'Ptprc', 
	                     'Cd3e', 'G0s2', 'Pecam1', 'Ccl21a', 'Myh11', 'Mpz', 'Mbp', 'Scn7a', 'Des', 'Dmd', 'Myl1', 'Tnnt3', 'Mylpf', 'Mlana', 'Pmel', 'Kit',
	                     'Hba-a1', 'Hbb-bs', 'percent.mt', 'percent.ribo', 'Lars2'), 
	        stack = T,
	        flip = T,
	        pt.size = 0,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(p1,filename = "Healthy Tissue Harmony2 Celltype VlnPlot2.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	table(tissue.integrated$SCT_snn_res.0.8 %in% c('46'), tissue.integrated$sampleid) # check '31','32','33',37','41','44','46'
	
	tissue.integrated <- subset(tissue.integrated, idents = c('31','32','33','37','41','44','46'), invert = T) # RBC, RBC-contaminated, or distinct doublet clusters
	
	table(tissue.integrated$scvi.clusters, tissue.integrated$sampleid) # check that individual samples do not dominate single clusters
	
	tissue.integrated$seurat_clusters <- droplevels(tissue.integrated$seurat_clusters)
	tissue.integrated$SCT_snn_res.0.8 <- droplevels(tissue.integrated$SCT_snn_res.0.8)
	
	tissue.integrated <- RenameIdents(tissue.integrated,
	                                  '1' =  'Immune',
	                                  '2' =  'Fibroblast',
	                                  '3' =  'Immune',
	                                  '4' =  'Keratinocyte',
	                                  '5' =  'Fibroblast',
	                                  '6' =  'Keratinocyte',
	                                  '7' =  'Keratinocyte',
	                                  '8' =  'Keratinocyte',
	                                  '9' =  'Keratinocyte',
	                                  '10' = 'Keratinocyte',
	                                  '11' = 'Fibroblast',
	                                  '12' = 'Keratinocyte',
	                                  '13' = 'Immune',
	                                  '14' = 'Fibroblast',
	                                  '15' = 'Keratinocyte',
	                                  '16' = 'Keratinocyte',
	                                  '17' = 'Vascular',
	                                  '18' = 'Salivary',
	                                  '19' = 'Immune',
	                                  '20' = 'Keratinocyte',
	                                  '21' = 'Keratinocyte',
	                                  '22' = 'Keratinocyte',
	                                  '23' = 'Keratinocyte',
	                                  '24' = 'Keratinocyte',
	                                  '25' = 'Keratinocyte',
	                                  '26' = 'Keratinocyte',
	                                  '27' = 'Keratinocyte',
	                                  '28' = 'Keratinocyte',
	                                  '29' = 'Keratinocyte',
	                                  '30' = 'Keratinocyte',
	                                  '34' = 'Salivary',
	                                  '35' = 'Salivary',
	                                  '36' = 'Fibroblast',
	                                  '38' = 'Vascular',
	                                  '39' = 'Neural',
	                                  '40' = 'Neural',
	                                  '42' = 'Keratinocyte',
	                                  '43' = 'Mesenchymal',
	                                  '45' = 'Mesenchymal',
	                                  '47' = 'Salivary',
	                                  '48' = 'Keratinocyte',
	                                  '49' = 'Salivary',
	                                  '50' = 'Salivary')
	
	levels(tissue.integrated) <- c("Keratinocyte", 
	                               "Immune",
	                               "Fibroblast",
	                               "Vascular",
	                               "Neural",
	                               'Mesenchymal',
	                               'Salivary')
	
	tissue.integrated$celltype = tissue.integrated@active.ident
	
	Idents(tissue.integrated) <- 'celltype'
	
	umap.plot <- DimPlot(tissue.integrated,
	                     label = TRUE,
	                     repel = TRUE,
	                     pt.size = 0.75,
	                     reduction = 'umap.harmony',
	                     label.size = 10,
	                     raster = FALSE) + NoLegend()
	ggsave2(umap.plot,
	        filename = "Healthy Tissue Harmony2 Celltype UMAP.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.integrated) <- 'seurat_clusters'
	
	umap.plot <- DimPlot(tissue.integrated,
	                     label = TRUE,
	                     repel = TRUE,
	                     pt.size = 0.75,
	                     reduction = 'umap.harmony',
	                     label.size = 8,
	                     raster = FALSE) + NoLegend()
	ggsave2(umap.plot,
	        filename = "Healthy Tissue Harmony2 Cluster UMAP.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	saveRDS(tissue.integrated,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony2.integrated.rds",
	        compress = FALSE)
	
	tissue.integrated <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony2.integrated.rds')
	
	Idents(tissue.integrated) <- 'celltype'
	
	tissue.kera <- subset(tissue.integrated,
	                       idents = "Keratinocyte")
	tissue.imm <- subset(tissue.integrated,
	                      idents = "Immune")
	tissue.fib <- subset(tissue.integrated,
	                      idents = "Fibroblast")
	tissue.vasc <- subset(tissue.integrated,
	                       idents = "Vascular")
	tissue.neur <- subset(tissue.integrated,
	                       idents = "Neural")
	tissue.mes <- subset(tissue.integrated,
	                      idents = "Mesenchymal")
	tissue.sal <- subset(tissue.integrated,
	                     idents = 'Salivary')
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	# Keratinocyte
	
	tissue.kera <- SCTransform(tissue.kera,
	                           method = 'glmGamPoi',
	                           vars.to.regress = 'CC.Diff',
	                           vst.flavor = 'v2',
	                           do.scale = T,
	                           ncells = length(colnames(tissue.kera))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.kera <- RunPCA(tissue.kera, 
	                      npcs = 100, 
	                      assay = 'SCT')
	
	ElbowPlot(tissue.kera, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.kera <- FindNeighbors(tissue.kera, 
	                             reduction = 'harmony', 
	                             dims = 1:100,
	                             assay = 'SCT')
	
	tissue.kera2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.kera, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'kera.res',
	                                                                   resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.kera)
	
	tissue.kera3 <- data.frame(as.numeric(as.character(tissue.kera2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.kera3[ ,i] <- data.frame(as.numeric(as.character(tissue.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.kera3) <- col.names
	rownames(tissue.kera3) <- row.names
	
	rm(tissue.kera2)
	
	gc()
	
	tissue.kera2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.kera, 
	                                                                    method = 'igraph',
	                                                                    algorithm = 4,
	                                                                    cluster.name = 'kera.res',
	                                                                    resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.kera)
	
	tissue.kera4 <- data.frame(as.numeric(as.character(tissue.kera2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.kera4[ ,i] <- data.frame(as.numeric(as.character(tissue.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.kera4) <- col.names
	rownames(tissue.kera4) <- row.names
	
	rm(tissue.kera2)
	
	gc()
	
	tissue.kera2 <- foreach(n = seq(0.7,0.9,0.1)) %dopar% FindClusters(tissue.kera, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'kera.res',
	                                                                   resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,0.9,0.1))
	row.names <- colnames(tissue.kera)
	
	tissue.kera5 <- data.frame(as.numeric(as.character(tissue.kera2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.kera5[ ,i] <- data.frame(as.numeric(as.character(tissue.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.kera5) <- col.names
	rownames(tissue.kera5) <- row.names
	
	tissue.kera.meta <- data.frame(tissue.kera3, tissue.kera4, tissue.kera5)
	
	tissue.kera <- AddMetaData(tissue.kera, 
	                            tissue.kera.meta)
	
	rm(tissue.kera2)
	
	gc()
	
	tissue.kera <- RunUMAP(tissue.kera, 
	                       reduction = 'harmony', 
	                       reduction.name = 'umap.harmony',
	                       dims = 1:100,
	                       assay = 'SCT')
	
	saveRDS(tissue.kera,
	         file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.kera, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Immune
	
	tissue.imm <- SCTransform(tissue.imm,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.imm))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.imm <- RunPCA(tissue.imm, 
	                      npcs = 100, 
	                      assay = 'SCT')
	
	ElbowPlot(tissue.imm, ndims = 100) #decided to do 80 PCs for this analysis
	
	tissue.imm <- FindNeighbors(tissue.imm, 
	                            reduction = 'harmony', 
	                            dims = 1:80,
	                            assay = 'SCT')
	
	tissue.imm2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.imm, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'imm.res',
	                                                                   resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.imm)
	
	tissue.imm3 <- data.frame(as.numeric(as.character(tissue.imm2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.imm3[ ,i] <- data.frame(as.numeric(as.character(tissue.imm2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.imm3) <- col.names
	rownames(tissue.imm3) <- row.names
	
	tissue.imm <- AddMetaData(tissue.imm, 
	                           tissue.imm3)
	
	rm(tissue.imm2, tissue.imm3)
	
	gc()
	
	tissue.imm <- RunUMAP(tissue.imm, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:80,
	                      assay = 'SCT')
	
	saveRDS(tissue.imm,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.imm, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Fibroblast
	
	tissue.fib <- SCTransform(tissue.fib,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.fib))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.fib <- RunPCA(tissue.fib, 
	                      npcs = 100, 
	                      assay = 'SCT')
	
	ElbowPlot(tissue.fib, ndims = 100) #decided to do 80 PCs for this analysis
	
	tissue.fib <- FindNeighbors(tissue.fib, 
	                            reduction = 'harmony', 
	                            dims = 1:80,
	                            assay = 'SCT')
	
	tissue.fib2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.fib, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'fib.res',
	                                                                   resolution = n)
	col.names <- paste0('fib.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.fib)
	
	tissue.fib3 <- data.frame(as.numeric(as.character(tissue.fib2[[1]]@meta.data$fib.res)))
	
	for (i in 1:10) {
	  tissue.fib3[ ,i] <- data.frame(as.numeric(as.character(tissue.fib2[[i]]@meta.data$fib.res)))
	}  
	
	colnames(tissue.fib3) <- col.names
	rownames(tissue.fib3) <- row.names
	
	tissue.fib <- AddMetaData(tissue.fib, 
	                          tissue.fib3)
	
	rm(tissue.fib2, tissue.fib3)
	
	gc()
	
	tissue.fib <- RunUMAP(tissue.fib, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:80,
	                      assay = 'SCT')
	
	saveRDS(tissue.fib,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.fib, prefix = 'fib.res.')
	ggsave2(filename = "Healthy Tissue Fibroblast Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Vascular
	
	tissue.vasc <- SCTransform(tissue.vasc,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.vasc))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.vasc <- RunPCA(tissue.vasc, 
	                      npcs = 100, 
	                      assay = 'SCT')
	
	ElbowPlot(tissue.vasc, ndims = 100) #decided to do 60 PCs for this analysis
	
	tissue.vasc <- FindNeighbors(tissue.vasc, 
	                             reduction = 'harmony', 
	                             dims = 1:60,
	                             assay = 'SCT')
	
	tissue.vasc2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.vasc, 
	                                                                    method = 'igraph',
	                                                                    algorithm = 4,
	                                                                    cluster.name = 'vasc.res',
	                                                                    resolution = n)
	col.names <- paste0('vasc.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.vasc)
	
	tissue.vasc3 <- data.frame(as.numeric(as.character(tissue.vasc2[[1]]@meta.data$vasc.res)))
	
	for (i in 1:10) {
	  tissue.vasc3[ ,i] <- data.frame(as.numeric(as.character(tissue.vasc2[[i]]@meta.data$vasc.res)))
	}  
	
	colnames(tissue.vasc3) <- col.names
	rownames(tissue.vasc3) <- row.names
	
	tissue.vasc <- AddMetaData(tissue.vasc, 
	                            tissue.vasc3)
	
	rm(tissue.vasc2, tissue.vasc3)
	
	gc()
	
	tissue.vasc <- RunUMAP(tissue.vasc, 
	                       reduction = 'harmony', 
	                       reduction.name = 'umap.harmony',
	                       dims = 1:60,
	                       assay = 'SCT')
	
	saveRDS(tissue.vasc,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.vasc, prefix = 'vasc.res.')
	ggsave2(filename = "Healthy Tissue Vascular Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Neural
	
	tissue.neur <- SCTransform(tissue.neur,
	                           method = 'glmGamPoi',
	                           vars.to.regress = 'CC.Diff',
	                           vst.flavor = 'v2',
	                           do.scale = T,
	                           ncells = length(colnames(tissue.neur))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.neur <- RunPCA(tissue.neur, 
	                      npcs = 100, 
	                      assay = 'SCT')
	
	ElbowPlot(tissue.neur, ndims = 100) #decided to do 50 PCs for this analysis
	
	tissue.neur <- FindNeighbors(tissue.neur, 
	                             reduction = 'harmony', 
	                             dims = 1:50,
	                             assay = 'SCT')
	
	tissue.neur2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.neur, 
	                                                                    method = 'igraph',
	                                                                    algorithm = 4,
	                                                                    cluster.name = 'neur.res',
	                                                                    resolution = n)
	col.names <- paste0('neur.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.neur)
	
	tissue.neur3 <- data.frame(as.numeric(as.character(tissue.neur2[[1]]@meta.data$neur.res)))
	
	for (i in 1:10) {
	  tissue.neur3[ ,i] <- data.frame(as.numeric(as.character(tissue.neur2[[i]]@meta.data$neur.res)))
	}  
	
	colnames(tissue.neur3) <- col.names
	rownames(tissue.neur3) <- row.names
	
	tissue.neur <- AddMetaData(tissue.neur, 
	                            tissue.neur3)
	
	rm(tissue.neur2, tissue.neur3)
	
	gc()
	
	tissue.neur <- RunUMAP(tissue.neur, 
	                       reduction = 'harmony', 
	                       reduction.name = 'umap.harmony',
	                       dims = 1:50,
	                       assay = 'SCT')
	
	saveRDS(tissue.neur,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.neur, prefix = 'neur.res.')
	ggsave2(filename = "Healthy Tissue Neural Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Mesenchymal
	
	tissue.mes <- SCTransform(tissue.mes,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.mes))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.mes <- RunPCA(tissue.mes, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.mes, ndims = 100) #decided to do 50 PCs for this analysis
	
	tissue.mes <- FindNeighbors(tissue.mes, 
	                            reduction = 'harmony', 
	                            dims = 1:50,
	                            assay = 'SCT')
	
	tissue.mes2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.mes, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'mes.res',
	                                                                  resolution = n)
	col.names <- paste0('mes.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.mes)
	
	tissue.mes3 <- data.frame(as.numeric(as.character(tissue.mes2[[1]]@meta.data$mes.res)))
	
	for (i in 1:10) {
	  tissue.mes3[ ,i] <- data.frame(as.numeric(as.character(tissue.mes2[[i]]@meta.data$mes.res)))
	}  
	
	colnames(tissue.mes3) <- col.names
	rownames(tissue.mes3) <- row.names
	
	tissue.mes <- AddMetaData(tissue.mes, 
	                          tissue.mes3)
	
	rm(tissue.mes2, tissue.mes3)
	
	gc()
	
	tissue.mes <- RunUMAP(tissue.mes, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:50,
	                      assay = 'SCT')
	
	saveRDS(tissue.mes,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.mes, prefix = 'mes.res.')
	ggsave2(filename = "Healthy Tissue Mesenchymal Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Salivary
	
	tissue.sal <- SCTransform(tissue.sal,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.sal))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.sal <- RunPCA(tissue.sal, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.sal, ndims = 100) #decided to do 60 PCs for this analysis
	
	tissue.sal <- FindNeighbors(tissue.sal, 
	                            reduction = 'harmony', 
	                            dims = 1:60,
	                            assay = 'SCT')
	
	tissue.sal2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.sal, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'sal.res',
	                                                                  resolution = n)
	col.names <- paste0('sal.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.sal)
	
	tissue.sal3 <- data.frame(as.numeric(as.character(tissue.sal2[[1]]@meta.data$sal.res)))
	
	for (i in 1:10) {
	  tissue.sal3[ ,i] <- data.frame(as.numeric(as.character(tissue.sal2[[i]]@meta.data$sal.res)))
	}  
	
	colnames(tissue.sal3) <- col.names
	rownames(tissue.sal3) <- row.names
	
	tissue.sal <- AddMetaData(tissue.sal, 
	                          tissue.sal3)
	
	rm(tissue.sal2, tissue.sal3)
	
	gc()
	
	tissue.sal <- RunUMAP(tissue.sal, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:60,
	                      assay = 'SCT')
	
	saveRDS(tissue.sal,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.sal, prefix = 'sal.res.')
	ggsave2(filename = "Healthy Tissue Salivary Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.kera <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds")
	tissue.fib <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds")
	tissue.imm <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds")
	tissue.vasc <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds")
	tissue.neur <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds")
	tissue.mes <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds")
	tissue.sal <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds")
	
	Idents(tissue.kera) <- "kera.res.0.9"
	Idents(tissue.imm) <-  "imm.res.0.4"
	Idents(tissue.fib) <-  "fib.res.0.6"
	Idents(tissue.vasc) <- "vasc.res.0.1"
	Idents(tissue.neur) <- "neur.res.0.1"
	Idents(tissue.mes) <-  "mes.res.0.1"
	Idents(tissue.sal) <- 'sal.res.0.4'
	
	umap.kera <- DimPlot(tissue.kera,
	                     label = TRUE,
	                     repel = TRUE,
	                     pt.size = 0.5,
	                     label.size = 10,
	                     raster = FALSE) + NoLegend()
	
	umap.kera
	
	umap.imm <- DimPlot(tissue.imm,
	                    label = TRUE,
	                    repel = TRUE,
	                    pt.size = 1,
	                    label.size = 10,
	                    raster = FALSE) + NoLegend()
	
	umap.imm
	
	umap.fib <- DimPlot(tissue.fib,
	                    label = TRUE,
	                    repel = TRUE,
	                    pt.size = 1,
	                    label.size = 10,
	                    raster = FALSE) + NoLegend()
	
	umap.fib
	
	umap.vasc <- DimPlot(tissue.vasc,
	                     label = TRUE,
	                     repel = TRUE,
	                     pt.size = 2,
	                     label.size = 10,
	                     raster = FALSE) + NoLegend()
	
	umap.vasc
	
	umap.neur <- DimPlot(tissue.neur,
	                     label = TRUE,
	                     repel = TRUE,
	                     pt.size = 4,
	                     label.size = 10,
	                     raster = FALSE) + NoLegend()
	
	umap.neur
	
	umap.mes <- DimPlot(tissue.mes,
	                    label = TRUE,
	                    repel = TRUE,
	                    pt.size = 4,
	                    label.size = 10,
	                    raster = FALSE) + NoLegend()
	
	umap.mes
	
	umap.sal <- DimPlot(tissue.sal,
	                    label = TRUE,
	                    repel = TRUE,
	                    pt.size = 2,
	                    label.size = 10,
	                    raster = FALSE) + NoLegend()
	
	umap.sal
	
	
	ggsave2(umap.kera,
	        filename = "Healthy Tissue Keratinocyte UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(umap.imm,
	        filename = "Healthy Tissue Immune UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(umap.fib,
	        filename = "Healthy Tissue Fibroblast UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(umap.vasc,
	        filename = "Healthy Tissue Vascular UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(umap.neur,
	        filename = "Healthy Tissue Neural UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(umap.sal,
	        filename = "Healthy Tissue Salivary UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.kera) <- "kera.res.0.9"
	Idents(tissue.imm) <-  "imm.res.0.4"
	Idents(tissue.fib) <-  "fib.res.0.6"
	Idents(tissue.vasc) <- "vasc.res.0.1"
	Idents(tissue.neur) <- "neur.res.0.1"
	Idents(tissue.mes) <-  "mes.res.0.1"
	Idents(tissue.sal) <- 'sal.res.0.4'
	
	tissue.kera.markers <- FindAllMarkers(tissue.kera,
	                                      only.pos = T, 
	                                      test.use = "MAST",
	                                      latent.vars = "sex",
	                                      min.pct = 0.50,
	                                      logfc.threshold = 1.00,
	                                      return.thresh = 0.05,
	                                      assay = "RNA",
	                                      densify = T) 
	
	tissue.kera.markers <- tissue.kera.markers[order(tissue.kera.markers$cluster,
	                                                 -tissue.kera.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.kera.cellcounts <- table(tissue.kera@meta.data$kera.res.0.9,
	                                tissue.kera@meta.data$state)
	
	tissue.kera.cellcounts2 <- table(tissue.kera@meta.data$kera.res.0.9,
	                                 tissue.kera@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.kera.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.kera.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.kera.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.imm.markers <- FindAllMarkers(tissue.imm,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.imm.markers <- tissue.imm.markers[order(tissue.imm.markers$cluster,
	                                               -tissue.imm.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.imm.cellcounts <- table(tissue.imm@meta.data$imm.res.0.4,
	                               tissue.imm@meta.data$state)
	
	tissue.imm.cellcounts2 <- table(tissue.imm@meta.data$imm.res.0.4,
	                                tissue.imm@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.imm.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.imm.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.imm.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.fib.markers <- FindAllMarkers(tissue.fib,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.fib.markers <- tissue.fib.markers[order(tissue.fib.markers$cluster,
	                                               -tissue.fib.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.fib.cellcounts <- table(tissue.fib@meta.data$fib.res.0.6,
	                               tissue.fib@meta.data$state)
	
	tissue.fib.cellcounts2 <- table(tissue.fib@meta.data$fib.res.0.6,
	                                tissue.fib@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.fib.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.fib.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.fib.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.vasc.markers <- FindAllMarkers(tissue.vasc,
	                                      only.pos = T, 
	                                      test.use = "MAST",
	                                      latent.vars = "sex",
	                                      min.pct = 0.50,
	                                      logfc.threshold = 1.00,
	                                      return.thresh = 0.05,
	                                      assay = "RNA",
	                                      densify = T) 
	
	tissue.vasc.markers <- tissue.vasc.markers[order(tissue.vasc.markers$cluster,
	                                                 -tissue.vasc.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.vasc.cellcounts <- table(tissue.vasc@meta.data$vasc.res.0.1,
	                                tissue.vasc@meta.data$state)
	
	tissue.vasc.cellcounts2 <- table(tissue.vasc@meta.data$vasc.res.0.1,
	                                 tissue.vasc@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.vasc.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Vascular Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.vasc.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Vascular Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.vasc.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Vascular Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.neur.markers <- FindAllMarkers(tissue.neur,
	                                      only.pos = T, 
	                                      test.use = "MAST",
	                                      latent.vars = "sex",
	                                      min.pct = 0.50,
	                                      logfc.threshold = 1.00,
	                                      return.thresh = 0.05,
	                                      assay = "RNA",
	                                      densify = T) 
	
	tissue.neur.markers <- tissue.neur.markers[order(tissue.neur.markers$cluster,
	                                                 -tissue.neur.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.neur.cellcounts <- table(tissue.neur@meta.data$neur.res.0.1,
	                                tissue.neur@meta.data$state)
	
	tissue.neur.cellcounts2 <- table(tissue.neur@meta.data$neur.res.0.1,
	                                 tissue.neur@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.neur.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Neural Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.neur.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Neural Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.neur.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Neural Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.mes.markers <- FindAllMarkers(tissue.mes,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.mes.markers <- tissue.mes.markers[order(tissue.mes.markers$cluster,
	                                               -tissue.mes.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.mes.cellcounts <- table(tissue.mes@meta.data$mes.res.0.1,
	                               tissue.mes@meta.data$state)
	
	tissue.mes.cellcounts2 <- table(tissue.mes@meta.data$mes.res.0.1,
	                                tissue.mes@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.mes.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Mesenchymal Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.mes.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Mesenchymal Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.mes.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Mesenchymal Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.sal.markers <- FindAllMarkers(tissue.sal,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.sal.markers <- tissue.sal.markers[order(tissue.sal.markers$cluster,
	                                               -tissue.sal.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.sal.cellcounts <- table(tissue.sal@meta.data$sal.res.0.4,
	                               tissue.sal@meta.data$state)
	
	tissue.sal.cellcounts2 <- table(tissue.sal@meta.data$sal.res.0.4,
	                                tissue.sal@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.sal.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Salivary Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.sal.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Salivary Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.sal.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Salivary Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.kera <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds")
	tissue.fib <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds")
	tissue.imm <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds")
	tissue.vasc <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds")
	tissue.neur <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds")
	tissue.mes <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds')
	tissue.sal <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds")
	
	Idents(tissue.kera) <- "kera.res.0.9"
	Idents(tissue.imm) <-  "imm.res.0.4"
	Idents(tissue.fib) <-  "fib.res.0.6"
	Idents(tissue.vasc) <- "vasc.res.0.1"
	Idents(tissue.neur) <- "neur.res.0.1"
	Idents(tissue.mes) <-  "mes.res.0.1"
	Idents(tissue.sal) <- 'sal.res.0.4'
	
	DimPlot(tissue.imm,
	        cells = WhichCells(tissue.sal, idents = '5'),
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 0.50,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	VlnPlot(tissue.fib,
	        features = c('Col1a1', 'Dcn', 'Lum', 'Krt5', 'Krt14', 'Krt1', 'Sbsn', 'Mucl2', 'Pip', 'Dcpp2', 'Dcpp3', 'Ly6d', 'Pde4d'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	# Subcluster and annotate Keratinocytes
	###########
	
	tissue.kera <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds")
	
	Idents(tissue.kera) <- "kera.res.0.9"
	
	VlnPlot(tissue.kera,
	        features = c('Top2a', 'Cenpf', 'CC.Diff', 'S.Score', 'G2M.Score', 'Ptn', 'Bnc2', 'Sox4', 'Ucp2'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	tissue.kera <- subset(tissue.kera, idents = c('28', '30', '32'), invert = T) 
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	options(future.globals.maxSize = 1048576000)
	
	tissue.kera <- SCTransform(tissue.kera,
	                           method = 'glmGamPoi',
	                           vars.to.regress = 'CC.Diff',
	                           vst.flavor = 'v2',
	                           do.scale = T,
	                           ncells = length(colnames(tissue.kera))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.kera <- RunPCA(tissue.kera, 
	                      npcs = 100, 
	                      assay = 'SCT')
	
	ElbowPlot(tissue.kera, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.kera <- FindNeighbors(tissue.kera, 
	                             reduction = 'harmony', 
	                             dims = 1:100,
	                             assay = 'SCT')
	
	gc()
	
	tissue.kera2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.kera, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'kera.res',
	                                                                   resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.kera)
	
	tissue.kera3 <- data.frame(as.numeric(as.character(tissue.kera2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.kera3[ ,i] <- data.frame(as.numeric(as.character(tissue.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.kera3) <- col.names
	rownames(tissue.kera3) <- row.names
	
	rm(tissue.kera2)
	
	gc()
	
	tissue.kera2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.kera, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'kera.res',
	                                                                   resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.kera)
	
	tissue.kera4 <- data.frame(as.numeric(as.character(tissue.kera2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.kera4[ ,i] <- data.frame(as.numeric(as.character(tissue.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.kera4) <- col.names
	rownames(tissue.kera4) <- row.names
	
	rm(tissue.kera2)
	
	gc()
	
	tissue.kera2 <- foreach(n = seq(0.7,0.9,0.1)) %dopar% FindClusters(tissue.kera, 
	                                                                   method = 'igraph',
	                                                                   algorithm = 4,
	                                                                   cluster.name = 'kera.res',
	                                                                   resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,0.9,0.1))
	row.names <- colnames(tissue.kera)
	
	tissue.kera5 <- data.frame(as.numeric(as.character(tissue.kera2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.kera5[ ,i] <- data.frame(as.numeric(as.character(tissue.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.kera5) <- col.names
	rownames(tissue.kera5) <- row.names
	
	tissue.kera.meta <- data.frame(tissue.kera3, tissue.kera4, tissue.kera5)
	
	tissue.kera <- AddMetaData(tissue.kera, 
	                           tissue.kera.meta)
	
	rm(tissue.kera2)
	
	gc()
	
	tissue.kera <- RunUMAP(tissue.kera, 
	                       reduction = 'harmony', 
	                       reduction.name = 'umap.harmony',
	                       dims = 1:100,
	                       assay = 'SCT')
	
	saveRDS(tissue.kera,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.kera, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.kera) <- 'kera.res.0.9'
	
	DimPlot(tissue.kera,
	        split.by = 'state',
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Keratinocyte Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.kera.markers <- FindAllMarkers(tissue.kera,
	                                      only.pos = T, 
	                                      test.use = "MAST",
	                                      latent.vars = "sex",
	                                      min.pct = 0.50,
	                                      logfc.threshold = 1.00,
	                                      return.thresh = 0.05,
	                                      assay = "RNA",
	                                      densify = T) 
	
	tissue.kera.markers <- tissue.kera.markers[order(tissue.kera.markers$cluster,
	                                                 -tissue.kera.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.kera.cellcounts <- table(tissue.kera@meta.data$kera.res.0.9,
	                                tissue.kera@meta.data$state)
	
	tissue.kera.cellcounts2 <- table(tissue.kera@meta.data$kera.res.0.9,
	                                 tissue.kera@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.kera.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.kera.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.kera.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	Idents(tissue.kera) <- 'kera.res.0.9'
	
	VlnPlot(tissue.kera,
	        features = c('Krt5', 'Trp63', 'Sbsn', 'Flg',
	                     'Krt77', 'Serpinb2', 'Igfbp3', 'Ifi202b', 'Lgr6',
	                     'Krt4', 'Krt13', 'Mt3', 'Ada', 
	                     'Krt79', 'Plet1', 'Defb6', 'Cst6', 'Klk7', 'Lgr5', 
	                     'Krt15', 'Cxcl14', 'Sox9', 'Postn', 'Tenm2', 'Ptn','Elovl6', 'Mgst1'),
	        flip = T,
	        stack = T,
	        assay = 'RNA') + NoLegend()
	
	tissue.kera <- RenameIdents(tissue.kera,
	                            '1' =  'IFE',
	                            '2' =  'HF',
	                            '3' =  'HF',
	                            '4' =  'HF',
	                            '5' =  'IFE',
	                            '6' =  'HF',
	                            '7' =  'HF',
	                            '8' =  'IFE',
	                            '9' =  'HF',
	                            '10' = 'HF',
	                            '11' = 'HF',
	                            '12' = 'IFE',
	                            '13' = 'HF',
	                            '14' = 'IFE',
	                            '15' = 'HF',
	                            '16' = 'HF',
	                            '17' = 'IFE',
	                            '18' = 'IFE',
	                            '19' = 'HF',
	                            '20' = 'IFE',
	                            '21' = 'HF',
	                            '22' = 'IFE',
	                            '23' = 'IFE',
	                            '24' = 'HF',
	                            '25' = 'IFE',
	                            '26' = 'HF',
	                            '27' = 'HF',
	                            '28' = 'HF',
	                            '29' = 'HF',
	                            '30' = 'IFE')
	
	levels(tissue.kera) <- c('IFE',
	                         'HF')
	
	tissue.kera$location = tissue.kera@active.ident
	
	DimPlot(tissue.kera,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        reduction = 'umap.harmony',
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	saveRDS(tissue.kera,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds",
	        compress = FALSE)
	
	tissue.ife <- subset(tissue.kera, idents = 'IFE')
	tissue.hf <- subset(tissue.kera, idents = 'HF')
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	# IFE
	
	tissue.ife <- SCTransform(tissue.ife,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.ife))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.ife <- RunPCA(tissue.ife, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.ife, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.ife <- FindNeighbors(tissue.ife, 
	                            reduction = 'harmony', 
	                            dims = 1:100,
	                            assay = 'SCT')
	
	gc()
	
	tissue.ife@meta.data[ ,24:32] <- NULL # strip old kera.res clusters
	
	tissue.ife2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.ife, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'kera.res',
	                                                                  resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.ife)
	
	tissue.ife3 <- data.frame(as.numeric(as.character(tissue.ife2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.ife3[ ,i] <- data.frame(as.numeric(as.character(tissue.ife2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.ife3) <- col.names
	rownames(tissue.ife3) <- row.names
	
	rm(tissue.ife2)
	
	gc()
	
	tissue.ife2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.ife, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'kera.res',
	                                                                  resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.ife)
	
	tissue.ife4 <- data.frame(as.numeric(as.character(tissue.ife2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.ife4[ ,i] <- data.frame(as.numeric(as.character(tissue.ife2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.ife4) <- col.names
	rownames(tissue.ife4) <- row.names
	
	rm(tissue.ife2)
	
	gc()
	
	tissue.ife2 <- foreach(n = seq(0.7,1.0,0.1)) %dopar% FindClusters(tissue.ife, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'kera.res',
	                                                                  resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,1.0,0.1))
	row.names <- colnames(tissue.ife)
	
	tissue.ife5 <- data.frame(as.numeric(as.character(tissue.ife2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:4) {
	  tissue.ife5[ ,i] <- data.frame(as.numeric(as.character(tissue.ife2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.ife5) <- col.names
	rownames(tissue.ife5) <- row.names
	
	tissue.ife.meta <- data.frame(tissue.ife3, tissue.ife4, tissue.ife5)
	
	tissue.ife <- AddMetaData(tissue.ife, 
	                          tissue.ife.meta)
	
	rm(tissue.ife2)
	
	gc()
	
	tissue.ife <- RunUMAP(tissue.ife, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:100,
	                      assay = 'SCT')
	
	saveRDS(tissue.ife,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.ife, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte IFE Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.ife) <- 'kera.res.0.7'
	
	DimPlot(tissue.ife,
	        #split.by = 'state',
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Keratinocyte IFE UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.ife.markers <- FindAllMarkers(tissue.ife,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.ife.markers <- tissue.ife.markers[order(tissue.ife.markers$cluster,
	                                               -tissue.ife.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.ife.cellcounts <- table(tissue.ife@meta.data$kera.res.0.7,
	                               tissue.ife@meta.data$state)
	
	tissue.ife.cellcounts2 <- table(tissue.ife@meta.data$kera.res.0.7,
	                                tissue.ife@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.ife.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte IFE Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.ife.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte IFE Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.ife.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte IFE Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# Update IFE
	
	Idents(tissue.ife) <- 'kera.res.0.7'
	
	tissue.ife <- subset(tissue.ife, idents = c('10', '17'), invert = T)
	
	tissue.ife <- SCTransform(tissue.ife,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.ife))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.ife <- RunPCA(tissue.ife, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.ife, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.ife <- FindNeighbors(tissue.ife, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	gc()
	
	tissue.ife@meta.data[ ,24:32] <- NULL # strip old kera.res clusters
	
	tissue.ife2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.ife, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'kera.res',
	                                                                  resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.ife)
	
	tissue.ife3 <- data.frame(as.numeric(as.character(tissue.ife2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.ife3[ ,i] <- data.frame(as.numeric(as.character(tissue.ife2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.ife3) <- col.names
	rownames(tissue.ife3) <- row.names
	
	rm(tissue.ife2)
	
	gc()
	
	tissue.ife2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.ife, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'kera.res',
	                                                                  resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.ife)
	
	tissue.ife4 <- data.frame(as.numeric(as.character(tissue.ife2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.ife4[ ,i] <- data.frame(as.numeric(as.character(tissue.ife2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.ife4) <- col.names
	rownames(tissue.ife4) <- row.names
	
	rm(tissue.ife2)
	
	gc()
	
	tissue.ife2 <- foreach(n = seq(0.7,1.0,0.1)) %dopar% FindClusters(tissue.ife, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'kera.res',
	                                                                  resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,1.0,0.1))
	row.names <- colnames(tissue.ife)
	
	tissue.ife5 <- data.frame(as.numeric(as.character(tissue.ife2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:4) {
	  tissue.ife5[ ,i] <- data.frame(as.numeric(as.character(tissue.ife2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.ife5) <- col.names
	rownames(tissue.ife5) <- row.names
	
	tissue.ife.meta <- data.frame(tissue.ife3, tissue.ife4, tissue.ife5)
	
	tissue.ife <- AddMetaData(tissue.ife, 
	                          tissue.ife.meta)
	
	rm(tissue.ife2)
	
	gc()
	
	tissue.ife <- RunUMAP(tissue.ife, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.ife,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.ife, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte IFE Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.ife) <- 'kera.res.0.8'
	
	DimPlot(tissue.ife,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Keratinocyte IFE Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.ife.markers <- FindAllMarkers(tissue.ife,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.ife.markers <- tissue.ife.markers[order(tissue.ife.markers$cluster,
	                                               -tissue.ife.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.ife.cellcounts <- table(tissue.ife@meta.data$kera.res.0.8,
	                               tissue.ife@meta.data$state)
	
	tissue.ife.cellcounts2 <- table(tissue.ife@meta.data$kera.res.0.8,
	                                tissue.ife@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.ife.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte IFE Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.ife.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte IFE Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.ife.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte IFE Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.ife,
	        features = c('S.Score', 'G2M.Score', 'CC.Diff', 'Top2a', 'Mki67'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	# HF
	
	tissue.hf <- SCTransform(tissue.hf,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(tissue.hf))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.hf <- RunPCA(tissue.hf, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(tissue.hf, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.hf <- FindNeighbors(tissue.hf, 
	                           reduction = 'harmony', 
	                           dims = 1:100,
	                           assay = 'SCT')
	
	gc()
	
	tissue.hf@meta.data[ ,24:32] <- NULL # strip old kera.res clusters
	
	tissue.hf2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf3 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.hf3[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf3) <- col.names
	rownames(tissue.hf3) <- row.names
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf4 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.hf4[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf4) <- col.names
	rownames(tissue.hf4) <- row.names
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf2 <- foreach(n = seq(0.7,1.0,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,1.0,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf5 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:4) {
	  tissue.hf5[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf5) <- col.names
	rownames(tissue.hf5) <- row.names
	
	tissue.hf.meta <- data.frame(tissue.hf3, tissue.hf4, tissue.hf5)
	
	tissue.hf <- AddMetaData(tissue.hf, 
	                         tissue.hf.meta)
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf <- RunUMAP(tissue.hf, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     dims = 1:100,
	                     assay = 'SCT')
	
	saveRDS(tissue.hf,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.hf.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.hf, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte HF Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.hf) <- 'kera.res.0.7'
	
	DimPlot(tissue.hf,
	        #split.by = 'state',
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Keratinocyte HF UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.hf.markers <- FindAllMarkers(tissue.hf,
	                                    only.pos = T, 
	                                    test.use = "MAST",
	                                    latent.vars = "sex",
	                                    min.pct = 0.50,
	                                    logfc.threshold = 1.00,
	                                    return.thresh = 0.05,
	                                    assay = "RNA",
	                                    densify = T) 
	
	tissue.hf.markers <- tissue.hf.markers[order(tissue.hf.markers$cluster,
	                                             -tissue.hf.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.hf.cellcounts <- table(tissue.hf@meta.data$kera.res.0.7,
	                              tissue.hf@meta.data$state)
	
	tissue.hf.cellcounts2 <- table(tissue.hf@meta.data$kera.res.0.7,
	                               tissue.hf@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.hf.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.hf.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.hf.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.hf,
	        features = c('Tslp', 'Epgn', 'Ptgs2', 'Gadd45a', 'Igfbp3', 'Serpinb2', 'Cldn4'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	c18 <- FindMarkers(tissue.hf,
	                   ident.1 = '18',
	                   assay = 'RNA',
	                   densify = T,
	                   test.use = 'MAST',
	                   latent.vars = 'sex',
	                   logfc.threshold = 0.5,
	                   min.pct = 0.33)
	
	# Updated HF
	
	Idents(tissue.hf) <- 'kera.res.0.7'
	
	tissue.hf <- subset(tissue.hf, idents = c('15', '18', '19', '21', '22', '23', '24'), invert = T)
	
	tissue.hf <- SCTransform(tissue.hf,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(tissue.hf))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.hf <- RunPCA(tissue.hf, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(tissue.hf, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.hf <- FindNeighbors(tissue.hf, 
	                           reduction = 'harmony', 
	                           dims = 1:100,
	                           assay = 'SCT')
	
	gc()
	
	tissue.hf@meta.data[ ,24:32] <- NULL # strip old kera.res clusters
	
	tissue.hf2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf3 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.hf3[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf3) <- col.names
	rownames(tissue.hf3) <- row.names
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf4 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.hf4[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf4) <- col.names
	rownames(tissue.hf4) <- row.names
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf2 <- foreach(n = seq(0.7,1.0,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,1.0,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf5 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:4) {
	  tissue.hf5[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf5) <- col.names
	rownames(tissue.hf5) <- row.names
	
	tissue.hf.meta <- data.frame(tissue.hf3, tissue.hf4, tissue.hf5)
	
	tissue.hf <- AddMetaData(tissue.hf, 
	                         tissue.hf.meta)
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf <- RunUMAP(tissue.hf, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     dims = 1:100,
	                     assay = 'SCT')
	
	saveRDS(tissue.hf,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.hf.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.hf, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte HF Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.hf) <- 'kera.res.0.7'
	
	DimPlot(tissue.hf,
	        #split.by = 'state',
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Keratinocyte HF Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.hf.markers <- FindAllMarkers(tissue.hf,
	                                    only.pos = T, 
	                                    test.use = "MAST",
	                                    latent.vars = "sex",
	                                    min.pct = 0.50,
	                                    logfc.threshold = 1.00,
	                                    return.thresh = 0.05,
	                                    assay = "RNA",
	                                    densify = T) 
	
	tissue.hf.markers <- tissue.hf.markers[order(tissue.hf.markers$cluster,
	                                             -tissue.hf.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.hf.cellcounts <- table(tissue.hf@meta.data$kera.res.0.7,
	                              tissue.hf@meta.data$state)
	
	tissue.hf.cellcounts2 <- table(tissue.hf@meta.data$kera.res.0.7,
	                               tissue.hf@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.hf.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.hf.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.hf.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.hf,
	        features = c('Tslp', 'Epgn', 'Ptgs2', 'Gadd45a', 'Igfbp3'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	c1.10 <- FindMarkers(tissue.hf,
	                    ident.1 = c('1','10'),
	                    assay = 'RNA',
	                    densify = T,
	                    test.use = 'MAST',
	                    latent.vars = 'sex',
	                    only.pos = T,
	                    min.pct = 0.33,
	                    logfc.threshold = 0.5)
	
	# Updated HF v2
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	options(future.globals.maxSize = 1048576000)
	
	Idents(tissue.hf) <- 'kera.res.0.7'
	
	tissue.hf <- subset(tissue.hf, idents = c('17'), invert = T)
	
	tissue.hf <- SCTransform(tissue.hf,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(tissue.hf))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.hf <- RunPCA(tissue.hf, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(tissue.hf, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.hf <- FindNeighbors(tissue.hf, 
	                           reduction = 'harmony', 
	                           dims = 1:100,
	                           assay = 'SCT')
	
	gc()
	
	tissue.hf@meta.data[ ,24:32] <- NULL # strip old kera.res clusters
	
	tissue.hf2 <- foreach(n = seq(0.1,0.3,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.1,0.3,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf3 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.hf3[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf3) <- col.names
	rownames(tissue.hf3) <- row.names
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf2 <- foreach(n = seq(0.4,0.6,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.4,0.6,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf4 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:3) {
	  tissue.hf4[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf4) <- col.names
	rownames(tissue.hf4) <- row.names
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf2 <- foreach(n = seq(0.7,1.0,0.1)) %dopar% FindClusters(tissue.hf, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'kera.res',
	                                                                 resolution = n)
	
	col.names <- paste0('kera.res.', seq(0.7,1.0,0.1))
	row.names <- colnames(tissue.hf)
	
	tissue.hf5 <- data.frame(as.numeric(as.character(tissue.hf2[[1]]@meta.data$kera.res))) # to initialize the data frame
	
	for (i in 1:4) {
	  tissue.hf5[ ,i] <- data.frame(as.numeric(as.character(tissue.hf2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(tissue.hf5) <- col.names
	rownames(tissue.hf5) <- row.names
	
	tissue.hf.meta <- data.frame(tissue.hf3, tissue.hf4, tissue.hf5)
	
	tissue.hf <- AddMetaData(tissue.hf, 
	                         tissue.hf.meta)
	
	rm(tissue.hf2)
	
	gc()
	
	tissue.hf <- RunUMAP(tissue.hf, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     dims = 1:100,
	                     assay = 'SCT')
	
	saveRDS(tissue.hf,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.hf.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.hf, prefix = 'kera.res.')
	ggsave2(filename = "Healthy Tissue Keratinocyte HF Updated v2 Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.hf) <- 'kera.res.0.5'
	
	DimPlot(tissue.hf,
	        #split.by = 'state',
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Keratinocyte HF Updated v2 UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.hf.markers <- FindAllMarkers(tissue.hf,
	                                    only.pos = T, 
	                                    test.use = "MAST",
	                                    latent.vars = "sex",
	                                    min.pct = 0.50,
	                                    logfc.threshold = 1.00,
	                                    return.thresh = 0.05,
	                                    assay = "RNA",
	                                    densify = T) 
	
	tissue.hf.markers <- tissue.hf.markers[order(tissue.hf.markers$cluster,
	                                             -tissue.hf.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.hf.cellcounts <- table(tissue.hf@meta.data$kera.res.0.5,
	                              tissue.hf@meta.data$state)
	
	tissue.hf.cellcounts2 <- table(tissue.hf@meta.data$kera.res.0.5,
	                               tissue.hf@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.hf.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.hf.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.hf.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Keratinocyte HF Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.hf,
	        features = c('Il11ra1', 'Sfrp1', 'Scrg1', 'Gnmt', 'Edn2', 'Tagln', 'Lgr5', 'Fgf5', 'Barx2', 'Krt5', 'Krt14', 'Pthlh', 'Wfdc18', 'Krt6a', 'Cst6', 'Krt79'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	c2 <- FindMarkers(tissue.hf,
	                  ident.1 = c('2'),
	                  assay = 'RNA',
	                  densify = T,
	                  test.use = 'MAST',
	                  latent.vars = 'sex',
	                  only.pos = T,
	                  min.pct = 0.33,
	                  logfc.threshold = 0.5)
	
	# add metadata to keratinocyte subsets; subset cells and add metadata to tissue.kera
	
	tissue.ife <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds')
	tissue.hf <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.hf.rds')
	
	# IFE
	
	Idents(tissue.ife) = 'kera.res.0.8'
	
	tissue.ife <- RenameIdents(tissue.ife,
	                           '1' = 'Oral_Suprabasal_1',
	                           '2' = 'Proliferating_Keratinocyte',
	                           '3' = 'Epithelial_Suprabasal_2',
	                           '4' = 'Epithelial_Suprabasal_2',
	                           '5' = 'Oral_Basal',
	                           '6' = 'Epithelial_Suprabasal_1',
	                           '7' = 'Oral_Suprabasal_2',
	                           '8' = 'Proliferating_Keratinocyte',
	                           '9' = 'Epithelial_Basal',
	                           '10' = 'Oral_Suprabasal_1',
	                           '11' = 'Epithelial_Suprabasal_1',
	                           '12' = 'Epithelial_Suprabasal_1',
	                           '13' = 'Oral_Suprabasal_2',
	                           '14' = 'Epithelial_Basal',
	                           '15' = 'Oral_Suprabasal_1',
	                           '16' = 'Epithelial_Suprabasal_1')
	
	levels(tissue.ife) <- c('Epithelial_Basal',
	                        'Oral_Basal',
	                        'Proliferating_Keratinocyte',
	                        'Epithelial_Suprabasal_1',
	                        'Epithelial_Suprabasal_2',
	                        'Oral_Suprabasal_1',
	                        'Oral_Suprabasal_2')
	
	tissue.ife$L1subtype = tissue.ife@active.ident
	
	Idents(tissue.ife) = 'kera.res.0.8'
	
	tissue.ife <- RenameIdents(tissue.ife,
	                           '1' = 'Oral_Suprabasal_1b',
	                           '2' = 'Proliferating_Keratinocyte_1',
	                           '3' = 'Epithelial_Suprabasal_2a',
	                           '4' = 'Epithelial_Suprabasal_2b',
	                           '5' = 'Oral_Basal',
	                           '6' = 'Epithelial_Suprabasal_1b',
	                           '7' = 'Oral_Suprabasal_2a',
	                           '8' = 'Proliferating_Keratinocyte_2',
	                           '9' = 'Epithelial_Basal',
	                           '10' = 'Oral_Suprabasal_1b',
	                           '11' = 'Epithelial_Suprabasal_1b',
	                           '12' = 'Epithelial_Suprabasal_1a',
	                           '13' = 'Oral_Suprabasal_2b',
	                           '14' = 'Epithelial_Basal',
	                           '15' = 'Oral_Suprabasal_1a',
	                           '16' = 'Epithelial_Suprabasal_1a')
	
	levels(tissue.ife) <- c('Epithelial_Basal',
	                        'Oral_Basal',
	                        'Proliferating_Keratinocyte_1',
	                        'Proliferating_Keratinocyte_2',
	                        'Epithelial_Suprabasal_1a',
	                        'Epithelial_Suprabasal_1b',
	                        'Epithelial_Suprabasal_2a',
	                        'Epithelial_Suprabasal_2b',
	                        'Oral_Suprabasal_1a',
	                        'Oral_Suprabasal_1b',
	                        'Oral_Suprabasal_2a',
	                        'Oral_Suprabasal_2b')
	
	tissue.ife$L2subtype = tissue.ife@active.ident
	
	saveRDS(tissue.ife,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds",
	        compress = FALSE)
	
	# HF
	
	Idents(tissue.hf) = 'kera.res.0.5'
	
	tissue.hf <- RenameIdents(tissue.hf,
	                          '1' = 'Upper_HF_Basal',
	                          '2' = 'Lower_HF',
	                          '3' = 'Lower_HF',
	                          '4' = 'Lower_HF',
	                          '5' = 'Stress_Basal',
	                          '6' = 'Upper_HF_Suprabasal',
	                          '7' = 'Lower_HF',
	                          '8' = 'Upper_HF_Suprabasal',
	                          '9' = 'HFSC',
	                          '10' = 'Lower_HF',
	                          '11' = 'Upper_HF_Suprabasal',
	                          '12' = 'Upper_HF_Basal',
	                          '13' = 'Sebaceous',
	                          '14' = 'Upper_HF_Basal',
	                          '15' = 'Proliferating_Keratinocyte')
	
	levels(tissue.hf) <- c('Proliferating_Keratinocyte',
	                       'HFSC',
	                       'Lower_HF',
	                       'Sebaceous',
	                       'Upper_HF_Basal',
	                       'Upper_HF_Suprabasal',
	                       'Stress_Basal')
	
	tissue.hf$L1subtype = tissue.hf@active.ident
	
	Idents(tissue.hf) = 'kera.res.0.5'
	
	tissue.hf <- RenameIdents(tissue.hf,
	                          '1' = 'Upper_HF_Basal_1',
	                          '2' = 'Lower_HF_2',
	                          '3' = 'Lower_HF_3',
	                          '4' = 'Lower_HF_3',
	                          '5' = 'Stress_Basal',
	                          '6' = 'Upper_HF_Suprabasal_1',
	                          '7' = 'Lower_HF_1',
	                          '8' = 'Upper_HF_Suprabasal_3',
	                          '9' = 'HFSC',
	                          '10' = 'Lower_HF_1',
	                          '11' = 'Upper_HF_Suprabasal_2',
	                          '12' = 'Upper_HF_Basal_2',
	                          '13' = 'Sebaceous',
	                          '14' = 'Upper_HF_Basal_3',
	                          '15' = 'Proliferating_Keratinocyte_3')
	
	levels(tissue.hf) <- c('Proliferating_Keratinocyte_3',
	                       'HFSC',
	                       'Lower_HF_1',
	                       'Lower_HF_2',
	                       'Lower_HF_3',
	                       'Sebaceous',
	                       'Upper_HF_Basal_1',
	                       'Upper_HF_Basal_2',
	                       'Upper_HF_Basal_3',
	                       'Upper_HF_Suprabasal_1',
	                       'Upper_HF_Suprabasal_2',
	                       'Upper_HF_Suprabasal_3',
	                       'Stress_Basal')
	
	tissue.hf$L2subtype = tissue.hf@active.ident
	
	saveRDS(tissue.hf,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.hf.rds",
	        compress = FALSE)
	
	# Transfer metadata to tissue.kera
	
	tissue.kera <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds')
	
	kera.cells <- tissue.ife$orig.ident %>% as.data.frame %>% rownames
	hf.cells <- tissue.hf$orig.ident %>% as.data.frame %>% rownames
	
	kera.cells <- c(kera.cells, hf.cells)
	
	tissue.kera <- subset(tissue.kera, cells = kera.cells)
	
	kera.cells <- tissue.ife$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	hf.cells <- tissue.hf$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, hf.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	kera.cells <- tissue.ife$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	hf.cells <- tissue.hf$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, hf.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	tissue.kera <- AddMetaData(tissue.kera, L1subtype, col.name = 'L1subtype')
	tissue.kera <- AddMetaData(tissue.kera, L2subtype, col.name = 'L2subtype')
	
	options(future.globals.maxSize = 1048576000)
	
	tissue.kera <- SCTransform(tissue.kera,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.kera))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.kera <- RunPCA(tissue.kera, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.kera, ndims = 100) #decided to do 100 PCs for this analysis
	
	tissue.kera <- FindNeighbors(tissue.kera, 
	                            reduction = 'harmony', 
	                            dims = 1:100,
	                            assay = 'SCT')
	
	tissue.kera <- RunUMAP(tissue.kera, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:100,
	                      assay = 'SCT')
	
	saveRDS(tissue.kera,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds",
	        compress = FALSE)
	
	Idents(tissue.kera) = 'L1subtype'
	
	DimPlot(tissue.kera,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 6,
	        pt.size = 1,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Keratinocyte L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	###########
	
	# Subcluster and annotate Immune cells
	###########
	
	# Update Immune
	
	tissue.imm <- subset(tissue.imm, idents = c('11', '12', '15'), invert = T)
	
	tissue.imm <- SCTransform(tissue.imm,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.imm))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.imm <- RunPCA(tissue.imm, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.imm, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.imm <- FindNeighbors(tissue.imm, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	tissue.imm2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.imm, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'imm.res',
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.imm)
	
	tissue.imm3 <- data.frame(as.numeric(as.character(tissue.imm2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.imm3[ ,i] <- data.frame(as.numeric(as.character(tissue.imm2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.imm3) <- col.names
	rownames(tissue.imm3) <- row.names
	
	tissue.imm <- AddMetaData(tissue.imm, 
	                          tissue.imm3)
	
	rm(tissue.imm2, tissue.imm3)
	
	gc()
	
	tissue.imm <- RunUMAP(tissue.imm, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.imm,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.imm, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.imm) <- 'imm.res.0.7'
	
	DimPlot(tissue.imm,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Immune Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.imm.markers <- FindAllMarkers(tissue.imm,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.imm.markers <- tissue.imm.markers[order(tissue.imm.markers$cluster,
	                                               -tissue.imm.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.imm.cellcounts <- table(tissue.imm@meta.data$imm.res.0.7,
	                               tissue.imm@meta.data$state)
	
	tissue.imm.cellcounts2 <- table(tissue.imm@meta.data$imm.res.0.7,
	                                tissue.imm@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.imm.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.imm.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.imm.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# Update Immune #2
	
	tissue.imm <- subset(tissue.imm, idents = c('12', '13', '16', '18', '19'), invert = T)
	
	tissue.imm <- SCTransform(tissue.imm,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.imm))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.imm <- RunPCA(tissue.imm, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	
	ElbowPlot(tissue.imm, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.imm <- FindNeighbors(tissue.imm, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	tissue.imm2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.imm, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'imm.res',
	                                                                  resolution = n)
	
	
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.imm)
	
	tissue.imm3 <- data.frame(as.numeric(as.character(tissue.imm2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.imm3[ ,i] <- data.frame(as.numeric(as.character(tissue.imm2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.imm3) <- col.names
	rownames(tissue.imm3) <- row.names
	
	tissue.imm <- AddMetaData(tissue.imm, 
	                          tissue.imm3)
	
	rm(tissue.imm2, tissue.imm3)
	
	gc()
	
	tissue.imm <- RunUMAP(tissue.imm, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.imm,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.imm, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Updated v2 Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.imm) <- 'imm.res.0.8'
	
	DimPlot(tissue.imm,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Immune Updated v2 UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.imm.markers <- FindAllMarkers(tissue.imm,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.imm.markers <- tissue.imm.markers[order(tissue.imm.markers$cluster,
	                                               -tissue.imm.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.imm.cellcounts <- table(tissue.imm@meta.data$imm.res.0.8,
	                               tissue.imm@meta.data$state)
	
	tissue.imm.cellcounts2 <- table(tissue.imm@meta.data$imm.res.0.8,
	                                tissue.imm@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.imm.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.imm.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.imm.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# split Immune cells into Myeloid & Lymphoid
	
	tissue.imm <- RenameIdents(tissue.imm,
	                           '1' = 'Myeloid',
	                           '2' = 'Myeloid',
	                           '3' = 'Myeloid',
	                           '4' = 'Myeloid',
	                           '5' = 'Lymphoid',
	                           '6' = 'Myeloid',
	                           '7' = 'Myeloid',
	                           '8' = 'Myeloid',
	                           '9' = 'Lymphoid',
	                           '10' = 'Myeloid',
	                           '11' = 'Myeloid',
	                           '12' = 'Myeloid',
	                           '13' = 'Myeloid',
	                           '14' = 'Myeloid')
	
	levels(tissue.imm) <- c("Myeloid",
	                        "Lymphoid")
	
	tissue.imm$lineage <- tissue.imm@active.ident
	
	saveRDS(tissue.imm,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds",
	        compress = FALSE)
	
	tissue.lym <- subset(tissue.imm, idents = 'Lymphoid')
	
	Idents(tissue.imm) <- 'imm.res.0.8'
	
	tissue.neu <- subset(tissue.imm, idents = c('6','11','12','14')) #neutrophils
	tissue.mac <- subset(tissue.imm, idents = c('1','2','3','4','7','8','10','13')) #macrophage/dendritic
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	# Lymphoid
	
	tissue.lym <- SCTransform(tissue.lym,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.lym))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.lym <- RunPCA(tissue.lym, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.lym, ndims = 100) #decided to do 60 PCs for this analysis
	
	tissue.lym <- FindNeighbors(tissue.lym, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:60)
	
	tissue.lym@meta.data[ ,24:33] <- NULL # strip old imm.res clusters
	
	tissue.lym2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.lym, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.lym)
	
	tissue.lym3 <- data.frame(as.numeric(as.character(tissue.lym2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.lym3[ ,i] <- data.frame(as.numeric(as.character(tissue.lym2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.lym3) <- col.names
	rownames(tissue.lym3) <- row.names
	
	tissue.lym <- AddMetaData(tissue.lym, 
	                          tissue.lym3)
	
	rm(tissue.lym2, tissue.lym3)
	
	tissue.lym <- RunUMAP(tissue.lym, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:60)
	
	saveRDS(tissue.lym,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.lym.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.lym, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Lymphoid Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Neutrophil
	
	tissue.neu <- SCTransform(tissue.neu,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.neu))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.neu <- RunPCA(tissue.neu, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.neu, ndims = 100) #decided to do 70 PCs for this analysis
	
	tissue.neu <- FindNeighbors(tissue.neu, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:70)
	
	tissue.neu@meta.data[ ,24:33] <- NULL # strip old clusters
	
	tissue.neu2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.neu, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.neu)
	
	tissue.neu3 <- data.frame(as.numeric(as.character(tissue.neu2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.neu3[ ,i] <- data.frame(as.numeric(as.character(tissue.neu2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.neu3) <- col.names
	rownames(tissue.neu3) <- row.names
	
	tissue.neu <- AddMetaData(tissue.neu, 
	                          tissue.neu3)
	
	rm(tissue.neu2, tissue.neu3)
	
	tissue.neu <- RunUMAP(tissue.neu, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	saveRDS(tissue.neu,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.neu, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Neutrophil Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Macrophage/Dendritic
	
	tissue.mac <- SCTransform(tissue.mac,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.mac))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.mac <- RunPCA(tissue.mac, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.mac, ndims = 100) #decided to do 80 PCs for this analysis
	
	tissue.mac <- FindNeighbors(tissue.mac, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:80)
	
	tissue.mac@meta.data[ ,24:33] <- NULL # strip old clusters
	
	tissue.mac2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.mac, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.mac)
	
	tissue.mac3 <- data.frame(as.numeric(as.character(tissue.mac2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.mac3[ ,i] <- data.frame(as.numeric(as.character(tissue.mac2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.mac3) <- col.names
	rownames(tissue.mac3) <- row.names
	
	tissue.mac <- AddMetaData(tissue.mac, 
	                          tissue.mac3)
	
	rm(tissue.mac2, tissue.mac3)
	
	tissue.mac <- RunUMAP(tissue.mac, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:80)
	
	saveRDS(tissue.mac,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mac.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.mac, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Macrophage & Dendritic Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.lym) = 'imm.res.0.3'
	Idents(tissue.neu) = 'imm.res.0.6'
	Idents(tissue.mac) = 'imm.res.0.5'
	
	DimPlot(tissue.lym,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 3,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Immune Lymphoid UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	DimPlot(tissue.neu,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Immune Neutrophil UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	DimPlot(tissue.mac,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1.5,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Immune Macrophage & Dendritic UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.lym.markers <- FindAllMarkers(tissue.lym,
	                                      only.pos = T, 
	                                      test.use = "MAST",
	                                      latent.vars = "sex",
	                                      min.pct = 0.25,
	                                      logfc.threshold = 1.00,
	                                      return.thresh = 0.05,
	                                      assay = "RNA",
	                                      densify = T) 
	
	tissue.lym.markers <- tissue.lym.markers[order(tissue.lym.markers$cluster,
	                                                 -tissue.lym.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.lym.cellcounts <- table(tissue.lym@meta.data$imm.res.0.3,
	                                tissue.lym@meta.data$state)
	
	tissue.lym.cellcounts2 <- table(tissue.lym@meta.data$imm.res.0.3,
	                                 tissue.lym@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.lym.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Lymphoid Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.lym.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Lymphoid Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.lym.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Lymphoid Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.neu.markers <- FindAllMarkers(tissue.neu,
	                                      only.pos = T, 
	                                      test.use = "MAST",
	                                      latent.vars = "sex",
	                                      min.pct = 0.25,
	                                      logfc.threshold = 1.00,
	                                      return.thresh = 0.05,
	                                      assay = "RNA",
	                                      densify = T) 
	
	tissue.neu.markers <- tissue.neu.markers[order(tissue.neu.markers$cluster,
	                                                 -tissue.neu.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.neu.cellcounts <- table(tissue.neu@meta.data$imm.res.0.6,
	                                tissue.neu@meta.data$state)
	
	tissue.neu.cellcounts2 <- table(tissue.neu@meta.data$imm.res.0.6,
	                                 tissue.neu@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.neu.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.neu.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.neu.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	tissue.mac.markers <- FindAllMarkers(tissue.mac,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.25,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.mac.markers <- tissue.mac.markers[order(tissue.mac.markers$cluster,
	                                               -tissue.mac.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.mac.cellcounts <- table(tissue.mac@meta.data$imm.res.0.5,
	                               tissue.mac@meta.data$state)
	
	tissue.mac.cellcounts2 <- table(tissue.mac@meta.data$imm.res.0.5,
	                                tissue.mac@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.mac.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Macrophage & Dendritic Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.mac.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Macrophage & Dendritic Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.mac.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Macrophage & Dendritic Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.neu,
	        features = c('S100a8', 'S100a9', 'Slpi', 'Lyz2', 'Ly6g'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	# Neutrophil Update
	
	tissue.neu <- subset(tissue.neu, idents = c('7', '9', '10'), invert = T)
	
	tissue.neu <- SCTransform(tissue.neu,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.neu))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.neu <- RunPCA(tissue.neu, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.neu, ndims = 100) #decided to do 70 PCs for this analysis
	
	tissue.neu <- FindNeighbors(tissue.neu, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:70)
	
	tissue.neu@meta.data[ ,24:33] <- NULL # strip old clusters
	
	tissue.neu2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.neu, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.neu)
	
	tissue.neu3 <- data.frame(as.numeric(as.character(tissue.neu2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.neu3[ ,i] <- data.frame(as.numeric(as.character(tissue.neu2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.neu3) <- col.names
	rownames(tissue.neu3) <- row.names
	
	tissue.neu <- AddMetaData(tissue.neu, 
	                          tissue.neu3)
	
	rm(tissue.neu2, tissue.neu3)
	
	tissue.neu <- RunUMAP(tissue.neu, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	saveRDS(tissue.neu,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.neu, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Neutrophil Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.neu) = 'imm.res.0.8'
	
	DimPlot(tissue.neu,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Immune Neutrophil Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.neu.markers <- FindAllMarkers(tissue.neu,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.25,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.neu.markers <- tissue.neu.markers[order(tissue.neu.markers$cluster,
	                                               -tissue.neu.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.neu.cellcounts <- table(tissue.neu@meta.data$imm.res.0.8,
	                               tissue.neu@meta.data$state)
	
	tissue.neu.cellcounts2 <- table(tissue.neu@meta.data$imm.res.0.8,
	                                tissue.neu@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.neu.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Updated Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.neu.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Updated Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.neu.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Updated Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# Neutrophil Update v2
	
	Idents(tissue.neu) = 'imm.res.0.8'
	
	tissue.neu <- subset(tissue.neu, idents = c('9', '11'), invert = T)
	
	tissue.neu <- SCTransform(tissue.neu,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.neu))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.neu <- RunPCA(tissue.neu, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.neu, ndims = 100) #decided to do 70 PCs for this analysis
	
	tissue.neu <- FindNeighbors(tissue.neu, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:70)
	
	tissue.neu@meta.data[ ,24:33] <- NULL # strip old clusters
	
	tissue.neu2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.neu, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.neu)
	
	tissue.neu3 <- data.frame(as.numeric(as.character(tissue.neu2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  tissue.neu3[ ,i] <- data.frame(as.numeric(as.character(tissue.neu2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(tissue.neu3) <- col.names
	rownames(tissue.neu3) <- row.names
	
	tissue.neu <- AddMetaData(tissue.neu, 
	                          tissue.neu3)
	
	rm(tissue.neu2, tissue.neu3)
	
	tissue.neu <- RunUMAP(tissue.neu, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	saveRDS(tissue.neu,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.neu, prefix = 'imm.res.')
	ggsave2(filename = "Healthy Tissue Immune Neutrophil Updated v2 Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.neu) = 'imm.res.0.6'
	
	DimPlot(tissue.neu,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Immune Neutrophil Updated v2 UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.neu.markers <- FindAllMarkers(tissue.neu,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.25,
	                                     logfc.threshold = 0.50,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.neu.markers <- tissue.neu.markers[order(tissue.neu.markers$cluster,
	                                               -tissue.neu.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.neu.cellcounts <- table(tissue.neu@meta.data$imm.res.0.6,
	                               tissue.neu@meta.data$state)
	
	tissue.neu.cellcounts2 <- table(tissue.neu@meta.data$imm.res.0.6,
	                                tissue.neu@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.neu.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Updated v2 Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.neu.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Updated v2 Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.neu.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Immune Neutrophil Updated v2 Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.neu,
	        features = c('H2-Eb1', 'Cd74', 'Plbd1', 'Cd74'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	fib.L1subtype.markers = FindMarkers(tissue.neu,
	                                    ident.1 = c('2','8'),
	                                    only.pos = TRUE, 
	                                    test.use = "MAST",
	                                    latent.vars = "sex",
	                                    min.pct = 0.33,
	                                    logfc.threshold = 1,
	                                    assay = "RNA",
	                                    densify = TRUE) 
	
	fib.L1subtype.markers <- fib.L1subtype.markers[order(fib.L1subtype.markers$cluster,
	                                                     -fib.L1subtype.markers$avg_log2FC), ]
	
	
	# add metadata to immune cell subset; subset cells and add metadata to tissue.imm
	
	tissue.imm <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds')
	
	# Lymphoid
	
	Idents(tissue.lym) = 'imm.res.0.3'
	
	tissue.lym <- RenameIdents(tissue.lym,
	                           '1' = 'NK_cell',
	                           '2' = 'T_cell',
	                           '3' = 'T_cell',
	                           '4' = 'NK_cell',
	                           '5' = 'Proliferating_T_cell')
	
	levels(tissue.lym) <- c('T_cell', 
	                        'Proliferating_T_cell',
	                        'NK_cell')
	
	tissue.lym$L1subtype = tissue.lym@active.ident
	
	Idents(tissue.lym) = 'imm.res.0.3'
	
	tissue.lym <- RenameIdents(tissue.lym,
	                           '1' = 'NK_cell',
	                           '2' = 'T_cell',
	                           '3' = 'T_cell',
	                           '4' = 'NK_cell',
	                           '5' = 'Proliferating_T_cell')
	
	levels(tissue.lym) <- c('T_cell', 
	                        'Proliferating_T_cell',
	                        'NK_cell')
	
	tissue.lym$L2subtype = tissue.lym@active.ident
	
	saveRDS(tissue.lym,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.lym.rds",
	        compress = FALSE)
	
	# Macrophage, Dendritic, Langerhans, Monocyte
	
	Idents(tissue.mac) = 'imm.res.0.5'
	
	tissue.mac <- RenameIdents(tissue.mac,
	                           '1' = 'Macrophage_1',
	                           '2' = 'Dendritic_1',
	                           '3' = 'Macrophage_1',
	                           '4' = 'Dendritic_1',
	                           '5' = 'Macrophage_1',
	                           '6' = 'Monocyte',
	                           '7' = 'Langerhans',
	                           '8' = 'Macrophage_2')
	
	levels(tissue.mac) <- c('Macrophage_1', 
	                        'Macrophage_2',
	                        'Dendritic_1',
	                        'Langerhans',
	                        'Monocyte')
	
	tissue.mac$L1subtype = tissue.mac@active.ident
	
	Idents(tissue.mac) = 'imm.res.0.5'
	
	tissue.mac <- RenameIdents(tissue.mac,
	                           '1' = 'Macrophage_1a',
	                           '2' = 'Dendritic_1',
	                           '3' = 'Macrophage_1a',
	                           '4' = 'Dendritic_1',
	                           '5' = 'Macrophage_1b',
	                           '6' = 'Monocyte',
	                           '7' = 'Langerhans',
	                           '8' = 'Macrophage_2')
	
	levels(tissue.mac) <- c('Macrophage_1a', 
	                        'Macrophage_1b',
	                        'Macrophage_2',
	                        'Dendritic_1',
	                        'Langerhans',
	                        'Monocyte')
	
	tissue.mac$L2subtype = tissue.mac@active.ident
	
	saveRDS(tissue.mac,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mac.rds",
	        compress = FALSE)
	
	# Neutrophil (and additional macrophages)
	
	Idents(tissue.neu) = 'imm.res.0.6'
	
	tissue.neu <- RenameIdents(tissue.neu,
	                           '1' = 'Neutrophil_1',
	                           '2' = 'Dendritic_2',
	                           '3' = 'Neutrophil_1',
	                           '4' = 'Macrophage_3',
	                           '5' = 'Neutrophil_1',
	                           '6' = 'Neutrophil_1',
	                           '7' = 'Neutrophil_2',
	                           '8' = 'Dendritic_2')
	
	levels(tissue.neu) <- c('Neutrophil_1',
	                        'Neutrophil_2',
	                        'Dendritic_2',
	                        'Macrophage_3')
	
	tissue.neu$L1subtype = tissue.neu@active.ident
	
	Idents(tissue.neu) = 'imm.res.0.6'
	
	tissue.neu <- RenameIdents(tissue.neu,
	                           '1' = 'Neutrophil_1',
	                           '2' = 'Dendritic_2',
	                           '3' = 'Neutrophil_1',
	                           '4' = 'Macrophage_3',
	                           '5' = 'Neutrophil_1',
	                           '6' = 'Neutrophil_1',
	                           '7' = 'Neutrophil_2',
	                           '8' = 'Dendritic_2')
	
	levels(tissue.neu) <- c('Neutrophil_1',
	                        'Neutrophil_2',
	                        'Dendritic_2',
	                        'Macrophage_3')
	
	tissue.neu$L2subtype = tissue.neu@active.ident
	
	saveRDS(tissue.neu,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds",
	        compress = FALSE)
	
	# Transfer metadata to tissue.imm
	
	tissue.imm <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds')
	
	tissue.lym <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.lym.rds')
	tissue.mac <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mac.rds')
	tissue.neu <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds')
	
	lym.cells <- tissue.lym$orig.ident %>% as.data.frame %>% rownames
	mac.cells <- tissue.mac$orig.ident %>% as.data.frame %>% rownames
	neu.cells <- tissue.neu$orig.ident %>% as.data.frame %>% rownames
	
	imm.cells <- c(lym.cells, mac.cells, neu.cells)
	
	tissue.imm <- subset(tissue.imm, cells = imm.cells)
	
	lym.cells <- tissue.lym$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mac.cells <- tissue.mac$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	neu.cells <- tissue.neu$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(lym.cells, mac.cells, neu.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	lym.cells <- tissue.lym$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mac.cells <- tissue.mac$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	neu.cells <- tissue.neu$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(lym.cells, mac.cells, neu.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	tissue.imm <- AddMetaData(tissue.imm, L1subtype, col.name = 'L1subtype')
	tissue.imm <- AddMetaData(tissue.imm, L2subtype, col.name = 'L2subtype')
	
	tissue.imm <- SCTransform(tissue.imm,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.imm))) # workaround for 'none of the requested variables to regress...' error
	
	tissue.imm <- RunPCA(tissue.imm, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	
	ElbowPlot(tissue.imm, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.imm <- FindNeighbors(tissue.imm, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	tissue.imm <- RunUMAP(tissue.imm, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.imm,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds",
	        compress = FALSE)
	
	Idents(tissue.imm) = 'L1subtype'
	
	DimPlot(tissue.imm,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 2,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Immune L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	###########
	
	# Annotate Fibroblasts
	###########
	
	tissue.fib <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds")
	
	# Update Fibroblast
	
	Idents(tissue.fib) <-  "fib.res.0.6"
	
	tissue.fib <- subset(tissue.fib, idents = c('11', '15', '16'), invert = T)
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	tissue.fib <- SCTransform(tissue.fib,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.fib))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.fib <- RunPCA(tissue.fib, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.fib, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.fib <- FindNeighbors(tissue.fib, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	tissue.fib2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.fib, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'fib.res',
	                                                                  resolution = n)
	col.names <- paste0('fib.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.fib)
	
	tissue.fib3 <- data.frame(as.numeric(as.character(tissue.fib2[[1]]@meta.data$fib.res)))
	
	for (i in 1:10) {
	  tissue.fib3[ ,i] <- data.frame(as.numeric(as.character(tissue.fib2[[i]]@meta.data$fib.res)))
	}  
	
	colnames(tissue.fib3) <- col.names
	rownames(tissue.fib3) <- row.names
	
	tissue.fib <- AddMetaData(tissue.fib, 
	                          tissue.fib3)
	
	rm(tissue.fib2, tissue.fib3)
	
	gc()
	
	tissue.fib <- RunUMAP(tissue.fib, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.fib,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.fib, prefix = 'fib.res.')
	ggsave2(filename = "Healthy Tissue Fibroblast Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.fib) <- 'fib.res.1'
	
	DimPlot(tissue.fib,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Fibroblast Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.fib.markers <- FindAllMarkers(tissue.fib,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.fib.markers <- tissue.fib.markers[order(tissue.fib.markers$cluster,
	                                               -tissue.fib.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.fib.cellcounts <- table(tissue.fib@meta.data$fib.res.1,
	                               tissue.fib@meta.data$state)
	
	tissue.fib.cellcounts2 <- table(tissue.fib@meta.data$fib.res.1,
	                                tissue.fib@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.fib.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.fib.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.fib.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# Update Fibroblast v2
	
	Idents(tissue.fib) <- 'fib.res.1'
	
	tissue.fib <- subset(tissue.fib, idents = c('16', '21'), invert = T)
	
	tissue.fib <- SCTransform(tissue.fib,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.fib))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.fib <- RunPCA(tissue.fib, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.fib, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.fib <- FindNeighbors(tissue.fib, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	tissue.fib2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.fib, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'fib.res',
	                                                                  resolution = n)
	col.names <- paste0('fib.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.fib)
	
	tissue.fib3 <- data.frame(as.numeric(as.character(tissue.fib2[[1]]@meta.data$fib.res)))
	
	for (i in 1:10) {
	  tissue.fib3[ ,i] <- data.frame(as.numeric(as.character(tissue.fib2[[i]]@meta.data$fib.res)))
	}  
	
	colnames(tissue.fib3) <- col.names
	rownames(tissue.fib3) <- row.names
	
	tissue.fib <- AddMetaData(tissue.fib, 
	                          tissue.fib3)
	
	rm(tissue.fib2, tissue.fib3)
	
	gc()
	
	tissue.fib <- RunUMAP(tissue.fib, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.fib,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.fib, prefix = 'fib.res.')
	ggsave2(filename = "Healthy Tissue Fibroblast Updated v2 Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.fib) <- 'fib.res.0.6'
	
	DimPlot(tissue.fib,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Fibroblast Updated v2 UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.fib.markers <- FindAllMarkers(tissue.fib,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.fib.markers <- tissue.fib.markers[order(tissue.fib.markers$cluster,
	                                               -tissue.fib.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.fib.cellcounts <- table(tissue.fib@meta.data$fib.res.0.6,
	                               tissue.fib@meta.data$state)
	
	tissue.fib.cellcounts2 <- table(tissue.fib@meta.data$fib.res.0.6,
	                                tissue.fib@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.fib.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.fib.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.fib.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated v2 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# Update Fibroblast v3
	
	Idents(tissue.fib) <- 'fib.res.0.6'
	
	tissue.fib <- subset(tissue.fib, idents = c('10'), invert = T)
	
	tissue.fib <- SCTransform(tissue.fib,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.fib))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.fib <- RunPCA(tissue.fib, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.fib, ndims = 100) #decided to do 90 PCs for this analysis
	
	tissue.fib <- FindNeighbors(tissue.fib, 
	                            reduction = 'harmony', 
	                            dims = 1:90,
	                            assay = 'SCT')
	
	gc()
	
	tissue.fib2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.fib, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'fib.res',
	                                                                  resolution = n)
	col.names <- paste0('fib.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.fib)
	
	tissue.fib3 <- data.frame(as.numeric(as.character(tissue.fib2[[1]]@meta.data$fib.res)))
	
	for (i in 1:10) {
	  tissue.fib3[ ,i] <- data.frame(as.numeric(as.character(tissue.fib2[[i]]@meta.data$fib.res)))
	}  
	
	colnames(tissue.fib3) <- col.names
	rownames(tissue.fib3) <- row.names
	
	tissue.fib <- AddMetaData(tissue.fib, 
	                          tissue.fib3)
	
	rm(tissue.fib2, tissue.fib3)
	
	gc()
	
	tissue.fib <- RunUMAP(tissue.fib, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:90,
	                      assay = 'SCT')
	
	saveRDS(tissue.fib,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.fib, prefix = 'fib.res.')
	ggsave2(filename = "Healthy Tissue Fibroblast Updated v3 Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.fib) <- 'fib.res.0.5'
	
	DimPlot(tissue.fib,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1.5,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Fibroblast Updated v3 UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.fib.markers <- FindAllMarkers(tissue.fib,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.fib.markers <- tissue.fib.markers[order(tissue.fib.markers$cluster,
	                                               -tissue.fib.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.fib.cellcounts <- table(tissue.fib@meta.data$fib.res.0.5,
	                               tissue.fib@meta.data$state)
	
	tissue.fib.cellcounts2 <- table(tissue.fib@meta.data$fib.res.0.5,
	                                tissue.fib@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.fib.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated v3 Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.fib.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated v3 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.fib.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Fibroblast Updated v3 Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.fib,
	        features = c('Col1a1', 'Sparc', 'Dcn', 'Lum', 'Apcdd1', 'Apoe', 'Cxcl3',
	                     'Cxcl12', 'Gpx3', 'Mfap5', 'Plac8', 'Postn', 'Aspn', 'Enpp2'),
	        #split.by = 'state',
	        flip = T,
	        stack = T,
	        assay = 'RNA',
	        pt.size = 0)
	
	Idents(tissue.fib) <- 'fib.res.0.5'
	
	tissue.fib <- RenameIdents(tissue.fib,
	                           '1' = 'Dermal_Fibroblast_1',
	                           '2' = 'Stromal_Fibroblast_1',
	                           '3' = 'Stromal_Fibroblast_1',
	                           '4' = 'Dermal_Fibroblast_2',
	                           '5' = 'Dermal_Fibroblast_1',
	                           '6' = 'Dermal_Sheath',
	                           '7' = 'Dermal_Papilla',
	                           '8' = 'Dermal_Fibroblast_1',
	                           '9' = 'Dermal_Fibroblast_1',
	                           '10' = 'Stromal_Fibroblast_1',
	                           '11' = 'Dermal_Fibroblast_2',
	                           '12' = 'Stromal_Fibroblast_1',
	                           '13' = 'Stromal_Fibroblast_2',
	                           '14' = 'Dermal_Fibroblast_2')
	
	levels(tissue.fib) <- c('Dermal_Fibroblast_1',
	                        'Dermal_Fibroblast_2',
	                        'Stromal_Fibroblast_1',
	                        'Stromal_Fibroblast_2',
	                        'Dermal_Papilla',
	                        'Dermal_Sheath')
	
	tissue.fib$L1subtype = tissue.fib@active.ident
	
	Idents(tissue.fib) = 'fib.res.0.5'
	
	tissue.fib <- RenameIdents(tissue.fib,
	                           '1' = 'Dermal_Fibroblast_1a',
	                           '2' = 'Stromal_Fibroblast_1a',
	                           '3' = 'Stromal_Fibroblast_1a',
	                           '4' = 'Dermal_Fibroblast_2a',
	                           '5' = 'Dermal_Fibroblast_1a',
	                           '6' = 'Dermal_Sheath',
	                           '7' = 'Dermal_Papilla',
	                           '8' = 'Dermal_Fibroblast_1a',
	                           '9' = 'Dermal_Fibroblast_1b',
	                           '10' = 'Stromal_Fibroblast_1b',
	                           '11' = 'Dermal_Fibroblast_2b',
	                           '12' = 'Stromal_Fibroblast_1a',
	                           '13' = 'Stromal_Fibroblast_2',
	                           '14' = 'Dermal_Fibroblast_2a')
	
	levels(tissue.fib) <- c('Dermal_Fibroblast_1a',
	                        'Dermal_Fibroblast_1b',
	                        'Dermal_Fibroblast_2a',
	                        'Dermal_Fibroblast_2b',
	                        'Stromal_Fibroblast_1a',
	                        'Stromal_Fibroblast_1b',
	                        'Stromal_Fibroblast_2',
	                        'Dermal_Papilla',
	                        'Dermal_Sheath')
	
	tissue.fib$L2subtype = tissue.fib@active.ident
	
	saveRDS(tissue.fib,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds",
	        compress = FALSE)
	
	Idents(tissue.fib) = 'L1subtype'
	
	DimPlot(tissue.fib,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 1.5,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Fibroblast L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	
	
	
	
	
	
	
	
	
	###########
	
	# Annotate Vascular cells
	###########
	
	tissue.vasc <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds")
	
	Idents(tissue.vasc) = 'vasc.res.0.1'
	
	tissue.vasc <- RenameIdents(tissue.vasc,
	                            '1' = 'Endothelial',
	                            '2' = 'Endothelial',
	                            '3' = 'Lymph_Vessel')
	
	levels(tissue.vasc) <- c('Endothelial', 
	                         'Lymph_Vessel')
	
	tissue.vasc$L1subtype = tissue.vasc@active.ident
	
	Idents(tissue.vasc) = 'vasc.res.0.1'
	
	tissue.vasc <- RenameIdents(tissue.vasc,
	                            '1' = 'Endothelial_1',
	                            '2' = 'Endothelial_2',
	                            '3' = 'Lymph_Vessel')
	
	levels(tissue.vasc) <- c('Endothelial_1',
	                         'Endothelial_2',
	                         'Lymph_Vessel')
	
	tissue.vasc$L2subtype = tissue.vasc@active.ident
	
	# Transfer metadata to tissue.vasc
	
	vasc.cells <- tissue.vasc$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(vasc.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	vasc.cells <- tissue.vasc$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(vasc.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	tissue.vasc <- AddMetaData(tissue.vasc, L1subtype, col.name = 'L1subtype')
	tissue.vasc <- AddMetaData(tissue.vasc, L2subtype, col.name = 'L2subtype')
	
	saveRDS(tissue.vasc,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds",
	        compress = FALSE)
	
	Idents(tissue.vasc) = 'L1subtype'
	
	DimPlot(tissue.vasc,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 2,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Vascular L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	
	###########
	
	# Annotate Neural cells & melanocytes
	###########
	
	tissue.neur <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds")
	
	Idents(tissue.neur) = 'neur.res.0.1'
	
	tissue.neur <- RenameIdents(tissue.neur,
	                            '1' = 'Neural',
	                            '2' = 'Neural',
	                            '3' = 'Neural',
	                            '4' = 'Melanocyte')
	
	levels(tissue.neur) <- c('Neural', 
	                         'Melanocyte')
	
	tissue.neur$celltype = tissue.neur@active.ident
	
	Idents(tissue.neur) = 'neur.res.0.1'
	
	tissue.neur <- RenameIdents(tissue.neur,
	                           '1' = 'Myelinating_Schwann',
	                           '2' = 'Non_Myelinating_Schwann',
	                           '3' = 'Myelinating_Schwann',
	                           '4' = 'Melanocyte')
	
	levels(tissue.neur) <- c('Myelinating_Schwann', 
	                         'Non_Myelinating_Schwann',
	                         'Melanocyte')
	
	tissue.neur$L1subtype = tissue.neur@active.ident
	
	Idents(tissue.neur) = 'neur.res.0.1'
	
	tissue.neur <- RenameIdents(tissue.neur,
	                            '1' = 'Myelinating_Schwann',
	                            '2' = 'Non_Myelinating_Schwann',
	                            '3' = 'Myelinating_Schwann',
	                            '4' = 'Melanocyte')
	
	levels(tissue.neur) <- c('Myelinating_Schwann', 
	                         'Non_Myelinating_Schwann',
	                         'Melanocyte')
	
	tissue.neur$L2subtype = tissue.neur@active.ident
	
	saveRDS(tissue.neur,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds",
	        compress = FALSE)
	
	# Transfer metadata to tissue.neur
	
	neur.cells <- tissue.neur$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(neur.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	neur.cells <- tissue.neur$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(neur.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	tissue.neur <- AddMetaData(tissue.neur, L1subtype, col.name = 'L1subtype')
	tissue.neur <- AddMetaData(tissue.neur, L2subtype, col.name = 'L2subtype')
	
	saveRDS(tissue.neur,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds",
	        compress = FALSE)
	
	Idents(tissue.neur) = 'L1subtype'
	
	DimPlot(tissue.neur,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 3,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Neural L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	###########
	
	# Annotate Mesenchymal cells
	###########
	
	tissue.mes <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds")
	
	Idents(tissue.mes) = 'mes.res.0.1'
	
	tissue.mes <- RenameIdents(tissue.mes,
	                           '1' = 'Vascular',
	                           '2' = 'Skeletal_Muscle') # effectivey converts Mesenchymal celltype into Skeletal Muscle
	
	levels(tissue.mes) <- c('Vascular', 
	                        'Skeletal_Muscle')
	
	tissue.mes$celltype = tissue.mes@active.ident
	
	Idents(tissue.mes) = 'mes.res.0.1'
	
	tissue.mes <- RenameIdents(tissue.mes,
	                           '1' = 'Vascular_Smooth_Muscle',
	                           '2' = 'Skeletal_Muscle')
	
	levels(tissue.mes) <- c('Vascular_Smooth_Muscle', 
	                        'Skeletal_Muscle')
	
	tissue.mes$L1subtype = tissue.mes@active.ident
	
	Idents(tissue.mes) = 'mes.res.0.1'
	
	tissue.mes <- RenameIdents(tissue.mes,
	                           '1' = 'Vascular_Smooth_Muscle',
	                           '2' = 'Skeletal_Muscle')
	
	levels(tissue.mes) <- c('Vascular_Smooth_Muscle', 
	                        'Skeletal_Muscle')
	
	tissue.mes$L2subtype = tissue.mes@active.ident
	
	saveRDS(tissue.mes,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds",
	        compress = FALSE)
	
	Idents(tissue.mes) = 'L1subtype'
	
	DimPlot(tissue.mes,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 4,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Mesenchymal L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	
	###########
	
	# Annotate Salivary gland cells
	###########
	
	tissue.sal <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds")
	
	Idents(tissue.sal) <- 'sal.res.0.4'
	
	DimPlot(tissue.sal,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	VlnPlot(tissue.sal,
	        features = c('Ptprc', 'Cd3e', 'Cd74', 'Krt5', 'Krt14', 'Sbsn', 'Mucl2', 'Muc5b', 'Pip', 'Tff2', 'Tcea3', 'Dcpp2', 'Dcpp3'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	tissue.sal <- subset(tissue.sal, idents = c('1', '3', '13'), invert = T) #1 & 3 are keratinocytes, but since they clustered with salivary, remove
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	tissue.sal <- SCTransform(tissue.sal,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(tissue.sal))) # workaround for 'none of hte requested variables to regress...' error
	
	tissue.sal <- RunPCA(tissue.sal, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(tissue.sal, ndims = 100) #decided to do 60 PCs for this analysis
	
	tissue.sal <- FindNeighbors(tissue.sal, 
	                            reduction = 'harmony', 
	                            dims = 1:60,
	                            assay = 'SCT')
	
	tissue.sal2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(tissue.sal, 
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  cluster.name = 'sal.res',
	                                                                  resolution = n)
	col.names <- paste0('sal.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(tissue.sal)
	
	tissue.sal3 <- data.frame(as.numeric(as.character(tissue.sal2[[1]]@meta.data$sal.res)))
	
	for (i in 1:10) {
	  tissue.sal3[ ,i] <- data.frame(as.numeric(as.character(tissue.sal2[[i]]@meta.data$sal.res)))
	}  
	
	colnames(tissue.sal3) <- col.names
	rownames(tissue.sal3) <- row.names
	
	tissue.sal <- AddMetaData(tissue.sal, 
	                          tissue.sal3)
	
	rm(tissue.sal2, tissue.sal3)
	
	gc()
	
	tissue.sal <- RunUMAP(tissue.sal, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:60,
	                      assay = 'SCT')
	
	saveRDS(tissue.sal,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(tissue.sal, prefix = 'sal.res.')
	ggsave2(filename = "Healthy Tissue Salivary Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.sal) <- 'sal.res.0.5'
	
	DimPlot(tissue.sal,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 3,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	
	ggsave2(filename = "Healthy Tissue Salivary Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Skin & Oral/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.sal.markers <- FindAllMarkers(tissue.sal,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.50,
	                                     logfc.threshold = 1.00,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	tissue.sal.markers <- tissue.sal.markers[order(tissue.sal.markers$cluster,
	                                               -tissue.sal.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	tissue.sal.cellcounts <- table(tissue.sal@meta.data$sal.res.0.5,
	                               tissue.sal@meta.data$state)
	
	tissue.sal.cellcounts2 <- table(tissue.sal@meta.data$sal.res.0.5,
	                                tissue.sal@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(tissue.sal.markers,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Salivary Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(tissue.sal.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Salivary Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(tissue.sal.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Salivary Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(tissue.sal,
	        features = c('Mucl2', 'Amy1', 'Smr3a', 'Bpifa2', 'Pip'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	VlnPlot(tissue.sal,
	        features = c('Tcea3', 'Smim31', 'Smim22', 'Tesc', 'Tff2'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	c3.4 <- FindMarkers(tissue.sal,
	                    ident.1 = c('3','4'),
	                    logfc.threshold = 0.5,
	                    min.pct = 0.33,
	                    assay = 'RNA',
	                    test.use = 'MAST',
	                    latent.vars = 'sex')
	
	Idents(tissue.sal) = 'sal.res.0.5'
	
	tissue.sal <- RenameIdents(tissue.sal,
	                           '1' = 'Salivary_Basal',
	                           '2' = 'Salivary_Ductal',
	                           '3' = 'Salivary_Acinar',
	                           '4' = 'Salivary_Acinar',
	                           '5' = 'Salivary_Ductal',
	                           '6' = 'Salivary_Ductal',
	                           '7' = 'Salivary_Acinar',
	                           '8' = 'Salivary_Basal',
	                           '9' = 'Salivary_Basal',
	                           '10' = 'Salivary_Acinar',
	                           '11' = 'Salivary_Ductal',
	                           '12' = 'Salivary_Acinar')
	
	levels(tissue.sal) <- c('Salivary_Basal', 
	                        'Salivary_Acinar',
	                        'Salivary_Ductal')
	
	tissue.sal$L1subtype = tissue.sal@active.ident
	
	Idents(tissue.sal) = 'sal.res.0.5'
	
	tissue.sal <- RenameIdents(tissue.sal,
	                           '1' = 'Salivary_Basal_1',
	                           '2' = 'Salivary_Ductal_1',
	                           '3' = 'Salivary_Acinar_1',
	                           '4' = 'Salivary_Acinar_1',
	                           '5' = 'Salivary_Ductal_2',
	                           '6' = 'Salivary_Ductal_3',
	                           '7' = 'Salivary_Acinar_2',
	                           '8' = 'Salivary_Basal_2',
	                           '9' = 'Salivary_Basal_3',
	                           '10' = 'Salivary_Acinar_3',
	                           '11' = 'Salivary_Ductal_4',
	                           '12' = 'Salivary_Acinar_4')
	
	levels(tissue.sal) <- c('Salivary_Basal_1', 
	                        'Salivary_Basal_2', 
	                        'Salivary_Basal_3', 
	                        'Salivary_Acinar_1',
	                        'Salivary_Acinar_2',
	                        'Salivary_Acinar_3',
	                        'Salivary_Acinar_4',
	                        'Salivary_Ductal_1',
	                        'Salivary_Ductal_2',
	                        'Salivary_Ductal_3',
	                        'Salivary_Ductal_4')
	
	tissue.sal$L2subtype = tissue.sal@active.ident
	
	saveRDS(tissue.sal,
	        file = "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds",
	        compress = FALSE)
	
	###########
	
	# Pull metadata from subsetted objects and add info back to tissue.integrated as 'L1' & 'L2' subtype
	
	tissue.integrated <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony2.integrated.rds')
	
	tissue.kera <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds")
	tissue.fib <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds")
	tissue.imm <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds")
	tissue.vasc <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds")
	tissue.neur <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds")
	tissue.mes <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds")
	tissue.sal <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds")
	
	# celltype
	
	kera.cells <- tissue.kera$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	imm.cells <- tissue.imm$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	fib.cells <- tissue.fib$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	vasc.cells <- tissue.vasc$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	neur.cells <- tissue.neur$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	mes.cells <- tissue.mes$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	sal.cells <- tissue.sal$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, imm.cells, fib.cells, vasc.cells, neur.cells, mes.cells, sal.cells)
	
	celltype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(celltype) <- celltype[ ,1]
	celltype[ ,1] <- NULL
	colnames(celltype) <- 'celltype'
	
	# L1subtype
	
	kera.cells <- tissue.kera$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	imm.cells <- tissue.imm$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	fib.cells <- tissue.fib$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	vasc.cells <- tissue.vasc$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	neur.cells <- tissue.neur$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mes.cells <- tissue.mes$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	sal.cells <- tissue.sal$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, imm.cells, fib.cells, vasc.cells, neur.cells, mes.cells, sal.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	# L2subtype
	
	kera.cells <- tissue.kera$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	imm.cells <- tissue.imm$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	fib.cells <- tissue.fib$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	vasc.cells <- tissue.vasc$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	neur.cells <- tissue.neur$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mes.cells <- tissue.mes$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	sal.cells <- tissue.sal$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, imm.cells, fib.cells, vasc.cells, neur.cells, mes.cells, sal.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	# Add metadata to tissue.integrated and subset out unlabeled cells
	
	barcodes <- rownames(L1subtype) %>% as.character
	
	tissue.integrated <- subset(tissue.integrated,
	                            cells = barcodes)
	
	tissue.integrated <- AddMetaData(tissue.integrated, celltype, col.name = 'celltype')
	tissue.integrated <- AddMetaData(tissue.integrated, L1subtype, col.name = 'L1subtype')
	tissue.integrated <- AddMetaData(tissue.integrated, L2subtype, col.name = 'L2subtype')
	
	saveRDS(tissue.integrated, 
	        '/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony.integrated.rds',
	        compress = F)
	
	tissue.integrated <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony.integrated.rds')
	
	options(future.globals.maxSize = 1572864000)
	
	tissue.integrated <- SCTransform(tissue.integrated,
	                                 method = 'glmGamPoi',
	                                 vars.to.regress = 'CC.Diff',
	                                 vst.flavor = 'v2',
	                                 do.scale = T,
	                                 ncells = length(colnames(tissue.integrated))) # workaround for 'none of hte requested variables to regress...' error
	
	gc()
	
	tissue.integrated <- RunPCA(tissue.integrated, 
	                            npcs = 100, 
	                            assay = 'SCT')
	
	tissue.integrated <- FindNeighbors(tissue.integrated, 
	                                   reduction = 'harmony', 
	                                   assay = 'SCT',
	                                   dims = 1:100)
	
	tissue.integrated <- RunUMAP(tissue.integrated, 
	                             reduction = 'harmony', 
	                             reduction.name = 'umap.harmony',
	                             assay = 'SCT',
	                             dims = 1:100)
	
	Idents(tissue.integrated) = "celltype"
	
	DimPlot(tissue.integrated,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 0.75,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue Celltype UMAP Final.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.integrated) = "L1subtype"
	
	DimPlot(tissue.integrated,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 0.75,
	        label.size = 5,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue L1 Subtype UMAP Final.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(tissue.integrated) = "L2subtype"
	
	DimPlot(tissue.integrated,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 0.75,
	        label.size = 3,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Healthy Tissue L2 Subtype UMAP Final.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	saveRDS(tissue.integrated, 
	        '/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony.integrated.rds',
	        compress = F)
	
	# Find DEGs for each comparison in celltype & L1 subtype
	
	# Dox vs Control
	##########
	#celltype
	Idents(tissue.integrated) = "celltype"
	celltype = as.vector(unique(tissue.integrated$celltype))
	
	celltype.markers = list()
	
	for (i in 1:length(celltype)) {
	  tryCatch({
	    celltype.markers[[i]] = FindMarkers(tissue.integrated,
	                                        test.use = "MAST",
	                                        latent.vars = "sex",
	                                        logfc.threshold = 1.00,
	                                        min.pct = 0.25,
	                                        only.pos = F,
	                                        assay = "RNA",
	                                        group.by = "state",
	                                        ident.1 = "skin_dox_8_week",
	                                        ident.2 = "skin_control_8_week",
	                                        subset.ident = celltype[[i]],
	                                        densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(celltype.markers) = celltype
	
	celltype.markers.sort = list()
	
	for (j in 1:length(celltype)) {
	  
	  celltype.markers.sort[[j]] = celltype.markers[[j]]
	  celltype.markers.sort[[j]]$gene = rownames(celltype.markers[[j]])
	  celltype.markers.sort[[j]] = celltype.markers.sort[[j]][order(-celltype.markers.sort[[j]]$avg_log2FC), ]
	  celltype.markers.sort[[j]] = subset(celltype.markers.sort[[j]], p_val_adj < 0.05)
	  
	}
	
	names(celltype.markers.sort) = celltype
	
	openxlsx::write.xlsx(celltype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Dox vs Control Celltype DEGs.xlsx")
	
	# L1 subtype
	Idents(tissue.integrated) = "L1subtype"
	L1subtype = as.vector(unique(tissue.integrated$L1subtype))
	
	L1subtype.markers = list()
	
	for (i in 1:length(L1subtype)) {
	  tryCatch({
	    L1subtype.markers[[i]] = FindMarkers(tissue.integrated,
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         logfc.threshold = 1.00,
	                                         min.pct = 0.25,
	                                         only.pos = F,
	                                         assay = "RNA",
	                                         group.by = "state",
	                                         ident.1 = "skin_dox_8_week",
	                                         ident.2 = "skin_control_8_week",
	                                         subset.ident = L1subtype[[i]],
	                                         densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L1subtype.markers) = L1subtype
	
	L1subtype.markers.sort = list()
	
	for (j in 1:length(L1subtype)) {
	  tryCatch({
	    L1subtype.markers.sort[[j]] = L1subtype.markers[[j]]
	    L1subtype.markers.sort[[j]]$gene = rownames(L1subtype.markers[[j]])
	    L1subtype.markers.sort[[j]] = L1subtype.markers.sort[[j]][order(-L1subtype.markers.sort[[j]]$avg_log2FC), ]
	    L1subtype.markers.sort[[j]] = subset(L1subtype.markers.sort[[j]], p_val_adj < 0.05)
	    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L1subtype.markers.sort) = L1subtype
	
	openxlsx::write.xlsx(L1subtype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Dox vs Control L1 Subtype DEGs.xlsx")
	
	
	#L2subtype
	Idents(tissue.integrated) = "L2subtype"
	L2subtype = as.vector(unique(tissue.integrated$L2subtype))
	
	L2subtype.markers = list()
	
	for (i in 1:length(L2subtype)) {
	  tryCatch({
	    L2subtype.markers[[i]] = FindMarkers(tissue.integrated,
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         logfc.threshold = 1.00,
	                                         min.pct = 0.25,
	                                         only.pos = F,
	                                         assay = "RNA",
	                                         group.by = "state",
	                                         ident.1 = "skin_dox_8_week",
	                                         ident.2 = "skin_control_8_week",
	                                         subset.ident = L2subtype[[i]],
	                                         densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L2subtype.markers) = L2subtype[1:60]
	
	L2subtype.markers.sort = list()
	
	for (j in 1:length(L2subtype)) {
	  tryCatch({
	    L2subtype.markers.sort[[j]] = L2subtype.markers[[j]]
	    L2subtype.markers.sort[[j]]$gene = rownames(L2subtype.markers[[j]])
	    L2subtype.markers.sort[[j]] = L2subtype.markers.sort[[j]][order(-L2subtype.markers.sort[[j]]$avg_log2FC), ]
	    L2subtype.markers.sort[[j]] = subset(L2subtype.markers.sort[[j]], p_val_adj < 0.05)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L2subtype.markers.sort) = L2subtype[1:60]
	
	openxlsx::write.xlsx(L2subtype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Dox vs Control L2 Subtype DEGs.xlsx")
	##########
	
	# Buccal vs Control
	##########
	#celltype
	Idents(tissue.integrated) = "celltype"
	celltype = as.vector(unique(tissue.integrated$celltype))
	
	celltype.markers = list()
	
	for (i in 1:length(celltype)) {
	  tryCatch({
	    celltype.markers[[i]] = FindMarkers(tissue.integrated,
	                                        test.use = "MAST",
	                                        latent.vars = "sex",
	                                        logfc.threshold = 1.00,
	                                        min.pct = 0.25,
	                                        only.pos = F,
	                                        assay = "RNA",
	                                        group.by = "state",
	                                        ident.1 = "control_buccal_mucosa",
	                                        ident.2 = "skin_control_8_week",
	                                        subset.ident = celltype[[i]],
	                                        densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(celltype.markers) = celltype
	
	celltype.markers.sort = list()
	
	for (j in 1:length(celltype)) {
	  
	  celltype.markers.sort[[j]] = celltype.markers[[j]]
	  celltype.markers.sort[[j]]$gene = rownames(celltype.markers[[j]])
	  celltype.markers.sort[[j]] = celltype.markers.sort[[j]][order(-celltype.markers.sort[[j]]$avg_log2FC), ]
	  celltype.markers.sort[[j]] = subset(celltype.markers.sort[[j]], p_val_adj < 0.05)
	  
	}
	
	names(celltype.markers.sort) = celltype
	
	openxlsx::write.xlsx(celltype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Buccal vs Control Celltype DEGs.xlsx")
	
	# L1 subtype
	Idents(tissue.integrated) = "L1subtype"
	L1subtype = as.vector(unique(tissue.integrated$L1subtype))
	
	L1subtype.markers = list()
	
	for (i in 1:length(L1subtype)) {
	  tryCatch({
	    L1subtype.markers[[i]] = FindMarkers(tissue.integrated,
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         logfc.threshold = 1.00,
	                                         min.pct = 0.25,
	                                         only.pos = F,
	                                         assay = "RNA",
	                                         group.by = "state",
	                                         ident.1 = "control_buccal_mucosa",
	                                         ident.2 = "skin_control_8_week",
	                                         subset.ident = L1subtype[[i]],
	                                         densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L1subtype.markers) = L1subtype
	
	L1subtype.markers.sort = list()
	
	for (j in 1:length(L1subtype)) {
	  tryCatch({
	    L1subtype.markers.sort[[j]] = L1subtype.markers[[j]]
	    L1subtype.markers.sort[[j]]$gene = rownames(L1subtype.markers[[j]])
	    L1subtype.markers.sort[[j]] = L1subtype.markers.sort[[j]][order(-L1subtype.markers.sort[[j]]$avg_log2FC), ]
	    L1subtype.markers.sort[[j]] = subset(L1subtype.markers.sort[[j]], p_val_adj < 0.05)
	    
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L1subtype.markers.sort) = L1subtype
	
	openxlsx::write.xlsx(L1subtype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Buccal vs Control L1 Subtype DEGs.xlsx")
	
	
	#L2subtype
	Idents(tissue.integrated) = "L2subtype"
	L2subtype = as.vector(unique(tissue.integrated$L2subtype))
	
	L2subtype.markers = list()
	
	for (i in 1:length(L2subtype)) {
	  tryCatch({
	    L2subtype.markers[[i]] = FindMarkers(tissue.integrated,
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         logfc.threshold = 1.00,
	                                         min.pct = 0.25,
	                                         only.pos = F,
	                                         assay = "RNA",
	                                         group.by = "state",
	                                         ident.1 = "control_buccal_mucosa",
	                                         ident.2 = "skin_control_8_week",
	                                         subset.ident = L2subtype[[i]],
	                                         densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L2subtype.markers) = L2subtype[1:60]
	
	L2subtype.markers.sort = list()
	
	for (j in 1:length(L2subtype)) {
	  tryCatch({
	    L2subtype.markers.sort[[j]] = L2subtype.markers[[j]]
	    L2subtype.markers.sort[[j]]$gene = rownames(L2subtype.markers[[j]])
	    L2subtype.markers.sort[[j]] = L2subtype.markers.sort[[j]][order(-L2subtype.markers.sort[[j]]$avg_log2FC), ]
	    L2subtype.markers.sort[[j]] = subset(L2subtype.markers.sort[[j]], p_val_adj < 0.05)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L2subtype.markers.sort) = L2subtype[1:60]
	
	openxlsx::write.xlsx(L2subtype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Skin & Oral/Healthy Tissue Buccal vs Control L2 Subtype DEGs.xlsx")
	##########
	
	# enrichR 
	
	#######################
	#pull gene names from celltype, L1subtype, and L2subtype; organize 
	
	#celltype
	celltype.dox.genes = list()
	celltype.buccal.genes = list()
	
	for (i in 1:length(celltype)) {
	  celltype.dox.genes[[i]] = rownames(celltype.dox.markers.sort[[i]])
	  celltype.buccal.genes[[i]] = rownames(celltype.buccal.markers.sort[[i]])
	}
	
	names(celltype.dox.genes) = celltype
	names(celltype.buccal.genes)= celltype
	
	celltype.genes.list = list(celltype.dox.genes, celltype.buccal.genes)
	
	names(celltype.genes.list) = c('dox', 'buccal')
	
	#simple L2subtype
	L1subtype.dox.genes = list()
	L1subtype.buccal.genes = list()
	
	for (i in 1:length(L1subtype)) {
	  L1subtype.dox.genes[[i]] = rownames(L1subtype.dox.markers.sort[[i]])
	  L1subtype.buccal.genes[[i]] = rownames(L1subtype.buccal.markers.sort[[i]])
	}
	
	names(L1subtype.dox.genes) = L1subtype
	names(L1subtype.buccal.genes)= L1subtype
	
	L1subtype.genes.list = list(L1subtype.dox.genes, L1subtype.buccal.genes)
	
	names(L1subtype.genes.list) = c('dox', 'buccal')
	
	#L2subtype
	L2subtype.dox.genes = list()
	L2subtype.buccal.genes = list()
	
	for (i in 1:length(L2subtype)) {
	  L2subtype.dox.genes[[i]] = rownames(L2subtype.dox.markers.sort[[i]])
	  L2subtype.buccal.genes[[i]] = rownames(L2subtype.buccal.markers.sort[[i]])
	}
	
	names(L2subtype.dox.genes) = L2subtype
	names(L2subtype.buccal.genes)= L2subtype
	
	L2subtype.genes.list = list(L2subtype.dox.genes, L2subtype.buccal.genes)
	
	names(L2subtype.genes.list) = c('dox', 'buccal')
	
	#run enrichR
	
	setEnrichrSite("Enrichr") # Human/mouse genes
	
	websiteLive = T
	
	all.dbs = listEnrichrDbs()
	
	dbs = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 
	        'Reactome_2022', 'Tabula_Muris', 'Mouse_Gene_Atlas', 'Jensen_TISSUES', 'KEGG_2019_Mouse', 'MSigDB_Hallmark_2020', 
	        'WikiPathways_2019_Mouse', 'Panther_2016', 'BioPlex_2017', 'ChEA_2022', 'ENCODE_and_ChEA_Consensus_TFs_from_CHIP-X', 
	        'ENCODE_TF_ChIP-seq_2015', 'Epigenomics_Roadmap_HM_ChIP-seq', 'Transcription_Factor_PPIs', 
	        'Tissue_Protein_Expression_from_ProteomicsDB', 'TRRUST_Transcription_Factors_2019', 'KOMP2_Mouse_Phenotypes_2022', 
	        'miRTarBase_2017')
	
	#celltype
	celltype.dox.enrichr = list()
	celltype.buccal.enrichr = list()
	
	for (i in 1:length(celltype)) {
	  if (websiteLive) {
	    celltype.dox.enrichr[[i]] = enrichr(celltype.dox.genes[[i]], dbs)
	    celltype.buccal.enrichr[[i]] = enrichr(celltype.buccal.genes[[i]], dbs)
	  }
	}
	
	names(celltype.dox.enrichr) = celltype
	names(celltype.buccal.enrichr)= celltype
	
	celltype.enrichr.list = list(celltype.dox.enrichr, celltype.buccal.enrichr)
	
	names(celltype.enrichr.list) = c('Dox vs Control', 'Buccal vs Control')
	
	#L1subtype
	L1subtype.dox.enrichr = list()
	L1subtype.buccal.enrichr = list()
	
	for (i in 1:length(L1subtype)) {
	  if (websiteLive) {
	    L1subtype.dox.enrichr[[i]] = enrichr(L1subtype.dox.genes[[i]], dbs)
	    L1subtype.buccal.enrichr[[i]] = enrichr(L1subtype.buccal.genes[[i]], dbs)
	  }
	}
	
	names(L1subtype.dox.enrichr) = L1subtype
	names(L1subtype.buccal.enrichr)= L1subtype
	
	L1subtype.enrichr.list = list(L1subtype.dox.enrichr, L1subtype.buccal.enrichr)
	
	names(L1subtype.enrichr.list) = c('Dox vs Control', 'Buccal vs Control')
	
	#L2subtype
	L2subtype.dox.enrichr = list()
	L2subtype.buccal.enrichr = list()
	
	for (i in 1:length(L2subtype)) {
	  if (websiteLive) {
	    L2subtype.dox.enrichr[[i]] = enrichr(L2subtype.dox.genes[[i]], dbs)
	    L2subtype.buccal.enrichr[[i]] = enrichr(L2subtype.buccal.genes[[i]], dbs)
	  }
	}
	
	names(L2subtype.dox.enrichr) = L2subtype
	names(L2subtype.buccal.enrichr)= L2subtype
	
	L2subtype.enrichr.list = list(L2subtype.dox.enrichr, L2subtype.buccal.enrichr)
	
	names(L2subtype.enrichr.list) = c('Dox vs Control', 'Buccal vs Control')
	
	enrichr.list = list(celltype.enrichr.list, L1subtype.enrichr.list, L2subtype.enrichr.list)
	names(enrichr.list) = c('Celltype', 'L1 L2subtype', 'L2 L2subtype')
	
	saveRDS(enrichr.list,
	        '/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.enrichrlist.rds',
	        compress = T)
	###########################
	 
	#Cellchat
	###########################
	tissue.integrated <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony.integrated.rds")
	
	# subtype
	Idents(tissue.integrated) = "state"
	
	skin.control <- subset(tissue.integrated, 
	                       idents = "skin_control_8_week")
	skin.dox <- subset(tissue.integrated, 
	                   idents = "skin_dox_8_week")
	buccal <- subset(tissue.integrated,
	                 idents = "control_buccal_mucosa")
	
	# drop unused metadata
	skin.control@meta.data$celltype <- droplevels(skin.control@meta.data$celltype, exclude = setdiff(levels(skin.control@meta.data$celltype), unique(skin.control@meta.data$celltype)))
	skin.dox@meta.data$celltype <- droplevels(skin.dox@meta.data$celltype, exclude = setdiff(levels(skin.dox@meta.data$celltype),unique(skin.dox@meta.data$celltype)))
	buccal@meta.data$celltype <- droplevels(buccal@meta.data$celltype, exclude = setdiff(levels(buccal@meta.data$celltype),unique(buccal@meta.data$celltype)))
	
	skin.control@meta.data$L1subtype <- droplevels(skin.control@meta.data$L1subtype, exclude = setdiff(levels(skin.control@meta.data$L1subtype),unique(skin.control@meta.data$L1subtype)))
	skin.dox@meta.data$L1subtype <- droplevels(skin.dox@meta.data$L1subtype, exclude = setdiff(levels(skin.dox@meta.data$L1subtype),unique(skin.dox@meta.data$L1subtype)))
	buccal@meta.data$L1subtype <- droplevels(buccal@meta.data$L1subtype, exclude = setdiff(levels(buccal@meta.data$L1subtype),unique(buccal@meta.data$L1subtype)))
	
	skin.control@meta.data$L2subtype <- droplevels(skin.control@meta.data$L2subtype, exclude = setdiff(levels(skin.control@meta.data$L2subtype),unique(skin.control@meta.data$L2subtype)))
	skin.dox@meta.data$L2subtype <- droplevels(skin.dox@meta.data$L2subtype, exclude = setdiff(levels(skin.dox@meta.data$L2subtype),unique(skin.dox@meta.data$L2subtype)))
	buccal@meta.data$L2subtype <- droplevels(buccal@meta.data$L2subtype, exclude = setdiff(levels(buccal@meta.data$L2subtype),unique(buccal@meta.data$L2subtype)))
	
	# L2 subtype
	skco.cellchat <- createCellChat(object = skin.control,
	                                assay = "RNA",
	                                group.by = "L2subtype")
	skdo.cellchat <- createCellChat(object = skin.dox,
	                                assay = "RNA",
	                                group.by = "L2subtype")
	orbu.cellchat <- createCellChat(object = buccal,
	                                assay = "RNA",
	                                group.by = "L2subtype")
	
	skco.cellchat@DB <- CellChatDB.mouse
	skdo.cellchat@DB <- CellChatDB.mouse
	orbu.cellchat@DB <- CellChatDB.mouse
	
	cellchat.list <- list(skco.cellchat, skdo.cellchat, orbu.cellchat)
	
	names(cellchat.list) <- c("Control 8 Week", "+Dox 8 Week", "Buccal")
	
	cellchat.list <- lapply(cellchat.list, function(x){
	  x <- subsetData(x)
	  x <- identifyOverExpressedGenes(x)
	  x <- identifyOverExpressedInteractions(x)
	  x <- computeCommunProb(x, raw.use = T)
	  x <- filterCommunication(x, min.cells = 10)
	  x <- computeCommunProbPathway(x)
	  x <- aggregateNet(x)
	  x <- netAnalysis_computeCentrality(x)
	})
	
	saveRDS(cellchat.list, 
	        "/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.cellchat.L2subtypelist.rds",
	        compress = F)
	
	# L1 subtype
	
	skco.cellchat <- createCellChat(object = skin.control,
	                                assay = "RNA",
	                                group.by = "L1subtype")
	skdo.cellchat <- createCellChat(object = skin.dox,
	                                assay = "RNA",
	                                group.by = "L1subtype")
	orbu.cellchat <- createCellChat(object = buccal,
	                                assay = "RNA",
	                                group.by = "L1subtype")
	
	skco.cellchat@DB <- CellChatDB.mouse
	skdo.cellchat@DB <- CellChatDB.mouse
	orbu.cellchat@DB <- CellChatDB.mouse
	
	cellchat.list <- list(skco.cellchat, skdo.cellchat, orbu.cellchat)
	
	names(cellchat.list) <- c("Control 8 Week", "+Dox 8 Week", "Buccal")
	
	cellchat.list = lapply(cellchat.list, function(x){
	  x <- subsetData(x)
	  x <- identifyOverExpressedGenes(x)
	  x <- identifyOverExpressedInteractions(x)
	  x <- computeCommunProb(x, raw.use = T)
	  x <- filterCommunication(x, min.cells = 10)
	  x <- computeCommunProbPathway(x)
	  x <- aggregateNet(x)
	  x <- netAnalysis_computeCentrality(x)
	})
	
	saveRDS(cellchat.list, 
	        "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/healthytissue.cellchat.L1subtypelist.rds",
	        compress = F)
	
	# celltype
	
	skco.cellchat <- createCellChat(object = skin.control,
	                                assay = "RNA",
	                                group.by = "celltype")
	skdo.cellchat <- createCellChat(object = skin.dox,
	                                assay = "RNA",
	                                group.by = "celltype")
	orbu.cellchat <- createCellChat(object = buccal,
	                                assay = "RNA",
	                                group.by = "celltype")
	
	skco.cellchat@DB <- CellChatDB.mouse
	skdo.cellchat@DB <- CellChatDB.mouse
	orbu.cellchat@DB <- CellChatDB.mouse
	
	cellchat.list <- list(skco.cellchat, skdo.cellchat, orbu.cellchat)
	
	names(cellchat.list) <- c("Control 8 Week", "+Dox 8 Week", "Buccal")
	
	cellchat.list <- lapply(cellchat.list, function(x){
	  x <- subsetData(x)
	  x <- identifyOverExpressedGenes(x)
	  x <- identifyOverExpressedInteractions(x)
	  x <- computeCommunProb(x, raw.use = T)
	  x <- filterCommunication(x, min.cells = 10)
	  x <- computeCommunProbPathway(x)
	  x <- aggregateNet(x)
	  x <- netAnalysis_computeCentrality(x)
	})
	
	saveRDS(cellchat.list, 
	        "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/healthytissue.cellchat.celltypelist.rds",
	        compress = F)
	
	#Celltype
	#####
	
	cellchat.celltype.list <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/CellChat/healthytissue.cellchat.celltypelist.rds')
	
	pathways.con <- cellchat.celltype.list$`Control 8 Week`@netP$pathways
	pathways.dox <- cellchat.celltype.list$`+Dox 8 Week`@netP$pathways
	pathways.buc <- cellchat.celltype.list$`Buccal`@netP$pathways
	
	pathways.skin = intersect(cellchat.celltype.list$`Control 8 Week`@netP$pathways,
	                          cellchat.celltype.list$`+Dox 8 Week`@netP$pathways)
	
	pathways = intersect(pathways.skin, pathways.buc)
	
	cellchat.celltype = mergeCellChat(cellchat.celltype.list, 
	                                  add.names = names(cellchat.celltype.list))
	
	#Number and Strength of Interactions in merged objects
	gg1 = compareInteractions(cellchat.celltype, show.legend = F, group = c(1:3))
	gg2 = compareInteractions(cellchat.celltype, show.legend = F, group = c(1:3), measure = "weight")
	
	gg1|gg2
	ggsave2(filename = "Healthy Tissue Celltype Interaction Number and Strength.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/",
	        width = 4000,
	        height = 2000,
	        units = "px")
	
	#Major sources and targets in 2D space
	
	num.link = sapply(cellchat.celltype.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	
	weight.MinMax = c(40,293)
	
	gg = list()
	for (i in 1:length(cellchat.celltype.list)) {
	  gg[[i]] = netAnalysis_signalingRole_scatter(cellchat.celltype.list[[i]], 
	                                              title = names(cellchat.celltype.list)[i], 
	                                              weight.MinMax = weight.MinMax)
	}
	
	patchwork::wrap_plots(plots = gg)
	
	ggsave2(filename = "Healthy Tissue Celltype Major Sources and Targets 2D.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/",
	        width = 6000,
	        height = 2000,
	        units = "px")
	
	
	#Compare overall info flow of each signaling pathway
	
	gg1 = rankNet(cellchat.celltype, measure = "weight", mode = "comparison", stacked = F, font.size = 12, comparison = c(1:3), do.stat = TRUE)
	gg2 = rankNet(cellchat.celltype, measure = "weight", mode = "comparison", stacked = T, font.size = 12, comparison = c(1:3), do.stat = TRUE)
	
	gg1|gg2
	
	ggsave2(filename = "Healthy Tissue Celltype Overall Signaling Pathway Information Flow Combined.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/",
	        width = 4000,
	        height = 3000,
	        units = "px",
	        limitsize = F)
	
	
	#Compare signaling associated with each cell population
	
	#All tissue
	
	pathway.union = unique(c(cellchat.celltype.list[[1]]@netP$pathways, 
	                         cellchat.celltype.list[[2]]@netP$pathways,
	                         cellchat.celltype.list[[3]]@netP$pathways))
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[1]],
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.celltype.list)[1], 
	                                        width = 10, 
	                                        height = 20)
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[2]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[2], 
	                                        width = 10, 
	                                        height = 20)
	
	ht3 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[3]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[3], 
	                                        width = 10, 
	                                        height = 20)
	
	draw(ht1 + ht2 + ht3, ht_gap = unit(1, "cm"), height = unit(20, "cm"))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/Healthy Tissue Celltype Outgoing Signal Comparison.svg",
	          width = 21,
	          height = 21)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[1]],
	                                        pattern = "incoming", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.celltype.list)[1], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "GnBu")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[2]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[2], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "GnBu")
	
	ht3 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[3]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[3], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "GnBu")
	
	draw(ht1 + ht2 + ht3, ht_gap = unit(1, "cm"), height = unit(20, "cm"))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/Healthy Tissue Celltype Incoming Signal Comparison.svg",
	          width = 21,
	          height = 21)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[1]],
	                                        pattern = "all", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.celltype.list)[1], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "OrRd")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[2]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[2], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "OrRd")
	
	ht3 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[3]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[3], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "OrRd")
	
	draw(ht1 + ht2 + ht3, ht_gap = unit(1, "cm"), height = unit(20, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/Healthy Tissue Celltype Overall Signal Comparison.svg",
	          width = 20,
	          height = 20)
	
	dev.off()
	
	#DEG analysis to identify up/downregulated L-R pairs
	
	cells.plot <- as.character(levels(cellchat.celltype.list$`Control 8 Week`@meta$celltype))
	
	cells.kera <- cells.plot[1] 
	cells.immu <- cells.plot[2]
	cells.fibr <- cells.plot[3]
	cells.vasc <- cells.plot[4]
	cells.neur <- cells.plot[5]
	cells.sal  <- cells.plot[6]
	
	plot.list <- list()
	
	plot.list <- list(cells.kera, cells.fibr, cells.immu, 
	                  cells.vasc, cells.neur, cells.sal)
	
	names(plot.list) <- c('Keratinocyte', 'Fibroblast', 'Immune',
	                      'Vascular', 'Neural', 'Salivary')
	
	# +Dox as comparative dataset (pos.dataset)
	pos.dataset = "+Dox 8 Week"
	features.name = pos.dataset
	
	cellchat.celltype = identifyOverExpressedGenes(cellchat.celltype, 
	                                               group.dataset = "datasets", 
	                                               pos.dataset = pos.dataset, 
	                                               features.name = features.name, 
	                                               only.pos = FALSE, 
	                                               thresh.pc = 0.10, 
	                                               thresh.fc = 0.25, 
	                                               thresh.p = 0.05)
	
	cellchat.celltype = identifyOverExpressedInteractions(cellchat.celltype)
	
	net = netMappingDEG(cellchat.celltype, features.name = features.name)
	
	net.con = subsetCommunication(cellchat.celltype, 
	                              net = net, 
	                              datasets = "Control 8 Week",
	                              ligand.logFC = -0.3, 
	                              receptor.logFC = NULL)
	
	net.dox = subsetCommunication(cellchat.celltype,
	                              net = net, 
	                              datasets = "+Dox 8 Week",
	                              ligand.logFC = 0.2, 
	                              receptor.logFC = NULL)
	
	net.buc = subsetCommunication(cellchat.celltype,
	                              net = net, 
	                              datasets = "Buccal",
	                              ligand.logFC = -0.3, 
	                              receptor.logFC = NULL)
	
	color.use = scPalette(nlevels(cellchat.celltype.list[[1]]@idents))
	names(color.use) <- levels(cellchat.celltype.list[[1]]@idents)
	
	par(mfrow = c(1,3), xpd=TRUE)
	
	#Control 
	netVisual_chord_gene(cellchat.celltype.list[[1]], 
	                     sources.use = c(1:5), 
	                     targets.use = c(1:5), #exclude Salivary
	                     slot.name = 'net', 
	                     net = net.con, 
	                     lab.cex = 2, 
	                     small.gap = 3.5, 
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.celltype.list)[1]))
	#Dox
	netVisual_chord_gene(cellchat.celltype.list[[2]], 
	                     sources.use = c(1:5), 
	                     targets.use = c(1:5), #exclude Salivary
	                     slot.name = 'net', 
	                     net = net.dox, 
	                     lab.cex = 2, 
	                     small.gap = 3.5, 
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.celltype.list)[2]))
	#Buccal
	netVisual_chord_gene(cellchat.celltype.list[[3]], 
	                     sources.use = c(1:5), 
	                     targets.use = c(1:5), #exclude Salivary
	                     slot.name = 'net', 
	                     net = net.buc, 
	                     lab.cex = 2, 
	                     small.gap = 3.5, 
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.celltype.list)[3]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/Healthy Tissue Celltype Upregulated Signaling.svg",
	          width = 60,
	          height = 20)
	dev.off()
	
	#Compare gene expression distribution
	
	cellchat.celltype@meta$datasets = factor(cellchat.celltype@meta$datasets, 
	                                         levels = c("Control 8 Week", "+Dox 8 Week",
	                                                    'Buccal', 'Hard Palate')) # set factor level
	
	pathways.plot <- pathway.union
	
	for (i in pathways.plot) {
	  
	  plotGeneExpression(cellchat.celltype,
	                     signaling = i, 
	                     split.by = "datasets", 
	                     type = 'violin',
	                     colors.ggplot = T)
	  
	  ggsave2(filename = paste0('Healthy Tissue Celltype ', i, ' Signaling Gene Expression Violin Plots.svg'),
	          path = "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/Gene Expression Violin Plots",
	          width = 4000,
	          height = 4000,
	          units = "px",
	          limitsize = F)
	}
	
	#####
	
	#L1 Subtype
	#####
	
	cellchat.L1subtype.list <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/CellChat/healthytissue.cellchat.L1subtypelist.rds')
	
	pathways.con <- cellchat.L1subtype.list$`Control 8 Week`@netP$pathways
	pathways.dox <- cellchat.L1subtype.list$`+Dox 8 Week`@netP$pathways
	pathways.buc <- cellchat.L1subtype.list$`Buccal`@netP$pathways
	
	pathways.skin = intersect(cellchat.L1subtype.list$`Control 8 Week`@netP$pathways,
	                          cellchat.L1subtype.list$`+Dox 8 Week`@netP$pathways)
	
	pathways = intersect(pathways.skin, cellchat.L1subtype.list$`Buccal`@netP$pathways)
	
	cellchat.L1subtype = mergeCellChat(cellchat.L1subtype.list, 
	                                   add.names = names(cellchat.L1subtype.list))
	
	cells.plot <- as.character(levels(cellchat.L1subtype.list$`+Dox 8 Week`@meta$L1subtype))
	cells.kera <- cells.plot[1:11] 
	cells.immu <- cells.plot[13:19]
	cells.fibr <- cells.plot[20:23]
	cells.vasc <- cells.plot[24:26]
	cells.neur <- cells.plot[27:28]
	cells.sal  <- cells.plot[c(12,29:32)]
	
	plot.list <- list(cells.kera, cells.immu, cells.fibr, 
	                  cells.vasc, cells.neur, cells.sal)
	names(plot.list) <- c('Keratinocyte', 'Immune', 'Fibroblast',
	                      'Vascular', 'Neural', 'Salivary')
	
	color.use = scPalette(nlevels(cellchat.L1subtype.list$`+Dox 8 Week`@idents))
	names(color.use) <- levels(cellchat.L1subtype.list$`+Dox 8 Week`@idents)
	
	#Number and Strength of Interactions in merged objects
	gg1 = compareInteractions(cellchat.L1subtype, show.legend = F, group = c(1:3))
	gg2 = compareInteractions(cellchat.L1subtype, show.legend = F, group = c(1:3), measure = "weight")
	
	gg1|gg2
	ggsave2(filename = "Healthy Tissue L1 Subtype Interaction Number and Strength.svg",
	        path = "../Skin & Oral/CellChat/",
	        width = 4000,
	        height = 2000,
	        units = "px")
	
	#Major sources and targets in 2D space
	
	num.link = sapply(cellchat.L1subtype.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	
	weight.MinMax = c(0,2108)
	
	gg = list()
	for (i in 1:length(cellchat.L1subtype.list)) {
	  gg[[i]] = netAnalysis_signalingRole_scatter(cellchat.L1subtype.list[[i]], 
	                                              title = names(cellchat.L1subtype.list)[i], 
	                                              weight.MinMax = weight.MinMax,
	                                              color.use = color.use)
	}
	
	patchwork::wrap_plots(plots = gg)
	
	ggsave2(filename = "Healthy Tissue L1 Subtype Major Sources and Targets 2D.svg",
	        path = "../Skin & Oral/CellChat/",
	        width = 9000,
	        height = 3000,
	        units = "px")
	
	
	#Compare overall info flow of each signaling pathway
	
	gg1 = rankNet(cellchat.L1subtype, measure = "weight", mode = "comparison", stacked = F, font.size = 12, comparison = c(1:3), do.stat = TRUE)
	gg2 = rankNet(cellchat.L1subtype, measure = "weight", mode = "comparison", stacked = T, font.size = 12, comparison = c(1:3), do.stat = TRUE)
	
	gg1|gg2
	
	ggsave2(filename = "Healthy Tissue L1 Subtype Overall Signaling Pathway Information Flow Combined.svg",
	        path = "../Skin & Oral/CellChat/",
	        width = 6000,
	        height = 6000,
	        units = "px",
	        limitsize = F)
	
	
	#Compare signaling associated with each cell population
	
	#All tissue
	
	pathway.union = unique(c(cellchat.L1subtype.list[[1]]@netP$pathways, 
	                         cellchat.L1subtype.list[[2]]@netP$pathways,
	                         cellchat.L1subtype.list[[3]]@netP$pathways))
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[1]],
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L1subtype.list)[1], 
	                                        width = 10, 
	                                        height = 30)
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[2]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[2], 
	                                        width = 10, 
	                                        height = 30)
	
	ht3 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[3]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[3], 
	                                        width = 10, 
	                                        height = 30)
	
	draw(ht1 + ht2 + ht3, ht_gap = unit(1, "cm"), height = unit(30, "cm"))
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Outgoing Signal Comparison.svg",
	          width = 21,
	          height = 32)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[1]],
	                                        pattern = "incoming", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L1subtype.list)[1], 
	                                        width = 10, 
	                                        height = 30, 
	                                        color.heatmap = "GnBu")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[2]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[2], 
	                                        width = 10, 
	                                        height = 30, 
	                                        color.heatmap = "GnBu")
	
	ht3 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[3]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[3], 
	                                        width = 10, 
	                                        height = 30, 
	                                        color.heatmap = "GnBu")
	
	draw(ht1 + ht2 + ht3, ht_gap = unit(1, "cm"), height = unit(30, "cm"))
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Incoming Signal Comparison.svg",
	          width = 21,
	          height = 32)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[1]],
	                                        pattern = "all", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L1subtype.list)[1], 
	                                        width = 10, 
	                                        height = 30, 
	                                        color.heatmap = "OrRd")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[2]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[2], 
	                                        width = 10, 
	                                        height = 30, 
	                                        color.heatmap = "OrRd")
	
	ht3 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[3]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[3], 
	                                        width = 10, 
	                                        height = 30, 
	                                        color.heatmap = "OrRd")
	
	draw(ht1 + ht2 + ht3, ht_gap = unit(1, "cm"), height = unit(30, "cm"))
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Overall Signal Comparison.svg",
	          width = 21,
	          height = 32)
	
	dev.off()
	
	#DEG analysis to identify up/downregulated L-R pairs
	
	#Control as comparative dataset (pos.dataset)
	pos.dataset = "+Dox 8 Week"
	features.name = pos.dataset
	
	tryCatch({
	  cellchat.L1subtype = identifyOverExpressedGenes(cellchat.L1subtype, 
	                                                       group.dataset = "datasets", 
	                                                       pos.dataset = pos.dataset, 
	                                                       features.name = features.name, 
	                                                       only.pos = FALSE, 
	                                                       thresh.pc = 0.10, 
	                                                       thresh.fc = 0.25, 
	                                                       thresh.p = 0.05)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	
	cellchat.L1subtype = identifyOverExpressedInteractions(cellchat.L1subtype)
	
	net = netMappingDEG(cellchat.L1subtype, features.name = features.name)
	
	net.con = subsetCommunication(cellchat.L1subtype, 
	                              net = net, 
	                              datasets = "Control 8 Week",
	                              ligand.logFC = -0.1, 
	                              receptor.logFC = NULL)
	
	net.dox = subsetCommunication(cellchat.L1subtype,
	                              net = net, 
	                              datasets = "+Dox 8 Week",
	                              ligand.logFC = 0.1, 
	                              receptor.logFC = NULL)
	
	net.buc = subsetCommunication(cellchat.L1subtype,
	                              net = net, 
	                              datasets = "Buccal",
	                              ligand.logFC = -0.1, 
	                              receptor.logFC = NULL)
	
	gene.con = extractGeneSubsetFromPair(net.con, cellchat.L1subtype)
	gene.dox = extractGeneSubsetFromPair(net.dox, cellchat.L1subtype)
	gene.buc = extractGeneSubsetFromPair(net.buc, cellchat.L1subtype)
	
	pairLR.use.con = net.con[, "interaction_name", drop = F]
	pairLR.use.dox = net.dox[, "interaction_name", drop = F]
	pairLR.use.buc = net.buc[, "interaction_name", drop = F]
	
	lgd <- ComplexHeatmap::Legend(at = names(color.use), 
	                              type = "grid", 
	                              grid_height = unit(10, 'mm'),
	                              grid_width = unit(10, 'mm'),
	                              labels_gp = gpar(fontsize = 20),
	                              legend_gp = grid::gpar(fill = color.use), 
	                              title = "Cell State")
	
	par(mfrow = c(1,1), xpd=T)
	
	ComplexHeatmap::draw(lgd,  
	                     just = c('centre'))
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Signaling Legend.svg",
	          width = 10,
	          height = 15)
	
	#Control 
	
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                       color.use = color.use,
	                       sources.use = cells.kera, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.con, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Keratinocyte to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[1])))
	}
	  
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Keratinocyte Sender Control Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                       color.use = color.use,
	                       sources.use = cells.fibr, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.con, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Fibroblast to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[1])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Fibroblast Sender Control Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                       color.use = color.use,
	                       sources.use = cells.immu, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.con, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Immune to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[1])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Immune Sender Control Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                       color.use = color.use,
	                       sources.use = cells.vasc, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.con, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Vascular to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[1])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Vascular Sender Control Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                       color.use = color.use,
	                       sources.use = cells.neur, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.con, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Neural to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[1])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Neural Sender Control Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	
	#Dox
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[2]], 
	                       color.use = color.use,
	                       sources.use = cells.kera, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.dox, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Keratinocyte to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[2])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Keratinocyte Sender +Dox Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[2]], 
	                       color.use = color.use,
	                       sources.use = cells.fibr, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.dox, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Fibroblast to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[2])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Fibroblast Sender +Dox Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[2]], 
	                       color.use = color.use,
	                       sources.use = cells.immu, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.dox, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Immune to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[2])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Immune Sender +Dox Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[2]], 
	                       color.use = color.use,
	                       sources.use = cells.vasc, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.dox, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Vascular to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[2])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Vascular Sender +Dox Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[2]], 
	                       color.use = color.use,
	                       sources.use = cells.neur, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.dox, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Neural to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[2])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Neural Sender +Dox Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	
	#Buccal
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[3]], 
	                       color.use = color.use,
	                       sources.use = cells.kera, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.buc, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Keratinocyte to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[3])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Keratinocyte Sender Buccal Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[3]], 
	                       color.use = color.use,
	                       sources.use = cells.fibr, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.buc, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Fibroblast to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[3])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Fibroblast Sender Buccal Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[3]], 
	                       color.use = color.use,
	                       sources.use = cells.immu, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.buc, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Immune to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[3])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Immune Sender Buccal Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[3]], 
	                       color.use = color.use,
	                       sources.use = cells.vasc, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.buc, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Vascular to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[3])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Vascular Sender Buccal Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[3]], 
	                       color.use = color.use,
	                       sources.use = cells.neur, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.buc, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Neural to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[3])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Neural Sender Buccal Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	
	par(mfrow = c(2,3), xpd=T)
	
	for (i in seq_along(plot.list)) {
	  netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                       color.use = color.use,
	                       sources.use = cells.sal, 
	                       targets.use = plot.list[[i]], 
	                       slot.name = 'net', 
	                       net = net.buc, 
	                       lab.cex = 1, 
	                       small.gap = 1, 
	                       big.gap = 6,
	                       show.legend = T,
	                       title.name = paste0("Up-regulated signaling from Salivary to ", 
	                                           names(plot.list)[i],
	                                           " in ",
	                                           names(cellchat.L1subtype.list[1])))
	}
	
	dev.print(svg,
	          "../Skin & Oral/CellChat/Healthy Tissue L1 Subtype Upregulated Salivary Sender Buccal Skin Signaling.svg",
	          width = 120,
	          height = 80)
	
	dev.off()
	
	
	#Compare gene expression distribution
	
	#Compare gene expression distribution
	
	cellchat.L1subtype@meta$datasets = factor(cellchat.L1subtype@meta$datasets, 
	                                         levels = c("Control 8 Week", "+Dox 8 Week",
	                                                    'Buccal')) # set factor level
	
	pathways.plot <- pathway.union
	
	for (i in pathways.plot) {
	  
	  plotGeneExpression(cellchat.L1subtype,
	                     signaling = i, 
	                     split.by = "datasets", 
	                     type = 'violin',
	                     colors.ggplot = T)
	  
	  ggsave2(filename = paste0('Healthy Tissue L1 Subtype ', i, ' Signaling Gene Expression Violin Plots.svg'),
	          path = "/data/overmilleram/scRNAseq/Skin & Oral/CellChat/Gene Expression Violin Plots/L1Subtype",
	          width = 4000,
	          height = 4000,
	          units = "px",
	          limitsize = F)
	}
	
	#####
	###########################
	
	# MultiNicheNet
	###########################
	
	tissue.integrated <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony.integrated.rds')
	
	# Celltype
	###################
	#Load LR matrix and other files manually
	lr_network = readRDS("/data/overmilleram/scRNAseq/lr_network_mouse.rds")
	ligand_target_matrix = readRDS("/data/overmilleram/scRNAseq/ligand_target_matrix_mouse.rds")
	lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
	colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
	rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
	
	#Convert Seurat object to SCE
	
	table(tissue.integrated$celltype, tissue.integrated$state)
	
	Idents(tissue.integrated) <- 'celltype' # remove Salivary, Melanocyte, & Skeletal Muscle cells from analysis
	
	tissue.integrated <- subset(tissue.integrated,
	                            idents = c('Salivary', 'Skeletal_Muscle', 'Melanocyte'),
	                            invert = T)
	
	tissue.integrated$celltype <- droplevels(tissue.integrated$celltype) # remove unused factor levels
	tissue.integrated$L1subtype <- droplevels(tissue.integrated$L1subtype) # remove unused factor levels
	tissue.integrated$L2subtype <- droplevels(tissue.integrated$L2subtype) # remove unused factor levels
	
	tissue.sce <- as.SingleCellExperiment(tissue.integrated, assay = "RNA")
	
	tissue.sce <- alias_to_symbol_SCE(tissue.sce, "mouse")
	
	#define metadata info
	sample_id = "sampleid"
	group_id = "state"
	celltype_id = "celltype"
	covariates = "sex"
	batches = NA
	
	#define senders/receivers
	senders_oi = SummarizedExperiment::colData(tissue.sce)[,celltype_id] %>% unique()
	receivers_oi = SummarizedExperiment::colData(tissue.sce)[,celltype_id] %>% unique()
	
	#define minimum # of cells to consider for analysis
	min_cells = 50
	
	#Define ligand activity analysis parameters
	
	logFC_threshold <- 0.50
	p_val_threshold <- 0.05
	fraction_cutoff <- 0.05
	p_val_adj <- TRUE 
	empirical_pval <- F
	
	top_n_target <-250
	verbose <- TRUE
	cores_system <- as.numeric(parallel::detectCores())
	n.cores <- min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type
	
	#define prioritization weights, prepare grouping objects
	
	prioritizing_weights_DE <- c("de_ligand" = 1,
	                             "de_receptor" = 1)
	prioritizing_weights_activity <- c("activity_scaled" = 2)
	
	prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
	                                                 "exprs_receptor" = 2)
	
	prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)
	
	prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
	                                              "abund_receiver" = 0)
	
	prioritizing_weights <- c(prioritizing_weights_DE, 
	                          prioritizing_weights_activity, 
	                          prioritizing_weights_expression_specificity,
	                          prioritizing_weights_expression_sufficiency, 
	                          prioritizing_weights_relative_abundance)
	
	#Define contrasts 
	contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	
	contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	
	abundance_expression_info = get_abundance_expression_info(sce = tissue.sce, 
	                                                          sample_id = sample_id, 
	                                                          group_id = group_id, 
	                                                          celltype_id = celltype_id,
	                                                          min_cells = min_cells, 
	                                                          senders_oi = senders_oi, 
	                                                          receivers_oi = receivers_oi, 
	                                                          lr_network = lr_network, 
	                                                          batches = batches)
	
	#differential abundance per group
	abundance_expression_info$abund_plot_group + NoLegend() + theme(axis.text.x=element_text(angle=30,hjust=1))
	ggsave2(filename = "Healthy Tissue Celltype Abundances by State.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#perform DE analysis for each cell type 
	DE_info <- get_DE_info(sce = tissue.sce, 
	                       sample_id = sample_id, 
	                       group_id = group_id, 
	                       celltype_id = celltype_id, 
	                       batches = batches, 
	                       covariates = covariates, 
	                       contrasts_oi = contrasts_oi, 
	                       min_cells = min_cells)
	
	#Combine DE info for ligand-senders and receptors-receivers
	celltype_de = DE_info$celltype_de$de_output_tidy
	
	sender_receiver_de <- combine_sender_receiver_de(sender_de = celltype_de,
	                                                 receiver_de = celltype_de,
	                                                 senders_oi = senders_oi,
	                                                 receivers_oi = receivers_oi,
	                                                 lr_network = lr_network)
	
	#Run NicheNet ligand activity analysis
	
	ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
	  receiver_de = celltype_de,
	  receivers_oi = receivers_oi,
	  ligand_target_matrix = ligand_target_matrix,
	  logFC_threshold = logFC_threshold,
	  p_val_threshold = p_val_threshold,
	  p_val_adj = p_val_adj,
	  top_n_target = top_n_target,
	  verbose = verbose, 
	  n.cores = n.cores
	)
	
	sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
	
	metadata_combined <- SummarizedExperiment::colData(tissue.sce) %>% tibble::as_tibble()
	
	grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
	
	colnames(grouping_tbl) <- c("sample","group")
	
	#run prioritization
	
	prioritization_tables <- generate_prioritization_tables(
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de = sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  contrast_tbl = contrast_tbl,
	  sender_receiver_tbl = sender_receiver_tbl,
	  grouping_tbl = grouping_tbl,
	  prioritizing_weights = prioritizing_weights,
	  fraction_cutoff = fraction_cutoff, 
	  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
	  abundance_data_sender = abundance_expression_info$abundance_data_sender
	)
	
	#Add info on prior knowledge and expresion correlation between ligand-receptor and target expression
	
	lr_target_prior_cor <- lr_target_prior_cor_inference(receivers_oi, 
	                                                     abundance_expression_info, 
	                                                     celltype_de, 
	                                                     grouping_tbl, 
	                                                     prioritization_tables, 
	                                                     ligand_target_matrix, 
	                                                     logFC_threshold = logFC_threshold, 
	                                                     p_val_threshold = p_val_threshold, 
	                                                     p_val_adj = p_val_adj)
	
	#save outputs
	multinichenet_celltype_output <- list(
	  celltype_info = abundance_expression_info$celltype_info,
	  celltype_de = celltype_de,
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de =  sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  prioritization_tables = prioritization_tables,
	  grouping_tbl = grouping_tbl,
	  lr_target_prior_cor = lr_target_prior_cor
	) %>% make_lite_output()
	
	saveRDS(multinichenet_celltype_output, 
	        "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_celltype_output.rds", 
	        compress = T)
	###################
	
	# L1 Subtype
	###################
	#Load LR matrix and other files manually
	lr_network = readRDS("/data/overmilleram/scRNAseq/lr_network_mouse.rds")
	ligand_target_matrix = readRDS("/data/overmilleram/scRNAseq/ligand_target_matrix_mouse.rds")
	lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
	colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
	rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
	
	#Convert Seurat object to SCE
	
	Idents(tissue.integrated) <- 'celltype' # remove Salivary, Melanocyte, & Skeletal Muscle cells from analysis
	
	tissue.integrated <- subset(tissue.integrated,
	                            idents = c('Salivary', 'Skeletal_Muscle', 'Melanocyte'),
	                            invert = T)
	
	tissue.integrated$celltype <- droplevels(tissue.integrated$celltype) # remove unused factor levels
	tissue.integrated$L1subtype <- droplevels(tissue.integrated$L1subtype) # remove unused factor levels
	tissue.integrated$L2subtype <- droplevels(tissue.integrated$L2subtype) # remove unused factor levels
	
	tissue.sce <- as.SingleCellExperiment(tissue.integrated, assay = "RNA")
	
	tissue.sce <- alias_to_symbol_SCE(tissue.sce, "mouse")
	
	#define metadata info
	sample_id = "sampleid"
	group_id = "state"
	celltype_id = "L1subtype"
	covariates = "sex"
	batches = NA
	
	#define senders/receivers
	senders_oi = SummarizedExperiment::colData(tissue.sce)[,celltype_id] %>% unique()
	receivers_oi = SummarizedExperiment::colData(tissue.sce)[,celltype_id] %>% unique()
	
	#define minimum # of cells to consider for analysis
	min_cells = 50
	
	#Define ligand activity analysis parameters
	
	logFC_threshold <- 0.50
	p_val_threshold <- 0.05
	fraction_cutoff <- 0.05
	p_val_adj <- TRUE 
	empirical_pval <- F
	
	top_n_target <-250
	verbose <- TRUE
	cores_system <- as.numeric(parallel::detectCores())
	n.cores <- min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type
	
	#define prioritization weights, prepare grouping objects
	
	prioritizing_weights_DE <- c("de_ligand" = 1,
	                             "de_receptor" = 1)
	prioritizing_weights_activity <- c("activity_scaled" = 2)
	
	prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
	                                                 "exprs_receptor" = 2)
	
	prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)
	
	prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
	                                              "abund_receiver" = 0)
	
	prioritizing_weights <- c(prioritizing_weights_DE, 
	                          prioritizing_weights_activity, 
	                          prioritizing_weights_expression_specificity,
	                          prioritizing_weights_expression_sufficiency, 
	                          prioritizing_weights_relative_abundance)
	
	#Define contrasts 
	contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	
	contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	
	abundance_expression_info = get_abundance_expression_info(sce = tissue.sce, 
	                                                          sample_id = sample_id, 
	                                                          group_id = group_id, 
	                                                          celltype_id = celltype_id,
	                                                          min_cells = min_cells, 
	                                                          senders_oi = senders_oi, 
	                                                          receivers_oi = receivers_oi, 
	                                                          lr_network = lr_network, 
	                                                          batches = batches)
	
	#differential abundance per group
	abundance_expression_info$abund_plot_group + NoLegend() + theme(axis.text.x=element_text(angle=30,hjust=1))
	ggsave2(filename = "Healthy Tissue L1subtype Abundances by State.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#perform DE analysis for each cell type 
	DE_info <- get_DE_info(sce = tissue.sce, 
	                       sample_id = sample_id, 
	                       group_id = group_id, 
	                       celltype_id = celltype_id, 
	                       batches = batches, 
	                       covariates = covariates, 
	                       contrasts_oi = contrasts_oi, 
	                       min_cells = min_cells)
	
	#Combine DE info for ligand-senders and receptors-receivers
	celltype_de = DE_info$celltype_de$de_output_tidy
	
	sender_receiver_de <- combine_sender_receiver_de(sender_de = celltype_de,
	                                                 receiver_de = celltype_de,
	                                                 senders_oi = senders_oi,
	                                                 receivers_oi = receivers_oi,
	                                                 lr_network = lr_network)
	
	#Run NicheNet ligand activity analysis
	
	ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
	  receiver_de = celltype_de,
	  receivers_oi = receivers_oi,
	  ligand_target_matrix = ligand_target_matrix,
	  logFC_threshold = logFC_threshold,
	  p_val_threshold = p_val_threshold,
	  p_val_adj = p_val_adj,
	  top_n_target = top_n_target,
	  verbose = verbose, 
	  n.cores = n.cores
	)
	
	sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
	
	metadata_combined <- SummarizedExperiment::colData(tissue.sce) %>% tibble::as_tibble()
	
	grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
	
	colnames(grouping_tbl) <- c("sample","group")
	
	#run prioritization
	
	prioritization_tables <- generate_prioritization_tables(
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de = sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  contrast_tbl = contrast_tbl,
	  sender_receiver_tbl = sender_receiver_tbl,
	  grouping_tbl = grouping_tbl,
	  prioritizing_weights = prioritizing_weights,
	  fraction_cutoff = fraction_cutoff, 
	  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
	  abundance_data_sender = abundance_expression_info$abundance_data_sender
	)
	
	#Add info on prior knowledge and expresion correlation between ligand-receptor and target expression
	
	lr_target_prior_cor <- lr_target_prior_cor_inference(receivers_oi, 
	                                                     abundance_expression_info, 
	                                                     celltype_de, 
	                                                     grouping_tbl, 
	                                                     prioritization_tables, 
	                                                     ligand_target_matrix, 
	                                                     logFC_threshold = logFC_threshold, 
	                                                     p_val_threshold = p_val_threshold, 
	                                                     p_val_adj = p_val_adj)
	
	#save outputs
	multinichenet_L1subtype_output <- list(
	  L1subtype_info = abundance_expression_info$L1subtype_info,
	  L1subtype_de = celltype_de,
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de =  sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  prioritization_tables = prioritization_tables,
	  grouping_tbl = grouping_tbl,
	  lr_target_prior_cor = lr_target_prior_cor
	) %>% make_lite_output()
	
	saveRDS(multinichenet_L1subtype_output, 
	        "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_L1subtype_output.rds", 
	        compress = F)
	
	# Neutrophil
	###################
	#Load LR matrix and other files manually
	lr_network = readRDS("/data/overmilleram/scRNAseq/lr_network_mouse.rds")
	ligand_target_matrix = readRDS("/data/overmilleram/scRNAseq/ligand_target_matrix_mouse.rds")
	lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
	colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
	rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
	
	#Convert Seurat object to SCE
	
	Idents(tissue.integrated) <- 'celltype' # remove Salivary, Melanocyte, & Skeletal Muscle cells from analysis
	
	tissue.integrated <- subset(tissue.integrated,
	                            idents = c('Salivary', 'Skeletal_Muscle', 'Melanocyte'),
	                            invert = T)
	
	tissue.integrated$celltype <- droplevels(tissue.integrated$celltype) # remove unused factor levels
	tissue.integrated$L1subtype <- droplevels(tissue.integrated$L1subtype) # remove unused factor levels
	tissue.integrated$L2subtype <- droplevels(tissue.integrated$L2subtype) # remove unused factor levels
	
	tissue.sce <- as.SingleCellExperiment(tissue.integrated, assay = "RNA")
	
	tissue.sce <- alias_to_symbol_SCE(tissue.sce, "mouse")
	
	#define metadata info
	sample_id = "sampleid"
	group_id = "state"
	celltype_id = "L1subtype"
	covariates = "sex"
	batches = NA
	
	#define senders/receivers
	senders_oi = c('Neutrophil_1', 'Neutrophil_2')
	receivers_oi = SummarizedExperiment::colData(tissue.sce)[,celltype_id] %>% unique()
	
	#define minimum # of cells to consider for analysis
	min_cells = 10
	
	#Define ligand activity analysis parameters
	
	logFC_threshold <- 0.50
	p_val_threshold <- 0.05
	fraction_cutoff <- 0.05
	p_val_adj <- TRUE 
	empirical_pval <- F
	
	top_n_target <-250
	verbose <- TRUE
	cores_system <- as.numeric(parallel::detectCores())
	n.cores <- 8 
	
	#define prioritization weights, prepare grouping objects
	
	prioritizing_weights_DE <- c("de_ligand" = 1,
	                             "de_receptor" = 1)
	prioritizing_weights_activity <- c("activity_scaled" = 2)
	
	prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
	                                                 "exprs_receptor" = 2)
	
	prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)
	
	prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
	                                              "abund_receiver" = 0)
	
	prioritizing_weights <- c(prioritizing_weights_DE, 
	                          prioritizing_weights_activity, 
	                          prioritizing_weights_expression_specificity,
	                          prioritizing_weights_expression_sufficiency, 
	                          prioritizing_weights_relative_abundance)
	
	#Define contrasts 
	contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	
	contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	
	abundance_expression_info = get_abundance_expression_info(sce = tissue.sce, 
	                                                          sample_id = sample_id, 
	                                                          group_id = group_id, 
	                                                          celltype_id = celltype_id,
	                                                          min_cells = min_cells, 
	                                                          senders_oi = senders_oi, 
	                                                          receivers_oi = receivers_oi, 
	                                                          lr_network = lr_network, 
	                                                          batches = batches)
	
	#differential abundance per group
	abundance_expression_info$abund_plot_group + NoLegend() + theme(axis.text.x=element_text(angle=30,hjust=1))
	
	#perform DE analysis for each cell type 
	DE_info <- get_DE_info(sce = tissue.sce, 
	                       sample_id = sample_id, 
	                       group_id = group_id, 
	                       celltype_id = celltype_id, 
	                       batches = batches, 
	                       covariates = covariates, 
	                       contrasts_oi = contrasts_oi, 
	                       min_cells = min_cells)
	
	#Combine DE info for ligand-senders and receptors-receivers
	celltype_de = DE_info$celltype_de$de_output_tidy
	
	sender_receiver_de <- combine_sender_receiver_de(sender_de = celltype_de,
	                                                 receiver_de = celltype_de,
	                                                 senders_oi = senders_oi,
	                                                 receivers_oi = receivers_oi,
	                                                 lr_network = lr_network)
	
	#Run NicheNet ligand activity analysis
	
	ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
	  receiver_de = celltype_de,
	  receivers_oi = receivers_oi,
	  ligand_target_matrix = ligand_target_matrix,
	  logFC_threshold = logFC_threshold,
	  p_val_threshold = p_val_threshold,
	  p_val_adj = p_val_adj,
	  top_n_target = top_n_target,
	  verbose = verbose, 
	  n.cores = n.cores
	)
	
	sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
	
	metadata_combined <- SummarizedExperiment::colData(tissue.sce) %>% tibble::as_tibble()
	
	grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
	
	colnames(grouping_tbl) <- c("sample","group")
	
	#run prioritization
	
	prioritization_tables <- generate_prioritization_tables(
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de = sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  contrast_tbl = contrast_tbl,
	  sender_receiver_tbl = sender_receiver_tbl,
	  grouping_tbl = grouping_tbl,
	  prioritizing_weights = prioritizing_weights,
	  fraction_cutoff = fraction_cutoff, 
	  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
	  abundance_data_sender = abundance_expression_info$abundance_data_sender
	)
	
	#Add info on prior knowledge and expresion correlation between ligand-receptor and target expression
	
	lr_target_prior_cor <- lr_target_prior_cor_inference(receivers_oi, 
	                                                     abundance_expression_info, 
	                                                     celltype_de, 
	                                                     grouping_tbl, 
	                                                     prioritization_tables, 
	                                                     ligand_target_matrix, 
	                                                     logFC_threshold = logFC_threshold, 
	                                                     p_val_threshold = p_val_threshold, 
	                                                     p_val_adj = p_val_adj)
	
	#save outputs
	multinichenet_neutrophil_output <- list(
	  L1subtype_info = abundance_expression_info$L1subtype_info,
	  L1subtype_de = celltype_de,
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de =  sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  prioritization_tables = prioritization_tables,
	  grouping_tbl = grouping_tbl,
	  lr_target_prior_cor = lr_target_prior_cor
	) %>% make_lite_output()
	
	saveRDS(multinichenet_neutrophil_output, 
	        "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_neutrophil_output.rds", 
	        compress = F)
	###################
	
	#Analyze MultiNicheNet data
	
	#Celltype
	###############
	multinichenet_celltype_output <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_celltype_output.rds')
	
	#Load LR matrix and other files manually
	lr_network = readRDS("/data/overmilleram/scRNAseq/lr_network_mouse.rds")
	ligand_target_matrix = readRDS("/data/overmilleram/scRNAseq/ligand_target_matrix_mouse.rds")
	lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
	colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
	rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
	
	#Define contrasts 
	contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	
	contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	
	#visualization of top 50 interactions across all groups
	prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 50, rank_per_group = FALSE)
	
	prioritized_tbl_oi = multinichenet_celltype_output$prioritization_tables$group_prioritization_tbl %>%
	  filter(id %in% prioritized_tbl_oi_all$id) %>%
	  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
	
	prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
	
	senders_receivers <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular', 'Neural')
	
	color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E15759', '#592B02FF')
	names(color.use) <- senders_receivers
	
	circos_list = make_circos_group_comparison(prioritized_tbl_oi, color.use, color.use)
	
	#top 30 skin_control_8_week
	
	prioritized_tbl_oi_skco_30 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                30, 
	                                                groups_oi = "skin_control_8_week")
	
	circos_skco = make_circos_one_group(prioritized_tbl_oi_skco_30, color.use, color.use)
	
	#top 30 skin_dox_8_week
	
	prioritized_tbl_oi_skdo_30 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                30, 
	                                                groups_oi = "skin_dox_8_week")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_30, color.use, color.use)
	
	#top 30 buccal_control_3_week
	
	prioritized_tbl_oi_orbu_30 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                30, 
	                                                groups_oi = "control_buccal_mucosa")
	
	circos_orbu = make_circos_one_group(prioritized_tbl_oi_orbu_30, color.use, color.use)
	
	
	circos_skco[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Control Top30 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue +Dox Top30 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_orbu[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Buccal Top30 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_orbu[[2]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Top30 Circos Legend.svg",
	          width = 2,
	          height = 4)
	dev.off()
	
	
	#top 50 skin_control_8_week
	
	prioritized_tbl_oi_skco_50 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "skin_control_8_week")
	
	circos_skco = make_circos_one_group(prioritized_tbl_oi_skco_50, color.use, color.use)
	
	#top 50 skin_dox_8_week
	
	prioritized_tbl_oi_skdo_50 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "skin_dox_8_week")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_50, color.use, color.use)
	
	#top 50 buccal_control_3_week
	
	prioritized_tbl_oi_orbu_50 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "control_buccal_mucosa")
	
	circos_orbu = make_circos_one_group(prioritized_tbl_oi_orbu_50, color.use, color.use)
	
	circos_skco[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Control Top50 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue +Dox Top50 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_orbu[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Buccal Top50 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_orbu[[2]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Top50 Circos Legend.svg",
	          width = 2,
	          height = 4)
	dev.off()
	
	group_oi = list("skin_control_8_week", "skin_dox_8_week","control_buccal_mucosa")
	
	prioritized_tbl_oi_100 = lapply(X = group_oi,
	                                FUN = function(x){
	                                  get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 100, groups_oi = x)
	                                })
	
	plot_oi_100 = lapply(X = prioritized_tbl_oi_100,
	                     FUN = function(x){
	                       make_sample_lr_prod_activity_plots(multinichenet_celltype_output$prioritization_tables, x)
	                     })
	
	ggsave2(plot_oi_100[[1]],
	        filename = "Healthy Tissue Top100 Control LR Products and Ligand Activity Celltype.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 4000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[2]],
	        filename = "Healthy Tissue Top100 Dox LR Products and Ligand Activity Celltype.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 4000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[3]],
	        filename = "Healthy Tissue Top100 Buccal LR Products and Ligand Activity Celltype.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 4000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	
	#visualize expression-correlated target genes of L-R pairs 
	
	receivers_oi <- unique(senders_receivers)
	
	top_n_target <- 250
	
	lr_filtered <- multinichenet_celltype_output$lr_target_prior_cor 
	
	test <- multinichenet_celltype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)
	
	lr_filtered <- lr_filtered %>% inner_join(multinichenet_celltype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) 
	
	lr_filtered <- lr_filtered %>% inner_join(contrast_tbl) 
	
	lr_filtered <- lr_filtered %>% rowwise() 
	
	receivers_oi <- as.list(senders_receivers) 
	
	names(receivers_oi) = senders_receivers
	
	#Control Skin
	skco = filter(lr_filtered, group == "skin_control_8_week")
	
	skco.up = skco %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skco.down = skco %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skconew = bind_rows(skco.up, skco.down)
	
	prioritized_tbl_oi_skco = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "skin_control_8_week", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skco = make_lr_target_correlation_plot(multinichenet_celltype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skco,  
	                                                          skconew, 
	                                                          multinichenet_celltype_output$grouping_tbl, 
	                                                          multinichenet_celltype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skco[[1]],
	        filename = "Healthy Tissue Control Skin Celltype LR Target Correlation Plot.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 25000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skco[[2]],
	        filename = "Healthy Tissue Control Skin Celltype LR Target Correlation Plot Legend.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	#Dox Skin
	skdo = filter(lr_filtered, group == "skin_dox_8_week")
	
	skdo.up = skdo %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skdo.down = skdo %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skdonew = bind_rows(skdo.up, skdo.down)
	
	prioritized_tbl_oi_skdo = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "skin_dox_8_week", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skdo = make_lr_target_correlation_plot(multinichenet_celltype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skdo,  
	                                                          skdonew, 
	                                                          multinichenet_celltype_output$grouping_tbl, 
	                                                          multinichenet_celltype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skdo[[1]],
	        filename = "Healthy Tissue +Dox Skin Celltype LR Target Correlation Plot.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 15000,
	        height =5000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skdo[[2]],
	        filename = "Healthy Tissue +Dox Skin Celltype LR Target Correlation Plot Legend.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	#Buccal
	orbu = filter(lr_filtered, group == "control_buccal_mucosa")
	
	orbu.up = orbu %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	orbu.down = orbu %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	orbunew = bind_rows(orbu.up, orbu.down)
	
	prioritized_tbl_oi_orbu = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "control_buccal_mucosa", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_orbu = make_lr_target_correlation_plot(multinichenet_celltype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_orbu,  
	                                                          orbunew, 
	                                                          multinichenet_celltype_output$grouping_tbl, 
	                                                          multinichenet_celltype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_orbu[[1]],
	        filename = "Healthy Tissue Buccal Celltype LR Target Correlation Plot.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 35000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_orbu[[2]],
	        filename = "Healthy Tissue Buccal Celltype LR Target Correlation Plot Legend.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	###############
	
	#L1 subtype
	###############
	multinichenet_L1subtype_output <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_L1subtype_output.rds')
	
	#Define contrasts 
	contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	
	contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	
	#visualization of top 50 interactions across all groups
	prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 50, rank_per_group = FALSE)
	
	prioritized_tbl_oi = multinichenet_L1subtype_output$prioritization_tables$group_prioritization_tbl %>%
	  filter(id %in% prioritized_tbl_oi_all$id) %>%
	  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
	
	prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
	
	senders_receivers <- unique(levels(tissue.sce$L1subtype))
	
	color.use = paletteer::paletteer_d("Polychrome::palette36", 27)
	names(color.use) <- senders_receivers
	
	circos_list = make_circos_group_comparison(prioritized_tbl_oi, color.use, color.use)
	
	#top 30 skin_control_8_week
	
	prioritized_tbl_oi_skco_50 = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "skin_control_8_week")
	
	circos_skco = make_circos_one_group(prioritized_tbl_oi_skco_50, color.use, color.use)
	
	#top 30 skin_dox_8_week
	
	prioritized_tbl_oi_skdo_50 = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "skin_dox_8_week")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_30, color.use, color.use)
	
	#top 30 buccal_control_3_week
	
	prioritized_tbl_oi_orbu_50 = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "control_buccal_mucosa")
	
	circos_orbu = make_circos_one_group(prioritized_tbl_oi_orbu_30, color.use, color.use)
	
	circos_skco[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Control Top30 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue +Dox Top30 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_orbu[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Buccal Top30 Celltype Circos.svg",
	          width = 7.5,
	          height = 7.5)
	dev.off()
	
	circos_orbu[[2]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Top30 Circos Legend.svg",
	          width = 2,
	          height = 4)
	dev.off()
	
	group_oi = list("skin_control_8_week", "skin_dox_8_week","control_buccal_mucosa")
	
	prioritized_tbl_oi_100 = lapply(X = group_oi,
	                                FUN = function(x){
	                                  get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 100, groups_oi = x)
	                                })
	
	plot_oi_100 = lapply(X = prioritized_tbl_oi_100,
	                     FUN = function(x){
	                       make_sample_lr_prod_activity_plots(multinichenet_L1subtype_output$prioritization_tables, x)
	                     })
	
	ggsave2(plot_oi_100[[1]],
	        filename = "Healthy Tissue Top100 Control LR Products and Ligand Activity L1 Subtype.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 5000,
	        height = 10000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[2]],
	        filename = "Healthy Tissue Top100 Dox LR Products and Ligand Activity L1 Subtype.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 5000,
	        height = 10000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[3]],
	        filename = "Healthy Tissue Top100 Buccal LR Products and Ligand Activity L1 Subtype.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 5000,
	        height = 10000,
	        units = "px",
	        limitsize = F)
	
	#visualize expression-correlated target genes of L-R pairs 
	
	receivers_oi <- unique(senders_receivers)
	
	top_n_target <- 250
	
	lr_filtered <- multinichenet_L1subtype_output$lr_target_prior_cor 
	
	test <- multinichenet_L1subtype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)
	
	lr_filtered <- lr_filtered %>% inner_join(multinichenet_L1subtype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) 
	
	lr_filtered <- lr_filtered %>% inner_join(contrast_tbl) 
	
	lr_filtered <- lr_filtered %>% rowwise() 
	
	receivers_oi <- as.list(senders_receivers) 
	
	names(receivers_oi) = senders_receivers
	
	#Control Skin
	skco = filter(lr_filtered, group == "skin_control_8_week")
	
	skco.up = skco %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skco.down = skco %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skconew = bind_rows(skco.up, skco.down)
	
	prioritized_tbl_oi_skco = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "skin_control_8_week", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skco = make_lr_target_correlation_plot(multinichenet_L1subtype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skco,  
	                                                          skconew, 
	                                                          multinichenet_L1subtype_output$grouping_tbl, 
	                                                          multinichenet_L1subtype_output$L1subtype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skco[[1]],
	        filename = "Healthy Tissue Control Skin Celltype LR Target Correlation Plot.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 25000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skco[[2]],
	        filename = "Healthy Tissue Control Skin Celltype LR Target Correlation Plot Legend.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	#Dox Skin
	skdo = filter(lr_filtered, group == "skin_dox_8_week")
	
	skdo.up = skdo %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skdo.down = skdo %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skdonew = bind_rows(skdo.up, skdo.down)
	
	prioritized_tbl_oi_skdo = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "skin_dox_8_week", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skdo = make_lr_target_correlation_plot(multinichenet_L1subtype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skdo,  
	                                                          skdonew, 
	                                                          multinichenet_L1subtype_output$grouping_tbl, 
	                                                          multinichenet_L1subtype_output$L1subtype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skdo[[1]],
	        filename = "Healthy Tissue +Dox Skin Celltype LR Target Correlation Plot.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 15000,
	        height =5000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skdo[[2]],
	        filename = "Healthy Tissue +Dox Skin Celltype LR Target Correlation Plot Legend.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	#Buccal
	orbu = filter(lr_filtered, group == "control_buccal_mucosa")
	
	orbu.up = orbu %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	orbu.down = orbu %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	orbunew = bind_rows(orbu.up, orbu.down)
	
	prioritized_tbl_oi_orbu = get_top_n_lr_pairs(multinichenet_L1subtype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "control_buccal_mucosa", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_orbu = make_lr_target_correlation_plot(multinichenet_L1subtype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_orbu,  
	                                                          orbunew, 
	                                                          multinichenet_L1subtype_output$grouping_tbl, 
	                                                          multinichenet_L1subtype_output$L1subtype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_orbu[[1]],
	        filename = "Healthy Tissue Buccal Celltype LR Target Correlation Plot.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 35000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_orbu[[2]],
	        filename = "Healthy Tissue Buccal Celltype LR Target Correlation Plot Legend.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/",
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	multinichenet_woundskin_celltype_output <- readRDS('/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin_multinichenet_celltype_output.rds')
	multinichenet_woundskin_L1subtype_output <- readRDS('/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin_multinichenet_L1subtype_output.rds')
	
	Idents(wound.integrated) <- 'celltype' # remove Melanocyte & Mesenchymal cells from analysis
	
	wound.sce <- subset(wound.integrated,
	                    idents = c('Melanocyte', 'Skeletal_Muscle'),
	                    invert = T)
	
	wound.sce$celltype <- droplevels(wound.sce$celltype) # remove unused factor levels
	wound.sce$L1subtype <- droplevels(wound.sce$L1subtype) # remove unused factor levels
	wound.sce$L2subtype <- droplevels(wound.sce$L2subtype) # remove unused factor levels
	
	#Define contrasts 
	
	contrasts_oi <- c("'wound_control_8_week-wound_dox_8_week','wound_dox_8_week-wound_control_8_week'")
	
	contrast_tbl <- tibble(contrast = c("wound_control_8_week-wound_dox_8_week","wound_dox_8_week-wound_control_8_week"),
	                       group = c("wound_control_8_week","wound_dox_8_week"))
	
	senders_receivers <- unique(levels(wound.sce$L1subtype))
	
	color.use <- c('#0142FEFF','#1CADE4FF','#324376FF', '#219089FF','#486048FF',
	               '#B0D8D0FF','#D0F0F0BA','#6088A0FF', '#787850FF','#CC7A88FF', '#B04030FF',
	               '#FFFF80FF','#FFD700FF','#73652DFF', '#C4A000FF','#E5B17EFF', '#99540FFF',
	               '#FCD47CFF','#F57206FF','#4838A8FF', '#97AD3DFF','#15CC31FF', '#415521FF',
	               '#CCF4CFFF','#009CCEFF','#0F6A81DA', '#75FB8AFF','#7FD4C1FF', '#CD3122FF','#F6C4E5FF',
	               '#FF0010FF','#FFB79FFF')
	
	# Epithelial_Basal, Epithelial_Suprabasal, Proliferating_Keratinocyte, Upper_HF_Suprabasal, Sebaceous, 
	# HFSC_1, HFSC_2, Lower_HF_1, Lower_HF_2, Wound_Activated, Wound_Migratory
	# Dermal_Fibroblast_1, Dermal_Fibroblast_2, Dermal_Sheath, Dermal_Papilla, Myofibroblast_1, Myofibroblast_2
	#'T_cell', 'NK_Cell', 'pDC', 'Macrophage_1', 'Macrophage_2', 'Macrophage_3' 
	#'Dendritic_1', 'Dendritic_2', 'Dendritic_3', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2'
	#'Endothelial, 'Lymph_vessel'
	
	#color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E83800FF')
	
	names(color.use) <- senders_receivers
	
	# comparison multinichenet
	
	# create sample-level data frame for these interactions
	sample_data_wound = multinichenet_woundskin_L1subtype_output$prioritization_tables$sample_prioritization_tbl %>% 
	  dplyr::filter(id %in% prioritized_tbl_oi_wound$id) %>% 
	  dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
	                lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
	  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>% 
	  dplyr::arrange(sender, .by_group = TRUE) 
	
	sample_data_wound = sample_data_wound %>% dplyr::mutate(sender_receiver = factor(sender_receiver, 
	                                                                                 levels = sample_data_wound$sender_receiver %>% 
	                                                                                   unique()))
	
	# define the time point and group and link it all together
	
	grouping_tbl2_wound = multinichenet_woundskin_L1subtype_output$grouping_tbl %>% 
	  dplyr::inner_join(multinichenet_woundskin_L1subtype_output$prioritization_tables$sample_prioritization_tbl %>% 
	                      dplyr::distinct(sample, keep_receiver, keep_sender))
	
	grouping_tbl2_wound = grouping_tbl2_wound %>% inner_join(tibble(group = c('wound_control_8_week', 'wound_dox_8_week'), 
	                                                                contrast = c('wound_control_8_week', 'wound_dox_8_week')))
	
	grouping_tbl2_wound$condition = "control"
	
	grouping_tbl2_wound$condition[grouping_tbl2_wound$group %in% c('wound_dox_8_week')] = "dox"
	
	sample_data_wound = sample_data_wound %>% ungroup() %>% 
	  mutate(sampleid = sample_data_wound$sample) %>% 
	  inner_join(grouping_tbl2_wound)
	
	sample_data_wound = sample_data_wound %>% filter(keep_sender & keep_receiver) %>% 
	  mutate(group = factor(group, levels = c('wound_control_8_week', 'wound_dox_8_week')), 
	         condition = factor(condition, levels = c('control', 'dox')))
	
	# bring it together
	
	ex.1 <- sample_data_wound %>% filter(keep_receiver == 1 & keep_sender == 1) %>% ungroup() %>% 
	  dplyr::select(id, 
	                condition, 
	                sampleid,
	                ligand_receptor_pb_prod) 
	
	aggregate.ex.1 <- aggregate(ligand_receptor_pb_prod ~ id + condition, ex.1, mean)
	aggregate.ex.2 <- aggregate.ex.1[aggregate.ex.1$condition == 'control', ] %>% tidyr::spread(condition, ligand_receptor_pb_prod)
	
	join.ex.1 <- inner_join(ex.1, aggregate.ex.2, by = c('id'))
	
	names(join.ex.1) <- c('id', 'condition', 'sampleid', 'ligand_receptor_pb_prod', 'avg_control_lrpp')
	
	join.ex.2 <- join.ex.1 %>%  
	  mutate(diff = ligand_receptor_pb_prod-avg_control_lrpp,
	         fc = ligand_receptor_pb_prod/avg_control_lrpp) %>% 
	  mutate(lfc = log(fc)) %>% 
	  arrange(-lfc)
	
	sample_data_wound_2 <- sample_data_wound %>% inner_join(join.ex.2)
	
	order_sampleid = sample_data_wound_2 %>% 
	  group_by(sampleid) %>% 
	  summarise(sum_diff = sum(diff, 
	                           na.rm = TRUE)) %>% 
	  arrange(-sum_diff) %>% 
	  pull(sampleid)
	
	# bubble plot
	
	max_diff = abs(sample_data_wound_2$diff) %>% max(na.rm = TRUE)
	
	custom_scale_color = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, 
	                                                                              name = "RdBu") %>%  rev(), 
	                                           values = c(0, 0.30, 0.425, 0.5, 0.575, 0.70, 1), 
	                                           limits = c(-1 * max_diff, max_diff))
	
	p_lr_prod_change_wound = sample_data_wound_2 %>% 
	  mutate(patient = factor(sampleid, 
	                          levels = order_sampleid)) %>%
	  ggplot(aes(sampleid, 
	             lr_interaction, 
	             color = diff)) +
	  geom_point(size = 5) +
	  facet_grid(sender_receiver~contrast, 
	             scales = "free",
	             space = "free", 
	             switch = "y") +
	  theme_light() +  
	  theme(axis.ticks = element_blank(), 
	        axis.title = element_blank(), 
	        axis.text.y = element_text(face = "bold.italic", 
	                                   size = 9), 
	        axis.text.x = element_text(size = 9, 
	                                   angle = 90, 
	                                   hjust = 0), 
	        panel.grid.major = element_blank(), 
	        panel.grid.minor = element_blank(), 
	        panel.spacing.x = unit(0.4, 
	                               "lines"), 
	        panel.spacing.y = unit(0.25, 
	                               "lines"), 
	        strip.text.x.top = element_text(size = 8, 
	                                        color = "black", 
	                                        face = "bold", 
	                                        angle = 0), 
	        strip.text.y.left = element_text(size = 8, 
	                                         color = "black", 
	                                         face = "bold", 
	                                         angle = 0), 
	        strip.background = element_rect(color = "darkgrey",
	                                        fill = "whitesmoke", 
	                                        size = 1.5, 
	                                        linetype = "solid")) + 
	  custom_scale_color +
	  xlab("") +
	  ylab("")
	
	p_lr_prod_change_wound
	
	###############
	
	#Neutrophil 
	###############
	multinichenet_neutrophil_output <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_neutrophil_output.rds')
	
	#Define contrasts 
	contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	
	contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	
	#visualization of top 50 interactions across all groups
	prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_neutrophil_output$prioritization_tables, 50, rank_per_group = FALSE)
	
	prioritized_tbl_oi = multinichenet_neutrophil_output$prioritization_tables$group_prioritization_tbl %>%
	  filter(id %in% prioritized_tbl_oi_all$id) %>%
	  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
	
	prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
	
	senders_receivers <- unique(levels(tissue.sce$L1subtype))
	
	color.use = paletteer::paletteer_d("Polychrome::palette36", 35)
	names(color.use) <- senders_receivers
	
	circos_list = make_circos_group_comparison(prioritized_tbl_oi, color.use, color.use)
	
	circos_list[[2]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Control Top50 Share LR Links Neutrophil.svg",
	          width = 5.5,
	          height = 5.5)
	dev.off()
	
	circos_list[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue +Dox Top50 Share LR Links Neutrophil.svg",
	          width = 5.5,
	          height = 5.5)
	dev.off()
	
	circos_list[[4]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/Healthy Tissue Legend Top50 Share LR Links Neutrophil.svg",
	          width = 5,
	          height = 10)
	dev.off()
	
	###########################
	
	#Pseudotime trajectory analysis on keratinocytes
	###########################
	
	tissue.ife <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds')
	
	Idents(tissue.ife) <- 'L2subtype'
	
	DimPlot(tissue.ife,
	        pt.size = 1.5,
	        label = T,
	        repel = T,
	        label.size = 7,
	        raster = F) + NoLegend()
	ggsave2(filename = "Healthy Keratinocyte IFE L2 Subtype UMAP Final.svg",
	        path = "/data/overmilleram/scRNAseq/Skin & Oral/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	tissue.ife[["RNA3"]] <- as(object = tissue.ife[["RNA"]], Class = "Assay") # convert to Seurat v3 object to play nice downstream # make sure counts matrix is there
	
	set.seed(42) #same one Seurat and others use
	
	# keep UMAP embeddings calculated from Seurat for pseudotime analysis
	
	ife.cds <- as.cell_data_set(tissue.ife, assay = 'RNA3')
	
	reducedDim(ife.cds, type = "PCA") <- tissue.ife@reductions$pca@cell.embeddings
	
	ife.cds@reduce_dim_aux$prop_var_expl <- tissue.ife@reductions$pca@stdev
	
	ife.cds@int_colData@listData$reducedDims$UMAP <- tissue.ife@reductions$umap.harmony@cell.embeddings
	
	ife.cds@reduce_dim_aux@listData[["UMAP"]] <- tissue.ife@reductions[["umap.harmony"]]@cell.embeddings
	
	ife.cds@reduce_dim_aux$gene_loadings <- tissue.ife@reductions[["pca"]]@feature.loadings
	
	ife.cds <- cluster_cells(cds = ife.cds, 
	                         reduction_method = "UMAP")
	
	ife.cds <- learn_graph(ife.cds, 
	                       use_partition = T)
	
	# Use helper function to acquire all nodes associated with Proliferating_Keratinocyte_1 L2subtype
	
	get_earliest_principal_node <- function(ife.cds, L2subtype = 'Proliferating_Keratinocyte_1'){
	  
	  cell_ids <- which(colData(ife.cds)[, "L2subtype"] == L2subtype)
	  
	  closest_vertex <- ife.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
	  
	  closest_vertex <- as.matrix(closest_vertex[colnames(ife.cds), ])
	  
	  root_pr_nodes <- igraph::V(principal_graph(ife.cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
	  
	  root_pr_nodes
	}
	
	ife.cds <- order_cells(ife.cds, 
	                       root_pr_nodes = get_earliest_principal_node(ife.cds))
	
	plot_cells(cds = ife.cds,
	           color_cells_by = "pseudotime",
	           label_principal_points = T,
	           label_roots = F,
	           label_branch_points = F,
	           label_leaves = F,
	           cell_size = 1,
	           show_trajectory_graph = F)
	
	ggsave2(filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Pseudotime Blank UMAP.svg",  
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	plot_cells(cds = ife.cds,
	           color_cells_by = "pseudotime",
	           label_roots = F,
	           label_branch_points = F,
	           label_leaves = F,
	           cell_size = 1,
	           show_trajectory_graph = T,
	           trajectory_graph_color = 'black',
	           trajectory_graph_segment_size = 1)
	
	ggsave2(filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Pseudotime Trajectory Line UMAP.svg",  
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	ife.cds <- estimate_size_factors(ife.cds)
	
	ife.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(tissue.ife[["RNA3"]])
	
	saveRDS(ife.cds, 
	        '/data/overmilleram/scRNAseq/Skin & Oral/ife.cds.rds',
	        compress = F)
	
	#Genes that change as a function of pseudotime
	#Took 1.5 hours to run
	
	tissue.ife <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/keraife.rds')
	ife.cds <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ife.cds.rds')
	
	ife_cds_pr_test_res <- graph_test(ife.cds, 
	                                  neighbor_graph = "principal_graph",
	                                  cores = 4)
	
	saveRDS(ife_cds_pr_test_res, 
	        "ifepseudotime.rds",
	        compress = FALSE)
	
	ife_cds_pr_test_res <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ifepseudotime.rds')
	
	pr_deg_ids <- row.names(subset(ife_cds_pr_test_res, q_value < 0.0001))
	
	ife.cds <- preprocess_cds(ife.cds) # need to do for find_gene_modules() to work
	 
	ife_gene_module_df <- find_gene_modules(ife.cds[pr_deg_ids,], 
	                                        resolution = 0.001) 
	
	ife_cell_group_df <- tibble(cell=row.names(colData(ife.cds)), 
	                            cell_group=colData(ife.cds)$L1subtype)
	
	ife_agg_mat <- aggregate_gene_expression(ife.cds, 
	                                         ife_gene_module_df, 
	                                         ife_cell_group_df)
	
	row.names(ife_agg_mat) <- stringr::str_c("Module ", 
	                                         row.names(ife_agg_mat))
	
	pheatmap <- pheatmap::pheatmap(ife_agg_mat,
	                               color = paletteer_c("scico::bam", 100),
	                               scale = "column", 
	                               cluster_rows = T,
	                               cluster_cols = T,
	                               fontsize = 8,
	                               clustering_method = "ward.D2")
	
	ggsave2(pheatmap,
	        filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte IFE Pseudotime Modules Heatmap.svg",  
	        width = 1500,
	        height = 6000,
	        units = "px")
	
	ife.pseudo.genes <- ife_cds_pr_test_res[order(-ife_cds_pr_test_res$morans_I), ]
	ife.pseudo.genes <- subset(ife.pseudo.genes, q_value < 0.0001) 
	
	ife.pseudo.module <- ife_gene_module_df[order(ife_gene_module_df$module,
	                                              decreasing = F), ]
	
	##Make Excel sheet with markers and cell numbers
	#write.xlsx(ife.pseudo.genes,
	#           file = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte IFE Pseudotime.xlsx",
	#           sheetName = "Differentiated Genes Along Pseudotime",
	#           append = F)
	#write.xlsx(ife.pseudo.module,
	#           file = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte IFE Pseudotime.xlsx",
	#           sheetName = "Pseudotime Modules",
	#           append = T)
	
	#Extract pseudotime values and add to tissue.ife Seurat object metadata
	
	tissue.ife <- AddMetaData(object = tissue.ife,
	                        metadata = ife.cds@principal_graph_aux@listData$UMAP$pseudotime,
	                        col.name = "pseudotime")
	
	#FeaturePlot(tissue.ife, 
	#            features = "pseudotime",
	#            pt.size = 0.5) & scale_color_viridis_c() #worked
	#
	##Extract genes from modules for each specific simple subtype, run enrichR and generate Seurat module score
	#
	#ife.pseudo.module <- read.xlsx('/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte IFE Pseudotime.xlsx',
	#                               sheetIndex = 2,
	#                               header = T)
	#
	#ife.pseudo.module <- ife.pseudo.module[-1]
	
	ife_agg_mat_fil <- ife_agg_mat
	
	ife_agg_mat_fil <- ife_agg_mat_fil %>% as.data.frame() %>% rownames_to_column() 
	
	ife_agg_mat_fil[ ,1] <- gsub('Module ', '', ife_agg_mat_fil[ ,1]) %>% as.numeric(ife_agg_mat_fil[ ,1])
	
	epi.basal.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Basal > 0.1]
	ora.basal.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Basal > 0.1]
	prolifera.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Proliferating_Keratinocyte > 0.1]
	epi.supra1.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Suprabasal_1 > 0.1]
	epi.supra2.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Suprabasal_2 > 0.1]
	ora.supra1.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Suprabasal_1 > 0.1]
	ora.supra2.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Suprabasal_2 > 0.1]
	epi.supra.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Suprabasal_1 > 0.1 | ife_agg_mat_fil$Epithelial_Suprabasal_2 > 0.1]
	ora.supra.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Suprabasal_1 > 0.1 | ife_agg_mat_fil$Oral_Suprabasal_2 > 0.1]
	
	epi.basal.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.basal.module, ]
	ora.basal.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.basal.module, ]
	prolifera.genes <- ife.pseudo.module[ife.pseudo.module$module == prolifera.module, ]
	epi.supra1.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.supra1.module, ]
	epi.supra2.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.supra2.module, ]
	ora.supra1.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.supra1.module, ]
	ora.supra2.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.supra2.module, ]
	epi.supra.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.supra.module, ]
	ora.supra.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.supra.module, ]
	
	epi.basal.genes <- as.character(epi.basal.genes$id)
	ora.basal.genes <- as.character(ora.basal.genes$id)
	prolifera.genes <- as.character(prolifera.genes$id)
	epi.supra1.genes <- as.character(epi.supra1.genes$id)
	epi.supra2.genes <- as.character(epi.supra2.genes$id)
	ora.supra1.genes <- as.character(ora.supra1.genes$id)
	ora.supra2.genes <- as.character(ora.supra2.genes$id)
	epi.supra.genes <- as.character(epi.supra.genes$id)
	ora.supra.genes <- as.character(ora.supra.genes$id)
	
	#run enrichR
	##############
	setEnrichrSite("Enrichr") # Human/mouse genes
	
	websiteLive = T
	
	all.dbs = listEnrichrDbs()
	
	dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", 
	        'Reactome_2022', 'Tabula_Muris', 'Mouse_Gene_Atlas', 'Jensen_TISSUES', 'KEGG_2019_Mouse', 'MSigDB_Hallmark_2020', 
	        'WikiPathways_2019_Mouse', 'Panther_2016', 'BioPlex_2017', 'ChEA_2022', 'ENCODE_and_ChEA_Consensus_TFs_from_CHIP-X', 
	        'ENCODE_TF_ChIP-seq_2015', 'Epigenomics_Roadmap_HM_ChIP-seq', 'Transcription_Factor_PPIs', 
	        'Tissue_Protein_Expression_from_ProteomicsDB', 'TRRUST_Transcription_Factors_2019', 'KOMP2_Mouse_Phenotypes_2022', 
	        'miRTarBase_2017')
	
	epi.basal.enrichr <- enrichr(epi.basal.genes, dbs)
	ora.basal.enrichr <- enrichr(ora.basal.genes, dbs)
	prolifera.enrichr <- enrichr(prolifera.genes, dbs)
	epi.supra.enrichr <- enrichr(epi.supra.genes, dbs)
	ora.supra.enrichr <- enrichr(ora.supra.genes, dbs)
	
	require(dplyr)
	epi.basal.enrichr.plots = list()
	
	for (i in (1:21)) {
	  tryCatch({
	    epi.basal.enrichr.plots[[i]] <- plotEnrich(epi.basal.enrichr[[i]],
	                                           showTerms = 20, 
	                                           numChar = 100,
	                                           y = 'Ratio', 
	                                           orderBy = 'Combined.Score',
	                                           title = paste0(dbs[i],' - epi.basal.Module'))
	    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(epi.basal.enrichr.plots) <- dbs
	
	(epi.basal.enrichr.plots[[1]] + epi.basal.enrichr.plots[[2]] + epi.basal.enrichr.plots[[3]])/
	  (epi.basal.enrichr.plots[[4]] + epi.basal.enrichr.plots[[5]] + epi.basal.enrichr.plots[[6]])/
	  (epi.basal.enrichr.plots[[7]] + epi.basal.enrichr.plots[[8]] + epi.basal.enrichr.plots[[9]])/
	  (epi.basal.enrichr.plots[[10]] + epi.basal.enrichr.plots[[11]] + epi.basal.enrichr.plots[[12]])/
	  (epi.basal.enrichr.plots[[13]] + epi.basal.enrichr.plots[[14]] + epi.basal.enrichr.plots[[15]])/
	  (epi.basal.enrichr.plots[[16]] + epi.basal.enrichr.plots[[17]] + epi.basal.enrichr.plots[[18]])/
	  (epi.basal.enrichr.plots[[19]] + epi.basal.enrichr.plots[[20]] + epi.basal.enrichr.plots[[21]]) 
	
	ggsave('/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Epithelial Basal Gene Module enrichR.svg',
	       device = 'svg',
	       units = 'px',
	       height = 7500,
	       width = 10000,
	       limitsize = F)
	
	ora.basal.enrichr.plots = list()
	
	for (i in (1:21)) {
	  tryCatch({
	    ora.basal.enrichr.plots[[i]] <- plotEnrich(ora.basal.enrichr[[i]],
	                                           showTerms = 20, 
	                                           numChar = 100,
	                                           y = 'Ratio', 
	                                           orderBy = 'Combined.Score',
	                                           title = paste0(dbs[i],' - ora.basal.Module'))
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(ora.basal.enrichr.plots) <- dbs
	
	(ora.basal.enrichr.plots[[1]] + ora.basal.enrichr.plots[[2]] + ora.basal.enrichr.plots[[3]])/
	  (ora.basal.enrichr.plots[[4]] + ora.basal.enrichr.plots[[5]] + ora.basal.enrichr.plots[[6]])/
	  (ora.basal.enrichr.plots[[7]] + ora.basal.enrichr.plots[[8]] + ora.basal.enrichr.plots[[9]])/
	  (ora.basal.enrichr.plots[[10]] + ora.basal.enrichr.plots[[11]] + ora.basal.enrichr.plots[[12]])/
	  (ora.basal.enrichr.plots[[13]] + ora.basal.enrichr.plots[[14]] + ora.basal.enrichr.plots[[15]])/
	  (ora.basal.enrichr.plots[[16]] + ora.basal.enrichr.plots[[17]] + ora.basal.enrichr.plots[[18]])/
	  (ora.basal.enrichr.plots[[19]] + ora.basal.enrichr.plots[[20]] + ora.basal.enrichr.plots[[21]]) 
	
	ggsave('/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Oral Basal Gene Module enrichR.svg',
	       device = 'svg',
	       units = 'px',
	       height = 7500,
	       width = 10000,
	       limitsize = F)
	
	prolifera.enrichr.plots = list()
	
	for (i in (1:21)) {
	  tryCatch({
	    prolifera.enrichr.plots[[i]] <- plotEnrich(prolifera.enrichr[[i]],
	                                           showTerms = 20, 
	                                           numChar = 100,
	                                           y = 'Ratio', 
	                                           orderBy = 'Combined.Score',
	                                           title = paste0(dbs[i],' - prolifera.Module'))
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(prolifera.enrichr.plots) <- dbs
	
	(prolifera.enrichr.plots[[1]] + prolifera.enrichr.plots[[2]] + prolifera.enrichr.plots[[3]])/
	  (prolifera.enrichr.plots[[4]] + prolifera.enrichr.plots[[5]] + prolifera.enrichr.plots[[6]])/
	  (prolifera.enrichr.plots[[7]] + prolifera.enrichr.plots[[8]] + prolifera.enrichr.plots[[9]])/
	  (prolifera.enrichr.plots[[10]] + prolifera.enrichr.plots[[11]] + prolifera.enrichr.plots[[12]])/
	  (prolifera.enrichr.plots[[13]] + prolifera.enrichr.plots[[14]] + prolifera.enrichr.plots[[15]])/
	  (prolifera.enrichr.plots[[16]] + prolifera.enrichr.plots[[17]] + prolifera.enrichr.plots[[18]])/
	  (prolifera.enrichr.plots[[19]] + prolifera.enrichr.plots[[20]] + prolifera.enrichr.plots[[21]]) 
	
	ggsave('/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Proliferating Gene Module enrichR.svg',
	       device = 'svg',
	       units = 'px',
	       height = 7500,
	       width = 10000,
	       limitsize = F)
	
	epi.supra.enrichr.plots = list()
	
	for (i in (1:21)) {
	  tryCatch({
	    epi.supra.enrichr.plots[[i]] <- plotEnrich(epi.supra.enrichr[[i]],
	                                               showTerms = 20, 
	                                               numChar = 100,
	                                               y = 'Ratio', 
	                                               orderBy = 'Combined.Score',
	                                               title = paste0(dbs[i],' - epi.supra.Module'))
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(epi.supra.enrichr.plots) <- dbs
	
	(epi.supra.enrichr.plots[[1]] + epi.supra.enrichr.plots[[2]] + epi.supra.enrichr.plots[[3]])/
	  (epi.supra.enrichr.plots[[4]] + epi.supra.enrichr.plots[[5]] + epi.supra.enrichr.plots[[6]])/
	  (epi.supra.enrichr.plots[[7]] + epi.supra.enrichr.plots[[8]] + epi.supra.enrichr.plots[[9]])/
	  (epi.supra.enrichr.plots[[10]] + epi.supra.enrichr.plots[[11]] + epi.supra.enrichr.plots[[12]])/
	  (epi.supra.enrichr.plots[[13]] + epi.supra.enrichr.plots[[14]] + epi.supra.enrichr.plots[[15]])/
	  (epi.supra.enrichr.plots[[16]] + epi.supra.enrichr.plots[[17]] + epi.supra.enrichr.plots[[18]])/
	  (epi.supra.enrichr.plots[[19]] + epi.supra.enrichr.plots[[20]] + epi.supra.enrichr.plots[[21]]) 
	ggsave('/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Epithelial Suprabasal Gene Module enrichR.svg',
	       device = 'svg',
	       units = 'px',
	       height = 7500,
	       width = 10000,
	       limitsize = F)
	
	ora.supra.enrichr.plots = list()
	
	for (i in (1:21)) {
	  tryCatch({
	    ora.supra.enrichr.plots[[i]] <- plotEnrich(ora.supra.enrichr[[i]],
	                                               showTerms = 20, 
	                                               numChar = 100,
	                                               y = 'Ratio', 
	                                               orderBy = 'Combined.Score',
	                                               title = paste0(dbs[i],' - ora.supra.Module'))
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(ora.supra.enrichr.plots) <- dbs
	
	(ora.supra.enrichr.plots[[1]] + ora.supra.enrichr.plots[[2]] + ora.supra.enrichr.plots[[3]])/
	  (ora.supra.enrichr.plots[[4]] + ora.supra.enrichr.plots[[5]] + ora.supra.enrichr.plots[[6]])/
	  (ora.supra.enrichr.plots[[7]] + ora.supra.enrichr.plots[[8]] + ora.supra.enrichr.plots[[9]])/
	  (ora.supra.enrichr.plots[[10]] + ora.supra.enrichr.plots[[11]] + ora.supra.enrichr.plots[[12]])/
	  (ora.supra.enrichr.plots[[13]] + ora.supra.enrichr.plots[[14]] + ora.supra.enrichr.plots[[15]])/
	  (ora.supra.enrichr.plots[[16]] + ora.supra.enrichr.plots[[17]] + ora.supra.enrichr.plots[[18]])/
	  (ora.supra.enrichr.plots[[19]] + ora.supra.enrichr.plots[[20]] + ora.supra.enrichr.plots[[21]]) 
	ggsave('/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Oral Suprabasal Gene Module enrichR.svg',
	       device = 'svg',
	       units = 'px',
	       height = 7500,
	       width = 10000,
	       limitsize = F)
	
	##############
	
	epi.basal.genes <- list(epi.basal.genes)
	ora.basal.genes <- list(ora.basal.genes)
	prolifera.genes <- list(prolifera.genes)
	epi.supra1.genes <- list(epi.supra1.genes)
	epi.supra2.genes <- list(epi.supra2.genes)
	ora.supra1.genes <- list(ora.supra1.genes) 
	ora.supra2.genes <- list(ora.supra2.genes) 
	epi.supra.genes <- list(epi.supra.genes)
	ora.supra.genes <- list(ora.supra.genes) # genes must be in a list for AddModuleScore to work
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = epi.basal.genes,
	                             ctrl = 100,
	                             name = 'Epi.Basal.Score',
	                             assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = ora.basal.genes,
	                             ctrl = 100,
	                             name = 'Oral.Basal.Score',
	                             assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                           features = prolifera.genes,
	                           ctrl = 100,
	                           name = 'Proliferating.Kera.Score',
	                           assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = epi.supra1.genes,
	                             ctrl = 100,
	                             name = 'Epi.Suprabasal_1.Score',
	                             assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = epi.supra2.genes,
	                             ctrl = 100,
	                             name = 'Epi.Suprabasal_2.Score',
	                             assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                           features = ora.supra1.genes,
	                           ctrl = 100,
	                           name = 'Oral.Suprabasal_1.Score',
	                           assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = ora.supra2.genes,
	                             ctrl = 100,
	                             name = 'Oral.Suprabasal_2.Score',
	                             assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = epi.supra.genes,
	                             ctrl = 100,
	                             name = 'Combined.Epi.Suprabasal.Score',
	                             assay = 'RNA')
	
	tissue.ife <- AddModuleScore(tissue.ife,
	                             features = ora.supra.genes,
	                             ctrl = 100,
	                             name = 'Combined.Oral.Suprabasal.Score',
	                             assay = 'RNA')
	
	saveRDS(tissue.ife,
	        '/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds',
	        compress = F)
	
	Idents(tissue.ife) <- 'L1subtype'
	DefaultAssay(tissue.ife) <- 'RNA'
	
	tissue.ife.basal <- subset(tissue.ife,
	                           idents = c('Epithelial_Basal', 'Oral_Basal', 'Proliferating_Keratinocyte')) 
	
	tissue.ife.supra <- subset(tissue.ife,
	                           idents = c('Epithelial_Suprabasal_1', 'Epithelial_Suprabasal_2', 'Oral_Suprabasal_1', 'Oral_Suprabasal_2')) 
	
	Idents(tissue.ife) <- 'celltype'
	Idents(tissue.ife.basal) <- 'celltype'
	Idents(tissue.ife.supra) <- 'celltype'
	
	VlnPlot(tissue.ife,
	        features = c('Oral.Suprabasal_1.Score1', 'Oral.Suprabasal_2.Score1', 'Epi.Suprabasal_1.Score1', 'Epi.Suprabasal_2.Score1', 
	                     'Combined.Oral.Suprabasal.Score1', 'Combined.Epi.Suprabasal.Score1',
	                     'Epi.Basal.Score1', 'Oral.Basal.Score1', 'Proliferating.Kera.Score1'),
	        pt.size = 0,
	        flip = T,
	        stack = T,
	        cols = c('#FB4142','blue','green'),
	        split.by = 'state',
	        raster = F) + NoLegend()
	
	ggsave2(filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte IFE Pseudotime Gene Module (All Modules) VlnPlots.svg",  
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	VlnPlot(tissue.ife.basal,
	        features = c('Epi.Basal.Score1', 'Oral.Basal.Score1', 'Proliferating.Kera.Score1'),
	        pt.size = 0,
	        flip = T,
	        stack = T,
	        cols = c('#FB4142','blue','green'),
	        split.by = 'state',
	        raster = F) #+ NoLegend()
	
	ggsave2(filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Basal IFE Pseudotime Gene Module (Basal Modules) VlnPlots.svg",  
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	VlnPlot(tissue.ife.supra,
	        features = c('Oral.Suprabasal_1.Score1', 'Oral.Suprabasal_2.Score1', 'Epi.Suprabasal_1.Score1', 'Epi.Suprabasal_2.Score1'),
	        pt.size = 0,
	        flip = T,
	        stack = T,
	        cols = c('#FB4142','blue','green'),
	        split.by = 'state',
	        raster = F) #+ NoLegend()
	
	ggsave2(filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Suprabasal IFE Pseudotime Gene Module (Suprabasal Modules) VlnPlots.svg",  
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	VlnPlot(tissue.ife.supra,
	        features = c('Combined.Oral.Suprabasal.Score1', 'Combined.Epi.Suprabasal.Score1'),
	        pt.size = 0,
	        flip = T,
	        stack = T,
	        cols = c('#FB4142','blue','green'),
	        split.by = 'state',
	        raster = F) #+ NoLegend()
	
	ggsave2(filename = "/data/overmilleram/scRNAseq/Skin & Oral/Keratinocyte Suprabasal IFE Pseudotime Gene Module (Combined Suprabasal Modules) VlnPlots.svg",  
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	
	# show select gene expression as a function of pseudotime
	
	pseudo_genes <- c('Tgm3', 'Krt77', 'Krt4', 'Krt6b', 'Sbsn', 'Krt5')
	
	skco_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('skin_control_8_week')]
	
	skdo_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('skin_dox_8_week')]
	
	orbu_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('control_buccal_mucosa')]
	
	skco <- plot_genes_in_pseudotime(skco_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	skdo <- plot_genes_in_pseudotime(skdo_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	orbu <- plot_genes_in_pseudotime(orbu_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	
	#extract module scores to perform Kruskal-Wallis
	
	module.score <- data.frame(EpiSuprabasalCombined = tissue.ife$Combined.Epi.Suprabasal.Score1, 
	                           OralSuprabasalCombined = tissue.ife$Combined.Oral.Suprabasal.Score1, 
	                           EpiSuprabasal1 = tissue.ife$Epi.Suprabasal_1.Score1,
	                           EpiSuprabasal2 = tissue.ife$Epi.Suprabasal_2.Score1,
	                           OralSuprabasal1 = tissue.ife$Oral.Suprabasal_1.Score1,
	                           OralSuprabasal2 = tissue.ife$Oral.Suprabasal_2.Score1,
	                           state = tissue.ife$state)
	
	ggpubr::ggboxplot(module.score, #check distrubution on a box plot
	                  x = "state", 
	                  y = "OralSuprabasal1", 
	                  color = "state", 
	                  palette = c('#3D79F3FF','#34A74BFF','#F9B90AFF'),
	                  order = c('skin_control_8_week', 'skin_dox_8_week', 'control_buccal_mucosa'),
	                  ylab = "Module Score", 
	                  xlab = F)
	
	pairwise.wilcox.test(module.score$EpiSuprabasal1, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-09 at least)
	pairwise.wilcox.test(module.score$EpiSuprabasal2, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-09 at least)
	pairwise.wilcox.test(module.score$OralSuprabasal1, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	pairwise.wilcox.test(module.score$OralSuprabasal2, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 6.4e-5 at least)
	
	pairwise.wilcox.test(module.score$EpiSuprabasalCombined, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-09 at least)
	pairwise.wilcox.test(module.score$OralSuprabasalCombined, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-09 at least)

### Wound Skin 

	library(dplyr)
	library(ggplot2)
	library(Matrix)
	library(patchwork)
	library(gdata)
	library(Seurat)
	library(SeuratObject)
	library(SeuratWrappers)
	library(BiocManager)
	library(metap)
	library(cowplot)
	library(sctransform)
	library(xlsx)
	library(glmGamPoi)
	library(clustree)
	library(biomaRt)
	library(monocle3)
	library(magrittr)
	library(slingshot)
	library(paletteer)
	library(nichenetr)
	library(multinichenetr)
	library(tidyverse)
	library(circlize)
	library(scales)
	library(kableExtra)
	library(knitr)
	library(SoupX)
	library(scDblFinder)
	library(viridis)
	library(CellChat)
	library(uwot)
	library(ComplexHeatmap)
	library(enrichR)
	
	#Set current working directory
	setwd("/data/overmilleram/scRNAseq/Wound/")
	################
	
	#Preprocess data
	
	#Use SoupX to remove ambient RNA
	toc.fcs1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_1/filtered_feature_bc_matrix/")
	tod.fcs1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_1/raw_feature_bc_matrix/")
	toc.fcs2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_2/filtered_feature_bc_matrix/")
	tod.fcs2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_2/raw_feature_bc_matrix/")
	toc.fcs3 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_3/filtered_feature_bc_matrix/")
	tod.fcs3 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_3/raw_feature_bc_matrix/")
	toc.fcs4 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_4/filtered_feature_bc_matrix/")
	tod.fcs4 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Control_D4_4/raw_feature_bc_matrix/")
	toc.mcs1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Control_D4_1/filtered_feature_bc_matrix/")
	tod.mcs1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Control_D4_1/raw_feature_bc_matrix/")
	toc.mcs2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Control_D4_2/filtered_feature_bc_matrix/")
	tod.mcs2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Control_D4_2/raw_feature_bc_matrix/")
	toc.fds1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_1/filtered_feature_bc_matrix/")
	tod.fds1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_1/raw_feature_bc_matrix/")
	toc.fds2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_2/filtered_feature_bc_matrix/")
	tod.fds2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_2/raw_feature_bc_matrix/")
	toc.fds3 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_3/filtered_feature_bc_matrix/")
	tod.fds3 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_3/raw_feature_bc_matrix/")
	toc.fds4 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_4/filtered_feature_bc_matrix/")
	tod.fds4 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Female_Dox_D4_4/raw_feature_bc_matrix/")
	toc.mds1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Dox_D4_1/filtered_feature_bc_matrix/")
	tod.mds1 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Dox_D4_1/raw_feature_bc_matrix/")
	toc.mds2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Dox_D4_2/filtered_feature_bc_matrix/")
	tod.mds2 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Dox_D4_2/raw_feature_bc_matrix/")
	toc.mds3 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Dox_D4_3/filtered_feature_bc_matrix/")
	tod.mds3 <- Seurat::Read10X("../../Sample Matrices/D4 Skin/8w_Male_Dox_D4_3/raw_feature_bc_matrix/")
	
	sc.fcs1 <- SoupChannel(tod.fcs1,toc.fcs1, calcSoupProfile = T)
	sc.fcs2 <- SoupChannel(tod.fcs2,toc.fcs2, calcSoupProfile = T)
	sc.fcs3 <- SoupChannel(tod.fcs3,toc.fcs3, calcSoupProfile = T)
	sc.fcs4 <- SoupChannel(tod.fcs4,toc.fcs4, calcSoupProfile = T)
	sc.mcs1 <- SoupChannel(tod.mcs1,toc.mcs1, calcSoupProfile = T)
	sc.mcs2 <- SoupChannel(tod.mcs2,toc.mcs2, calcSoupProfile = T)
	sc.fds1 <- SoupChannel(tod.fds1,toc.fds1, calcSoupProfile = T)
	sc.fds2 <- SoupChannel(tod.fds2,toc.fds2, calcSoupProfile = T)
	sc.fds3 <- SoupChannel(tod.fds3,toc.fds3, calcSoupProfile = T)
	sc.fds4 <- SoupChannel(tod.fds4,toc.fds4, calcSoupProfile = T)
	sc.mds1 <- SoupChannel(tod.mds1,toc.mds1, calcSoupProfile = T)
	sc.mds2 <- SoupChannel(tod.mds2,toc.mds2, calcSoupProfile = T)
	sc.mds3 <- SoupChannel(tod.mds3,toc.mds3, calcSoupProfile = T)
	
	seur.fcs1 <- CreateSeuratObject(sc.fcs1$toc)
	seur.fcs2 <- CreateSeuratObject(sc.fcs2$toc)
	seur.fcs3 <- CreateSeuratObject(sc.fcs3$toc)
	seur.fcs4 <- CreateSeuratObject(sc.fcs4$toc)
	seur.mcs1 <- CreateSeuratObject(sc.mcs1$toc)
	seur.mcs2 <- CreateSeuratObject(sc.mcs2$toc)
	seur.fds1 <- CreateSeuratObject(sc.fds1$toc)
	seur.fds2 <- CreateSeuratObject(sc.fds2$toc)
	seur.fds3 <- CreateSeuratObject(sc.fds3$toc)
	seur.fds4 <- CreateSeuratObject(sc.fds4$toc)
	seur.mds1 <- CreateSeuratObject(sc.mds1$toc)
	seur.mds2 <- CreateSeuratObject(sc.mds2$toc)
	seur.mds3 <- CreateSeuratObject(sc.mds3$toc)
	
	seur.list <- list(seur.fcs1=seur.fcs1, seur.fcs2=seur.fcs2, seur.fcs3=seur.fcs3, seur.fcs4=seur.fcs4,
	                  seur.mcs1=seur.mcs1, seur.mcs2=seur.mcs2,
	                  seur.fds1=seur.fds1, seur.fds2=seur.fds2, seur.fds3=seur.fds3, seur.fds4=seur.fds4, 
	                  seur.mds1=seur.mds1, seur.mds2=seur.mds2, seur.mds3=seur.mds3)
	
	seur.list <- lapply(X = seur.list, 
	                        FUN = function(x) {
	                          x = NormalizeData(x)
	                          x = FindVariableFeatures(x, 
	                                                    selection.method = "vst",
	                                                    nfeatures = 5000)
	                          x = ScaleData(x, verbose = TRUE)
	                          x = RunPCA(x, 
	                                      npcs = 50, verbose = TRUE)
	                          x = FindNeighbors(x, 
	                                             reduction = "pca", 
	                                             dims = 1:50)
	                          x = FindClusters(x,
	                                            algorithm = 3,
	                                            resolution = 1)
	                        })
	
	list2env(seur.list, 
	         .GlobalEnv)
	
	sc.fcs1 <- setClusters(sc.fcs1, seur.fcs1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fcs2 <- setClusters(sc.fcs2, seur.fcs2$seurat_clusters) %>% autoEstCont() %>% adjustCounts() 
	sc.fcs3 <- setClusters(sc.fcs3, seur.fcs3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fcs4 <- setClusters(sc.fcs4, seur.fcs4$seurat_clusters) %>% autoEstCont() %>% adjustCounts() 
	sc.mcs1 <- setClusters(sc.mcs1, seur.mcs1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mcs2 <- setClusters(sc.mcs2, seur.mcs2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds1 <- setClusters(sc.fds1, seur.fds1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds2 <- setClusters(sc.fds2, seur.fds2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds3 <- setClusters(sc.fds3, seur.fds3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.fds4 <- setClusters(sc.fds4, seur.fds4$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mds1 <- setClusters(sc.mds1, seur.mds1$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mds2 <- setClusters(sc.mds2, seur.mds2$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	sc.mds3 <- setClusters(sc.mds3, seur.mds3$seurat_clusters) %>% autoEstCont() %>% adjustCounts()
	
	fcs1 <- CreateSeuratObject(sc.fcs1, min.cells = 3, min.features = 200)
	fcs2 <- CreateSeuratObject(sc.fcs2, min.cells = 3, min.features = 200)
	fcs3 <- CreateSeuratObject(sc.fcs3, min.cells = 3, min.features = 200)
	fcs4 <- CreateSeuratObject(sc.fcs4, min.cells = 3, min.features = 200)
	mcs1 <- CreateSeuratObject(sc.mcs1, min.cells = 3, min.features = 200)
	mcs2 <- CreateSeuratObject(sc.mcs2, min.cells = 3, min.features = 200)
	fds1 <- CreateSeuratObject(sc.fds1, min.cells = 3, min.features = 200)
	fds2 <- CreateSeuratObject(sc.fds2, min.cells = 3, min.features = 200)
	fds3 <- CreateSeuratObject(sc.fds3, min.cells = 3, min.features = 200)
	fds4 <- CreateSeuratObject(sc.fds4, min.cells = 3, min.features = 200)
	mds1 <- CreateSeuratObject(sc.mds1, min.cells = 3, min.features = 200)
	mds2 <- CreateSeuratObject(sc.mds2, min.cells = 3, min.features = 200)
	mds3 <- CreateSeuratObject(sc.mds3, min.cells = 3, min.features = 200)
	
	#Use scDblFinder to determine # of doublets
	
	set.seed(42)
	
	fcs1sce <- scDblFinder(GetAssayData(fcs1, slot = "counts"))
	fcs2sce <- scDblFinder(GetAssayData(fcs2, slot = "counts"))
	fcs3sce <- scDblFinder(GetAssayData(fcs3, slot = "counts"))
	fcs4sce <- scDblFinder(GetAssayData(fcs4, slot = "counts"))
	mcs1sce <- scDblFinder(GetAssayData(mcs1, slot = "counts"))
	mcs2sce <- scDblFinder(GetAssayData(mcs2, slot = "counts"))
	fds1sce <- scDblFinder(GetAssayData(fds1, slot = "counts"))
	fds2sce <- scDblFinder(GetAssayData(fds2, slot = "counts"))
	fds3sce <- scDblFinder(GetAssayData(fds3, slot = "counts"))
	fds4sce <- scDblFinder(GetAssayData(fds4, slot = "counts"))
	mds1sce <- scDblFinder(GetAssayData(mds1, slot = "counts"))
	mds2sce <- scDblFinder(GetAssayData(mds2, slot = "counts"))
	mds3sce <- scDblFinder(GetAssayData(mds3, slot = "counts"))
	
	fcs1$scDblFinder.class <- fcs1sce$scDblFinder.class
	fcs2$scDblFinder.class <- fcs2sce$scDblFinder.class
	fcs3$scDblFinder.class <- fcs3sce$scDblFinder.class
	fcs4$scDblFinder.class <- fcs4sce$scDblFinder.class
	mcs1$scDblFinder.class <- mcs1sce$scDblFinder.class
	mcs2$scDblFinder.class <- mcs2sce$scDblFinder.class
	fds1$scDblFinder.class <- fds1sce$scDblFinder.class
	fds2$scDblFinder.class <- fds2sce$scDblFinder.class
	fds3$scDblFinder.class <- fds3sce$scDblFinder.class
	fds4$scDblFinder.class <- fds4sce$scDblFinder.class
	mds1$scDblFinder.class <- mds1sce$scDblFinder.class
	mds2$scDblFinder.class <- mds2sce$scDblFinder.class
	mds3$scDblFinder.class <- mds3sce$scDblFinder.class
	
	Idents(fcs1) <- "scDblFinder.class"
	Idents(fcs2) <- "scDblFinder.class"
	Idents(fcs3) <- "scDblFinder.class"
	Idents(fcs4) <- "scDblFinder.class"
	Idents(mcs1) <- "scDblFinder.class"
	Idents(mcs2) <- "scDblFinder.class"
	Idents(fds1) <- "scDblFinder.class"
	Idents(fds2) <- "scDblFinder.class"
	Idents(fds3) <- "scDblFinder.class"
	Idents(fds4) <- "scDblFinder.class"
	Idents(mds1) <- "scDblFinder.class"
	Idents(mds2) <- "scDblFinder.class"
	Idents(mds3) <- "scDblFinder.class"
	
	fcs1 <- subset(fcs1, idents = "singlet")
	fcs2 <- subset(fcs2, idents = "singlet")
	fcs3 <- subset(fcs3, idents = "singlet")
	fcs4 <- subset(fcs4, idents = "singlet")
	mcs1 <- subset(mcs1, idents = "singlet")
	mcs2 <- subset(mcs2, idents = "singlet")
	fds1 <- subset(fds1, idents = "singlet")
	fds2 <- subset(fds2, idents = "singlet")
	fds3 <- subset(fds3, idents = "singlet")
	fds4 <- subset(fds4, idents = "singlet")
	mds1 <- subset(mds1, idents = "singlet")
	mds2 <- subset(mds2, idents = "singlet")
	mds3 <- subset(mds3, idents = "singlet")
	
	#File cleanup
	rm(fcs1sce, fcs2sce, fcs3sce, fcs4sce, fcs5sce, fcs6sce,
	   mcs1sce, mcs2sce, mcs3sce, mcs4sce, mcs5sce,
	   fds1sce, fds2sce, fds3sce, fds4sce, fds5sce, fds6sce,
	   mds1sce, mds2sce, mds3sce, mds4sce, mds5sce, mds6sce,
	   sc.fcs1, sc.fcs2, sc.fcs3, sc.fcs4, sc.fcs5, sc.fcs6,
	   sc.mcs1, sc.mcs2, sc.mcs3, sc.mcs4, sc.mcs5,
	   sc.fds1, sc.fds2, sc.fds3, sc.fds4, sc.fds5, sc.fds6,
	   sc.mds1, sc.mds2, sc.mds3, sc.mds4, sc.mds5, sc.mds6,
	   seur.fcs1, seur.fcs2, seur.fcs3, seur.fcs4, seur.fcs5, seur.fcs6,
	   seur.mcs1, seur.mcs2, seur.mcs3, seur.mcs4, seur.mcs5,
	   seur.fds1, seur.fds2, seur.fds3, seur.fds4, seur.fds5, seur.fds6,
	   seur.mds1, seur.mds2, seur.mds3, seur.mds4, seur.mds5, seur.mds6,
	   tod.fcs1, tod.fcs2, tod.fcs3, tod.fcs4, tod.fcs5, tod.fcs6,
	   tod.mcs1, tod.mcs2, tod.mcs3, tod.mcs4, tod.mcs5,
	   tod.fds1, tod.fds2, tod.fds3, tod.fds4, tod.fds5, tod.fds6,
	   tod.mds1, tod.mds2, tod.mds3, tod.mds4, tod.mds5, tod.mds6,
	   toc.fcs1, toc.fcs2, toc.fcs3, toc.fcs4, toc.fcs5, toc.fcs6,
	   toc.mcs1, toc.mcs2, toc.mcs3, toc.mcs4, toc.mcs5,
	   toc.fds1, toc.fds2, toc.fds3, toc.fds4, toc.fds5, toc.fds6,
	   toc.mds1, toc.mds2, toc.mds3, toc.mds4, toc.mds5, toc.mds6,
	   fcs1sce, fcs2sce, fcs3sce, fcs4sce, fcs5sce, fcs6sce,
	   mcs1sce, mcs2sce, mcs3sce, mcs4sce, mcs5sce,
	   fds1sce, fds2sce, fds3sce, fds4sce, fds5sce, fds6sce,
	   mds1sce, mds2sce, mds3sce, mds4sce, mds5sce, mds6sce,
	   fcbm1sce, fcbm2sce, fchp1sce, fchp2sce, mcbm1sce, mchp1sce,
	   sc.fcs1, sc.fcs2, sc.fcs3, sc.fcs4, sc.fcs5, sc.fcs6,
	   sc.mcs1, sc.mcs2, sc.mcs3, sc.mcs4, sc.mcs5,
	   sc.fds1, sc.fds2, sc.fds3, sc.fds4, sc.fds5, sc.fds6,
	   sc.mds1, sc.mds2, sc.mds3, sc.mds4, sc.mds5, sc.mds6,
	   sc.fcbm1, sc.fcbm2, sc.fchp1, sc.fchp2, sc.mcbm1, sc.mchp1,
	   seur.fcs1, seur.fcs2, seur.fcs3, seur.fcs4, seur.fcs5, seur.fcs6,
	   seur.mcs1, seur.mcs2, seur.mcs3, seur.mcs4, seur.mcs5,
	   seur.fds1, seur.fds2, seur.fds3, seur.fds4, seur.fds5, seur.fds6,
	   seur.mds1, seur.mds2, seur.mds3, seur.mds4, seur.mds5, seur.mds6,
	   seur.fcbm1, seur.fcbm2, seur.fchp1, seur.fchp2, seur.mcbm1, seur.mchp1,
	   tod.fcs1, tod.fcs2, tod.fcs3, tod.fcs4, tod.fcs5, tod.fcs6,
	   tod.mcs1, tod.mcs2, tod.mcs3, tod.mcs4, tod.mcs5,
	   tod.fds1, tod.fds2, tod.fds3, tod.fds4, tod.fds5, tod.fds6,
	   tod.mds1, tod.mds2, tod.mds3, tod.mds4, tod.mds5, tod.mds6,
	   tod.fcbm1, tod.fcbm2, tod.fchp1, tod.fchp2, tod.mcbm1, tod.mchp1,
	   toc.fcs1, toc.fcs2, toc.fcs3, toc.fcs4, toc.fcs5, toc.fcs6,
	   toc.mcs1, toc.mcs2, toc.mcs3, toc.mcs4, toc.mcs5,
	   toc.fds1, toc.fds2, toc.fds3, toc.fds4, toc.fds5, toc.fds6,
	   toc.mds1, toc.mds2, toc.mds3, toc.mds4, toc.mds5, toc.mds6, 
	   toc.fcbm1, toc.fcbm2, toc.fchp1, toc.fchp2, toc.mcbm1, toc.mchp1, 
	   fcbm1sce, fcbm2sce, fcbm3sce, fchp1sce, fchp2sce, fchp3sce, mcbm1sce, mcbm2sce, mchp1sce, mchp2sce,
	   fcbm1.desouped, fcbm2.desouped, fcbm3.desouped, fchp1.desouped, fchp2.desouped, fchp3.desouped, mcbm1.desouped, mcbm2.desouped, mchp1.desouped, mchp2.desouped,
	   sc.fcbm1, sc.fcbm2, sc.fcbm3, sc.fchp1, sc.fchp2, sc.fchp3, sc.mcbm1, sc.mcbm2, sc.mchp1, sc.mchp2,
	   seur.fcbm1, seur.fcbm2, seur.fcbm3, seur.fchp1, seur.fchp2, seur.fchp3, seur.mcbm1, seur.mcbm2, seur.mchp1, seur.mchp2,
	   tod.fcbm1, tod.fcbm2, tod.fcbm3, tod.fchp1, tod.fchp2, tod.fchp3, tod.mcbm1, tod.mcbm2, tod.mchp1, tod.mchp2,
	   toc.fcbm1, toc.fcbm2, toc.fcbm3, toc.fchp1, toc.fchp2, toc.fchp3, toc.mcbm1, toc.mcbm2, toc.mchp1, toc.mchp2)
	
	#Add condition identity
	fcs1$condition <- "control"
	fcs2$condition <- "control"
	fcs3$condition <- "control"
	fcs4$condition <- "control"
	mcs1$condition <- "control"
	mcs2$condition <- "control"
	fds1$condition <- "dox"
	fds2$condition <- "dox"
	fds3$condition <- "dox"
	fds4$condition <- "dox"
	mds1$condition <- "dox"
	mds2$condition <- "dox"
	mds3$condition <- "dox"
	
	#Add tissue ID
	fcs1$tissue <- "wound"
	fcs2$tissue <- "wound"
	fcs3$tissue <- "wound"
	fcs4$tissue <- "wound"
	mcs1$tissue <- "wound"
	mcs2$tissue <- "wound"
	fds1$tissue <- "wound"
	fds2$tissue <- "wound"
	fds3$tissue <- "wound"
	fds4$tissue <- "wound"
	mds1$tissue <- "wound"
	mds2$tissue <- "wound"
	mds3$tissue <- "wound"
	
	#Add Pitx1 induction identity
	fcs1$induction <- "8_week_skin"
	fcs2$induction <- "8_week_skin"
	fcs3$induction <- "8_week_skin"
	fcs4$induction <- "8_week_skin"
	mcs1$induction <- "8_week_skin"
	mcs2$induction <- "8_week_skin"
	fds1$induction <- "8_week_skin"
	fds2$induction <- "8_week_skin"
	fds3$induction <- "8_week_skin"
	fds4$induction <- "8_week_skin"
	mds1$induction <- "8_week_skin"
	mds2$induction <- "8_week_skin"
	mds3$induction <- "8_week_skin"
	
	#Add skin state
	fcs1$state <- "wound_control_8_week"
	fcs2$state <- "wound_control_8_week"
	fcs3$state <- "wound_control_8_week"
	fcs4$state <- "wound_control_8_week"
	mcs1$state <- "wound_control_8_week"
	mcs2$state <- "wound_control_8_week"
	fds1$state <- "wound_dox_8_week"
	fds2$state <- "wound_dox_8_week"
	fds3$state <- "wound_dox_8_week"
	fds4$state <- "wound_dox_8_week"
	mds1$state <- "wound_dox_8_week"
	mds2$state <- "wound_dox_8_week"
	mds3$state <- "wound_dox_8_week"
	
	#Add sex identity
	fcs1$sex <- "female"
	fcs2$sex <- "female"
	fcs3$sex <- "female"
	fcs4$sex <- "female"
	mcs1$sex <- "male"
	mcs2$sex <- "male"
	fds1$sex <- "female"
	fds2$sex <- "female"
	fds3$sex <- "female"
	fds4$sex <- "female"
	mds1$sex <- "male"
	mds2$sex <- "male"
	mds3$sex <- "male"
	
	#Add batch #
	fcs1$batch <- "1"
	fcs2$batch <- "1"
	fcs3$batch <- "2"
	fcs4$batch <- "2"
	mcs1$batch <- "1"
	mcs2$batch <- "2"
	fds1$batch <- "1"
	fds2$batch <- "1"
	fds3$batch <- "2"
	fds4$batch <- "2"
	mds1$batch <- "1"
	mds2$batch <- "2"
	mds3$batch <- "2"
	
	#Add ID#
	fcs1$sampleid <- "fcs1"
	fcs2$sampleid <- "fcs2"
	fcs3$sampleid <- "fcs3"
	fcs4$sampleid <- "fcs4"
	mcs1$sampleid <- "mcs1"
	mcs2$sampleid <- "mcs2"
	fds1$sampleid <- "fds1"
	fds2$sampleid <- "fds2"
	fds3$sampleid <- "fds3"
	fds4$sampleid <- "fds4"
	mds1$sampleid <- "mds1"
	mds2$sampleid <- "mds2"
	mds3$sampleid <- "mds3"
	
	#Add %mitochondrial gene to metadata, filter based on nFeature_RNA & percent.mt, normalize data, 
	fcs1[["percent.mt"]] <- PercentageFeatureSet(fcs1, pattern = "^mt-")
	fcs2[["percent.mt"]] <- PercentageFeatureSet(fcs2, pattern = "^mt-")
	fcs3[["percent.mt"]] <- PercentageFeatureSet(fcs3, pattern = "^mt-")
	fcs4[["percent.mt"]] <- PercentageFeatureSet(fcs4, pattern = "^mt-")
	mcs1[["percent.mt"]] <- PercentageFeatureSet(mcs1, pattern = "^mt-")
	mcs2[["percent.mt"]] <- PercentageFeatureSet(mcs2, pattern = "^mt-")
	fds1[["percent.mt"]] <- PercentageFeatureSet(fds1, pattern = "^mt-")
	fds2[["percent.mt"]] <- PercentageFeatureSet(fds2, pattern = "^mt-")
	fds3[["percent.mt"]] <- PercentageFeatureSet(fds3, pattern = "^mt-")
	fds4[["percent.mt"]] <- PercentageFeatureSet(fds4, pattern = "^mt-")
	mds1[["percent.mt"]] <- PercentageFeatureSet(mds1, pattern = "^mt-")
	mds2[["percent.mt"]] <- PercentageFeatureSet(mds2, pattern = "^mt-")
	mds3[["percent.mt"]] <- PercentageFeatureSet(mds3, pattern = "^mt-")
	
	fcs1[["percent.ribo"]] <- PercentageFeatureSet(fcs1, pattern = "Rp[sl]")
	fcs2[["percent.ribo"]] <- PercentageFeatureSet(fcs2, pattern = "Rp[sl]")
	fcs3[["percent.ribo"]] <- PercentageFeatureSet(fcs3, pattern = "Rp[sl]")
	fcs4[["percent.ribo"]] <- PercentageFeatureSet(fcs4, pattern = "Rp[sl]")
	mcs1[["percent.ribo"]] <- PercentageFeatureSet(mcs1, pattern = "Rp[sl]")
	mcs2[["percent.ribo"]] <- PercentageFeatureSet(mcs2, pattern = "Rp[sl]")
	fds1[["percent.ribo"]] <- PercentageFeatureSet(fds1, pattern = "Rp[sl]")
	fds2[["percent.ribo"]] <- PercentageFeatureSet(fds2, pattern = "Rp[sl]")
	fds3[["percent.ribo"]] <- PercentageFeatureSet(fds3, pattern = "Rp[sl]")
	fds4[["percent.ribo"]] <- PercentageFeatureSet(fds4, pattern = "Rp[sl]")
	mds1[["percent.ribo"]] <- PercentageFeatureSet(mds1, pattern = "Rp[sl]")
	mds2[["percent.ribo"]] <- PercentageFeatureSet(mds2, pattern = "Rp[sl]")
	mds3[["percent.ribo"]] <- PercentageFeatureSet(mds3, pattern = "Rp[sl]")
	
	fcs1 <- subset(fcs1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fcs2 <- subset(fcs2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5)
	fcs3 <- subset(fcs3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5)
	fcs4 <- subset(fcs4, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mcs1 <- subset(mcs1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mcs2 <- subset(mcs2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds1 <- subset(fds1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds2 <- subset(fds2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds3 <- subset(fds3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	fds4 <- subset(fds4, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds1 <- subset(mds1, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds2 <- subset(mds2, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	mds3 <- subset(mds3, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 5) 
	
	wound.list <- list(fcs1 = fcs1, fcs2 = fcs2, fcs3 = fcs3, fcs4 = fcs4,
	                     mcs1 = mcs1, mcs2 = mcs2,
	                     fds1 = fds1, fds2 = fds2, fds3 = fds3, fds4 = fds4,
	                     mds1 = mds1, mds2 = mds2, mds3 = mds3)
	
	saveRDS(wound.list,
	        file = "woundskinlist.rds",
	        compress = T)
	
	#File cleanup
	rm(fcs1, fcs2, fcs3, fcs4,
	   mcs1, mcs2,
	   fds1, fds2, fds3, fds4,
	   mds1, mds2, mds3, seur.list)
	
	# RPCA integration
	########
	
	wound.list <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskinlist.rds')
	
	m.s.genes <- cc.genes$s.genes %>% convert_human_to_mouse_symbols
	m.g2m.genes <- cc.genes$g2m.genes %>% convert_human_to_mouse_symbols
	
	fcs1 <- wound.list$fcs1 
	fcs2 <- wound.list$fcs2
	fcs3 <- wound.list$fcs3
	fcs4 <- wound.list$fcs4
	mcs1 <- wound.list$mcs1
	mcs2 <- wound.list$mcs2
	fds1 <- wound.list$fds1
	fds2 <- wound.list$fds2
	fds3 <- wound.list$fds3
	fds4 <- wound.list$fds4
	mds1 <- wound.list$mds1
	mds2 <- wound.list$mds2
	mds3 <- wound.list$mds3
	
	fcs1[['RNA']] <- as(object = fcs1[['RNA']], Class = 'Assay')
	fcs2[['RNA']] <- as(object = fcs2[['RNA']], Class = 'Assay')
	fcs3[['RNA']] <- as(object = fcs3[['RNA']], Class = 'Assay')
	fcs4[['RNA']] <- as(object = fcs4[['RNA']], Class = 'Assay')
	mcs1[['RNA']] <- as(object = mcs1[['RNA']], Class = 'Assay')
	mcs2[['RNA']] <- as(object = mcs2[['RNA']], Class = 'Assay')
	fds1[['RNA']] <- as(object = fds1[['RNA']], Class = 'Assay')
	fds2[['RNA']] <- as(object = fds2[['RNA']], Class = 'Assay')
	fds3[['RNA']] <- as(object = fds3[['RNA']], Class = 'Assay')
	fds4[['RNA']] <- as(object = fds4[['RNA']], Class = 'Assay')
	mds1[['RNA']] <- as(object = mds1[['RNA']], Class = 'Assay')
	mds2[['RNA']] <- as(object = mds2[['RNA']], Class = 'Assay')
	mds3[['RNA']] <- as(object = mds3[['RNA']], Class = 'Assay')
	
	wound.list <- list(fcs1 = fcs1, fcs2 = fcs2, fcs3 = fcs3, fcs4 = fcs4,
	                   mcs1 = mcs1, mcs2 = mcs2,
	                   fds1 = fds1, fds2 = fds2, fds3 = fds3, fds4 = fds4,
	                   mds1 = mds1, mds2 = mds2, mds3 = mds3)
	
	#File cleanup
	rm(fcs1, fcs2, fcs3, fcs4,
	   mcs1, mcs2,
	   fds1, fds2, fds3, fds4,
	   mds1, mds2, mds3)
	
	# Harmony integration
	
	wound.list <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskinlist.rds')
	
	#Merge Seurat objects to proceed with CellCycleScoring
	wound.merge <- merge(x = wound.list[[1]],
	                     y = wound.list[-1],
	                     add.cell.ids = names(wound.list),
	                     merge.data = F,
	                     project = 'wound_skin')
	
	wound.merge <- NormalizeData(wound.merge)
	
	wound.merge <- JoinLayers(wound.merge) # need to join layers as CellCycleScoring can't do separate layers (GetAssayData() doesn't work for multliple layers in v5 assay)
	
	m.s.genes <- cc.genes$s.genes %>% convert_human_to_mouse_symbols
	m.g2m.genes <- cc.genes$g2m.genes %>% convert_human_to_mouse_symbols
	
	wound.merge <- CellCycleScoring(wound.merge, 
	                                s.features = m.s.genes,
	                                g2m.features = m.g2m.genes, 
	                                set.ident = T)
	
	wound.merge[['CC.Diff']] <- wound.merge[['S.Score']] - wound.merge[['G2M.Score']]
	
	wound.merge[['RNA']] <- split(wound.merge[['RNA']], f = wound.merge$sampleid)
	
	wound.merge <- FindVariableFeatures(wound.merge, 
	                                    selection.method = 'vst', 
	                                    nfeatures = 5000)
	
	wound.merge <- SCTransform(wound.merge,
	                           method = 'glmGamPoi',
	                           vars.to.regress = 'CC.Diff',
	                           vst.flavor = 'v2',
	                           ncells = length(colnames(wound.merge)))
	
	wound.merge <- RunPCA(wound.merge, npcs = 100)
	
	ElbowPlot(wound.merge, ndims = 100) # choose 90 dimensions
	
	saveRDS(wound.merge,
	        '/data/overmilleram/scRNAseq/Wound/woundskin.merge.rds',
	        compress = F)
	
	wound.integrated <- IntegrateLayers(wound.merge, 
	                                    method = HarmonyIntegration, 
	                                    orig.reduction = 'pca',
	                                    new.reduction = 'harmony',
	                                    assay = 'SCT',
	                                    theta = 4, # variable that increasing diversity of clusters as increased (i.e., higher theta = more, unique clusters)
	                                    sigma = 0.1,
	                                    npcs = 90, 
	                                    verbose = T)
	
	wound.integrated <- FindNeighbors(wound.integrated, 
	                                  reduction = 'harmony', 
	                                  assay = 'SCT',
	                                  dims = 1:90)
	
	wound.integrated <- RunUMAP(wound.integrated, 
	                            reduction = 'harmony', 
	                            reduction.name = 'umap.harmony',
	                            assay = 'SCT',
	                            dims = 1:90)
	
	wound.integrated <- FindClusters(wound.integrated, 
	                                 method = 'igraph',
	                                 algorithm = 4,
	                                 resolution = 1)
	
	DimPlot(wound.integrated,
	        #split.by = 'sampleid',
	        #ncol = 3,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 0.75,
	        reduction = 'umap.harmony',
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Harmony UMAP.svg",
	        path = "/data/overmilleram/scRNAseq/Wound/",
	        width = 5000,
	        height =5000,
	        units = "px")
	
	wound.integrated <- JoinLayers(wound.integrated, assay = 'RNA')
	
	saveRDS(wound.integrated,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds",
	        compress = FALSE)
	
	wound.integrated.markers <- FindAllMarkers(wound.integrated,
	                                           only.pos = T, 
	                                           test.use = "MAST",
	                                           latent.vars = "sex",
	                                           min.pct = 0.50,
	                                           logfc.threshold = 2.00,
	                                           return.thresh = 0.05,
	                                           assay = "RNA",
	                                           densify = T) 
	
	wound.integrated.markers[ ,1] <- wound.integrated.markers[ ,7]
	
	wound.integrated.markers <- wound.integrated.markers[order(wound.integrated.markers$cluster,
	                                                         -wound.integrated.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.integrated.cellcounts <- table(wound.integrated@meta.data$wound.res.1,
	                                    wound.integrated@meta.data$state)
	
	wound.integrated.cellcounts2 <- table(wound.integrated@meta.data$wound.res.1,
	                                     wound.integrated@meta.data$sampleid)
	
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.integrated.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           row.names = F,
	           append = FALSE)
	write.xlsx(wound.integrated.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.integrated.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(wound.integrated, 
	        features = c('Krt5', 'Sbsn', 'Krt79', 'Krt4', 'Tgm3', 'Muc19', 'Bpifb2', 'Col1a1', 'Dcn', 'Lum', 'Ptprc', 
	                     'Cd3e', 'G0s2', 'Pecam1', 'Ccl21a', 'Myh11', 'Mpz', 'Mbp', 'Scn7a', 'Des', 'Dmd', 'Myl1', 'Tnnt3', 'Mylpf', 'Mlana', 'Pmel', 'Kit',
	                     'Hba-a1', 'Hbb-bs', 'percent.mt', 'percent.ribo', 'Lars2'), 
	        stack = T,
	        flip = T,
	        pt.size = 0,
	        assay = 'RNA') + NoLegend()
	
	wound.integrated <- subset(wound.integrated, idents = c('20', '36'), invert = T) # remove Lars2 junk cluster 20 and RBC cluster 36
	
	wound.integrated <- RenameIdents(wound.integrated,
	                                 "1" =  'Immune',
	                                 "2" =  'Immune',
	                                 "3" =  'Immune',
	                                 "4" =  'Keratinocyte',
	                                 "5" =  'Immune',
	                                 "6" =  'Keratinocyte',
	                                 "7" =  'Keratinocyte',
	                                 "8" =  'Fibroblast',
	                                 "9" =  'Immune',
	                                 "10" = 'Immune',
	                                 "11" = 'Keratinocyte',
	                                 "12" = 'Immune',
	                                 "13" = 'Keratinocyte',
	                                 "14" = 'Immune',
	                                 "15" = 'Immune',
	                                 "16" = 'Immune',
	                                 "17" = 'Keratinocyte',
	                                 "18" = 'Keratinocyte',
	                                 "19" = 'Keratinocyte',
	                                 "21" = 'Vascular',
	                                 "22" = 'Immune',
	                                 "23" = 'Immune',
	                                 "24" = 'Immune',
	                                 "25" = 'Immune',
	                                 "26" = 'Keratinocyte',
	                                 "27" = 'Keratinocyte',
	                                 "28" = 'Immune',
	                                 "29" = 'Mesenchymal',
	                                 '30' = 'Melanocyte',
	                                 '31' = 'Immune',
	                                 '32' = 'Fibroblast',
	                                 '33' = 'Immune',
	                                 '34' = 'Immune',
	                                 '35' = 'Immune')
	
	levels(wound.integrated) <- c('Keratinocyte', 
	                              'Immune',
	                              'Fibroblast',
	                              'Vascular',
	                              'Melanocyte',
	                              'Mesenchymal')
	
	wound.integrated$celltype <- wound.integrated@active.ident
	Idents(wound.integrated) <- "celltype"
	
	DimPlot(wound.integrated,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 1,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Cell Type UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 4000,
	        height =4000,
	        units = "px")
	
	DimPlot(wound.integrated,
	        label = F,
	        pt.size = 1,
	        label.size = 10,
	        split.by = 'state',
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Cell Type Split State UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 8000,
	        height = 4000,
	        units = "px")
	
	saveRDS(wound.integrated,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds",
	        compress = FALSE)
	
	Idents(wound.integrated) <- "celltype"
	
	wound.kera <- subset(wound.integrated,
	                     idents = "Keratinocyte")
	wound.imm <- subset(wound.integrated,
	                    idents = "Immune")
	wound.fib <- subset(wound.integrated,
	                    idents = "Fibroblast")
	wound.vasc <- subset(wound.integrated,
	                     idents = "Vascular")
	wound.mel <- subset(wound.integrated,
	                    idents = "Melanocyte")
	wound.mes <- subset(wound.integrated,
	                    idents = "Mesenchymal")
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	# Keratinocyte
	
	wound.kera <- SCTransform(wound.kera,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(wound.kera))) # workaround for 'none of the requested variables to regress...' error
	
	wound.kera <- RunPCA(wound.kera, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(wound.kera, ndims = 100) #decided to do 100 PCs for this analysis
	
	wound.kera <- FindNeighbors(wound.kera, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:100)
	
	wound.kera2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.kera, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'kera.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('kera.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.kera)
	
	wound.kera3 <- data.frame(as.numeric(as.character(wound.kera2[[1]]@meta.data$kera.res)))
	
	for (i in 1:10) {
	  wound.kera3[ ,i] <- data.frame(as.numeric(as.character(wound.kera2[[i]]@meta.data$kera.res)))
	}  
	
	colnames(wound.kera3) <- col.names
	rownames(wound.kera3) <- row.names
	
	wound.kera <- AddMetaData(wound.kera, 
	                         wound.kera3)
	
	rm(wound.kera2, wound.kera3)
	
	gc()
	
	wound.kera <- RunUMAP(wound.kera, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     assay = 'SCT',
	                     dims = 1:100)
	
	saveRDS(wound.kera,
	        file = "woundskin.kera.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.kera, prefix = 'kera.res.')
	ggsave2(filename = "Wound Skin Keratinocyte Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Immune
	
	wound.imm <- SCTransform(wound.imm,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(wound.imm))) # workaround for 'none of the requested variables to regress...' error
	
	wound.imm <- RunPCA(wound.imm, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(wound.imm, ndims = 100) #decided to do 100 PCs for this analysis
	
	wound.imm <- FindNeighbors(wound.imm, 
	                           assay = 'SCT',
	                           reduction = 'harmony', 
	                           dims = 1:100)
	
	wound.imm2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.imm, 
	                                                                 cluster.name = 'imm.res',
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.imm)
	
	wound.imm3 <- data.frame(as.numeric(as.character(wound.imm2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  wound.imm3[ ,i] <- data.frame(as.numeric(as.character(wound.imm2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(wound.imm3) <- col.names
	rownames(wound.imm3) <- row.names
	
	wound.imm <- AddMetaData(wound.imm, 
	                          wound.imm3)
	
	rm(wound.imm2, wound.imm3)
	
	gc()
	
	wound.imm <- RunUMAP(wound.imm, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     assay = 'SCT',
	                     dims = 1:100)
	
	saveRDS(wound.imm,
	        file = "woundskin.imm.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.imm, prefix = 'imm.res.')
	ggsave2(filename = "Wound Skin Immune Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Fibroblast
	
	wound.fib <- SCTransform(wound.fib,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(wound.fib))) # workaround for 'none of the requested variables to regress...' error
	
	wound.fib <- RunPCA(wound.fib, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(wound.fib, ndims = 100) #decided to do 60 PCs for this analysis
	
	wound.fib <- FindNeighbors(wound.fib, 
	                           assay = 'SCT',
	                           reduction = 'harmony', 
	                           dims = 1:70)
	
	wound.fib2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.fib, 
	                                                                 method = 'igraph',
	                                                                 cluster.name = 'fib.res',
	                                                                 algorithm = 4,
	                                                                 resolution = n)
	col.names <- paste0('fib.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.fib)
	
	wound.fib3 <- data.frame(as.numeric(as.character(wound.fib2[[1]]@meta.data$fib.res)))
	
	for (i in 1:10) {
	  wound.fib3[ ,i] <- data.frame(as.numeric(as.character(wound.fib2[[i]]@meta.data$fib.res)))
	}  
	
	colnames(wound.fib3) <- col.names
	rownames(wound.fib3) <- row.names
	
	wound.fib <- AddMetaData(wound.fib, 
	                         wound.fib3)
	
	rm(wound.fib2, wound.fib3)
	
	gc()
	
	wound.fib <- RunUMAP(wound.fib, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     assay = 'SCT',
	                     dims = 1:70)
	
	saveRDS(wound.fib,
	        file = "woundskin.fib.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.fib, prefix = 'fib.res.')
	ggsave2(filename = "Wound Skin Fibroblast Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Vascular
	
	wound.vasc <- SCTransform(wound.vasc,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(wound.vasc))) # workaround for 'none of the requested variables to regress...' error
	
	wound.vasc <- RunPCA(wound.vasc, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(wound.vasc, ndims = 100) #decided to do 80 PCs for this analysis
	
	wound.vasc <- FindNeighbors(wound.vasc, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:80)
	
	wound.vasc2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.vasc, 
	                                                                  cluster.name = 'vasc.res',
	                                                                  method = 'igraph',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('vasc.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.vasc)
	
	wound.vasc3 <- data.frame(as.numeric(as.character(wound.vasc2[[1]]@meta.data$vasc.res)))
	
	for (i in 1:10) {
	  wound.vasc3[ ,i] <- data.frame(as.numeric(as.character(wound.vasc2[[i]]@meta.data$vasc.res)))
	}  
	
	colnames(wound.vasc3) <- col.names
	rownames(wound.vasc3) <- row.names
	
	wound.vasc <- AddMetaData(wound.vasc, 
	                         wound.vasc3)
	
	rm(wound.vasc2, wound.vasc3)
	
	gc()
	
	wound.vasc <- RunUMAP(wound.vasc, 
	                      assay = 'SCT',
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      dims = 1:80)
	
	saveRDS(wound.vasc,
	        file = "woundskin.vasc.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.vasc, prefix = 'vasc.res.')
	ggsave2(filename = "Wound Skin Vascular Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Melanocyte
	
	wound.mel <- SCTransform(wound.mel,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(wound.mel))) # workaround for 'none of the requested variables to regress...' error
	
	wound.mel <- RunPCA(wound.mel, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(wound.mel, ndims = 100) #decided to do 80 PCs for this analysis
	
	wound.mel <- FindNeighbors(wound.mel, 
	                           assay = 'SCT',
	                           reduction = 'harmony', 
	                           dims = 1:80)
	
	wound.mel2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.mel, 
	                                                                 cluster.name = 'mel.res',
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 resolution = n)
	col.names <- paste0('mel.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.mel)
	
	wound.mel3 <- data.frame(as.numeric(as.character(wound.mel2[[1]]@meta.data$mel.res)))
	
	for (i in 1:10) {
	  wound.mel3[ ,i] <- data.frame(as.numeric(as.character(wound.mel2[[i]]@meta.data$mel.res)))
	}  
	
	colnames(wound.mel3) <- col.names
	rownames(wound.mel3) <- row.names
	
	wound.mel <- AddMetaData(wound.mel, 
	                         wound.mel3)
	
	rm(wound.mel2, wound.mel3)
	
	gc()
	
	wound.mel <- RunUMAP(wound.mel, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     assay = 'SCT',
	                     dims = 1:80)
	
	saveRDS(wound.mel,
	        file = "woundskin.mel.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.mel, prefix = 'mel.res.')
	ggsave2(filename = "Wound Skin Melanocyte Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Mesenchymal
	
	wound.mes <- SCTransform(wound.mes,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(wound.mes))) # workaround for 'none of the requested variables to regress...' error
	
	wound.mes <- RunPCA(wound.mes, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(wound.mes, ndims = 100) #decided to do 100 PCs for this analysis
	
	wound.mes <- FindNeighbors(wound.mes, 
	                           assay = 'SCT',
	                           reduction = 'harmony', 
	                           dims = 1:100)
	
	wound.mes2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.mes, 
	                                                                 cluster.name = 'mes.res',
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 resolution = n)
	col.names <- paste0('mes.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.mes)
	
	wound.mes3 <- data.frame(as.numeric(as.character(wound.mes2[[1]]@meta.data$mes.res)))
	
	for (i in 1:10) {
	  wound.mes3[ ,i] <- data.frame(as.numeric(as.character(wound.mes2[[i]]@meta.data$mes.res)))
	}  
	
	colnames(wound.mes3) <- col.names
	rownames(wound.mes3) <- row.names
	
	wound.mes <- AddMetaData(wound.mes, 
	                         wound.mes3)
	
	rm(wound.mes2, wound.mes3)
	
	gc()
	
	wound.mes <- RunUMAP(wound.mes, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     assay = 'SCT',
	                     dims = 1:100)
	
	saveRDS(wound.mes,
	        file = "woundskin.mes.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.mes, prefix = 'mes.res.')
	ggsave2(filename = "Wound Skin Mesenchymal Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	wound.kera <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.kera.rds")
	wound.fib <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.fib.rds")
	wound.imm <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds")
	wound.vasc <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.vasc.rds")
	wound.mel <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.mel.rds")
	wound.mes <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.mes.rds")
	
	Idents(wound.kera) = 'kera.res.0.8'
	Idents(wound.imm) =  'imm.res.0.2' # will split into myeloid & lymphoid and recluster
	Idents(wound.fib) =  'fib.res.0.5' 
	Idents(wound.vasc) = 'vasc.res.0.3' 
	Idents(wound.mel) = 'mel.res.0.1'
	Idents(wound.mes) =  'mes.res.0.3' 
	
	DimPlot(wound.kera,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Keratinocyte UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height =5000,
	        units = "px")
	
	
	DimPlot(wound.imm,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Immune UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height =5000,
	        units = "px")
	
	
	DimPlot(wound.fib,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 4,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Fibroblast UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height =5000,
	        units = "px")
	
	
	DimPlot(wound.vasc,
	         label = TRUE,
	         repel = TRUE,
	         pt.size = 7,
	         label.size = 10,
	         raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Vascular UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height =5000,
	        units = "px")
	
	DimPlot(wound.mel,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 8,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Melanocyte UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height =5000,
	        units = "px")
	
	DimPlot(wound.mes,
	         label = TRUE,
	         repel = TRUE,
	         pt.size = 8,
	         label.size = 10,
	         raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Mesenchymal UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height =5000,
	        units = "px")
	
	Idents(wound.kera) = 'kera.res.0.8'
	Idents(wound.imm) =  'imm.res.0.2' # will split into myeloid & lymphoid and recluster
	Idents(wound.fib) =  'fib.res.0.5' 
	Idents(wound.vasc) = 'vasc.res.0.3' 
	Idents(wound.mel) = 'mel.res.0.1'
	Idents(wound.mes) =  'mes.res.0.3'
	
	wound.kera.harmony.markers <- FindAllMarkers(wound.kera,
	                                          only.pos = T, 
	                                          test.use = "MAST",
	                                          latent.vars = "sex",
	                                          min.pct = 0.33,
	                                          logfc.threshold = 1.00,
	                                          return.thresh = 0.05,
	                                          assay = "RNA",
	                                          densify = T) 
	
	wound.kera.harmony.markers <- wound.kera.harmony.markers[order(wound.kera.harmony.markers$cluster,
	                                                         -wound.kera.harmony.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.kera.harmony.cellcounts <- table(wound.kera@meta.data$kera.res.0.8,
	                                    wound.kera@meta.data$state)
	
	wound.kera.harmony.cellcounts2 <- table(wound.kera@meta.data$kera.res.0.8,
	                                     wound.kera@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.kera.harmony.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Keratinocyte Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.kera.harmony.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Keratinocyte Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.kera.harmony.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Keratinocyte Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.imm.harmony.markers <- FindAllMarkers(wound.imm,
	                                         only.pos = T, 
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         min.pct = 0.33,
	                                         logfc.threshold = 1.00,
	                                         return.thresh = 0.05,
	                                         assay = "RNA",
	                                         densify = T) 
	
	wound.imm.harmony.markers <- wound.imm.harmony.markers[order(wound.imm.harmony.markers$cluster,
	                                                         -wound.imm.harmony.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.imm.harmony.cellcounts <- table(wound.imm@meta.data$imm.res.0.2,
	                                    wound.imm@meta.data$state)
	
	wound.imm.harmony.cellcounts2 <- table(wound.imm@meta.data$imm.res.0.2,
	                                     wound.imm@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.imm.harmony.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Immune Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.imm.harmony.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Immune Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.imm.harmony.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Immune Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.fib.harmony.markers <- FindAllMarkers(wound.fib,
	                                          only.pos = T, 
	                                          test.use = "MAST",
	                                          latent.vars = "sex",
	                                          min.pct = 0.33,
	                                          logfc.threshold = 1.00,
	                                          return.thresh = 0.05,
	                                          assay = "RNA",
	                                          densify = T) 
	
	wound.fib.harmony.markers <- wound.fib.harmony.markers[order(wound.fib.harmony.markers$cluster,
	                                                         -wound.fib.harmony.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.fib.harmony.cellcounts <- table(wound.fib@meta.data$fib.res.0.5,
	                                    wound.fib@meta.data$state)
	
	wound.fib.harmony.cellcounts2 <- table(wound.fib@meta.data$fib.res.0.5,
	                                     wound.fib@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.fib.harmony.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Fibroblast Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.fib.harmony.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Fibroblast Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.fib.harmony.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Fibroblast Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.vasc.harmony.markers <- FindAllMarkers(wound.vasc,
	                                          only.pos = T, 
	                                          test.use = "MAST",
	                                          latent.vars = "sex",
	                                          min.pct = 0.33,
	                                          logfc.threshold = 1.00,
	                                          return.thresh = 0.05,
	                                          assay = "RNA",
	                                          densify = T) 
	
	wound.vasc.harmony.markers <- wound.vasc.harmony.markers[order(wound.vasc.harmony.markers$cluster,
	                                                         -wound.vasc.harmony.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.vasc.harmony.cellcounts <- table(wound.vasc@meta.data$vasc.res.0.3,
	                                    wound.vasc@meta.data$state)
	
	wound.vasc.harmony.cellcounts2 <- table(wound.vasc@meta.data$vasc.res.0.3,
	                                     wound.vasc@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.vasc.harmony.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Vascular Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.vasc.harmony.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Vascular Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.vasc.harmony.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Vascular Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.mel.harmony.markers <- FindAllMarkers(wound.mel,
	                                            only.pos = T, 
	                                            test.use = "MAST",
	                                            latent.vars = "sex",
	                                            min.pct = 0.33,
	                                            logfc.threshold = 1.5,
	                                            return.thresh = 0.05,
	                                            assay = "RNA",
	                                            densify = T) 
	
	wound.mel.harmony.markers <- wound.mel.harmony.markers[order(wound.mel.harmony.markers$cluster,
	                                                             -wound.mel.harmony.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.mel.harmony.cellcounts <- table(wound.mel@meta.data$mel.res.0.1,
	                                      wound.mel@meta.data$state)
	
	wound.mel.harmony.cellcounts2 <- table(wound.mel@meta.data$mel.res.0.1,
	                                       wound.mel@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.mel.harmony.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Melanocyte Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.mel.harmony.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Melanocyte Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.mel.harmony.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Melanocyte Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.mes.harmony.markers <- FindAllMarkers(wound.mes,
	                                          only.pos = T, 
	                                          test.use = "MAST",
	                                          latent.vars = "sex",
	                                          min.pct = 0.33,
	                                          logfc.threshold = 0.50,
	                                          return.thresh = 0.05,
	                                          assay = "RNA",
	                                          densify = T) 
	
	wound.mes.harmony.markers <- wound.mes.harmony.markers[order(wound.mes.harmony.markers$cluster,
	                                                         -wound.mes.harmony.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.mes.harmony.cellcounts <- table(wound.mes@meta.data$mes.res.0.3,
	                                    wound.mes@meta.data$state)
	
	wound.mes.harmony.cellcounts2 <- table(wound.mes@meta.data$mes.res.0.3,
	                                     wound.mes@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.mes.harmony.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Mesenchymal Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.mes.harmony.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Mesenchymal Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.mes.harmony.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Mesenchymal Harmony Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	# Annotate Keratinocytes
	###########
	
	wound.kera <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.kera.rds")
	
	Idents(wound.kera) <- "kera.res.0.8"
	
	VlnPlot(wound.kera,
	features = c('Krt14', 'Col17a1', 'Trp63', 'Sbsn', 'Krt10', 'Tgm1', 'Skint5', 'Krt6a', 'Krt16', 'Krt17'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	c6.9 <- FindMarkers(wound.kera,
	                    ident.1 = c('6','9'),
	                    assay = 'RNA',
	                    only.pos = T,
	                    test.use = 'MAST',
	                    latent.vars = 'sex',
	                    logfc.threshold = 0.75,
	                    min.pct = 0.33,
	                    densify = T)
	
	wound.kera <- RenameIdents(wound.kera,
	                           '1' =  'Epithelial_Basal',
	                           '2' =  'Wound_Activated',
	                           '3' =  'HFSC_1',
	                           '4' =  'Epithelial_Suprabasal',
	                           '5' =  'Upper_HF_Suprabasal',
	                           '6' =  'Lower_HF_1',
	                           '7' =  'Wound_Migratory',
	                           '8' =  'Lower_HF_2',
	                           '9' =  'Lower_HF_1',
	                           '10' = 'HFSC_1',
	                           '11' = 'HFSC_2',
	                           '12' = 'HFSC_1',
	                           '13' = 'Proliferating_Keratinocyte',
	                           '14' = 'HFSC_2',
	                           '15' = 'Sebaceous',
	                           '16' = 'Proliferating_Keratinocyte')
	
	levels(wound.kera) <- c('Epithelial_Basal',
	                        'Epithelial_Suprabasal',
	                        'Proliferating_Keratinocyte',
	                        'Upper_HF_Suprabasal',
	                        'Sebaceous',
	                        'HFSC_1',
	                        'HFSC_2',
	                        'Lower_HF_1',
	                        'Lower_HF_2',
	                        'Wound_Activated',
	                        'Wound_Migratory')
	
	wound.kera$L1subtype = wound.kera@active.ident
	
	Idents(wound.kera) <- "kera.res.0.8"
	
	wound.kera <- RenameIdents(wound.kera,
	                           '1' =  'Epithelial_Basal',
	                           '2' =  'Wound_Activated',
	                           '3' =  'HFSC_1a',
	                           '4' =  'Epithelial_Suprabasal',
	                           '5' =  'Upper_HF_Suprabasal',
	                           '6' =  'Lower_HF_1',
	                           '7' =  'Wound_Migratory',
	                           '8' =  'Lower_HF_2',
	                           '9' =  'Lower_HF_1',
	                           '10' = 'HFSC_1b',
	                           '11' = 'HFSC_2a',
	                           '12' = 'HFSC_1c',
	                           '13' = 'Proliferating_Keratinocyte_2',
	                           '14' = 'HFSC_2b',
	                           '15' = 'Sebaceous',
	                           '16' = 'Proliferating_Keratinocyte_1')
	
	levels(wound.kera) <- c('Epithelial_Basal',
	                        'Epithelial_Suprabasal',
	                        'Proliferating_Keratinocyte_1',
	                        'Proliferating_Keratinocyte_2',
	                        'Upper_HF_Suprabasal',
	                        'Sebaceous',
	                        'HFSC_1a',
	                        'HFSC_1b',
	                        'HFSC_1c',
	                        'HFSC_2a',
	                        'HFSC_2b',
	                        'Lower_HF_1',
	                        'Lower_HF_2',
	                        'Wound_Activated',
	                        'Wound_Migratory')
	
	wound.kera$L2subtype = wound.kera@active.ident
	
	saveRDS(wound.kera,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.kera.rds",
	        compress = FALSE)
	
	Idents(wound.kera) = 'L1subtype'
	
	DimPlot(wound.kera,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 4,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Keratinocyte L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	###########
	
	# Subcluster and annotate immune cells
	###########
	
	Idents(wound.imm) =  'imm.res.0.2'
	
	VlnPlot(wound.imm, 
	        features = c('C1qa', 'Folr2', 'Ccl8', 'H2-Eb1', 'Plbd1', 'Cd74', 'Cd207', 'Mfge8', 'Fn1', 'Plac8', 'Thbs1'), 
	        stack = T,
	        flip = T,
	        pt.size = 0,
	        assay = 'RNA') + NoLegend()
	
	# split Immune cells into Myeloid, Granulocytes, & Lymphoid
	
	wound.imm <- RenameIdents(wound.imm,
	                           '1' = 'Myeloid',
	                           '2' = 'Myeloid',
	                           '3' = 'Granulocyte',
	                           '4' = 'Granulocyte',
	                           '5' = 'Myeloid',
	                           '6' = 'Lymphoid',
	                           '7' = 'Myeloid',
	                           '8' = 'Myeloid',
	                           '9' = 'Myeloid',
	                           '10' = 'Granulocyte',
	                           '11' = 'Myeloid',
	                           '12' = 'Myeloid',
	                           '13' = 'Myeloid')
	
	levels(wound.imm) <- c("Myeloid",
	                       'Granulocyte',
	                       "Lymphoid")
	
	wound.imm$lineage <- wound.imm@active.ident
	
	saveRDS(wound.imm,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds",
	        compress = FALSE)
	
	wound.lym <- subset(wound.imm, idents = 'Lymphoid') #t-cells/nk-cells
	wound.gra <- subset(wound.imm, idents = 'Granulocyte') #neutrophils/basophils
	wound.mye <- subset(wound.imm, idents = 'Myeloid') #macrophage/dendritic/langerhans/monocyte
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	# Lymphoid
	
	wound.lym <- SCTransform(wound.lym,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(wound.lym))) # workaround for 'none of the requested variables to regress...' error
	
	wound.lym <- RunPCA(wound.lym, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(wound.lym, ndims = 100) #decided to do 70 PCs for this analysis
	
	wound.lym <- FindNeighbors(wound.lym, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:70)
	
	wound.lym@meta.data[ ,24:33] <- NULL # strip old imm.res clusters
	
	wound.lym2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.lym, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.lym)
	
	wound.lym3 <- data.frame(as.numeric(as.character(wound.lym2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  wound.lym3[ ,i] <- data.frame(as.numeric(as.character(wound.lym2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(wound.lym3) <- col.names
	rownames(wound.lym3) <- row.names
	
	wound.lym <- AddMetaData(wound.lym, 
	                          wound.lym3)
	
	rm(wound.lym2, wound.lym3)
	
	wound.lym <- RunUMAP(wound.lym, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	saveRDS(wound.lym,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.lym.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.lym, prefix = 'imm.res.')
	ggsave2(filename = "Wound Skin Immune Lymphoid Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Granulocytes
	
	wound.gra <- SCTransform(wound.gra,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(wound.gra))) # workaround for 'none of the requested variables to regress...' error
	
	wound.gra <- RunPCA(wound.gra, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(wound.gra, ndims = 100) #decided to do 70 PCs for this analysis
	
	wound.gra <- FindNeighbors(wound.gra, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:70)
	
	wound.gra@meta.data[ ,24:33] <- NULL # strip old clusters
	
	wound.gra2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.gra, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.gra)
	
	wound.gra3 <- data.frame(as.numeric(as.character(wound.gra2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  wound.gra3[ ,i] <- data.frame(as.numeric(as.character(wound.gra2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(wound.gra3) <- col.names
	rownames(wound.gra3) <- row.names
	
	wound.gra <- AddMetaData(wound.gra, 
	                          wound.gra3)
	
	rm(wound.gra2, wound.gra3)
	
	wound.gra <- RunUMAP(wound.gra, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	saveRDS(wound.gra,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.gra.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.gra, prefix = 'imm.res.')
	ggsave2(filename = "Wound Skin Immune Granulocyte Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	# Myeloid (Macrophage/Dendritic/Langerhans/Monocyte)
	
	wound.mye <- SCTransform(wound.mye,
	                          method = 'glmGamPoi',
	                          vars.to.regress = 'CC.Diff',
	                          vst.flavor = 'v2',
	                          do.scale = T,
	                          ncells = length(colnames(wound.mye))) # workaround for 'none of the requested variables to regress...' error
	
	wound.mye <- RunPCA(wound.mye, 
	                     npcs = 100, 
	                     assay = 'SCT')
	
	ElbowPlot(wound.mye, ndims = 100) #decided to do 90 PCs for this analysis
	
	wound.mye <- FindNeighbors(wound.mye, 
	                            assay = 'SCT',
	                            reduction = 'harmony', 
	                            dims = 1:90)
	
	wound.mye@meta.data[ ,24:33] <- NULL # strip old clusters
	
	wound.mye2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.mye, 
	                                                                  method = 'igraph',
	                                                                  cluster.name = 'imm.res',
	                                                                  algorithm = 4,
	                                                                  resolution = n)
	col.names <- paste0('imm.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.mye)
	
	wound.mye3 <- data.frame(as.numeric(as.character(wound.mye2[[1]]@meta.data$imm.res)))
	
	for (i in 1:10) {
	  wound.mye3[ ,i] <- data.frame(as.numeric(as.character(wound.mye2[[i]]@meta.data$imm.res)))
	}  
	
	colnames(wound.mye3) <- col.names
	rownames(wound.mye3) <- row.names
	
	wound.mye <- AddMetaData(wound.mye, 
	                          wound.mye3)
	
	rm(wound.mye2, wound.mye3)
	
	wound.mye <- RunUMAP(wound.mye, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:90)
	
	saveRDS(wound.mye,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.mye.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.mye, prefix = 'imm.res.')
	ggsave2(filename = "Wound Skin Immune Myeloid Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(wound.lym) = 'imm.res.1'
	Idents(wound.gra) = 'imm.res.0.2'
	Idents(wound.mye) = 'imm.res.0.8'
	
	DimPlot(wound.lym,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 7,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Immune Lymphoid UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	DimPlot(wound.gra,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 3,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Immune Granulocyte UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	DimPlot(wound.mye,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 2,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Immune Myeloid UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	wound.lym.markers <- FindAllMarkers(wound.lym,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.33,
	                                     logfc.threshold = 0.50,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	wound.lym.markers <- wound.lym.markers[order(wound.lym.markers$cluster,
	                                               -wound.lym.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.lym.cellcounts <- table(wound.lym@meta.data$imm.res.1,
	                               wound.lym@meta.data$state)
	
	wound.lym.cellcounts2 <- table(wound.lym@meta.data$imm.res.1,
	                                wound.lym@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.lym.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Lymphoid Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.lym.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Lymphoid Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.lym.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Lymphoid Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.gra.markers <- FindAllMarkers(wound.gra,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.33,
	                                     logfc.threshold = 0.50,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	wound.gra.markers <- wound.gra.markers[order(wound.gra.markers$cluster,
	                                               -wound.gra.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.gra.cellcounts <- table(wound.gra@meta.data$imm.res.0.2,
	                               wound.gra@meta.data$state)
	
	wound.gra.cellcounts2 <- table(wound.gra@meta.data$imm.res.0.2,
	                                wound.gra@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.gra.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Granulocyte Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.gra.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Granulocyte Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.gra.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Granulocyte Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.mye.markers <- FindAllMarkers(wound.mye,
	                                     only.pos = T, 
	                                     test.use = "MAST",
	                                     latent.vars = "sex",
	                                     min.pct = 0.33,
	                                     logfc.threshold = 0.50,
	                                     return.thresh = 0.05,
	                                     assay = "RNA",
	                                     densify = T) 
	
	wound.mye.markers <- wound.mye.markers[order(wound.mye.markers$cluster,
	                                               -wound.mye.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.mye.cellcounts <- table(wound.mye@meta.data$imm.res.0.8,
	                               wound.mye@meta.data$state)
	
	wound.mye.cellcounts2 <- table(wound.mye@meta.data$imm.res.0.8,
	                                wound.mye@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.mye.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Myeloid Harmony Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.mye.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Myeloid Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.mye.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Immune Myeloid Harmony Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	VlnPlot(wound.mye,
	        features = c('Lyst', 'Retreg1', 'Lrmda', 'Irak2', 'Lpl', 'Mmp12', 'Cd36', 'Rnf128', 'Ctsk', 'Itgax'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	FeaturePlot(wound.mye,
	            features = 'Foxp3',
	            pt.size = 3,
	            label = F,
	            order = T,
	            cols = c('grey', 'red'))
	
	c.3 <- FindMarkers(wound.mye,
	                   ident.1 = c('3','14'),
	                   assay = 'RNA',
	                   test.use = 'MAST',
	                   latent.vars = 'sex',
	                   only.pos = T,
	                   min.pct = 0.33,
	                   logfc.threshold = 0.75)
	
	Idents(wound.mye) <- 'imm.res.1'
	
	c.3a <- FindMarkers(wound.mye,
	                    ident.1 = '6',
	                    assay = 'RNA',
	                    test.use = 'MAST',
	                    latent.vars = 'sex',
	                    only.pos = T,
	                    min.pct = 0.33,
	                    logfc.threshold = 0.75)
	
	c.3b <- FindMarkers(wound.mye,
	                    ident.1 = '15',
	                    assay = 'RNA',
	                    test.use = 'MAST',
	                    latent.vars = 'sex',
	                    only.pos = T,
	                    min.pct = 0.33,
	                    logfc.threshold = 0.75)
	  
	# add metadata to immune cell subset; subset cells and add metadata to wound.imm
	
	wound.imm <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds')
	wound.lym <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.lym.rds')
	wound.gra <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.gra.rds')
	wound.mye <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.mye.rds')
	
	Idents(wound.lym) = 'imm.res.1'
	Idents(wound.gra) = 'imm.res.0.2'
	Idents(wound.mye) = 'imm.res.0.8'
	
	# Lymphoid
	
	Idents(wound.lym) = 'imm.res.1'
	
	# update lineage
	wound.lym <- RenameIdents(wound.lym,
	                          '1' = 'Lymphoid',
	                          '2' = 'Lymphoid',
	                          '3' = 'Lymphoid',
	                          '4' = 'Lymphoid',
	                          '5' = 'Lymphoid',
	                          '6' = 'Lymphoid',
	                          '7' = 'Lymphoid',
	                          '8' = 'Myeloid',
	                          '9' = 'Lymphoid',
	                          '10' = 'Lymphoid',
	                          '11' = 'Lymphoid')
	
	levels(wound.lym) <- c('Lymphoid', 
	                       'Myeloid')
	
	wound.lym$lineage = wound.lym@active.ident
	
	Idents(wound.lym) = 'imm.res.1'
	
	wound.lym <- RenameIdents(wound.lym,
	                          '1' = 'NK_cell',
	                          '2' = 'T_cell',
	                          '3' = 'T_cell',
	                          '4' = 'NK_cell',
	                          '5' = 'T_cell',
	                          '6' = 'T_cell',
	                          '7' = 'T_cell',
	                          '8' = 'pDC',
	                          '9' = 'T_cell',
	                          '10' = 'T_cell',
	                          '11' = 'NK_cell')
	
	levels(wound.lym) <- c('T_cell', 
	                       'NK_cell',
	                       'pDC')
	
	wound.lym$L1subtype = wound.lym@active.ident
	
	Idents(wound.lym) = 'imm.res.1'
	
	wound.lym <- RenameIdents(wound.lym,
	                          '1' = 'NK_cell',
	                          '2' = 'T_cell',
	                          '3' = 'T_cell',
	                          '4' = 'NK_cell',
	                          '5' = 'T_cell',
	                          '6' = 'T_cell',
	                          '7' = 'T_cell',
	                          '8' = 'pDC',
	                          '9' = 'T_cell',
	                          '10' = 'T_cell',
	                          '11' = 'NK_cell')
	
	levels(wound.lym) <- c('T_cell', 
	                       'NK_cell',
	                       'pDC')
	
	wound.lym$L2subtype = wound.lym@active.ident
	
	saveRDS(wound.lym,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.lym.rds",
	        compress = FALSE)
	
	# Macrophage, Dendritic, Langerhans, Monocyte
	
	Idents(wound.mye) = 'imm.res.0.8'
	
	wound.mye <- RenameIdents(wound.mye,
	                          '1' = 'Dendritic_1',
	                          '2' = 'Macrophage_1',
	                          '3' = 'Macrophage_3',
	                          '4' = 'Dendritic_1',
	                          '5' = 'Macrophage_1',
	                          '6' = 'Macrophage_2',
	                          '7' = 'Monocyte',
	                          '8' = 'Macrophage_2',
	                          '9' = 'Macrophage_3',
	                          '10' = 'Macrophage_3',
	                          '11' = 'Monocyte',
	                          '12' = 'Macrophage_3',
	                          '13' = 'Langerhans',
	                          '14' = 'Macrophage_3',
	                          '15' = 'Dendritic_2',
	                          '16' = 'Macrophage_3',
	                          '17' = 'Dendritic_3')
	
	levels(wound.mye) <- c('Macrophage_1', 
	                       'Macrophage_2',
	                       'Macrophage_3',
	                       'Dendritic_1',
	                       'Dendritic_2',
	                       'Dendritic_3',
	                       'Langerhans',
	                       'Monocyte')
	
	wound.mye$L1subtype = wound.mye@active.ident
	
	Idents(wound.mye) = 'imm.res.0.8'
	
	wound.mye <- RenameIdents(wound.mye,
	                          '1' = 'Dendritic_1a',
	                          '2' = 'Macrophage_1a',
	                          '3' = 'Macrophage_3a',
	                          '4' = 'Dendritic_1b',
	                          '5' = 'Macrophage_1b',
	                          '6' = 'Macrophage_2',
	                          '7' = 'Monocyte_1',
	                          '8' = 'Macrophage_2',
	                          '9' = 'Macrophage_3b',
	                          '10' = 'Macrophage_3c',
	                          '11' = 'Monocyte_2',
	                          '12' = 'Macrophage_3d',
	                          '13' = 'Langerhans',
	                          '14' = 'Macrophage_3e',
	                          '15' = 'Dendritic_2',
	                          '16' = 'Macrophage_3f',
	                          '17' = 'Dendritic_3')
	
	levels(wound.mye) <- c('Macrophage_1a',
	                       'Macrophage_1b',
	                       'Macrophage_2',
	                       'Macrophage_3a',
	                       'Macrophage_3b',
	                       'Macrophage_3c',
	                       'Macrophage_3d',
	                       'Macrophage_3e',
	                       'Macrophage_3f',
	                       'Dendritic_1a',
	                       'Dendritic_1b',
	                       'Dendritic_2',
	                       'Dendritic_3',
	                       'Langerhans',
	                       'Monocyte_1',
	                       'Monocyte_2')
	
	wound.mye$L2subtype = wound.mye@active.ident
	
	saveRDS(wound.mye,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.mye.rds",
	        compress = FALSE)
	
	# Granulocytes
	
	Idents(wound.gra) = 'imm.res.0.2'
	
	wound.gra <- RenameIdents(wound.gra,
	                          '1' = 'Neutrophil_1',
	                          '2' = 'Neutrophil_1',
	                          '3' = 'Neutrophil_1',
	                          '4' = 'Neutrophil_2')
	
	levels(wound.gra) <- c('Neutrophil_1', 
	                       'Neutrophil_2')
	
	wound.gra$L1subtype = wound.gra@active.ident
	
	Idents(wound.gra) = 'imm.res.0.2'
	
	wound.gra <- RenameIdents(wound.gra,
	                          '1' = 'Neutrophil_1',
	                          '2' = 'Neutrophil_1',
	                          '3' = 'Neutrophil_1',
	                          '4' = 'Neutrophil_2')
	
	levels(wound.gra) <- c('Neutrophil_1', 
	                       'Neutrophil_2')
	
	wound.gra$L2subtype = wound.gra@active.ident
	
	saveRDS(wound.gra,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.gra.rds",
	        compress = FALSE)
	
	# Transfer metadata to wound.imm
	
	lym.cells <- wound.lym$lineage %>% as.data.frame %>% rownames_to_column('barcode')
	mye.cells <- wound.mye$lineage %>% as.data.frame %>% rownames_to_column('barcode')
	gra.cells <- wound.gra$lineage %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(lym.cells, mye.cells, gra.cells)
	
	lineage <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(lineage) <- lineage[ ,1]
	lineage[ ,1] <- NULL
	colnames(lineage) <- 'lineage'
	
	lym.cells <- wound.lym$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mye.cells <- wound.mye$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	gra.cells <- wound.gra$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(lym.cells, mye.cells, gra.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	lym.cells <- wound.lym$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mye.cells <- wound.mye$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	gra.cells <- wound.gra$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(lym.cells, mye.cells, gra.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	wound.imm <- AddMetaData(wound.imm, lineage, col.name = 'lineage')
	wound.imm <- AddMetaData(wound.imm, L1subtype, col.name = 'L1subtype')
	wound.imm <- AddMetaData(wound.imm, L2subtype, col.name = 'L2subtype')
	
	saveRDS(wound.imm,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds",
	        compress = FALSE)
	
	Idents(wound.imm) = 'L1subtype'
	
	DimPlot(wound.imm,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 2,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Immune L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	###########
	
	# Annotate Fibroblasts
	###########
	
	wound.fib <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.fib.rds")
	
	Idents(wound.fib) =  'fib.res.0.5' 
	
	VlnPlot(wound.fib,
	        features = c('Col1a1', 'Col1a2', 'Dcn', 'Lum', 'Ptprc'),
	        flip = T,
	        stack = T,
	        assay = 'RNA')
	
	wound.fib <- RenameIdents(wound.fib,
	                          '1' =  'Dermal_Fibroblast_1',
	                          '2' =  'Myofibroblast_1',
	                          '3' =  'Myofibroblast_1',
	                          '4' =  'Dermal_Sheath',
	                          '5' =  'Dermal_Papilla',
	                          '6' =  'Dermal_Fibroblast_2',
	                          '7' =  'Myofibroblast_2',
	                          '8' =  'Dermal_Fibroblast_1')
	
	levels(wound.fib) <- c('Dermal_Fibroblast_1',
	                        'Dermal_Fibroblast_2',
	                        'Dermal_Sheath',
	                        'Dermal_Papilla',
	                        'Myofibroblast_1',
	                        'Myofibroblast_2')
	
	wound.fib$L1subtype = wound.fib@active.ident
	
	Idents(wound.fib) <- "fib.res.0.5"
	
	wound.fib <- RenameIdents(wound.fib,
	                          '1' =  'Dermal_Fibroblast_1a',
	                          '2' =  'Myofibroblast_1a',
	                          '3' =  'Myofibroblast_1b',
	                          '4' =  'Dermal_Sheath',
	                          '5' =  'Dermal_Papilla',
	                          '6' =  'Dermal_Fibroblast_2',
	                          '7' =  'Myofibroblast_2',
	                          '8' =  'Dermal_Fibroblast_1b')
	
	levels(wound.fib) <- c('Dermal_Fibroblast_1a',
	                       'Dermal_Fibroblast_1b',
	                       'Dermal_Fibroblast_2',
	                       'Dermal_Sheath',
	                       'Dermal_Papilla',
	                       'Myofibroblast_1a',
	                       'Myofibroblast_1b',
	                       'Myofibroblast_2')
	
	wound.fib$L2subtype = wound.fib@active.ident
	
	saveRDS(wound.fib,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.fib.rds",
	        compress = FALSE)
	
	Idents(wound.fib) = 'L1subtype'
	
	DimPlot(wound.fib,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 2,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Fibroblast L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	
	###########
	
	# Annotate Vascular cells
	###########
	
	wound.vasc <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.vasc.rds")
	
	Idents(wound.vasc) = 'vasc.res.0.3'
	
	wound.vasc <- RenameIdents(wound.vasc,
	                            '1' = 'Endothelial',
	                            '2' = 'Endothelial',
	                            '3' = 'Lymph_Vessel')
	
	levels(wound.vasc) <- c('Endothelial', 
	                        'Lymph_Vessel')
	
	wound.vasc$L1subtype = wound.vasc@active.ident
	
	Idents(wound.vasc) = 'vasc.res.0.3'
	
	wound.vasc <- RenameIdents(wound.vasc,
	                            '1' = 'Endothelial_1',
	                            '2' = 'Endothelial_2',
	                            '3' = 'Lymph_Vessel')
	
	levels(wound.vasc) <- c('Endothelial_1',
	                         'Endothelial_2',
	                         'Lymph_Vessel')
	
	wound.vasc$L2subtype = wound.vasc@active.ident
	
	saveRDS(wound.vasc,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.vasc.rds",
	        compress = FALSE)
	
	Idents(wound.vasc) = 'L1subtype'
	
	DimPlot(wound.vasc,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 4,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Vascular L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	
	###########
	
	# Annotate Melanocytes
	###########
	
	wound.mel <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.mel.rds")
	
	Idents(wound.mel) = 'mel.res.0.1'
	
	wound.mel <- RenameIdents(wound.mel,
	                          '1' = 'Skeletal_Muscle',
	                          '2' = 'Melanocyte')
	
	levels(wound.mel) <- c('Skeletal_Muscle', 
	                       'Melanocyte')
	
	wound.mel$celltype = wound.mel@active.ident
	
	Idents(wound.mel) = 'mel.res.0.1'
	
	wound.mel <- RenameIdents(wound.mel,
	                            '1' = 'Satellite_Cells',
	                            '2' = 'Melanocyte')
	
	levels(wound.mel) <- c('Satellite_Cells', 
	                       'Melanocyte')
	
	wound.mel$L1subtype = wound.mel@active.ident
	
	Idents(wound.mel) = 'mel.res.0.1'
	
	wound.mel <- RenameIdents(wound.mel,
	                          '1' = 'Satellite_Cells',
	                          '2' = 'Melanocyte')
	
	levels(wound.mel) <- c('Satellite_Cells', 
	                       'Melanocyte')
	
	wound.mel$L2subtype = wound.mel@active.ident
	
	saveRDS(wound.mel,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.mel.rds",
	        compress = FALSE)
	
	DimPlot(wound.mel,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 3,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Melanocyte L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	###########
	
	# Annotate Mesenchymal cells
	###########
	
	# Update Mesenchymal
	
	Idents(wound.mes) <-  "mes.res.0.3"
	
	wound.mes <- subset(wound.mes, idents = c('2'), invert = T)
	
	library(foreach)
	library(doParallel)
	library(doMC)
	library(rbenchmark)
	
	registerDoMC(cores = future::availableCores())
	
	wound.mes <- SCTransform(wound.mes,
	                         method = 'glmGamPoi',
	                         vars.to.regress = 'CC.Diff',
	                         vst.flavor = 'v2',
	                         do.scale = T,
	                         ncells = length(colnames(wound.mes))) # workaround for 'none of hte requested variables to regress...' error
	
	wound.mes <- RunPCA(wound.mes, 
	                    npcs = 100, 
	                    assay = 'SCT')
	
	ElbowPlot(wound.mes, ndims = 100) #decided to do 90 PCs for this analysis
	
	wound.mes <- FindNeighbors(wound.mes, 
	                           reduction = 'harmony', 
	                           dims = 1:90,
	                           assay = 'SCT')
	
	wound.mes2 <- foreach(n = seq(0.1,1.0,0.1)) %dopar% FindClusters(wound.mes, 
	                                                                 method = 'igraph',
	                                                                 algorithm = 4,
	                                                                 cluster.name = 'mes.res',
	                                                                 resolution = n)
	col.names <- paste0('mes.res.', seq(0.1,1.0,0.1))
	row.names <- colnames(wound.mes)
	
	wound.mes3 <- data.frame(as.numeric(as.character(wound.mes2[[1]]@meta.data$mes.res)))
	
	for (i in 1:10) {
	  wound.mes3[ ,i] <- data.frame(as.numeric(as.character(wound.mes2[[i]]@meta.data$mes.res)))
	}  
	
	colnames(wound.mes3) <- col.names
	rownames(wound.mes3) <- row.names
	
	wound.mes <- AddMetaData(wound.mes, 
	                         wound.mes3)
	
	rm(wound.mes2, wound.mes3)
	
	gc()
	
	wound.mes <- RunUMAP(wound.mes, 
	                     reduction = 'harmony', 
	                     reduction.name = 'umap.harmony',
	                     dims = 1:90,
	                     assay = 'SCT')
	
	saveRDS(wound.mes,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.mes.rds",
	        compress = FALSE)
	
	#Clustree
	clustree(wound.mes, prefix = 'mes.res.')
	ggsave2(filename = "Wound Skin Mesenchymal Updated Clustree.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(wound.mes) <- 'mes.res.0.1'
	
	DimPlot(wound.mes,
	        label = TRUE,
	        repel = TRUE,
	        pt.size = 8,
	        label.size = 10,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Mesenchymal Updated UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	wound.mes.markers <- FindAllMarkers(wound.mes,
	                                    only.pos = T, 
	                                    test.use = "MAST",
	                                    latent.vars = "sex",
	                                    min.pct = 0.50,
	                                    logfc.threshold = 1.00,
	                                    return.thresh = 0.05,
	                                    assay = "RNA",
	                                    densify = T) 
	
	wound.mes.markers <- wound.mes.markers[order(wound.mes.markers$cluster,
	                                             -wound.mes.markers$avg_log2FC), ]
	
	#Identify how many cells are present in each cluster and by induction time course in each dataset
	
	wound.mes.cellcounts <- table(wound.mes@meta.data$mes.res.0.1,
	                              wound.mes@meta.data$state)
	
	wound.mes.cellcounts2 <- table(wound.mes@meta.data$mes.res.0.1,
	                               wound.mes@meta.data$sampleid)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.mes.markers,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Mesenchymal Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	write.xlsx(wound.mes.cellcounts,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Mesenchymal Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-State",
	           append = TRUE)
	write.xlsx(wound.mes.cellcounts2,
	           file = '/data/overmilleram/scRNAseq/Wound/Wound Skin Mesenchymal Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Cell Counts-Sample ID",
	           append = TRUE)
	
	wound.mes <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.mes.rds")
	
	Idents(wound.mes) = 'mes.res.0.1'
	
	wound.mes <- RenameIdents(wound.mes,
	                          '1' = 'Skeletal_Muscle',
	                          '2' = 'Skeletal_Muscle')
	
	levels(wound.mes) <- c('Skeletal_Muscle')
	
	wound.mes$celltype = wound.mes@active.ident
	
	Idents(wound.mes) = 'mes.res.0.1'
	
	wound.mes <- RenameIdents(wound.mes,
	                           '1' = 'Skeletal_Muscle',
	                           '2' = 'Skeletal_Muscle')
	
	levels(wound.mes) <- c('Skeletal_Muscle')
	
	wound.mes$L1subtype = wound.mes@active.ident
	
	Idents(wound.mes) = 'mes.res.0.1'
	
	wound.mes <- RenameIdents(wound.mes,
	                           '1' = 'Skeletal_Muscle_1',
	                           '2' = 'Skeletal_Muscle_2')
	
	levels(wound.mes) <- c('Skeletal_Muscle_1', 
	                       'Skeletal_Muscle_2')
	
	wound.mes$L2subtype = wound.mes@active.ident
	
	saveRDS(wound.mes,
	        file = "/data/overmilleram/scRNAseq/Wound/woundskin.mes.rds",
	        compress = FALSE)
	
	Idents(wound.mes) = 'L1subtype'
	
	DimPlot(wound.mes,
	        label = TRUE,
	        repel = TRUE,
	        label.size = 7,
	        pt.size = 4,
	        raster = FALSE) + NoLegend()
	ggsave2(filename = "Wound Skin Mesenchymal L1 Subtype UMAP.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	###########
	
	# Pull metadata from subsetted objects and add info back to tissue.integrated as 'L1' & 'L2' subtype
	
	wound.integrated <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds')
	
	wound.kera <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.kera.rds")
	wound.fib <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.fib.rds")
	wound.imm <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds")
	wound.vasc <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.vasc.rds")
	wound.mel <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.mel.rds")
	wound.mes <- readRDS("/data/overmilleram/scRNAseq/Wound/woundskin.mes.rds")
	
	# celltype
	kera.cells <- wound.kera$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	fib.cells <- wound.fib$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	imm.cells <- wound.imm$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	vasc.cells <- wound.vasc$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	mel.cells <- wound.mel$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	mes.cells <- wound.mes$celltype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, fib.cells, imm.cells, vasc.cells, mel.cells, mes.cells)
	
	celltype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(celltype) <- celltype[ ,1]
	celltype[ ,1] <- NULL
	colnames(celltype) <- 'celltype'
	
	# L1subtype
	kera.cells <- wound.kera$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	fib.cells <- wound.fib$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	imm.cells <- wound.imm$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	vasc.cells <- wound.vasc$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mel.cells <- wound.mel$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mes.cells <- wound.mes$L1subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, fib.cells, imm.cells, vasc.cells, mel.cells, mes.cells)
	
	L1subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L1subtype) <- L1subtype[ ,1]
	L1subtype[ ,1] <- NULL
	colnames(L1subtype) <- 'L1subtype'
	
	# L2subtype
	kera.cells <- wound.kera$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	fib.cells <- wound.fib$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	imm.cells <- wound.imm$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	vasc.cells <- wound.vasc$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mel.cells <- wound.mel$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	mes.cells <- wound.mes$L2subtype %>% as.data.frame %>% rownames_to_column('barcode')
	
	cell.list <- list(kera.cells, fib.cells, imm.cells, vasc.cells, mel.cells, mes.cells)
	
	L2subtype <- Reduce(function(x, y) full_join(x, y), cell.list)
	
	rownames(L2subtype) <- L2subtype[ ,1]
	L2subtype[ ,1] <- NULL
	colnames(L2subtype) <- 'L2subtype'
	
	# Add metadata to wound.integrated and subset out unlabeled cells
	
	barcodes <- rownames(L1subtype) %>% as.character
	
	wound.integrated <- subset(wound.integrated,
	                           cells = barcodes)
	
	wound.integrated <- AddMetaData(wound.integrated, celltype, col.name = 'celltype')
	wound.integrated <- AddMetaData(wound.integrated, L1subtype, col.name = 'L1subtype')
	wound.integrated <- AddMetaData(wound.integrated, L2subtype, col.name = 'L2subtype')
	
	wound.integrated <- RunPCA(wound.integrated, 
	                           npcs = 100, 
	                           assay = 'SCT')
	
	wound.integrated <- FindNeighbors(wound.integrated, 
	                                  reduction = 'harmony', 
	                                  assay = 'RNA',
	                                  dims = 1:90)
	
	wound.integrated <- RunUMAP(wound.integrated, 
	                            reduction = 'harmony', 
	                            reduction.name = 'umap.harmony',
	                            assay = 'RNA',
	                            dims = 1:90)
	
	Idents(wound.integrated) <- 'celltype'
	
	umap.wound.celltype <- DimPlot(wound.integrated,
	                              label = TRUE,
	                              repel = TRUE,
	                              pt.size = 1,
	                              label.size = 10,
	                              raster = FALSE) + NoLegend()
	ggsave2(umap.wound.celltype,
	        filename = "Wound Skin Celltype UMAP Final.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(wound.integrated) <- 'L1subtype'
	
	umap.wound.L1 <- DimPlot(wound.integrated,
	                        label = TRUE,
	                        repel = TRUE,
	                        pt.size = 1,
	                        label.size = 6,
	                        raster = FALSE) + NoLegend()
	
	ggsave2(umap.wound.L1,
	        filename = "Wound Skin L1 Subtype UMAP Final.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	Idents(wound.integrated) <- 'L2subtype'
	
	umap.wound.L2 <- DimPlot(wound.integrated,
	                        label = TRUE,
	                        repel = TRUE,
	                        pt.size = 1,
	                        label.size = 4,
	                        raster = FALSE) + NoLegend()
	
	ggsave2(umap.wound.L2,
	        filename = "Wound Skin L2 Subtype UMAP Final.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	saveRDS(wound.integrated,
	        '/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds',
	        compress = F)
	
	# Find DEGs for each comparison in celltype & L1 subtype
	
	wound.integrated <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds')
	
	# Dox vs Control
	##########
	#celltype
	Idents(wound.integrated) = "celltype"
	celltype = as.vector(unique(wound.integrated$celltype))
	
	celltype.markers = list()
	
	for (i in 1:length(celltype)) {
	  tryCatch({
	    celltype.markers[[i]] = FindMarkers(wound.integrated,
	                                        test.use = "MAST",
	                                        latent.vars = "sex",
	                                        logfc.threshold = 1.00,
	                                        min.pct = 0.25,
	                                        only.pos = F,
	                                        assay = "RNA",
	                                        group.by = "condition",
	                                        ident.1 = "dox",
	                                        ident.2 = "control",
	                                        subset.ident = celltype[[i]],
	                                        densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(celltype.markers) = celltype
	
	celltype.markers.sort = list()
	
	for (j in 1:length(celltype)) {
	  
	  celltype.markers.sort[[j]] = celltype.markers[[j]]
	  celltype.markers.sort[[j]]$gene = rownames(celltype.markers[[j]])
	  celltype.markers.sort[[j]] = celltype.markers.sort[[j]][order(-celltype.markers.sort[[j]]$avg_log2FC), ]
	  celltype.markers.sort[[j]] = subset(celltype.markers.sort[[j]], p_val_adj < 0.05)
	  
	  
	}
	
	names(celltype.markers.sort) = celltype
	
	openxlsx::write.xlsx(celltype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Wound/Wound Skin Dox vs Control Celltype DEGs.xlsx")
	
	# L1 subtype
	Idents(wound.integrated) = "L1subtype"
	L1subtype = as.vector(unique(wound.integrated$L1subtype))
	
	L1subtype.markers = list()
	
	for (i in 1:length(L1subtype)) {
	  tryCatch({
	    L1subtype.markers[[i]] = FindMarkers(wound.integrated,
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         logfc.threshold = 1.00,
	                                         min.pct = 0.25,
	                                         only.pos = F,
	                                         assay = "RNA",
	                                         group.by = "condition",
	                                         ident.1 = "dox",
	                                         ident.2 = "control",
	                                         subset.ident = L1subtype[[i]],
	                                         densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L1subtype.markers) = L1subtype
	
	L1subtype.markers.sort = list()
	
	for (j in 1:length(L1subtype)) {
	  tryCatch({
	    L1subtype.markers.sort[[j]] = L1subtype.markers[[j]]
	    L1subtype.markers.sort[[j]]$gene = rownames(L1subtype.markers[[j]])
	    L1subtype.markers.sort[[j]] = L1subtype.markers.sort[[j]][order(-L1subtype.markers.sort[[j]]$avg_log2FC), ]
	    L1subtype.markers.sort[[j]] = subset(L1subtype.markers.sort[[j]], p_val_adj < 0.05)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L1subtype.markers.sort) = L1subtype
	
	openxlsx::write.xlsx(L1subtype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Wound/Wound Skin Dox vs Control L1 Subtype DEGs.xlsx")
	
	
	#L2subtype
	Idents(wound.integrated) = "L2subtype"
	L2subtype = as.vector(unique(wound.integrated$L2subtype))
	
	L2subtype.markers = list()
	
	for (i in 1:length(L2subtype)) {
	  tryCatch({
	    L2subtype.markers[[i]] = FindMarkers(wound.integrated,
	                                         test.use = "MAST",
	                                         latent.vars = "sex",
	                                         logfc.threshold = 1.00,
	                                         min.pct = 0.25,
	                                         only.pos = F,
	                                         assay = "RNA",
	                                         group.by = "condition",
	                                         ident.1 = "dox",
	                                         ident.2 = "control",
	                                         subset.ident = L2subtype[[i]],
	                                         densify = T)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L2subtype.markers) = L2subtype
	
	L2subtype.markers.sort = list()
	
	for (j in 1:length(L2subtype)) {
	  tryCatch({
	    L2subtype.markers.sort[[j]] = L2subtype.markers[[j]]
	    L2subtype.markers.sort[[j]]$gene = rownames(L2subtype.markers[[j]])
	    L2subtype.markers.sort[[j]] = L2subtype.markers.sort[[j]][order(-L2subtype.markers.sort[[j]]$avg_log2FC), ]
	    L2subtype.markers.sort[[j]] = subset(L2subtype.markers.sort[[j]], p_val_adj < 0.05)
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	
	names(L2subtype.markers.sort) = L2subtype
	
	openxlsx::write.xlsx(L2subtype.markers.sort, 
	                     file = "/data/overmilleram/scRNAseq/Wound/Wound Skin Dox vs Control L2 Subtype DEGs.xlsx")
	##########
	
	# enrichR 
	#######################
	#pull gene names from celltype, simplesubtype, and subtype; organize 
	
	# celltype
	celltype.dox.genes = list()
	
	for (i in 1:length(celltype)) {
	  celltype.dox.genes[[i]] = rownames(celltype.markers.sort[[i]])
	}
	
	names(celltype.dox.genes) = celltype
	
	celltype.genes.list = list(celltype.dox.genes)
	
	names(celltype.genes.list) = c('dox')
	
	# L1 subtype
	L1subtype.dox.genes = list()
	
	for (i in 1:length(L1subtype)) {
	  L1subtype.dox.genes[[i]] = rownames(L1subtype.markers.sort[[i]])
	}
	
	names(L1subtype.dox.genes) = L1subtype
	
	L1subtype.genes.list = list(L1subtype.dox.genes)
	
	names(L1subtype.genes.list) = c('dox')
	
	#L2subtype
	L2subtype.dox.genes = list()
	
	for (i in 1:length(L2subtype)) {
	  L2subtype.dox.genes[[i]] = rownames(L2subtype.markers.sort[[i]])
	}
	
	names(L2subtype.dox.genes) = L2subtype
	
	L2subtype.genes.list = list(L2subtype.dox.genes)
	
	names(L2subtype.genes.list) = c('dox')
	
	# run enrichR
	
	setEnrichrSite("Enrichr") # Human/mouse genes
	
	websiteLive = T
	
	all.dbs = listEnrichrDbs()
	
	dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", 
	        'Reactome_2022', 'Tabula_Muris', 'Mouse_Gene_Atlas', 'Jensen_TISSUES', 'KEGG_2019_Mouse', 'MSigDB_Hallmark_2020', 
	        'WikiPathways_2019_Mouse', 'Panther_2016', 'BioPlex_2017', 'ChEA_2022', 'ENCODE_and_ChEA_Consensus_TFs_from_CHIP-X', 
	        'ENCODE_TF_ChIP-seq_2015', 'Epigenomics_Roadmap_HM_ChIP-seq', 'Transcription_Factor_PPIs', 
	        'Tissue_Protein_Expression_from_ProteomicsDB', 'TRRUST_Transcription_Factors_2019', 'KOMP2_Mouse_Phenotypes_2022', 
	        'miRTarBase_2017')
	
	# celltype
	celltype.dox.enrichr = list()
	
	for (i in 1:length(celltype)) {
	  if (websiteLive) {
	    celltype.dox.enrichr[[i]] = enrichr(celltype.dox.genes[[i]], dbs)
	  }
	}
	
	names(celltype.dox.enrichr) = celltype
	
	celltype.enrichr.list = list(celltype.dox.enrichr)
	
	names(celltype.enrichr.list) = c('Dox vs Control')
	
	# L1subtype
	L1subtype.dox.enrichr = list()
	
	for (i in 1:length(L1subtype)) {
	  if (websiteLive) {
	    L1subtype.dox.enrichr[[i]] = enrichr(L1subtype.dox.genes[[i]], dbs)
	  }
	}
	
	names(L1subtype.dox.enrichr) = L1subtype
	
	L1subtype.enrichr.list = list(L1subtype.dox.enrichr)
	
	names(L1subtype.enrichr.list) = c('Dox vs Control')
	
	#L2subtype
	L2subtype.dox.enrichr = list()
	
	for (i in 1:length(L2subtype)) {
	  if (websiteLive) {
	    L2subtype.dox.enrichr[[i]] = enrichr(L2subtype.dox.genes[[i]], dbs)
	  }
	}
	
	names(L2subtype.dox.enrichr) = L2subtype
	
	L2subtype.enrichr.list = list(L2subtype.dox.enrichr)
	
	names(L2subtype.enrichr.list) = c('Dox vs Control')
	
	enrichr.list = list(celltype.enrichr.list, L1subtype.enrichr.list, L2subtype.enrichr.list)
	names(enrichr.list) = c('Celltype', 'L1subtype', 'L2subtype')
	
	saveRDS(enrichr.list,
	        '/data/overmilleram/scRNAseq/Wound/woundskin.enrichrlist.rds',
	        compress = T)
	###########################
	
	#Cellchat
	###################
	
	wound.integrated <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.final.rds')
	
	#celltype
	Idents(wound.integrated) = "state"
	
	skin.control = subset(wound.integrated, 
	                      idents = "wound_control_8_week")
	skin.dox = subset(wound.integrated, 
	                  idents = "wound_dox_8_week")
	
	skco.cellchat = createCellChat(object = skin.control,
	                               assay = "RNA",
	                               group.by = "celltype")
	skdo.cellchat = createCellChat(object = skin.dox,
	                               assay = "RNA",
	                               group.by = "celltype")
	
	skco.cellchat@DB = CellChatDB.mouse
	skdo.cellchat@DB = CellChatDB.mouse
	
	cellchat.celltype.list = list(skco.cellchat, skdo.cellchat)
	
	names(cellchat.celltype.list) = c("Control 8 Week", "+Dox 8 Week")
	
	cellchat.celltype.list = lapply(cellchat.celltype.list, function(x){
	  x = subsetData(x)
	  x = identifyOverExpressedGenes(x)
	  x = identifyOverExpressedInteractions(x)
	  x = computeCommunProb(x, raw.use = T)
	  x = filterCommunication(x, min.cells = 10)
	  x = computeCommunProbPathway(x)
	  x = aggregateNet(x)
	  x = netAnalysis_computeCentrality(x)
	})
	
	saveRDS(cellchat.celltype.list, 
	        '/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.celltype.list.rds',
	        compress = F)
	
	#L1subtype
	Idents(wound.integrated) = "state"
	
	skin.control = subset(wound.integrated, 
	                      idents = "wound_control_8_week")
	skin.dox = subset(wound.integrated, 
	                  idents = "wound_dox_8_week")
	
	skco.cellchat = createCellChat(object = skin.control,
	                               assay = "RNA",
	                               group.by = "L1subtype")
	skdo.cellchat = createCellChat(object = skin.dox,
	                               assay = "RNA",
	                               group.by = "L1subtype")
	
	skco.cellchat@DB = CellChatDB.mouse
	skdo.cellchat@DB = CellChatDB.mouse
	
	cellchat.L1subtype.list = list(skco.cellchat, skdo.cellchat)
	
	names(cellchat.L1subtype.list) = c("Control 8 Week", "+Dox 8 Week")
	
	cellchat.L1subtype.list = lapply(cellchat.L1subtype.list, function(x){
	  x = subsetData(x)
	  x = identifyOverExpressedGenes(x)
	  x = identifyOverExpressedInteractions(x)
	  x = computeCommunProb(x, raw.use = T)
	  x = filterCommunication(x, min.cells = 10)
	  x = computeCommunProbPathway(x)
	  x = aggregateNet(x)
	  x = netAnalysis_computeCentrality(x)
	})
	
	saveRDS(cellchat.L1subtype.list, 
	        '/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.L1subtype.list.rds',
	        compress = F)
	
	#L2subtype
	Idents(wound.integrated) = "state"
	
	skin.control = subset(wound.integrated, 
	                      idents = "wound_control_8_week")
	skin.dox = subset(wound.integrated, 
	                  idents = "wound_dox_8_week")
	
	skco.cellchat = createCellChat(object = skin.control,
	                               assay = "RNA",
	                               group.by = "L2subtype")
	skdo.cellchat = createCellChat(object = skin.dox,
	                               assay = "RNA",
	                               group.by = "L2subtype")
	
	skco.cellchat@DB = CellChatDB.mouse
	skdo.cellchat@DB = CellChatDB.mouse
	
	cellchat.L2subtype.list = list(skco.cellchat, skdo.cellchat)
	
	names(cellchat.L2subtype.list) = c("Control 8 Week", "+Dox 8 Week")
	
	cellchat.L2subtype.list = lapply(cellchat.L2subtype.list, function(x){
	  x = subsetData(x)
	  x = identifyOverExpressedGenes(x)
	  x = identifyOverExpressedInteractions(x)
	  x = computeCommunProb(x, raw.use = T)
	  x = filterCommunication(x, min.cells = 10)
	  x = computeCommunProbPathway(x)
	  x = aggregateNet(x)
	  x = netAnalysis_computeCentrality(x)
	})
	
	saveRDS(cellchat.L2subtype.list, 
	        '/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.L2subtype.list.rds',
	        compress = F)
	
	# analyze CellChat objects
	
	# celltype
	
	cellchat.celltype = mergeCellChat(cellchat.celltype.list, add.names = names(cellchat.celltype.list))
	
	#Number and Strength of Interactions in merged objects
	gg1 = compareInteractions(cellchat.celltype, show.legend = F, group = c(1:2))
	gg2 = compareInteractions(cellchat.celltype, show.legend = F, group = c(1:2), measure = "weight")
	
	gg1|gg2
	ggsave2(filename = "Wound Skin Celltype Interaction Number and Strength.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 4000,
	        height = 2000,
	        units = "px")
	
	#Major sources and targets in 2D space
	
	num.link = sapply(cellchat.celltype.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	
	weight.MinMax = c(min(num.link), max(num.link))
	
	gg = list()
	for (i in 1:length(cellchat.celltype.list)) {
	  gg[[i]] = netAnalysis_signalingRole_scatter(cellchat.celltype.list[[i]], 
	                                              title = names(cellchat.celltype.list)[i], 
	                                              weight.MinMax = weight.MinMax)
	}
	
	patchwork::wrap_plots(plots = gg)
	
	ggsave2(filename = "Wound Skin Celltype Major Sources and Targets 2D.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 6000,
	        height = 3000,
	        units = "px")
	
	#Compare overall info flow of each signaling pathway
	
	#Dox vs Control
	gg1 = rankNet(cellchat.celltype, measure = "weight", mode = "comparison", stacked = T, font.size = 12, comparison = c(1,2), do.stat = TRUE)
	gg2 = rankNet(cellchat.celltype, measure = "weight", mode = "comparison", stacked = F, font.size = 12, comparison = c(1,2), do.stat = TRUE)
	
	gg1/gg2
	
	ggsave2(filename = "Wound Skin Celltype Overall Signaling Pathway Information Flow Combined.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 4000,
	        height = 10000,
	        units = "px",
	        limitsize = F)
	
	#Compare signaling associated with each cell population
	
	pathway.union = unique(c(cellchat.celltype.list[[1]]@netP$pathways, 
	                         cellchat.celltype.list[[2]]@netP$pathways))
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[1]],
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.celltype.list)[1], 
	                                        width = 10, 
	                                        height = 20)
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[2]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[2], 
	                                        width = 10, 
	                                        height = 20)
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(20, "cm"))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype Outgoing Signal Comparison.svg",
	          width = 20,
	          height = 20)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[1]],
	                                        pattern = "incoming", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.celltype.list)[1], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "GnBu")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[2]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[2], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "GnBu")
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(20, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype Incoming Signal Comparison.svg",
	          width = 20,
	          height = 20)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[1]],
	                                        pattern = "all", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.celltype.list)[1], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "OrRd")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.celltype.list[[2]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.celltype.list)[2], 
	                                        width = 10, 
	                                        height = 20, 
	                                        color.heatmap = "OrRd")
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(20, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype Overall Signal Comparison.svg",
	          width = 20,
	          height = 20)
	
	dev.off()
	
	#Cell-cell communication network circos
	
	#Control
	groupSize <- as.numeric(table(cellchat.celltype@idents$`Control 8 Week`))
	mat <- cellchat.celltype@net$`Control 8 Week`$weight
	par(mfrow = c(2,3), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype Control Cell-Cell Network.svg",
	          width = 20,
	          height = 10)
	
	dev.off()
	
	#Dox
	groupSize <- as.numeric(table(cellchat.celltype@idents$`+Dox 8 Week`))
	mat <- cellchat.celltype@net$`+Dox 8 Week`$weight
	par(mfrow = c(2,3), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype +Dox Cell-Cell Network.svg",
	          width = 20,
	          height = 10)
	
	dev.off()
	
	# DEG analysis to identify up/downregulated L-R pairs
	
	#+Dox as comparative dataset (pos.dataset)
	pos.dataset = "+Dox 8 Week"
	features.name = pos.dataset
	
	cellchat.celltype = identifyOverExpressedGenes(cellchat.celltype, 
	                                               group.dataset = "datasets", 
	                                               pos.dataset = pos.dataset, 
	                                               features.name = features.name, 
	                                               only.pos = FALSE, 
	                                               thresh.pc = 0.25, 
	                                               thresh.fc = 0.75, 
	                                               thresh.p = 0.05)
	
	cellchat.celltype = identifyOverExpressedInteractions(cellchat.celltype)
	
	net = netMappingDEG(cellchat.celltype, features.name = features.name)
	
	net.up = subsetCommunication(cellchat.celltype, 
	                             net = net, 
	                             datasets = "+Dox 8 Week",
	                             ligand.logFC = 0.2, 
	                             receptor.logFC = NULL)
	
	net.down = subsetCommunication(cellchat.celltype,
	                               net = net, 
	                               datasets = "Control 8 Week",
	                               ligand.logFC = -0.2, 
	                               receptor.logFC = NULL)
	
	gene.up = extractGeneSubsetFromPair(net.up, cellchat.celltype)
	gene.down = extractGeneSubsetFromPair(net.down, cellchat.celltype)
	
	pairLR.use.up = net.up[, "interaction_name", drop = F]
	pairLR.use.down = net.down[, "interaction_name", drop = F]
	
	par(mfrow = c(1,1), xpd=TRUE)
	
	#Control
	# Downregulated 
	netVisual_chord_gene(cellchat.celltype.list[[1]], 
	                     sources.use = c(1:6), 
	                     targets.use = c(1:6), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.8, 
	                     small.gap = 3.5, 
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.celltype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype +Dox Downregulated Signaling.svg",
	          width = 15,
	          height = 15)
	dev.off()
	
	#Upregulated
	netVisual_chord_gene(cellchat.celltype.list[[2]], 
	                     sources.use = c(1:6), 
	                     targets.use = c(1:6), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.8, 
	                     small.gap = 3.5, 
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.celltype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin Celltype +Dox Upregulated Signaling.svg",
	          width = 15,
	          height = 15)
	dev.off()
	
	# Compare cell-cell communication - SKIPPED
	#####
	# Chord diagram
	#pathways.show <- c("CSPG4") 
	#
	#par(mfrow = c(1,2), xpd=TRUE)
	#
	#for (i in 1:length(cellchat.celltype.list)) {
	#  netVisual_aggregate(cellchat.celltype.list[[i]], 
	#                      signaling = pathways.show, 
	#                      layout = "chord", 
	#                      signaling.name = paste(pathways.show, 
	#                                             names(cellchat.celltype.list)[i]))
	#}
	#
	#dev.print(svg,
	#          "../Healthy Tissue 8 Week/CellChat/Healthy Tissue EPHB Signaling Network.svg",
	#          width = 20,
	#          height = 10)
	#dev.off()
	#
	#control.pathways.show <- list("MK", "CSPG4", "RELN", "IL1", "SEMA7") 
	#
	#par(mfrow = c(2,3), xpd=TRUE)
	#
	#for (i in 1:length(control.pathways.show)) {
	#  netVisual_aggregate(cellchat.celltype.list[[1]], 
	#                      signaling = control.pathways.show[[i]], 
	#                      layout = "chord", 
	#                      signaling.name = paste(control.pathways.show[[i]], 
	#                                             names(cellchat.celltype.list)[1]))
	#}
	#
	#dev.print(svg,
	#          "../Healthy Tissue 8 Week/CellChat/Healthy Tissue Control-Specific Signaling Networks.svg",
	#          width = 30,
	#          height = 20)
	#dev.off()
	#
	#dox.pathways.show <- list("IL10", "NT", "VISFATIN") 
	#
	#par(mfrow = c(1,3), xpd=TRUE)
	#
	#for (i in 1:length(dox.pathways.show)) {
	#  netVisual_aggregate(cellchat.celltype.list[[2]], 
	#                      signaling = dox.pathways.show[[i]], 
	#                      layout = "chord", 
	#                      signaling.name = paste(dox.pathways.show[[i]], 
	#                                             names(cellchat.celltype.list)[2]))
	#}
	#
	#dev.print(svg,
	#          "../Healthy Tissue 8 Week/CellChat/Healthy Tissue Dox-Specific Signaling Networks.svg",
	#          width = 30,
	#          height = 20)
	##dev.off()
	#
	## Compare gene expression distribution
	#
	#cellchat.celltype@meta$datasets = factor(cellchat.celltype@meta$datasets, 
	#                                        levels = c("Control 8 Week", "+Dox 8 Week")) # set factor level
	#
	#pathways.plot = unique(cellchat.celltype@netP$`Control 8 Week`$pathways)
	#
	## didn't work due to a Assay5 issue (FetchData.Assay5(): layer data is not found in the object)
	#
	#for (i in pathways.plot) {
	#  
	#  plotGeneExpression(cellchat.celltype,
	#                     signaling = i, 
	#                     split.by = "datasets", 
	#                     colors.ggplot = T)
	#  
	#  ggsave2(filename = paste0('Wound Skin ', i, ' Celltype Signaling Gene Expression Violin Plots.svg'),
	#          path = "/data/overmilleram/scRNAseq/Wound/CellChat/Celltype NetP Violin Plots",
	#          width = 6000,
	#          height = 4000,
	#          units = "px",
	#          limitsize = F)
	#} 
	#####
	
	# L1subtype
	
	cellchat.L1subtype.list <- readRDS('/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.L1subtype.list.rds')
	
	cellchat.L1subtype = mergeCellChat(cellchat.L1subtype.list, add.names = names(cellchat.L1subtype.list))
	
	#Number and Strength of Interactions in merged objects
	gg1 = compareInteractions(cellchat.L1subtype, show.legend = F, group = c(1:2))
	gg2 = compareInteractions(cellchat.L1subtype, show.legend = F, group = c(1:2), measure = "weight")
	
	gg1|gg2
	ggsave2(filename = "Wound Skin L1 Subtype Interaction Number and Strength.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 4000,
	        height = 2000,
	        units = "px")
	
	#Major sources and targets in 2D space
	
	num.link = sapply(cellchat.L1subtype.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	
	weight.MinMax = c(min(num.link), max(num.link))
	
	gg = list()
	for (i in 1:length(cellchat.L1subtype.list)) {
	  gg[[i]] = netAnalysis_signalingRole_scatter(cellchat.L1subtype.list[[i]], 
	                                              title = names(cellchat.L1subtype.list)[i], 
	                                              weight.MinMax = weight.MinMax)
	}
	
	patchwork::wrap_plots(plots = gg)
	
	ggsave2(filename = "Wound Skin L1 Subtype Major Sources and Targets 2D.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 6000,
	        height = 3000,
	        units = "px")
	
	#Compare overall info flow of each signaling pathway
	
	#Dox vs Control
	gg1 = rankNet(cellchat.L1subtype, measure = "weight", mode = "comparison", stacked = T, font.size = 12, comparison = c(1,2), do.stat = TRUE)
	gg2 = rankNet(cellchat.L1subtype, measure = "weight", mode = "comparison", stacked = F, font.size = 12, comparison = c(1,2), do.stat = TRUE)
	
	gg1/gg2
	
	ggsave2(filename = "Wound Skin L1 Subtype Overall Signaling Pathway Information Flow Combined.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 4000,
	        height = 10000,
	        units = "px",
	        limitsize = F)
	
	#Compare signaling associated with each cell population
	
	pathway.union = unique(c(cellchat.L1subtype.list[[1]]@netP$pathways, 
	                         cellchat.L1subtype.list[[2]]@netP$pathways))
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[1]],
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L1subtype.list)[1], 
	                                        width = 10, 
	                                        height = 25)
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[2]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[2], 
	                                        width = 10, 
	                                        height = 25)
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(25, "cm"))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype Outgoing Signal Comparison.svg",
	          width = 20,
	          height = 25)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[1]],
	                                        pattern = "incoming", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L1subtype.list)[1], 
	                                        width = 10, 
	                                        height = 25, 
	                                        color.heatmap = "GnBu")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[2]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[2], 
	                                        width = 10, 
	                                        height = 25, 
	                                        color.heatmap = "GnBu")
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(25, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype Incoming Signal Comparison.svg",
	          width = 20,
	          height = 25)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[1]],
	                                        pattern = "all", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L1subtype.list)[1], 
	                                        width = 10, 
	                                        height = 25, 
	                                        color.heatmap = "OrRd")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L1subtype.list[[2]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L1subtype.list)[2], 
	                                        width = 10, 
	                                        height = 25, 
	                                        color.heatmap = "OrRd")
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(25, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype Overall Signal Comparison.svg",
	          width = 20,
	          height = 25)
	
	dev.off()
	
	#Cell-cell communication network circos
	
	#Control
	groupSize <- as.numeric(table(cellchat.L1subtype@idents$`Control 8 Week`))
	mat <- cellchat.L1subtype@net$`Control 8 Week`$weight
	par(mfrow = c(5,6), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype Control Cell-Cell Network.svg",
	          width = 25,
	          height = 20)
	
	dev.off()
	
	#Dox
	groupSize <- as.numeric(table(cellchat.L1subtype@idents$`+Dox 8 Week`))
	mat <- cellchat.L1subtype@net$`+Dox 8 Week`$weight
	par(mfrow = c(5,6), xpd=TRUE)
	for (i in 1:nrow(mat)) {
	  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	  mat2[i, ] <- mat[i, ]
	  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	}
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype +Dox Cell-Cell Network.svg",
	          width = 25,
	          height = 20)
	
	dev.off()
	
	# DEG analysis to identify up/downregulated L-R pairs
	
	#+Dox as comparative dataset (pos.dataset)
	pos.dataset = "+Dox 8 Week"
	features.name = pos.dataset
	
	cellchat.L1subtype = identifyOverExpressedGenes(cellchat.L1subtype, 
	                                               group.dataset = "datasets", 
	                                               pos.dataset = pos.dataset, 
	                                               features.name = features.name, 
	                                               only.pos = FALSE, 
	                                               thresh.pc = 0.25, 
	                                               thresh.fc = 0.75, 
	                                               thresh.p = 0.05)
	
	cellchat.L1subtype = identifyOverExpressedInteractions(cellchat.L1subtype)
	
	net = netMappingDEG(cellchat.L1subtype, features.name = features.name)
	
	net.up = subsetCommunication(cellchat.L1subtype, 
	                             net = net, 
	                             datasets = "+Dox 8 Week",
	                             ligand.logFC = 0.2, 
	                             receptor.logFC = NULL)
	
	net.down = subsetCommunication(cellchat.L1subtype,
	                               net = net, 
	                               datasets = "Control 8 Week",
	                               ligand.logFC = -0.2, 
	                               receptor.logFC = NULL)
	
	gene.up = extractGeneSubsetFromPair(net.up, cellchat.L1subtype)
	gene.down = extractGeneSubsetFromPair(net.down, cellchat.L1subtype)
	
	pairLR.use.up = net.up[, "interaction_name", drop = F]
	pairLR.use.down = net.down[, "interaction_name", drop = F]
	
	par(mfrow = c(1,1), xpd=TRUE)
	
	#Control
	# Downregulated 
	netVisual_chord_gene(cellchat.L1subtype.list[[1]], 
	                     sources.use = c(14:26), 
	                     targets.use = c(1:26), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L1subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype +Dox Downregulated Signaling 2.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	#Upregulated
	netVisual_chord_gene(cellchat.L1subtype.list[[2]], 
	                     sources.use = c(14:26), 
	                     targets.use = c(1:26), 
	                     slot.name = 'net', 
	                     net = net.up,  
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5, 
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L1subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L1 Subtype +Dox Upregulated Signaling 2.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	# L2subtype
	
	cellchat.L2subtype.list <- readRDS('/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.L2subtype.list.rds')
	
	cellchat.L2subtype = mergeCellChat(cellchat.L2subtype.list, add.names = names(cellchat.L2subtype.list))
	
	#Number and Strength of Interactions in merged objects
	gg1 = compareInteractions(cellchat.L2subtype, show.legend = F, group = c(1:2))
	gg2 = compareInteractions(cellchat.L2subtype, show.legend = F, group = c(1:2), measure = "weight")
	
	gg1|gg2
	ggsave2(filename = "Wound Skin L2 Subtype Interaction Number and Strength.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 4000,
	        height = 2000,
	        units = "px")
	
	#Major sources and targets in 2D space
	
	num.link = sapply(cellchat.L2subtype.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
	
	weight.MinMax = c(min(num.link), max(num.link))
	
	gg = list()
	for (i in 1:length(cellchat.L2subtype.list)) {
	  gg[[i]] = netAnalysis_signalingRole_scatter(cellchat.L2subtype.list[[i]], 
	                                              title = names(cellchat.L2subtype.list)[i], 
	                                              weight.MinMax = weight.MinMax)
	}
	
	patchwork::wrap_plots(plots = gg)
	
	ggsave2(filename = "Wound Skin L2 Subtype Major Sources and Targets 2D.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 6000,
	        height = 3000,
	        units = "px")
	
	#Compare overall info flow of each signaling pathway
	
	#Dox vs Control
	gg1 = rankNet(cellchat.L2subtype, measure = "weight", mode = "comparison", stacked = T, font.size = 12, comparison = c(1,2), do.stat = TRUE)
	gg2 = rankNet(cellchat.L2subtype, measure = "weight", mode = "comparison", stacked = F, font.size = 12, comparison = c(1,2), do.stat = TRUE)
	
	gg1/gg2
	
	ggsave2(filename = "Wound Skin L2 Subtype Overall Signaling Pathway Information Flow Combined.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/CellChat',
	        width = 4000,
	        height = 15000,
	        units = "px",
	        limitsize = F)
	
	#Compare signaling associated with each cell population
	
	pathway.union = unique(c(cellchat.L2subtype.list[[1]]@netP$pathways, 
	                         cellchat.L2subtype.list[[2]]@netP$pathways))
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L2subtype.list[[1]],
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L2subtype.list)[1], 
	                                        width = 15, 
	                                        height = 40)
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L2subtype.list[[2]], 
	                                        pattern = "outgoing", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L2subtype.list)[2], 
	                                        width = 15, 
	                                        height = 40)
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(40, "cm"))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype Outgoing Signal Comparison.svg",
	          width = 30,
	          height = 40)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L2subtype.list[[1]],
	                                        pattern = "incoming", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L2subtype.list)[1], 
	                                        width = 15, 
	                                        height = 40, 
	                                        color.heatmap = "GnBu")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L2subtype.list[[2]], 
	                                        pattern = "incoming", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L2subtype.list)[2], 
	                                        width = 15, 
	                                        height = 40, 
	                                        color.heatmap = "GnBu")
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(40, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype Incoming Signal Comparison.svg",
	          width = 30,
	          height = 40)
	
	dev.off()
	
	ht1 = netAnalysis_signalingRole_heatmap(cellchat.L2subtype.list[[1]],
	                                        pattern = "all", 
	                                        signaling = pathway.union, 
	                                        title = names(cellchat.L2subtype.list)[1], 
	                                        width = 15, 
	                                        height = 40, 
	                                        color.heatmap = "OrRd")
	
	ht2 = netAnalysis_signalingRole_heatmap(cellchat.L2subtype.list[[2]], 
	                                        pattern = "all", 
	                                        signaling = pathway.union,
	                                        title = names(cellchat.L2subtype.list)[2], 
	                                        width = 15, 
	                                        height = 40, 
	                                        color.heatmap = "OrRd")
	
	draw(ht1+ht2, ht_gap = unit(1, "cm"), height = unit(40, "cm"))
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype Overall Signal Comparison.svg",
	          width = 30,
	          height = 40)
	
	dev.off()
	
	#Cell-cell communication network circos
	
	#Control
	#groupSize <- as.numeric(table(cellchat.L2subtype@idents$`Control 8 Week`))
	#mat <- cellchat.L2subtype@net$`Control 8 Week`$weight
	#par(mfrow = c(8,7), xpd=TRUE)
	#for (i in 1:nrow(mat)) {
	#  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	#  mat2[i, ] <- mat[i, ]
	#  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	#}
	#
	#dev.print(svg,
	#          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype Control Cell-Cell Network.svg",
	#          width = 30,
	#          height = 35)
	#
	#dev.off()
	#
	##Dox
	#groupSize <- as.numeric(table(cellchat.L2subtype@idents$`+Dox 8 Week`))
	#mat <- cellchat.L2subtype@net$`+Dox 8 Week`$weight
	#par(mfrow = c(8,7), xpd=TRUE)
	#for (i in 1:nrow(mat)) {
	#  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
	#  mat2[i, ] <- mat[i, ]
	#  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
	#}
	#
	#dev.print(svg,
	#          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Cell-Cell Network.svg",
	#          width = 30,
	#          height = 35)
	#
	#dev.off()
	
	# DEG analysis to identify up/downregulated L-R pairs
	
	# +Dox as comparative dataset (pos.dataset)
	pos.dataset = "+Dox 8 Week"
	features.name = pos.dataset
	
	cellchat.L2subtype = identifyOverExpressedGenes(cellchat.L2subtype, 
	                                                group.dataset = "datasets", 
	                                                pos.dataset = pos.dataset, 
	                                                features.name = features.name, 
	                                                only.pos = FALSE, 
	                                                thresh.pc = 0.25, 
	                                                thresh.fc = 0.75, 
	                                                thresh.p = 0.05)
	
	cellchat.L2subtype = identifyOverExpressedInteractions(cellchat.L2subtype)
	
	net = netMappingDEG(cellchat.L2subtype, features.name = features.name)
	
	net.up = subsetCommunication(cellchat.L2subtype, 
	                             net = net, 
	                             datasets = "+Dox 8 Week",
	                             ligand.logFC = 0.2, 
	                             receptor.logFC = NULL)
	
	net.down = subsetCommunication(cellchat.L2subtype,
	                               net = net, 
	                               datasets = "Control 8 Week",
	                               ligand.logFC = -0.2, 
	                               receptor.logFC = NULL)
	
	gene.up = extractGeneSubsetFromPair(net.up, cellchat.L2subtype)
	gene.down = extractGeneSubsetFromPair(net.down, cellchat.L2subtype)
	
	pairLR.use.up = net.up[, "interaction_name", drop = F]
	pairLR.use.down = net.down[, "interaction_name", drop = F]
	
	par(mfrow = c(1,1), xpd=TRUE)
	
	#Control
	# Downregulated 
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(1:11), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down,  
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Keratinocyte Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(12:17), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Fibroblast Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(19:22,24), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5, 
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Lymhpoid Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(23,25:30), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Myeloid Sender 1.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(31:32), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Myeloid Sender 2.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(33), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Myeloid Sender 3.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(34), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Myeloid Sender 4.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(35:39), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Myeloid Sender 5.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(18,40:43), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.down, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Vascular & Mesenchymal Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	#Upregulated
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(1:11), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Keratinocyte Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(12:17), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up,
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Fibroblast Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(19:22,24), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Lymhpoid Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(23,25:30), 
	                     targets.use = c(1:56), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Myeloid Sender 1.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(31:32), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Myeloid Sender 2.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(33), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Myeloid Sender 3.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[2]], 
	                     sources.use = c(34), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5, 
	                     title.name = paste0("Up-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Upregulated Signaling Myeloid Sender 4.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(35:39), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Myeloid Sender 5.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	netVisual_chord_gene(cellchat.L2subtype.list[[1]], 
	                     sources.use = c(18,40:43), 
	                     targets.use = c(1:43), 
	                     slot.name = 'net', 
	                     net = net.up, 
	                     lab.cex = 0.4, 
	                     small.gap = 1, 
	                     big.gap = 5,
	                     title.name = paste0("Down-regulated signaling in ", 
	                                         names(cellchat.L2subtype.list)[2]))
	
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/CellChat/Wound Skin L2 Subtype +Dox Downregulated Signaling Vascular & Mesenchymal Sender.svg",
	          width = 20,
	          height = 20)
	dev.off()
	
	#Compare gene expression distribution
	
	cellchat.celltype.list <- readRDS('/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.celltype.list.rds')
	cellchat.celltype = mergeCellChat(cellchat.celltype.list, add.names = names(cellchat.celltype.list))
	
	cellchat.L1subtype.list <- readRDS('/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.L1subtype.list.rds')
	cellchat.L1subtype = mergeCellChat(cellchat.L1subtype.list, add.names = names(cellchat.L1subtype.list))
	
	cellchat.L2subtype.list <- readRDS('/data/overmilleram/scRNAseq/Wound/CellChat/woundskin.cellchat.L2subtype.list.rds')
	cellchat.L2subtype = mergeCellChat(cellchat.L2subtype.list, add.names = names(cellchat.L2subtype.list))
	
	
	cellchat.celltype@meta$datasets = factor(cellchat.celltype@meta$datasets, 
	                                         levels = c("Control 8 Week", "+Dox 8 Week",
	                                                    'Buccal', 'Hard Palate')) # set factor level
	
	cellchat.L1subtype@meta$datasets = factor(cellchat.L1subtype@meta$datasets, 
	                                         levels = c("Control 8 Week", "+Dox 8 Week",
	                                                    'Buccal', 'Hard Palate')) # set factor level
	
	cellchat.L2subtype@meta$datasets = factor(cellchat.L2subtype@meta$datasets, 
	                                         levels = c("Control 8 Week", "+Dox 8 Week",
	                                                    'Buccal', 'Hard Palate')) # set factor level
	
	celltype.pathways.plot = c('MK', 'OSM', 'ncWNT', 'CSPG4', 'EDN', 'DESMOSOME', 'IL1', 'PDGF', 'SEMA7', 'ESAM', 'HSPG',
	                           'MPZ', 'SELPLG', 'THY1', 'ICAM', 'CDH5', 'WNT', 
	                           'THBS', 'CCL', 'FN1', 'VEGF', 'CD45', 'IL6')
	
	for (i in celltype.pathways.plot) {
	  
	  plotGeneExpression(cellchat.celltype,
	                     signaling = i, 
	                     split.by = "datasets", 
	                     colors.ggplot = T)
	  
	  ggsave2(filename = paste0('Wound Skin Celltype ', i, ' Signaling Gene Expression Violin Plots.svg'),
	          path = "/data/overmilleram/scRNAseq/Wound/CellChat/Gene Expression Violin Plots",
	          width = 4000,
	          height = 4000,
	          units = "px",
	          limitsize = F)
	}
	
	L1subtype.pathways.plot = c('MK', 'NPNT', 'ncWNT', 'MPZ', 'DESMOSOME', 'ESAM', 'SPP1', 'TNF', 'PDGF', 'PTPRM', 'ACTIVIN', 'JAM',
	                            'SELE', 'EPHA', 'LAMININ', 'TGFb', 'THY1', 'PERIOSTIN', 'EPHB', 'FN1', 'BMP', 'ANGPTL', 'VEGF', 'IGF', 
	                            'PECAM1', 'EGF', 'CDH', 'NOTCH', 'IL1', 'CCL', 'CXCL', 'PTN', 'KIT', 'THBS', 'FGF', 'CDH1', 'OSM', 'IL6', 
	                            'IL2', 'HGF', 'GDF', 'IL4')
	
	for (i in L1subtype.pathways.plot) {
	  
	  plotGeneExpression(cellchat.L1subtype,
	                     signaling = i, 
	                     split.by = "datasets", 
	                     colors.ggplot = T)
	  
	  ggsave2(filename = paste0('Wound Skin L1 Subtype ', i, ' Signaling Gene Expression Violin Plots.svg'),
	          path = "/data/overmilleram/scRNAseq/Wound/CellChat/Gene Expression Violin Plots",
	          width = 4000,
	          height = 4000,
	          units = "px",
	          limitsize = F)
	}
	
	
	
	
	###################
	
	# MultiNichetR
	###################
	
	wound.integrated <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds')
	
	#Load LR matrix and other files manually
	lr_network = readRDS("/data/overmilleram/scRNAseq/lr_network_mouse.rds")
	ligand_target_matrix = readRDS("/data/overmilleram/scRNAseq/ligand_target_matrix_mouse.rds")
	lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
	colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
	rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
	
	#Convert Seurat object to SCE
	
	table(wound.integrated$celltype, wound.integrated$state)
	
	Idents(wound.integrated) <- 'celltype' # remove Melanocyte & Mesenchymal cells from analysis
	
	wound.sce <- subset(wound.integrated,
	                    idents = c('Melanocyte', 'Mesenchymal'),
	                    invert = T)
	
	wound.sce$celltype <- droplevels(wound.sce$celltype) # remove unused factor levels
	wound.sce$L1subtype <- droplevels(wound.sce$L1subtype) # remove unused factor levels
	wound.sce$L2subtype <- droplevels(wound.sce$L2subtype) # remove unused factor levels
	
	wound.sce <- as.SingleCellExperiment(wound.sce, assay = "RNA")
	
	wound.sce <- alias_to_symbol_SCE(wound.sce, "mouse")
	
	#define metadata info
	sample_id = "sampleid"
	group_id = "state"
	celltype_id = "celltype"
	covariates = "sex"
	batches = NA
	
	#define senders/receivers
	senders_oi = SummarizedExperiment::colData(wound.sce)[,celltype_id] %>% unique()
	receivers_oi = SummarizedExperiment::colData(wound.sce)[,celltype_id] %>% unique()
	
	#define minimum # of cells to consider for analysis
	min_cells = 20
	
	#Define ligand activity analysis parameters
	
	logFC_threshold <- 0.50
	p_val_threshold <- 0.05
	fraction_cutoff <- 0.05
	p_val_adj <- F 
	empirical_pval <- F
	
	top_n_target <-250
	verbose <- TRUE
	cores_system <- as.numeric(parallel::detectCores())
	n.cores <- min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type
	
	#define prioritization weights, prepare grouping objects
	
	prioritizing_weights_DE <- c("de_ligand" = 1,
	                             "de_receptor" = 1)
	prioritizing_weights_activity <- c("activity_scaled" = 2)
	
	prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
	                                                 "exprs_receptor" = 2)
	
	prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)
	
	prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
	                                              "abund_receiver" = 0)
	
	prioritizing_weights <- c(prioritizing_weights_DE, 
	                          prioritizing_weights_activity, 
	                          prioritizing_weights_expression_specificity,
	                          prioritizing_weights_expression_sufficiency, 
	                          prioritizing_weights_relative_abundance)
	
	#Define contrasts 
	contrasts_oi <- c("'wound_control_8_week-wound_dox_8_week','wound_dox_8_week-wound_control_8_week'")
	
	contrast_tbl <- tibble(contrast = c("wound_control_8_week-wound_dox_8_week","wound_dox_8_week-wound_control_8_week"),
	                       group = c("wound_control_8_week","wound_dox_8_week")) 
	
	abundance_expression_info = get_abundance_expression_info(sce = wound.sce, 
	                                                          sample_id = sample_id, 
	                                                          group_id = group_id, 
	                                                          celltype_id = celltype_id,
	                                                          min_cells = min_cells, 
	                                                          senders_oi = senders_oi, 
	                                                          receivers_oi = receivers_oi, 
	                                                          lr_network = lr_network, 
	                                                          batches = batches)
	
	#differential abundance per group
	abundance_expression_info$abund_plot_group + NoLegend() + theme(axis.text.x=element_text(angle=30,hjust=1))
	ggsave2(filename = "Wound Skin Celltype Abundances by State.svg",
	        path = "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#perform DE analysis for each cell type 
	DE_info <- get_DE_info(sce = wound.sce, 
	                       sample_id = sample_id, 
	                       group_id = group_id, 
	                       celltype_id = celltype_id, 
	                       batches = batches, 
	                       covariates = covariates, 
	                       contrasts_oi = contrasts_oi, 
	                       min_cells = min_cells)
	
	#Combine DE info for ligand-senders and receptors-receivers
	celltype_de = DE_info$celltype_de$de_output_tidy
	
	sender_receiver_de <- combine_sender_receiver_de(sender_de = celltype_de,
	                                                 receiver_de = celltype_de,
	                                                 senders_oi = senders_oi,
	                                                 receivers_oi = receivers_oi,
	                                                 lr_network = lr_network)
	
	#Run NicheNet ligand activity analysis
	
	ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
	  receiver_de = celltype_de,
	  receivers_oi = receivers_oi,
	  ligand_target_matrix = ligand_target_matrix,
	  logFC_threshold = logFC_threshold,
	  p_val_threshold = p_val_threshold,
	  p_val_adj = p_val_adj,
	  top_n_target = top_n_target,
	  verbose = verbose, 
	  n.cores = n.cores
	)
	
	sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
	
	metadata_combined <- SummarizedExperiment::colData(wound.sce) %>% tibble::as_tibble()
	
	grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
	
	colnames(grouping_tbl) <- c("sample","group")
	
	#run prioritization
	
	prioritization_tables <- generate_prioritization_tables(
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de = sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  contrast_tbl = contrast_tbl,
	  sender_receiver_tbl = sender_receiver_tbl,
	  grouping_tbl = grouping_tbl,
	  prioritizing_weights = prioritizing_weights,
	  fraction_cutoff = fraction_cutoff, 
	  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
	  abundance_data_sender = abundance_expression_info$abundance_data_sender
	)
	
	#Add info on prior knowledge and expresion correlation between ligand-receptor and target expression
	
	lr_target_prior_cor <- lr_target_prior_cor_inference(receivers_oi, 
	                                                     abundance_expression_info, 
	                                                     celltype_de, 
	                                                     grouping_tbl, 
	                                                     prioritization_tables, 
	                                                     ligand_target_matrix, 
	                                                     logFC_threshold = logFC_threshold, 
	                                                     p_val_threshold = p_val_threshold, 
	                                                     p_val_adj = p_val_adj)
	
	#save outputs
	multinichenet_celltype_output <- list(
	  celltype_info = abundance_expression_info$celltype_info,
	  celltype_de = celltype_de,
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de =  sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  prioritization_tables = prioritization_tables,
	  grouping_tbl = grouping_tbl,
	  lr_target_prior_cor = lr_target_prior_cor
	) %>% make_lite_output()
	
	saveRDS(multinichenet_celltype_output, 
	        "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin_multinichenet_celltype_output.rds", 
	        compress = T)
	
	##########
	
	#L1 Subtype
	##########
	
	wound.integrated <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds')
	
	#Load LR matrix and other files manually
	lr_network = readRDS("/data/overmilleram/scRNAseq/lr_network_mouse.rds")
	ligand_target_matrix = readRDS("/data/overmilleram/scRNAseq/ligand_target_matrix_mouse.rds")
	lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
	colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
	rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
	
	#Convert Seurat object to SCE
	
	table(wound.integrated$celltype, wound.integrated$state)
	
	Idents(wound.integrated) <- 'celltype' # remove Melanocyte & Mesenchymal cells from analysis
	
	wound.sce <- subset(wound.integrated,
	                    idents = c('Melanocyte', 'Mesenchymal'),
	                    invert = T)
	
	wound.sce$celltype <- droplevels(wound.sce$celltype) # remove unused factor levels
	wound.sce$L1subtype <- droplevels(wound.sce$L1subtype) # remove unused factor levels
	wound.sce$L2subtype <- droplevels(wound.sce$L2subtype) # remove unused factor levels
	
	wound.sce <- as.SingleCellExperiment(wound.sce, assay = "RNA")
	
	wound.sce <- alias_to_symbol_SCE(wound.sce, "mouse")
	
	#define metadata info
	sample_id = "sampleid"
	group_id = "state"
	celltype_id = "L1subtype"
	covariates = "sex"
	batches = NA
	
	#define senders/receivers
	senders_oi = SummarizedExperiment::colData(wound.sce)[,celltype_id] %>% unique()
	receivers_oi = SummarizedExperiment::colData(wound.sce)[,celltype_id] %>% unique()
	
	#define minimum # of cells to consider for analysis
	min_cells = 10
	
	#Define ligand activity analysis parameters
	
	logFC_threshold <- 0.50
	p_val_threshold <- 0.05
	fraction_cutoff <- 0.05
	p_val_adj <- F 
	empirical_pval <- F
	
	top_n_target <-250
	verbose <- TRUE
	cores_system <- as.numeric(parallel::detectCores())
	n.cores <- min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type
	
	#define prioritization weights, prepare grouping objects
	
	prioritizing_weights_DE <- c("de_ligand" = 1,
	                             "de_receptor" = 1)
	prioritizing_weights_activity <- c("activity_scaled" = 2)
	
	prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
	                                                 "exprs_receptor" = 2)
	
	prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)
	
	prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
	                                              "abund_receiver" = 0)
	
	prioritizing_weights <- c(prioritizing_weights_DE, 
	                          prioritizing_weights_activity, 
	                          prioritizing_weights_expression_specificity,
	                          prioritizing_weights_expression_sufficiency, 
	                          prioritizing_weights_relative_abundance)
	
	#Define contrasts 
	contrasts_oi <- c("'wound_control_8_week-wound_dox_8_week','wound_dox_8_week-wound_control_8_week'")
	
	contrast_tbl <- tibble(contrast = c("wound_control_8_week-wound_dox_8_week","wound_dox_8_week-wound_control_8_week"),
	                       group = c("wound_control_8_week","wound_dox_8_week")) 
	
	abundance_expression_info = get_abundance_expression_info(sce = wound.sce, 
	                                                          sample_id = sample_id, 
	                                                          group_id = group_id, 
	                                                          celltype_id = celltype_id,
	                                                          min_cells = min_cells, 
	                                                          senders_oi = senders_oi, 
	                                                          receivers_oi = receivers_oi, 
	                                                          lr_network = lr_network, 
	                                                          batches = batches)
	
	#differential abundance per group
	abundance_expression_info$abund_plot_group + NoLegend() + theme(axis.text.x=element_text(angle=30,hjust=1))
	ggsave2(filename = "Wound Skin L1 Subtype Abundances by State.svg",
	        path = "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#perform DE analysis for each cell type 
	DE_info <- get_DE_info(sce = wound.sce, 
	                       sample_id = sample_id, 
	                       group_id = group_id, 
	                       celltype_id = celltype_id, 
	                       batches = batches, 
	                       covariates = covariates, 
	                       contrasts_oi = contrasts_oi, 
	                       min_cells = min_cells)
	
	#Combine DE info for ligand-senders and receptors-receivers
	celltype_de = DE_info$celltype_de$de_output_tidy
	
	sender_receiver_de <- combine_sender_receiver_de(sender_de = celltype_de,
	                                                 receiver_de = celltype_de,
	                                                 senders_oi = senders_oi,
	                                                 receivers_oi = receivers_oi,
	                                                 lr_network = lr_network)
	
	#Run NicheNet ligand activity analysis
	
	ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
	  receiver_de = celltype_de,
	  receivers_oi = receivers_oi,
	  ligand_target_matrix = ligand_target_matrix,
	  logFC_threshold = logFC_threshold,
	  p_val_threshold = p_val_threshold,
	  p_val_adj = p_val_adj,
	  top_n_target = top_n_target,
	  verbose = verbose, 
	  n.cores = n.cores
	)
	
	sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
	
	metadata_combined <- SummarizedExperiment::colData(wound.sce) %>% tibble::as_tibble()
	
	grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
	
	colnames(grouping_tbl) <- c("sample","group")
	
	#run prioritization
	
	prioritization_tables <- generate_prioritization_tables(
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de = sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  contrast_tbl = contrast_tbl,
	  sender_receiver_tbl = sender_receiver_tbl,
	  grouping_tbl = grouping_tbl,
	  prioritizing_weights = prioritizing_weights,
	  fraction_cutoff = fraction_cutoff, 
	  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
	  abundance_data_sender = abundance_expression_info$abundance_data_sender
	)
	
	#Add info on prior knowledge and expresion correlation between ligand-receptor and target expression
	
	lr_target_prior_cor <- lr_target_prior_cor_inference(receivers_oi, 
	                                                     abundance_expression_info, 
	                                                     celltype_de, 
	                                                     grouping_tbl, 
	                                                     prioritization_tables, 
	                                                     ligand_target_matrix, 
	                                                     logFC_threshold = logFC_threshold, 
	                                                     p_val_threshold = p_val_threshold, 
	                                                     p_val_adj = p_val_adj)
	
	#save outputs
	multinichenet_celltype_output <- list(
	  celltype_info = abundance_expression_info$celltype_info,
	  celltype_de = celltype_de,
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de =  sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  prioritization_tables = prioritization_tables,
	  grouping_tbl = grouping_tbl,
	  lr_target_prior_cor = lr_target_prior_cor
	) %>% make_lite_output()
	
	saveRDS(multinichenet_celltype_output, 
	        "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin_multinichenet_L1subtype_output.rds", 
	        compress = T)
	
	##########
	
	#L2 Subtype
	##########
	
	#Convert Seurat object to SCE
	
	table(wound.integrated$condition, wound.integrated$L2subtype)
	
	wound.sce <- wound.integrated
	
	wound.sce@meta.data$L2subtype <- droplevels(wound.sce@meta.data$L2subtype, 
	                                            exclude = setdiff(levels(wound.sce@meta.data$L2subtype),
	                                                              unique(wound.sce@meta.data$L2subtype))) #must drop unused factor levels
	
	skin.sce <- as.SingleCellExperiment(wound.sce, assay = "RNA")
	
	skin.sce <- alias_to_symbol_SCE(skin.sce, "mouse")
	
	colData(skin.sce)$sampleid <- colData(skin.sce)$sampleid %>% make.names()
	colData(skin.sce)$condition <- colData(skin.sce)$condition %>% make.names()
	colData(skin.sce)$L2subtype <- colData(skin.sce)$L2subtype %>% make.names()
	colData(skin.sce)$sex <- colData(skin.sce)$sex %>% make.names()
	colData(skin.sce)$batch <- colData(skin.sce)$batch %>% make.names()
	
	#define metadata info
	sample_id = "sampleid"
	group_id = "condition"
	celltype_id = "L2subtype"
	covariates = 'sex'
	batches = NA
	
	#define senders/receivers
	senders_oi <- SummarizedExperiment::colData(skin.sce)[,celltype_id] %>% unique()
	receivers_oi <- SummarizedExperiment::colData(skin.sce)[,celltype_id] %>% unique()
	
	#define minimum # of cells to consider for analysis
	min_cells <- 10
	
	#Define ligand activity analysis parameters
	
	logFC_threshold <- 0.50
	p_val_threshold <- 0.05
	fraction_cutoff <- 0.05
	p_val_adj <- TRUE 
	empirical_pval <- F
	
	top_n_target <-250
	verbose <- TRUE
	cores_system <- as.numeric(parallel::detectCores())
	n.cores <- min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type
	
	#define prioritization weights, prepare grouping objects
	
	prioritizing_weights_DE <- c("de_ligand" = 1,
	                             "de_receptor" = 1)
	prioritizing_weights_activity <- c("activity_scaled" = 2)
	
	prioritizing_weights_expression_specificity <- c("exprs_ligand" = 2,
	                                                 "exprs_receptor" = 2)
	
	prioritizing_weights_expression_sufficiency <- c("frac_exprs_ligand_receptor" = 1)
	
	prioritizing_weights_relative_abundance <- c( "abund_sender" = 0,
	                                              "abund_receiver" = 0)
	
	prioritizing_weights <- c(prioritizing_weights_DE, 
	                          prioritizing_weights_activity, 
	                          prioritizing_weights_expression_specificity,
	                          prioritizing_weights_expression_sufficiency, 
	                          prioritizing_weights_relative_abundance)
	
	#Define contrasts 
	contrasts_oi <- c("'control-dox','dox-control'")
	
	contrast_tbl <- tibble(contrast = c("control-dox","dox-control"),
	                       group = c("control","dox")) 
	
	abundance_expression_info <- get_abundance_expression_info(sce = skin.sce, 
	                                                           sample_id = sample_id, 
	                                                           group_id = group_id, 
	                                                           celltype_id = celltype_id,
	                                                           min_cells = min_cells, 
	                                                           senders_oi = senders_oi, 
	                                                           receivers_oi = receivers_oi, 
	                                                           lr_network = lr_network, 
	                                                           batches = batches)
	
	#differential abundance per group
	abundance_expression_info$abund_plot_group + NoLegend() + theme(axis.text.x=element_text(angle=30,hjust=1))
	ggsave2(filename = "Wound Skin L2 Subtype Abundances by Condition.svg",
	        path = "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/",
	        width = 6000,
	        height =6000,
	        units = "px")
	
	#perform DE analysis for each cell type 
	DE_info <- get_DE_info(sce = skin.sce, 
	                       sample_id = sample_id, 
	                       group_id = group_id, 
	                       celltype_id = celltype_id, 
	                       batches = batches, 
	                       covariates = covariates, 
	                       contrasts_oi = contrasts_oi, 
	                       min_cells = min_cells)
	
	#Combine DE info for ligand-senders and receptors-receivers
	celltype_de <- DE_info$celltype_de$de_output_tidy
	
	sender_receiver_de <- combine_sender_receiver_de(sender_de = celltype_de,
	                                                 receiver_de = celltype_de,
	                                                 senders_oi = senders_oi,
	                                                 receivers_oi = receivers_oi,
	                                                 lr_network = lr_network)
	
	
	
	#Run NicheNet ligand activity analysis
	
	ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
	  receiver_de = celltype_de,
	  receivers_oi = receivers_oi,
	  ligand_target_matrix = ligand_target_matrix,
	  logFC_threshold = logFC_threshold,
	  p_val_threshold = p_val_threshold,
	  p_val_adj = p_val_adj,
	  top_n_target = top_n_target,
	  verbose = verbose, 
	  n.cores = n.cores
	)
	
	sender_receiver_tbl <- sender_receiver_de %>% dplyr::distinct(sender, receiver)
	
	metadata_combined <- SummarizedExperiment::colData(skin.sce) %>% tibble::as_tibble()
	
	grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
	
	colnames(grouping_tbl) <- c("sample","group")
	
	#run prioritization
	
	prioritization_tables <- generate_prioritization_tables(
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de = sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  contrast_tbl = contrast_tbl,
	  sender_receiver_tbl = sender_receiver_tbl,
	  grouping_tbl = grouping_tbl,
	  prioritizing_weights = prioritizing_weights,
	  fraction_cutoff = fraction_cutoff, 
	  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
	  abundance_data_sender = abundance_expression_info$abundance_data_sender
	)
	
	#Add info on prior knowledge and expresion correlation between ligand-receptor and target expression
	
	lr_target_prior_cor <- lr_target_prior_cor_inference(receivers_oi, 
	                                                     abundance_expression_info, 
	                                                     celltype_de, 
	                                                     grouping_tbl, 
	                                                     prioritization_tables, 
	                                                     ligand_target_matrix, 
	                                                     logFC_threshold = logFC_threshold, 
	                                                     p_val_threshold = p_val_threshold, 
	                                                     p_val_adj = p_val_adj)
	
	#save outputs
	multinichenet_skin_L2subtype_output <- list(
	  celltype_info = abundance_expression_info$celltype_info,
	  celltype_de = celltype_de,
	  sender_receiver_info = abundance_expression_info$sender_receiver_info,
	  sender_receiver_de =  sender_receiver_de,
	  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
	  prioritization_tables = prioritization_tables,
	  grouping_tbl = grouping_tbl,
	  lr_target_prior_cor = lr_target_prior_cor
	) %>% make_lite_output()
	
	saveRDS(multinichenet_skin_L2subtype_output, 
	        "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/healthyoral.multinichenet.L2subtype.rds", 
	        compress = T)
	
	##########
	#######################
	
	example <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_celltype_output.rds')
	
	multinichenet_woundskin_celltype_output <- readRDS('/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin.multinichenet.celltype.rds')
	multinichenet_woundskin_L1subtype_output <- readRDS('/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin.multinichenet.L1subtype.rds')
	
	#Analyze MultiNicheNet data
	
	#Healthy Tissue (all states)
	#######################
	
	#Celltype
	###############
	
	#Define contrasts 
	
	contrasts_oi <- c("'wound_control_8_week-wound_dox_8_week','wound_dox_8_week-wound_control_8_week'")
	
	contrast_tbl <- tibble(contrast = c("wound_control_8_week-wound_dox_8_week","wound_dox_8_week-wound_control_8_week"),
	                       group = c("wound_control_8_week","wound_dox_8_week"))
	
	#visualization of top 50 interactions across all groups
	prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 50, rank_per_group = FALSE)
	
	prioritized_tbl_oi = multinichenet_celltype_output$prioritization_tables$group_prioritization_tbl %>%
	  filter(id %in% prioritized_tbl_oi_all$id) %>%
	  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
	
	prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
	
	senders_receivers <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular')
	
	color.use = c('#A7BCD3', '#59A14F', '#B07AA1', '#E83800FF')
	names(color.use) <- senders_receivers
	
	circos_list = make_circos_group_comparison(prioritized_tbl_oi, color.use, color.use)
	
	# control
	
	prioritized_tbl_oi_skco_50 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "wound_control_8_week")
	
	circos_skco = make_circos_one_group(prioritized_tbl_oi_skco_50, color.use, color.use)
	
	# dox
	
	prioritized_tbl_oi_skdo_50 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                50, 
	                                                groups_oi = "wound_dox_8_week")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_50, color.use, color.use)
	
	circos_skco[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/Wound Skin Control Top50 Celltype Circos.svg",
	          width = 10,
	          height = 10)
	dev.off()
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/Wound Skin +Dox Top50 Celltype Circos.svg",
	          width = 10,
	          height = 10)
	dev.off()
	
	circos_skdo[[2]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/Wound Skin Top30 Circos Legend.svg",
	          width = 6,
	          height = 10)
	dev.off()
	
	group_oi = list("wound_control_8_week", "wound_dox_8_week")
	
	prioritized_tbl_oi_100 = lapply(X = group_oi,
	                                FUN = function(x){
	                                  get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 100, groups_oi = x)
	                                })
	
	plot_oi_100 = lapply(X = prioritized_tbl_oi_100,
	                     FUN = function(x){
	                       make_sample_lr_prod_activity_plots(multinichenet_celltype_output$prioritization_tables, x)
	                     })
	
	ggsave2(plot_oi_100[[1]],
	        filename = "Wound Skin Control Top100 Celltype LR Products and Ligand Activity.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 3500,
	        height =5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[2]],
	        filename = "Wound Skin Dox Top100 Celltype LR Products and Ligand Activity.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 3500,
	        height =5000,
	        units = "px",
	        limitsize = F)
	
	prioritized_tbl_oi_100 = lapply(X = group_oi,
	                                FUN = function(x){
	                                  get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 100, groups_oi = x)
	                                })
	
	plot_oi_100 = lapply(X = prioritized_tbl_oi_100,
	                     FUN = function(x){
	                       make_sample_lr_prod_activity_plots(multinichenet_celltype_output$prioritization_tables, x)
	                     })
	
	ggsave2(plot_oi_100[[1]],
	        filename = "Wound Skin Control Top100 Celltype LR Products and Ligand Activity.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 3500,
	        height =5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[2]],
	        filename = "Wound Skin Dox Top100 Celltype LR Products and Ligand Activity.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 3500,
	        height =5000,
	        units = "px",
	        limitsize = F)
	
	#visualize expression-correlated target genes of L-R pairs 
	
	receivers_oi <- unique(senders_receivers)
	
	top_n_target <- 250
	
	lr_filtered <- multinichenet_celltype_output$lr_target_prior_cor 
	
	test <- multinichenet_celltype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)
	
	lr_filtered <- lr_filtered %>% inner_join(multinichenet_celltype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) 
	
	lr_filtered <- lr_filtered %>% inner_join(contrast_tbl) 
	
	lr_filtered <- lr_filtered %>% rowwise() 
	
	receivers_oi <- as.list(senders_receivers) 
	
	names(receivers_oi) = senders_receivers
	
	#Control Skin
	skco = filter(lr_filtered, group == "wound_control_8_week")
	
	skco.up = skco %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skco.down = skco %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skconew = bind_rows(skco.up, skco.down)
	
	prioritized_tbl_oi_skco = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "wound_control_8_week", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skco = make_lr_target_correlation_plot(multinichenet_celltype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skco,  
	                                                          skconew, 
	                                                          multinichenet_celltype_output$grouping_tbl, 
	                                                          multinichenet_celltype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skco[[1]],
	        filename = "Wound Skin Control Celltype LR Target Correlation Plot.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 4000,
	        height =4000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skco[[2]],
	        filename = "Wound Skin Control Celltype LR Target Correlation Plot Legend.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	#Dox Skin
	skdo = filter(lr_filtered, group == "wound_dox_8_week")
	
	skdo.up = skdo %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skdo.down = skdo %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skdonew = bind_rows(skdo.up, skdo.down)
	
	prioritized_tbl_oi_skdo = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "wound_dox_8_week", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skdo = make_lr_target_correlation_plot(multinichenet_celltype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skdo,  
	                                                          skdonew, 
	                                                          multinichenet_celltype_output$grouping_tbl, 
	                                                          multinichenet_celltype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skdo[[1]],
	        filename = "Wound Skin Dox Celltype LR Target Correlation Plot.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 5000,
	        height =4000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skdo[[2]],
	        filename = "Wound Skin Dox Celltype LR Target Correlation Plot Legend.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	###############
	
	#L1 Subtype
	###############
	
	#Define contrasts 
	
	contrasts_oi <- c("'wound_control_8_week-wound_dox_8_week','wound_dox_8_week-wound_control_8_week'")
	
	contrast_tbl <- tibble(contrast = c("wound_control_8_week-wound_dox_8_week","wound_dox_8_week-wound_control_8_week"),
	                       group = c("wound_control_8_week","wound_dox_8_week"))
	
	#visualization of top 50 interactions across all groups
	prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 50, rank_per_group = FALSE)
	
	prioritized_tbl_oi = multinichenet_woundskin_L1subtype_output$prioritization_tables$group_prioritization_tbl %>%
	  filter(id %in% prioritized_tbl_oi_all$sampleid) %>%
	  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
	
	prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
	
	senders_receivers <- c('Keratinocyte', 'Fibroblast', 'Immune', 'Vascular', 'Melanocyte', 'Mesenchymal')
	
	color.use = c('#A7BCD3', '#B07AA1', '#59A14F', '#E15759', '#6A452B', '#E6C2A2FF')
	names(color.use) <- senders_receivers
	
	circos_list = make_circos_group_comparison(prioritized_tbl_oi, color.use, color.use)
	
	# control
	
	prioritized_tbl_oi_skco_30 = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                                30, 
	                                                groups_oi = "control")
	
	circos_skco = make_circos_one_group(prioritized_tbl_oi_skco_30, color.use, color.use)
	
	# dox
	
	prioritized_tbl_oi_skdo_30 = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                                30, 
	                                                groups_oi = "dox")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_30, color.use, color.use)
	
	circos_skco[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/Wound Skin Control Top30 L1 Subtype Circos.svg",
	          width = 10,
	          height = 10)
	dev.off()
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/Wound Skin Dox Top30 L1 Subtype Circos.svg",
	          width = 10,
	          height = 10)
	dev.off()
	
	circos_skdo[[2]]
	dev.print(svg,
	          "/data/overmilleram/scRNAseq/Wound/MultiNicheNet/Wound Skin Top30 Circos Legend.svg",
	          width = 6,
	          height = 10)
	dev.off()
	
	group_oi = list("control", "dox")
	
	prioritized_tbl_oi_100 = lapply(X = group_oi,
	                                FUN = function(x){
	                                  get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 100, groups_oi = x)
	                                })
	
	plot_oi_100 = lapply(X = prioritized_tbl_oi_100,
	                     FUN = function(x){
	                       make_sample_lr_prod_activity_plots(multinichenet_woundskin_L1subtype_output$prioritization_tables, x)
	                     })
	
	ggsave2(plot_oi_100[[1]],
	        filename = "Wound Skin Control Top100 L1 Subtype LR Products and Ligand Activity.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 3750,
	        height =5500,
	        units = "px",
	        limitsize = F)
	
	ggsave2(plot_oi_100[[2]],
	        filename = "Wound Skin Dox Top100 L1 Subtype LR Products and Ligand Activity.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 4000,
	        height =5500,
	        units = "px",
	        limitsize = F)
	
	#visualize expression-correlated target genes of L-R pairs 
	
	receivers_oi <- unique(senders_receivers)
	
	top_n_target <- 250
	
	lr_filtered <- multinichenet_woundskin_L1subtype_output$lr_target_prior_cor 
	
	lr_filtered <- lr_filtered %>% inner_join(multinichenet_woundskin_L1subtype_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) 
	
	lr_filtered <- lr_filtered %>% inner_join(contrast_tbl) 
	
	lr_filtered <- lr_filtered %>% rowwise() 
	
	receivers_oi <- as.list(senders_receivers) 
	
	names(receivers_oi) = senders_receivers
	
	#Control Skin
	skco = filter(lr_filtered, group == "control")
	
	skco.up = skco %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skco.down = skco %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skconew = bind_rows(skco.up, skco.down)
	
	prioritized_tbl_oi_skco = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "control", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skco = make_lr_target_correlation_plot(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skco,  
	                                                          skconew, 
	                                                          multinichenet_woundskin_L1subtype_output$grouping_tbl, 
	                                                          multinichenet_woundskin_L1subtype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skco[[1]],
	        filename = "Wound Skin Control Skin L1 Subtype LR Target Correlation Plot.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 4000,
	        height =4000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skco[[2]],
	        filename = "Wound Skin Control Skin L1 Subtype LR Target Correlation Plot Legend.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	#Dox Skin
	skdo = filter(lr_filtered, group == "dox")
	
	skdo.up = skdo %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
	skdo.down = skdo %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson > -0.50 | spearman > -0.50))
	skdonew = bind_rows(skdo.up, skdo.down)
	
	prioritized_tbl_oi_skdo = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                             50, 
	                                             groups_oi = "dox", 
	                                             receivers_oi = receivers_oi)
	
	lr_target_cor_plot_skdo = make_lr_target_correlation_plot(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                                          prioritized_tbl_oi_skdo,  
	                                                          skdonew, 
	                                                          multinichenet_woundskin_L1subtype_output$grouping_tbl, 
	                                                          multinichenet_woundskin_L1subtype_output$celltype_info, 
	                                                          receivers_oi,
	                                                          plot_legend = FALSE)
	
	ggsave2(lr_target_cor_plot_skdo[[1]],
	        filename = "Wound Skin +Dox Skin L1 Subtype LR Target Correlation Plot.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 5000,
	        height =4000,
	        units = "px",
	        limitsize = F)
	ggsave2(lr_target_cor_plot_skdo[[2]],
	        filename = "Wound Skin +Dox Skin L1 Subtype LR Target Correlation Plot Legend.svg",
	        path = '/data/overmilleram/scRNAseq/Wound/MultiNicheNet',
	        width = 2000,
	        height =2000,
	        units = "px",
	        limitsize = F)
	
	# comparison multinichenet
	
	prioritized_tbl_oi_wound = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables, 
	                                              10, 
	                                              rank_per_group = TRUE, 
	                                              groups_oi = "wound_dox_8_week")
	
	# create sample-level data frame for these interactions
	sample_data_wound = multinichenet_woundskin_L1subtype_output$prioritization_tables$sample_prioritization_tbl %>% 
	  dplyr::filter(id %in% prioritized_tbl_oi_wound$id) %>% 
	  dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
	                lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
	  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>% 
	  dplyr::arrange(sender, .by_group = TRUE) 
	
	sample_data_wound = sample_data_wound %>% dplyr::mutate(sender_receiver = factor(sender_receiver, 
	                                                                                 levels = sample_data_wound$sender_receiver %>% 
	                                                                                   unique()))
	
	# define the time point and group and link it all together
	
	grouping_tbl2_wound = multinichenet_woundskin_L1subtype_output$grouping_tbl %>% 
	  dplyr::inner_join(multinichenet_woundskin_L1subtype_output$prioritization_tables$sample_prioritization_tbl %>% 
	                      dplyr::distinct(sample, keep_receiver, keep_sender))
	
	grouping_tbl2_wound = grouping_tbl2_wound %>% inner_join(tibble(group = c('wound_control_8_week', 'wound_dox_8_week'), 
	                                                                contrast = c('wound_control_8_week', 'wound_dox_8_week')))
	
	grouping_tbl2_wound$condition = "control"
	
	grouping_tbl2_wound$condition[grouping_tbl2_wound$group %in% c('wound_dox_8_week')] = "dox"
	
	sample_data_wound = sample_data_wound %>% ungroup() %>% 
	  mutate(sampleid = sample_data_wound$sample) %>% 
	  inner_join(grouping_tbl2_wound)
	
	sample_data_wound = sample_data_wound %>% filter(keep_sender & keep_receiver) %>% 
	  mutate(group = factor(group, levels = c('wound_control_8_week', 'wound_dox_8_week')), 
	         condition = factor(condition, levels = c('control', 'dox')))
	
	# bring it together
	
	ex.1 <- sample_data_wound %>% filter(keep_receiver == 1 & keep_sender == 1) %>% ungroup() %>% 
	  dplyr::select(id, 
	                condition, 
	                sampleid,
	                ligand_receptor_pb_prod) 
	
	aggregate.ex.1 <- aggregate(ligand_receptor_pb_prod ~ id + condition, ex.1, mean)
	aggregate.ex.2 <- aggregate.ex.1[aggregate.ex.1$condition == 'control', ] %>% tidyr::spread(condition, ligand_receptor_pb_prod)
	
	join.ex.1 <- inner_join(ex.1, aggregate.ex.2, by = c('id'))
	
	names(join.ex.1) <- c('id', 'condition', 'sampleid', 'ligand_receptor_pb_prod', 'avg_control_lrpp')
	
	join.ex.2 <- join.ex.1 %>%  
	  mutate(diff = ligand_receptor_pb_prod-avg_control_lrpp,
	         fc = ligand_receptor_pb_prod/avg_control_lrpp) %>% 
	  mutate(lfc = log(fc)) %>% 
	  arrange(-lfc)
	
	sample_data_wound_2 <- sample_data_wound %>% inner_join(join.ex.2)
	
	order_sampleid = sample_data_wound_2 %>% 
	  group_by(sampleid) %>% 
	  summarise(sum_diff = sum(diff, 
	                           na.rm = TRUE)) %>% 
	  arrange(-sum_diff) %>% 
	  pull(sampleid)
	
	# boxplot
	
	p_lr_prod_change_boxplot_wound = sample_data_wound_2 %>% 
	  mutate(sampleid = factor(sampleid, 
	                           levels = order_sampleid)) %>% 
	  ggplot(aes(x = contrast, 
	             y = ligand_receptor_pb_prod, 
	             fill = condition,
	             group = group)) +
	  geom_boxplot() + 
	  facet_wrap(id~.) +
	  theme_bw() +
	  xlab("") + ylab("") 
	
	p_lr_prod_change_boxplot_wound
	
	# bubble plot
	
	max_diff = abs(sample_data_wound_2$diff) %>% max(na.rm = TRUE)
	
	custom_scale_color = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>%  rev(), 
	                                           values = c(0, 0.30, 0.425, 0.5, 0.575, 0.70, 1), 
	                                           limits = c(-1 * max_diff, max_diff))
	
	p_lr_prod_change_wound = sample_data_wound_2 %>% 
	  mutate(patient = factor(sampleid, 
	                          levels = order_sampleid)) %>%
	  ggplot(aes(sampleid, 
	             lr_interaction, 
	             color = diff)) +
	  geom_point(size = 5) +
	  facet_grid(sender_receiver~contrast, 
	             scales = "free",
	             space = "free", 
	             switch = "y") +
	  theme_light() +  
	  theme(axis.ticks = element_blank(), 
	        axis.title = element_blank(), 
	        axis.text.y = element_text(face = "bold.italic", 
	                                   size = 9), 
	        axis.text.x = element_text(size = 9, 
	                                   angle = 90, 
	                                   hjust = 0), 
	        panel.grid.major = element_blank(), 
	        panel.grid.minor = element_blank(), 
	        panel.spacing.x = unit(0.4, 
	                               "lines"), 
	        panel.spacing.y = unit(0.25, 
	                               "lines"), 
	        strip.text.x.top = element_text(size = 10, 
	                                        color = "black", 
	                                        face = "bold", 
	                                        angle = 0), 
	        strip.text.y.left = element_text(size = 9, 
	                                         color = "black", 
	                                         face = "bold", 
	                                         angle = 0), 
	        strip.background = element_rect(color = "darkgrey",
	                                        fill = "whitesmoke", 
	                                        size = 1.5, 
	                                        linetype = "solid")) + 
	  custom_scale_color +
	  xlab("") +
	  ylab("")
	
	p_lr_prod_change_wound
	
	
	library(easyalluvial)
	
	alluvial.table <- data.frame(paste(prioritized_tbl_oi_imm_wound$sender, prioritized_tbl_oi_imm_wound$ligand, sep = '_'), 
	                             paste(prioritized_tbl_oi_imm_wound$receiver, prioritized_tbl_oi_imm_wound$receptor, sep = '_'))
	names(alluvial.table) <- c('sender_ligand', 'receiver_receptor')
	alluvial.table <- alluvial.table[order(alluvial.table$sender_ligand), ]
	
	flow_color <- c(rep('#CD3122FF', times = 10), rep('#F6C4E5FF', times = 2))
	
	target_color <- c('#CD3122FF','#CCF4CFFF','#CCF4CFFF','#FFFF80FF','#FFD700FF',
	                  '#FFD700FF','#C4A000FF','#FF0010FF','#6088A0FF','#FFB79FFF',
	                  '#FFB79FFF','#97AD3DFF','#15CC31FF','#15CC31FF','#15CC31FF',
	                  '#15CC31FF','#15CC31FF','#415521FF','#415521FF','#E5B17EFF',
	                  '#99540FFF','#99540FFF','#CD3122FF','#F6C4E5FF','#F6C4E5FF',
	                  '#F6C4E5FF')
	
	alluvial_wide(alluvial.table, 
	              max_variables = 2,
	              fill_by = 'last_variable',
	              col_vector_flow = target_color,
	              col_vector_value = flow_color,
	              stratum_labels = T,
	              stratum_label_size = 2) + theme_void()
	
	ggsave2(filename = "Fig 6e Wound Skin Pitx1+ Top50 Neutrophil-L1subtype Alluvial Labeled.svg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "Fig 6e Wound Skin Pitx1+ Top50 Neutrophil-L1subtype Alluvial Labeled.jpeg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	alluvial_wide(alluvial.table, 
	              max_variables = 2,
	              fill_by = 'last_variable',
	              col_vector_flow = target_color,
	              col_vector_value = flow_color,
	              stratum_labels = F) + theme_void()
	
	ggsave2(filename = "Fig 6e Wound Skin Pitx1+ Top50 Neutrophil-L1subtype Alluvial Unlabeled.svg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "Fig 6e Wound Skin Pitx1+ Top50 Neutrophil-L1subtype Alluvial Unlabeled.jpeg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")

## Xenium in situ analysis

### Healthy Skin and Buccal Mucosa

	library(dplyr)
	library(tidyverse)
	library(ggplot2)
	library(Matrix)
	library(BiocManager)
	library(xlsx)
	library(readxl)
	library(paletteer)
	library(future)
	library(Giotto)
	library(ggplot2)
	library(cowplot)
	library(qs)
	
	setwd('/data/overmilleram/Xenium/Giotto/')
	
	#installGiottoEnvironment(force_environment = T) #install Giotto Python environment -- takes a few minutes
	
	results.folder = '/data/overmilleram/Xenium/Giotto/'
	
	healthy.xenium <- c('Male Skin 1', 'Male Skin 2', 'Buccal')
	
	xenium.folder = '/data/overmilleram/Xenium/Datasets/'
	
	instrs <- createGiottoInstructions(save_dir = results.folder,
	                                   save_plot = F,
	                                   show_plot = F,
	                                   return_plot = T)
	
	healthy.giotto <- list()
	
	for (i in healthy.xenium) {
	  healthy.giotto[[i]] <- createGiottoXeniumObject(xenium_dir = paste0(xenium.folder, i),
	                                                  data_to_use = 'subcellular',
	                                                  bounds_to_load = c('cell', 'nucleus'),
	                                                  qv_threshold = 20,
	                                                  h5_expression = F,
	                                                  instructions = instrs,
	                                                  cores = 8) # set number of cores to use
	}
	
	merge.healthy.giotto <- joinGiottoObjects(healthy.giotto,
	                                          gobject_names = healthy.xenium,
	                                          join_method = 'shift',
	                                          x_shift = c(0, 6000, 12000),
	                                          verbose = T)
	
	# plot centroid info
	spatPlot2D(merge.healthy.giotto,
	           spat_unit = 'cell',
	           point_shape = 'no_border',
	           point_size = 0.4,
	           point_alpha = 0.4)
	
	ggsave(filename = "Healthy Skin Giotto Skin & Oral Spatial Plot.svg",
	       path = '/data/overmilleram/Xenium/Giotto/',
	       width = 15000,
	       height = 20000,
	       units = "px",
	       limitsize = F)
	
	# generate aggregated expression based on feature & boundary (polygon) info
	merge.healthy.giotto <- calculateOverlapRaster(merge.healthy.giotto,
	                                               spatial_info = 'cell',
	                                               feat_info = 'rna')
	
	# assign polygon overlaps info to expression matrix
	merge.healthy.giotto <- overlapToMatrix(merge.healthy.giotto,
	                                        poly_info = 'cell',
	                                        feat_info = 'rna',
	                                        name = 'raw')
	
	# append feature metadata
	
	healthy.meta <- data.frame('feat_ID' = rownames(merge.healthy.giotto@expression[["cell"]][["rna"]][["raw"]]@exprMat))
	
	merge.healthy.giotto <- addFeatMetadata(gobject = merge.healthy.giotto,
	                                        feat_type = 'rna',
	                                        spat_unit = 'cell',
	                                        new_metadata = healthy.meta,
	                                        by_column = TRUE,
	                                        column_feat_ID = 'feat_ID')
	
	# data filtering
	
	filterCombinations(merge.healthy.giotto, 
	                   expression_thresholds = 1,
	                   feat_det_in_min_cells = c(10,10,10,10,10),
	                   min_det_feats_per_cell = c(1,2,3,4,5)) # chose 1,10,3
	
	merge.healthy.giotto <- filterGiotto(gobject = merge.healthy.giotto,
	                                     spat_unit = 'cell',
	                                     poly_info = 'cell',
	                                     expression_threshold = 1,
	                                     feat_det_in_min_cells = 10,
	                                     min_det_feats_per_cell = 3)
	
	# add data statistics
	merge.healthy.giotto <- addStatistics(merge.healthy.giotto, expression_values = 'raw')
	
	# normalize expression
	merge.healthy.giotto <- normalizeGiotto(gobject = merge.healthy.giotto,
	                                        spat_unit = 'cell',
	                                        scalefactor = 5000,
	                                        verbose = T)
	
	#calculate highly variable features
	merge.healthy.giotto <- calculateHVF(gobject = merge.healthy.giotto,
	                                     spat_unit = 'cell')
	
	#PCA
	merge.healthy.giotto <- runPCA(gobject = merge.healthy.giotto,
	                               spat_unit = 'cell',
	                               expression_values = 'scaled',
	                               feats_to_use = NULL,
	                               scale_unit = F,
	                               center = F)
	
	# Visualize Screeplot and PCA
	screePlot(merge.healthy.giotto, # 100 PCs accounts for ~100% of cumulative variance
	          ncp = 100)
	
	# Save Giotto object
	saveGiotto(merge.healthy.giotto,
	           dir = '/data/overmilleram/Xenium/Giotto',
	           foldername = 'Healthy Skin & Buccal Merged Giotto Object',
	           method = 'qs',
	           overwrite = T)
	
	# Load Giotto object
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	
	# UMAP
	merge.healthy.giotto <- runUMAP(merge.healthy.giotto,
	                                dimensions_to_use = 1:100,
	                                n_components = 3,
	                                spat_unit = 'cell')
	
	# sNN & Leiden clustering
	merge.healthy.giotto <- createNearestNetwork(merge.healthy.giotto,
	                                             dimensions_to_use = 1:100,
	                                             k = 30,
	                                             spat_unit = 'cell')
	
	merge.healthy.giotto <- doLeidenCluster(merge.healthy.giotto,
	                                        resolution = 2,
	                                        n_iterations = 200,
	                                        spat_unit = 'cell')
	
	# Save Giotto object
	saveGiotto(merge.healthy.giotto,
	           dir = '/data/overmilleram/Xenium/Giotto',
	           foldername = 'Healthy Skin & Buccal Merged Giotto Object',
	           overwrite = T,
	           method = 'qs')
	
	# visualize UMAP cluster results
	plotUMAP(gobject = merge.healthy.giotto,
	         spat_unit = 'cell',
	         cell_color = 'leiden_clus',
	         show_legend = FALSE,
	         point_size = 0.1,
	         point_shape = 'no_border')
	
	ggsave(filename = "Healthy Skin & Buccal Merged Giotto UMAP.svg",
	       path = '/data/overmilleram/Xenium/Giotto/',
	       width = 4000,
	       height = 4000,
	       units = "px",
	       limitsize = F)
	
	healthy.umap.3d <- plotUMAP_3D(gobject = merge.healthy.giotto, 
	                               cell_color = 'leiden_clus', 
	                               point_size = 4,
	                               show_center_label = F)
	
	htmlwidgets::saveWidget(plotly::as_widget(healthy.umap.3d), "Healthy Skin_Oral Merged Giotto 3D UMAP.html")
	
	showClusterHeatmap(gobject = merge.healthy.giotto, 
	                   cluster_column = 'leiden_clus',
	                   save_plot = T,
	                   save_param = list(save_dir = '/data/overmilleram/Xenium/Giotto',
	                                     save_name = 'Healthy Skin & Buccal Merged Giotto Cluster Comparison Heatmap',
	                                     save_format = 'svg'))
	
	# Find markers
	
	healthy.markers = findMarkers_one_vs_all(gobject = merge.healthy.giotto,
	                                         method = 'gini',
	                                         expression_values = 'normalized',
	                                         cluster_column = 'leiden_clus',
	                                         rank_score = 2)
	
	#Make Excel sheet with markers and cell numbers
	
	write.xlsx(healthy.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	# load in Giotto objects
	
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	#export cluster annotations for Xenium Explorer
	clustering <- merge.healthy.giotto@cell_metadata$cell$rna@metaDT %>% as.data.frame
	
	clustering[ ,2:5] <- NULL
	
	clustering.ms1 <- clustering[1:177796, ]
	clustering.ms2 <- clustering[177797:265190, ]
	clustering.buc <- clustering[265191:331586, ]
	
	clustering.ms1[ ,1] <- gsub('Male Skin 1-', '', clustering.ms1[ ,1])
	clustering.ms2[ ,1] <- gsub('Male Skin 2-', '', clustering.ms2[ ,1])
	clustering.buc[ ,1] <- gsub('Buccal-', '', clustering.buc[ ,1])
	
	colnames(clustering.ms1) <- c('cell_id', 'group')
	colnames(clustering.ms2) <- c('cell_id', 'group')
	colnames(clustering.buc) <- c('cell_id', 'group')
	
	clustering.ms1 <- clustering.ms1[order(clustering.ms1$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	clustering.ms2 <- clustering.ms2[order(clustering.ms2$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	clustering.buc <- clustering.buc[order(clustering.buc$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	
	
	write.csv(clustering.ms1,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Oral Merged Giotto Male Skin 1 Cell ID-Group Updated for Explorer.csv',
	          row.names = F)
	
	write.csv(clustering.ms2,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Oral Merged Giotto Male Skin 2 Cell ID-Group Updated for Explorer.csv',
	          row.names = F)
	
	write.csv(clustering.buc,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Oral Merged Giotto Buccal Cell ID-Group for Updated Explorer.csv',
	          row.names = F)
	
	violinPlot(merge.healthy.giotto, 
	           feats = c('Krt5', 'Ccl27a', 'Col1a1', 'Saa3', 'Cd74', 'Foxp3', 'Slpi', 'Pecam1', 'Lyve1', 'Vwf', 'Mpz', 'Ncmap'), 
	           cluster_column = 'leiden_clus',
	           strip_position = 'right',
	           save_plot = F)
	
	# add celltype & subtype metadata
	
	tissue.annotations <- read_xlsx('/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Updated Cluster Markers.xlsx',
	                                sheet = 1,
	                                range = 'M1:O64',
	                                col_names = T)
	
	data.files = list.files(path = '/data/overmilleram/Xenium/Cell Coordinates/Healthy Skin/',
	                        pattern = "*.xlsx")
	
	data <- lapply(data.files, 
	               function(x) read_xlsx(paste0('/data/overmilleram/Xenium/Cell Coordinates/Healthy Skin/', x), 
	                                     sheet = 1,
	                                     range = cell_rows(c(4,NA)),
	                                     col_names = F))
	
	file.names <- gsub('_cell_stats.xlsx', '', data.files) # extract ending sequence
	file.names <- gsub(' ', '_', file.names) # replace spaces with underscores
	
	names(data) <- file.names
	
	#tissue.cells <- lapply(data, 
	#                       function(x) x$...1 %>% data.frame %>% set_names('cell_ID')) # extract only cell_ID's into a new list
	#
	#for(i in 1:14){
	#  tissue.cells[[i]] <- cbind(tissue.cells[[i]], sample = names(tissue.cells[i]))
	#} # add a new column including the name of the sample
	#
	#tissue.cells <- reduce(tissue.cells, full_join) # full_join all listed df's into one df for addition to metadata
	
	tissue.meta <- pDataDT(merge.healthy.giotto) # pull existing metadata
	
	#tissue.meta$cell_ID <- stringi::stri_replace_all_regex(tissue.meta$cell_ID,
	#                                                       pattern = c('Male Skin 1-', 'Male Skin 2-', 'Buccal-'),
	#                                                       replacement = '',
	#                                                       vectorize = F)
	
	tissue.meta <- left_join(tissue.meta, tissue.annotations, join_by(leiden_clus == Cluster)) %>% select(cell_ID, Celltype, L1Subtype)
	
	#tissue.meta <- inner_join(tissue.meta, tissue.cells, by = 'cell_ID') 
	
	tissue.meta <- column_to_rownames(tissue.meta, 'cell_ID')
	
	merge.healthy.giotto <- addCellMetadata(merge.healthy.giotto2,
	                                        spat_unit = 'cell',
	                                        feat_type = 'rna',
	                                        tissue.meta,
	                                        by_column = F)
	
	head(pDataDT(merge.healthy.giotto), 10) # check if metadata correct
	
	# Find celltype markers & make Excel sheet
	
	healthy.markers = findMarkers_one_vs_all(gobject = merge.healthy.giotto,
	                                         method = 'gini',
	                                         expression_values = 'normalized',
	                                         cluster_column = 'Celltype',
	                                         rank_score = 2)
	
	write.xlsx(healthy.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Celltype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	# plot top cell markers in heatmap
	
	topgini_genes = unique(healthy.markers[, head(.SD, 5), by = 'cluster']$feats)
	
	plotMetaDataHeatmap(merge.healthy.giotto, 
	                    custom_cluster_order = c('Fibroblast', 'Immune', 'Vascular', 'Neural', 'Mesenchymal', 'Keratinocyte', 'Salivary'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'Celltype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = topgini_genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 3b Xenium Celltype Marker Heatmap.svg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 3b Xenium Celltype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	# Find L1Subtype markers & make Excel sheet
	
	healthy.markers = findMarkers_one_vs_all(gobject = merge.healthy.giotto,
	                                         method = 'gini',
	                                         expression_values = 'normalized',
	                                         cluster_column = 'L1Subtype',
	                                         rank_score = 2)
	
	write.xlsx(healthy.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto L1Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	# plot top cell markers in heatmap
	
	topgini_genes = unique(healthy.markers[, head(.SD, 5), by = 'cluster']$feats)
	
	plotMetaDataHeatmap(merge.healthy.giotto, 
	                    #custom_cluster_order = c('Fibroblast', 'Immune', 'Vascular', 'Neural', 'Mesenchymal', 'Keratinocyte', 'Salivary'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = topgini_genes) + theme_minimal()
	
	#ggsave2(filename = "SupFig 3b Xenium L1Subtype Marker Heatmap.svg",
	#        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	#        width = 5000,
	#        height = 5000,
	#        units = "px")
	#ggsave2(filename = "SupFig 3b Xenium L1Subtype Marker Heatmap.jpeg",
	#        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	#        width = 5,
	#        height = 5,
	#        dpi = 600,
	#        units = "in")
	
	# Save Giotto object
	saveGiotto(merge.healthy.giotto,
	           dir = '/data/overmilleram/Xenium/Giotto',
	           foldername = 'Healthy Skin & Buccal Merged Giotto Object',
	           overwrite = T,
	           method = 'RDS')
	
	# export additional cluster annotations for Xenium Explorer
	clustering <- merge.healthy.giotto@cell_metadata$cell$rna@metaDT %>% as.data.frame
	clustering.celltype <- clustering[ ,c(1:2,7)]
	clustering.L1subtype <- clustering[ ,c(1:2,8)]
	
	clustering.celltype.ms1 <- subset(clustering.celltype, list_ID == 'Male Skin 1')
	clustering.L1subtype.ms1 <- subset(clustering.L1subtype, list_ID == 'Male Skin 1')
	clustering.celltype.ms2 <- subset(clustering.celltype, list_ID == 'Male Skin 2')
	clustering.L1subtype.ms2 <- subset(clustering.L1subtype, list_ID == 'Male Skin 2')
	clustering.celltype.buc <- subset(clustering.celltype, list_ID == 'Buccal')
	clustering.L1subtype.buc <- subset(clustering.L1subtype, list_ID == 'Buccal')
	
	clustering.celltype.ms1[ ,2] <- NULL
	clustering.L1subtype.ms1[ ,2] <- NULL
	clustering.celltype.ms2[ ,2] <- NULL
	clustering.L1subtype.ms2[ ,2] <- NULL
	clustering.celltype.buc[ ,2] <- NULL
	clustering.L1subtype.buc[ ,2] <- NULL
	
	clustering.celltype.ms1[ ,1] <- gsub('Male Skin 1-', '', clustering.celltype.ms1[ ,1])
	clustering.L1subtype.ms1[ ,1] <- gsub('Male Skin 1-', '', clustering.L1subtype.ms1[ ,1])
	clustering.celltype.ms2[ ,1] <- gsub('Male Skin 2-', '', clustering.celltype.ms2[ ,1])
	clustering.L1subtype.ms2[ ,1] <- gsub('Male Skin 2-', '', clustering.L1subtype.ms2[ ,1])
	clustering.celltype.buc[ ,1] <- gsub('Buccal-', '', clustering.celltype.buc[ ,1])
	clustering.L1subtype.buc[ ,1] <- gsub('Buccal-', '', clustering.L1subtype.buc[ ,1])
	
	colnames(clustering.celltype.ms1) <- c('cell_id', 'group')
	colnames(clustering.L1subtype.ms1) <- c('cell_id', 'group')
	colnames(clustering.celltype.ms2) <- c('cell_id', 'group')
	colnames(clustering.L1subtype.ms2) <- c('cell_id', 'group')
	colnames(clustering.celltype.buc) <- c('cell_id', 'group')
	colnames(clustering.L1subtype.buc) <- c('cell_id', 'group')
	
	clustering.celltype.ms1 <- clustering.celltype.ms1[order(clustering.celltype.ms1$group), ] #order column alphabetically by 'group' so that it makes sense in Explorer
	clustering.L1subtype.ms1 <- clustering.L1subtype.ms1[order(clustering.L1subtype.ms1$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	clustering.celltype.ms2 <- clustering.celltype.ms2[order(clustering.celltype.ms2$group), ] #order column alphabetically by 'group' so that it makes sense in Explorer
	clustering.L1subtype.ms2 <- clustering.L1subtype.ms2[order(clustering.L1subtype.ms2$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	clustering.celltype.buc <- clustering.celltype.buc[order(clustering.celltype.buc$group), ] #order column alphabetically by 'group' so that it makes sense in Explorer
	clustering.L1subtype.buc <- clustering.L1subtype.buc[order(clustering.L1subtype.buc$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	
	write.csv(clustering.celltype.ms1,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Buccal Merged Giotto Male Skin 1 Cell ID-Celltype Updated for Explorer.csv',
	          row.names = F)
	write.csv(clustering.L1subtype.ms1,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Buccal Merged Giotto Male Skin 1 Cell ID-L1Subtype Updated for Explorer.csv',
	          row.names = F)
	write.csv(clustering.celltype.ms2,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Buccal Merged Giotto Male Skin 2 Cell ID-Celltype Updated for Explorer.csv',
	          row.names = F)
	write.csv(clustering.L1subtype.ms2,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Buccal Merged Giotto Male Skin 2 Cell ID-L2Subtype Updated for Explorer.csv',
	          row.names = F)
	write.csv(clustering.celltype.buc,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Buccal Merged Giotto Buccal Cell ID-Celltype Updated for Explorer.csv',
	          row.names = F)
	write.csv(clustering.L1subtype.buc,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Healthy Skin & Buccal Merged Giotto Buccal Cell ID-L1Subtype Updated for Explorer.csv',
	          row.names = F)

### Wound Skin

	library(dplyr)
	library(tidyverse)
	library(ggplot2)
	library(Matrix)
	library(BiocManager)
	library(xlsx)
	library(readxl)
	library(paletteer)
	library(future)
	library(Giotto)
	library(ggplot2)
	library(cowplot)
	library(qs)
	
	setwd('/data/overmilleram/Xenium/Giotto/')
	
	#installGiottoEnvironment(force_environment = T) #install Giotto Python environment -- takes a few minutes
	
	results.folder = '/data/overmilleram/Xenium/Giotto/'
	
	instrs <- createGiottoInstructions(save_dir = results.folder,
	                                   save_plot = F,
	                                   show_plot = F,
	                                   return_plot = T)
	
	# Read in and save cell barcodes to be analyzed for wounds (remove eschar)
	
	data.files = list.files(path = '/data/overmilleram/Xenium/Cell Coordinates/Wound Skin/',
	                        pattern = "*.xlsx")
	
	data <- lapply(data.files, 
	               function(x) read_xlsx(paste0('/data/overmilleram/Xenium/Cell Coordinates/Wound Skin/', x), 
	                                     sheet = 1,
	                                     range = cell_rows(c(4,NA)),
	                                     col_names = F))
	
	file.names <- gsub('_cells_stats.xlsx', '', data.files)
	
	names(data) <- file.names
	
	#wound.fc1 <- merge(data[[1]][ ,1], data[[2]][ ,1], all = T) %>% pull
	#wound.fc2 <- merge(data[[3]][ ,1], data[[4]][ ,1], all = T) %>% pull
	wound.mc1 <- merge(data[[11]][ ,1], data[[12]][ ,1], all = T ) %>% pull
	
	#wound.fd1 <- merge(data[[5]][ ,1], data[[6]][ ,1], all = T) %>% pull
	#wound.fd2.a <- merge(data[[7]][ ,1], data[[8]][ ,1], all = T) 
	#wound.fd2.b <- merge(data[[9]][ ,1], data[[10]][ ,1], all = T) 
	#wound.fd2 <- merge(wound.fd2.a, wound.fd2.b, all = T) %>% pull
	wound.md1.a <- merge(data[[13]][ ,1], data[[14]][ ,1], all = T) 
	wound.md1.b <- merge(data[[15]][ ,1], data[[16]][ ,1], all = T)
	wound.md1 <- merge(wound.md1.a, wound.md1.b, all = T) %>% pull
	wound.md2.a <- merge(data[[17]][ ,1], data[[18]][ ,1], all = T)
	wound.md2.b <- merge(data[[19]][ ,1], data[[20]][ ,1], all = T)
	wound.md2 <- merge(wound.md2.a, wound.md2.b, all = T) %>% pull
	
	rm(wound.fd2.a, wound.fd2.b, wound.md1.a, wound.md1.b, wound.md2.a, wound.md2.b) # cleanup
	
	cells.male.wound <- c(wound.mc1, wound.md1, wound.md2)
	#cells.female.wound <- c(wound.fc1, wound.fc2, wound.fd1, wound.fd2)
	
	# load data into Giotto object
	
	wound.xenium <- c('Male Wound 1')
	
	xenium.folder = '/data/overmilleram/Xenium/Datasets/'
	
	wound.giotto <- createGiottoXeniumObject(xenium_dir = '/data/overmilleram/Xenium/Datasets/Male Wound 1/',
	                                         data_to_use = 'subcellular',
	                                         bounds_to_load = c('cell', 'nucleus'),
	                                         qv_threshold = 20,
	                                         h5_expression = F,
	                                         instructions = instrs,
	                                         cores = 8) # set number of cores to use
	
	#for (i in wound.xenium) {
	#  wound.giotto[[i]] <- createGiottoXeniumObject(xenium_dir = paste0(xenium.folder, i),
	#                                                data_to_use = 'subcellular',
	#                                                bounds_to_load = c('cell', 'nucleus'),
	#                                                qv_threshold = 20,
	#                                                h5_expression = F,
	#                                                instructions = instrs,
	#                                                cores = 8) # set number of cores to use
	#}
	
	#merge.healthy.giotto <- joinGiottoObjects(healthy.giotto,
	#                                          gobject_names = healthy.xenium,
	#                                          join_method = 'shift',
	#                                          x_shift = c(0, 6000, 12000),
	#                                          verbose = T)
	
	wound.giotto <- subsetGiotto(wound.giotto,
	                             cell_ids = cells.male.wound)
	
	# plot centroid info
	spatPlot2D(wound.giotto,
	           spat_unit = 'cell',
	           point_shape = 'no_border',
	           point_size = 0.4,
	           point_alpha = 0.4)
	
	ggsave(filename = "Wound Skin Giotto Spatial Plot.svg",
	       path = '/data/overmilleram/Xenium/Giotto/',
	       width = 15000,
	       height = 20000,
	       units = "px",
	       limitsize = F)
	
	# generate aggregated expression based on feature & boundary (polygon) info
	wound.giotto <- calculateOverlapRaster(wound.giotto,
	                                       spatial_info = 'cell',
	                                       feat_info = 'rna')
	
	# assign polygon overlaps info to expression matrix
	wound.giotto <- overlapToMatrix(wound.giotto,
	                                poly_info = 'cell',
	                                feat_info = 'rna',
	                                name = 'raw')
	
	# append feature metadata
	
	wound.meta <- data.frame('feat_ID' = rownames(wound.giotto@expression[["cell"]][["rna"]][["raw"]]@exprMat))
	
	wound.giotto <- addFeatMetadata(gobject = wound.giotto,
	                                feat_type = 'rna',
	                                spat_unit = 'cell',
	                                new_metadata = wound.meta,
	                                by_column = TRUE,
	                                column_feat_ID = 'feat_ID')
	
	# data filtering
	
	filterCombinations(wound.giotto, 
	                   expression_thresholds = 1,
	                   feat_det_in_min_cells = c(10,10,10,10,10),
	                   min_det_feats_per_cell = c(1,2,3,4,5)) # chose 1,10,3
	
	wound.giotto <- filterGiotto(gobject = wound.giotto,
	                             spat_unit = 'cell',
	                             poly_info = 'cell',
	                             expression_threshold = 1,
	                             feat_det_in_min_cells = 10,
	                             min_det_feats_per_cell = 3)
	
	# add data statistics
	wound.giotto <- addStatistics(wound.giotto, expression_values = 'raw')
	
	# normalize expression
	wound.giotto <- normalizeGiotto(gobject = wound.giotto,
	                                spat_unit = 'cell',
	                                scalefactor = 5000,
	                                verbose = T)
	
	#calculate highly variable features
	wound.giotto <- calculateHVF(gobject = wound.giotto,
	                             spat_unit = 'cell')
	
	#PCA
	wound.giotto <- runPCA(gobject = wound.giotto,
	                       spat_unit = 'cell',
	                       expression_values = 'scaled',
	                       feats_to_use = NULL,
	                       scale_unit = F,
	                       center = F)
	
	# Visualize Screeplot and PCA
	screePlot(wound.giotto, # 100 PCs accounts for ~100% of cumulative variance
	          ncp = 100)
	
	# Save Giotto object
	saveGiotto(wound.giotto,
	           dir = '/data/overmilleram/Xenium/Giotto',
	           foldername = 'Wound Skin Giotto Object',
	           method = 'qs',
	           overwrite = T)
	
	# Load Giotto object
	wound.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Wound Skin Merged Giotto Object/')
	
	
	# UMAP
	wound.giotto <- runUMAP(wound.giotto,
	                         dimensions_to_use = 1:100,
	                         n_components = 3,
	                         spat_unit = 'cell')
	
	# sNN & Leiden clustering
	wound.giotto <- createNearestNetwork(wound.giotto,
	                                     dimensions_to_use = 1:100,
	                                     k = 30,
	                                     spat_unit = 'cell')
	
	wound.giotto <- doLeidenCluster(wound.giotto,
	                                resolution = 2,
	                                n_iterations = 200,
	                                spat_unit = 'cell')
	
	# Save Giotto object
	saveGiotto(wound.giotto,
	           dir = '/data/overmilleram/Xenium/Giotto',
	           foldername = 'Wound Skin Giotto Object',
	           overwrite = T,
	           method = 'qs')
	
	# visualize UMAP cluster results
	plotUMAP(gobject = wound.giotto,
	         spat_unit = 'cell',
	         cell_color = 'leiden_clus',
	         show_legend = FALSE,
	         point_size = 0.1,
	         point_shape = 'no_border')
	
	ggsave(filename = "Wound Skin Merged Giotto UMAP.svg",
	       path = '/data/overmilleram/Xenium/Giotto/',
	       width = 4000,
	       height = 4000,
	       units = "px",
	       limitsize = F)
	
	wound.umap.3d <- plotUMAP_3D(gobject = wound.giotto, 
	                               cell_color = 'leiden_clus', 
	                               point_size = 4,
	                               show_center_label = F)
	
	htmlwidgets::saveWidget(plotly::as_widget(wound.umap.3d), "Wound Skin Giotto 3D UMAP.html")
	
	showClusterHeatmap(gobject = wound.giotto, 
	                   cluster_column = 'leiden_clus',
	                   save_plot = T,
	                   save_param = list(save_dir = '/data/overmilleram/Xenium/Giotto',
	                                     save_name = 'Wound Skin Giotto Cluster Comparison Heatmap',
	                                     save_format = 'svg'))
	
	# Find markers
	
	wound.markers = findMarkers_one_vs_all(gobject = wound.giotto,
	                                         method = 'gini',
	                                         expression_values = 'normalized',
	                                         cluster_column = 'leiden_clus',
	                                         rank_score = 2)
	
	#Make Excel sheet with markers and cell numbers
	write.xlsx(wound.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Updated Cluster Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	# load in Giotto objects
	
	wound.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Object/')
	
	# export cluster annotations for Xenium Explorer
	clustering <- wound.giotto@cell_metadata$cell$rna@metaDT %>% as.data.frame
	
	clustering[ ,2:4] <- NULL
	
	colnames(clustering) <- c('cell_id', 'group')
	
	clustering <- clustering[order(clustering$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	
	write.csv(clustering,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Wound Skin Giotto Male Wound 1 Cell ID-Group Updated for Explorer.csv',
	          row.names = F)
	
	# add celltype metadata into Giotto object
	
	wound.meta <- pDataDT(wound.giotto)
	
	# add celltype & subtype metadata
	
	wound.annotations <- read_xlsx('/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Updated Cluster Markers.xlsx',
	                                sheet = 1,
	                                range = 'M1:O44',
	                                col_names = T)
	
	#data.files = list.files(path = '/data/overmilleram/Xenium/Cell Coordinates/Wound Skin/',
	#                        pattern = "*.xlsx")
	#
	#data <- lapply(data.files, 
	#               function(x) read_xlsx(paste0('/data/overmilleram/Xenium/Cell Coordinates/Healthy Skin/', x), 
	#                                     sheet = 1,
	#                                     range = cell_rows(c(4,NA)),
	#                                     col_names = F))
	#
	#file.names <- gsub('_cell_stats.xlsx', '', data.files) # extract ending sequence
	#file.names <- gsub(' ', '_', file.names) # replace spaces with underscores
	#
	#names(data) <- file.names
	#
	#wound.cells <- lapply(data, 
	#                       function(x) x$...1 %>% data.frame %>% set_names('cell_ID')) # extract only cell_ID's into a new list
	#
	#for(i in 1:14){
	#  wound.cells[[i]] <- cbind(wound.cells[[i]], sample = names(wound.cells[i]))
	#} # add a new column including the name of the sample
	#
	#wound.cells <- reduce(wound.cells, full_join) # full_join all listed df's into one df for addition to metadata
	
	wound.meta <- left_join(wound.meta, wound.annotations, join_by(leiden_clus == Cluster)) %>% select(cell_ID, Celltype, L1Subtype)
	
	#wound.meta <- inner_join(wound.meta, wound.cells, by = 'cell_ID') 
	
	wound.meta <- column_to_rownames(wound.meta, 'cell_ID')
	
	wound.giotto <- addCellMetadata(wound.giotto,
	                                spat_unit = 'cell',
	                                feat_type = 'rna',
	                                wound.meta,
	                                by_column = F)
	
	head(pDataDT(wound.giotto), 10) # check if metadata correct
	
	# Find celltype markers & make Excel sheet
	
	wound.markers = findMarkers_one_vs_all(gobject = wound.giotto,
	                                       method = 'gini',
	                                       expression_values = 'normalized',
	                                       cluster_column = 'Celltype',
	                                       rank_score = 2)
	
	write.xlsx(wound.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Celltype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	# plot top cell markers in heatmap
	
	topgini_genes = unique(wound.markers[, head(.SD, 5), by = 'cluster']$feats)
	
	plotMetaDataHeatmap(wound.giotto, 
	                    custom_cluster_order = c('Mesenchymal', 'Vascular', 'Immune', 'Fibroblast', 'Neural', 'Keratinocyte'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'Celltype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = topgini_genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 9f Wound Skin Xenium Celltype Marker Heatmap.svg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 9f Xenium Wound Skin Celltype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	# Find L1subtype markers & make Excel sheet
	
	wound.markers = findMarkers_one_vs_all(gobject = wound.giotto,
	                                       method = 'gini',
	                                       expression_values = 'normalized',
	                                       cluster_column = 'L1Subtype',
	                                       rank_score = 2)
	
	write.xlsx(wound.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto L1Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	# plot top cell markers in heatmap
	
	topgini_genes = unique(wound.markers[, head(.SD, 5), by = 'cluster']$feats)
	
	plotMetaDataHeatmap(wound.giotto, 
	                    #custom_cluster_order = c('Mesenchymal', 'Vascular', 'Immune', 'Fibroblast', 'Neural', 'Keratinocyte'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = topgini_genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 9f Wound Skin Xenium L1Subtype Marker Heatmap.svg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 9f Xenium Wound Skin L1Subtype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/scRNAseq/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	# Save Giotto object
	saveGiotto(wound.giotto,
	           dir = '/data/overmilleram/Xenium/Giotto',
	           foldername = 'Wound Skin Giotto Object',
	           overwrite = T,
	           method = 'RDS')
	
	# export additional cluster annotations for Xenium Explorer
	clustering <- wound.giotto@cell_metadata$cell$rna@metaDT %>% as.data.frame
	clustering.celltype <- clustering
	clustering.L1subtype <- clustering
	
	clustering.celltype[ ,c(2:5,7)] <- NULL
	clustering.L1subtype[ ,c(2:5,6)] <- NULL
	
	colnames(clustering.celltype) <- c('cell_id', 'group')
	colnames(clustering.L1subtype) <- c('cell_id', 'group')
	
	clustering.celltype <- clustering.celltype[order(clustering.celltype$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	clustering.L1subtype <- clustering.L1subtype[order(clustering.L1subtype$group), ] #order column numerically by 'group' so that it makes sense in Explorer
	
	write.csv(clustering.celltype,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Wound Skin Giotto Male Wound 1 Cell ID-Celltype Updated for Explorer.csv',
	          row.names = F)
	write.csv(clustering.L1subtype,
	          file = '/data/overmilleram/Xenium/Explorer Cluster Masks/Wound Skin Giotto Male Wound 1 Cell ID-L1Subtype Updated for Explorer.csv',
	          row.names = F)

## Manuscript Figures

	library(dplyr)
	library(ggplot2)
	library(Matrix)
	library(patchwork)
	library(gdata)
	library(Seurat)
	library(SeuratObject)
	library(SeuratWrappers)
	library(BiocManager)
	library(metap)
	library(cowplot)
	library(sctransform)
	library(xlsx)
	library(glmGamPoi)
	library(clustree)
	library(biomaRt)
	library(monocle3)
	library(magrittr)
	library(slingshot)
	library(paletteer)
	library(nichenetr)
	library(multinichenetr)
	library(tidyverse)
	library(circlize)
	library(scales)
	library(kableExtra)
	library(knitr)
	library(SoupX)
	library(scDblFinder)
	library(viridis)
	library(uwot)
	library(ComplexHeatmap)
	library(enrichR)
	library(Giotto)
	library(qs)
	
	#Set current working directory
	setwd("/data/overmilleram/Manuscript Figures/")
	
	tissue.integrated <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.harmony.integrated.rds')
	wound.integrated <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.harmony.integrated.rds')
	
	#Fig1
	#######
	
	# 1c
	color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E83800FF', '#D088E0FF', '#FFAE34', '#805840FF', '#808898FF')
	names(color.use) <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular', 'Neural', 'Salivary', 'Melanocyte', 'Skeletal_Muscle')
	
	Idents(tissue.integrated) <- 'celltype'
	
	DimPlot(tissue.integrated,
	        pt.size = 0.25,
	        label = F,
	        cols = color.use,
	        split.by = 'state',
	        ncol = 3,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 1c Celltype UMAP Split State.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(filename = "Fig 1c Celltype UMAP Split State.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 1200,
	        units = "in",
	        limitsize = F)
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	tissue.integrated.test <- sc_utils(tissue.integrated)
	
	prop_test <- permutation_test(tissue.integrated.test, 
	                              cluster_identity = "celltype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "control_buccal_mucosa",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5))
	
	prop_test <- permutation_test(tissue.integrated.test, 
	                              cluster_identity = "celltype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "skin_dox_8_week",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5))
	
	table(tissue.integrated$state, tissue.integrated$celltype) #completed graph in Excel
	
	# 1d
	
	celltypes <- c('keratinocyte', 'immune', 'fibroblast', 'vascular', 'neural', 'mesenchymal', 'salivary')
	xenium.control <- c(16772, 11546, 11783, 3488, 976, 5784, 67) 
	xenium.pitx1 <- c(92803, 46792, 53052, 11320, 1721, 8847, 232)
	xenium.buccal <- c(20800, 6526, 10878, 10976, 2607, 11162, 3440)
	
	xenium.prop.cvp <- rbind(xenium.control, xenium.pitx1)
	xenium.prop.cvb <- rbind(xenium.control, xenium.buccal)
	
	colnames(xenium.prop.cvp) <- celltypes
	colnames(xenium.prop.cvb) <- celltypes
	
	#control vs pitx1+
	#p-value (all are >>>0.05)
	prop.test(c(xenium.prop.cvp[1,1], xenium.prop.cvp[2,1]), c(sum(xenium.prop.cvp[ ,1]), sum(xenium.prop.cvp[ ,1])), cor = F, conf.level = 0.975)$p.value #keratinocyte
	prop.test(c(xenium.prop.cvp[1,2], xenium.prop.cvp[2,2]), c(sum(xenium.prop.cvp[ ,2]), sum(xenium.prop.cvp[ ,2])), cor = F, conf.level = 0.975)$p.value #immune
	prop.test(c(xenium.prop.cvp[1,3], xenium.prop.cvp[2,3]), c(sum(xenium.prop.cvp[ ,3]), sum(xenium.prop.cvp[ ,3])), cor = F, conf.level = 0.975)$p.value #fibroblast
	prop.test(c(xenium.prop.cvp[1,4], xenium.prop.cvp[2,4]), c(sum(xenium.prop.cvp[ ,4]), sum(xenium.prop.cvp[ ,4])), cor = F, conf.level = 0.975)$p.value #vascular
	prop.test(c(xenium.prop.cvp[1,5], xenium.prop.cvp[2,5]), c(sum(xenium.prop.cvp[ ,5]), sum(xenium.prop.cvp[ ,5])), cor = F, conf.level = 0.975)$p.value #neural
	prop.test(c(xenium.prop.cvp[1,6], xenium.prop.cvp[2,6]), c(sum(xenium.prop.cvp[ ,6]), sum(xenium.prop.cvp[ ,6])), cor = F, conf.level = 0.975)$p.value #mesenchymal
	prop.test(c(xenium.prop.cvp[1,7], xenium.prop.cvp[2,7]), c(sum(xenium.prop.cvp[ ,7]), sum(xenium.prop.cvp[ ,7])), cor = F, conf.level = 0.975)$p.value #salivary
	
	#log2fc
	log2((xenium.prop.cvp[2,1]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,1]/sum(xenium.prop.cvp[1,]))) #keratinocyte, 0.3773
	log2((xenium.prop.cvp[2,2]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,2]/sum(xenium.prop.cvp[1,]))) #immune,      -0.0720
	log2((xenium.prop.cvp[2,3]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,3]/sum(xenium.prop.cvp[1,]))) #fibroblast,   0.0799
	log2((xenium.prop.cvp[2,4]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,4]/sum(xenium.prop.cvp[1,]))) #vascular,    -0.3924
	log2((xenium.prop.cvp[2,5]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,5]/sum(xenium.prop.cvp[1,]))) #neural,      -1.2725
	log2((xenium.prop.cvp[2,6]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,6]/sum(xenium.prop.cvp[1,]))) #mesenchymal, -1.4777
	log2((xenium.prop.cvp[2,7]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,7]/sum(xenium.prop.cvp[1,]))) #salivary,    -0.2990
	
	#control vs buccal
	#p-value (all are >>>0.05)
	prop.test(c(xenium.prop.cvb[1,1], xenium.prop.cvb[2,1]), c(sum(xenium.prop.cvb[ ,1]), sum(xenium.prop.cvb[ ,1])), cor = F, conf.level = 0.975)$p.value #keratinocyte
	prop.test(c(xenium.prop.cvb[1,2], xenium.prop.cvb[2,2]), c(sum(xenium.prop.cvb[ ,2]), sum(xenium.prop.cvb[ ,2])), cor = F, conf.level = 0.975)$p.value #immune
	prop.test(c(xenium.prop.cvb[1,3], xenium.prop.cvb[2,3]), c(sum(xenium.prop.cvb[ ,3]), sum(xenium.prop.cvb[ ,3])), cor = F, conf.level = 0.975)$p.value #fibroblast
	prop.test(c(xenium.prop.cvb[1,4], xenium.prop.cvb[2,4]), c(sum(xenium.prop.cvb[ ,4]), sum(xenium.prop.cvb[ ,4])), cor = F, conf.level = 0.975)$p.value #vascular
	prop.test(c(xenium.prop.cvb[1,5], xenium.prop.cvb[2,5]), c(sum(xenium.prop.cvb[ ,5]), sum(xenium.prop.cvb[ ,5])), cor = F, conf.level = 0.975)$p.value #neural
	prop.test(c(xenium.prop.cvb[1,6], xenium.prop.cvb[2,6]), c(sum(xenium.prop.cvb[ ,6]), sum(xenium.prop.cvb[ ,6])), cor = F, conf.level = 0.975)$p.value #mesenchymal
	prop.test(c(xenium.prop.cvb[1,7], xenium.prop.cvb[2,7]), c(sum(xenium.prop.cvb[ ,7]), sum(xenium.prop.cvb[ ,7])), cor = F, conf.level = 0.975)$p.value #salivary
	
	#log2fc
	log2((xenium.prop.cvb[2,1]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,1]/sum(xenium.prop.cvb[1,]))) #keratinocyte,-0.0865
	log2((xenium.prop.cvb[2,2]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,2]/sum(xenium.prop.cvb[1,]))) #immune,      -1.2202
	log2((xenium.prop.cvb[2,3]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,3]/sum(xenium.prop.cvb[1,]))) #fibroblast,  -0.5124
	log2((xenium.prop.cvb[2,4]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,4]/sum(xenium.prop.cvb[1,]))) #vascular,     1.2568
	log2((xenium.prop.cvb[2,5]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,5]/sum(xenium.prop.cvb[1,]))) #neural,       1.0204
	log2((xenium.prop.cvb[2,6]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,6]/sum(xenium.prop.cvb[1,]))) #mesenchymal,  0.5514
	log2((xenium.prop.cvb[2,7]/sum(xenium.prop.cvb[2,]))/(xenium.prop.cvb[1,7]/sum(xenium.prop.cvb[1,]))) #salivary,     5.2860
	
	#######
	
	#Supplementary Fig 2
	#######
	# Sup2a
	Idents(tissue.integrated) <- 'sampleid'
	
	VlnPlot(tissue.integrated,
	        pt.size = 0,
	        features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
	        split.by = 'sampleid',
	        cols = paletteer::paletteer_d("ggsci::default_ucscgb"),
	        stack = T,
	        flip = T) + NoLegend()
	
	ggsave2(filename = "SupFig 2a Quality Control VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 4000,
	        height = 4000,
	        units = "px")
	
	ggsave2(filename = "SupFig 2a Quality Control VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	# Sup2b
	
	celltype.cells.by.sample <- data.frame(table(tissue.integrated$sampleid, tissue.integrated$celltype))
	
	colnames(celltype.cells.by.sample) <- c('sampleid', 'celltype', 'cells')
	
	write.xlsx(celltype.cells.by.sample,
	           '/data/overmilleram/Manuscript Figures/Sup Fig 2b Cell Counts.xlsx')
	
	# Sup2c
	Idents(tissue.integrated) <- 'celltype'
	
	color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E83800FF', '#D088E0FF', '#FFAE34', '#805840FF', '#808898FF')
	names(color.use) <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular', 'Neural', 'Salivary', 'Melanocyte', 'Skeletal_Muscle')
	
	DimPlot(tissue.integrated,
	        pt.size = 0.5,
	        label = F,
	        cols = color.use, 
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 2c Celltype UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	ggsave2(filename = "SupFig 2c Celltype UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	# Sup2d
	
	Idents(tissue.integrated) <- 'celltype'
	DefaultAssay(tissue.integrated) <- 'RNA'
	
	celltype.markers = FindAllMarkers(tissue.integrated,
	                                  only.pos = TRUE, 
	                                  test.use = "MAST",
	                                  latent.vars = "sex",
	                                  min.pct = 0.75,
	                                  logfc.threshold = 1.5,
	                                  assay = "RNA",
	                                  densify = TRUE) 
	
	VlnPlot(tissue.integrated,
	        pt.size = 0,
	        cols = c("#7DB3E2FF","#7DB3E2FF",
	                 "#59A14F","#59A14F", 
	                 "#EED58CFF","#EED58CFF", 
	                 "#E83800FF","#E83800FF", 
	                 "#D088E0FF","#D088E0FF",
	                 '#FFAE34', '#FFAE34',
	                 '#805840FF','#805840FF',
	                 '#808898FF','#808898FF'),
	        features = c('Krt5', 'Perp', 
	                     'Ptprc', 'Fcer1g', 
	                     'Col1a1', 'Dcn',
	                     'Pecam1', 'Gng11',
	                     'Mpz', 'Ncmap',
	                     'Muc19', 'Bpifb2',
	                     'Pmel', 'Mlana', 
	                     'Tnnc2', 'Tnni2'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 2d Celltype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	ggsave2(filename = "SupFig 2d Celltype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	#######
	
	#Supplementary Fig 3
	#######
	
	# Sup3b
	
	# load in Giotto objects
	
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	# Find celltype markers
	
	healthy.markers = findMarkers_one_vs_all(gobject = merge.healthy.giotto,
	                                         method = 'gini',
	                                         expression_values = 'normalized',
	                                         cluster_column = 'Celltype',
	                                         rank_score = 2)
	
	#Make Excel sheet with markers and cell numbers
	
	write.xlsx(healthy.markers,
	           # file = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Celltype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	healthy.markers <- readxl::read_xlsx('/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Celltype Markers.xlsx', 
	                                     sheet = 1,
	                                     col_names = T)
	
	healthy.markers[ ,1] <- NULL
	
	healthy.markers <- healthy.markers[order(healthy.markers$cluster,
	                                         -healthy.markers$expression), ]
	
	# plot top cell markers in heatmap
	
	marker.genes <- healthy.markers %>% slice_max(expression_gini, n=5, by=cluster)
	marker.genes <- unique(marker.genes[ ,1]) %>% unlist
	
	plotMetaDataHeatmap(merge.healthy.giotto, 
	                    custom_cluster_order = c('Fibroblast', 'Immune', 'Vascular', 'Neural', 'Mesenchymal', 'Keratinocyte', 'Salivary'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'Celltype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 3b Xenium Celltype Marker Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 3b Xenium Celltype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	#######
	
	#Fig 2
	#######
	
	tissue.kera <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds')
	kera.ife <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.ife.rds')
	kera.hf <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.hf.rds')
	
	# 2a
	Idents(kera.ife) <- 'L1subtype'
	
	DimPlot(kera.ife,
	        pt.size = 1,
	        split.by = 'state',
	        ncol = 3,
	        label = F,
	        cols = c('#0142FEFF','#802880FF','#324376FF', '#1CADE4FF', '#CCEEFFFF', '#C70E7BC1', '#F8D0F8FF'),
	        #Epithelial_Basal, Oral_Basal, Proliferating_Keratinocyte, Epithelial_Suprabasal_1, Epithelial_Suprabasal_2, Oral_Suprabasal_1, Oral_Suprabasal_2
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 2a IFE Keratinocyte L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(filename = "Fig 2a IFE Keratinocyte L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	kera.ife.test <- sc_utils(kera.ife)
	
	prop_test <- permutation_test(kera.ife.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "control_buccal_mucosa",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5))
	
	prop_test <- permutation_test(kera.ife.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "skin_dox_8_week",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5))
	
	table(kera.ife$state, kera.ife$L1subtype) #completed graph in Excel
	
	# 2b
	
	set.seed(42)
	
	ife.cds <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ife.cds.rds')
	
	plot_cells(cds = ife.cds,
	           color_cells_by = "pseudotime",
	           label_principal_points = T,
	           label_roots = F,
	           label_branch_points = F,
	           label_leaves = F,
	           cell_size = 1,
	           show_trajectory_graph = F) + NoLegend()
	
	ggsave2(filename = 'Fig 2b IFE Keratinocyte Pseudotime (No Trajectory Lines) UMAP.svg',
	        path = '/data/overmilleram/Manuscript Figures/',  
	        width = 5000,
	        height =5000,
	        units = "px")
	
	ggsave2(filename = "Fig 2b IFE Keratinocyte Pseudotime (No Trajectory Lines) UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	# 2c 
	
	set.seed(42)
	
	ife.cds <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ife.cds.rds')
	
	ife_cds_pr_test_res <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ifepseudotime.rds')
	
	pr_deg_ids <- row.names(subset(ife_cds_pr_test_res, q_value < 0.0001))
	
	ife.cds <- preprocess_cds(ife.cds) # need to do for find_gene_modules() to work
	
	ife_gene_module_df <- find_gene_modules(ife.cds[pr_deg_ids,], 
	                                        resolution = 0.001) 
	
	ife_cell_group_df <- tibble(cell=row.names(colData(ife.cds)), 
	                            cell_group=colData(ife.cds)$L1subtype)
	
	ife_agg_mat <- aggregate_gene_expression(ife.cds, 
	                                         ife_gene_module_df, 
	                                         ife_cell_group_df)
	
	row.names(ife_agg_mat) <- stringr::str_c("Module ", 
	                                         row.names(ife_agg_mat))
	
	ife.pseudo.genes <- ife_cds_pr_test_res[order(-ife_cds_pr_test_res$morans_I), ]
	ife.pseudo.genes <- subset(ife.pseudo.genes, q_value < 0.0001) 
	
	ife.pseudo.module <- ife_gene_module_df[order(ife_gene_module_df$module,
	                                              decreasing = F), ]
	
	pseudo_genes <- c('Krt77', 'Tgm3', 'Crabp2', 'Krt4', 'Krt6a', 'Krt6b', 'Sprr2a3', 'Dsc2', 'Cnfn', 'Sbsn', 'Krt5')
	
	skco_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('skin_control_8_week')]
	
	skdo_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('skin_dox_8_week')]
	
	orbu_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('control_buccal_mucosa')]
	
	skco <- plot_genes_in_pseudotime(skco_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	skdo <- plot_genes_in_pseudotime(skdo_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	orbu <- plot_genes_in_pseudotime(orbu_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	ggsave2(skco + NoLegend(),
	        filename = "Fig 2c IFE Keratinocyte Control Differentiation Genes over Pseudotime.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 4,
	        height = 12,
	        dpi = 1200,
	        units = "in")
	
	ggsave2(skdo + NoLegend(),
	        filename = "Fig 2c IFE Keratinocyte +Dox Differentiation Genes over Pseudotime.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 4,
	        height = 12,
	        dpi = 1200,
	        units = "in")
	
	ggsave2(orbu + NoLegend(),
	        filename = "Fig 2c IFE Keratinocyte Buccal Differentiation Genes over Pseudotime.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 4,
	        height = 12,
	        dpi = 1200,
	        units = "in")
	
	#Fig 2d 
	
	#Extract genes from positively correlated modules for each specific simplesubtype, run enrichR and generate Seurat module score
	
	ife_agg_mat_fil <- ife_agg_mat
	
	ife_agg_mat_fil <- ife_agg_mat_fil %>% as.data.frame() %>% rownames_to_column() 
	
	ife_agg_mat_fil[ ,1] <- gsub('Module ', '', ife_agg_mat_fil[ ,1]) %>% as.numeric(ife_agg_mat_fil[ ,1])
	
	epi.basal.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Basal > 0.1]
	ora.basal.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Basal > 0.1]
	prolifera.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Proliferating_Keratinocyte > 0.1]
	epi.supra1.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Suprabasal_1 > 0.1]
	epi.supra2.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Suprabasal_2 > 0.1]
	ora.supra1.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Suprabasal_1 > 0.1]
	ora.supra2.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Suprabasal_2 > 0.1]
	epi.supra.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Epithelial_Suprabasal_1 > 0.1 | ife_agg_mat_fil$Epithelial_Suprabasal_2 > 0.1]
	ora.supra.module <- ife_agg_mat_fil$rowname[ife_agg_mat_fil$Oral_Suprabasal_1 > 0.1 | ife_agg_mat_fil$Oral_Suprabasal_2 > 0.1]
	
	epi.basal.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.basal.module, ]
	ora.basal.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.basal.module, ]
	prolifera.genes <- ife.pseudo.module[ife.pseudo.module$module == prolifera.module, ]
	epi.supra1.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.supra1.module, ]
	epi.supra2.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.supra2.module, ]
	ora.supra1.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.supra1.module, ]
	ora.supra2.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.supra2.module, ]
	epi.supra.genes <- ife.pseudo.module[ife.pseudo.module$module == epi.supra.module, ]
	ora.supra.genes <- ife.pseudo.module[ife.pseudo.module$module == ora.supra.module, ]
	
	epi.basal.genes <- as.character(epi.basal.genes$id)
	ora.basal.genes <- as.character(ora.basal.genes$id)
	prolifera.genes <- as.character(prolifera.genes$id)
	epi.supra1.genes <- as.character(epi.supra1.genes$id)
	epi.supra2.genes <- as.character(epi.supra2.genes$id)
	ora.supra1.genes <- as.character(ora.supra1.genes$id)
	ora.supra2.genes <- as.character(ora.supra2.genes$id)
	epi.supra.genes <- as.character(epi.supra.genes$id)
	ora.supra.genes <- as.character(ora.supra.genes$id)
	
	# genes must be in a list for AddModuleScore to work
	epi.basal.genes <- list(epi.basal.genes)
	ora.basal.genes <- list(ora.basal.genes)
	prolifera.genes <- list(prolifera.genes)
	epi.supra1.genes <- list(epi.supra1.genes)
	epi.supra2.genes <- list(epi.supra2.genes)
	ora.supra1.genes <- list(ora.supra1.genes) 
	ora.supra2.genes <- list(ora.supra2.genes) 
	epi.supra.genes <- list(epi.supra.genes)
	ora.supra.genes <- list(ora.supra.genes) 
	
	kera.ife <- AddModuleScore(kera.ife,
	                           features = epi.basal.genes,
	                           ctrl = 100,
	                           name = 'Epi.Basal.Score',
	                           assay = 'RNA')
	
	kera.ife <- AddModuleScore(kera.ife,
	                          features = ora.basal.genes,
	                          ctrl = 100,
	                          name = 'Oral.Basal.Score',
	                          assay = 'RNA')
	
	kera.ife <- AddModuleScore(kera.ife,
	                           features = prolifera.genes,
	                           ctrl = 100,
	                           name = 'Proliferating.Kera.Score',
	                           assay = 'RNA')
	
	kera.ife <- AddModuleScore(kera.ife,
	                             features = epi.supra1.genes,
	                             ctrl = 100,
	                             name = 'Epi.Suprabasal_1.Score',
	                             assay = 'RNA')
	
	kera.ife <- AddModuleScore(kera.ife,
	                             features = epi.supra2.genes,
	                             ctrl = 100,
	                             name = 'Epi.Suprabasal_2.Score',
	                             assay = 'RNA')
	
	kera.ife <- AddModuleScore(kera.ife,
	                             features = ora.supra1.genes,
	                             ctrl = 100,
	                             name = 'Oral.Suprabasal_1.Score',
	                             assay = 'RNA')
	
	kera.ife <- AddModuleScore(kera.ife,
	                             features = ora.supra2.genes,
	                             ctrl = 100,
	                             name = 'Oral.Suprabasal_2.Score',
	                             assay = 'RNA')
	
	Idents(kera.ife) <- 'L1subtype'
	DefaultAssay(kera.ife) <- 'RNA'
	
	kera.ife.supra <- subset(kera.ife,
	                           idents = c('Epithelial_Suprabasal_1', 'Epithelial_Suprabasal_2', 'Oral_Suprabasal_1', 'Oral_Suprabasal_2')) 
	
	Idents(kera.ife.supra) <- 'celltype'
	
	VlnPlot(kera.ife.supra,
	        features = c('Oral.Suprabasal_1.Score1', 'Oral.Suprabasal_2.Score1', 'Epi.Suprabasal_1.Score1', 'Epi.Suprabasal_2.Score1'),
	        pt.size = 0,
	        flip = T,
	        stack = T,
	        cols = c('#F9B90AFF','#3D79F3FF','#34A74BFF'),
	        split.by = 'state',
	        assay = 'RNA',
	        raster = F) + NoLegend()
	ggsave2(filename = "Fig 2d Keratinocyte Suprabasal IFE Pseudotime Gene Module (Suprabasal Modules) VlnPlots.svg",  
	        path = '/data/overmilleram/Manuscript Figures/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#extract module scores to perform Kruskal-Wallis
	
	module.score <- data.frame(EpiSuprabasal1 = kera.ife.supra$Epi.Suprabasal_1.Score1,
	                           EpiSuprabasal2 = kera.ife.supra$Epi.Suprabasal_2.Score1,
	                           OralSuprabasal1 = kera.ife.supra$Oral.Suprabasal_1.Score1,
	                           OralSuprabasal2 = kera.ife.supra$Oral.Suprabasal_2.Score1,
	                           state = kera.ife.supra$state)
	
	ggpubr::ggboxplot(module.score, #check distrubution on a box plot
	          x = "state", 
	          y = "OralSuprabasal2", 
	          color = "state", 
	          palette = c('#3D79F3FF','#34A74BFF','#F9B90AFF'),
	          order = c('skin_control_8_week', 'skin_dox_8_week', 'control_buccal_mucosa'),
	          ylab = "Module Score", 
	          xlab = F)
	
	pairwise.wilcox.test(module.score$EpiSuprabasal1, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	pairwise.wilcox.test(module.score$EpiSuprabasal2, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	
	pairwise.wilcox.test(module.score$OralSuprabasal1, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	pairwise.wilcox.test(module.score$OralSuprabasal2, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	
	#######
	
	#Supplementary Fig 4
	#######
	
	tissue.kera <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.kera.rds')
	
	cols = c('#0142FEFF','#802880FF','#324376FF', '#1CADE4FF', '#CCEEFFFF', '#C70E7BC1', '#F8D0F8FF',
	         '#B0D8D0FF', '#6088A0FF', '#486048FF', '#5860B8FF', '#219089FF', '#A45D6BFF')
	
	#Epithelial_Basal, Oral_Basal, Proliferating_Keratinocyte, Epithelial_Suprabasal_1, Epithelial_Suprabasal_2, Oral_Suprabasal_1, Oral_Suprabasal_2,
	#HFSC, Lower_HF, Sebaceous, Upper_HF_Basal, Upper_HF_Suprabasal, Stress_Basal
	
	
	#Sup4a
	Idents(tissue.kera) <- 'L1subtype'
	
	DimPlot(tissue.kera,
	        cols = cols,
	        pt.size = 1,
	        label = F,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 4a Keratinocyte L1 Subtype UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	ggsave2(filename = "SupFig 4a Keratinocyte L1 Subtype UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	#Sup4b
	
	Idents(tissue.kera) <- 'L1subtype'
	
	keraL1subtype.markers = FindAllMarkers(tissue.kera,
	                                       only.pos = TRUE, 
	                                       test.use = "poisson",
	                                       latent.vars = "sex",
	                                       min.pct = 0.50,
	                                       logfc.threshold = 2,
	                                       assay = "RNA",
	                                       densify = TRUE) 
	
	keraL1subtype.markers <- keraL1subtype.markers[order(keraL1subtype.markers$cluster,
	                                                     -keraL1subtype.markers$avg_log2FC), ]
	
	write.xlsx(keraL1subtype.markers,
	           file = '/data/overmilleram/Manuscript Figures/SupFig 4b Keratinocyte L1 Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	Lower_HF_markers <- FindMarkers(tissue.kera,
	                                ident.1 = 'Lower_HF',
	                                only.pos = T,
	                                assay = 'RNA',
	                                test.use = 'MAST',
	                                latent.vars = 'sex',
	                                min.pct = 0.25,
	                                logfc.threshold = 1,
	                                densify = T)
	
	#Epithelial_Basal, Oral_Basal, Proliferating_Keratinocyte, Epithelial_Suprabasal_1, Epithelial_Suprabasal_2, Oral_Suprabasal_1, Oral_Suprabasal_2,
	#HFSC, Lower_HF, Sebaceous, Upper_HF_Basal, Upper_HF_Suprabasal, Stress_Basal
	
	cols = c('#0142FEFF','#802880FF','#324376FF', '#1CADE4FF', '#CCEEFFFF', '#C70E7BC1', '#F8D0F8FF',
	         '#B0D8D0FF', '#6088A0FF', '#486048FF', '#5860B8FF', '#219089FF', '#A45D6BFF')
	
	VlnPlot(tissue.kera,
	        pt.size = 0,
	        cols = c(rep('#0142FEFF', times = 2),
	                 rep('#802880FF', times = 2),
	                 rep('#324376FF', times = 2),
	                 rep('#1CADE4FF', times = 2),
	                 rep('#CCEEFFFF', times = 2),
	                 rep('#C70E7BC1', times = 2),
	                 rep('#F8D0F8FF', times = 2),
	                 rep('#B0D8D0FF', times = 2),
	                 rep('#6088A0FF', times = 2),
	                 rep('#486048FF', times = 2),
	                 rep('#5860B8FF', times = 2),
	                 rep('#219089FF', times = 2),
	                 rep('#A45D6BFF', times = 2)),
	       features = c('Krt14', 'Dsc2', 
	                    'Sox6', 'Igfbp2', 
	                    'Hist1h2ap', 'Mki67', 
	                    'Krt10', 'Clca3a2', 
	                    'Hal', 'Ivl', 
	                    'Selenbp1', 'Gsta4',
	                    'Serpinb3a', 'Ada',
	                    'S100a4', 'Postn', 
	                    'Sox4', 'Sox5',
	                    'Cidea', 'Plin2', 
	                    'Thbs1', 'Itga2',
	                    'Fst', 'Krt17', 
	                    'Fos', 'Atf3'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 4b Keratinocyte L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 4b Keratinocyte L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#Sup4c
	
	Idents(tissue.kera) <- 'L1subtype'
	
	DimPlot(tissue.kera,
	        pt.size = 0.75,
	        split.by = "state",
	        label = F,
	        cols = cols,
	        ncol = 3,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 4c Keratinocyte L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	
	ggsave2(filename = "SupFig 4c Keratinocyte L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	table(tissue.kera$state, tissue.kera$L1subtype) #completed graph in Excel
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	skin.kera.test = sc_utils(tissue.kera)
	
	prop_test = permutation_test(skin.kera.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "skin_dox_8_week",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(2)) 
	
	table(tissue.kera$state)
	
	table(tissue.kera$state, tissue.kera$L1subtype) #completed in Excel
	
	#Sup4d
	Idents(kera.hf) <- 'state'
	
	kera.hf <- subset(kera.hf, idents = c('skin_control_8_week', 'skin_dox_8_week'))
	
	kera.hf$L1subtype <- droplevels(kera.hf$L1subtype)
	
	Idents(kera.hf) <- 'L1subtype'
	
	cols = c('#324376FF', '#B0D8D0FF', '#6088A0FF', '#486048FF', '#5860B8FF', '#219089FF', '#A45D6BFF')
	#Proliferating_Keratinocyte, HFSC, Lower_HF, Sebaceous, Upper_HF_Basal, Upper_HF_Suprabasal, Stress_Basal
	
	DimPlot(kera.hf,
	        pt.size = 0.75,
	        label = F,
	        split.by = 'state',
	        ncol = 2,
	        cols = cols,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 4d Keratinocyte HF L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 4d Keratinocyte HF L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	kera.hf.test <- sc_utils(kera.hf)
	
	prop_test <- permutation_test(kera.hf.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "skin_dox_8_week",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(2))
	
	table(kera.hf$state, kera.hf$L1subtype) #completed graph in Excel
	
	#Sup4f
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	metadata = pDataDT(merge.healthy.giotto)
	
	subset_cell_IDs = metadata[L1Subtype %in% c('Dermal Papilla', 'Epithelial Basal 1', 'HFSC', 'Lower HF 1', 
	                                            'Lower HF 2', 'Lower HF Inflamm', 'Sebaceous', 'Upper HF Suprabasal 1')]$cell_ID
	
	hf.giotto = subsetGiotto(merge.healthy.giotto, cell_ids = subset_cell_IDs)
	
	hf.meta <- pDataDT(hf.giotto) # pull existing metadata
	
	hf.meta$L1Subtype <- gsub(c('Lower HF 1|Lower HF 2|Lower HF Inflamm'), 'Lower HF', hf.meta$L1Subtype) # combine Lower HF populations into one cluster
	
	hf.giotto <- addCellMetadata(hf.giotto,
	                             spat_unit = 'cell',
	                             feat_type = 'rna',
	                             hf.meta,
	                             by_column = F)
	
	unique(pDataDT(hf.giotto)$L1Subtype) # check if metadata correct
	
	marker.genes <- c('Krt15', 'Lgr5', 'Sox9', 'Dlx3', 'Dsg4', 'Ptch1', 'Fabp5', 'Elovl6', 'Cidea', 
	                  'Krt5', 'Col17a1', 'Krt77', 'Krt17', 'Sprr1a', 'Plet1', 'Wif1', 'Vim', 'Cyp26b1')
	
	plotMetaDataHeatmap(hf.giotto, 
	                    custom_cluster_order = c('HFSC', 'Lower HF', 'Sebaceous', 'Epithelial Basal 1', 'Upper HF Suprabasal 1', 'Dermal Papilla'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes,
	                    custom_feat_order = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 4f Xenium HF L1Subtype Marker Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 4f Xenium HF L1Subtype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	#######
	
	#Supplementary Fig 5
	#######
	
	#Sup5a
	
	set.seed(42)
	
	ife.cds <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ife.cds.rds')
	
	ife_cds_pr_test_res <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ifepseudotime.rds')
	
	pr_deg_ids <- row.names(subset(ife_cds_pr_test_res, q_value < 0.0001))
	
	ife.cds <- preprocess_cds(ife.cds) # need to do for find_gene_modules() to work
	
	ife_gene_module_df <- find_gene_modules(ife.cds[pr_deg_ids,], 
	                                        resolution = 0.001) 
	
	ife_cell_group_df <- tibble(cell=row.names(colData(ife.cds)), 
	                            cell_group=colData(ife.cds)$L1subtype)
	
	ife_agg_mat <- aggregate_gene_expression(ife.cds, 
	                                         ife_gene_module_df, 
	                                         ife_cell_group_df)
	
	row.names(ife_agg_mat) <- stringr::str_c("Module ", 
	                                         row.names(ife_agg_mat))
	
	ife.pseudo.genes <- ife_cds_pr_test_res[order(-ife_cds_pr_test_res$morans_I), ]
	ife.pseudo.genes <- subset(ife.pseudo.genes, q_value < 0.0001) 
	
	ife.pseudo.module <- ife_gene_module_df[order(ife_gene_module_df$module,
	                                              decreasing = F), ]
	
	pseudo_genes <- c('Pitx1', 'Pitx2', 'Pax9', 'Sox2', 'Trp63')
	
	skco_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('skin_control_8_week')]
	
	skdo_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('skin_dox_8_week')]
	
	orbu_pseudo_genes_cds <- ife.cds[rowData(ife.cds)$gene_short_name %in% pseudo_genes,
	                                 colData(ife.cds)$state %in% c('control_buccal_mucosa')]
	
	skco <- plot_genes_in_pseudotime(skco_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	skdo <- plot_genes_in_pseudotime(skdo_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	orbu <- plot_genes_in_pseudotime(orbu_pseudo_genes_cds,
	                                 color_cells_by = 'pseudotime',
	                                 min_expr = 0.5,
	                                 cell_size = 0.4,
	                                 panel_order = pseudo_genes,
	                                 label_by_short_name = T)
	
	ggsave2(skco + NoLegend(),
	        filename = "SupFig 5a IFE Keratinocyte Control Differentiation Transcription Factor Genes over Pseudotime.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 2,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	ggsave2(skdo + NoLegend(),
	        filename = "SupFig 5a IFE Keratinocyte Pitx1+ Differentiation Transcription Factor Genes over Pseudotime.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 2,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	ggsave2(orbu + NoLegend(),
	        filename = "SupFig 5a IFE Keratinocyte Buccal Differentiation Transcription Factor Genes over Pseudotime.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 2,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	#Sup5b
	
	set.seed(42)
	
	ife_cds_pr_test_res <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/ifepseudotime.rds')
	
	pr_deg_ids <- row.names(subset(ife_cds_pr_test_res, q_value < 0.0001))
	
	ife.cds <- preprocess_cds(ife.cds) # need to do for find_gene_modules() to work
	
	ife_gene_module_df <- find_gene_modules(ife.cds[pr_deg_ids,], 
	                                        resolution = 0.001) 
	
	ife_cell_group_df <- tibble(cell=row.names(colData(ife.cds)), 
	                            cell_group=colData(ife.cds)$L1subtype)
	
	ife_agg_mat <- aggregate_gene_expression(ife.cds, 
	                                         ife_gene_module_df, 
	                                         ife_cell_group_df)
	
	row.names(ife_agg_mat) <- stringr::str_c("Module ", 
	                                         row.names(ife_agg_mat))
	
	ife.pseudo.genes <- ife_cds_pr_test_res[order(-ife_cds_pr_test_res$morans_I), ]
	ife.pseudo.genes <- subset(ife.pseudo.genes, q_value < 0.0001) 
	
	ife.pseudo.module <- ife_gene_module_df[order(ife_gene_module_df$module,
	                                              decreasing = F), ]
	
	pheatmap <- pheatmap::pheatmap(ife_agg_mat,
	                               color = paletteer_c("scico::bam", 100),
	                               scale = "column", 
	                               cluster_rows = T,
	                               cluster_cols = T,
	                               fontsize = 8,
	                               clustering_method = "ward.D2")
	
	ggsave2(pheatmap,
	        filename = 'SupFig 5b IFE Keratinocyte Pseudotime Modules Heatmap.svg', 
	        path = '/data/overmilleram/Manuscript Figures/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#Sup5c
	
	Idents(kera.ife) <- 'L1subtype'
	DefaultAssay(kera.ife) <- 'RNA'
	
	kera.ife.basal <- subset(kera.ife,
	                         idents = c('Epithelial_Basal', 'Oral_Basal', 'Proliferating_Keratinocyte')) 
	
	Idents(kera.ife.basal) <- 'celltype'
	
	VlnPlot(kera.ife.basal,
	        features = c('Epi.Basal.Score1', 'Oral.Basal.Score1', 'Proliferating.Kera.Score1'),
	        pt.size = 0,
	        flip = T,
	        stack = T,
	        cols = c('#F9B90AFF','#3D79F3FF','#34A74BFF'),
	        split.by = 'state',
	        assay = 'RNA',
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 5c Keratinocyte Basal IFE Pseudotime Gene Module (Basal Modules) VlnPlots.svg",  
	        path = '/data/overmilleram/Manuscript Figures/',
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	#extract module scores to perform Kruskal-Wallis
	
	module.score <- data.frame(EpiBasal = kera.ife.basal$Epi.Basal.Score1, 
	                           OralBasal = kera.ife.basal$Oral.Basal.Score1,
	                           Proliferative = kera.ife.basal$Proliferating.Kera.Score1,
	                           state = kera.ife.basal$state)
	
	ggpubr::ggboxplot(module.score, #check distrubution on a box plot
	                  x = "state", 
	                  y = "Proliferative", 
	                  color = "state", 
	                  palette = c('#3D79F3FF','#34A74BFF','#F9B90AFF'),
	                  order = c('skin_control_8_week', 'skin_dox_8_week', 'control_buccal_mucosa'),
	                  ylab = "Module Score", 
	                  xlab = F)
	
	pairwise.wilcox.test(module.score$EpiBasal, module.score$state, p.adjust.method = 'BH') #control vs buccal p < 2e-16 (control vs dox p = 0.054)
	pairwise.wilcox.test(module.score$OralBasal, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	pairwise.wilcox.test(module.score$Proliferative, module.score$state, p.adjust.method = 'BH') #all comparisons are significant (p < 2e-16 at least)
	
	#######
	
	#Fig 3
	#######
	
	tissue.imm <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds')
	tissue.fib <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds')
	tissue.vasc <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds')
	tissue.neur <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds')
	
	#3a
	Idents(tissue.fib) <- 'L1subtype'
	
	DimPlot(tissue.fib,
	        pt.size = 1.5,
	        label = F,
	        split.by = "state",
	        ncol = 3,
	        cols = c('#FFFF80FF','#FFD700FF','#F19425FF', '#BE5C00FF', '#C4A000FF', '#73652DFF'),
	        #Dermal_Fibroblast_1, Dermal_Fibroblast_2, Stromal_Fibroblast_1, Stromal_Fibroblast_2, Dermal_Papilla, Dermal_Sheath
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 3a Fibroblast L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "Fig 3a Fibroblast L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	tissue.fib.test <- sc_utils(tissue.fib)
	
	prop_test <- permutation_test(tissue.fib.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "control_buccal_mucosa",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4))
	
	prop_test <- permutation_test(tissue.fib.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "skin_dox_8_week",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4))
	
	table(tissue.fib$state, tissue.fib$L1subtype) #completed graph in Excel
	
	#3b
	Idents(tissue.imm) <- 'L1subtype'
	
	DimPlot(tissue.imm,
	        pt.size = 2,
	        label = F,
	        split.by = "state",
	        ncol = 3,
	        cols = c('#FCD47CFF','#AC570FFF','#F57206FF', '#97AD3DFF','#15CC31FF', 
	                 '#CCF4CFFF','#75FB8AFF','#7FD4C1FF', '#CD3122FF','#E7A6CDFF', '#008D98FF','#415521FF'),
	        #'T_cell', 'Proliferating_T_Cell', 'NK_Cell', 'Macrophage_1', 'Macrophage_2', 
	        #'Dendritic_1', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2', 'Dendritic_2', 'Macrophage_3'
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 3b Immune L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "Fig 3b Immune L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	tissue.imm.test <- sc_utils(tissue.imm)
	
	prop_test <- permutation_test(tissue.imm.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "control_buccal_mucosa",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4))
	
	prop_test <- permutation_test(tissue.imm.test, 
	                              cluster_identity = "L1subtype",
	                              sample_1 = "skin_control_8_week", 
	                              sample_2 = "skin_dox_8_week",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4))
	
	table(tissue.imm$state, tissue.imm$L1subtype) #completed graph in Excel
	
	#3d
	
	tissue.neu <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds')
	
	Idents(tissue.neu) <- 'L1subtype'
	
	tissue.neu <- subset(tissue.neu, idents = c('Neutrophil_1', 'Neutrophil_2'))
	
	tissue.neu$L1subtype <- droplevels(tissue.neu$L1subtype)
	tissue.neu$L2subtype <- droplevels(tissue.neu$L2subtype)
	
	tissue.neu <- RunUMAP(tissue.neu, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	DimPlot(tissue.neu,
	        pt.size = 4,
	        label = F,
	        split.by = "state",
	        ncol = 3,
	        cols = c('#CD3122FF','#E7A6CDFF'),
	        #Neutrophil_1, Neutrophil_2
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 3d Neutrophil Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "Fig 3d Neutrophil Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	#3e
	
	multinichenet_neutrophil_output <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_neutrophil_output.rds')
	
	Idents(tissue.integrated) <- 'celltype' # remove Salivary, Melanocyte, & Skeletal Muscle cells from analysis
	
	tissue.integrated2 <- subset(tissue.integrated,
	                            idents = c('Salivary', 'Skeletal_Muscle', 'Melanocyte'),
	                            invert = T)
	
	tissue.integrated2$celltype <- droplevels(tissue.integrated2$celltype) # remove unused factor levels
	tissue.integrated2$L1subtype <- droplevels(tissue.integrated2$L1subtype) # remove unused factor levels
	tissue.integrated2$L2subtype <- droplevels(tissue.integrated2$L2subtype) # remove unused factor levels
	
	tissue.sce <- as.SingleCellExperiment(tissue.integrated2, assay = "RNA")
	
	tissue.sce <- alias_to_symbol_SCE(tissue.sce, "mouse")
	
	rm(tissue.integrated2)
	
	senders_receivers <- unique(levels(tissue.sce$L1subtype))
	
	color.use = c('#0142FEFF','#802880FF','#324376FF', '#1CADE4FF', '#CCEEFFFF', '#C70E7BC1', '#F8D0F8FF',
	              '#B0D8D0FF', '#6088A0FF', '#486048FF', '#5860B8FF', '#219089FF', '#A45D6BFF',
	              '#FCD47CFF','#AC570FFF','#F57206FF', '#97AD3DFF','#15CC31FF', 
	              '#CCF4CFFF','#75FB8AFF','#7FD4C1FF', '#CD3122FF','#E7A6CDFF', '#008D98FF','#415521FF',
	              '#FFFF80FF','#FFD700FF','#F19425FF', '#BE5C00FF', '#C4A000FF', '#73652DFF',
	              '#FF0010FF', '#FFB79FFF', '#D7CADEFF', '#9632B8FF','#990000FF')
	
	names(color.use) <- senders_receivers
	
	# Epithelial_Basal, Oral_Basal, Proliferating_Keratinocyte, Epithelial_Suprabasal_1, Epithelial_Suprabasal_2, Oral_Suprabasal_1, Oral_Suprabasal_2,
	# HFSC, Lower_HF, Sebaceous, Upper_HF_Basal, Upper_HF_Suprabasal, Stress_Basal
	# 'T_cell', 'Proliferating_T_Cell', 'NK_Cell', 'Macrophage_1', 'Macrophage_2', 
	# 'Dendritic_1', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2', 'Dendritic_2', 'Macrophage_3'
	# Dermal_Fibroblast_1, Dermal_Fibroblast_2, Stromal_Fibroblast_1, Stromal_Fibroblast_2, Dermal_Papilla, Dermal_Sheath
	# Endothelial, Lymph_Vessel, Myelinating Schwann, Non-myelinating Schwann, Vascular Smooth Muscle
	
	
	##Define contrasts 
	#
	#contrasts_oi <- c("'skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2','skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2','control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2'")
	#
	#contrast_tbl <- tibble(contrast = c("skin_control_8_week-(skin_dox_8_week+control_buccal_mucosa)/2","skin_dox_8_week-(skin_control_8_week+control_buccal_mucosa)/2","control_buccal_mucosa-(skin_control_8_week+skin_dox_8_week)/2"),
	#                       group = c("skin_control_8_week","skin_dox_8_week","control_buccal_mucosa")) 
	#
	##visualization of top 50 interactions across all groups
	#prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_neutrophil_output$prioritization_tables, 50, rank_per_group = FALSE)
	#
	#prioritized_tbl_oi = multinichenet_neutrophil_output$prioritization_tables$group_prioritization_tbl %>%
	#  filter(id %in% prioritized_tbl_oi_all$id) %>%
	#  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
	#
	#prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0
	#
	#circos_list = make_circos_group_comparison(prioritized_tbl_oi, color.use, color.use)
	#
	#circos_list[[2]]
	#dev.print(svg,
	#          "/data/overmilleram/Manuscript Figures/Fig 3e Healthy Tissue Control Top50 Shared LR Links Neutrophil Circos.svg",
	#          width = 5.5,
	#          height = 5.5)
	#dev.off()
	#
	#circos_list[[1]]
	#dev.print(svg,
	#          "/data/overmilleram/Manuscript Figures/Fig 3e Healthy Tissue +Dox Top50 Shared LR Links Neutrophil Circos.svg",
	#          width = 5.5,
	#          height = 5.5)
	#dev.off()
	#
	#circos_list[[4]]
	#dev.print(svg,
	#          "/data/overmilleram/Manuscript Figures/Fig 3e Healthy Tissue Legend Top50 Shared LR Links Neutrophil.svg",
	#          width = 5,
	#          height = 10)
	#dev.off()
	
	#top 50 skin_dox_8_week
	
	prioritized_tbl_oi_skdo_40 = get_top_n_lr_pairs(multinichenet_neutrophil_output$prioritization_tables, 
	                                                40, 
	                                                groups_oi = "skin_dox_8_week")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_40, color.use, color.use)
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/Fig 3e Healthy Tissue Pitx1+ Top40 Neutrophil-L1subtype Circos.svg",
	          width = 5.6,
	          height = 5.6)
	dev.off()
	
	circos_skdo[[2]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/Fig 3e Healthy Tissue Pitx1+ Top40 Neutrophil-L1subtype Circos Legend.svg",
	          width = 5,
	          height = 10)
	dev.off()
	
	#######
	
	#Supplementary Fig 7
	#######
	
	tissue.fib <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.fib.rds")
	tissue.imm <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.imm.rds")
	tissue.vasc <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds")
	tissue.neur <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds")
	tissue.mes <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds")
	tissue.sal <- readRDS("/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.sal.rds")
	
	#Sup6a
	
	cols = c('#FFFF80FF','#FFD700FF','#F19425FF', '#BE5C00FF', '#C4A000FF', '#73652DFF')
	#Dermal_Fibroblast_1, Dermal_Fibroblast_2, Stromal_Fibroblast_1, Stromal_Fibroblast_2, Dermal_Papilla, Dermal_Sheath
	
	Idents(fib.integrated) <- 'L1subtype'
	
	fib.L1subtype.markers = FindAllMarkers(fib.integrated,
	                                       only.pos = TRUE, 
	                                       test.use = "MAST",
	                                       latent.vars = "sex",
	                                       min.pct = 0.5,
	                                       logfc.threshold = 2,
	                                       assay = "RNA",
	                                       densify = TRUE) 
	
	fib.L1subtype.markers <- fib.L1subtype.markers[order(fib.L1subtype.markers$cluster,
	                                                     -fib.L1subtype.markers$avg_log2FC), ]
	
	write.xlsx(fib.L1subtype.markers,
	           file = '/data/overmilleram/Manuscript Figures/SupFig 6a Fibroblast L1 Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	Derm2_markers <- FindMarkers(fib.integrated,
	                             ident.1 = 'Dermal_Fibroblast_2',
	                             only.pos = T,
	                             assay = 'RNA',
	                             test.use = 'MAST',
	                             latent.vars = 'sex',
	                             min.pct = 0.25,
	                             logfc.threshold = 1,
	                             densify = T)
	
	VlnPlot(fib.integrated,
	        pt.size = 0,
	        cols = c(rep('#FFFF80FF', times = 3),
	                 rep('#FFD700FF', times = 3),
	                 rep('#F19425FF', times = 3),
	                 rep('#BE5C00FF', times = 3),
	                 rep('#C4A000FF', times = 3),
	                 rep('#73652DFF', times = 3)),
	        features = c('Col1a1', 'Lum', 'Tgfbi',
	                     'Angptl1', 'Spon2', 'S100a4',
	                     'Cxcl1', 'Cxcl2', 'Plac8',
	                     'Col6a5', 'Dpp4', 'Smpd3',
	                     'Bcl2', 'Dkk2', 'Spon1', 
	                     'Mylk', 'Tshz3', 'Slit2'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 6a Fibroblast L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 6a Fibroblast L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#Sup6b
	
	Idents(tissue.imm) <- 'L1subtype'
	
	imm.L1subtype.markers = FindAllMarkers(tissue.imm,
	                                       only.pos = TRUE, 
	                                       test.use = "MAST",
	                                       latent.vars = "sex",
	                                       min.pct = 0.5,
	                                       logfc.threshold = 2,
	                                       assay = "RNA",
	                                       densify = TRUE) 
	
	imm.L1subtype.markers <- imm.L1subtype.markers[order(imm.L1subtype.markers$cluster,
	                                                     -imm.L1subtype.markers$avg_log2FC), ]
	
	neu1.v.neu2.markers <- FindMarkers(tissue.imm,
	                                   ident.1 = 'Neutrophil_1',
	                                   ident.2 = 'Neutrophil_2', 
	                                   assay = 'RNA',
	                                   test.use = 'MAST',
	                                   latent.vars = 'sex',
	                                   logfc.threshold = 0.5,
	                                   min.pct = 0.25,
	                                   densify = T)
	
	
	write.xlsx(imm.L1subtype.markers,
	           file = '/data/overmilleram/Manuscript Figures/SupFig 6b Immune L1 Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	tissue.imm[['L1subtype']] <- factor(x = tissue.imm@meta.data$L1subtype, 
	                                    levels = c('T_cell', 'Proliferating_T_Cell', 'NK_Cell', 'Macrophage_1', 'Macrophage_2', 'Macrophage_3',
	                                               'Dendritic_1', 'Dendritic_2', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2'))
	
	
	cols = c('#FCD47CFF','#AC570FFF','#F57206FF','#97AD3DFF','#15CC31FF','#415521FF',
	         '#CCF4CFFF','#008D98FF','#75FB8AFF','#7FD4C1FF','#CD3122FF','#E7A6CDFF')
	#'T_cell', 'Proliferating_T_Cell', 'NK_Cell', 'Macrophage_1', 'Macrophage_2', Macrophage_3,
	#'Dendritic_1', 'Dendritic_2',  'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2'
	
	names(cols) <- factor_level
	
	VlnPlot(tissue.imm,
	        pt.size = 0,
	        cols = c(rep('#FCD47CFF', times = 2),
	                 rep('#AC570FFF', times = 2),
	                 rep('#F57206FF', times = 2),
	                 rep('#97AD3DFF', times = 2),
	                 rep('#15CC31FF', times = 2),
	                 rep('#415521FF', times = 2),
	                 rep('#CCF4CFFF', times = 2),
	                 rep('#008D98FF', times = 2),
	                 rep('#75FB8AFF', times = 2),
	                 rep('#7FD4C1FF', times = 2),
	                 rep('#CD3122FF', times = 2),
	                 rep('#E7A6CDFF', times = 2)),
	        features = c('Icos', 'Camk4',
	                     'Hist1h1b', 'Top2a',
	                     'Trdc', 'Xcl1',
	                     'Ccl7', 'Folr2',
	                     'Spp1', 'Gpnmb',
	                     'Lrmda', 'Msrb1',
	                     'H2-Eb1', 'Cd74',
	                     'Clec4e', 'Cxcl2',
	                     'Cd207', 'Mfge8',
	                     'Plac8', 'Fn1',
	                     'Ifitm6', 'Retnlg',
	                     'S100a8', 'Cstb'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 6b Immune L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 6b Immune L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#Sup6c
	
	tissue.vasc <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.vasc.rds')
	tissue.neur <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neur.rds')
	tissue.mes <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.mes.rds')
	
	Idents(tissue.vasc) <- 'L1subtype'
	
	DimPlot(tissue.vasc,
	        pt.size = 3,
	        label = F,
	        split.by = "state",
	        ncol = 3,
	        cols = c('#FF0010FF', '#FFB79FFF'),
	        #Endothelial, Lymph_Vessel, Vascular Smooth Muscle
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 6c Vascular L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "SupFig 6c Vascular L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	vasc.test = sc_utils(tissue.vasc)
	
	prop_test = permutation_test(vasc.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "skin_dox_8_week",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5)) 
	
	prop_test = permutation_test(vasc.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "control_buccal_mucosa",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5)) 
	
	table(tissue.vasc$state, tissue.vasc$L1subtype) #completed graph in Excel 
	
	cols = c('#FF0010FF', '#FFB79FFF')
	#Endothelial, Lymph_Vessel
	
	Idents(tissue.vasc) <- 'L1subtype'
	
	vasc.L1subtype.markers = FindAllMarkers(tissue.vasc,
	                                       only.pos = TRUE, 
	                                       test.use = "poisson",
	                                       latent.vars = "sex",
	                                       min.pct = 0.5,
	                                       logfc.threshold = 2,
	                                       assay = "RNA",
	                                       densify = TRUE) 
	
	VlnPlot(tissue.vasc,
	        pt.size = 0,
	        cols = c('#FF0010FF','#FF0010FF','#FF0010FF',
	                 '#FFB79FFF','#FFB79FFF','#FFB79FFF'),
	        features = c('Flt1', 'Cd36', 'Pecam1',
	                     'Ccl21a', 'Lyve1', 'Mmrn1'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 6c Vascular L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 6c Vascular L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#Sup6d
	
	Idents(tissue.neur) <- 'L1subtype'
	
	DimPlot(tissue.neur,
	        pt.size = 4,
	        label = F,
	        split.by = "state",
	        ncol = 3,
	        cols = c('#D7CADEFF', '#9632B8FF', '#805840FF'),
	        #Myelinating Schwann, Non-myelinating Schwann
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 6d Neural L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "SupFig 6d Neural L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	neur.test = sc_utils(tissue.neur)
	
	prop_test = permutation_test(neur.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "skin_dox_8_week",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(2)) 
	
	prop_test = permutation_test(neur.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "control_buccal_mucosa",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5)) 
	
	table(tissue.neur$state, tissue.neur$L1subtype) #completed graph in Excel 
	
	cols = c('#D7CADEFF', '#9632B8FF', '#805840FF')
	#Myelinating Schwann, Non-myelinating Schwann
	
	Idents(tissue.neur) <- 'L1subtype'
	
	neur.L1subtype.markers = FindAllMarkers(tissue.neur,
	                                        only.pos = TRUE, 
	                                        test.use = "poisson",
	                                        latent.vars = "sex",
	                                        min.pct = 0.5,
	                                        logfc.threshold = 2,
	                                        assay = "RNA",
	                                        densify = TRUE) 
	
	DefaultAssay(tissue.neur) <- 'RNA'
	
	ana.hf <- FindMarkers(tissue.neur,
	                      ident.1 = 'Anagen_HF',
	                      only.pos = T,
	                      test.use = "poisson",
	                      latent.vars = "sex",
	                      min.pct = 0.33,
	                      logfc.threshold = 0.50,
	                      assay = "RNA",
	                      densify = TRUE)
	
	VlnPlot(tissue.neur,
	        pt.size = 0,
	        cols = c('#D7CADEFF','#D7CADEFF','#D7CADEFF',
	                 '#9632B8FF','#9632B8FF','#9632B8FF',
	                 '#805840FF','#805840FF','#805840FF'),
	        features = c('Mpz', 'Ncmap', 'Ctnna3',
	                     'Csmd1', 'Kcna2', 'Scn7a',
	                     'Mlana','Pmel','Kit'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 6d Neural L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 6d Neural L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#Sup6e
	
	Idents(tissue.mes) <- 'L1subtype'
	
	DimPlot(tissue.mes,
	        pt.size = 6,
	        label = F,
	        split.by = "state",
	        ncol = 3,
	        cols = c('#990000FF', '#808898FF'),
	        #Vascular_Smooth_Muscle, Skeletal_Muscle
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 6e Mesenchymal L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "SupFig 6e Mesenchymal L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	mes.test = sc_utils(tissue.mes)
	
	prop_test = permutation_test(mes.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "skin_dox_8_week",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5)) # no significant changes
	
	prop_test = permutation_test(mes.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "skin_control_8_week", 
	                             sample_2 = "control_buccal_mucosa",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(1.5)) # no significant changes
	
	table(tissue.mes$state, tissue.mes$L1subtype) #completed graph in Excel 
	
	#Myelinating Schwann, Non-myelinating Schwann
	
	VlnPlot(tissue.mes,
	        pt.size = 0,
	        cols = c('#990000FF','#990000FF','#990000FF',
	                 '#808898FF','#808898FF','#808898FF'),
	        features = c('Myh11', 'Tagln', 'Acta2',
	                     'Mylpf', 'Tnnc2', 'Tnni2'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 6e Mesenchymal L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 6e Mesenchymal L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#######
	
	#Supplementary Fig 8
	#######
	
	#Sup8c
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	metadata = pDataDT(merge.healthy.giotto)
	
	subset_cell_IDs = metadata[L1Subtype %in% c('Epithelial Basal 1', 'Epithelial Basal 2', 'Epithelial Suprabasal 1', 'Epithelial Suprabasal 2',
	                                            'HFSC', 'Upper HF Suprabasal 1', 'Lower HF Inflamm', 'Lower HF 1', 'Lower HF 2',
	                                            'Sebaceous', 'Neutrophil')]$cell_ID
	
	fig8c.giotto = subsetGiotto(merge.healthy.giotto, cell_ids = subset_cell_IDs)
	
	fig8c.meta <- pDataDT(fig8c.giotto) # pull existing metadata
	
	fig8c.meta$L1Subtype <- gsub(c('HFSC|Upper HF Suprabasal 1|Lower HF 1|Lower HF 2|Lower HF Inflamm'), 'Hair Follicle', fig8c.meta$L1Subtype) # combine HF populations into one cluster
	fig8c.meta$L1Subtype <- gsub(c('Epithelial Basal 1|Epithelial Basal 2|Epithelial Suprabasal 1|Epithelial Suprabasal 2'), 'Interfollicular Epi', fig8c.meta$L1Subtype) # combine IFE populations into one cluster
	
	fig8c.giotto <- addCellMetadata(fig8c.giotto,
	                                spat_unit = 'cell',
	                                feat_type = 'rna',
	                                fig8c.meta,
	                                by_column = F)
	
	unique(pDataDT(fig8c.giotto)$L1Subtype) # check if metadata correct
	
	marker.genes <- c('Il1rap', 'Il1b', 'Ccr1', 'Hdc', 'S100a8', 
	                  'Col17a1', 'Lor', 'Krt77', 'Krt78', 'Krt5',
	                  'Ptch1', 'Dsg4', 'Lgr5', 'Krt75', 'Sox9',
	                  'Elovl6', 'Cidea', 'Fabp5', 'Rora', 'Pparg')
	
	plotMetaDataHeatmap(fig8c.giotto, 
	                    custom_cluster_order = c('Neutrophil', 'Interfollicular Epi', 'Hair Follicle', 'Sebaceous'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes,
	                    custom_feat_order = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 8c Xenium HF L1Subtype Marker Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 8c Xenium HF L1Subtype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	# Sup8d
	
	tissue.neu <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/healthytissue.neu.rds')
	
	Idents(tissue.neu) <- 'L1subtype'
	
	tissue.neu <- subset(tissue.neu, idents = c('Neutrophil_1', 'Neutrophil_2'))
	
	tissue.neu$L1subtype <- droplevels(tissue.neu$L1subtype)
	tissue.neu$L2subtype <- droplevels(tissue.neu$L2subtype)
	
	tissue.neu <- RunUMAP(tissue.neu, 
	                      reduction = 'harmony', 
	                      reduction.name = 'umap.harmony',
	                      assay = 'SCT',
	                      dims = 1:70)
	
	Idents(tissue.neu) <- 'L1subtype'
	
	neu.L1subtype.markers = FindAllMarkers(tissue.neu,
	                                       only.pos = TRUE, 
	                                       test.use = "MAST",
	                                       latent.vars = "sex",
	                                       min.pct = 0.5,
	                                       logfc.threshold = 0.5,
	                                       assay = "RNA",
	                                       densify = TRUE) 
	
	neu.L1subtype.markers <- neu.L1subtype.markers[order(neu.L1subtype.markers$cluster,
	                                                     -neu.L1subtype.markers$avg_log2FC), ]
	
	write.xlsx(neu.L1subtype.markers,
	           file = '/data/overmilleram/Manuscript Figures/SupFig 8d Neutrophil L1 Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	cols = c('#CD3122FF','#E7A6CDFF')
	#Neutrophil_1, Neutrophil_2
	
	VlnPlot(tissue.neu,
	        pt.size = 0,
	        cols = c(rep('#CD3122FF', times = 4),
	                 rep('#E7A6CDFF', times = 4)),
	        features = c('Retnlg', 'Lcn2', 'S100a8', 'S100a9',
	                     'Cstb', 'Ccl3', 'Mif', 'Atox1'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 8d Neutrophil L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "SupFig 8d Neutrophil L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#######
	
	#Fig 5
	#######
	
	#5a/c
	multinichenet_celltype_output <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_celltype_output.rds')
	
	senders_receivers <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular', 'Neural')
	
	color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E15759', '#592B02FF')
	names(color.use) <- senders_receivers
	
	#top 40 skin_dox_8_week
	
	prioritized_tbl_oi_skdo_40 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                40, 
	                                                groups_oi = "skin_dox_8_week")
	
	circos_skdo = make_circos_one_group(prioritized_tbl_oi_skdo_40, color.use, color.use)
	
	#top 40 control_buccal_mucosa
	
	prioritized_tbl_oi_orbu_40 = get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 
	                                                40, 
	                                                groups_oi = "control_buccal_mucosa")
	
	circos_orbu = make_circos_one_group(prioritized_tbl_oi_orbu_40, color.use, color.use)
	
	
	circos_skdo[[1]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/Fig 5a Healthy Tissue Pitx1+ Top40 Celltype Circos.svg",
	          width = 5.6,
	          height = 5.6)
	dev.off()
	
	circos_orbu[[1]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/Fig 5c Healthy Tissue Buccal Top40 Celltype Circos.svg",
	          width = 5.6,
	          height = 5.6)
	dev.off()
	
	circos_orbu[[2]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/Fig 5ac Healthy Tissue Top40 Circos Legend.svg",
	          width = 2,
	          height = 4)
	dev.off()
	
	#######
	
	#Supplementary Fig 9
	#######
	multinichenet_celltype_output <- readRDS('/data/overmilleram/scRNAseq/Skin & Oral/MultiNicheNet/healthytissue_multinichenet_celltype_output.rds')
	
	group_oi = list("skin_dox_8_week","control_buccal_mucosa")
	
	prioritized_tbl_oi_100 = lapply(X = group_oi,
	                                FUN = function(x){
	                                  get_top_n_lr_pairs(multinichenet_celltype_output$prioritization_tables, 100, groups_oi = x)
	                                })
	
	plot_oi_100 = lapply(X = prioritized_tbl_oi_100,
	                     FUN = function(x){
	                       make_sample_lr_prod_activity_plots(multinichenet_celltype_output$prioritization_tables, x)
	                     })
	
	
	ggsave2(plot_oi_100[[1]],
	        filename = "SupFig 9 Top100 Dox LR Products and Ligand Activity Celltype.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 18,
	        units = "in",
	        dpi = 300,
	        limitsize = F)
	
	ggsave2(plot_oi_100[[2]],
	        filename = "SupFig 9 Top100 Buccal LR Products and Ligand Activity Celltype.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 15,
	        height = 18,
	        units = "in",
	        dpi = 300,
	        limitsize = F)
	
	#######
	
	#Fig 6
	#######
	
	# 6d
	
	color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E83800FF', '#805840FF', '#808898FF')
	names(color.use) <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular', 'Melanocyte', 'Skeletal_Muscle')
	
	Idents(wound.integrated) <- 'celltype'
	
	DimPlot(wound.integrated,
	        pt.size = 1,
	        label = F,
	        cols = color.use,
	        split.by = 'state',
	        ncol = 2,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 6d Wound Skin Celltype UMAP Split State.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "Fig 6d Wound Skin Celltype UMAP Split State.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 5,
	        dpi = 1200,
	        units = "in",
	        limitsize = F)
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	wound.integrated.test <- sc_utils(wound.integrated)
	
	prop_test <- permutation_test(wound.integrated.test, 
	                              cluster_identity = "celltype",
	                              sample_1 = "wound_control_8_week", 
	                              sample_2 = "wound_dox_8_week",
	                              sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4))
	
	table(wound.integrated$state, wound.integrated$celltype) #completed graph in Excel
	
	#######
	
	#Supplementary Fig 10
	#######
	
	#Sup10a
	
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	metadata = pDataDT(merge.healthy.giotto)
	
	subset_cell_IDs = metadata[L1Subtype %in% c('Epithelial Basal 1', 'Epithelial Basal 2', 'Epithelial Suprabasal 1', 'Epithelial Suprabasal 2',
	                                            'Lower HF 1', 'Lower HF 2', 'Lower HF Inflamm', 'HFSC', 'Upper HF Suprabasal 1', 'Sebaceous', 
	                                            'Dermal Fibroblast', 'Stromal Fibroblast', 'Inflammatory Fibroblast', 'Dermal Sheath', 'Dermal Papilla',
	                                            'T-cell', 'Plasma/B-cell', 'Macrophage', 'Dendritic', 'Langerhans', 'Neutrophil')]$cell_ID
	
	cell.giotto = subsetGiotto(merge.healthy.giotto, cell_ids = subset_cell_IDs)
	
	unique(pDataDT(cell.giotto)$L1Subtype) # check if metadata correct
	
	marker.genes <- c('Shh', 'Ptch1', 'Ptch2', 'Ccl3', 'Ccr1', 'Ccr5', 'Il1b', 'Il1r1', 'Il1rap')
	
	cluster.order <- c('Epithelial Basal 1', 'Epithelial Basal 2', 'Epithelial Suprabasal 1', 'Epithelial Suprabasal 2',
	                   'Lower HF 1', 'Lower HF 2', 'Lower HF Inflamm', 'HFSC', 'Upper HF Suprabasal 1', 'Sebaceous', 
	                   'Dermal Fibroblast', 'Stromal Fibroblast', 'Inflammatory Fibroblast', 'Dermal Sheath', 'Dermal Papilla',
	                   'T-cell', 'Plasma/B-cell', 'Macrophage', 'Dendritic', 'Langerhans', 'Neutrophil')
	
	plotMetaDataHeatmap(cell.giotto, 
	                    custom_cluster_order = cluster.order,
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes,
	                    custom_feat_order = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 10a Xenium Pitx1+ MultiNicheNet Celltype Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 10a Xenium Pitx1+ MultiNicheNet Celltype Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	#Sup10b
	
	subset_cell_IDs = metadata[L1Subtype %in% c('Epithelial Basal 1', 'Epithelial Basal 2', 'Epithelial Suprabasal 1', 'Epithelial Suprabasal 2',
	                                            'Lower HF 1', 'Lower HF 2', 'Lower HF Inflamm', 'HFSC', 'Upper HF Suprabasal 1', 'Sebaceous', 
	                                            'Oral Suprabasal 1', 'Oral Suprabasal 2', 'Oral Suprabasal 3', 'Oral Suprabasal 4',
	                                            'Dermal Fibroblast', 'Stromal Fibroblast', 'Inflammatory Fibroblast', 'Dermal Sheath', 'Dermal Papilla')]$cell_ID
	
	cell.giotto = subsetGiotto(merge.healthy.giotto, cell_ids = subset_cell_IDs)
	
	unique(pDataDT(cell.giotto)$L1Subtype) # check if metadata correct
	
	marker.genes <- c('Dsg1b', 'Dsg2', 'Dsc2', 'Dsc3', 'Dsg4', 'Wnt4', 'Wnt10a', 'Wnt16', 'Igf1', 'Itgb4')
	
	cluster.order <- c('Epithelial Basal 1', 'Epithelial Basal 2', 'Epithelial Suprabasal 1', 'Epithelial Suprabasal 2',
	                   'Lower HF 1', 'Lower HF 2', 'Lower HF Inflamm', 'HFSC', 'Upper HF Suprabasal 1', 'Sebaceous', 
	                   'Oral Suprabasal 1', 'Oral Suprabasal 2', 'Oral Suprabasal 3', 'Oral Suprabasal 4',
	                   'Dermal Fibroblast', 'Stromal Fibroblast', 'Inflammatory Fibroblast', 'Dermal Sheath', 'Dermal Papilla')
	
	plotMetaDataHeatmap(cell.giotto, 
	                    custom_cluster_order = cluster.order,
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes,
	                    custom_feat_order = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 10b Xenium Buccal MultiNicheNet Celltype Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 10b Xenium Buccal MultiNicheNet Celltype Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	#######
	
	#Supplementary Fig 11
	#######
	
	# Sup11a
	
	color.use = c('#7DB3E2FF','#59A14F', '#EED58CFF', '#E83800FF', '#805840FF', '#808898FF')
	names(color.use) <- c('Keratinocyte', 'Immune', 'Fibroblast', 'Vascular', 'Melanocyte', 'Skeletal_Muscle')
	
	Idents(wound.integrated) <- 'celltype'
	
	DimPlot(wound.integrated,
	        pt.size = 1,
	        label = F,
	        cols = color.use,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "SupFig 11a Wound Skin Celltype UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px",
	        limitsize = F)
	ggsave2(filename = "SupFig 11a Wound Skin Celltype UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in",
	        limitsize = F)
	
	# Sup11b
	
	Idents(wound.integrated) <- 'sampleid'
	
	VlnPlot(wound.integrated,
	        pt.size = 0,
	        features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
	        split.by = 'sampleid',
	        cols = paletteer::paletteer_d("ggsci::default_ucscgb"),
	        stack = T,
	        flip = T) + NoLegend()
	
	ggsave2(filename = "SupFig 11b Wound Skin Quality Control VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	
	ggsave2(filename = "SupFig 11b Wound Skin Quality Control VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	# Sup11c
	
	celltype.cells.by.sample <- data.frame(table(wound.integrated$sampleid, wound.integrated$celltype))
	
	colnames(celltype.cells.by.sample) <- c('sampleid', 'celltype', 'cells') # completed in Excel
	
	write.xlsx(celltype.cells.by.sample,
	           '/data/overmilleram/Manuscript Figures/Sup Fig 11b Cell Counts table.xlsx')
	
	#######
	
	#Supplementary Fig 12
	#######
	
	# Sup12a
	
	VlnPlot(wound.integrated,
	        pt.size = 0,
	        cols = c(rep("#7DB3E2FF", times = 2),
	                 rep("#59A14F", times = 2),
	                 rep("#EED58CFF", times = 2), 
	                 rep("#E83800FF", times = 2),
	                 rep("#805840FF", times = 2),
	                 rep('#808898FF', times = 2)),
	        features = c('Krt5', 'Perp', 
	                     'Ptprc', 'Fcer1g', 
	                     'Col1a1', 'Dcn',
	                     'Pecam1', 'Cdh5',
	                     'Mlana', 'Pmel',
	                     'Tnnt3', 'Des'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 12a Celltype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 12a Celltype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 1200,
	        units = "in")
	
	# Sup12b
	
	# load in Giotto objects
	
	wound.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Object/')
	
	wound.markers = findMarkers_one_vs_all(gobject = wound.giotto,
	                                       method = 'gini',
	                                       expression_values = 'normalized',
	                                       cluster_column = 'Celltype',
	                                       rank_score = 2)
	
	write.xlsx(wound.markers,
	           file = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Celltype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	wound.markers <- readxl::read_xlsx('/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Celltype Markers.xl', 
	                                   sheet = 1,
	                                   col_names = T)
	wound.markers[ ,1] <- NULL
	wound.markers <- wound.markers[order(wound.markers$cluster,
	                                     -wound.markers$expression), ]
	
	# plot top cell markers in heatmap
	
	marker.genes <- wound.markers %>% slice_max(expression_gini, n=5, by=cluster)
	marker.genes <- unique(marker.genes[ ,1]) %>% unlist
	
	plotMetaDataHeatmap(wound.giotto, 
	                    custom_cluster_order = c('Mesenchymal', 'Vascular', 'Immune', 'Fibroblast', 'Neural', 'Keratinocyte'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'Celltype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 12b Wound Skin Xenium Celltype Marker Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 12b Xenium Wound Skin Celltype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	# Sup12c
	
	celltypes <- c('keratinocyte', 'immune', 'fibroblast', 'vascular', 'neural', 'mesenchymal')
	xenium.control <- c(5210, 4059, 2643, 728, 175, 962) 
	xenium.pitx1 <- c(25217, 42829, 19395, 4106, 1050, 3785)
	
	xenium.prop.cvp <- rbind(xenium.control, xenium.pitx1)
	
	colnames(xenium.prop.cvp) <- celltypes
	
	#control vs pitx1+
	#p-value (all are >>>0.05)
	prop.test(c(xenium.prop.cvp[1,1], xenium.prop.cvp[2,1]), c(sum(xenium.prop.cvp[ ,1]), sum(xenium.prop.cvp[ ,1])), cor = F, conf.level = 0.975)$p.value #keratinocyte
	prop.test(c(xenium.prop.cvp[1,2], xenium.prop.cvp[2,2]), c(sum(xenium.prop.cvp[ ,2]), sum(xenium.prop.cvp[ ,2])), cor = F, conf.level = 0.975)$p.value #immune
	prop.test(c(xenium.prop.cvp[1,3], xenium.prop.cvp[2,3]), c(sum(xenium.prop.cvp[ ,3]), sum(xenium.prop.cvp[ ,3])), cor = F, conf.level = 0.975)$p.value #fibroblast
	prop.test(c(xenium.prop.cvp[1,4], xenium.prop.cvp[2,4]), c(sum(xenium.prop.cvp[ ,4]), sum(xenium.prop.cvp[ ,4])), cor = F, conf.level = 0.975)$p.value #vascular
	prop.test(c(xenium.prop.cvp[1,5], xenium.prop.cvp[2,5]), c(sum(xenium.prop.cvp[ ,5]), sum(xenium.prop.cvp[ ,5])), cor = F, conf.level = 0.975)$p.value #neural
	prop.test(c(xenium.prop.cvp[1,6], xenium.prop.cvp[2,6]), c(sum(xenium.prop.cvp[ ,6]), sum(xenium.prop.cvp[ ,6])), cor = F, conf.level = 0.975)$p.value #mesenchymal
	
	#log2fc
	log2((xenium.prop.cvp[2,1]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,1]/sum(xenium.prop.cvp[1,]))) #keratinocyte, -0.5314
	log2((xenium.prop.cvp[2,2]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,2]/sum(xenium.prop.cvp[1,]))) #immune,        0.5929
	log2((xenium.prop.cvp[2,3]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,3]/sum(xenium.prop.cvp[1,]))) #fibroblast,    0.0689
	log2((xenium.prop.cvp[2,4]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,4]/sum(xenium.prop.cvp[1,]))) #vascular,     -0.3108
	log2((xenium.prop.cvp[2,5]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,5]/sum(xenium.prop.cvp[1,]))) #neural,       -0.2215
	log2((xenium.prop.cvp[2,6]/sum(xenium.prop.cvp[2,]))/(xenium.prop.cvp[1,6]/sum(xenium.prop.cvp[1,]))) #mesenchymal,  -0.8303
	
	#######
	
	#Fig6
	#######
	
	# 6a
	
	wound.kera <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.kera.rds')
	
	cols = c('#0142FEFF','#1CADE4FF','#324376FF', '#219089FF', '#486048FF',
	         '#B0D8D0FF', '#D0F0F0BA', '#6088A0FF', '#787850FF', '#F88080FF', '#B04030FF')
	
	# Epithelial_Basal, Epithelial_Suprabasal, Proliferating_Keratinocyte, Upper_HF_Suprabasal, Sebaceous, 
	# HFSC_1, HFSC_2, Lower_HF_1, Lower_HF_2, Wound_Activated, Wound_Migratory  
	
	Idents(wound.kera) <- 'L1subtype'
	
	DimPlot(wound.kera,
	        split.by = 'state',
	        cols = cols,
	        ncol = 2,
	        pt.size = 2,
	        label = F,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 6a Wound Keratinocyte L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "Fig 6a Wound Keratinocyte L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	skin.kera.test = sc_utils(wound.kera)
	
	prop_test = permutation_test(skin.kera.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "wound_control_8_week", 
	                             sample_2 = "wound_dox_8_week",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4)) 
	
	table(wound.kera$state, wound.kera$L1subtype) #completed in Excel
	
	# 6c
	
	wound.imm <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds')
	
	cols = c('#FCD47CFF', '#F57206FF', '#4838A8FF', '#97AD3DFF','#15CC31FF', '#415521FF',
	         '#CCF4CFFF', '#009CCEFF', '#0F6A81DA',  '#75FB8AFF','#7FD4C1FF', '#CD3122FF','#E7A6CDFF')
	#'T_cell', 'NK_Cell', 'pDC', 'Macrophage_1', 'Macrophage_2', 'Macrophage_3' 
	#'Dendritic_1', 'Dendritic_2', 'Dendritic_3', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2'
	
	Idents(wound.imm) <- 'L1subtype'
	
	DimPlot(wound.imm,
	        split.by = 'state',
	        cols = cols,
	        ncol = 2,
	        pt.size = 2,
	        label = F,
	        raster = F) + NoLegend()
	
	ggsave2(filename = "Fig 6c Wound Immune L1 Subtype Split State UMAP.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "Fig 6c Wound Immune L1 Subtype Split State UMAP.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	library(scProportionTest) #check if proportions of cells are statistically different between Control & +Dox Skin
	
	skin.imm.test = sc_utils(wound.imm)
	
	prop_test = permutation_test(skin.imm.test, 
	                             cluster_identity = "L1subtype",
	                             sample_1 = "wound_control_8_week", 
	                             sample_2 = "wound_dox_8_week",
	                             sample_identity = "state")
	
	permutation_plot(prop_test,
	                 log2FD_threshold = log2(4)) 
	
	table(wound.imm$state, wound.imm$L1subtype) #completed in Excel
	
	#######
	
	#Supplementary Fig 13
	#######
	# Sup13a
	
	wound.kera <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.kera.rds')
	
	Idents(wound.kera) <- 'L1subtype'
	
	woundkera.L1subtype.markers = FindAllMarkers(wound.kera,
	                                             only.pos = TRUE, 
	                                             test.use = "poisson",
	                                             latent.vars = "sex",
	                                             min.pct = 0.50,
	                                             logfc.threshold = 2,
	                                             assay = "RNA",
	                                             densify = TRUE) 
	
	woundkera.L1subtype.markers <- woundkera.L1subtype.markers[order(woundkera.L1subtype.markers$cluster,
	                                                                 -woundkera.L1subtype.markers$avg_log2FC), ]
	
	write.xlsx(woundkera.L1subtype.markers,
	           file = '/data/overmilleram/Manuscript Figures/SupFig 13a Wound Keratinocyte L1 Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	Wound_HFSC1_markers <- FindMarkers(wound.kera,
	                                   ident.1 = 'HFSC_1',
	                                   only.pos = T,
	                                   assay = 'RNA',
	                                   test.use = 'MAST',
	                                   latent.vars = 'sex',
	                                   min.pct = 0.25,
	                                   logfc.threshold = 1,
	                                   densify = T)
	
	cols = c('#0142FEFF','#1CADE4FF','#324376FF', '#219089FF', '#486048FF',
	         '#B0D8D0FF', '#D0F0F0BA', '#6088A0FF', '#787850FF', '#F88080FF', '#B04030FF')
	
	# Epithelial_Basal, Epithelial_Suprabasal, Proliferating_Keratinocyte, Upper_HF_Suprabasal, Sebaceous, 
	# HFSC_1, HFSC_2, Lower_HF_1, Lower_HF_2, Wound_Activated, Wound_Migratory  
	
	Idents(wound.kera) <- 'L1subtype'
	
	DimPlot(wound.kera,
	        cols = cols,
	        pt.size = 2,
	        label = T,
	        raster = F) + NoLegend()
	
	VlnPlot(wound.kera,
	        pt.size = 0,
	        cols = c(rep('#0142FEFF', times = 2),
	                 rep('#1CADE4FF', times = 2),
	                 rep('#324376FF', times = 2),
	                 rep('#219089FF', times = 2),
	                 rep('#486048FF', times = 2),
	                 rep('#B0D8D0FF', times = 2),
	                 rep('#D0F0F0BA', times = 2),
	                 rep('#6088A0FF', times = 2),
	                 rep('#787850FF', times = 2),
	                 rep('#F88080FF', times = 2),
	                 rep('#B04030FF', times = 2)),
	        features = c('Krt14', 'Col23a1', 
	                     'Krt10', 'Serpinb3b', 
	                     'Hist1h2ap', 'Cenpf', 
	                     'Cst6', 'Defb6', 
	                     'Scd1', 'Mgst1', 
	                     'Angptl4', 'Fst',
	                     'Cd34', 'Tgm5',
	                     'Sox5', 'Ptch1', 
	                     'Lef1', 'Pik3r3',
	                     'Stfa1', 'Sprr2a3', 
	                     'Igfbp3', 'Itga2'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 13a Wound Keratinocyte L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 13a Wound Keratinocyte L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	
	# Sup12b
	
	wound.imm <- readRDS('/data/overmilleram/scRNAseq/Wound/woundskin.imm.rds')
	
	Idents(wound.imm) <- 'L1subtype'
	
	woundimm.L1subtype.markers = FindAllMarkers(wound.imm,
	                                            only.pos = TRUE, 
	                                            test.use = "poisson",
	                                            latent.vars = "sex",
	                                            min.pct = 0.50,
	                                            logfc.threshold = 2,
	                                            assay = "RNA",
	                                            densify = TRUE) 
	
	woundimm.L1subtype.markers <- woundimm.L1subtype.markers[order(woundimm.L1subtype.markers$cluster,
	                                                               -woundimm.L1subtype.markers$avg_log2FC), ]
	
	write.xlsx(woundimm.L1subtype.markers,
	           file = '/data/overmilleram/Manuscript Figures/SupFig 13a Wound Immune L1 Subtype Markers.xlsx',
	           sheetName = "Cluster Markers",
	           append = FALSE)
	
	Wound_Mac3_markers <- FindMarkers(wound.imm,
	                                  ident.1 = 'Macrophage_3',
	                                  only.pos = T,
	                                  assay = 'RNA',
	                                  test.use = 'MAST',
	                                  latent.vars = 'sex',
	                                  min.pct = 0.25,
	                                  logfc.threshold = 1,
	                                  densify = T)
	
	cols = c('#FCD47CFF', '#F57206FF', '#4838A8FF', '#97AD3DFF','#15CC31FF', '#415521FF',
	         '#CCF4CFFF', '#009CCEFF', '#0F6A81DA',  '#75FB8AFF','#7FD4C1FF', '#CD3122FF','#E7A6CDFF')
	#'T_cell', 'NK_Cell', 'pDC', 'Macrophage_1', 'Macrophage_2', 'Macrophage_3' 
	#'Dendritic_1', 'Dendritic_2', 'Dendritic_3', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2'
	
	VlnPlot(wound.imm,
	        pt.size = 0,
	        cols = c(rep('#FCD47CFF', times = 2),
	                 rep('#F57206FF', times = 2),
	                 rep('#4838A8FF', times = 2),
	                 rep('#97AD3DFF', times = 2),
	                 rep('#15CC31FF', times = 2),
	                 rep('#415521FF', times = 2),
	                 rep('#CCF4CFFF', times = 2),
	                 rep('#009CCEFF', times = 2),
	                 rep('#0F6A81DA', times = 2),
	                 rep('#75FB8AFF', times = 2),
	                 rep('#7FD4C1FF', times = 2),
	                 rep('#CD3122FF', times = 2),
	                 rep('#E7A6CDFF', times = 2)),
	        features = c('Icos', 'St6galnac3',
	                     'Trdc', 'Xcl1',
	                     'Siglech', 'Ccr9', 
	                     'Ccl8', 'Cbr2',
	                     'Arg1', 'Fn1',
	                     'Cd36', 'Spp1',
	                     'H2-Ab1', 'Tmem176b', 
	                     'Cacnb3', 'Ccl22',
	                     'Cxcl10', 'Ifit1', 
	                     'Cd207', 'Mfge8',  
	                     'Plac8', 'Vcan',
	                     'Retnlg', 'G0s2',
	                     'S100a8', 'Ccl4'),
	        stack = T,
	        flip = T,
	        assay = 'RNA') + NoLegend()
	
	ggsave2(filename = "SupFig 13b Wound Immune L1 Subtype Marker VlnPlot.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 13b Wound Immune L1 Subtype Marker VlnPlot.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 10,
	        height = 10,
	        dpi = 600,
	        units = "in")
	
	#Sup12d
	wound.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Wound Skin Giotto Object/')
	merge.healthy.giotto <- loadGiotto(path_to_folder = '/data/overmilleram/Xenium/Giotto/Healthy Skin & Buccal Merged Giotto Object/')
	
	metadata = pDataDT(wound.giotto)
	
	subset_cell_IDs = metadata[L1Subtype %in% c('Activated HF', 'Activated Suprabasal 1', 'Activated Suprabasal 2',
	                                            'Epithelial Basal', 'Epithelial Suprabasal')]$cell_ID
	
	fig11d.wound.giotto = subsetGiotto(wound.giotto, cell_ids = subset_cell_IDs)
	
	marker.genes <- c('Krt6a', 'Krt16', 'S100a8', 'Dsc2', 'Aldh1a3', 'Il1b', 'Il1rn', 'Slpi',
	                  'Krt1', 'Krt4', 'Tgm3', 'Flg', 'Hrnr', 'Cnfn', 'Col17a1', 'Hist1h1b', 'Krt77', 'Egfr')
	
	plotMetaDataHeatmap(fig11d.wound.giotto, 
	                    custom_cluster_order = c('Activated HF', 'Activated Suprabasal 1', 'Activated Suprabasal 2',
	                                             'Epithelial Basal', 'Epithelial Suprabasal'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes,
	                    custom_feat_order = marker.genes) + theme_minimal() 
	
	ggsave2(filename = "SupFig 13d Xenium Wound L1Subtype Marker Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 13d Xenium Wound L1Subtype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	metadata = pDataDT(merge.healthy.giotto)
	
	subset_cell_IDs = metadata[L1Subtype %in% c('Oral Suprabasal 1', 'Oral Suprabasal 2', 'Oral Suprabasal 3', 'Oral Suprabasal 4')]$cell_ID
	
	fig11d.buccal.giotto = subsetGiotto(merge.healthy.giotto, cell_ids = subset_cell_IDs)
	
	marker.genes <- c('Dsg2', 'Col17a1', 'Krt6a', 'Krt5',
	                  'Mt3', 'Dsg3', 'Barx2', 'Krt4',
	                  'Dsg1b', 'Dsc1', 'Klf4', 'Krt16', 
	                  'Sprr3', 'Lor', 'Hrnr', 'Cidea')
	
	plotMetaDataHeatmap(fig11d.buccal.giotto, 
	                    custom_cluster_order = c('Oral Suprabasal 1', 'Oral Suprabasal 2', 'Oral Suprabasal 3', 'Oral Suprabasal 4'),
	                    expression_values = 'scaled',
	                    metadata_cols = 'L1Subtype',
	                    gradient_color = c('#172869FF','#FFFFFFFF','#A50021FF'),
	                    selected_feats = marker.genes,
	                    custom_feat_order = marker.genes) + theme_minimal()
	
	ggsave2(filename = "SupFig 13d Xenium Buccal L1Subtype Marker Heatmap.svg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5000,
	        height = 5000,
	        units = "px")
	ggsave2(filename = "SupFig 13d Xenium Buccal L1Subtype Marker Heatmap.jpeg",
	        path = "/data/overmilleram/Manuscript Figures/",
	        width = 5,
	        height = 5,
	        dpi = 600,
	        units = "in")
	
	# Sup12e
	
	multinichenet_woundskin_L1subtype_output <- readRDS('/data/overmilleram/scRNAseq/Wound/MultiNicheNet/woundskin_multinichenet_L1subtype_output.rds')
	
	contrasts_oi <- c("'wound_control_8_week-wound_dox_8_week','wound_dox_8_week-wound_control_8_week'")
	
	contrast_tbl <- tibble(contrast = c("wound_control_8_week-wound_dox_8_week","wound_dox_8_week-wound_control_8_week"),
	                       group = c("wound_control_8_week","wound_dox_8_week"))
	
	senders_receivers <- unique(levels(wound.integrated$L1subtype))
	senders_receivers <- senders_receivers[-c(33:35)] # remove muscle cells and melanocytes
	
	color.use <- c('#0142FEFF','#1CADE4FF','#324376FF', '#219089FF','#486048FF',
	               '#B0D8D0FF','#D0F0F0BA','#6088A0FF', '#787850FF','#CC7A88FF', '#B04030FF',
	               '#FFFF80FF','#FFD700FF','#73652DFF', '#C4A000FF','#E5B17EFF', '#99540FFF',
	               '#FCD47CFF','#F57206FF','#4838A8FF', '#97AD3DFF','#15CC31FF', '#415521FF',
	               '#CCF4CFFF','#009CCEFF','#0F6A81DA', '#75FB8AFF','#7FD4C1FF', '#CD3122FF','#F6C4E5FF',
	               '#FF0010FF','#FFB79FFF')
	
	# Epithelial_Basal, Epithelial_Suprabasal, Proliferating_Keratinocyte, Upper_HF_Suprabasal, Sebaceous, 
	# HFSC_1, HFSC_2, Lower_HF_1, Lower_HF_2, Wound_Activated, Wound_Migratory
	# Dermal_Fibroblast_1, Dermal_Fibroblast_2, Dermal_Sheath, Dermal_Papilla, Myofibroblast_1, Myofibroblast_2
	#'T_cell', 'NK_Cell', 'pDC', 'Macrophage_1', 'Macrophage_2', 'Macrophage_3' 
	#'Dendritic_1', 'Dendritic_2', 'Dendritic_3', 'Langerhans', 'Monocyte', 'Neutrophil_1', 'Neutrophil_2'
	#'Endothelial, 'Lymph_vessel'
	
	names(color.use) <- senders_receivers
	
	# comparison multinichenet
	
	prioritized_tbl_oi_imm_wound = get_top_n_lr_pairs(multinichenet_woundskin_L1subtype_output$prioritization_tables,
	                                                  30,
	                                                  groups_oi = 'wound_dox_8_week',
	                                                  senders_oi = c('Neutrophil_1', 'Neutrophil_2'))
	
	circos_imm = make_circos_one_group(prioritized_tbl_oi_imm_wound, color.use, color.use)
	
	circos_imm[[1]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/SupFig 13e Wound Skin Pitx1+ Top50 Neutrophil-L1subtype Circos.svg",
	          width = 5.6,
	          height = 5.6)
	dev.off()
	
	circos_imm[[2]]
	dev.print(svg,
	          "/data/overmilleram/Manuscript Figures/SupFig 13e Wound Skin Pitx1+ Top50 Neutrophil-L1subtype Circos Legend.svg",
	          width = 5,
	          height = 10)
	dev.off()


