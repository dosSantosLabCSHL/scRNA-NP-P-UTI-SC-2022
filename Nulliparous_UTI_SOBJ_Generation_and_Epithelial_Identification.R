######
# Preprocessing of scRNA-seq data from 
# No-UTI Mice (Bioproject : PRJNA677888) - 2 Samples
# UTI Mice (Bioproject : PRJNA855880) - 2 Samples
#Steps Included:
# 1. Load in Data
######

#Required Libraries
library(dplyr)
library(Seurat)
library(clustree)
library(ggplot2)

#1. Load in Libraries 
#Nulliparous UTI Infected Samples
UTI_sobj_1 <- Read10X(data.dir ="filepathway/Null_UTI_S1/filtered_feature_bc_matrix")
UTI_sobj_1 <-CreateSeuratObject(UTI_sobj_1, 
                                project = "UTI_sobj_1", 
                                min.cells = 3)
UTI_sobj_1 #View the number of cells and features in your sobj before QC

UTI_sobj_2 <- Read10X(data.dir ="filepathway/Null_UTI_S2/filtered_feature_bc_matrix")
UTI_sobj_2 <-CreateSeuratObject(UTI_sobj_2, 
                                project = "UTI_sobj_2", 
                                min.cells = 3)
UTI_sobj_2 #View the number of cells and features in your sobj before QC

#Nulliparous NO-UTI Infected Samples
NoUTI_sobj_1 <- Read10X(data.dir ="filepathway/Nulliparous_WT_S1/filtered_feature_bc_matrix")
NoUTI_sobj_1 <-CreateSeuratObject(NoUTI_sobj_1, 
                                  project = "NoUTI_sobj_1", 
                                  min.cells = 3)
NoUTI_sobj_1 #View the number of cells and features in your sobj before QC

NoUTI_sobj_2 <- Read10X(data.dir ="filepathway/Nulliparous_WT_S2/filtered_feature_bc_matrix")
NoUTI_sobj_2 <-CreateSeuratObject(NoUTI_sobj_2, 
                                  project = "NoUTI_sobj_2", 
                                  min.cells = 3)
NoUTI_sobj_2 #View the number of cells and features in your sobj before QC

#2. Run QC
#Define %MT content
UTI_sobj_1[["percent.mt"]] <- PercentageFeatureSet(UTI_sobj_1,
                                                   pattern = c("^mt-"))
UTI_sobj_2[["percent.mt"]] <- PercentageFeatureSet(UTI_sobj_2, 
                                                   pattern = c("^mt-"))
NoUTI_sobj_1[["percent.mt"]] <- PercentageFeatureSet(NoUTI_sobj_1, 
                                                     pattern = c("^mt-"))
NoUTI_sobj_2[["percent.mt"]] <- PercentageFeatureSet(NoUTI_sobj_2, 
                                                     pattern = c("^mt-"))

#Conduct QC in UTI_sobj_1
#Visualize the number of genes, the amount of RNA and the percentage of mitochodrial genes in each cell
VlnPlot(UTI_sobj_1, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#Determine the metrics of these features that would be considered outliers and adjust the line below accordingly
UTI_sobj_1 <- subset(UTI_sobj_1, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
#Examine these features in Vlnplots again to examine appropriate cells were removed
VlnPlot(UTI_sobj_1, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#View sobj to examine how the number of cells and features changed after following QC protocol
UTI_sobj_1

#Conduct QC in UTI_sobj_2
#Visualize the number of genes, the amount of RNA and the percentage of mitochodrial genes in each cell
VlnPlot(UTI_sobj_2, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#Determine the metrics of these features that would be considered outliers and adjust the line below accordingly
UTI_sobj_2 <- subset(UTI_sobj_2, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
#Examine these features in Vlnplots again to examine appropriate cells were removed
VlnPlot(UTI_sobj_2, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#View sobj to examine how the number of cells and features changed after following QC protocol
UTI_sobj_2

#Conduct QC in NoUTI_sobj_1
#Visualize the number of genes, the amount of RNA and the percentage of mitochodrial genes in each cell
VlnPlot(NoUTI_sobj_1, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#Determine the metrics of these features that would be considered outliers and adjust the line below accordingly
NoUTI_sobj_1 <- subset(NoUTI_sobj_1, 
                       subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
#Examine these features in Vlnplots again to examine appropriate cells were removed
VlnPlot(NoUTI_sobj_1, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#View sobj to examine how the number of cells and features changed after following QC protocol
NoUTI_sobj_1

#Conduct QC in NoUTI_sobj_2
#Visualize the number of genes, the amount of RNA and the percentage of mitochodrial genes in each cell
VlnPlot(NoUTI_sobj_2, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#Determine the metrics of these features that would be considered outliers and adjust the line below accordingly
NoUTI_sobj_2 <- subset(NoUTI_sobj_2, 
                       subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
#Examine these features in Vlnplots again to examine appropriate cells were removed
VlnPlot(NoUTI_sobj_2, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"),
        group.by = 'orig.ident')
#View sobj to examine how the number of cells and features changed after following QC protocol
NoUTI_sobj_2

#3. Merge and Integrate SOBJs
#Merge the sobjs
Nulliparous_UTI_sobj.merge <-merge(UTI_sobj_1, 
                                   c(UTI_sobj_2, 
                                     NoUTI_sobj_1,
                                     NoUTI_sobj_2))

#Split merged object by sample - Typically stored in orig.ident
Nulliparous_UTI_sobj.list <- SplitObject(Nulliparous_UTI_sobj.merge, 
                                         split.by = "orig.ident")

#Normalize and Find Highly Variable Features
for (i in 1:length(Nulliparous_UTI_sobj.list)) {
  Nulliparous_UTI_sobj.list[[i]] <- NormalizeData(Nulliparous_UTI_sobj.list[[i]], 
                                                  verbose = FALSE)
  Nulliparous_UTI_sobj.list[[i]] <- FindVariableFeatures(Nulliparous_UTI_sobj.list[[i]], 
                                                         selection.method = "vst",
                                                         nfeatures = 2000, 
                                                         verbose = FALSE)
}

#Find anchors and integrate data
reference.list <- Nulliparous_UTI_sobj.list[c("UTI_sobj_1", 
                                              "UTI_sobj_2",
                                              "NoUTI_sobj_1",
                                              "NoUTI_sobj_2")]
Nulliparous_UTI_sobj.anchors <- FindIntegrationAnchors(object.list = reference.list , 
                                                       dims = 1:30)
Nulliparous_UTI_sobj.integrated <- IntegrateData(anchorset = Nulliparous_UTI_sobj.anchors, 
                                                 dims = 1:30)

#Scale Data
DefaultAssay(Nulliparous_UTI_sobj.integrated) <- "integrated"
Nulliparous_UTI_sobj.integrated <- ScaleData(Nulliparous_UTI_sobj.integrated)

#4.Perform Linear Regression
#Run PCA
Nulliparous_UTI_sobj.integrated <- RunPCA(Nulliparous_UTI_sobj.integrated,
                                          features = VariableFeatures(object = Nulliparous_UTI_sobj.integrated))
#Visualize PCA reduction
#DimPlot(Nulliparous_UTI_sobj.integrated, 
#        reduction = "pca")

#Determine Dimensionality of the dataset
ElbowPlot(Nulliparous_UTI_sobj.integrated) #Choose Value were the value on the y axis remains about constant

#Determine Resolution of the dataset
Nulliparous_UTI_sobj.integrated <- FindNeighbors(Nulliparous_UTI_sobj.integrated, 
                                                 dims = 1:10) #change depending on the dimensionality used
Nulliparous_UTI_sobj.integrated <- FindClusters(Nulliparous_UTI_sobj.integrated, 
                                                resolution = 0.2)
Nulliparous_UTI_sobj.integrated <- FindClusters(Nulliparous_UTI_sobj.integrated, 
                                                resolution = 0.4)
Nulliparous_UTI_sobj.integrated <- FindClusters(Nulliparous_UTI_sobj.integrated, 
                                                resolution = 0.5)
Nulliparous_UTI_sobj.integrated <- FindClusters(Nulliparous_UTI_sobj.integrated, 
                                                resolution = 0.8)
Nulliparous_UTI_sobj.integrated <- FindClusters(Nulliparous_UTI_sobj.integrated, 
                                                resolution = 1)
clustree(Nulliparous_UTI_sobj.integrated, prefix = "integrated_snn_res.")

#Cluster Cells: Adjust Parameters in the following functions based on elbowplot and clustree
Nulliparous_UTI_sobj.integrated <- FindNeighbors(Nulliparous_UTI_sobj.integrated, 
                                                 dims = 1:10)
Nulliparous_UTI_sobj.integrated <- FindClusters(Nulliparous_UTI_sobj.integrated,
                                                resolution = 0.5)

#Visualize Clusters Using UMAP
Nulliparous_UTI_sobj.integrated <- RunUMAP(Nulliparous_UTI_sobj.integrated, 
                                           dims = 1:10) #changes dims as needed
DimPlot(Nulliparous_UTI_sobj.integrated, reduction = "umap")

#Save Seurat Object
saveRDS(Nulliparous_UTI_sobj.integrated, file = "filepathway/Nulliparous_UTI_sobj.integrated.rds")

#6. Define Epithelial Cells 
#Read in File if Neccessary 
Nulliparous_UTI_sobj.integrated <- readRDS("filepathway/Nulliparous_UTI_sobj.integrated.rds")

#Load in Gene Markers : Derived from DOI: https://doi.org/10.1007/s10911-021-09486-3
SC_Identification_Markers <- read_excel("filepathway/SC_Identification_Markers.xlsx")

#Make sure in RNA
DefaultAssay(Nulliparous_UTI_sobj.integrated) <- "RNA"
Nulliparous_UTI_sobj.integrated <- ScaleData(Nulliparous_UTI_sobj.integrated)

#Visualize Markers Across Clusters
DotPlot(Nulliparous_UTI_sobj.integrated, 
        features = SC_Identification_Markers) + RotatedAxis()

#Examine Specific Genes in FP
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Epcam")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Cd3e")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Ms4a1")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Cd19")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Tyrobp")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Sparc")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Fabp4")
FeaturePlot(Nulliparous_UTI_sobj.integrated, 
            features = "Vcam1")

#Subset Epithelial Cells
Nulliparous_UTI_Epith <- subset(Nulliparous_UTI_sobj.integrated, 
                                idents = c("2", "4", "5", "7", "10", "16"))
#Examine proper clusters were subsetted
DimPlot(Nulliparous_UTI_Epith, 
        reduction = "umap", 
        label = TRUE)

#7. Recluster Epithelial Cell Populations 
#make sure in RNA assay for initial step
DefaultAssay(Nulliparous_UTI_Epith) <- "RNA"
ScaleData(Nulliparous_UTI_Epith)

#Find Highly Variable Features
Nulliparous_UTI_Epith <- FindVariableFeatures(Nulliparous_UTI_Epith, 
                                              selection.method = "vst", 
                                              nfeatures = 2000)

#return to integrated assay
all.genes <- rownames(Nulliparous_UTI_Epith)
DefaultAssay(Nulliparous_UTI_Epith) <- "integrated"
Nulliparous_UTI_Epith <- ScaleData(Nulliparous_UTI_Epith, 
                                   features = all.genes)

#Cluster Cells
Nulliparous_UTI_Epith <- RunPCA(Nulliparous_UTI_Epith, 
                                features = VariableFeatures(object = Nulliparous_UTI_Epith))
#Determine dimensionality
ElbowPlot(Nulliparous_UTI_Epith)
#Determine Resolution
Nulliparous_UTI_Epith <- FindNeighbors(Nulliparous_UTI_Epith, 
                                       dims = 1:10)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.2)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.4)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.5)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.8)
clustree(Nulliparous_UTI_Epith, prefix = "integrated_snn_res.")
#Visualize Clusters in UMAP
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.5)
Nulliparous_UTI_Epith <- RunUMAP(Nulliparous_UTI_Epith, 
                                 reduction = 'pca', 
                                 dims = 1:10) 
DimPlot(Nulliparous_UTI_Epith, 
        reduction = "umap", 
        label = TRUE)

#Save First Epithelail Re-cluster
saveRDS(Nulliparous_UTI_Epith, file = "filepathway/Nulliparous_UTI_Epith_R1.rds") 
#Examine Quality of Clusters 
#If Working with an Integrted Sobj Make sure in RNA Assay
DefaultAssay(Nulliparous_UTI_Epith) <- "RNA"
Nulliparous_UTI_Epith<- ScaleData(Nulliparous_UTI_Epith)

#Determine Percent Ribosomal Genes
Nulliparous_UTI_Epith[["percent.ribo"]] <- PercentageFeatureSet(Nulliparous_UTI_Epith, 
                                                                pattern = c("^Rp[sl]"))

#Visualize Features and Markers
VlnPlot(Nulliparous_UTI_Epith, features = c("nFeature_RNA", "percent.mt", "percent.ribo"))
DotPlot(Nulliparous_UTI_Epith, features = SC_Identification_Markers) +RotatedAxis()

#8. Second Epithelial Re-clustering 
#subset high quality clustes (removed those of low quality or of non-epithelail origin)
Nulliparous_UTI_Epith <- subset(Nulliparous_UTI_Epith, 
                                idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "10"))

DimPlot(Nulliparous_UTI_Epith, 
        reduction = "umap", 
        label = TRUE)

#If working with an Integrated sobj, make sure in RNA assay for initial step
DefaultAssay(Nulliparous_UTI_Epith) <- "RNA"
ScaleData(Nulliparous_UTI_Epith)

#Find Highly Variable Features
Nulliparous_UTI_Epith <- FindVariableFeatures(Nulliparous_UTI_Epith, 
                                              selection.method = "vst", 
                                              nfeatures = 2000)

#If working with an Integrated sobj, return to integrated assay
all.genes <- rownames(Nulliparous_UTI_Epith)
DefaultAssay(Nulliparous_UTI_Epith) <- "integrated"
Nulliparous_UTI_Epith <- ScaleData(Nulliparous_UTI_Epith, 
                                   features = all.genes)

#Cluster Cells
Nulliparous_UTI_Epith <- RunPCA(Nulliparous_UTI_Epith, 
                                features = VariableFeatures(object = Nulliparous_UTI_Epith))
#Determine dimensionality
ElbowPlot(Nulliparous_UTI_Epith)
#Determine Resolution
Nulliparous_UTI_Epith <- FindNeighbors(Nulliparous_UTI_Epith, 
                                       dims = 1:8)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.2)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.4)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.5)
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.8)
clustree(Nulliparous_UTI_Epith, prefix = "integrated_snn_res.")
#Visualize Clusters in UMAP
Nulliparous_UTI_Epith <- FindClusters(Nulliparous_UTI_Epith, 
                                      resolution = 0.5)
Nulliparous_UTI_Epith <- RunUMAP(Nulliparous_UTI_Epith, 
                                 reduction = 'pca', 
                                 dims = 1:8) 
DimPlot(Nulliparous_UTI_Epith, 
        reduction = "umap", 
        label = TRUE)

#Go back to RNA Assay
DefaultAssay(Nulliparous_UTI_Epith) <- "RNA"
Nulliparous_UTI_Epith<- ScaleData(Nulliparous_UTI_Epith)

#Determine Percent Ribosomal Genes
Nulliparous_UTI_Epith[["percent.ribo"]] <- PercentageFeatureSet(Nulliparous_UTI_Epith, 
                                                                pattern = c("^Rp[sl]"))

#Visualize Features and Markers
VlnPlot(Nulliparous_UTI_Epith, features = c("nFeature_RNA", "percent.mt", "percent.ribo"))
DotPlot(Nulliparous_UTI_Epith, features = General_Epith_NonEpith_Markers) +RotatedAxis()


#Save Epithelail Re-cluster
saveRDS(Nulliparous_UTI_Epith, file = "filepathway/Nulliparous_UTI_Epith.rds")

#9. Identification of Epithelial Clusters 
#Load in Gene Lists
#Load in Epithelial Signature
MG_Epith_Markers <- read_excel("filepathway/MG_Epith_Markers.xlsx")
MG_Epith_Markers <- MG_Epith_Markers [,1:1]

#Order clusters based on Dendogram
#Make sure in Integrated Assay
DefaultAssay(Nulliparous_UTI_Epith) <- "integrated"
Nulliparous_UTI_Epith<- ScaleData(Nulliparous_UTI_Epith)
#Dendrogram
Nulliparous_UTI_Epith <- BuildClusterTree(Nulliparous_UTI_Epith, reorder = TRUE)
plot(Nulliparous_UTI_Epith@tools$BuildClusterTree, show.tip.label = TRUE)

#Make Sure in RNA assay
DefaultAssay(Nulliparous_UTI_Epith) <- "RNA"
Nulliparous_UTI_Epith <- ScaleData(Nulliparous_UTI_Epith)

#Examine Markers across clusters
DotPlot(MG_Epith_Markers, features = General_Epith_Markers) +RotatedAxis()

#10. Adding Labels to Meta Data
#Label Cells based on Condition
temp_sobj <- Nulliparous_UTI_Epith
temp_sobj <- SetIdent(temp_sobj, 
                      value = "orig.ident")
DimPlot(temp_sobj)
#Create New Names
Condition <- c("UTI", "UTI",
               "Non-UTI", "Non-UTI")
names(Condition) <- levels(temp_sobj)
temp_sobj <- RenameIdents(temp_sobj, Condition)
#Add this data to a new meta data in your sobj
Nulliparous_UTI_Epith@meta.data$Condition <- temp_sobj@active.ident
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith,
                                  value = "Condition")
DimPlot(Nulliparous_UTI_Epith)

#Label Cells According to General Cell Identity 
temp_sobj <- Nulliparous_UTI_Epith
General_Cell_Identity <- c("BM", "LHS", "LASP", 
                           "LHS", "LASP", "LASP", 
                           "LASP", "NA", "LHS", 
                           "NA", "NA")
names(General_Cell_Identity) <- levels(temp_sobj)
temp_sobj <- RenameIdents(temp_sobj, General_Cell_Identity)
#Add this data to a new meta data in your sobj
Nulliparous_UTI_Epith@meta.data$General_Cell_Identity <- temp@active.ident
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith, value = "General_Cell_Identity")
DimPlot(Nulliparous_UTI_Epith)
Nulliparous_UTI_Epith <- subset(Nulliparous_UTI_Epith, idents = c("BM", "LHS", "LASP")) # Used in some analysis

#Combine Conditions and General Cell Identity
#want to combine idents of seurat.clusters and Conditions
test <- Nulliparous_UTI_Epith
Idents(object = test) <- paste(Idents(object = test), 
                               test$Condition, 
                               sep = '_')
Nulliparous_UTI_Epith@meta.data$General_Cell_Identity_Conditions <- test@active.ident
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith, 
                                  value = "General_Cell_Identity_Conditions")
Nulliparous_UTI_Epith$General_Cell_Identity_Conditions <- factor(Nulliparous_UTI_Epith$General_Cell_Identity_Conditions,
                                                                 levels = c("BM_Non-UTI", "BM_UTI",
                                                                            "LHS_Non-UTI", "LHS_UTI",
                                                                            "LASP_Non-UTI", "LASP_UTI"))
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith, 
                                  value = "General_Cell_Identity_Conditions")

#Save Seurat Object with Updated Meta Data
saveRDS(Nulliparous_UTI_Epith, file = "filepathway/Nulliparous_UTI_Epith.rds")