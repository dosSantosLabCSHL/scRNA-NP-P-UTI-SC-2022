######
# Preprocessing of scRNA-seq data from 
# PBS-PLI Mice (Bioproject : PRJNA855880) - 2 Samples
# UTI-PLI Mice (Bioproject : PRJNA855880) - 2 Samples
#Steps Included:
# 1. Load in Data
# 2. Run QC
# 3. Perform Linear Regression
# 4. Determine Parameters (dims and res)
# 5. Define Epithelial Cells 
# 6. Re-cluster Epithelial Cell Populations 
# 7. Identification of Epithelial Clusters 
# 8. Adding Labels to Meta Data
######

#Required Libraries
library(dplyr)
library(Seurat)
library(clustree)
library(ggplot2)

#used code
SL04.data <- Read10X(data.dir = "/Users/slewis/Desktop/smlscrnaseq/SL04/filtered_feature_bc_matrix/")
SL04_Sobj <- CreateSeuratObject(counts = SL04.data, project = "SL04_scRNA", min.cells = 3, min.features = 200)
View(SL04_Sobj@meta.data) #8,368  cells

#1. Load in Libraries 
#PLI PBS and UTI Samples
PLI_sobj <- Read10X(data.dir ="filepathway/PLI/filtered_feature_bc_matrix")
PLI_sobj <-CreateSeuratObject(PLI_sobj_1, 
                                project = "PLI_sobj_1", 
                                min.cells = 3)
PLI_sobj #View the number of cells and features in your sobj before QC
#20,973 cells

# Add column in metadata to allow further subsetting of samples
PLI_sobj@meta.data$temp1<-row.names(PLI_sobj@meta.data)
PLI_sobj@meta.data$Condition<-ifelse(grepl("-5",PLI_sobj@meta.data$temp1),'PLI_PBS',PLI_sobj@meta.data$Condition)
PLI_sobj@meta.data$Condition<-ifelse(grepl("-6",PLI_sobj@meta.data$temp1),'PLI_PBS',PLI_sobj@meta.data$Condition)
PLI_sobj@meta.data$Condition<-ifelse(grepl("-7",PLI_sobj@meta.data$temp1),'PLI_UTI',PLI_sobj@meta.data$Condition)
PLI_sobj@meta.data$Condition<-ifelse(grepl("-8",PLI_sobj@meta.data$temp1),'PLI_UTI',PLI_sobj@meta.data$Condition)
PLI_sobj@meta.data$temp1<-NULL # Remove temporary column

#2. Run QC
#Define %MT content
PLI_sobj[["percent.mt"]] <- PercentageFeatureSet(PLI_sobj,pattern = c("^mt-"))

#Conduct QC 
#Visualize the number of genes, the amount of RNA and the percentage of mitochodrial genes in each cell
VlnPlot(PLI_sobj, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"))
#normalize data
PLI_sobj <- NormalizeData(PLI_sobj, normalization.method = "LogNormalize", scale.factor = 10000)
PLI_sobj <- FindVariableFeatures(PLI_sobj, selection.method = "vst", nfeatures = 2000)

#Determine the metrics of these features that would be considered outliers and adjust the line below accordingly
PLI_sobj <- subset(PLI_sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
#Examine these features in Vlnplots again to examine appropriate cells were removed
VlnPlot(PLI_sobj, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"))
#View sobj to examine how the number of cells and features changed after following QC protocol
View(PLI_sobj@meta.data)

#3.Perform Linear Regression
#Run PCA
PLI_sobj <- RunPCA(PLI_sobj,features = VariableFeatures(object = PLI_sobj))
#Visualize PCA reduction
#DimPlot(PLI_sobj, reduction = "pca")

# 4. Determine Parameters (dims and res)
ElbowPlot(PLI_sobj) #Choose Value were the value on the y axis remains about constant

#Determine Resolution of the dataset
PLI_sobj <- FindNeighbors(PLI_sobj, dims = 1:10) #change depending on the dimensionality used
PLI_sobj <- FindClusters(PLI_sobj, resolution = 0.2)
PLI_sobj <- FindClusters(PLI_sobj, resolution = 0.4)
PLI_sobj <- FindClusters(PLI_sobj, resolution = 0.5)
PLI_sobj <- FindClusters(PLI_sobj, resolution = 0.8)
PLI_sobj <- FindClusters(PLI_sobj, resolution = 1)
clustree(PLI_sobj, prefix = "integrated_snn_res.")

#Cluster Cells: Adjust Parameters in the following functions based on elbowplot and clustree
PLI_sobj <- FindNeighbors(PLI_sobj, dims = 1:12)
PLI_sobj <- FindClusters(PLI_sobj,resolution = 0.5)

#Visualize Clusters Using UMAP
PLI_sobj <- RunUMAP(PLI_sobj, dims = 1:12) #changes dims as needed
DimPlot(PLI_sobj, reduction = "umap")
#initial umap has 18 clusters

#Examine these features in Vlnplots again to examine appropriate cells were removed
VlnPlot(PLI_sobj, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"))
#remove cluster 16 due to high mito RNA and/or low features and counts

#Save Seurat Object
saveRDS(PLI_sobj, file = "filepathway/PLI_sobj.rds")

#4. Define Epithelial Cells 
#Read in File if Neccessary 
PLI_sobj <- readRDS("filepathway/PLI_sobj.rds")

#Load in Gene Markers : Derived from DOI: https://doi.org/10.1007/s10911-021-09486-3
SC_Identification_Markers <- read_excel("filepathway/SC_Identification_Markers.xlsx")

#Visualize Markers Across Clusters
DotPlot(PLI_sobj, 
        features = SC_Identification_Markers) + RotatedAxis()

#Examine Specific Genes in FP
FeaturePlot(PLI_sobj, 
            features = "Epcam")
FeaturePlot(PLI_sobj, 
            features = "Cd3e")
FeaturePlot(PLI_sobj, 
            features = "Ms4a1")
FeaturePlot(PLI_sobj, 
            features = "Cd19")
FeaturePlot(PLI_sobj, 
            features = "Tyrobp")
FeaturePlot(PLI_sobj, 
            features = "Sparc")
FeaturePlot(PLI_sobj, 
            features = "Fabp4")
FeaturePlot(PLI_sobj, 
            features = "Vcam1")

#Subset Epithelial Cells
#0, 3, 4, 6, 11, 15, 17
PLI_sobj_Epith <- subset(PLI_sobj_Epith, 
                                idents = c("0", "3", "4", "6", "11", "15", "17"))
#Examine proper clusters were subsetted
DimPlot(PLI_sobj_Epith, 
        reduction = "umap", 
        label = TRUE)

#5. Recluster Epithelial Cell Populations 
ScaleData(PLI_sobj_Epith)

#Find Highly Variable Features
PLI_sobj_Epith <- FindVariableFeatures(PLI_sobj_Epith, 
                                              selection.method = "vst", 
                                              nfeatures = 2000)


#Cluster Cells
PLI_sobj_Epith <- RunPCA(PLI_sobj_Epith, 
                                features = VariableFeatures(object = PLI_sobj_Epith))
#Determine dimensionality
ElbowPlot(PLI_sobj_Epith)
#Determine Resolution
PLI_sobj_Epith <- FindNeighbors(PLI_sobj_Epith, 
                                       dims = 1:10)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.2)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.4)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.5)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.8)
clustree(PLI_sobj_Epith, prefix = "integrated_snn_res.")
#Visualize Clusters in UMAP
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.5)
PLI_sobj_Epith <- RunUMAP(PLI_sobj_Epith, 
                                 reduction = 'pca', 
                                 dims = 1:10) 
DimPlot(PLI_sobj_Epith, 
        reduction = "umap", 
        label = TRUE)
# have 12 clusters

#Save First Epithelail Re-cluster
saveRDS(PLI_sobj_Epith, file = "filepathway/PLI_sobj_Epith.rds") 
#Examine Quality of Clusters 

#Determine Percent Ribosomal Genes
PLI_sobj_Epith[["percent.ribo"]] <- PercentageFeatureSet(PLI_sobj_Epith, 
                                                                pattern = c("^Rp[sl]"))

#Visualize Features and Markers
VlnPlot(PLI_sobj_Epith, features = c("nFeature_RNA", "percent.mt", "percent.ribo"))
DotPlot(PLI_sobj_Epith, features = SC_Identification_Markers) +RotatedAxis()

#6. Second Epithelial Re-clustering 
#subset high quality clustes (removed those of low quality or of non-epithelail origin - 6,10)
PLI_sobj_Epith <- subset(PLI_sobj_Epith, 
                                idents = c("0", "1", "2", "3", "4", "5", "7", "8", "9"))

DimPlot(PLI_sobj_Epith, 
        reduction = "umap", 
        label = TRUE)

#Find Highly Variable Features
PLI_sobj_Epith <- FindVariableFeatures(PLI_sobj_Epith, 
                                              selection.method = "vst", 
                                              nfeatures = 2000)

#Cluster Cells
PLI_sobj_Epith <- RunPCA(PLI_sobj_Epith, 
                                features = VariableFeatures(object = PLI_sobj_Epith))
#Determine dimensionality
ElbowPlot(PLI_sobj_Epith)
#Determine Resolution
PLI_sobj_Epith <- FindNeighbors(PLI_sobj_Epith, 
                                       dims = 1:10)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.2)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.4)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.5)
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.8)
clustree(PLI_sobj_Epith, prefix = "integrated_snn_res.")
#Visualize Clusters in UMAP
PLI_sobj_Epith <- FindClusters(PLI_sobj_Epith, 
                                      resolution = 0.5)
PLI_sobj_Epith <- RunUMAP(PLI_sobj_Epith, 
                                 reduction = 'pca', 
                                 dims = 1:10) 
DimPlot(PLI_sobj_Epith, 
        reduction = "umap", 
        label = TRUE)

#Determine Percent Ribosomal Genes
PLI_sobj_Epith[["percent.ribo"]] <- PercentageFeatureSet(PLI_sobj_Epith, 
                                                                pattern = c("^Rp[sl]"))

#Visualize Features and Markers
VlnPlot(PLI_sobj_Epith, features = c("nFeature_RNA", "percent.mt", "percent.ribo"))
DotPlot(PLI_sobj_Epith, features = General_Epith_NonEpith_Markers) +RotatedAxis()


#Save Epithelail Re-cluster
saveRDS(PLI_sobj_Epith, file = "filepathway/PLI_sobj_Epith.rds")

#7. Identification of Epithelial Clusters 
#Load in Gene Lists
#Load in Epithelial Signature
MG_Epith_Markers <- read_excel("filepathway/MG_Epith_Markers.xlsx")
MG_Epith_Markers <- MG_Epith_Markers [,1:1]

#Order clusters based on Dendogram

#Dendrogram
PLI_sobj_Epith <- BuildClusterTree(PLI_sobj_Epith, reorder = TRUE)
plot(PLI_sobj_Epith@tools$BuildClusterTree, show.tip.label = TRUE)

#Examine Markers across clusters
DotPlot(MG_Epith_Markers, features = General_Epith_Markers) +RotatedAxis()

#8. Adding Labels to Meta Data
#Label Cells based on Condition
temp_sobj <- PLI_sobj_Epith
temp_sobj <- SetIdent(temp_sobj, 
                      value = "orig.ident")
DimPlot(temp_sobj)


#Label Cells According to General Cell Identity 
temp_sobj <- PLI_sobj_Epith
General_Cell_Identity <- c("BM1", "BM2",
                           "BL1", "BL2",
                           "LHS1",
                           "LASP1", "LASP2", "LASP3",
                           "LASP4", "LASP5")
names(General_Cell_Identity) <- levels(temp_sobj)
temp_sobj <- RenameIdents(temp_sobj, General_Cell_Identity)
#Add this data to a new meta data in your sobj
PLI_sobj_Epith@meta.data$General_Cell_Identity <- temp@active.ident
PLI_sobj_Epith <- SetIdent(PLI_sobj_Epith, value = "General_Cell_Identity")
DimPlot(PLI_sobj_Epith)
PLI_sobj_Epith <- subset(PLI_sobj_Epith, idents = c("BM", "LHS", "LASP")) # Used in some analysis

#Save Seurat Object with Updated Meta Data
saveRDS(PLI_sobj_Epith, file = "filepathway/PLI_sobj.rds")
