#####
#Analysis for figures done in Nulliparous No-UTI and UTI samples in just Epithelial Clusters
#Present in Manuscript:
#Host response during unresolved urinary tract infection alters female mammary tissue homeostasis through collagen deposition and TIMP1 
#Generation of Figures: 2A-C, 4A
#Generation of Supplementary Figures: 3A-B
#Generation of Supplementary Data: 1-3
#####
#Load in Libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(UCell) #Generate Signature Score
library(ggpubr) #Generate Stats for Desired Comparisions 
library(readxl) #Read in excel files
library(writexl) #Write out excel files

#Figure 2A: Split UMAP 
DimPlot(Nulliparous_UTI_Epith, 
        split.by = "Condition")

#Figure 2B: Lineage Identification, Visualization DotPlot
DotPlot(Nulliparous_UTI_Epith, 
        features = c("Epcam","Krt5","Krt14",
                     "Trp63","Pdpn","Krt8",
                     "Krt18","Prom1","Esr1",
                     "Prlr","Pgr","Ly6a",
                     "Cd14","Kit","Csn1s1",
                     "Csn2","Csn3",
                     "Mki67","Ube2c","Top2a"), 
        cols = "RdBu")+ 
  RotatedAxis() 

#Figure 2C: Vlnplot of YAP Signaling Signature
#YAP Signaling Signature Derived from: doi: 10.1016/j.cell.2011.09.048.
YAP_Signaling_Signature <- read_excel("/filepathway/YAP_Signaling_Signature.xlsx")
#Add Score to Metadata
Nulliparous_UTI_Epith <- AddModuleScore_UCell(Nulliparous_UTI_Epith, 
                                        features = YAP_Signaling_Signature)
#Visualize
VlnPlot(Nulliparous_UTI_Epith, 
        features = "YAP_Signaling_Signature_UCell", 
        split.by = "Condition",
        cols = c("darkblue", "darkred"))
#Generate Stats
#Change ident and quantify Stats
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith, value = "General_Cell_Identity_Conditions")
#Define Comparisons
my_comparisons <- list(c("BM_Non-UTI", "BM_UTI"), 
                       c("LHS_Non-UTI", "LHS_UTI"),
                       c("LASP_Non-UTI", "LASP_UTI"))
VlnPlot(nullwildtype_Epith, 
        features = "YAP_Signaling_Signature_UCell") + 
  geom_boxplot(position=position_dodge(1), color="black") +
  stat_compare_means(comparisons = my_comparisons) +
  ylim(0, 1)

#Supplmentary Data 1: YAP signaling signature in NP Epithelial Clusters (related to Fig. 2C)
YAP_Signaling_Signature <- YAP_Signaling_Signature$YAP_Signaling_Signature
YAP_Signaling_Signature.markers <- FindAllMarkers(Nulliparous_UTI_Epith, 
                                                  features = intersect(rownames(Nulliparous_UTI_Epith),
                                                                       YAP_Signaling_Signature),
                                                  logfc.thresh
                                                  min.cells.feature = 0,
                                                  min.cells.group = 0,
                                                  return.thresh = 1)
write_xlsx(YAP_Signaling_Signature.markers, "filepathway/Null_YAP_Signaling_Signature_Expression.xlsx")

#Supplementary Figure 2A: Vlnplot of MEC-ECM Communication Signature
#MEC_ECM_CommunicationSignaling Signature Derived from: KEGG
MEC_ECM_Communication_Signature <- read_excel("/filepathway/MEC_ECM_Communication_Signature.xlsx")
#Add Score to Metadata
Nulliparous_UTI_Epith <- AddModuleScore_UCell(Nulliparous_UTI_Epith, 
                                              features = MEC_ECM_Communication_Signature)
#Visualize
VlnPlot(Nulliparous_UTI_Epith, 
        features = "MEC_ECM_Communication_Signature_UCell", 
        split.by = "Condition",
        cols = c("darkblue", "darkred"))
#Generate Stats
#Change ident and quantify Stats
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith, value = "General_Cell_Identity_Conditions")
#Define Comparisons
my_comparisons <- list(c("BM_Non-UTI", "BM_UTI"), 
                       c("LHS_Non-UTI", "LHS_UTI"),
                       c("LASP_Non-UTI", "LASP_UTI"))
VlnPlot(nullwildtype_Epith, 
        features = "MEC_ECM_Communication_Signature_UCell") + 
  geom_boxplot(position=position_dodge(1), color="black") +
  stat_compare_means(comparisons = my_comparisons) +
  ylim(0, 1)

#Supplmentary Data 2: MEC-ECM Communication Signature in NP Epithelial Clusters (related to Fig. S3A)
MEC_ECM_Communication_Signature <- MEC_ECM_Communication_Signature$MEC_ECM_Communication_Signature
MEC_ECM_Communication_Signature.markers <- FindAllMarkers(Nulliparous_UTI_Epith, 
                                                  features = intersect(rownames(Nulliparous_UTI_Epith),
                                                                       MEC_ECM_Communication_Signature),
                                                  logfc.thresh
                                                  min.cells.feature = 0,
                                                  min.cells.group = 0,
                                                  return.thresh = 1)
write_xlsx(MEC_ECM_Communication_Signature.markers, "filepathway/Null_MEC_ECM_Communication_Signature_Expression.xlsx")

#Supplementary Figure 3B: Vlnplot of Mechano-sensing Signaling Signature
#YAP Signaling Signature Derived from: doi: 10.1016/j.cell.2011.09.048.
Mechano-sensing_Signature <- read_excel("/filepathway/Mechano-sensing_Signature.xlsx")
#Add Score to Metadata
Nulliparous_UTI_Epith <- AddModuleScore_UCell(Nulliparous_UTI_Epith, 
                                              features = Mechano-sensing_Signature)
#Visualize
VlnPlot(Nulliparous_UTI_Epith, 
        features = "Mechano-sensing_Signature_UCell", 
        split.by = "Condition",
        cols = c("darkblue", "darkred"))
#Generate Stats
#Change ident and quantify Stats
Nulliparous_UTI_Epith <- SetIdent(Nulliparous_UTI_Epith, value = "General_Cell_Identity_Conditions")
#Define Comparisons
my_comparisons <- list(c("BM_Non-UTI", "BM_UTI"), 
                       c("LHS_Non-UTI", "LHS_UTI"),
                       c("LASP_Non-UTI", "LASP_UTI"))
VlnPlot(nullwildtype_Epith, 
        features = "Mechano-sensing_Signature_UCell") + 
  geom_boxplot(position=position_dodge(1), color="black") +
  stat_compare_means(comparisons = my_comparisons) +
  ylim(0, 1)

#Supplmentary Data 3: mechano-signaling signature in NP Epithelial Clusters (related to Fig. S3B)
Mechano-sensing_Signature <- Mechano-sensing_Signature$Mechano-sensing_Signature
Mechano-sensing_Signature.markers <- FindAllMarkers(Nulliparous_UTI_Epith, 
                                                  features = intersect(rownames(Nulliparous_UTI_Epith),
                                                                       Mechano-sensing_Signature),
                                                  logfc.thresh
                                                  min.cells.feature = 0,
                                                  min.cells.group = 0,
                                                  return.thresh = 1)
write_xlsx(Mechano-sensing_Signature.markers, "filepathway/Null_Mechano-sensing_Signature_Expression.xlsx")

