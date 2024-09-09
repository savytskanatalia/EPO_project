# Summary of results
CRWvsPlacebo_RiboTag <- read.delim("/data/natalia/EPO/results/CRWvsPlacebo_RiboTag_genes_DEresults_filt.txt")
EPOvsPlacebo_RiboTag <- read.delim("/data/natalia/EPO/results/EPOvsPlacebo_RiboTag_genes_DEresults_filt.txt")
HypoxiavsPlacebo_RiboTag <- read.delim("/data/natalia/EPO/results/HypoxiavsPlacebo_RiboTag_genes_DEresults_filt.txt")

CRWvsPlacebo_RiboTag$Group<-"CRWvsPlacebo"
EPOvsPlacebo_RiboTag$Group<-"EPOvsPlacebo"
HypoxiavsPlacebo_RiboTag$Group<-"HypoxiavsPlacebo"

RiboTag_analysis<-rbind(CRWvsPlacebo_RiboTag,EPOvsPlacebo_RiboTag,HypoxiavsPlacebo_RiboTag)
RiboTag_analysis$FDR_joint<-p.adjust(RiboTag_analysis$pvalue,method="fdr")

# stringent
RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.05 & abs(RiboTag_analysis$log2FoldChange)>=1,]
# 6 hits for CRW+Hypoxia
RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.1 & abs(RiboTag_analysis$log2FoldChange)>=0.5,]
# only 1 hit for EPO, ~100 hits for CRW and Hypoxia

# DE in CRW
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & abs(RiboTag_analysis$log2FoldChange)>=0.5) & RiboTag_analysis$Group=="CRWvsPlacebo",])
# [1] 73
# UP in CRW
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & RiboTag_analysis$log2FoldChange>=0.5) & RiboTag_analysis$Group=="CRWvsPlacebo",])
# [1] 28
# DE in Hypoxia
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & abs(RiboTag_analysis$log2FoldChange)>=0.5) & RiboTag_analysis$Group=="HypoxiavsPlacebo",])
# [1] 41
# UP in hypoxia
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & RiboTag_analysis$log2FoldChange>=0.5) & RiboTag_analysis$Group=="HypoxiavsPlacebo",])
# [1] 13

# checking overlap for UP
up_hyp<-RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & RiboTag_analysis$log2FoldChange>=0.5) & RiboTag_analysis$Group=="HypoxiavsPlacebo","X"]
up_crw<-RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & RiboTag_analysis$log2FoldChange>=0.5) & RiboTag_analysis$Group=="CRWvsPlacebo","X"]
up_hyp[up_hyp %in% up_crw]
# [1] "ENSMUSG00000003545.4"  "ENSMUSG00000037868.16" "ENSMUSG00000085609.3"  "ENSMUSG00000021250.14" "ENSMUSG00000045903.10"
# checking overlap for DOWN
down_hyp<-RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & RiboTag_analysis$log2FoldChange<(-0.5)) & RiboTag_analysis$Group=="HypoxiavsPlacebo","X"]
down_crw<-RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.1 & RiboTag_analysis$log2FoldChange<(-0.5)) & RiboTag_analysis$Group=="CRWvsPlacebo","X"]
down_hyp[down_hyp %in% down_crw]
# [1] "ENSMUSG00000027375.15" "ENSMUSG00000044576.7"  "ENSMUSG00000039683.17" "ENSMUSG00000032667.18" "ENSMUSG00000036634.16"
# [6] "ENSMUSG00000069919.8"  "ENSMUSG00000020932.15" "ENSMUSG00000054889.11" "ENSMUSG00000008393.10" "ENSMUSG00000046160.7" 
# [11] "ENSMUSG00000031425.16"
# get the IDs for mildly DE; plot a heatmap for their expression between the groups - with clustering by both rows and columns?
mildly_de<-unique(RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.1 & abs(RiboTag_analysis$log2FoldChange)>=0.5,"X"])

# get normalized expression for these genes; t, z0score, plot
# dds2_nc 
# scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))
# full_meta_epo[full_meta_epo$i.Tag=="RiboTag_IAA",]

dds2_nc<-as.data.frame(dds2_nc)
dds2_nc<-dds2_nc[mildly_de,]

# heatmap
library(ComplexHeatmap)
library(pheatmap)
library(circlize)

# Scale rows (Z-score)
dt_zs<-dds2_nc %>% pheatmap:::scale_rows()


# need as matrix

# For -+ gradient
col_fun = colorRamp2(c(-2.5, 0, 2.5), c("dodgerblue4", "white", "red"))


htl<-Heatmap(as.matrix(dt_zs),cluster_columns  = T, col=col_fun,
             row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 10),
             cluster_rows  = TRUE,
             heatmap_legend_param = list(
               title = "Z-score", at = c(-2.5,-1,0,1,2.5), 
               labels = c(-2.5,-1,0,1,2.5),direction = "horizontal"))

draw(htl,heatmap_legend_side = "bottom")

# col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")
c("#154360", "#FF5733", "#FFC300","grey") 
# Example with column annotation
ha_col = HeatmapAnnotation(Group = full_meta_epo[full_meta_epo$i.Tag=="RiboTag_IAA","group"],
                           col = list(Group = c("CRW" = "#154360", "EPO" = "#FF5733", "Hypoxia" = "#FFC300",
                                                "Placebo" = "grey")),
                                      show_annotation_name = c(Group = FALSE),
                           annotation_legend_param = list(
                             Group = list(direction = "horizontal",nrow = 1)))


htl<-Heatmap(as.matrix(dt_zs),cluster_columns  = T, col=col_fun,
             row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 15),
             cluster_rows  = TRUE,
             heatmap_legend_param = list(
               title = "Z-score", at = c(-2.5,-1,0,1,2.5),
               labels = c(-2.5,-1,0,1,2.5),direction = "horizontal"), top_annotation = ha_col)

draw(htl,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


# examining results from total
CRWvsPlacebo_TotalRNA <- read.delim("/data/natalia/EPO/results/CRWvsPlacebo_TotalRNA_genes_DEresults_filt.txt")
EPOvsPlacebo_TotalRNA <- read.delim("/data/natalia/EPO/results/EPOvsPlacebo_TotalRNA_genes_DEresults_filt.txt")
HypoxiavsPlacebo_TotalRNA <- read.delim("/data/natalia/EPO/results/HypoxiavsPlacebo_TotalRNA_genes_DEresults_filt.txt")

CRWvsPlacebo_TotalRNA$Group<-"CRWvsPlacebo"
EPOvsPlacebo_TotalRNA$Group<-"EPOvsPlacebo"
HypoxiavsPlacebo_TotalRNA$Group<-"HypoxiavsPlacebo"
Total_analysis<-rbind(CRWvsPlacebo_TotalRNA,EPOvsPlacebo_TotalRNA,HypoxiavsPlacebo_TotalRNA)

Total_analysis$FDR_joint<-p.adjust(Total_analysis$pvalue,method="fdr")
# get the list of DE (both mildly and strong); for whole list - make a heatmap; for strong list - check the genes
# check overlaps

totalrna_de<-Total_analysis[abs(Total_analysis$log2FoldChange)>=0.5 & Total_analysis$FDR_joint<=0.05,"X"]
dds3_nc<-as.data.frame(dds3_nc)
dds3_nc<-dds3_nc[totalrna_de,]


dt_zs2<-dds3_nc %>% pheatmap:::scale_rows()

# remove accidental duplicates of rows
# "ENSMUSG00000025270.14.1", "ENSMUSG00000073940.4.2", "ENSMUSG00000073940.4.1", "ENSMUSG00000052305.7.2", "ENSMUSG00000052305.7.1","ENSMUSG00000069919.8.2","ENSMUSG00000069919.8.1", "ENSMUSG00000069917.8.2", "ENSMUSG00000069917.8.1", "ENSMUSG00000037868.16.1", "ENSMUSG00000045903.10.1", "ENSMUSG00000003545.4.1", "ENSMUSG00000041324.15.1", "ENSMUSG00000021250.14.1"
dt_zs2<-dt_zs2[!(rownames(dt_zs2) %in% c("ENSMUSG00000025270.14.1", "ENSMUSG00000073940.4.2", "ENSMUSG00000073940.4.1", "ENSMUSG00000052305.7.2", "ENSMUSG00000052305.7.1","ENSMUSG00000069919.8.2","ENSMUSG00000069919.8.1", "ENSMUSG00000069917.8.2", "ENSMUSG00000069917.8.1", "ENSMUSG00000037868.16.1", "ENSMUSG00000045903.10.1", "ENSMUSG00000003545.4.1", "ENSMUSG00000041324.15.1", "ENSMUSG00000021250.14.1")),]

# need as matrix

# For -+ gradient
col_fun2 = colorRamp2(c(-2.5, 0, 4), c("dodgerblue4", "white", "red"))



# col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")
c("#154360", "#FF5733", "#FFC300","grey") 
# Example with column annotation
ha_col2 = HeatmapAnnotation(Group = full_meta_epo[full_meta_epo$i.Tag!="RiboTag_IAA","group"],
                           col = list(Group = c("CRW" = "#154360", "EPO" = "#FF5733", "Hypoxia" = "#FFC300",
                                                "Placebo" = "grey")),
                           show_annotation_name = c(Group = FALSE),
                           annotation_legend_param = list(
                             Group = list(direction = "horizontal",nrow = 1)))


htl2<-Heatmap(as.matrix(dt_zs2),cluster_columns  = T, col=col_fun2,
             row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 15),
             cluster_rows  = TRUE,
             heatmap_legend_param = list(
               title = "Z-score", at = c(-2.5,-1,0,1,2.5),
               labels = c(-2.5,-1,0,1,2.5),direction = "horizontal"), top_annotation = ha_col2)

draw(htl2,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
# 23067a054L2_01

total_ribo_overlap<-totalrna_de$X[totalrna_de$X %in% c(up_hyp,down_hyp,up_crw,down_crw)]


# To generate:
# - volcano plots per group comparison for Ribo FDR<0.05
# - Heatmap for prev
# Heatmap for cellular markers across all samples


# cell markers 

# Mouse cell markers
# Microglia: Iba1, Cx3cr1, Csf1r (ENSMUSG00000024397, ENSMUSG00000052336, ENSMUSG00000024621)
# Neuron: NeuN, Map2, Dcx, Gad1, Stmn2 (ENSMUSG00000025576 , ENSMUSG00000025576 , ENSMUSG00000027500 )
# Oligodendrocytes: Olig2, Sox10, Ptgds (ENSMUSG00000039830,ENSMUSG00000033006 )
# Endothelial cell: Cldn5, Flt1, CD31 ENSMUSG00000041378 

# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000374
# PV: Slc1a6, Atp2a3, Ppp1r17, Car8, Homer3, Pcp2
# GABA PV+: Slc32a1,Slc6a1, Gad1, Pvalb, Gad2, Lhx6
# GABA PV-: Gpr6, Adora2a,Vip, Sst,Penk, Htr3a
# Glut: Slc17a7,Neurod6,Dkk3,Slc17a6,Stx1a,Satb2
# Astro: Apoe,Aqp4,Fabq7,Slc6a11,S100b,Slc1a3

# add prev for other cell types
# load table with markers
mouse_cell_type_markers <- read.delim("/data/natalia/EPO/mouse_cell_type_markers.txt")
mouse_cell_type_markers$FullName<-""
for (i in mouse_cell_type_markers$Ensembl){
  mouse_cell_type_markers[mouse_cell_type_markers$Ensembl==i,"FullName"]<-rownames(dds_nc)[rownames(dds_nc) %like% i]
}

# subset normalized expression 
dds_nc_cmarkers<-dds_nc[mouse_cell_type_markers$FullName,]
identical(rownames(dds_nc_cmarkers),mouse_cell_type_markers$FullName)
# [1] TRUE
rownames(dds_nc_cmarkers)<-mouse_cell_type_markers$Gene_Symbol
#  construct heatmap object



dds_nc_cmarkers_zs<-dds_nc_cmarkers %>% pheatmap:::scale_rows()
# For -+ gradient
col_fun_m = colorRamp2(c(-2, 0, 4), c("dodgerblue4", "white", "red"))



# col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")
c("#154360", "#FF5733", "#FFC300","grey") 
# Example with column annotation
ha_col_m = HeatmapAnnotation(Group = full_meta_epo$group, Experiment = full_meta_epo$i.Tag,
                            col = list(Group = c("CRW" = "#154360", "EPO" = "#FF5733", "Hypoxia" = "#FFC300",
                                                 "Placebo" = "grey"),
                                       Experiment = c("TotalRNA" = "#822681FF", "RiboTag_IAA" = "#B63679FF")),
                            show_annotation_name = c(Group = FALSE, Experiment = FALSE),
                            annotation_legend_param = list(
                              Group = list(direction = "horizontal",nrow = 1), Experiment = list(direction = "horizontal",nrow = 1)))

library(scales)
library(viridis)
# show_col(viridis_pal(option = "viridis")(10))
c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")



ha_row_m = rowAnnotation(CellType = mouse_cell_type_markers$CellType, BroadType = mouse_cell_type_markers$Group,
                              col = list(CellType = c("Astroglia" = "#1F9E89FF", "Microglia" = "#35B779FF", "Oligodendrocytes" = "#6DCD59FF", 
                                                      "Neuron (general)" = "#440154FF", "PV (neu)" = "#482878FF", "GABA PV+" = "#3E4A89FF",  "GABA PV-" = "#31688EFF", "Glutamatergic" = "#26828EFF",  
                                                       "Endothelial" =  "#B4DE2CFF", "Erythroid" = "#FDE725FF"),
                                         BroadType = c("Glia" ="#35B779FF" ,"Neuron" = "#31688EFF" ,"Other" = "#FDE725FF")),
                              show_annotation_name = c(CellType = FALSE, BroadType = FALSE),
                              annotation_legend_param = list(
                                CellType = list(direction = "horizontal",nrow = 2), BroadType = list(direction = "horizontal",nrow = 1)))



htl_m<-Heatmap(as.matrix(dds_nc_cmarkers_zs),cluster_columns  = T, col=col_fun_m,
              row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 15),
              cluster_rows  = F,
              heatmap_legend_param = list(
                title = "Z-score", at = c(-2.5,-1,0,1,2.5),
                labels = c(-2.5,-1,0,1,2.5),direction = "horizontal"), top_annotation = ha_col_m, right_annotation = ha_row_m)

draw(htl_m,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

# get numbers for ribo w/o log2FC cutoff
nrow(RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.05 & RiboTag_analysis$Group=="CRWvsPlacebo",])
# [1] 381
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.05 & RiboTag_analysis$Group=="CRWvsPlacebo") & RiboTag_analysis$log2FoldChange>0,])
# [1] 224
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.05 & RiboTag_analysis$Group=="HypoxiavsPlacebo") & RiboTag_analysis$log2FoldChange>0,])
# [1] 68
nrow(RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.05 & RiboTag_analysis$Group=="HypoxiavsPlacebo",])
# [1] 203
nrow(RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.05 & RiboTag_analysis$Group=="EPOvsPlacebo",])
# [1] 109
nrow(RiboTag_analysis[(RiboTag_analysis$FDR_joint<=0.05 & RiboTag_analysis$Group=="EPOvsPlacebo") & RiboTag_analysis$log2FoldChange>0,])
# [1] 44

# Volcano plots per comparison
library(EnhancedVolcano)
# Example of data frame
#                            baseMean   log2FoldChange    lfcSE       stat        pvalue          padj
# Gene1   22.301338      -2.073119 0.08390326 -24.708447 8.677517e-135 2.394500e-133

EnhancedVolcano(RiboTag_analysis[RiboTag_analysis$Group=="CRWvsPlacebo",],
                lab = RiboTag_analysis[RiboTag_analysis$Group=="CRWvsPlacebo","X"],
                x = 'log2FoldChange', #colname to be plotted on X
                y = 'FDR_joint',           # colname to be plotted on Y
                title = 'CRW vs Placebo',
                pCutoff = 0.05,       # Cutoff for Y
                FCcutoff = 0.1, ylim = 0.0000005)       # Cutoff for X


EnhancedVolcano(RiboTag_analysis[RiboTag_analysis$Group=="HypoxiavsPlacebo",],
                lab = RiboTag_analysis[RiboTag_analysis$Group=="HypoxiavsPlacebo","X"],
                x = 'log2FoldChange', #colname to be plotted on X
                y = 'FDR_joint',           # colname to be plotted on Y
                title = 'Hypoxia vs Placebo',
                pCutoff = 0.05,       # Cutoff for Y
                FCcutoff = 0.1, ylim = 0.0000005) 

EnhancedVolcano(RiboTag_analysis[RiboTag_analysis$Group=="EPOvsPlacebo",],
                lab = RiboTag_analysis[RiboTag_analysis$Group=="EPOvsPlacebo","X"],
                x = 'log2FoldChange', #colname to be plotted on X
                y = 'FDR_joint',           # colname to be plotted on Y
                title = 'EPO vs Placebo',
                pCutoff = 0.05,       # Cutoff for Y
                FCcutoff = 0.1, ylim = 0.0000005, xlim=c(-1, 1)) 


# Heatmap
ribo_softfilter<-unique(RiboTag_analysis[RiboTag_analysis$FDR_joint<=0.05,"X"])
# get normalized cnts and subset to genes of interest
ribotag_genes_MedRat_norm <- read.delim("/data/natalia/EPO/results/ribotag_genes_MedRat_norm.txt", row.names=1)
colnames(ribotag_genes_MedRat_norm)<-str_remove(colnames(ribotag_genes_MedRat_norm),"X")
ribotag_genes_MedRat_norm<-ribotag_genes_MedRat_norm[ribo_softfilter,]
# plot



dds_zs<-ribotag_genes_MedRat_norm %>% pheatmap:::scale_rows()
# For -+ gradient
col_fun_s = colorRamp2(c(-3.5, 0, 3.5), c("dodgerblue4", "white", "red"))


# Example with column annotation
ha_col_s = HeatmapAnnotation(Group = full_meta_epo[full_meta_epo$i.Tag=="RiboTag_IAA","group"],
                             col = list(Group = c("CRW" = "#154360", "EPO" = "#FF5733", "Hypoxia" = "#FFC300",
                                                  "Placebo" = "grey")),
                             show_annotation_name = c(Group = FALSE),
                             annotation_legend_param = list(
                               Group = list(direction = "horizontal",nrow = 1)))



htl_s<-Heatmap(as.matrix(dds_zs),cluster_columns  = T, col=col_fun_s,
               row_names_gp = gpar(fontsize = 0),column_names_gp = gpar(fontsize = 15),
               cluster_rows  = T,
               heatmap_legend_param = list(
                 title = "Z-score", at = c(-3.5,-1,0,1,3.5),
                 labels = c(-3.5,-1,0,1,3.5),direction = "horizontal"), top_annotation = ha_col_s)

draw(htl_s,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

# get GO annotations per DE set of genes
# get all expressed genes to contrast against
all_genes<-unique(RiboTag_analysis[RiboTag_analysis$baseMean>=3,"X"])

# Edit all genes and tested genes to drop suffix after .
# gsub("\\..*","",ipd_dif_exp_isos_down)
all_genes<-gsub("\\..*","",all_genes)

Hypo_genes<-gsub("\\..*","",RiboTag_analysis[RiboTag_analysis$Group=="HypoxiavsPlacebo" & RiboTag_analysis$FDR_joint<=0.05,"X"])
CRW_genes<-gsub("\\..*","",RiboTag_analysis[RiboTag_analysis$Group=="CRWvsPlacebo" & RiboTag_analysis$FDR_joint<=0.05,"X"])
EPO_genes<-gsub("\\..*","",RiboTag_analysis[RiboTag_analysis$Group=="EPOvsPlacebo" & RiboTag_analysis$FDR_joint<=0.05,"X"])


Hypo_tested <- factor(as.integer(all_genes %in% Hypo_genes))
names(Hypo_tested) <- all_genes

CRW_tested <- factor(as.integer(all_genes %in% CRW_genes))
names(CRW_tested) <- all_genes

EPO_tested <- factor(as.integer(all_genes %in% EPO_genes))
names(EPO_tested) <- all_genes

# GO

Hypo_allRes_tested<-myGO(Hypo_tested)
fwrite(Hypo_allRes_tested,file="/data/natalia/EPO/results/Ribo_HypovsPlacebo_above3baseMean_softfilt_GO.txt",sep="\t",row.names=TRUE)
CRW_allRes_tested<-myGO(CRW_tested)
fwrite(CRW_allRes_tested,file="/data/natalia/EPO/results/Ribo_CRWvsPlacebo_above3baseMean_softfilt_GO.txt",sep="\t",row.names=TRUE)
EPO_allRes_tested<-myGO(EPO_tested)
fwrite(EPO_allRes_tested,file="/data/natalia/EPO/results/Ribo_EPOvsPlacebo_above3baseMean_softfilt_GO.txt",sep="\t",row.names=TRUE)

ggplot(CRW_allRes_tested, aes(x=Term, y=Significant, fill=Type)) +
  geom_bar(stat="identity")+coord_flip()+ 
  scale_fill_viridis(discrete=TRUE, direction=1)+
  labs(fill="", y="Number of genes", x="Significantly Enriched Terms")+theme_light()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# remove factoring, do it out of function manually as we have identical shortening of related but different terms

myGO<-function(antisense){
  # GO analysis for Biological Process
  GOdata_as <- new("topGOdata",
                   ontology = "BP",
                   allGenes = antisense,
                   nodeSize = 5,
                   annot = annFUN.org,
                   mapping = "org.Mm.eg.db",
                   ID = "ensembl")
  
  resultFisher_as <- runTest(GOdata_as, algorithm = "classic", statistic = "fisher")
  allRes_as <- GenTable(GOdata_as, classicFisher = resultFisher_as)
  allRes_as$p.adj = p.adjust(allRes_as$classicFisher, "fdr")
  
  # GO analysis for Cellular Component
  GOdata_as_cc <- new("topGOdata",
                      ontology = "CC",
                      allGenes = antisense,
                      nodeSize = 5,
                      annot = annFUN.org,
                      mapping = "org.Mm.eg.db",
                      ID = "ensembl")
  
  resultFisher_as_cc <- runTest(GOdata_as_cc, algorithm = "classic", statistic = "fisher")
  allRes_as_cc <- GenTable(GOdata_as_cc, classicFisher = resultFisher_as_cc)
  allRes_as_cc$p.adj = p.adjust(allRes_as_cc$classicFisher, "fdr")
  
  # GO for Molecular Function
  GOdata_as_mf <- new("topGOdata",
                      ontology = "MF",
                      allGenes = antisense,
                      nodeSize = 5,
                      annot = annFUN.org,
                      mapping = "org.Mm.eg.db",
                      ID = "ensembl")
  
  resultFisher_as_mf <- runTest(GOdata_as_mf, algorithm = "classic", statistic = "fisher")
  allRes_as_mf <- GenTable(GOdata_as_mf, classicFisher = resultFisher_as_mf)
  allRes_as_mf$p.adj = p.adjust(allRes_as_mf$classicFisher, "fdr")
  print("GO itself done")
  
  # merge tables
  allRes_as$Type<-"Biological Process"
  allRes_as_cc$Type<-"Cellular Component"
  allRes_as_mf$Type<-"Molecular Function"
  allRes_as_full<-rbind(allRes_as,allRes_as_cc,allRes_as_mf)
  allRes_as_full$FDR<-p.adjust(allRes_as_full$classicFisher, "fdr")
  print("fdr adjusted")
  return(allRes_as_full)
}

#   allRes_as_full$Term<-as.factor(allRes_as_full$Term)
#   allRes_as_full$Term<-ordered(allRes_as_full$Term, levels=as.character(allRes_as_full$Term))

Hypo_allRes_tested<-myGO(Hypo_tested)
fwrite(Hypo_allRes_tested,file="/data/natalia/EPO/results/Ribo_HypovsPlacebo_above3baseMean_softfilt_GO.txt",sep="\t",row.names=TRUE)
# 27 DNA-binding transcription activator activity, RNA polymerase II-specific
# 30 DNA-binding transcription activator activity
Hypo_allRes_tested$Term[c(27,30)]<-c("DNA-bind.transcr.activator activity, RNA PolII","DNA-bind.transcription activator activity")
Hypo_allRes_tested$Term<-as.factor(Hypo_allRes_tested$Term)
Hypo_allRes_tested$Term<-ordered(Hypo_allRes_tested$Term, levels=as.character(Hypo_allRes_tested$Term))

EPO_allRes_tested<-myGO(EPO_tested)
fwrite(EPO_allRes_tested,file="/data/natalia/EPO/results/Ribo_EPOvsPlacebo_above3baseMean_softfilt_GO.txt",sep="\t",row.names=TRUE)

# 27 transmitter-gated monoatomic ion channel activity involved in regulation of postsynaptic membrane potential
# 30 transmitter-gated monoatomic ion channel activity
EPO_allRes_tested$Term[c(27,30)]<-c("transmitter-gated m.ion ch.activity,regulation of postsynaptic membrane","transmitter-gated monoatomic ion channel activity")
EPO_allRes_tested$Term<-as.factor(EPO_allRes_tested$Term)
EPO_allRes_tested$Term<-ordered(EPO_allRes_tested$Term, levels=as.character(EPO_allRes_tested$Term))

ggplot(Hypo_allRes_tested, aes(x=Term, y=Significant, fill=Type)) +
  geom_bar(stat="identity")+coord_flip()+ 
  scale_fill_viridis(discrete=TRUE, direction=1)+
  labs(fill="", y="Number of genes", x="Significantly Enriched Terms")+theme_light()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


ggplot(EPO_allRes_tested, aes(x=Term, y=Significant, fill=Type)) +
  geom_bar(stat="identity")+coord_flip()+ 
  scale_fill_viridis(discrete=TRUE, direction=1)+
  labs(fill="", y="Number of genes", x="Significantly Enriched Terms")+theme_light()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


nrow(Total_analysis[Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="CRWvsPlacebo",]) # [1] 2208
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="CRWvsPlacebo") & Total_analysis$log2FoldChange>0,]) # [1] 1196
nrow(Total_analysis[Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="HypoxiavsPlacebo",]) # [1] 2999
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="HypoxiavsPlacebo") & Total_analysis$log2FoldChange>0,]) # [1] 1671
nrow(Total_analysis[Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="EPOvsPlacebo",]) # [1] 1616
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="EPOvsPlacebo") & Total_analysis$log2FoldChange>0,]) # [1] 963



nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="CRWvsPlacebo") & abs(Total_analysis$log2FoldChange)>=0.5,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="CRWvsPlacebo") & Total_analysis$log2FoldChange>=0.5,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="HypoxiavsPlacebo") & abs(Total_analysis$log2FoldChange)>=0.5,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="HypoxiavsPlacebo") & Total_analysis$log2FoldChange>=0.5,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="EPOvsPlacebo") & abs(Total_analysis$log2FoldChange)>=0.5,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="EPOvsPlacebo") & Total_analysis$log2FoldChange>=0.5,]) 


nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="CRWvsPlacebo") & abs(Total_analysis$log2FoldChange)>=1,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="CRWvsPlacebo") & Total_analysis$log2FoldChange>=1,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="HypoxiavsPlacebo") & abs(Total_analysis$log2FoldChange)>=1,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="HypoxiavsPlacebo") & Total_analysis$log2FoldChange>=1,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="EPOvsPlacebo") & abs(Total_analysis$log2FoldChange)>=1,]) 
nrow(Total_analysis[(Total_analysis$FDR_joint<=0.05 & Total_analysis$Group=="EPOvsPlacebo") & Total_analysis$log2FoldChange>=1,]) 
