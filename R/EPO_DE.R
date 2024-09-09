# EPO polyA
# direct cDNAseq from Manuela-Chiara project, bulk, human fcx

# dir.create("/data/natalia/R")
myPaths <- .libPaths()

myPaths <- c("/data/natalia/R", myPaths)

.libPaths(myPaths)
library(Rsubread)
library(data.table)
library(stringr)
library(tidyr)
library(purrr)
library(dplyr)
library(topGO)
library(RSQLite)
# "org.Mm.eg.db"
library(org.Mm.eg.db)
library(viridis)
library(ggplot2)
library(ggsignif)
library(ggpubr)



# GO analysis
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
  allRes_as_full$Term<-as.factor(allRes_as_full$Term)
  allRes_as_full$Term<-ordered(allRes_as_full$Term, levels=as.character(allRes_as_full$Term))
  return(allRes_as_full)
}



# PCA function
pca_dt<-function(vst_tab,ntopvar=20000){
  rv <- rowVars(assay(vst_tab))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntopvar, length(rv)))]
  pca2 <- prcomp(t(assay(vst_tab)[select, ]))
  percentVar <- pca2$sdev^2/sum(pca2$sdev^2)

  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d2 <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2], PC3=pca2$x[, 3], PC4 = pca2$x[, 4],
                   PC5 = pca2$x[, 5], PC6 = pca2$x[, 6],  
                   PC7 = pca2$x[, 7], PC8 = pca2$x[, 8], PC9=pca2$x[, 9], PC10 = pca2$x[, 10],
                   PC11 = pca2$x[, 11], PC12 = pca2$x[, 12], PC13=pca2$x[, 13], PC14 = pca2$x[, 14])
  attr(d2, "percentVar") <- percentVar[1:14]
  return(d2)
}


# function to load based on succesfull strategy
# generate single count table per bam file
load_stringtie_cnt<-function(bam_file, gtf_file,outdir){
  df<-featureCounts(files = bam_file,annot.ext=gtf_file,
                    GTF.featureType="exon", isGTFAnnotationFile=TRUE,GTF.attrType="gene_id",useMetaFeatures=TRUE,
                    strandSpecific=2,isPairedEnd=TRUE,countMultiMappingReads=FALSE)
  df_cnt<- as.data.frame(df$counts) 
  df_cnt$Gene<-rownames(df_cnt)
  print(paste("Completed sample",colnames(df_cnt)[1]))
  fwrite(df_cnt,paste0(outdir,colnames(df_cnt)[1],"_raw_counts.txt"),sep="\t",quote=F,row.names=T)
  return(df_cnt)
}

# function to run through all of them one by one

bam_load<-function(mypath, gtf_file,outdir){
  temp<-list.files(path=mypath, pattern="*.bam")
  temp<-paste0(mypath,temp)
  myfiles<-lapply(temp,load_stringtie_cnt, gtf_file,outdir)
  gene_tab<-myfiles %>% purrr::reduce(full_join, by="Gene")
  fwrite(gene_tab,paste0(outdir,"Joint_gene_raw_cnts.txt"),sep="\t",quote=F,row.names=F)
  return(gene_tab)
}


dir.create("/data/natalia/EPO/counts/")

polyA_raw_cnts<-bam_load("/data/natalia/EPO/bam/", "/data/natalia/EPO/gencode.vM35.annotation.gtf", "/data/natalia/EPO/counts/")

rownames(polyA_raw_cnts)<-polyA_raw_cnts$Gene
polyA_raw_cnts<-polyA_raw_cnts[,-2]
colnames(polyA_raw_cnts)<-str_remove(colnames(polyA_raw_cnts),"_merged_MM_100Aligned.sortedByCoord.out_deduplicated.bam")
colnames(polyA_raw_cnts)<-str_remove( colnames(polyA_raw_cnts), "_S[0-9]*")
# meta table 

full_meta_epo <- read.delim("/data/natalia/EPO/full_meta_epo.txt")
full_meta_epo<-full_meta_epo[full_meta_epo$Sequencing=="totalRNA",c("id", "qbic","secondaryName", "Tag", "i.id", "mouse_id",
                                                                    "eppendorf_code", "group","i.Tag")]

full_meta_epo<-full_meta_epo[match(colnames(polyA_raw_cnts),full_meta_epo$id),]


rownames(full_meta_epo)<-full_meta_epo$id
identical(rownames(full_meta_epo),colnames(polyA_raw_cnts))
# [1] TRUE
# remove chromosome Y genes
gencode_vM35 <- read.delim("/data/natalia/EPO/gencode.vM35.annotation.gtf", header=FALSE, comment.char="#")
gencode_vM35<-gencode_vM35[gencode_vM35$V1=="chrY",]
gencode_vM35<-separate(gencode_vM35,V9,sep="; ",into=c("gene","rest"),extra="merge")
gencode_vM35$gene<-str_remove(gencode_vM35$gene,"gene_id ")

polyA_raw_cnts_woY<-polyA_raw_cnts[!(rownames(polyA_raw_cnts) %in% gencode_vM35$gene),]
full_meta_epo$test<-paste(full_meta_epo$Tag,full_meta_epo$group,"_")
library(DESeq2)

dds<-DESeqDataSetFromMatrix(countData=polyA_raw_cnts_woY,
                            colData=full_meta_epo,
                            design=~test)
dds<-estimateSizeFactors(dds)
sizeFactors(dds)
idx.nz <- apply(counts(dds), 1, function(x) { all(x > 0)})
# 23067a001_01   23067a002_01   23067a003_01   23067a004_01   23067a005_01   23067a006_01   23067a007_01   23067a008_01 
# 0.8192041      0.4807162      0.5578605      0.7372265      0.7130288      0.7474391      0.8845816      0.6669549 
# 23067a009_01   23067a010_01   23067a011_01   23067a012_01   23067a013_01   23067a014_01   23067a015_01   23067a016_01 
# 0.8406966      0.8060564      0.5654573      0.4942186      0.7080679      0.4923836      0.5747959      0.5390389 
# 23067a017_01   23067a018_01   23067a019_01   23067a020_01   23067a021_01   23067a022_01   23067a023_01   23067a024_01 
# 0.7880003      0.4808060      0.8523397      0.9106477      0.8639203      0.6656781      0.6278267      0.7215727 
# 23067a048L2_01 23067a049L2_01 23067a050L2_01 23067a051L2_01 23067a052L2_01 23067a053L2_01 23067a054L2_01 23067a055L2_01 
# 1.8659482      1.2511899      1.6517732      1.4389504      1.6209988      1.0000107      1.5247734      1.4775343 
# 23067a056L2_01 23067a057L2_01 23067a058L2_01 23067a059L2_01 23067a060L2_01 23067a061L2_01 23067a062L2_01 23067a063L2_01 
# 1.5542609      1.5443767      1.7139064      1.6024711      1.4010599      1.2815565      1.2030293      1.3276769 
# 23067a064L2_01 23067a065L2_01 23067a066L2_01 23067a067L2_01 23067a068L2_01 23067a069L2_01 23067a070L2_01   23067a094_01 
# 1.5322270      1.4092229      1.7049953      1.6476738      1.6650151      1.6673013      1.5949443      1.5539000


dir.create("/data/natalia/EPO/results")
library(geneplotter)
png("/data/natalia/EPO/results/genes_ECDF.png", width = 900, height = 1200)
multiecdf( counts(dds, normalized = T)[idx.nz ,],
           xlab="mean counts", xlim=c(0, 1000))
dev.off()


png("/data/natalia/EPO/results/genes_density.png", width = 900, height = 1200)
multidensity( counts(dds, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000))
dev.off()

pdf("/data/natalia/EPO/results/genes_pairwiseMAs.pdf")
MA.idx = t(combn(1:48, 2))
for( i in 1:48){
  plotMD(counts(dds, normalized = T)[idx.nz ,],c(MA.idx[i,1],MA.idx[i,2]),
         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3) )
}
dev.off()




# rlog takes too much time and resources, thus VST

vst_dds<-vst(dds, blind=TRUE, nsub=8000)


pca_dt<-function(vst_tab,ntopvar=20000){
  rv <- rowVars(assay(vst_tab))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntopvar, length(rv)))]
  pca2 <- prcomp(t(assay(vst_tab)[select, ]))
  percentVar <- pca2$sdev^2/sum(pca2$sdev^2)
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d2 <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2], PC3=pca2$x[, 3], PC4 = pca2$x[, 4],
                   PC5 = pca2$x[, 5], PC6 = pca2$x[, 6],  
                   PC7 = pca2$x[, 7], PC8 = pca2$x[, 8], PC9=pca2$x[, 9], PC10 = pca2$x[, 10],
                   PC11 = pca2$x[, 11], PC12 = pca2$x[, 12], PC13=pca2$x[, 13], PC14 = pca2$x[, 14])
  attr(d2, "percentVar") <- percentVar[1:14]
  return(d2)
}

pca_dds<-pca_dt(vst_dds, ntopvar=8000)
pca_dds[,15:16]<-full_meta_epo[,c("Tag","group")]
pca_dds$itag<-full_meta_epo$i.Tag
fwrite(pca_dds,"/data/natalia/EPO/results/genes_top8000_PCA.txt",sep="\t",row.names=TRUE)

library(cowplot)
library(ggplot2)
library(umap)

library(data.table)
library(tidyr)
library(stringr)
# "#154360", "#FF5733", "#FFC300"
p1<-ggplot(pca_dds, aes(x=PC1,y=PC2)) +theme_minimal() +
  geom_point(aes(fill=group,color=group,shape=itag),size = 8, alpha = 0.85)  +
  xlab(paste0("PC1",": ", round(attributes(pca_dds)$percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2",": ", round(attributes(pca_dds)$percentVar[2] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))


p2<-ggplot(pca_dds, aes(x=PC3,y=PC4)) +theme_minimal() +
  geom_point(aes(fill=group,color=group,shape=itag),size = 8, alpha = 0.85)  +
  xlab(paste0("PC3",": ", round(attributes(pca_dds)$percentVar[3] * 100), "% variance")) +
  ylab(paste0("PC4",": ", round(attributes(pca_dds)$percentVar[4] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p3<-ggplot(pca_dds, aes(x=PC5,y=PC6)) +theme_minimal() +
  geom_point(aes(fill=group,color=group,shape=itag),size = 8, alpha = 0.85)  +
  xlab(paste0("PC5",": ", round(attributes(pca_dds)$percentVar[5] * 100), "% variance")) +
  ylab(paste0("PC6",": ", round(attributes(pca_dds)$percentVar[6] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

# change design formula to have tag (i.tag) and group separately; test by i.tag 
dds$group<-factor(dds$group,levels=c("Placebo","EPO","Hypoxia","CRW"))
dds$i.Tag<-factor(dds$i.Tag,levels=c("TotalRNA","RiboTag_IAA"))
design(dds)<- ~ group + i.Tag

# Getting ribotag (neuronal) markers
dds<-estimateDispersions(dds)

png("/data/natalia/EPO/results/genes_dispersion.png", width = 1300, height = 1000)
plotDispEsts(dds)
dev.off()
dds_nc<-counts(dds, normalized=TRUE)
fwrite(as.data.frame(dds_nc),"/data/natalia/EPO/results/genes_MedRat_norm.txt", row.names = TRUE, sep="\t")

dds <- nbinomWaldTest(dds, maxit=500)
# total as control; ribo as test
dds_res<-results(dds, contrast=,c("i.Tag","RiboTag_IAA","TotalRNA"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds_res),"/data/natalia/EPO/results/Total_vs_RiboTag_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/Total_vs_RiboTag_genes_unfilt_MA.png", width = 1300, height = 1000)
plotMA( dds_res, ylim = c(-1, 1) )
dev.off()

png("/data/natalia/EPO/results/Total_vs_RiboTag_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds_res_2 <- dds_res[ !is.na(dds_res$padj), ]
dds_res_2 <- dds_res_2[ !is.na(dds_res_2$pvalue), ]
png("/data/natalia/EPO/results/Total_vs_RiboTag_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds_res_2$pvalue, breaks=20, col="grey" )
dev.off()

png("/data/natalia/EPO/results/Total_vs_RiboTag_genes_filt_MA.png", width = 1300, height = 1000)
plotMA( dds_res_2, ylim = c(-1, 1) )
dev.off()
dds_res_2[,"FDR"]<-p.adjust(dds_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds_res_2),file="/data/natalia/EPO/results/Total_vs_RiboTag_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)
dds_res_2<-as.data.frame(dds_res_2)
# get IDs for "markers" (sign.upregulated) and characterize functionally the markers?
neuro_up_genes<-rownames(dds_res_2[(dds_res_2$log2FoldChange>=1 & dds_res_2$FDR<0.05) & dds_res_2$baseMean>=3,])
length(neuro_up_genes)
# [1] 1146
neuro_up_genes_ext<-rownames(dds_res_2[dds_res_2$log2FoldChange>=1 & dds_res_2$FDR<0.05,])
length(neuro_up_genes_ext)
# [1] 1980
# GO analysis against all genes on the list
all_genes<-rownames(dds_res_2[dds_res_2$baseMean>=3,])
length(all_genes)
# [1] 18446
all_genes_ext<-rownames(dds_res_2)
# [1] 26215
##############################################################################
##############################################################################

# GO setup


# Edit all genes and tested genes to drop suffix after .
# gsub("\\..*","",ipd_dif_exp_isos_down)
all_genes<-gsub("\\..*","",all_genes)
all_genes_ext<-gsub("\\..*","",all_genes_ext)
neuro_up_genes<-gsub("\\..*","",neuro_up_genes)
neuro_up_genes_ext<-gsub("\\..*","",neuro_up_genes_ext)

# create vector for testing

tested <- factor(as.integer(all_genes %in% neuro_up_genes))
names(tested) <- all_genes

tested_ext <- factor(as.integer(all_genes_ext %in% neuro_up_genes_ext))
names(tested_ext) <- all_genes_ext

# GO

allRes_tested<-myGO(tested)
fwrite(allRes_tested,file="/data/natalia/EPO/results/NeuronsUp_above3baseMean_GO.txt",sep="\t",row.names=TRUE)

allextRes_tested<-myGO(tested_ext)
fwrite(allextRes_tested,file="/data/natalia/EPO/results/NeuronsUp_extended_GO.txt",sep="\t",row.names=TRUE)


# plot
ggplot(allRes_tested, aes(x=Term, y=Significant, fill=Type)) +
  geom_bar(stat="identity")+coord_flip()+ 
  scale_fill_viridis(discrete=TRUE, direction=1)+
  labs(fill="", y="Number of genes", x="Significantly Enriched Terms")+theme_light()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())





allextRes_tested$Term<-as.factor(allextRes_tested$Term)
allextRes_tested<-allextRes_tested[1:29,]
allextRes_tested$Term<-ordered(allextRes_tested$Term, levels=as.character(allextRes_tested$Term))

ggplot(allextRes_tested, aes(x=Term, y=Significant, fill=Type)) +
  geom_bar(stat="identity")+coord_flip()+ 
  scale_fill_viridis(discrete=TRUE, direction=1)+
  labs(fill="", y="Number of genes", x="Significantly Enriched Terms")+theme_light()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# now get full DE lists and plot; then we can run DE per tag-group, pairwise comparisons,correction for multiple testings

neuro_de_genes<-rownames(dds_res_2[(abs(dds_res_2$log2FoldChange)>=1 & dds_res_2$FDR<0.05) & dds_res_2$baseMean>=3,])
neuro_de_genes<-gsub("\\..*","",neuro_de_genes)
tested_de <- factor(as.integer(all_genes %in% neuro_de_genes))
names(tested_de) <- all_genes

allDERes_tested<-myGO(tested_de)
fwrite(allDERes_tested,file="/data/natalia/EPO/results/NeuronsDE_above3baseMean_GO.txt",sep="\t",row.names=TRUE)

ggplot(allDERes_tested, aes(x=Term, y=Significant, fill=Type)) +
  geom_bar(stat="identity")+coord_flip()+ 
  scale_fill_viridis(discrete=TRUE, direction=1)+
  labs(fill="", y="Number of genes", x="Significantly Enriched Terms")+theme_light()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# neuron, glia, etc markers (from single cell) plot violin plots per total vs neuron groups

# Mouse cell markers
# Microglia: Iba1, Cx3cr1, Csf1r (ENSMUSG00000024397, ENSMUSG00000052336, ENSMUSG00000024621)
# Astrocyte: Gfap, Aqp4, S100b, Sox9 (?), Slc1a3 (ENSMUSG00000020932 )
# Neuron: NeuN, Map2, Dcx, Gad1, Stmn2 (ENSMUSG00000025576 , ENSMUSG00000025576 , ENSMUSG00000027500 )
# Oligodendrocytes: Olig2, Sox10, Ptgds (ENSMUSG00000039830,ENSMUSG00000033006 )
# Vessel cell: Cd31 ENSMUSG00000020717
# Endothelial cell: Cldn5, Flt1, CD31 ENSMUSG00000041378 

# dds_nc
# "ENSMUSG00000024397.15","ENSMUSG00000052336.9", "ENSMUSG00000020932.15", "ENSMUSG00000025576.18", 
# "ENSMUSG00000027500.11", "ENSMUSG00000015222.19", "ENSMUSG00000039830.10", "ENSMUSG00000033006.11",
# "ENSMUSG00000020717.20", "ENSMUSG00000041378.3"

dds_nc_marker<-dds_nc[rownames(dds_nc) %in% c("ENSMUSG00000024397.15","ENSMUSG00000052336.9", "ENSMUSG00000020932.15", "ENSMUSG00000025576.18", 
                                              "ENSMUSG00000027500.11", "ENSMUSG00000015222.19", "ENSMUSG00000039830.10", "ENSMUSG00000033006.11",
                                              "ENSMUSG00000020717.20", "ENSMUSG00000041378.3"),]
dds_nc_marker<-as.data.frame(t(dds_nc_marker))
dds_nc_marker[,11:12]<-full_meta_epo[,c("i.Tag","group")]
# plots


ggboxplot(dds_nc_marker, x = "i.Tag", y = "ENSMUSG00000025576.18",
          color = "i.Tag", palette = "jco")+ 
  stat_compare_means(comparisons = list(c("RiboTag_IAA", "TotalRNA")),label = "p.signif", size = 8)

for (i in c("ENSMUSG00000024397.15","ENSMUSG00000052336.9", "ENSMUSG00000020932.15", "ENSMUSG00000025576.18", 
            "ENSMUSG00000027500.11", "ENSMUSG00000015222.19", "ENSMUSG00000039830.10", "ENSMUSG00000033006.11",
            "ENSMUSG00000020717.20", "ENSMUSG00000041378.3")){
  p<-ggboxplot(dds_nc_marker, x = "i.Tag", y = i,
            color = "i.Tag", palette = "jco")+ 
    stat_compare_means(comparisons = list(c("RiboTag_IAA", "TotalRNA")),label = "p.signif", size = 8)
  png(paste0("/data/natalia/EPO/results/",i,"_Total_vs_RiboTag.png"), width = 700, height = 700)
  print(p)
  dev.off()
}

# Run comparisons

dds2<-DESeqDataSetFromMatrix(countData=polyA_raw_cnts_woY[,full_meta_epo[full_meta_epo$i.Tag=="RiboTag_IAA","id"]],
                            colData=full_meta_epo[full_meta_epo$i.Tag=="RiboTag_IAA",],
                            design=~group)
dds2<-estimateSizeFactors(dds2)
sizeFactors(dds2)
# 23067a001_01 23067a002_01 23067a003_01 23067a004_01 23067a005_01 23067a006_01 23067a007_01 23067a008_01 23067a009_01 
# 1.2299819    0.7191061    0.8311900    1.0935059    1.0602534    1.1047516    1.3123012    1.0013435    1.2520161 
# 23067a010_01 23067a011_01 23067a012_01 23067a013_01 23067a014_01 23067a015_01 23067a016_01 23067a017_01 23067a018_01 
# 1.1962560    0.8378231    0.7373819    1.0584348    0.7375861    0.8598548    0.8103913    1.1703713    0.7207587 
# 23067a019_01 23067a020_01 23067a021_01 23067a022_01 23067a023_01 23067a024_01 
# 1.2522482    1.3494013    1.2731530    0.9838075    0.9241044    1.0786069 
vst_dds2<-vst(dds2, blind=TRUE, nsub=8000)

pca_dds2<-pca_dt(vst_dds2, ntopvar=8000)
pca_dds2[,15:16]<-full_meta_epo[full_meta_epo$i.Tag=="RiboTag_IAA",c("i.Tag","group")]
p1<-ggplot(pca_dds2, aes(x=PC1,y=PC2)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC1",": ", round(attributes(pca_dds2)$percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2",": ", round(attributes(pca_dds2)$percentVar[2] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p2<-ggplot(pca_dds2, aes(x=PC3,y=PC4)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC3",": ", round(attributes(pca_dds2)$percentVar[3] * 100), "% variance")) +
  ylab(paste0("PC4",": ", round(attributes(pca_dds2)$percentVar[4] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p3<-ggplot(pca_dds2, aes(x=PC5,y=PC6)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC5",": ", round(attributes(pca_dds2)$percentVar[5] * 100), "% variance")) +
  ylab(paste0("PC6",": ", round(attributes(pca_dds2)$percentVar[6] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p4<-ggplot(pca_dds2, aes(x=PC7,y=PC8)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC7",": ", round(attributes(pca_dds2)$percentVar[7] * 100), "% variance")) +
  ylab(paste0("PC8",": ", round(attributes(pca_dds2)$percentVar[8] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

png("/data/natalia/EPO/results/PC12_RiboTag.png", width = 700, height = 700)
print(p1)
dev.off()

png("/data/natalia/EPO/results/PC34_RiboTag.png", width = 700, height = 700)
print(p2)
dev.off()

png("/data/natalia/EPO/results/PC56_RiboTag.png", width = 700, height = 700)
print(p3)
dev.off()

png("/data/natalia/EPO/results/PC78_RiboTag.png", width = 700, height = 700)
print(p4)
dev.off()

# de
fwrite(pca_dds2,"/data/natalia/EPO/results/ribotag_genes_top8000_PCA.txt",sep="\t",row.names=TRUE)
dds2<-estimateDispersions(dds2)
dds2_nc<-counts(dds2, normalized=TRUE)
fwrite(as.data.frame(dds2_nc),"/data/natalia/EPO/results/ribotag_genes_MedRat_norm.txt", row.names = TRUE, sep="\t")

dds2 <- nbinomWaldTest(dds2, maxit=500)

# comparisons epo vs placebo, crw vs placebo, hypoxia vs placebo

# placebo as control; epo as test
dds2_res<-results(dds2, contrast=,c("group","EPO","Placebo"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds2_res),"/data/natalia/EPO/results/EPOvsPlacebo_RiboTag_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/EPOvsPlacebo_RiboTag_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds2_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds2_res_2 <- dds2_res[ !is.na(dds2_res$padj), ]
dds2_res_2 <- dds2_res_2[ !is.na(dds2_res_2$pvalue), ]
png("/data/natalia/EPO/results/EPOvsPlacebo_RiboTag_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds2_res_2$pvalue, breaks=20, col="grey" )
dev.off()

dds2_res_2[,"FDR"]<-p.adjust(dds2_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds2_res_2),file="/data/natalia/EPO/results/EPOvsPlacebo_RiboTag_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)



# placebo as control; crw as test
dds2_res<-results(dds2, contrast=,c("group","CRW","Placebo"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds2_res),"/data/natalia/EPO/results/CRWvsPlacebo_RiboTag_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/CRWvsPlacebo_RiboTag_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds2_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds2_res_2 <- dds2_res[ !is.na(dds2_res$padj), ]
dds2_res_2 <- dds2_res_2[ !is.na(dds2_res_2$pvalue), ]
png("/data/natalia/EPO/results/CRWvsPlacebo_RiboTag_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds2_res_2$pvalue, breaks=20, col="grey" )
dev.off()

dds2_res_2[,"FDR"]<-p.adjust(dds2_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds2_res_2),file="/data/natalia/EPO/results/CRWvsPlacebo_RiboTag_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)

# placebo as control; hypoxia as test
dds2_res<-results(dds2, contrast=,c("group","Hypoxia","Placebo"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds2_res),"/data/natalia/EPO/results/HypoxiavsPlacebo_RiboTag_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/HypoxiavsPlacebo_RiboTag_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds2_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds2_res_2 <- dds2_res[ !is.na(dds2_res$padj), ]
dds2_res_2 <- dds2_res_2[ !is.na(dds2_res_2$pvalue), ]
png("/data/natalia/EPO/results/HypoxiavsPlacebo_RiboTag_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds2_res_2$pvalue, breaks=20, col="grey" )
dev.off()

dds2_res_2[,"FDR"]<-p.adjust(dds2_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds2_res_2),file="/data/natalia/EPO/results/HypoxiavsPlacebo_RiboTag_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)


########################################################
# TOTAL

dds3<-DESeqDataSetFromMatrix(countData=polyA_raw_cnts_woY[,full_meta_epo[full_meta_epo$i.Tag=="TotalRNA","id"]],
                             colData=full_meta_epo[full_meta_epo$i.Tag=="TotalRNA",],
                             design=~group)
dds3<-estimateSizeFactors(dds3)
sizeFactors(dds3)

vst_dds3<-vst(dds3, blind=TRUE, nsub=8000)

pca_dds3<-pca_dt(vst_dds3, ntopvar=8000)
pca_dds3[,15:16]<-full_meta_epo[full_meta_epo$i.Tag=="TotalRNA",c("i.Tag","group")]
p1<-ggplot(pca_dds3, aes(x=PC1,y=PC2)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC1",": ", round(attributes(pca_dds3)$percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2",": ", round(attributes(pca_dds3)$percentVar[2] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p2<-ggplot(pca_dds3, aes(x=PC3,y=PC4)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC3",": ", round(attributes(pca_dds3)$percentVar[3] * 100), "% variance")) +
  ylab(paste0("PC4",": ", round(attributes(pca_dds3)$percentVar[4] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p3<-ggplot(pca_dds3, aes(x=PC5,y=PC6)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC5",": ", round(attributes(pca_dds3)$percentVar[5] * 100), "% variance")) +
  ylab(paste0("PC6",": ", round(attributes(pca_dds3)$percentVar[6] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

p4<-ggplot(pca_dds3, aes(x=PC7,y=PC8)) +theme_minimal() +
  geom_point(aes(fill=group,color=group),size = 8, alpha = 0.85)  +
  xlab(paste0("PC7",": ", round(attributes(pca_dds3)$percentVar[7] * 100), "% variance")) +
  ylab(paste0("PC8",": ", round(attributes(pca_dds3)$percentVar[8] * 100), "% variance")) + coord_fixed() + 
  theme(
    text=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position="bottom") + 
  scale_fill_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))+ 
  scale_color_manual(values=c("#154360", "#FF5733", "#FFC300","grey"))

png("/data/natalia/EPO/results/PC12_TotalRNA.png", width = 700, height = 700)
print(p1)
dev.off()

png("/data/natalia/EPO/results/PC34_TotalRNA.png", width = 700, height = 700)
print(p2)
dev.off()

png("/data/natalia/EPO/results/PC56_TotalRNA.png", width = 700, height = 700)
print(p3)
dev.off()

png("/data/natalia/EPO/results/PC78_TotalRNA.png", width = 700, height = 700)
print(p4)
dev.off()

# de
fwrite(pca_dds3,"/data/natalia/EPO/results/TotalRNA_genes_top8000_PCA.txt",sep="\t",row.names=TRUE)
dds3<-estimateDispersions(dds3)
dds3_nc<-counts(dds3, normalized=TRUE)
fwrite(as.data.frame(dds3_nc),"/data/natalia/EPO/results/TotalRNA_genes_MedRat_norm.txt", row.names = TRUE, sep="\t")

dds3 <- nbinomWaldTest(dds3, maxit=500)

# comparisons epo vs placebo, crw vs placebo, hypoxia vs placebo

# placebo as control; epo as test
dds3_res<-results(dds3, contrast=,c("group","EPO","Placebo"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds3_res),"/data/natalia/EPO/results/EPOvsPlacebo_TotalRNA_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/EPOvsPlacebo_TotalRNA_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds3_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds3_res_2 <- dds3_res[ !is.na(dds3_res$padj), ]
dds3_res_2 <- dds3_res_2[ !is.na(dds3_res_2$pvalue), ]
png("/data/natalia/EPO/results/EPOvsPlacebo_TotalRNA_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds3_res_2$pvalue, breaks=20, col="grey" )
dev.off()

dds3_res_2[,"FDR"]<-p.adjust(dds3_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds3_res_2),file="/data/natalia/EPO/results/EPOvsPlacebo_TotalRNA_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)



# placebo as control; crw as test
dds3_res<-results(dds3, contrast=,c("group","CRW","Placebo"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds3_res),"/data/natalia/EPO/results/CRWvsPlacebo_TotalRNA_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/CRWvsPlacebo_TotalRNA_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds3_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds3_res_2 <- dds3_res[ !is.na(dds3_res$padj), ]
dds3_res_2 <- dds3_res_2[ !is.na(dds3_res_2$pvalue), ]
png("/data/natalia/EPO/results/CRWvsPlacebo_TotalRNA_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds3_res_2$pvalue, breaks=20, col="grey" )
dev.off()

dds3_res_2[,"FDR"]<-p.adjust(dds3_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds3_res_2),file="/data/natalia/EPO/results/CRWvsPlacebo_TotalRNA_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)

# placebo as control; hypoxia as test
dds3_res<-results(dds3, contrast=,c("group","Hypoxia","Placebo"), pAdjustMethod = "BH")

fwrite(as.data.frame(dds3_res),"/data/natalia/EPO/results/HypoxiavsPlacebo_TotalRNA_genes_DEresults_full.txt",sep="\t",row.names=TRUE)

png("/data/natalia/EPO/results/HypoxiavsPlacebo_TotalRNA_genes_unfilt_hist_pval.png", width = 1300, height = 1000)
hist( dds3_res$pvalue, breaks=20, col="grey" )
dev.off()
# "bimodal", but disappears when I remove NA padj (very few counts genes driven)
dds3_res_2 <- dds3_res[ !is.na(dds3_res$padj), ]
dds3_res_2 <- dds3_res_2[ !is.na(dds3_res_2$pvalue), ]
png("/data/natalia/EPO/results/HypoxiavsPlacebo_TotalRNA_genes_filt_hist_pval.png", width = 1300, height = 1000)
hist( dds3_res_2$pvalue, breaks=20, col="grey" )
dev.off()

dds3_res_2[,"FDR"]<-p.adjust(dds3_res_2$pvalue, method = "BH")
fwrite(as.data.frame(dds3_res_2),file="/data/natalia/EPO/results/HypoxiavsPlacebo_TotalRNA_genes_DEresults_filt.txt",sep="\t",row.names=TRUE)
