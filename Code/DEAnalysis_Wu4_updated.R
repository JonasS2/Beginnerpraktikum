# creating DESeq2 dataset with appropriate design and MA plot the data

#to run first to lines in Rstudio: own_ss <- readRDS("./Desktop/Analysis-master/rds/own_sample_sheet.rds")
# data <- readRDS("./Desktop/Analysis-master/rds/Freichel.rds")
own_ss <- readRDS("rds/own_sample_sheet.rds")
data <- readRDS("rds/Freichel.rds")
library(DESeq2)
library(openxlsx)

# loading and preparing gene alignement
library(biomaRt)
ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
fck <- getBM(attributes = c("ensembl_gene_id","description"), filters="ensembl_gene_id",values=data$gene_id,mart = ensembl)
colnames(fck)[1] <- "gene_id"

#renaming columns
x <- own_ss$data_Names
own_ss[1] <- rownames(own_ss)
rownames(own_ss) <- x

#only the Wu4 subset is used with treatments "saline" = control and "isopproterenol"

o <- grepl("Wu4",own_ss$Genotype)
u <- grepl("p",own_ss$p_t)
s <- grepl("saline",own_ss$Treatment)
q <- grepl("Isoproterenol",own_ss$Treatment)
w <- q | s
m <- u & o & w
helper <-  logical(1) 
helper[2:97] <- m

#factors are combined in treatment column to prepare DE
own_ss$Treatment <- factor(paste0(own_ss$Treatment,own_ss$Genotype))
##############

##############

deseqMat <-
  DESeqDataSetFromMatrix(
    countData = data[helper],
    colData = own_ss[m, ],
    design = ~ Treatment
  )
ppp <- DESeq(deseqMat)

res <- results(ppp)

plotMA(res, main="DESeq2", ylim=c(-2,2))

##############

#resultsNames(ppp)

#################


# Effect of Wu4_KO (WT_saline vs Wu4Ko_saline)
res <- results(ppp, contrast = c("Treatment", "salineWu4_control", "salineWu4_KO"))

png('./Figures/MA_Plot_updated/Wu4/WT_salinevsWu4KO_saline.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files_updated/Wu4/WT_salinevsWu4KO_saline.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



# Effect of Isoprot in WT animals (WT_saline vs WT_Isoprot)
res <- results(ppp, contrast = c("Treatment", "salineWu4_control", "IsoproterenolWu4_control"))

png('./Figures/MA_Plot_updated/Wu4/WT_salinevsWT_Isoprot.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files_updated/Wu4/WT_salinevsWT_Isoprot.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



# Effect of Wu4_KO (WT_isopoterenol vs Wu4Ko_isoproterenol)
res <- results(ppp, contrast = c("Treatment", "IsoproterenolWu4_control", "IsoproterenolWu4_KO"))

png('./Figures/MA_Plot_updated/Wu4/WT_IsoproterenolvsWu4KO_Isoproterenol.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files_updated/Wu4/WT_IsoproterenolvsWu4KO_Isoproterenol.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")



# Effect of Isoprot in Wu4KO animals (Wu4KO_saline vs Wu4KO_Isoprot)
res <- results(ppp, contrast = c("Treatment", "salineWu4_KO", "IsoproterenolWu4_KO"))

png('./Figures/MA_Plot_updated/Wu4/Wu4KO_salinevsWu4KO_Isoproterenol.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files_updated/Wu4/Wu4KO_salinevsWu4KO_Isoproterenol.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")