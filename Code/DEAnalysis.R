# creating DESeq2 dataset with appropriate design and MA plot the data

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

#only the stim subset is used

deseqMat <-
  DESeqDataSetFromMatrix(
    countData = data[2:17],
    colData = own_ss[1:16, ],
    design = ~ Genotype + Treatment + Genotype:Treatment
  )

ppp <- DESeq(deseqMat)

res <- results(ppp)

plotMA(res, main="DESeq2", ylim=c(-2,2))

results(ppp,
list(
"Genotype_Stim1.2_KO_vs_Stim1.2_control",
"GenotypeStim1.2_KO.Treatmentsaline")
)
############### comments and notes

#resultsNames(ppp)

################################
# Control Treatment effect
res <- results(ppp, contrast = c("Treatment", "saline", "Angiotensin"))

png('./Figures/MA_Plot/Control_treatmenteffect.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()

#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Control_Treatmenteffect.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")

#KO treatment effect
res <- results(ppp,list( contrast = c("Treatment_saline_vs_Angiotensin", "GenotypeStim1.2_KO.Treatmentsaline")))
png('./Figures/MA_Plot/KO_treatmenteffect.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/KO_Treatmenteffect.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")


#difference between treatment effect on genotypes
res <- results(ppp, name = "GenotypeStim1.2_KO.Treatmentsaline")
png('./Figures/MA_Plot/Genotype_treatmenteffectdiff.png',width = 40, height = 20, units = "cm", res=90)
plotMA(res, main="DESeq2")#, ylim=c(-2,2))
dev.off()
#merge with biomaRT and save as excel file
trans <- as.data.frame(res)
trans$"gene_id" <- data$gene_id
trans <- trans[,c("gene_id","baseMean","log2FoldChange","lfcSE", "stat","pvalue","padj")]
total <- merge(trans,fck,by="gene_id")
write.xlsx(total, "./Excel_files/Genotype_Treatmenteffectdiff.xlsx", sheetName="Sheet1",colNames = TRUE, borders = "columns")

