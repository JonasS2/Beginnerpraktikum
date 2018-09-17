# reading in the data matrix and convert to an DGEList 
library(edgeR)


# d <- read.csv(file = "/Volumes/prj/internship_jschimmer/stringtie/gene_count_matrix.csv")

# data <- readRDS("/Users/jonas/Desktop/Analysis-master/rds/Freichel.rds")

data <- readRDS("rds/Freichel.rds")

library(biomaRt)
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="dec2016.archive.ensembl.org")
resMArt=getBM(attributes=c("ensembl_gene_id","description","mgi_symbol"),mart=mart)
resMArt2=getBM(attributes=c("ensembl_gene_id","gene_biotype"),mart=mart)

#/prj/micropeptides20/micropeptides/tables
RNASeqSFB= readRDS("rds/Freichel.rds")
own_ss <- readRDS("rds/own_sample_sheet.rds")

##############################
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
masterTable <- RNASeqSFB[helper]


rownames(masterTable)=RNASeqSFB[,1]


samp=read.delim("sample_sheets/polyA_sample_sheet.txt",header=F) # use of this additional file is quicker
# samp=read.delim("./Desktop/Analysis-master/sample_sheets/polyA_sample_sheet.txt",header=F)

cond=rep("Control",ncol(masterTable))
cond[grep("KO",samp[,1])]<-"KO"
cond=factor(cond)

treat=rep("saline",ncol(masterTable))
treat[grep("renol",samp[,1])]<-"Iso"
treat=factor(treat)

ens2symb=resMArt[,3]
names(ens2symb)=resMArt[,1]

selID=c(1:ncol(masterTable)) #shortcut 

prefix="Wu4_experiment_";

masterTable=masterTable[,selID]                   # does nothing
cond=droplevels(cond[selID])                      # drops unused levels in factor cond (nothing in this script)
treat=relevel(droplevels(treat[selID]),"saline")  # drops unused levels first and set "saline" then as base level

require(edgeR)

z<-DGEList(counts=masterTable, remove.zeros=T)
z<-calcNormFactors(z) 

fourFactor=factor(paste(cond,treat,sep="_"))


design <- model.matrix(~0+cond*treat)     # has interaction term and is used only for it later
design2 <- model.matrix(~0+fourFactor)    # four factor without interaction term
z2<-z;


#design 
z <- estimateGLMCommonDisp(z,design)  #estimates common dispersion 
z <- estimateGLMTrendedDisp(z,design) #estimates abundance dispersion trend
z <- estimateGLMTagwiseDisp(z,design) #estimates bayes tagwise dispersion
fit <- glmQLFit(z,design)             #fits the estimated dispersions to the count data

#design2
z2 <- estimateGLMCommonDisp(z2,design2)
z2 <- estimateGLMTrendedDisp(z2,design2)
z2 <- estimateGLMTagwiseDisp(z2,design2) 
fit2 <- glmQLFit(z2,design2)          #fits the estimated dispersions to the count data

#first contrast "conti" compares KO vs Control
Conti<-makeContrasts((fourFactorKO_Iso+fourFactorKO_saline)-(fourFactorControl_Iso+fourFactorControl_saline),levels=design2)
#second contrast "conti2" compares Iso vs Saline (i swapped saline and iso because saline is control treatment)
Conti2<-makeContrasts((fourFactorControl_saline+fourFactorKO_saline)-(fourFactorControl_Iso+fourFactorKO_Iso),levels=design2)

################################






v <- voom(z, design, plot=TRUE) 


vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, coef=4) #TODO
###check if the contrast is correct!!

efit <- eBayes(vfit)

## /Users/jonas/Desktop/Analysis-master
png('./Figures_edgeR/voom_Plot_interaction.png',width = 40, height = 20, units = "cm", res=90)
plotSA(efit)
dev.off()

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)



# Examination of the number of DE Genes
summary(decideTests(efit))

png('./Figures_edgeR/vennDiagramm_interaction.png',width = 40, height = 20, units = "cm", res=90)
vennDiagram(dt[,1], circle.col=c("turquoise", "salmon"))
dev.off()

de.common <- which(dt[,1]!=0)
length(de.common)
write.fit(tfit, dt, file="results.txt")

# Examination of individual DE Genes fro top to bottom:
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))

library(Glimma)
#anno = z$counts


#TODO save the html somwhere
glMDPlot(tfit, anno = z$counts, coef=1, status=dt, main=colnames(tfit)[1],
         id.column="ENTREZID", counts=z$counts ,groups = design, 
         path = "./Figures_edgeR", folder = "html",
         html = "glMDPlot_interaction", launch=TRUE)
