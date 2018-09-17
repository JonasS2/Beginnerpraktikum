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


# write sheet
lrt <- glmQLFTest(fit, coef=4) # coef= 4 gets the interaction term

lrt2 <- glmQLFTest(fit2, contrast=Conti) # does the actual calculations for KO vs Control 
lrt3 <- glmQLFTest(fit2, contrast=Conti2)# does the actual calculations for Iso vs Saline

dfOutGenotype<-topTags(lrt2,nrow(fit2$coefficients),p.value=1) # extracts the top DE tags
dfOutGenotype=data.frame(ID=rownames(dfOutGenotype),dfOutGenotype)

dfOutTreatment<-topTags(lrt3,nrow(fit2$coefficients),p.value=1) # why is it named Ribo old script cp
dfOutTreatment=data.frame(ID=rownames(dfOutTreatment),dfOutTreatment)

dfOut<-topTags(lrt,nrow(fit$coefficients),p.value=1)
dfOut=data.frame(ID=rownames(dfOut),dfOut)

#add all results in one table, KOvsControl, Iso vs Saline and interaction term
dfOut=merge(dfOut,dfOutGenotype[,c(1,2,6)],by.x=1,by.y=1,all.x=T)
colnames(dfOut)[ncol(dfOut)-1]="KOvsWTlog2fc"
colnames(dfOut)[ncol(dfOut)]="KOvsWT.FDR"


dfOut=merge(dfOut,dfOutTreatment[,c(1,2,6)],by.x=1,by.y=1,all.x=T) 
colnames(dfOut)[ncol(dfOut)-1]="SalvsIsolog2fc"
colnames(dfOut)[ncol(dfOut)]="SalvsIso.FDR"

dfOut=merge(dfOut,resMArt,by.x=1,by.y=1,all.x=T)
dfOut=merge(dfOut,data.frame(ID=rownames(cpm(z)),cpm(z)),by.x=1,by.y=1)
dfOut=dfOut[order(dfOut[,"FDR.x"]),]

dfOutGenotype=merge(dfOutGenotype,resMArt,by.x=1,by.y=1,all.x=T)
dfOutGenotype=dfOutGenotype[order(dfOutGenotype[,"FDR"]),]

dfOutTreatment=merge(dfOutTreatment,resMArt,by.x=1,by.y=1,all.x=T)
dfOutTreatment=dfOutTreatment[order(dfOutTreatment[,"FDR"]),]

#exit(0);
#Fold change RNA-seq / RIbo-seq und Interaction

require(openxlsx);

wb <- createWorkbook()

addWorksheet(wb, sheetName = paste0(prefix,"_Int"));
writeDataTable(wb, sheet = 1, x= dfOut);
addWorksheet(wb, sheetName = paste0(prefix,"_KO"));
writeDataTable(wb, sheet = 2, x= dfOutGenotype);
addWorksheet(wb, sheetName = paste0(prefix,"_Iso"));
writeDataTable(wb, sheet = 3, x= dfOutTreatment);

saveWorkbook(wb, paste0(prefix,"_May2018.xlsx"), overwrite = TRUE)

pdf(paste0(prefix,"May2018.pdf"))
#smoothScatter(dfOut[,c("RNAlog2fc","Ribolog2fc")],main=paste0(prefix," FDR 1%"))

#abline(a=0,b=1,col="yellow")
#points(dfOut[which(dfOut[,"FDR.x"]<0.01),c("RNAlog2fc","Ribolog2fc")],col="red")
#text(dfOut[which(dfOut[,"FDR.x"]<0.01),c("RNAlog2fc","Ribolog2fc")],dfOut[which(dfOut[,"FDR.x"]<0.01),"mgi_symbol"],cex=0.7)
dev.off()

