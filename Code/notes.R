strsplit(names(d)[1],split = "_")
[[1]]
[1] "gene" "id"  

> strsplit(names(d)[2],split = "_")
[[1]]
[1] "ATII.26" ""        "E01"     "p"       "S1"      "L001"   

> pste(strsplit(names(d)[1],split = "_"))
Error: could not find function "pste"
> paste(strsplit(names(d)[1],split = "_"))
[1] "c(\"gene\", \"id\")"
> paste(strsplit(names(d)[1],split = "_")[[1]])
[1] "gene" "id"  
> paste(strsplit(names(d)[2],split = "_")[[1]])
[1] "ATII.26" ""        "E01"     "p"       "S1"      "L001"   
> paste(strsplit(names(d)[2],split = "_")[[1]], collapse="__")
[1] "ATII.26____E01__p__S1__L001"
> paste(strsplit(names(d)[2],split = "_")[[1]], collapse="_")
[1] "ATII.26__E01_p_S1_L001"
> paste(strsplit(names(d)[2],split = "_")[[1]][1:4], collapse="_")
[1] "ATII.26__E01_p"
merge()


tmp <- grep("ATII", d,value = TRUE) # get values of grep

f[[16]][3]


###################### own stuff



getBM(attributes = "description", filters="ensembl_gene_id",values=data$gene_id[1],mart = ensembl)


getBM(attributes = c("ensembl_gene_id",description"), filters="ensembl_gene_id",values=data$gene_id[1],mart = ensembl)
Error: unexpected string constant in "getBM(attributes = c("ensembl_gene_id",description"), filters=""
> getBM(attributes = c("ensembl_gene_id","description"), filters="ensembl_gene_id",values=data$gene_id[1],mart = ensembl)
     ensembl_gene_id                                                                      description



View(own_sample_sheet)
> d <- own_sample_sheet$namework
> grep(ATII,d)
Error in grep(ATII, d) : object 'ATII' not found
> grep(d,ATII)
Error in grep(d, ATII) : object 'ATII' not found
> ?grep
> tmp <- grep("ATII", d)
> tmp
[1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
> tmp2 <- names(tmp)
> tmp2
NULL
> tmp <- names(grep("ATII", d))
> tmp
NULL
> View(sample_sheet)
> View(own_sample_sheet)
> ?grep
> tmp <- grep("ATII", d,value = TRUE)
> tmp



###################
ERROR: configuration failed for package ‘XML’
* removing ‘/home/user/R/x86_64-pc-linux-gnu-library/3.2/XML’
ERROR: dependency ‘XML’ is not available for package ‘annotate’
* removing ‘/home/user/R/x86_64-pc-linux-gnu-library/3.2/annotate’
ERROR: dependency ‘annotate’ is not available for package ‘genefilter’
* removing ‘/home/user/R/x86_64-pc-linux-gnu-library/3.2/genefilter’
ERROR: dependency ‘annotate’ is not available for package ‘geneplotter’
* removing ‘/home/user/R/x86_64-pc-linux-gnu-library/3.2/geneplotter’
ERROR: dependencies ‘genefilter’, ‘geneplotter’ are not available for package ‘DESeq2’
* removing ‘/home/user/R/x86_64-pc-linux-gnu-library/3.2/DESeq2’




######################
o <- grepl("Stim",own_ss$Info_treatment_group)
plotPCA(rld[,o],"Info_treatment_group")
relevel(own_ss$Info_treatment_group, "TRPC1/TRPC4-KO-sham")


png('./Figures/PCA.png',width = 40, height = 20, units = "cm", res=90)
colData(rld)$Contr <- "Treat"
colData(rld)$Contr[grepl("Contr",colData(rld)$Info_treatment_group)] <- "Contr"

plotPCA(rld,"Info_treatment_group")

dev.off()

plotPCA(rld[,own_ss$p_t == "p"],"p_t")

plotPCA(rld,"Contr")
#############################################
> q + geom_text(aes(label=colData(rld)$V3[x]))
> q + geom_text(aes(label=colData(rld)$V3[x])) + facet_wrap(~colData(rld)$Treatment)
Error in `$<-.data.frame`(`*tmp*`, "PANEL", value = c(3L, 3L, 3L, 3L,  : 
                                                        replacement has 96 rows, data has 80
> q + geom_text(aes(label=colData(rld)$V3[x])) + facet_wrap(~colData(rld)$Treatment[x])
> install.packages("plotly")

ggplot(df[o,], aes(
  x = PC1,
  y = PC2,
  colour = colData(rld[,o])$Treatment
))  + geom_point() + xlab("PC1") + ylab("PC2")

pr <- prcomp(log2(1 + data[2:97]))

> plot(pr$x[,1], pr$x[,2], col=factor(colData(rld)$Treatment))
> counts(ddsMat)

##########################
df <- as_data_frame(q$rotation)

ggplot(df,aes(x=df$PC1,df$PC2))  + geom_point()
##########################


q <- prcomp(t(log2(1 + (data[,o]))))
df <- as.data.frame(q$x)

ggplot(data = df[,]) + 
  geom_point(aes(x=PC1,y=PC2))   
+ xlab("PC1") + ylab("PC2") 

################
res$padj < 1e-5

# 4 figures arranged in 2 rows and 2 columns
par(mfrow=c(1,3))
plot(main="Figures/MA_Plot/KO_treatmenteffect")
plot(main="Figures/MA_Plot/Control_treatmenteffect")
hist(main="Figures/MA_Plot/Genotype_treatmenteffectdiff")



