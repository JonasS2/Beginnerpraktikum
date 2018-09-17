
rds/Freichel.rds:
	Rscript ./Code/Load.R > ./log/load.log
	
hcluster: rds/Freichel.rds SampleSheet
	Rscript ./Code/hcluster.R > ./log/$@.txt
	
SampleSheet: rds/Freichel.rds
	Rscript ./Code/SampleSheet.R > ./log/sheet.log
	
Heatmap: SampleSheet
	Rscript ./Code/QC.R > ./log/QC.log
	
PCA: SampleSheet
	Rscript ./Code/PCA.R > ./log/PCA.log
	Rscript ./Code/ggplot.R > ./log/ggplotPCA.log
	
DEAnalysis: SampleSheet
	Rscript ./Code/DEAnalysis.R > ./log/DEAnalysis.log
	
DEAnalysis_Wu4: SampleSheet
	Rscript ./Code/DEAnalysis_Wu4.R > ./log/DEAnalysis_Wu4.log
	
DE_edgeR_Wu4: SampleSheet
	Rscript ./Code/DE_edgeR_Wu4.R > ./log/DE_edgeR_Wu4.log
	
Gene_ontology_KOvsCo_Wu4: SampleSheet
	Rscript ./Code/Gene_ontology_KOvsCo_Wu4.R > ./log/Gene_ontology_KOvsCo_Wu4.log
	
Gene_ontology_SalvsIso_Wu4: SampleSheet
	Rscript ./Code/Gene_ontology_SalvsIso_Wu4.R > ./log/Gene_ontology_SalvsIso_Wu4.log
	
Gene_ontology_interaction_Wu4: SampleSheet
	Rscript ./Code/Gene_ontology_interaction_Wu4.R > ./log/Gene_ontology_interaction_Wu4.log