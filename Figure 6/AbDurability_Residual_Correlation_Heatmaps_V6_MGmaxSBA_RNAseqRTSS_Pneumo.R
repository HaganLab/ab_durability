#Antibody response longevity correlation analysis - cross-vaccine
#Original data
#6/11/23
library(tidyverse)
library(data.table)
library(meta)
library(Biobase)
library(R.matlab)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(fgsea)
library(org.Hs.eg.db)
library(RColorBrewer)
library(DESeq2)

rm(list=ls())
setwd("/Users/tlhagan/Documents/Stanford/Antibody Durability/")
#Parameters
Vax010_tp=28 #TP to use for Vax010
Immunity_tp=7 #TP to use for Immunity
CHI_tp=28 #TP to use for CHI
Pfizer_tp=28 #TP to use for Pfizer data
YF_tp=7 #TP to use for YF studies
MG_tp=7 #TP to use for MG studies
malaria_tp=62 #TP to us for RTSS
collapse_option='meanFC'
p_cut=0.3
corrMethod='pearson'
#Load BTMs
#load('~/Documents/Emory/BTM_files/BTM_for_GSEA_20131008_geneID.RData') #Entrez
BTM_list=readRDS("/Users/tlhagan/Documents/Stanford/HIPC 2/Signatures Project/Analysis/Virtual Study/btm_list_2020_12_23.rds") #Symbol
#Load BTM groups
BTM_groups=read.delim('/Users/tlhagan/Documents/Emory/BTM_files/BTM functional groups V3.txt', na.strings=c(''))
#Match names with list
#BTM_groups=BTM_groups[match(str_remove(names(BTM_list), ' -.*'),BTM_groups$BTM),] #Entrez
BTM_groups=BTM_groups[match(sub('.*\\(S','S',sub('.*\\(M','M',names(BTM_list))),sub('.*\\(S','S',sub('.*\\(M','M',BTM_groups$NAME))),] #Symbol
names(BTM_list)=BTM_groups$NAME
# #Optional:Remove BTMs without subgroup
# BTM_list=BTM_list[sapply(BTM_groups$SUBGROUP, function(x) nchar(x)>0)]
# BTM_groups=BTM_groups[sapply(BTM_groups$SUBGROUP, function(x) nchar(x)>0),]
#Set BTM group as factor and manually set order
BTM_groups$SUBGROUP=factor(BTM_groups$SUBGROUP,
                           levels=c('ANTIGEN PRESENTATION', 'INFLAMMATORY/TLR/CHEMOKINES','INTERFERON/ANTIVIRAL SENSING',
                                    'MONOCYTES', 'DC ACTIVATION', 'NEUTROPHILS', 'SIGNAL TRANSDUCTION', 'ECM AND MIGRATION',
                                    'PLATELETS', 'ENERGY METABOLISM', 'CELL CYCLE', 'NK CELLS', 'T CELLS', 'B CELLS', 'PLASMA CELLS'))
#Sort BTM list by subgroup
BTM_list=BTM_list[order(BTM_groups$SUBGROUP)]
BTM_groups=BTM_groups[order(BTM_groups$SUBGROUP),]
rownames(BTM_groups)=BTM_groups$NAME
#Set BTM colors
cols=colorRampPalette(brewer.pal(length(levels(BTM_groups$SUBGROUP)), "Set1"))
mycolors_BTM=setNames(cols(length(levels(BTM_groups$SUBGROUP))), levels(BTM_groups$SUBGROUP))

#Load Data
#Vax010
#Load GE data
gene_exp=read.table('/Users/tlhagan/Documents/Emory/Vax010/Datasets/Vax010/Vax010.RMA.all.samples.txt', header=T)
rownames(gene_exp)=str_remove(rownames(gene_exp), '_PM')
#Collapse to gene symbol (remove probes matching multiple gene, take probe with highest mean exp per gene)
#Load lookup
probe_symbol_lookup=read.delim('/Users/tlhagan/Documents/Emory/Scripts and tables/AffyU133plus2_Mar2016_probe_entrez_symbol_desc.txt', sep='\t')
probe_symbol_lookup=probe_symbol_lookup[!is.na(probe_symbol_lookup$Symbol),]
#Add symbol and mean exp columns
gene_exp$mean=rowMeans(gene_exp)
gene_exp$symbol=probe_symbol_lookup$Symbol[match(rownames(gene_exp),probe_symbol_lookup$Probe)]
gene_exp=gene_exp[!is.na(gene_exp$symbol),]
#Collapse by highest mean per gene symbol
gene_exp=gene_exp %>% group_by(symbol) %>% dplyr::slice(which.max(mean))
gene_exp=data.frame(gene_exp)
rownames(gene_exp)=gene_exp$symbol
gene_exp=subset(gene_exp, select=-c(symbol,mean))
#Create metadata
metadata=data.frame(Subject=strtrim(colnames(gene_exp),7), Day=as.numeric(str_remove(gsub(".*[.]([^.]+)[.].*", "\\1", colnames(gene_exp)),'Day')))
rownames(metadata)=colnames(gene_exp)
#Create eset
eset=ExpressionSet(as.matrix(gene_exp), AnnotatedDataFrame(metadata))
rm(gene_exp,metadata)
#Compute D0/21 normalized FC (log2)
tp=sort(setdiff(unique(eset$Day),c(0,21)))
ind=lapply(tp, function(x) which(eset$Day==x))
exp_FC=vector("list",length(tp))
for (i in 1:(length(tp))) {
  if (tp[i]<21){
    ind_D0=which(0==eset$Day)
  }
  else {
    ind_D0=which(21==eset$Day)
  }
  common=intersect(eset$Subject[ind[[i]]],eset$Subject[ind_D0])
  ia=na.omit(match(common,eset$Subject[ind[[i]]]))
  ib=na.omit(match(common,eset$Subject[ind_D0]))
  exp_FC[[i]]=ExpressionSet(exprs(eset[,ind[[i]][ia]])-exprs(eset[,ind_D0[ib]]), eset@phenoData[ind[[i]][ia],])
}
#Keep only tp of interest
Vax010_exp_FC=exp_FC[[which(tp==Vax010_tp)]]
#Get HAI Residuals/FC
Vax010_D100_resid=read.table('/Users/tlhagan/Documents/Emory/Vax010/Analysis/Longevity/Vax010_Adj_D100D42_HAI_Residual_w_D100_D42_HAI_MN_select.txt', header=TRUE, sep='\t')
rownames(Vax010_D100_resid)=Vax010_D100_resid$Subject
#Find matching subjects
common=intersect(Vax010_exp_FC$Subject, Vax010_D100_resid$Subject)
ind_exp=na.omit(match(common,Vax010_exp_FC$Subject))
ind_resid=na.omit(match(common,Vax010_D100_resid$Subject))
Vax010_exp_FC=Vax010_exp_FC[,ind_exp]
Vax010_D100_resid=Vax010_D100_resid[ind_resid,]
rownames(Vax010_D100_resid)=colnames(Vax010_exp_FC)
#Create expressionset
Vax010_exp_FC=ExpressionSet(exprs(Vax010_exp_FC),phenoData = AnnotatedDataFrame(Vax010_D100_resid))
rm(eset,tp,ind,exp_FC,ind_D0,ia,ib,common,Vax010_D100_resid,ind_exp,ind_resid)

#Immunity
#Load GE data
gene_exp=read.table('/Users/tlhagan/Documents/Emory/Vax010/Datasets/Immunity/Normalized_Immunity.txt', header=T, row.names=1)
#Collapse to gene symbol (remove probes matching multiple gene, take probe with highest mean exp per gene)
#Load lookup
probe_symbol_lookup=read.delim('/Users/tlhagan/Documents/Emory/Scripts and tables/AffyU133plus2_Mar2016_probe_entrez_symbol_desc.txt', sep='\t')
probe_symbol_lookup=probe_symbol_lookup[!is.na(probe_symbol_lookup$Symbol),]
#Add symbol and mean exp columns
gene_exp$mean=rowMeans(gene_exp)
gene_exp$symbol=probe_symbol_lookup$Symbol[match(rownames(gene_exp),probe_symbol_lookup$Probe)]
gene_exp=gene_exp[!is.na(gene_exp$symbol),]
#Collapse by highest mean per gene symbol
gene_exp=gene_exp %>% group_by(symbol) %>% dplyr::slice(which.max(mean))
gene_exp=data.frame(gene_exp)
rownames(gene_exp)=gene_exp$symbol
gene_exp=subset(gene_exp, select=-c(symbol,mean))
#Create metadata
metadata=data.frame(Subject=str_remove(str_remove(colnames(gene_exp),'_.*'),'X'),
                    Day=as.numeric(gsub(".*_D([^.]+)", "\\1", colnames(gene_exp))))
rownames(metadata)=colnames(gene_exp)
#Create eset
eset=ExpressionSet(as.matrix(gene_exp), AnnotatedDataFrame(metadata))
rm(gene_exp,metadata)
#Compute D0 normalized FC (log2)
tp=sort(setdiff(unique(eset$Day),0))
ind=lapply(tp, function(x) which(eset$Day==x))
exp_FC=vector("list",length(tp))
for (i in 1:(length(tp))) {
  ind_D0=which(0==eset$Day)
  common=intersect(eset$Subject[ind[[i]]],eset$Subject[ind_D0])
  ia=na.omit(match(common,eset$Subject[ind[[i]]]))
  ib=na.omit(match(common,eset$Subject[ind_D0]))
  exp_FC[[i]]=ExpressionSet(exprs(eset[,ind[[i]][ia]])-exprs(eset[,ind_D0[ib]]), eset@phenoData[ind[[i]][ia],])
}
#Keep only tp of interest
Immunity_exp_FC=exp_FC[[which(tp==Immunity_tp)]]
#Load HAI data
Immunity_HAI_raw=readMat('/Users/tlhagan/Documents/Emory/Vax010/Datasets/Immunity/Immunity_D180_D30_MaxFC_submatch.mat')
Immunity_HAI_raw$D180.sub.labels=unlist(Immunity_HAI_raw$D180.sub.labels)
Immunity_HAI_raw$D30.D180.Max.FC.labels=unlist(Immunity_HAI_raw$D30.D180.Max.FC.labels)
Immunity_HAI=data.frame('pt_id'=Immunity_HAI_raw$D180.sub.labels,
                        'D30_MaxFC'=Immunity_HAI_raw$D30.D180.Max.FC[1,],
                        'D180_FC_D30MaxStrain'=Immunity_HAI_raw$D30.D180.Max.FC[2,])
rownames(Immunity_HAI)=Immunity_HAI$pt_id
#Compute HAI residual
ab_lm=lm(D180_FC_D30MaxStrain ~ D30_MaxFC, data=Immunity_HAI)
ab_resid=resid(ab_lm)
ab_resid=data.frame(pt_id=names(ab_resid), HAI_D180_resid=ab_resid)
ab_titers=left_join(Immunity_HAI,ab_resid)
rownames(ab_titers)=ab_titers$pt_id
#Find matching subjects
common=intersect(Immunity_exp_FC$Subject, ab_titers$pt_id)
ind_exp=na.omit(match(common,Immunity_exp_FC$Subject))
ind_resid=na.omit(match(common,ab_titers$pt_id))
Immunity_exp_FC=Immunity_exp_FC[,ind_exp]
ab_titers=ab_titers[ind_resid,]
rownames(ab_titers)=colnames(Immunity_exp_FC)
#Create expressionset
Immunity_exp_FC=ExpressionSet(exprs(Immunity_exp_FC), AnnotatedDataFrame(ab_titers))
rm(Immunity_HAI,Immunity_HAI_raw)

#CHI
#Load GE data
CHI=readMat('/Users/tlhagan/Documents/Emory/Vax010/Datasets/CHI/CHI_gexp_genes_FC.mat')
#Get D28D21 FC data
CHI_exp_FC=CHI$exp.FC[[4]][[1]]
dimnames(CHI_exp_FC)=list(unlist(CHI$probelist),str_replace(unlist(CHI$sub.labels.FC[[4]])[seq(2, ncol(CHI_exp_FC)*3+1, 3)],'_','.'))
#Load HAI data
CHI_HAI=read.delim('/Users/tlhagan/Documents/Emory/Vax010/Datasets/CHI/H5N1_CHI_HAI_Adj_median_edit.txt')
#Convert to log2
CHI_HAI$D42_HAI=log2(CHI_HAI$D42_HAI)
CHI_HAI$D100_HAI=log2(CHI_HAI$D100_HAI)
#Compute HAI residuals
fit=lm(D100_HAI ~ D42_HAI, data=CHI_HAI)
CHI_HAI$HAI_resid=fit$residuals
rownames(CHI_HAI)=CHI_HAI$Subject
#Find matching subjects
common=intersect(colnames(CHI_exp_FC), CHI_HAI$Subject)
ind_exp=na.omit(match(common,colnames(CHI_exp_FC)))
ind_resid=na.omit(match(common,CHI_HAI$Subject))
#Create expressionset
CHI_exp_FC=ExpressionSet(as.matrix(CHI_exp_FC[,ind_exp]), AnnotatedDataFrame(CHI_HAI[ind_resid,]))
rm(CHI, CHI_HAI, fit, ind_exp, ind_resid)

#Pfizer
#Load data
RNAseq=readRDS('/Users/tlhagan/Documents/Stanford/COVID/Vaccine Trial/RNAseq/dat_mat_v2_corrected.rds')
pheno=readRDS('/Users/tlhagan/Documents/Stanford/COVID/Vaccine Trial/RNAseq/pheno_v6_CORRECTED.RDS')
#Create eset
pheno$fullid=pheno$fullid_new
rownames(pheno)=pheno$fullid
eset=ExpressionSet(as.matrix(RNAseq[,-1]), AnnotatedDataFrame(pheno[match(colnames(RNAseq[,-1]),pheno$fullid),]))
rm(RNAseq, pheno)
#Create numeric day column
eset$day_num=as.numeric(str_replace_all(eset$day_agg, c('Day'='','BL'='0')))
#Exclude subject with SNP genotype mismatch
bad_sample_ID=c('PID2054_Day21_12.Jan')
eset=eset[,-match(bad_sample_ID,colnames(eset))]
#Normalize with DESeq
mat=DESeqDataSetFromMatrix(exprs(eset), pData(eset), ~1)
mat=estimateSizeFactors(mat)
mat_norm=counts(mat, normalized=TRUE)
#Log2 normalize
mat_norm=log2(mat_norm+1)
#Return to eset
Pfizer=ExpressionSet(as.matrix(mat_norm), eset@phenoData)
rm(eset,mat,mat_norm)
#Create numeric day column
Pfizer$day_num=as.numeric(str_replace_all(Pfizer$day_agg, c('Day'='','22.23'='22','BL'='0')))
#Load D210 Ab data, add to pData
ab_titers=read.delim('/Users/tlhagan/Documents/Stanford/COVID/Vaccine Trial/RNAseq/Longevity/PfizerCOVID_nAb_D0_D120_D210.txt')
ab_titers$pt_id=as.character(ab_titers$pt_id)
rownames(ab_titers)=ab_titers$pt_id
#Convert to log2
ab_titers[,2:6]=log2(ab_titers[,2:6])
#Fit D210/D42 titers with regression, compute residual
ab_lm=lm(nAb_D210 ~ nAb_D42, data=ab_titers)
ab_resid=resid(ab_lm)
ab_resid=data.frame(pt_id=names(ab_resid), nAb_D210_resid=ab_resid)
ab_titers=left_join(ab_titers,ab_resid)
ab_titers$pt_id=paste0('PID',ab_titers$pt_id)
temp=left_join(pData(Pfizer), ab_titers)
rownames(temp)=colnames(Pfizer)
pData(Pfizer)=temp

#Compute D0/21 normalized FC (log2)
tp=sort(setdiff(unique(Pfizer$day_num),c(0,21)))
ind=lapply(tp, function(x) which(Pfizer$day_num==x))
exp_FC=vector("list",length(tp))
for (i in 1:(length(tp))) {
  if (tp[i]<21){
    ind_D0=which(Pfizer$day_num==0)
  } else {
    ind_D0=which(Pfizer$day_num==21)
  }
  common=intersect(Pfizer$pt_id[ind[[i]]],Pfizer$pt_id[ind_D0])
  ia=na.omit(match(common,Pfizer$pt_id[ind[[i]]]))
  ib=na.omit(match(common,Pfizer$pt_id[ind_D0]))
  exp_FC[[i]]=ExpressionSet(exprs(Pfizer[,ind[[i]][ia]])-exprs(Pfizer[,ind_D0[ib]]), Pfizer@phenoData[ind[[i]][ia],])
}
#Keep only D7 (prime or boost)
Pfizer_exp_FC=exp_FC[[which(tp==Pfizer_tp)]]
rm(Pfizer,ab_lm,ab_resid,ab_titers,temp,tp,ind,exp_FC,ia,ib)

#Load YF studies
eset=readRDS("/Users/tlhagan/Documents/Stanford/HIPC 2/Signatures Project/Analysis/Virtual Study/2021_03_08/2021_03_08_young_noNorm_eset.rds")
colnames(eset)=eset$uid
#Create combined vaccine type/pathogen column
eset$pt=paste(eset$pathogen," (",eset$vaccine_type,")", sep='')
#Remove genes with NA
eset=eset[complete.cases(exprs(eset)),]
#Load YF studies ab data, convert titer to log2
peak_long_tp=list(sdy_1289=c(28,180), sdy_1294=c(28,84))
ab_sdy1289=read.delim('/Users/tlhagan/Documents/Stanford/HIPC 2/Signatures Project/Analysis/Virtual Study/Longevity/nab_SDY1289.tsv')
ab_sdy1289$Value.Preferred=log2(ab_sdy1289$Value.Preferred)
ab_sdy1294=read.delim('/Users/tlhagan/Documents/Stanford/HIPC 2/Signatures Project/Analysis/Virtual Study/Longevity/nab_SDY1294.tsv')
ab_sdy1294$Value.Preferred=log2(ab_sdy1294$Value.Preferred)
#Prune data frames
ab=list(ab_sdy1289, ab_sdy1294)
ab=lapply(1:length(ab), function(x) {y=ab[[x]][,c('Participant.ID','Study.Time.Collected','Value.Preferred')]; y})
rm(ab_sdy1289, ab_sdy1294)
ab_FC=vector("list",length(ab))
for (i in 1:length(ab)) {
  #Compute peak to long term FC
  ab_peak=ab[[i]] %>% dplyr::filter(Study.Time.Collected == peak_long_tp[[i]][1])
  ab_long=ab[[i]] %>% dplyr::filter(Study.Time.Collected == peak_long_tp[[i]][2])
  ab_FC[[i]]=inner_join(ab_peak,ab_long, by='Participant.ID')
  ab_FC[[i]]$FC=ab_FC[[i]]$Value.Preferred.y-ab_FC[[i]]$Value.Preferred.x
  #Compute peak to long term residual
  ab_lm=lm(Value.Preferred.y ~ Value.Preferred.x, data=ab_FC[[i]], na.action=na.exclude)
  ab_FC[[i]]$resid=resid(ab_lm)
}
#Merge and add to pData
ab_FC=do.call(rbind, ab_FC)
ab_FC=ab_FC[,c('Participant.ID','FC','resid')]
colnames(ab_FC)=c('participant_id','ab_peak_long_FC', 'ab_resid')
temp=pData(eset)
temp=left_join(temp,ab_FC, by='participant_id')
rownames(temp)=temp$uid
pData(eset)=temp
rm(temp)
#Remove subjects without residual
eset=eset[,!is.na(eset$ab_resid)]
#Find samples from timepoints of interest
tp_int=c(0,3,7,14)
#tp_int=unique(eset$study_time_collected[which(eset$study_time_collected>=0)]) #Alternate: use all timepoints (>=0)
ind=lapply(tp_int, function(x) which(eset$study_time_collected==x))
#Combine indices of all timepoints of interest
ind_all=Reduce(union,ind)
#Retain only samples from timepoints of interest
eset=eset[,ind_all]
#Recompute timepoint indices after removing extraneous timepoints
ind=lapply(tp_int, function(x) which(eset$study_time_collected==x))
#Create combined SDY/pathogen/vaccine type column
eset$SDY_pt=paste(eset$study,eset$pt)
#Create unique list of studies
matrix_uni=unique(eset$matrix)
#Compute D0 normalized FC
ind_D0=which(0==eset$study_time_collected)
common=lapply(2:length(ind),function(x) intersect(eset$participant_id[ind[[x]]],eset$participant_id[ind_D0]))
ia=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$participant_id[ind[[x]]])))
ib=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$participant_id[ind_D0])))
exp_FC=lapply(2:length(ind),function(x) eset[,ind[[x]][ia[[x-1]]]])
exp_FC=lapply(2:length(ind),function(x) {exprs(exp_FC[[x-1]])=exprs(exp_FC[[x-1]])-exprs(eset[,ind_D0[ib[[x-1]]]]); exp_FC[[x-1]]})
#Create separate esets for each HIPC study
SDY1289_exp_FC=exp_FC[[which(tp_int==YF_tp)-1]][,exp_FC[[which(tp_int==YF_tp)-1]]$matrix=='SDY1289_WholeBlood_MontrealCohort_Geo']
SDY1294_exp_FC=exp_FC[[which(tp_int==YF_tp)-1]][,exp_FC[[which(tp_int==YF_tp)-1]]$matrix=='SDY1294_PBMC_ChineseCohort_Geo']

#Load MG data
gene_exp=read.table('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Meni_SDY1260/Expression/gene_expression.txt', header=T)
rownames(gene_exp)=str_remove(rownames(gene_exp), '_PM')
#Collapse to gene symbol (remove probes matching multiple gene, take probe with highest mean exp per gene)
#Load lookup
probe_symbol_lookup=read.delim('/Users/tlhagan/Documents/Emory/Scripts and tables/AffyU133plus2_Mar2016_probe_entrez_symbol_desc.txt', sep='\t')
probe_symbol_lookup=probe_symbol_lookup[!is.na(probe_symbol_lookup$Symbol),]
#Add symbol and mean exp columns
gene_exp$mean=rowMeans(gene_exp)
gene_exp$symbol=probe_symbol_lookup$Symbol[match(rownames(gene_exp),probe_symbol_lookup$Probe)]
gene_exp=gene_exp[!is.na(gene_exp$symbol),]
#Collapse by highest mean per gene symbol
gene_exp=gene_exp %>% group_by(symbol) %>% dplyr::slice(which.max(mean))
gene_exp=data.frame(gene_exp)
rownames(gene_exp)=gene_exp$symbol
gene_exp=subset(gene_exp, select=-c(symbol,mean))
#Load metadata
metadata=read.table('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Meni_SDY1260/Samples/Annotation.txt', header=T)
rownames(metadata)=metadata$File
metadata=metadata[match(colnames(gene_exp),metadata$File),]
#Create eset
eset=ExpressionSet(as.matrix(gene_exp), AnnotatedDataFrame(metadata))
rm(gene_exp,metadata)
#Add titer information
ab_sdy1260=read.table('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Meni_SDY1260/Serology/SBA_serology_combined.txt', header=T)
days=unique(ab_sdy1260$Day)
#Remove GMT
ab_sdy1260$GMT_A_C=NULL
#Convert to log2
ab_sdy1260[,c('Titer_A','Titer_C')]=log2(ab_sdy1260[,c('Titer_A','Titer_C')])
#Compute FC
ab_sdy1260=pivot_wider(ab_sdy1260, id_cols = c('Subject'), values_from = c('Titer_A','Titer_C'), names_from=c('Day'))
ab_sdy1260[,grep('A',colnames(ab_sdy1260))]=ab_sdy1260[,grep('A',colnames(ab_sdy1260))]-ab_sdy1260$Titer_A_0
ab_sdy1260[,grep('C',colnames(ab_sdy1260))]=ab_sdy1260[,grep('C',colnames(ab_sdy1260))]-ab_sdy1260$Titer_C_0
#Remove D30 NAs
ab_sdy1260=ab_sdy1260[!(is.na(ab_sdy1260$Titer_A_30)|is.na(ab_sdy1260$Titer_C_30)),]
#Keep only data from strain with max FC
ab_titers=lapply(1:nrow(ab_sdy1260), function(x)
  if (ab_sdy1260$Titer_A_30[x]>ab_sdy1260$Titer_C_30[x]) {
    as.numeric(ab_sdy1260[x,grep('A',colnames(ab_sdy1260))])
  } else {
    as.numeric(ab_sdy1260[x,grep('C',colnames(ab_sdy1260))])
  })
ab_titers=as.data.frame(do.call(rbind, ab_titers))
colnames(ab_titers)=paste0('D',days)
ab_titers$Subject=ab_sdy1260$Subject
#Compute residual (each vaccine separately)
for (j in 1:length(unique(eset$Group))) {
  ab_lm=lm(as.formula('D180 ~ D30'),
           data=ab_titers[ab_titers$Subject %in% eset$Subject[eset$Group==unique(eset$Group)[j]],], na.action=na.exclude)
  ab_titers[which(ab_titers$Subject %in% eset$Subject[eset$Group==unique(eset$Group)[j]]),'ab_resid']=resid(ab_lm)
}
#Add to eset
temp=left_join(pData(eset), ab_titers, by='Subject')
rownames(temp)=temp$File
pData(eset)=temp
#Compute D0 normalized FC
ind_D0=which(0==eset$Day)
tp_int=unique(eset$Day)
ind=lapply(tp_int, function(x) which(eset$Day==x))
common=lapply(2:length(ind),function(x) intersect(eset$Subject[ind[[x]]],eset$Subject[ind_D0]))
ia=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$Subject[ind[[x]]])))
ib=lapply(2:length(ind),function(x) na.omit(match(common[[x-1]],eset$Subject[ind_D0])))
exp_FC=lapply(2:length(ind),function(x) eset[,ind[[x]][ia[[x-1]]]])
exp_FC=lapply(2:length(ind),function(x) {exprs(exp_FC[[x-1]])=exprs(exp_FC[[x-1]])-exprs(eset[,ind_D0[ib[[x-1]]]]); exp_FC[[x-1]]})
#Create separate esets for each vaccine
SDY1260_MCV4_exp_FC=exp_FC[[which(tp_int==MG_tp)-1]][,exp_FC[[which(tp_int==MG_tp)-1]]$Group=='MCV4']
SDY1260_MPSV4_exp_FC=exp_FC[[which(tp_int==MG_tp)-1]][,exp_FC[[which(tp_int==MG_tp)-1]]$Group=='MPSV4']

#Add malaria RRR
RTSS_exp_FC=readRDS('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/RTSS/RRR.RNAseq.respective.baselines.RDS')
tp_int=c(1,2,6,14,29,34,57,62)
RTSS_exp_FC=RTSS_exp_FC[[which(tp_int==malaria_tp)]]

#Add pneumococcal data
Pneumovax_exp_FC=readRDS('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Pneumococcal/Pneumovax.RDS')
Pneumovax_exp_FC=Pneumovax_exp_FC$Pneumovax.D7
Prevnar_exp_FC=readRDS('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Pneumococcal/Prevnar.RDS')
Prevnar_exp_FC=Prevnar_exp_FC$Prevnar.D7
# #Add pneumococcal data (top 20% serotypes residual)
# Pneumovax_exp_FC=readRDS('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Pneumococcal/Pneumovax.percentiles.RDS')
# Pneumovax_exp_FC=Pneumovax_exp_FC$Pneumovax.D7.20
# Prevnar_exp_FC=readRDS('/Users/tlhagan/Documents/Stanford/Antibody Durability/Human Datasets/Pneumococcal/Prevnar.percentiles.RDS')
# Prevnar_exp_FC=Prevnar_exp_FC$Prevnar.D7.20


#Merge into list
exp_FC=setNames(list(Vax010_exp_FC, Immunity_exp_FC, CHI_exp_FC, Pfizer_exp_FC,
                     SDY1289_exp_FC, SDY1294_exp_FC,
                     SDY1260_MCV4_exp_FC,
                     SDY1260_MPSV4_exp_FC,
                     RTSS_exp_FC,
#                     Pneumovax_exp_FC,
                     Prevnar_exp_FC),
                c('Vax010_D28','Immunity_D7','CHI_D28',paste0('Pfizer_D',Pfizer_tp),
                  paste0('SDY1289_YF_D',YF_tp),paste0('SDY1294_YF_D',YF_tp),
                  paste0('SDY1260_MCV4_SBAmaxAC_D',MG_tp),
                  paste0('SDY1260_MPSV4_SBAmaxAC_D',MG_tp),
                  paste0('RTSS_D236resid_D',malaria_tp),
#                  'Pneumovax_D7',
                  'Prevnar_D7'
                ))
#Standardize residuals
resid_list=c('D100_HAI_Residual','HAI_D180_resid', 'HAI_resid', 'nAb_D210_resid',
             'ab_resid', 'ab_resid',
             'ab_resid',
             'ab_resid',
             'Residual_D0_adj',
#             'resid',
             'resid')
for (i in 1:length(exp_FC)) {
  exp_FC[[i]]$ab_resid=pData(exp_FC[[i]])[,resid_list[i]]
}

#Subset datasets to common genes
common_gene=Reduce(intersect, lapply(exp_FC, function(x) rownames(x)))
# for (i in 1:length(exp_FC)) {
#   exp_FC[[i]]=ExpressionSet(exprs(exp_FC[[i]][na.omit(match(common_gene,rownames(exp_FC[[i]]))),]), exp_FC[[i]]@phenoData)
# }

#Remove BTMs with genes<cutoff in list
BTM_gene_cutoff=2
BTM_groups=BTM_groups[!sapply(BTM_list, function(x) length(intersect(x,common_gene)))<BTM_gene_cutoff,]
BTM_list=BTM_list[!sapply(BTM_list, function(x) length(intersect(x,common_gene)))<BTM_gene_cutoff]

# #Remove BTMs with no genes in list, order
# BTM_groups=BTM_groups[!sapply(BTM_list, function(x) is_empty(intersect(x,common_gene))),]
# BTM_list=BTM_list[!sapply(BTM_list, function(x) is_empty(intersect(x,common_gene)))]

#Collapse to BTM expression using meanFC or ssGSEA (depending on option)
if (collapse_option=='meanFC') {
  for (i in 1:length(exp_FC)) {
    exp_FC[[i]]=ExpressionSet(do.call(rbind, lapply(BTM_list, function(x) colMeans(exprs(exp_FC[[i]][na.omit(match(x,rownames(exp_FC[[i]]))),]),na.rm=TRUE))),
                              exp_FC[[i]]@phenoData)
  }
} else if (collapse_option=='ssGSEA') {
  for (i in 1:length(exp_FC)) {
    exp_FC[[i]]=ExpressionSet({z=apply(exprs(exp_FC[[i]]), 2, function(x) {y=fgsea(BTM_list, setNames(x,rownames(exp_FC[[i]]))); setNames(y$NES, y$pathway)});
    z}, exp_FC[[i]]@phenoData)
  }
}

#Find correlating features
ab_R=sapply(exp_FC, function(y) apply(exprs(y), 1, function(x) cor.test(x, y$ab_resid, method = corrMethod)$statistic))

#Plot BTMs related to platelets / cell adhesion
BTM_groups_int=c('PLATELETS', 'ECM AND MIGRATION')
for (i in 1:length(BTM_groups_int)) {
  df=ab_R[which(BTM_groups$SUBGROUP==BTM_groups_int[i]),]
  #df=ab_R[which(BTM_groups$SUBGROUP==BTM_groups_int[i]),c(1:4,7)]
  #ph_lim=ceiling(max(abs(df), na.rm = TRUE))
  ph_lim=2
  colorLS=colorRampPalette(colors = c("blue", "white", "red"))(n = 100)
  breakLS=seq(from = -ph_lim, to = ph_lim, length.out = 100)
  #Plot
  result=pheatmap::pheatmap(mat = df,
                  color = colorLS,
                  breaks = breakLS,
                  cluster_cols = TRUE,
                  cluster_rows = TRUE,
                  show_colnames = TRUE,
                  show_rownames = TRUE,
                  #annotation_col = df_anno,
                  #annotation_colors = adj_path_vt_colors,
                  fontsize = 10,
                  #cellheight = 10,
                  #cellwidth = 10,
                  filename=paste0('Residual_R_heatmap_',BTM_groups_int[i],'_genesymbol_MGmaxSBA_RNAseqRTSS_prevnarMAX_allgene_lim2.pdf'))
}

#Create scatterplots for 1 BTM across vaccines
BTM_plot='M85'
df=lapply(1:length(exp_FC), function(x) data.frame(FC=exprs(exp_FC[[x]])[which(BTM_groups$BTM==BTM_plot),],
                                                   ab_resid=exp_FC[[x]]$ab_resid,
                                                   vaccine=names(exp_FC)[x]))
df=Reduce(rbind, df)
#Plot
p=ggplot(df, aes(x=FC, y=ab_resid))+
  geom_point(shape=19)+
  xlab('log2FC')+
  ylab('Residual')+
  stat_cor()+
  theme_bw(base_size = 16)+
  ggtitle(BTM_groups$NAME[which(BTM_groups$BTM==BTM_plot)])+
  #theme(legend.position="none")+
  geom_smooth(method="lm", se=F, col='black')+
  facet_wrap(~vaccine, scales='free')
ggsave(paste0('BTM_AbResid_Scatterplot_perVaccine_',BTM_plot,'.pdf'), width=21, height=14)
