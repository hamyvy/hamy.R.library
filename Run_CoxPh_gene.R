library(dplyr)
library(data.table)
library("survival")
options(width=150)
library("optparse")


option_list <- list(
  make_option(c("-g", "--chromosome"), type="character", default=NULL, help="load score file name", metavar="character"),
  make_option(c("-f", "--scoreFolder"), type="character", default=NULL, help="directory where results files should be stored", metavar="character")
 )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


source("/hamy/wes/load/flib.R")

WORKING_DIR = '/hamy/wes/load/'
PHENO_DIR = '/hamy/phenotype/'

#------ load pheno -------------------------------------------------------------------------------------

dr <- fread(paste0(PHENO_DIR,'death_record_participant.csv'),sep=",",data.table=F)
#dc <- fread(paste0(PHENO_DIR,'death_record_death_cause.csv'),sep=",",data.table=F)
pheno <- fread(paste0(WORKING_DIR,'UKB_extracted_pheno.txt'),data.table=F)
pheno <- pheno[,c(1:5,9:21,52:57)]

pheno <- left_join(pheno,dr[,1:2],by=c("feid"="eid"))
pheno$livetime <- pheno$livetime +2
pheno[!is.na(pheno$p40007_i0),]$death <- 1
pheno[!is.na(pheno$p40007_i0),]$livetime <- pheno[!is.na(pheno$p40007_i0),]$p40007_i0
pheno$p40007_i0 <- NULL

pheno$centre <- factor(pheno$centre)
pheno$batch <- factor(pheno$batch)
pheno$sex <- factor(pheno$sex)

anc <- fread("UKBB_pan_ancestry.txt",data.table=F)
pheno <- inner_join(pheno,anc[,2:4],by=c('feid'='IID'))

pheno <- pheno[!is.na(pheno$sex),]
pheno <- pheno[!is.na(pheno$unrelated_anc),]
#-------------------------------------------------------------------------------------------------------
# CoxPH

PCs = paste0('PC',c(1:10))
covs = paste(c('age_at_baseline','sex','batch',PCs),collapse='+')


gene_nlof_snp <- fread(paste0(WORKING_DIR,'load_score/',opt$scoreFolder,'/snp_c',opt$chromosome,'.txt'),data.table=F)
gene_nlof_indel <- fread(paste0(WORKING_DIR,'load_score/',opt$scoreFolder,'/indel_c',opt$chromosome,'.txt'),data.table=F)
gene_nlof <- data.frame(sum_snp_indel(gene_nlof_snp,gene_nlof_indel))


res1 <- CoxPH(gene_nlof_snp,colnames(gene_nlof_snp)[2:ncol(gene_nlof_snp)],covs,pheno[pheno$unrelated_anc=="EUR",])
res1$varset <- 'snp'
res2 <- CoxPH(gene_nlof_indel,colnames(gene_nlof_indel)[2:ncol(gene_nlof_indel)],covs,pheno[pheno$unrelated_anc=="EUR",])
res2$varset <- 'indel'
res3 <- CoxPH(gene_nlof,colnames(gene_nlof)[2:ncol(gene_nlof)],covs,pheno[pheno$unrelated_anc=="EUR",])
res3$varset <- 'all'
res <- rbind(res1,res2,res3)
res$chrom <- opt$chromosome


write.table(res,paste0('CoxPH_EUR_bygene_',opt$scoreFolder,'_',opt$chromosome,'.txt'),quote=F,row.names=F,sep='\t')

