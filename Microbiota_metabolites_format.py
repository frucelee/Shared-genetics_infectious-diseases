#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000

###format microbiota datassets
cat microbiota_config | while read id
do
gunzip ${id}.gz
awk '$3<=1E-05' ${id} >123.subset.txt
Rscript --vanilla test.R 123.subset.txt
awk '{ print $1,$7}' 123.txt|sed '1d' >T2.txt
sed -i '1i SNP P' T2.txt
plink --bfile /workdir/shifang/scdrs/g1000_eur --clump T2.txt  --clump-r2 0.001 --out MetaGWAS_CAD_clumped_r0.05 --clump-kb 10000 --clump-p1 1E-05 --clump-p2 1 --threads 8
awk '{print $1,$3}' MetaGWAS_CAD_clumped_r0.05.clumped > SNP.valid
awk '{print $2}' SNP.valid >HO_clumped.txt
#head -n -2 HO_clumped.txt > new_file.txt
grep -w -F -f HO_clumped.txt 123.txt >/workdir/shifang/MR/microbiota/sig/${id}
rm -r 123.subset.txt HO_clumped.txt 123.txt T2.txt new_file.txt SNP.valid MetaGWAS_CAD_clumped_r0.05*
done

>>>R
##test.R
library(data.table)
library(tibble)
library(data.table)
library(magrittr)
library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
data0<-fread(exp_path,header=T)
colnames(data0)<-c("chromosome","base_pair_location","p_value","effect_allele","other_allele","beta","standard_error", "odds_ratio","ci_lower", "ci_upper","effect_allele_frequency", "variant_id")
ID<-fread("/workdir/shifang/scdrs/g1000_eur.bim",header=F)
data0$ID<-paste(data0$chromosome,data0$base_pair_location,sep="_")
ID$ID<-paste(ID$V1,ID$V4,sep="_")
bb<-intersect(ID$ID,data0$ID)
data0<-data0[data0$ID%in%bb,]
library(tibble)
ID<-data.frame(column_to_rownames(ID, var = "ID"))
ID$ID<-rownames(ID)
#rownames(ID)<-ID$ID
ID<-ID[data0$ID,]
data0$base_pair_location<-ID$V4
data0$SNP<-ID$V2
data0$chromosome<-ID$V1
# 2. 读取 1000 Genomes MAF 文件（SNP, A1, A2, MAF）
maf_table <- fread("/workdir/shifang/MR/merged.freq", header = TRUE)
# 3. 合并：按 SNP 对应
merged <- merge(data0, maf_table, by = "SNP")
# 4. 计算 EAF（effect_allele frequency）
# 假设 maf_table 中的 A1/A2 是两个真实的等位基因（未注明哪个是 minor）
merged$EAF <- mapply(function(effect_allele, a1, a2, maf) {
  # 判断哪个是 minor allele（假设字母靠前为 minor）
  alleles <- sort(c(a1, a2))
  minor <- alleles[1]
  major <- alleles[2]

  if (effect_allele == minor) {
    return(maf)
  } else if (effect_allele == major) {
    return(1 - maf)
  } else {
    return(NA)  # effect_allele 不在 A1/A2 中，返回 NA
  }
}, merged$effect_allele, merged$A1, merged$A2, merged$MAF)
data0<-merged
#data0$V2 <- sub('SNP', '', data0data0$V2)
#data0$V2 = paste0('rs', data0$V2)
#test<-fread("/scratch/users/s/h/shifang/QTL/pQTL_all_SNP_2807.txt")
#test=test[!duplicated(test$SNP),]
#data0 <- data0[data0$hm_rsid%in% test $SNP,]
#data0$hm_beta<-log(data0$odds_ratio)
#data0$standard_error=(data0$ci_upper-data0$ci_lower)/3.92
#data0$standard_error=sqrt(((data0$hm_beta)^2)/qchisq(data0$p_value,1,lower.tail=F))
bmi_exp_dat<-data0[,c("SNP","chromosome","base_pair_location","other_allele","effect_allele","beta","p_value","standard_error","EAF")]
head (bmi_exp_dat)
colnames(bmi_exp_dat)=c("SNP","Chromozes","pos","A2","A1","beta","pval","se","eaf")
#bmi_exp_dat<-subset(bmi_exp_dat,pval.outcome<=5e-08)
#bmi_exp_dat$pval.outcome<-10^((-1)*bmi_exp_dat$pval.outcome)
bmi_exp_dat$n<-c(rep("8956",dim(bmi_exp_dat)[1]))
bmi_exp_dat=bmi_exp_dat[!duplicated(bmi_exp_dat$SNP),]
bmi_exp_dat<-na.omit(bmi_exp_dat)
bmi_exp_dat$z<-bmi_exp_dat$beta/bmi_exp_dat$se
head (bmi_exp_dat)
colnames(bmi_exp_dat)=c("SNP","hg18","bq","A2","A1","beta","pval","se","eaf","n","z")
mhc <- bmi_exp_dat %>%
  dplyr::arrange(hg18, bq) %>%
  dplyr::filter(hg18 == 6) %>%
  dplyr::filter(bq >= 28477797 & bq <= 33448354)
mhcsnp <- mhc %>% dplyr::select(SNP)
# #
# # # remove MHC's snp from the full GWAS
bmi_exp_dat <- bmi_exp_dat %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))
write.table(bmi_exp_dat, '123.txt', quote = FALSE,row.names = FALSE)


###format metabolites datassets
cat metabolites_config | while read id
do
gunzip ${id}.gz
awk '$8<=5E-08' ${id} >123.subset.txt
Rscript --vanilla test.R 123.subset.txt
awk '{ print $1,$7}' 123.txt|sed '1d' >T2.txt
sed -i '1i SNP P' T2.txt
plink --bfile /workdir/shifang/scdrs/g1000_eur --clump T2.txt  --clump-r2 0.001 --out MetaGWAS_CAD_clumped_r0.05 --clump-kb 10000 --clump-p1 5E-08 --clump-p2 1 --threads 8
awk '{print $1,$3}' MetaGWAS_CAD_clumped_r0.05.clumped > SNP.valid
awk '{print $2}' SNP.valid >HO_clumped.txt
#head -n -2 HO_clumped.txt > new_file.txt
grep -w -F -f HO_clumped.txt 123.txt >/workdir/shifang/MR/metabolites/sig/${id}
rm -r 123.subset.txt HO_clumped.txt 123.txt T2.txt new_file.txt SNP.valid MetaGWAS_CAD_clumped_r0.05*
done

>>>R
##test.R
library(data.table)
library(tibble)
library(magrittr)
library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
data0<-fread(exp_path,header=T)
#colnames(data0)<-c("chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","effect_allele_frequency", "p_value", "variant_id","hm_coordinate_conversion", "hm_code", "rsid")
colnames(data0)<-c("chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","effect_allele_frequency", "p_value", "variant_id","rsid", "hm_code", "rsid1")
ID<-fread("/workdir/shifang/scdrs/g1000_eur.bim",header=F)
ID<-data.frame(column_to_rownames(ID, var = "V2"))
ID$V2<-rownames(ID)
bb<-intersect(data0$rsid,ID$V2)
data0<-data0[data0$rsid%in%bb,]
ID<-ID[data0$rsid,]
data0$base_pair_location<-ID$V4
data0$chromosome<-ID$V1
head (data0)
bmi_exp_dat<-data0[,c("rsid","chromosome","base_pair_location","other_allele","effect_allele","beta","p_value","standard_error","effect_allele_frequency")]
colnames(bmi_exp_dat)=c("SNP","Chromozes","pos","A2","A1","beta","pval","se","eaf")
bmi_exp_dat$n<-c(rep("8239",dim(bmi_exp_dat)[1]))
bmi_exp_dat$z<-bmi_exp_dat$beta/bmi_exp_dat$se
bmi_exp_dat=bmi_exp_dat[!duplicated(bmi_exp_dat$SNP),]
bmi_exp_dat$SNP[bmi_exp_dat$SNP==""] <- NA
bmi_exp_dat<-na.omit(bmi_exp_dat)
colnames(bmi_exp_dat)=c("SNP","hg18","bq","A2","A1","beta","pval","se","eaf","n","z")
mhc <- bmi_exp_dat %>%
  dplyr::arrange(hg18, bq) %>%
  dplyr::filter(hg18 == 6) %>%
  dplyr::filter(bq >= 28477797 & bq <= 33448354)
mhcsnp <- mhc %>% dplyr::select(SNP)
# #
# # # remove MHC's snp from the full GWAS
bmi_exp_dat <- bmi_exp_dat %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))
write.table(bmi_exp_dat, '123.txt', quote = FALSE,row.names = FALSE)

