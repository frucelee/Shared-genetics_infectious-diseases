#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

cat ID | while read PARAM; do
    # 
    P1=$(echo "$PARAM" | awk '{ print $1 }')
    P2=$(echo "$PARAM" | awk '{ print $3 }')
    P3=$(echo "$PARAM" | awk '{ print $5 }')
    P4=$(echo "$PARAM" | awk '{ print $6 }')
    P5=$(echo "$PARAM" | awk '{ print $7 }')

    echo "Processing: P1=$P1, P2=$P2, P3=$P3, P4=$P4, N=$P5"

    OUT_PREFIX="${P1}_${P2}_${P3}_${P4}"

    # 
    awk -v p2="$P2" -v start="$P3" -v end="$P4" '$2 == p2 && $3 >= start && $3 <= end { print $1, $11 }' \
        /scratch/users/s/h/shifang/ldsc/data/used/$P1 > 456.txt

    # 
    awk '{print $1}' 456.txt > SNP.txt

    # LD
    plink -bfile /scratch/users/s/h/shifang/ldsc/MAGMA/g1000_eur \
          --keep-allele-order \
          --r square \
          --extract SNP.txt \
          --out sig_locus_mt

    # fine-mapping
    Rscript --vanilla susie.R "$P5"

    # 
    sed -i "s/shifang/${OUT_PREFIX}/g" PIP_pip.csv

    # 
    cp PIP_pip.csv /scratch/users/s/h/shifang/ldsc/data/finemapping/${OUT_PREFIX}_pip.csv

    # 
    rm -f PIP_pip.csv 456.txt SNP.txt sig_locus_mt*
done

>>>R
##susie.R
library(susieR)
library(data.table)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library(tidyr)
library("tibble")
args <- commandArgs(trailingOnly=TRUE)
#exp_path <- args[1]
#N<-exp_path
N <- as.numeric(args[1])
#exposure
#args1 <- commandArgs(trailingOnly=TRUE)
#exp_path1 <- args1[2]
#exposure
ld = as.matrix(fread("sig_locus_mt.ld"))
st = read.table("456.txt", header=F)
rownames(ld)<-st$V1
colnames(ld)<-st$V1
#z_scores = st$V6/st$V7
#fitted_rss = susie_rss(z_scores, R=ld, L=10, refine=F)
fitted_rss = susie_rss(z = st$V2, R=ld, L=10, refine=F,n=N)
Susiedf <- data.frame("pip" = fitted_rss$pip,  "CS" = NA)%>%
    rownames_to_column("SNP")%>%
    separate(SNP, into =c("CHR","BP","REF","ALT","rsID"),remove = F) ##
  for (x in names(fitted_rss$sets$cs)){
    Susiedf$CS[fitted_rss$sets$cs[[x]]] <- x
  }
Susiedf<-subset(Susiedf,pip>0)
Susiedf$IDD<-c(rep("shifang",dim(Susiedf)[1]))
write.csv(Susiedf,"PIP_pip.csv",quote=F,row.names=F)
