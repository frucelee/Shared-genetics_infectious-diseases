#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

while read -r P1 P2 P3 P4 P5 P6 P7 P8 P9; do
    echo "Processing: $P1 $P2 $P3 $P4"

    OUT_PREFIX="${P1}_${P2}_${P3}"

    # 
    START=$((P5 - 500000))
    END=$((P5 + 500000))

    # 
    awk -v chr="$P4" -v start="$START" -v end="$END" '$2 == chr && $3 >= start && $3 <= end' \
        "/scratch/users/s/h/shifang/ldsc/data/used/$P2" > 123.txt

    awk -v chr="$P4" -v start="$START" -v end="$END" '$2 == chr && $3 >= start && $3 <= end' \
        "/scratch/users/s/h/shifang/ldsc/data/used/$P3" > 456.txt

    # fine-mapping
    Rscript --vanilla test.R2 123.txt 456.txt "$P6" "$P7" "$P8" "$P9"

    # 
    if [[ -f "123.csv" ]]; then
        mv 123.csv "${OUT_PREFIX}.csv"
    else
        echo "Warning: 123.csv not found for $OUT_PREFIX"
    fi

    # 
    rm -f 123.txt 456.txt

done <ccc

>>>R
#test.R2
library(data.table)
#library(magrittr) # needs to be run every time you start R and want to use %>%
#library(dplyr)
#library(tidyr)
#library("tibble")
#args <- commandArgs(trailingOnly=TRUE)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
data1<-fread(exp_path,header=T)
exp_path <- args[2]
data2<-fread(exp_path,header=T)
library(data.table)
N1 <- as.numeric(args[3])
N2<- as.numeric(args[4])
N3 <- as.numeric(args[5])
N4 <- as.numeric(args[6])
colnames(data1)=c("SNP","hg18","bq","A2","A1", "beta","pval", "se", "eaf", "n", "z" )
colnames(data2)=c("SNP","hg18","bq","A2","A1", "beta","pval", "se", "eaf", "n", "z" )
data1$varbeta <- (data1$SE)^2
data2$varbeta <- (data2$SE)^2
input <- merge(data1, data2, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
library("coloc")
#names(input )[names(input )=="V5"]="rs_id"
result <- coloc.abf(dataset1=list(pvalues=input$pval_eqtl, snp=input$SNP, type="cc", s=N2, N=N1,MAF=input$eaf_eqtl), dataset2=list(pvalues=input$pval_gwas,MAF=input$eaf_gwas, snp=input$SNP, type="cc", s=N4, N=N3))
#library(dplyr)
dd<-data.frame(t(data.frame(print(result[[1]]))))
write.csv(dd,"123.csv",quote=F,row.names=F)
