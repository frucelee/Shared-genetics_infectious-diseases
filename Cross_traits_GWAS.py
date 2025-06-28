#!/bin/bash
#
#SBATCH --output=10_mpi.txt
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000

##MTAG analysis
cat /scratch/users/s/h/shifang/ldsc/mtag/123.txt | while read PARAM
do
# PARAM contains the expression values for this gene
P1=`echo $PARAM | awk '{ print $1 }'`
P2=`echo $PARAM | awk '{ print $2 }'`
outfile="/scratch/users/s/h/shifang/ldsc/mtag/data/result/${P1}_${P2}"
python /scratch/users/s/h/shifang/ldsc/mtag/mtag.py --sumstats /scratch/users/s/h/shifang/ldsc/mtag/data/$P1.txt,/scratch/users/s/h/shifang/ldsc/mtag/data/$P2 --out ${outfile} --n_min 0.0 --stream_stdout --no_overlap
python /scratch/users/s/h/shifang/ldsc/mtag/mtag.py --sumstats /scratch/users/s/h/shifang/ldsc/mtag/data/$P1.txt,/scratch/users/s/h/shifang/ldsc/mtag/data/$P2 --out ${outfile} --stream_stdout --fdr --stream_stdout
done


##prepare the data for CPASSOC analysis
while read -r PARAM; do
  # 
  P1=$(echo "$PARAM" | awk '{ print $1 }')
  P2=$(echo "$PARAM" | awk '{ print $2 }')

  echo "Running Rscript for: $P1 $P2"

  # 
  Rscript --vanilla test.R10 "$P1.txt" "$P2"

  # 123.txt
  if [[ -f 123.txt ]]; then
    mv 123.txt "${P1}_${P2}.txt"
  else
    echo "Warning: 123.txt not found after running for $P1 $P2"
  fi

done < /scratch/users/s/h/shifang/ldsc/mtag/123.txt

>>>R
##test.R10
library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
data1<-fread(exp_path,header=T)
library(data.table)
exp_path1 <- args[2]
data2<-fread(exp_path1,header=T)
data1<-data1[-1,c(1,8)]
data2<-data2[-1,c(1,8)]
data3<-merge(data1,data2,by="snpid",all=T)
head (data3)
data3<-na.omit(data3)
head (data3)
write.table(data3, '123.txt', quote = FALSE,row.names = FALSE,sep="\t")

### CPASSOC with cross-trait statistical heterogeneity (SHet) analysis
while read -r line; do
  # 
  P1=$(echo "$line" | awk '{print $1}')
  P2=$(echo "$line" | awk '{print $2}')
  P3=$(echo "$line" | awk '{print $3}')
  P4=$(echo "$line" | awk '{print $5}')
  # 
  fname="${P1}_${P2}.txt"

  echo "Running Rscript with: $fname $P1 $P2"

  #
  Rscript --vanilla test.R21 "$fname" "$P3" "$P4"

  # 
  if [[ -f 456.txt ]]; then
    mv 456.txt "CAPSSOC_${P1}_${P2}.txt"
  else
    echo "Warning: 456.txt not found for ${P1} ${P2}"
  fi

done < /scratch/users/s/h/shifang/ldsc/data/xaa


