##MR analysis
#step 1: clump the datasets
cat config | while read id
do
Rscript --vanilla test.R1 ${id}
plink --bfile /scratch/users/s/h/shifang/ldsc/MAGMA/g1000_eur --clump test.txt0  --clump-r2 0.001 --out MetaGWAS_CAD_clumped_r0.05 --clump-kb 10000 --clump-p1 1 --clump-p2 1 --threads 8
awk '{print $1,$3}' MetaGWAS_CAD_clumped_r0.05.clumped > SNP.valid
awk '{print $2}' SNP.valid >HO_clumped.txt
rm -r test.txt0 *.log
rm -r MetaGWAS_CAD_clumped_r0.05.clumped
rm -r SNP.valid
Rscript --vanilla test.R2 ${id}
rm -r test.txt0 HO_clumped.txt
done

>>>R
#test.R1
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
#exposure
data0<-fread(exp_path,header=T)
data0<-subset(data0,pval<=5e-08)
test.txt0<-data0[,c("SNP","pval")]
colnames(test.txt0)<-c("SNP","P")
write.table(test.txt0, 'test.txt0', quote = FALSE,row.names = FALSE) 
#test.R2
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
#exposure
data0<-fread(exp_path,header=T)
data1<-fread("HO_clumped.txt",header=T)
data0<-data0[data0$SNP%in%data1$SNP,]
write.csv(data0,paste(exp_path,"csv",sep = "_clump."), quote = FALSE)