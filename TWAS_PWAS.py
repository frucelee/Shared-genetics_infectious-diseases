#!/bin/bash
#
#SBATCH --job-name=shifang
#SBATCH --output=45_mpi.txt
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=40000

#PWAS analysis
cat PWAS_ID | while read id
do
cp /scratch/users/s/h/shifang/ldsc/data/format/${id}.sumstats.gz /scratch/users/s/h/shifang/TWAS/${id}.sumstats.gz
gunzip ${id}.sumstats.gz
awk 'NF==5 {print $0 "\tNA"}' ${id}.sumstats | awk 'NR==1{$6="CHISQ"} 1' > ${id}_output.txt
for i in {1..23}; do
#until [ $CHR -lt 1 ]
#do
Rscript /scratch/users/s/h/shifang/TWAS/scripts/PWAS.assoc_test.R \
--sumstats ${id}_output.txt \
--weights /scratch/users/s/h/shifang/TWAS/PWAS_EA/Plasma_Protein_EA_hg19.pos \
--weights_dir /scratch/users/s/h/shifang/TWAS/PWAS_EA/Plasma_Protein_weights_EA/ \
--ref_ld_chr /scratch/users/s/h/shifang/TWAS/LDref/EUR/chr \
--force_model enet \
--chr ${i} \
--out /scratch/users/s/h/shifang/TWAS/PWAS_results/${id}_chr${i}.out
#
done
rm -r ${id}_output.txt ${id}.sumstats
done

##TWAS analysis
cat /scratch/users/s/h/shifang/TWAS/config | while read id
do
cp /scratch/users/s/h/shifang/ldsc/data/format/${id}.sumstats.gz /scratch/users/s/h/shifang/TWAS/TWAS_model/${id}.sumstats.gz
gunzip ${id}.sumstats.gz
awk 'NF==5 {print $0 "\tNA"}' ${id}.sumstats | awk 'NR==1{$6="CHISQ"} 1' > ${id}
done

ls *gz|awk '{print substr($0,1,length($0)-7)}' | sort | uniq > TWAS_config

cat TWAS_config | while read id
do
tar -xzvf ${id}.tar.gz
#awk 'NF==5 {print $0 "\tNA"}' ${id}.sumstats | awk 'NR==1{$6="CHISQ"} 1' > ${id}_output.txt
for i in {1..23}; do
Rscript /scratch/users/s/h/shifang/TWAS/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /scratch/users/s/h/shifang/TWAS/TWAS_model/finngen_R12_AB1_INTESTINAL_INFECTIONS.txt \
--weights /scratch/users/s/h/shifang/TWAS/TWAS_model/model/GTExv8.EUR.${id}.pos \
--weights_dir /scratch/users/s/h/shifang/TWAS/TWAS_model/model \
--ref_ld_chr /scratch/users/s/h/shifang/TWAS/LDref/EUR/chr \
--chr ${i} \
--out /scratch/users/s/h/shifang/TWAS/TWAS_model/results/GCST90475667.h.tsv_${id}_${i}.dat
done
rm -r ${id}*
done