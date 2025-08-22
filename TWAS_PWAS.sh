#!/bin/bash
#
#SBATCH --job-name=shifang
#SBATCH --output=45_mpi.txt
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=40000

#Build the model for plasma proteins
cat PWAS_ID_build | while read PARAM
do
# PARAM contains the expression values for this gene
GNAME=`echo $PARAM | awk '{ print $5 }'`
CHR=`echo $PARAM | awk '{ print $2 }'`
P0=`echo $PARAM | awk '{ p=$3 - 500e3; if(p<0) p=0; print p; }'`
P1=`echo $PARAM | awk '{ print $3 + 500e3 }'`
ID=`echo $PARAM | awk '{ print $1 }'`

# Convert the expression matrix to a PLINK format phenotype
plink --silent --bfile /globalscratch/users/s/h/shifang/eQTL/monkey/muscle_0yr --chr $CHR --from-bp $P0 --to-bp $P1 --make-bed --out OUT --pheno gene.txt --mpheno $ID --maf 0.0001 --allow-no-sex
# Run FUSION
Rscript FUSION.compute_weights.R --tmp /globalscratch/users/s/h/shifang/fusion/results/$GNAME.tmp --bfile OUT --PATH_plink /globalscratch/users/s/h/shifang/fusion/plink1.9/plink --out /globalscratch/users/s/h/shifang/fusion/results/$GNAME --verbose 0 --save_hsq  --models lasso,top1,enet --covar cov.txt --hsq_p 1.0 
rm OUT.*
done

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

###only keep the risk genes
> keep.tmp
while read base; do
  echo "${base}.hsq" >> keep.tmp
  echo "${base}.wgt.RDat" >> keep.tmp
done </scratch/users/s/h/shifang/TWAS/TWAS_model/model/raw/ID

##Begin TWAS analysis
cat TWAS_config | while read id
do
tar -xzvf ${id}.tar.gz

grep -w -F -f /scratch/users/s/h/shifang/TWAS/TWAS_model/model/raw/ID ${id}.nofilter.pos |awk '{ print $3}' >${id}_123

# 
> ${id}_keep.tmp
while read base; do
  echo "${base}.hsq" >> ${id}_keep.tmp
  echo "${base}.wgt.RDat" >> ${id}_keep.tmp
done < ${id}_123   # 

# 
#for f in *; do
#  if ! grep -qx "$f" ${id}_keep.tmp; then
#    echo "Deleting $f"
#    rm "$f"
#  fi
#done

# 
src_dir="${id}"
dst_dir="${id}_selected"
mkdir -p "$dst_dir"
# 
while read base; do
  for ext in ".hsq" ".wgt.RDat"; do
    file="${base}${ext}"
    if [ -f "$src_dir/$file" ]; then
      mv "$src_dir/$file" "$dst_dir/"
      echo "Moved $file"
    else
      echo "Warning: $file not found in $src_dir"
    fi
  done
done <"${id}_123"
rm -r ${id}
mv ${id}_selected ${id}
rm ${id}_keep.tmp

#awk 'NF==5 {print $0 "\tNA"}' ${id}.sumstats | awk 'NR==1{$6="CHISQ"} 1' > ${id}_output.txt
grep -w -F -f ${id}_123 ${id}.pos > ${id}.pos_tmp
mv ${id}.pos_tmp ${id}.pos
sed  -i '1i PANEL WGT ID CHR P0 P1 N' ${id}.pos

for i in {1..23}; do
Rscript /scratch/users/s/h/shifang/TWAS/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /scratch/users/s/h/shifang/TWAS/TWAS_model/finngen_R12_AB1_INTESTINAL_INFECTIONS.txt \
--weights /scratch/users/s/h/shifang/TWAS/TWAS_model/model/${id}.pos \
--weights_dir /scratch/users/s/h/shifang/TWAS/TWAS_model/model \
--ref_ld_chr /scratch/users/s/h/shifang/TWAS/LDREF/1000G.EUR. \
--chr ${i} \
--out /scratch/users/s/h/shifang/TWAS/TWAS_model/results/finngen_R12_AB1_INTESTINAL_INFECTIONS.txt_${id}_${i}.dat

done
rm -r ${id}* ${id}_123 
awk 'FNR==1 && NR!=1 {next} 1' /scratch/users/s/h/shifang/TWAS/TWAS_model/results/finngen_R12_AB1_INTESTINAL_INFECTIONS.txt_${id}* > /scratch/users/s/h/shifang/TWAS/TWAS_model/results/all/finngen_R12_AB1_INTESTINAL_INFECTIONS.txt_${id}.txt
done
