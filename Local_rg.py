#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

##GNOVA analysis
cat /scratch/users/s/h/shifang/ldsc/data/xaa | while read PARAM
do
  # PARAM contains the expression values for this gene
  P1=$(echo "$PARAM" | awk '{ print $2 }')
  P2=$(echo "$PARAM" | awk '{ print $3 }')
  P3=$(echo "$PARAM" | awk '{ print $5 }')
  ID=$(echo "$PARAM" | awk '{ print $1 }')

  # 
  BASENAME=$(basename "$P1" .sumstats.gz)  # 

  ./gnova.py \
    /scratch/users/s/h/shifang/ldsc/data/format/${ID}.txt.sumstats.gz \
    /scratch/users/s/h/shifang/ldsc/data/format/$P1.sumstats.gz \
    --N1 $P2 \
    --N2 $P3 \
    --bfile /scratch/users/s/h/shifang/ldsc/MAGMA/1000G_maf0.05/1000G_maf0.05_chr@ \
    --out /scratch/users/s/h/shifang/ldsc/GNOVA-master/${ID}_${BASENAME}.txt
done

##Local Rg analysis
input_file="/scratch/users/s/h/shifang/ldsc/data/xaa"
bfile_prefix="/scratch/users/s/h/shifang/ldsc/SUPERGNOVA/example/eur_chr@_SNPmaf5"
partition_prefix="/scratch/users/s/h/shifang/ldsc/SUPERGNOVA/partition/eur_chr@.bed"
output_dir="/scratch/users/s/h/shifang/ldsc/SUPERGNOVA/results"

while read -r line; do
    ID=$(echo "$line" | awk '{print $1}' | tr -d '\r\n')
    P1=$(echo "$line" | awk '{print $2}' | tr -d '\r\n')
    P2=$(echo "$line" | awk '{print $3}' | tr -d '\r\n')
    P3=$(echo "$line" | awk '{print $5}' | tr -d '\r\n')

    sumstats_file="/scratch/users/s/h/shifang/ldsc/data/format/${ID}.txt.sumstats.gz"
    ref_file="/scratch/users/s/h/shifang/ldsc/data/format/${P1}.sumstats.gz"
    out_prefix="${output_dir}/${ID}_${P1}"  # 

    echo "Running SUPERGNOVA: ID=$ID, Trait=$P1, N1=$P2, N2=$P3"

    python3 supergnova.py "$sumstats_file" "$ref_file" \
        --N1 "$P2" \
        --N2 "$P3" \
        --thread 8 \
        --bfile "$bfile_prefix" \
        --partition "$partition_prefix" \
        --out "$out_prefix"
done < "$input_file"
