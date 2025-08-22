#!/bin/bash
#
#SBATCH --output=66_mpi.txt
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000

cat ID | while read id
do
##for LDSC analysis
awk '{ print $1, $10, $11, $5, $4, $7}' ${id}>bbb_1
sed '1s/\<n\>/N/; 1s/\<z\>/Z/; 1s/\<pval\>/P/' bbb_1 > bbb_10
munge_sumstats.py --sumstats bbb_10 --out /scratch/users/s/h/shifang/ldsc/data/format/${id} --merge-alleles /scratch/users/s/h/shifang/ldsc/MAGMA/w_hm3.snplist
rm -r bbb_1 bbb_10
##for MTAG analysis
awk '{ print $1, $2, $3, $5, $4, $9,$11,$7,$10}' ${id} | sed 1d>/scratch/users/s/h/shifang/ldsc/mtag/data/${id}
sed  -i '1i snpid chr bpos a1 a2 freq z pval n' /scratch/users/s/h/shifang/ldsc/mtag/data/${id}
done
