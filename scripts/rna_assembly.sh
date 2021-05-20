#!/bin/bash -l

#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -J align_rna
#SBATCH --mail-type=ALL
#SBATCH --mail-user matilda.maxwell.5390@student.uu.se


# Load modules
module load bioinfo-tools bowtie samtools trinity
# commands

Trinity --seqType fq --max_memory 50G --left ~/genome_analysis/genome_analysis/out_files/trimmomatic_rna_seq/sel3_SRR1719266_1P.fq.gz,~/genome_analysis/genome_analysis/out_files/trimmomatic_rna_seq/sel3_SRR1719266_1U.fq.gz --right ~/genome_analysis/genome_analysis/out_files/trimmomatic_rna_seq/sel3_SRR1719266_2P.fq.gz,~/genome_analysis/genome_analysis/out_files/trimmomatic_rna_seq/sel3_SRR1719266_2U.fq.gz --CPU 2 --output /home/mama5390/genome_analysis/genome_analysis/out_files/trinity_rna_seq


