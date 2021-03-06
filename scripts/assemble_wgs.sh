#!/bin/bash -l

#SBATCH -A g2021012
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 16:00:00
#SBATCH -J wgs_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user matilda.maxwell.5390@student.uu.se

# Load modules
module load bioinfo-tools
module load soapdenovo

# Assembly
/sw/apps/bioinfo/SOAPdenovo/2.04-r240/rackham/bin/SOAPdenovo-63mer all -s config_file -K 49 -R -o graph_prefix 1>ass.log 2>ass.er
