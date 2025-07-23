#!/bin/bash -l

#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:20:00
#SBATCH -J Run_database 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user dima.mohsin.1749@student.uu.se
#SBATCH --mem=16G
 
# Load modules
module load bioinfo-tools
module load python/3.9.5
module load sqlite/3.34.0

#Run scripts

#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/database_setup.py
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/data_fetcher.py
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/parse_omim.py
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/gene_mapper.py

#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type gene-to-ncrna --gene BRCA1 --export csv
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type ncrna-to-gene --ncrna hsa-miR-21 --export csv
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type ncrna-to-gene --ncrna hsa-miR-125a-3p --export csv
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type disease --disease "breast cancer"
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type stats --export csv

