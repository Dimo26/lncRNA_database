#!/bin/bash -l

#SBATCH -A sens2023005
#SBATCH -M bianca
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:20:00
#SBATCH -J Run_database 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user dima.mohsin.1749@student.uu.se
#SBATCH --mem=16G
 
# Load modules
module load bioinfo-tools
module load python/3.9.5
module load sqlite/3.34.0
#pip install --user xlrd==2.0.1 openpyxl
#Run scripts

#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/database_setup.py
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/data_fetcher.py
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/parse_omim.py
#python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/gene_mapper.py
 

python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type gene-to-mirna --gene TP53
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type gene-to-lncrna --gene BRCA1
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type mirna-to-gene --mirna hsa-miR-21 
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type lncrna-to-gene --lncrna HOTAIR --csv
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type lncrna-to-gene --lncrna MALAT1 --evlncrnas-threshold 0.8
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type gene-to-ncrna --gene EGFR
python /proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/search_functions.py --search-type disease --disease "breast cancer"

