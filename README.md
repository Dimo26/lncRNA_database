# Database design and programming for analysis of non-coding RNA in relation to disease causing genes.  

## Introduction
The purpose of this project is to create a database which combines information about lncRNA and the adjacent genes of which are related to diseases. Non-coding RNAs (ncRNAs) do not code proteins but have important regulatory functions in cells and tissues. NcRNAs can be involved in different disease processes (Nemeth et al, 2024). The pupose of this database is to collect all this information in one place for easier and more flexible access.

## Functions and applications in Bianca
The database has been designed to be compatable with Bianca UPPMAX. The database has various applications. It can be used for a double ended search where the user can find ncRNA of which regulates certain genes and vice versa. It is also possible to gain information about regulationtypes such as uppregulation or downregulation.

### Search types

The user can search gene to ncRNA and vice versa. The user can also search disease association. The results are presented in csv folder where all the related information from the soruces are combined giving valuable information for genes or for ncRNAs, lncRNAs combining all the sources populating the tables of the database. 


## Database python files
A general description of the codes provided in each file and the functionality. 
There is a wide potential with future applications of this code. 
We can add more sources and integrate these to gain more information or adjust the functionality of the code. 

Since the statistics of the information extracted from the different sources are abundant, we can modify and expand this database to extract more information. The downside of using Bianca cluster is that the files extracted from the sources might not potentially be updated to the latest version with the latest research finds as Bianca is an offline cluster. In order to prevent this and for regular updates, I have linked the original sources of which I have extracted my files for populating the database. This can help users to update these files according to the latest versions for more accurate and expanded results. Once the database is created and the data is fetched one can use search functions to search genes, disease or ncRNA.

### Database_setup.py
As the file name describes this python code contains the database setup and creation using sqlite. There are 6 tables (entities) of which the database schema is made of each with attributes of which will be populated using various sources. The tables describe the entities seen in image 1 describing the schema and its relations. 




<img width="1576" height="1108" alt="image" src="https://github.com/user-attachments/assets/6e2752c8-d847-4931-b20c-bc64d6359faa" />

Image 1. ER diagram of the database schema with the tables NC-RNA genes, Regulations, Regulation details, Genes, HPO terms, OMIM entries and Gene disease associations. (Created with ERDplus.com, 2023)

### config.yaml 
In the configuration file you can find all the sources that contribute to populating the database. Some sources include both fasta and gtf formats. However, the code uses GTF files only. Some files come from sources as zipped or text files of which are taken care of in the code itself to access the content. GenCode was used to extract information about long noncoding RNAs specifically annotations, lncipedia was used for human transcripts. Disease associations and phenotype information are from HPO and OMIM. Regulation and prediction sources include Noncode, Mirdb and Targetscan. The database can be improved with additional sources of which can be parsed in the data fetcher code. Since the database needs to be applied for clinical and reasearch use in Bianca all the soruces have been downloaded to a specific directory allowing the database to function properly offline. The urls for the different databases are included in the readme to provide information about the sources and versions of the files downloaded. 

#### Database sources 

General_genes_mart: 

Gene_info: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz   
  
Gene2refseq: https://ftp.ncbi.nlm.nih.gov/gene/DATA/Gene2refseq

Gencode gtf: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.long_noncoding_RNAs.gtf.gz

Gencode lncrna fasta: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.lncRNA_transcripts.fa.gz

Lncipedia fasta url: https://lncipedia.org/downloads/lncipedia_5_0/full-database/lncipedia_5_0.fasta
Lncipedia gtf_url: https://lncipedia.org/downloads/lncipedia_5_0/full-database/lncipedia_5_0_hg38.gtf

OMIM mart file found through directory downloaded from ensembl biomart: https://www.ensembl.org/biomart/martview/d32f188e084d91493eadd4a96965194d

HPO obo: https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo

Noncode fasta: http://www.noncode.org/datadownload/NONCODEv6_human.fa.gz
Noncode gtf: http://www.noncode.org/datadownload/NONCODEv6_human_hg38_lncRNA.gtf.gz

Mirdb: http://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz

Targetscan: https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip

EVLncRNAs: https://request.sdklab-biophysics-dzu.net/uploads/files/EVLncRNAs3_human_data.zip
LncRNAdb: File from supervior lnctard2.0.zip

### data_fetcher.py

The data fetcher uses the sources in the config.yaml and populates the databases. It fetches miRNA-target regulations fro the miRDB file. It also fetches the lncRA genes from gencode gtf files with quality score of above 0.7. LncRNA egulations are created with high quality from literature. HPO and omim terms are fetched from obo and mart file respectively. 

### search_functions.py
This code has NCRNADatabase which is a class thar contains the search function. It goes through the data from the sources and builds the query based on the search of the user ie if its a gene to ncRNA, ncRNA to gene or disease to find regulations, targets and associations. It esentially integrates and combines the tables of the database and the fetched data for search result. In this code the option to export the result in csv format is available to help the user extract the information to local machine. 

### parse_omim_data.py
It parses information from OMIM and combines it in the database schema for linking disease associations with the genes in the database.  

### gene_mapper.py
Due to the difference in each and every database where biomart, ensembl, mirdb all define the genes differently. I wanted the user to search genes using gene names rather than ensembl ID etc. In order to do so, the gene_mapper maps the names in ensmebl with the information about the corresponding ncRNA by using ncbi gene info. 

### run database.sh
In this file one can submit the searches. I have provided with examples that can be modified. 


## References 
Nemeth, K., Bayraktar, R., Ferracin, M. et al. 2024. Non-coding RNAs in disease: from mechanisms to therapeutics. Nat Rev Genet 25, 211–232.

Kinsella, R. J., Kähäri, A., Haider, S., Zamora, J., Proctor, G., Spudich, G., Almeida-King, J., Staines, D., Derwent, P., Kerhornou, A., Kersey, P., & Flicek, P. 2011. Ensembl BioMarts: a hub for data retrieval across taxonomic space. Database: the journal of biological databases and curation, 2011, bar030.

Garth R. Brown, Vichet Hem, Kenneth S. Katz, Michael Ovetsky, Craig Wallin, Olga Ermolaeva, Igor Tolstoy, Tatiana Tatusova, Kim D. Pruitt, Donna R. Maglott, Terence D. Murphy. 2015. Gene: a gene-centered information resource at NCBI, Nucleic Acids Research. Volume 43. Issue D1. p. D36–D42.

Frankish A, Diekhans M, Ferreira AM, Johnson R, Jungreis I, Loveland J, Mudge JM, Sisu C, Wright J, Armstrong J, Barnes I, Berry A, Bignell A, Carbonell Sala S, Chrast J, Cunningham F, Di Domenico T, Donaldson S, Fiddes IT, García Girón C, Gonzalez JM, Grego T, Hardy M, Hourlier T, Hunt T, Izuogu OG, Lagarde J, Martin FJ, Martínez L, Mohanan S, Muir P, Navarro FCP, Parker A, Pei B, Pozo F, Ruffier M, Schmitt BM, Stapleton E, Suner MM, Sycheva I, Uszczynska-Ratajczak B, Xu J, Yates A, Zerbino D, Zhang Y, Aken B, Choudhary JS, Gerstein M, Guigó R, Hubbard TJP, Kellis M, Paten B, Reymond A, Tress ML, Flicek P. 2019. GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Research: 47(D1): D766-D773.

Volders PJ, Verheggen K, Menschaert G, Vandepoele K, Martens L, Vandesompele J, Mestdagh P. 2015. An update on LNCipedia: a database for annotated human lncRNA sequences. Nucleic Acids research 43: D174-D180. 

Hamosh, A., Scott, A. F., Amberger, J. S., Bocchini, C. A., & McKusick, V. A. 2005. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. Nucleic acids research, 33 (Database issue), D514–D517. 

Robinson, P. N., Köhler, S., Bauer, S., Seelow, D., Horn, D., & Mundlos, S. 2008. The Human Phenotype Ontology: a tool for annotating and analyzing human hereditary disease. American journal of human genetics, 83(5), 610–615.
  
Liu, C., Bai, B., Skogerbø, G., Cai, L., Deng, W., Zhang, Y., Bu, D., Zhao, Y., & Chen, R. 2005. NONCODE: an integrated knowledge database of non-coding RNAs. Nucleic acids research, 33(Database issue), D112–D115. 

Yuhao Chen, Xiaowei Wang. 2020 miRDB: an online database for prediction of functional microRNA targets, Nucleic Acids Research: 48, D1:D127–D131.

McGeary, S. E., Lin, K. S., Shi, C. Y., Pham, T. M., Bisaria, N., Kelley, G. M., & Bartel, D. P. 2019. The biochemical basis of microRNA targeting efficacy. Science (New York, N.Y.), 366(6472).

Zhou, B., Ji, B., Shen, C., Zhang, X., Yu, X., Huang, P., Yu, R., Zhang, H., Dou, X., Chen, Q., Zeng, Q., Wang, X., Cao, Z., Hu, G., Xu, S., Zhao, H., Yang, Y., Zhou, Y., & Wang, J. (2024). EVLncRNAs 3.0: an updated comprehensive database for manually curated functional long non-coding RNAs validated by low-throughput experiments. Nucleic acids research, 52(D1), D98–D106. 

Zhao, H., Yin, X., Xu, H., Liu, K., Liu, W., Wang, L., Zhang, C., Bo, L., Lan, X., Lin, S., Feng, K., Ning, S., Zhang, Y., & Wang, L. 2023. LncTarD 2.0: an updated comprehensive database for experimentally-supported functional lncRNA-target regulations in human diseases. Nucleic acids research, 51(D1), D199–D207. 


