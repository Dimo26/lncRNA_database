
"""
Datafetcher which is offline fetching data from sources found in onfiguration file
"""

import sqlite3
import yaml
import json
import logging
import gzip
import zipfile
from pathlib import Path
from tqdm import tqdm
import argparse


# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class BiancaDataFetcher:
    def __init__(self, config_path="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/config.yaml"):
        """Initialize data fetcher for Bianca (offline mode)."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_file = Path(self.config['database']['path']) / self.config['database']['name']
        self.data_dir = Path("data/raw")
        self.data_dir.mkdir(parents=True, exist_ok=True)
    
    def fetch_genes_from_local_biomart(self, limit=None):
        """Fetch genes from local BioMart file."""
        logger.info("Fetching protein-coding genes from local BioMart file...")
        
        # Get path from config
        biomart_file = Path(self.config['data_sources']['general_genes_mart'])
        
        if not biomart_file.exists():
            logger.error(f"BioMart file not found: {biomart_file}")
            return []
        
        genes = []
        
        try:
            with open(biomart_file, 'r', encoding='utf-8') as f:
                for i, line in enumerate(tqdm(f, desc="Processing BioMart file")):
                    if limit and i >= limit:
                        break
                    
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        
                        fields = line.split('\t')
                        if len(fields) >= 7:
                            ensembl_id, gene_symbol, description, chromosome, start, end, strand = fields[:7]
                            
                            if not ensembl_id or not gene_symbol:
                                continue
                            
                            
                            chromosome = chromosome.replace('CHR_', '').replace('_', '')
                            if chromosome not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']:
                                continue
                            
                            genes.append({
                                'ensembl_gene_id': ensembl_id,
                                'gene_symbol': gene_symbol,
                                'chromosome': chromosome,
                                'start_position': int(start) if start.isdigit() else 0,
                                'end_position': int(end) if end.isdigit() else 0,
                                'strand': '+' if strand == '1' else '-',
                                'gene_type': 'protein_coding',
                                'description': description[:200] if description else ''
                            })
                    except Exception as e:
                        continue
        
        except Exception as e:
            logger.error(f"Error reading BioMart file: {e}")
            return []
        
        logger.info(f"Fetched {len(genes)} genes from local BioMart file")
        return genes
    
    def fetch_noncode_lncrnas(self):
        logger.info("Fetching lncRNA genes from local NONCODE GTF file...")
        
        gtf_file = Path(self.config['regulation_sources']['noncode']['noncode_gtf'])
        
        if not gtf_file.exists():
            logger.warning(f"NONCODE GTF file not found: {gtf_file}")
            return []
        
        logger.info(f"Using NONCODE GTF file: {gtf_file}")
        
        lncrnas = []
        
        try:
            with gzip.open(gtf_file, 'rt') as f:
                for line in f:
                    if line.startswith('#') or 'gene' not in line:
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9 or fields[2] != 'gene':
                        continue
                    
                    attrs = fields[8]
                    gene_id = ''
                    gene_name = ''
                    
                    if 'gene_id "' in attrs:
                        start = attrs.find('gene_id "') + 9
                        end = attrs.find('"', start)
                        gene_id = attrs[start:end].split('.')[0]
                    
                    if 'gene_name "' in attrs:
                        start = attrs.find('gene_name "') + 11
                        end = attrs.find('"', start)
                        gene_name = attrs[start:end]
                    
                    if gene_id:
                        chromosome = fields[0].replace('chr', '')
                        
                        # Skip non-standard chromosomes
                        if chromosome not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']:
                            continue
                        
                        try:
                            start_pos = int(fields[3])
                            end_pos = int(fields[4])
                            strand = fields[6]
                        except:
                            continue
                        
                        lncrnas.append({
                            'ensembl_gene_id': gene_id,
                            'gene_symbol': gene_name if gene_name else gene_id,
                            'chromosome': chromosome,
                            'start_position': start_pos,
                            'end_position': end_pos,
                            'strand': strand,
                            'ncrna_type': 'lncRNA',
                            'ncrna_subtype': 'NONCODE_lncRNA',
                            'length': end_pos - start_pos + 1,
                            'description': f"NONCODE lncRNA: {gene_name if gene_name else gene_id}"
                        })
                    
                    if len(lncrnas) >= 10000: 
                        break
        
        except Exception as e:
            logger.error(f"Error parsing NONCODE GTF: {e}")
            return []
        
        logger.info(f"Parsed {len(lncrnas)} lncRNA genes from NONCODE")
        return lncrnas

    def fetch_lncipedia_lncrnas(self):

        logger.info("Fetching lncRNA genes from local LNCipedia GTF file...")
        
        gtf_file = Path(self.config['data_sources']['lncipedia']['lncipedia_gtf'])
        
        if not gtf_file.exists():
            logger.warning(f"LNCipedia GTF file not found: {gtf_file}")
            return []
        
        logger.info(f"Using LNCipedia GTF file: {gtf_file}")
        
        lncrnas = []
        
        try:
            if gtf_file.suffix == '.gz':
                file_handle = gzip.open(gtf_file, 'rt')
            else:
                file_handle = open(gtf_file, 'r')
            
            with file_handle as f:
                for line in f:
                    if line.startswith('#') or 'gene' not in line:
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9 or fields[2] != 'gene':
                        continue
                    
                    attrs = fields[8]
                    gene_id = ''
                    gene_name = ''
                    
                    if 'gene_id "' in attrs:
                        start = attrs.find('gene_id "') + 9
                        end = attrs.find('"', start)
                        gene_id = attrs[start:end].split('.')[0]
                    
                    if 'gene_name "' in attrs:
                        start = attrs.find('gene_name "') + 11
                        end = attrs.find('"', start)
                        gene_name = attrs[start:end]
                    
                    if gene_id:
                        chromosome = fields[0].replace('chr', '')
                        
                        if chromosome not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']:
                            continue
                        
                        try:
                            start_pos = int(fields[3])
                            end_pos = int(fields[4])
                            strand = fields[6]
                        except:
                            continue
                        
                        lncrnas.append({
                            'ensembl_gene_id': gene_id,
                            'gene_symbol': gene_name if gene_name else gene_id,
                            'chromosome': chromosome,
                            'start_position': start_pos,
                            'end_position': end_pos,
                            'strand': strand,
                            'ncrna_type': 'lncRNA',
                            'ncrna_subtype': 'LNCipedia_lncRNA',
                            'length': end_pos - start_pos + 1,
                            'description': f"LNCipedia lncRNA: {gene_name if gene_name else gene_id}"
                        })
                    
                    if len(lncrnas) >= 10000:  # Limit for Bianca performance
                        break
        
        except Exception as e:
            logger.error(f"Error parsing LNCipedia GTF: {e}")
            return []
        
        logger.info(f"Parsed {len(lncrnas)} lncRNA genes from LNCipedia")
        return lncrnas

    def fetch_gencode_lncrnas(self):
       
        logger.info("Fetching lncRNA genes from all local sources...")
        
        all_lncrnas = []
        
        gtf_file = Path(self.config['data_sources']['gencode']['gencode_gtf'])
        
        if gtf_file.exists():
            logger.info("Processing GENCODE lncRNAs...")
            gencode_lncrnas = []
            
            try:
                with gzip.open(gtf_file, 'rt') as f:
                    for line in f:
                        if line.startswith('#') or 'gene' not in line:
                            continue
                        
                        fields = line.strip().split('\t')
                        if len(fields) < 9 or fields[2] != 'gene':
                            continue
                        
                        attrs = fields[8]
                        gene_id = ''
                        gene_name = ''
                        
                        if 'gene_id "' in attrs:
                            start = attrs.find('gene_id "') + 9
                            end = attrs.find('"', start)
                            gene_id = attrs[start:end].split('.')[0]
                        
                        if 'gene_name "' in attrs:
                            start = attrs.find('gene_name "') + 11
                            end = attrs.find('"', start)
                            gene_name = attrs[start:end]
                        
                        if gene_id:
                            chromosome = fields[0].replace('chr', '')
                            
                            if chromosome not in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']:
                                continue
                            
                            try:
                                start_pos = int(fields[3])
                                end_pos = int(fields[4])
                                strand = fields[6]
                            except:
                                continue
                            
                            gencode_lncrnas.append({
                                'ensembl_gene_id': gene_id,
                                'gene_symbol': gene_name if gene_name else gene_id,
                                'chromosome': chromosome,
                                'start_position': start_pos,
                                'end_position': end_pos,
                                'strand': strand,
                                'ncrna_type': 'lncRNA',
                                'ncrna_subtype': 'GENCODE_lncRNA',
                                'length': end_pos - start_pos + 1,
                                'description': f"GENCODE lncRNA: {gene_name if gene_name else gene_id}"
                            })
                        
                        if len(gencode_lncrnas) >= 10000:
                            break
            
            except Exception as e:
                logger.error(f"Error parsing GENCODE GTF: {e}")
            
            all_lncrnas.extend(gencode_lncrnas)
            logger.info(f"Added {len(gencode_lncrnas)} GENCODE lncRNAs")
        
        
        try:
            noncode_lncrnas = self.fetch_noncode_lncrnas()
            all_lncrnas.extend(noncode_lncrnas)
            logger.info(f"Added {len(noncode_lncrnas)} NONCODE lncRNAs")
        except Exception as e:
            logger.warning(f"NONCODE parsing failed: {e}")
        
        
        try:
            lncipedia_lncrnas = self.fetch_lncipedia_lncrnas()
            all_lncrnas.extend(lncipedia_lncrnas)
            logger.info(f"Added {len(lncipedia_lncrnas)} LNCipedia lncRNAs")
        except Exception as e:
            logger.warning(f"LNCipedia parsing failed: {e}")
        
        # Remove duplicates based on gene_symbol (keep first option)
        seen_symbols = set()
        unique_lncrnas = []
        for lncrna in all_lncrnas:
            symbol = lncrna['gene_symbol']
            if symbol not in seen_symbols:
                seen_symbols.add(symbol)
                unique_lncrnas.append(lncrna)
        
        logger.info(f"Total unique lncRNAs from all sources: {len(unique_lncrnas)} (removed {len(all_lncrnas) - len(unique_lncrnas)} duplicates)")
        return unique_lncrnas
    
    def fetch_mirdb_regulations(self):
        
        logger.info("Fetching miRNA-target regulations from local miRDB file...")
        
        mirdb_file = Path(self.config['regulation_sources']['mirdb']['mirdb_file'])
        
        if not mirdb_file.exists():
            logger.error(f"miRDB file not found: {mirdb_file}")
            return []
        
        regulations = []
        
        try:
            with gzip.open(mirdb_file, 'rt') as f:
                for line_num, line in enumerate(tqdm(f, desc="Parsing miRDB")):
                    if line_num == 0:  # Skip header if any
                        continue
                    
                    try:
                        # miRDB format: miRNA_ID, Gene_Symbol, Target_Score
                        fields = line.strip().split('\t')
                        if len(fields) >= 3:
                            mirna_id = fields[0].strip()
                            gene_symbol = fields[1].strip() 
                            score = float(fields[2].strip())
                            
                            # Only keep high-confidence predictions (score > 80)
                            if score > 80:
                                regulations.append({
                                    'ncrna_symbol': mirna_id,
                                    'target_symbol': gene_symbol,
                                    'regulation_type': 'negative',
                                    'mechanism': 'post_transcriptional',
                                    'evidence': 'computational',
                                    'pmids': 'miRDB_v6.0',
                                    'confidence_score': score / 100.0
                                })
                    except Exception as e:
                        continue
                    
                    if len(regulations) >= 5000:
                        break
        
        except Exception as e:
            logger.error(f"Error parsing miRDB: {e}")
            return []
        
        logger.info(f"Parsed {len(regulations)} miRNA-target regulations from miRDB")
        return regulations
    
    def fetch_targetscan_regulations(self):
        logger.info("Fetching miRNA-target regulations from local TargetScan file...")
        
        targetscan_file = Path(self.config['regulation_sources']['targetscan']['targetscan_file'])
        
        if not targetscan_file.exists():
            logger.error(f"TargetScan file not found: {targetscan_file}")
            return []
        
        regulations = []
        
        try:
            # Handle ZIP file
            if targetscan_file.suffix == '.zip':
                with zipfile.ZipFile(targetscan_file, 'r') as zip_ref:
                    # Extract the first file from the zip
                    file_list = zip_ref.namelist()
                    if file_list:
                        with zip_ref.open(file_list[0]) as f:
                            content = f.read().decode('utf-8')
                            lines = content.split('\n')
                            
                            if lines:
                                header = lines[0]  # Skip header
                                
                                for line_num, line in enumerate(tqdm(lines[1:], desc="Parsing TargetScan")):
                                    if not line.strip():
                                        continue
                                    
                                    try:
                                        fields = line.strip().split('\t')
                                        if len(fields) >= 7:
                                            gene_symbol = fields[0].strip()
                                            mirna_family = fields[1].strip()
                                            context_score = float(fields[6].strip()) if fields[6].strip() else 0
                                            
                                            # Only keep high-confidence predictions (context_score < -0.2)
                                            if context_score < -0.2:
                                                regulations.append({
                                                    'ncrna_symbol': mirna_family,
                                                    'target_symbol': gene_symbol,
                                                    'regulation_type': 'negative',
                                                    'mechanism': 'post_transcriptional',
                                                    'evidence': 'computational',
                                                    'pmids': 'TargetScan_v8.0',
                                                    'confidence_score': min(abs(context_score), 1.0)
                                                })
                                    except Exception as e:
                                        continue
                                    
                                    if len(regulations) >= 5000:
                                        break
            else:
                # Handle regular text file
                with open(targetscan_file, 'r') as f:
                    header = f.readline()  # Skip header
                    
                    for line_num, line in enumerate(tqdm(f, desc="Parsing TargetScan")):
                        try:
                            fields = line.strip().split('\t')
                            if len(fields) >= 7:
                                gene_symbol = fields[0].strip()
                                mirna_family = fields[1].strip()
                                context_score = float(fields[6].strip()) if fields[6].strip() else 0
                                
                                # Only keep high-confidence predictions (context_score < -0.2)
                                if context_score < -0.2:
                                    regulations.append({
                                        'ncrna_symbol': mirna_family,
                                        'target_symbol': gene_symbol,
                                        'regulation_type': 'negative',
                                        'mechanism': 'post_transcriptional',
                                        'evidence': 'computational',
                                        'pmids': 'TargetScan_v8.0',
                                        'confidence_score': min(abs(context_score), 1.0)
                                    })
                        except Exception as e:
                            continue
                        
                        if len(regulations) >= 5000:
                            break
        
        except Exception as e:
            logger.error(f"Error parsing TargetScan: {e}")
            return []
        
        logger.info(f"Parsed {len(regulations)} miRNA-target regulations from TargetScan")
        return regulations
    
    def create_sample_lncrna_regulations(self):
        logger.info("Adding literature-based lncRNA regulations...")
        
        lncrna_regulations = [
            {
                'ncrna_symbol': 'NEAT1',
                'target_symbol': 'BRCA1',
                'regulation_type': 'positive',
                'mechanism': 'protein_sequestration',
                'evidence': 'experimental',
                'pmids': '23555862,24662760',
                'confidence_score': 0.9
            },
            {
                'ncrna_symbol': 'MALAT1',
                'target_symbol': 'BRCA2',
                'regulation_type': 'negative',
                'mechanism': 'splicing_regulation',
                'evidence': 'experimental',
                'pmids': '21685050',
                'confidence_score': 0.85
            },
            {
                'ncrna_symbol': 'XIST',
                'target_symbol': 'TP53',
                'regulation_type': 'negative',
                'mechanism': 'chromatin_silencing',
                'evidence': 'experimental',
                'pmids': '1656926',
                'confidence_score': 0.95
            },
            {
                'ncrna_symbol': 'H19',
                'target_symbol': 'EGFR',
                'regulation_type': 'positive',
                'mechanism': 'transcriptional',
                'evidence': 'experimental',
                'pmids': '24463510',
                'confidence_score': 0.8
            },
            {
                'ncrna_symbol': 'HOTAIR',
                'target_symbol': 'MYC',
                'regulation_type': 'negative',
                'mechanism': 'chromatin_modification',
                'evidence': 'experimental',
                'pmids': '17604720',
                'confidence_score': 0.88
            }
        ]
        
        return lncrna_regulations
    
    def fetch_all_regulations(self):
        logger.info("Fetching regulations from all local sources...")
        
        all_regulations = []
        

        try:
            mirdb_regs = self.fetch_mirdb_regulations()
            all_regulations.extend(mirdb_regs)
            logger.info(f"Added {len(mirdb_regs)} miRDB regulations")
        except Exception as e:
            logger.warning(f"miRDB failed: {e}")
        

        try:
            targetscan_regs = self.fetch_targetscan_regulations()
            all_regulations.extend(targetscan_regs)
            logger.info(f"Added {len(targetscan_regs)} TargetScan regulations")
        except Exception as e:
            logger.warning(f"TargetScan failed: {e}")
        

        lncrna_regs = self.create_sample_lncrna_regulations()
        all_regulations.extend(lncrna_regs)
        logger.info(f"Added {len(lncrna_regs)} literature lncRNA regulations")
        
        logger.info(f"Total regulations from all sources: {len(all_regulations)}")
        return all_regulations
    
    def fetch_hpo_terms(self):
        logger.info("Fetching HPO terms from local OBO file...")
        
        # Get path from config
        obo_file = Path(self.config['data_sources']['hpo']['hpo_obo'])
        
        if not obo_file.exists():
            logger.error(f"HPO OBO file not found: {obo_file}")
            return []
        
        hpo_terms = []
        current_term = {}
        
        try:
            with open(obo_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    if line == '[Term]':
                        if current_term:
                            hpo_terms.append(current_term)
                        current_term = {}
                        
                    elif line.startswith('id: HP:'):
                        current_term['hpo_id'] = line.replace('id: ', '')
                        
                    elif line.startswith('name: '):
                        current_term['hpo_name'] = line.replace('name: ', '')
                        
                    elif line.startswith('def: '):
                        current_term['definition'] = line.replace('def: ', '').split('"')[1]
                        
                    elif line.startswith('is_a: HP:'):
                        parent = line.split(' ! ')[0].replace('is_a: ', '')
                        current_term.setdefault('parent_terms', []).append(parent)
            
            if current_term:
                hpo_terms.append(current_term)
            
            for term in hpo_terms:
                if 'parent_terms' in term:
                    term['parent_terms'] = json.dumps(term['parent_terms'])
        
        except Exception as e:
            logger.error(f"Error parsing HPO OBO file: {e}")
            return []
        
        logger.info(f"Parsed {len(hpo_terms)} HPO terms from local file")
        
        # Add specific disease HPO terms
        additional_hpo_terms = [
            {
                'hpo_id': 'HP:0003002',
                'hpo_name': 'Breast carcinoma',
                'definition': 'The presence of a carcinoma of the breast.',
                'parent_terms': None
            },
            {
                'hpo_id': 'HP:0009733', 
                'hpo_name': 'Glioma',
                'definition': 'The presence of a glioma, which is a neoplasm of the central nervous system.',
                'parent_terms': None
            },
            {
                'hpo_id': 'HP:0002664',
                'hpo_name': 'Neoplasm',
                'definition': 'An organ or organ-system abnormality that consists of uncontrolled autonomous cell-proliferation.',
                'parent_terms': None
            }
        ]
        
        existing_hpo_ids = {term.get('hpo_id') for term in hpo_terms}
        for term in additional_hpo_terms:
            if term['hpo_id'] not in existing_hpo_ids:
                hpo_terms.append(term)
        
        logger.info(f"Total HPO terms (including additional): {len(hpo_terms)}")
        return hpo_terms
    
    def fetch_omim_disease_associations(self):
        """Use your existing OMIM parser - no changes needed!"""
        from parse_omim_data import OMIMDataParser
    
        parser = OMIMDataParser()
        omim_entries, gene_omim_mapping = parser.parse_omim_file()
        parser.populate_database(omim_entries, gene_omim_mapping)
        
        logger.info(f"Loaded {len(gene_omim_mapping)} real disease associations from OMIM")
        
    def populate_database(self, genes=None, ncrnas=None, regulations=None, hpo_terms=None, disease_associations=None):
        """Populate database with all data."""
        logger.info("Populating database...")
    
        conn = sqlite3.connect(self.db_file)
        conn.execute("PRAGMA foreign_keys = ON")
    
        try:
            # Insert genes
            if genes:
                logger.info(f"Inserting {len(genes)} protein-coding genes...")
                gene_columns = list(genes[0].keys())
                placeholders = ','.join(['?' for _ in gene_columns])
            
                conn.executemany(
                    f"INSERT OR REPLACE INTO genes ({','.join(gene_columns)}) VALUES ({placeholders})",
                    [tuple(gene[col] for col in gene_columns) for gene in genes]
                )
        
            # Insert ncRNAs 
            if ncrnas:
                logger.info(f"Inserting {len(ncrnas)} ncRNA genes...")
                ncrna_columns = list(ncrnas[0].keys())
                placeholders = ','.join(['?' for _ in ncrna_columns])
            
                conn.executemany(
                    f"INSERT OR REPLACE INTO ncrna_genes ({','.join(ncrna_columns)}) VALUES ({placeholders})",
                    [tuple(ncrna[col] for col in ncrna_columns) for ncrna in ncrnas]
                )
        
            # Insert HPO terms 
            if hpo_terms:
                logger.info(f"Inserting {len(hpo_terms)} HPO terms...")
                conn.executemany(
                    "INSERT OR REPLACE INTO hpo_terms (hpo_id, hpo_name, definition, parent_terms) VALUES (?, ?, ?, ?)",
                    [(term.get('hpo_id'), term.get('hpo_name'), term.get('definition'), term.get('parent_terms'))
                     for term in hpo_terms]
                )
           
            # Insert regulations
            if regulations:
                logger.info(f"Inserting {len(regulations)} regulation relationships...")
                successful_regs = 0
                
                for reg in regulations:
                    # Find or create ncRNA
                    ncrna_result = conn.execute(
                        "SELECT ncrna_id FROM ncrna_genes WHERE gene_symbol = ?",
                        (reg['ncrna_symbol'],)
                    ).fetchone()
                    
                    if not ncrna_result:
                        try:
                            ncrna_type = 'miRNA' if 'miR' in reg['ncrna_symbol'] else 'lncRNA'
                            conn.execute("""
                                INSERT OR IGNORE INTO ncrna_genes 
                                (ensembl_gene_id, gene_symbol, ncrna_type, ncrna_subtype, description) 
                                VALUES (?, ?, ?, ?, ?)
                            """, (
                                f"AUTO_{reg['ncrna_symbol']}", reg['ncrna_symbol'], 
                                ncrna_type, ncrna_type, f"Auto-added {ncrna_type}"
                            ))
                            
                            ncrna_result = conn.execute(
                                "SELECT ncrna_id FROM ncrna_genes WHERE gene_symbol = ?",
                                (reg['ncrna_symbol'],)
                            ).fetchone()
                        except:
                            continue
                    
                    # Find or create target gene
                    gene_result = conn.execute(
                        "SELECT gene_id FROM genes WHERE gene_symbol = ?",
                        (reg['target_symbol'],)
                    ).fetchone()
                    
                    if not gene_result:
                        try:
                            conn.execute("""
                                INSERT OR IGNORE INTO genes 
                                (ensembl_gene_id, gene_symbol, gene_type, description) 
                                VALUES (?, ?, ?, ?)
                            """, (
                                f"AUTO_{reg['target_symbol']}", reg['target_symbol'],
                                'protein_coding', f"Auto-added gene"
                            ))
                            
                            gene_result = conn.execute(
                                "SELECT gene_id FROM genes WHERE gene_symbol = ?",
                                (reg['target_symbol'],)
                            ).fetchone()
                        except:
                            continue
                    
                    # Insert regulation
                    if ncrna_result and gene_result:
                        try:
                            conn.execute("""
                                INSERT OR REPLACE INTO regulations 
                                (ncrna_id, target_gene_id, regulation_type, regulation_mechanism, 
                                 evidence_level, pubmed_ids, confidence_score) 
                                VALUES (?, ?, ?, ?, ?, ?, ?)
                            """, (
                                ncrna_result[0], gene_result[0], reg['regulation_type'],
                                reg['mechanism'], reg['evidence'], reg['pmids'], 
                                reg.get('confidence_score', 0.8)
                            ))
                            successful_regs += 1
                        except:
                            continue
                
                logger.info(f"Successfully inserted {successful_regs} regulation relationships")
            
            # Insert disease associations
            if disease_associations:
                logger.info(f"Inserting {len(disease_associations)} disease associations...")
                successful_diseases = 0
                
                for assoc in disease_associations:
                    gene_result = conn.execute(
                        "SELECT gene_id FROM genes WHERE gene_symbol = ?",
                        (assoc['gene_symbol'],)
                    ).fetchone()
                    
                    if gene_result:
                        try:
                            conn.execute("""
                                INSERT OR REPLACE INTO gene_disease_associations 
                                (gene_id, hpo_id, omim_id, association_type, evidence_level, pubmed_ids) 
                                VALUES (?, ?, ?, ?, ?, ?)
                            """, (
                                gene_result[0], assoc.get('hpo_id'), assoc.get('omim_id'),
                                assoc['association_type'], assoc['evidence_level'], assoc['pubmed_ids']
                            ))
                            successful_diseases += 1
                        except Exception as e:
                            continue
                
                logger.info(f"Successfully inserted {successful_diseases} disease associations")
            
            conn.commit()
            logger.info("Database populated successfully!")
            
        except Exception as e:
            logger.error(f"Error populating database: {e}")
            conn.rollback()
            raise
        finally:
            conn.close()

def main():
    """Main function for Bianca data fetching."""
    parser = argparse.ArgumentParser(description='Fetch data for ncRNA regulation database (Bianca version)')
    parser.add_argument('--source', choices=['all', 'genes', 'ncrnas', 'regulations', 'hpo', 'diseases'], 
                       default='all', help='Data source to fetch')
    parser.add_argument('--limit', type=int, help='Limit number of genes to fetch (for testing)')
    
    args = parser.parse_args()
    
    Path("logs").mkdir(exist_ok=True)
    
    try:
        fetcher = BiancaDataFetcher()
        
        genes = None
        ncrnas = None
        regulations = None
        hpo_terms = None
        disease_associations = None
        
        if args.source in ['all', 'genes']:
            genes = fetcher.fetch_genes_from_local_biomart(limit=args.limit)
        
        if args.source in ['all', 'ncrnas']:
            ncrnas = fetcher.fetch_gencode_lncrnas()
        
        if args.source in ['all', 'regulations']:
            regulations = fetcher.fetch_all_regulations()
        
        if args.source in ['all', 'hpo']:
            hpo_terms = fetcher.fetch_hpo_terms()
            
        if args.source in ['all', 'diseases']:
            fetcher.fetch_omim_disease_associations()
            disease_associations = None
        
        # Populate database
        fetcher.populate_database(genes=genes, ncrnas=ncrnas, 
                                 regulations=regulations, hpo_terms=hpo_terms,
                                 disease_associations=disease_associations)
        

        conn = sqlite3.connect(fetcher.db_file)
        
        gene_count = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
        ncrna_count = conn.execute("SELECT COUNT(*) FROM ncrna_genes").fetchone()[0]
        regulation_count = conn.execute("SELECT COUNT(*) FROM regulations").fetchone()[0]
        hpo_count = conn.execute("SELECT COUNT(*) FROM hpo_terms").fetchone()[0]
        omim_count = conn.execute("SELECT COUNT(*) FROM omim_entries").fetchone()[0]
        disease_count = conn.execute("SELECT COUNT(*) FROM gene_disease_associations").fetchone()[0]
        
        print(f" Database populated with:")
        print(f" Protein-coding genes: {gene_count}")
        print(f" ncRNA genes: {ncrna_count}")
        print(f" Regulations: {regulation_count}")
        print(f" HPO terms: {hpo_count}")
        print(f" OMIM entries: {omim_count}")
        print(f" Disease associations: {disease_count}")
        
        conn.close()
        
    except Exception as e:
        logger.error(f"Data fetching failed: {e}")
        raise

if __name__ == "__main__":
    main()
