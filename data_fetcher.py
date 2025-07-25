"""
Data fetcher with fast batch processing for regulations
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
from collections import defaultdict


from gene_mapper import PersistentGeneMapper

# Simple logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class BiancaDataFetcher:
    def __init__(self, config_path="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/config.yaml"):
        self.config_path = config_path
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_file = Path(self.config['database']['path']) / self.config['database']['name']
        self.data_dir = Path("data/raw")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Load thresholds from config
        self.thresholds = self.config.get('thresholds', {})
        self.limits = self.config.get('limits', {})
        
        # Initialize  gene mapper
        logger.info("Initializing enhanced gene mapper...")
        self.gene_mapper = None
        self.gene_lookup = {}
        self.ncrna_lookup = {}
    
    def _initialize_gene_mapper(self):
        """Initialize gene mapper and build fast lookup dictionaries"""
        if self.gene_mapper is None:
            try:
                self.gene_mapper = PersistentGeneMapper(str(self.db_file), self.config_path)
                # Build fresh mappings to include newly populated data
                self.gene_mapper.build_and_save_mappings(force_rebuild=True)
                logger.info("Enhanced gene mapper initialized successfully")
                
                # Build fast lookup dictionaries
                self._build_lookup_dictionaries()
                
            except Exception as e:
                logger.warning("Could not initialize gene mapper: {}. Using basic matching.".format(e))
                self.gene_mapper = None
    
    def _build_lookup_dictionaries(self):
        """Build fast lookup dictionaries for genes and ncRNAs - OPTIMIZED VERSION"""
        logger.info("Building fast lookup dictionaries...")
        
        max_gene_mappings = self.limits.get('max_gene_mappings', 50000)  # Reduced default
        max_ncrna_mappings = self.limits.get('max_ncrna_mappings', 10000)  # Reduced default
        max_variations = self.limits.get('max_lookup_variations', 3)  # Reduced default
        
        conn = sqlite3.connect(self.db_file)
        try:
            # Build gene lookup with CHUNKED processing
            logger.info("Processing gene mappings...")
            cursor = conn.execute("SELECT gene_id, gene_symbol FROM genes LIMIT {}".format(max_gene_mappings))
            
            gene_count = 0
            batch_size = 1000
            batch_genes = []
            
            for gene_id, gene_symbol in cursor:
                batch_genes.append((gene_id, gene_symbol))
                
                if len(batch_genes) >= batch_size:
                    self._process_gene_batch(batch_genes, max_variations)
                    gene_count += len(batch_genes)
                    batch_genes = []
                    
                    if gene_count % 5000 == 0:
                        logger.info("Processed {} genes...".format(gene_count))
                    
                    if gene_count >= max_gene_mappings:
                        break
            
            # Process remaining batch
            if batch_genes:
                self._process_gene_batch(batch_genes, max_variations)
                gene_count += len(batch_genes)
            
            # Build ncRNA lookup with CHUNKED processing  
            logger.info("Processing ncRNA mappings...")
            cursor = conn.execute("SELECT ncrna_id, gene_symbol FROM ncrna_genes LIMIT {}".format(max_ncrna_mappings))
            
            ncrna_count = 0
            batch_ncrnas = []
            
            for ncrna_id, gene_symbol in cursor:
                batch_ncrnas.append((ncrna_id, gene_symbol))
                
                if len(batch_ncrnas) >= batch_size:
                    self._process_ncrna_batch(batch_ncrnas, max_variations)
                    ncrna_count += len(batch_ncrnas)
                    batch_ncrnas = []
                    
                    if ncrna_count % 2000 == 0:
                        logger.info("Processed {} ncRNAs...".format(ncrna_count))
                    
                    if ncrna_count >= max_ncrna_mappings:
                        break
            
            # Process remaining batch
            if batch_ncrnas:
                self._process_ncrna_batch(batch_ncrnas, max_variations)
                ncrna_count += len(batch_ncrnas)
            
            logger.info("Built lookup dictionaries: {} gene variations, {} ncRNA variations".format(
                len(self.gene_lookup), len(self.ncrna_lookup)))
                
        finally:
            conn.close()

    def _process_gene_batch(self, batch_genes, max_variations):
        """Process a batch of genes efficiently"""
        for gene_id, gene_symbol in batch_genes:
            if not gene_symbol:
                continue
                
            # Create basic variations without expensive mapper calls
            variations = [
                gene_symbol,
                gene_symbol.upper(), 
                gene_symbol.lower(),
                gene_symbol.replace('-', ''),
                gene_symbol.replace('_', '')
            ]
            
            # Only call mapper if it exists and we haven't exceeded limits
            if self.gene_mapper and len(variations) < max_variations:
                try:
                    # Quick lookup - avoid expensive operations
                    standard_symbol = self.gene_mapper.find_gene_symbol(gene_symbol)
                    if standard_symbol and standard_symbol != gene_symbol:
                        variations.append(standard_symbol)
                except:
                    pass
            
            # Add to lookup (limit variations to prevent memory bloat)
            for var in variations[:max_variations]:
                if var and var not in self.gene_lookup:
                    self.gene_lookup[var] = gene_id

    def _process_ncrna_batch(self, batch_ncrnas, max_variations):
        """Process a batch of ncRNAs efficiently"""
        for ncrna_id, gene_symbol in batch_ncrnas:
            if not gene_symbol:
                continue
                
            # Create basic variations without expensive mapper calls
            variations = [
                gene_symbol,
                gene_symbol.upper(),
                gene_symbol.lower()
            ]
            
            # Add miRNA-specific variations if applicable
            if 'mir' in gene_symbol.lower():
                variations.extend([
                    gene_symbol.replace('mir', 'miR'),
                    gene_symbol.replace('miR', 'mir'),
                    'hsa-' + gene_symbol.replace('hsa-', '')
                ])
            
            # Only call mapper if it exists and we haven't exceeded limits  
            if self.gene_mapper and len(variations) < max_variations:
                try:
                    standard_symbol = self.gene_mapper.find_ncrna_symbol(gene_symbol)
                    if standard_symbol and standard_symbol != gene_symbol:
                        variations.append(standard_symbol)
                except:
                    pass
            
            # Add to lookup (limit variations)
            for var in variations[:max_variations]:
                if var and var not in self.ncrna_lookup:
                    self.ncrna_lookup[var] = ncrna_id
    
    def fetch_genes_from_local_biomart(self, limit=None):
        """Fetch genes from local BioMart file."""
        logger.info("Fetching protein-coding genes from local BioMart file...")
        
        biomart_file = Path(self.config['data_sources']['general_genes_mart'])
        
        if not biomart_file.exists():
            logger.error("BioMart file not found: {}".format(biomart_file))
            return []
        
        genes = []
        max_genes = limit or self.limits.get('max_genes_per_source', 50000)
        
        try:
            with open(biomart_file, 'r', encoding='utf-8') as f:
                for i, line in enumerate(tqdm(f, desc="Processing BioMart file")):
                    if max_genes > 0 and i >= max_genes:
                        break
                    
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        # Expected format: ensembl_id, gene_symbol, description, chromosome, start, end, strand
                        fields = line.split('\t')
                        if len(fields) >= 7:
                            ensembl_id, gene_symbol, description, chromosome, start, end, strand = fields[:7]
                            
                            if not ensembl_id or not gene_symbol:
                                continue
                            
                            # Clean chromosome name
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
            logger.error("Error reading BioMart file: {}".format(e))
            return []
        
        logger.info("Fetched {} genes from local BioMart file".format(len(genes)))
        return genes
    
    def fetch_gencode_lncrnas_combined(self):
        """Fetch lncRNA genes from all local sources: GENCODE + NONCODE + LNCipedia."""
        logger.info("Fetching lncRNA genes from all local sources...")
        
        all_lncrnas = []
        max_ncrnas = self.limits.get('max_ncrnas_per_source', 25000)
        
        # Process GENCODE
        gtf_file = Path(self.config['data_sources']['gencode']['gencode_gtf'])
        if gtf_file.exists():
            logger.info("Processing GENCODE lncRNAs...")
            try:
                with gzip.open(gtf_file, 'rt') as f:
                    for line in f:
                        if line.startswith('#') or 'gene' not in line:
                            continue
                        
                        if max_ncrnas > 0 and len(all_lncrnas) >= max_ncrnas:
                            break
                        
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
                            if chromosome in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']:
                                try:
                                    all_lncrnas.append({
                                        'ensembl_gene_id': gene_id,
                                        'gene_symbol': gene_name if gene_name else gene_id,
                                        'chromosome': chromosome,
                                        'start_position': int(fields[3]),
                                        'end_position': int(fields[4]),
                                        'strand': fields[6],
                                        'ncrna_type': 'lncRNA',
                                        'ncrna_subtype': 'GENCODE_lncRNA',
                                        'length': int(fields[4]) - int(fields[3]) + 1,
                                        'description': "GENCODE lncRNA: {}".format(gene_name if gene_name else gene_id)
                                    })
                                except:
                                    continue
            except Exception as e:
                logger.error("Error parsing GENCODE GTF: {}".format(e))
            
            logger.info("Added {} GENCODE lncRNAs".format(len(all_lncrnas)))
        
        # Remove duplicates
        seen_symbols = set()
        unique_lncrnas = []
        for lncrna in all_lncrnas:
            symbol = lncrna['gene_symbol']
            if symbol not in seen_symbols:
                seen_symbols.add(symbol)
                unique_lncrnas.append(lncrna)
        
        logger.info("Total unique lncRNAs: {} (removed {} duplicates)".format(
            len(unique_lncrnas), len(all_lncrnas) - len(unique_lncrnas)))
        return unique_lncrnas
        
    def fetch_mirdb_regulations(self, score_threshold=None):
        """Fetch miRNA-target regulations from local miRDB file."""
        logger.info("Fetching miRNA-target regulations from local miRDB file...")
        
        mirdb_file = Path(self.config['regulation_sources']['mirdb']['mirdb_file'])
        
        if not mirdb_file.exists():
            logger.error("miRDB file not found: {}".format(mirdb_file))
            return []
        
        # Use provided threshold or config threshold
        if score_threshold is None:
            score_threshold = self.thresholds.get('mirdb_score_threshold', 80)
        
        regulations = []
        max_regs = self.limits.get('max_regulations_per_source', 100000)
        
        try:
            with gzip.open(mirdb_file, 'rt') as f:
                for line_num, line in enumerate(tqdm(f, desc="Parsing miRDB")):
                    if line_num == 0:  # Skip header if any
                        continue
                    
                    if max_regs > 0 and len(regulations) >= max_regs:
                        break
                    
                    try:
                        fields = line.strip().split('\t')
                        if len(fields) >= 3:
                            mirna_id = fields[0].strip()
                            gene_symbol = fields[1].strip() 
                            score = float(fields[2].strip())
                            
                            if score > score_threshold:
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
        
        except Exception as e:
            logger.error("Error parsing miRDB: {}".format(e))
            return []
        
        logger.info("Parsed {} miRNA-target regulations from miRDB (threshold: {})".format(
            len(regulations), score_threshold))
        return regulations
    
    def fetch_targetscan_regulations(self, context_threshold=None):
        """Fetch miRNA-target regulations from local TargetScan file."""
        logger.info("Fetching miRNA-target regulations from local TargetScan file...")

        targetscan_file = Path(self.config['regulation_sources']['targetscan']['targetscan_file'])
        
        if not targetscan_file.exists():
            logger.error("TargetScan file not found: {}".format(targetscan_file))
            return []

        if context_threshold is None:
            context_threshold = self.thresholds.get('targetscan_context_score_threshold', -0.2)
        
        regulations = []
        max_regs = self.limits.get('max_regulations_per_source', 100000)
        
        try:
            # Handle ZIP file
            if targetscan_file.suffix == '.zip':
                with zipfile.ZipFile(targetscan_file, 'r') as zip_ref:
                    file_list = zip_ref.namelist()
                    if file_list:
                        with zip_ref.open(file_list[0]) as f:
                            content = f.read().decode('utf-8')
                            lines = content.split('\n')
                            
                            if lines:
                                for line_num, line in enumerate(tqdm(lines[1:], desc="Parsing TargetScan")):
                                    if not line.strip():
                                        continue
                                    
                                    if max_regs > 0 and len(regulations) >= max_regs:
                                        break
                                    
                                    try:
                                        fields = line.strip().split('\t')
                                        if len(fields) >= 7:
                                            gene_symbol = fields[0].strip()
                                            mirna_family = fields[1].strip()
                                            context_score = float(fields[6].strip()) if fields[6].strip() else 0
                                            
                                            if context_score < context_threshold:
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
        
        except Exception as e:
            logger.error("Error parsing TargetScan: {}".format(e))
            return []
        
        logger.info("Parsed {} miRNA-target regulations from TargetScan (threshold: {})".format(
            len(regulations), context_threshold))
        return regulations
    
    def fetch_evlncrnas_regulations(self, confidence_threshold=None):
        """Fetch lncRNA-target regulations from EVLncRNAs3 data."""
        logger.info("Fetching lncRNA-target regulations from EVLncRNAs3...")
        
        evlncrnas_file = Path(self.config['regulation_sources']['evlncrnas']['evlncrnas_file'])
        
        if not evlncrnas_file.exists():
            logger.error("EVLncRNAs file not found: {}".format(evlncrnas_file))
            return []
        
        if confidence_threshold is None:
            confidence_threshold = self.thresholds.get('evlncrnas_confidence_threshold', 0.5)
        
        regulations = []
        max_regs = self.limits.get('max_regulations_per_source', 100000)
        
        try:
            import pandas as pd
            
            with zipfile.ZipFile(evlncrnas_file, 'r') as zip_ref:
                file_list = zip_ref.namelist()
                logger.info("Files in EVLncRNAs ZIP: {}".format(file_list))
                
                # Try Excel files first, then text files
                excel_files = [f for f in file_list if f.endswith(('.xls', '.xlsx')) 
                              and not f.startswith('__MACOSX') and not '._' in f]
                text_files = [f for f in file_list if f.endswith(('.csv', '.tsv', '.txt', '.dat')) 
                             and not f.startswith('__MACOSX') and not '._' in f]
                
                all_files = excel_files + text_files
                
                if not all_files:
                    logger.warning("No readable data files found in EVLncRNAs ZIP")
                    return []
                
                # Try each file until we find one that works
                for data_file in all_files:
                    logger.info("Trying EVLncRNAs file: {}".format(data_file))
                    
                    try:
                        if data_file.endswith(('.xls', '.xlsx')):
                            # Read Excel file with error handling
                            try:
                                with zip_ref.open(data_file) as f:
                                    # Try pandas with different engines
                                    try:
                                        df = pd.read_excel(f, engine='openpyxl' if data_file.endswith('.xlsx') else 'xlrd')
                                    except ImportError as e:
                                        logger.warning("Excel library not available ({}), skipping Excel file: {}".format(e, data_file))
                                        continue
                                    except Exception as e:
                                        logger.warning("Could not read Excel file {}: {}".format(data_file, e))
                                        continue
                                    
                                    if df.empty or len(df.columns) < 2:
                                        logger.info("Excel file {} is empty or has insufficient columns".format(data_file))
                                        continue
                                    
                                    headers = [h.strip().lower() for h in df.columns]
                                    logger.info("Excel headers found: {}".format(headers[:10]))
                                    
                                    # Find relevant columns (more flexible matching)
                                    lncrna_col = None
                                    target_col = None
                                    regulation_col = None
                                    
                                    for i, header in enumerate(headers):
                                        if any(keyword in header for keyword in ['lncrna', 'longnoncoding', 'rna_name', 'regulator', 'symbol']):
                                            if lncrna_col is None:  # Take first match
                                                lncrna_col = i
                                        elif any(keyword in header for keyword in ['target', 'gene']) and 'lncrna' not in header:
                                            if target_col is None:  # Take first match
                                                target_col = i
                                        elif any(keyword in header for keyword in ['regulation', 'effect', 'direction', 'function']):
                                            if regulation_col is None:  # Take first match
                                                regulation_col = i
                                    
                                    # If we can't find proper columns, try first two columns as defaults
                                    if lncrna_col is None and len(df.columns) >= 1:
                                        lncrna_col = 0
                                        logger.info("Using first column as lncRNA column")
                                    if target_col is None and len(df.columns) >= 2:
                                        target_col = 1  
                                        logger.info("Using second column as target column")
                                    
                                    if lncrna_col is None or target_col is None:
                                        logger.info("Could not identify required columns in {}".format(data_file))
                                        continue
                                    
                                    logger.info("Using Excel columns - lncRNA: {}, target: {}, regulation: {}".format(
                                        lncrna_col, target_col, regulation_col))
                                    
                                    # Process rows with better error handling
                                    for idx, row in df.iterrows():
                                        if max_regs > 0 and len(regulations) >= max_regs:
                                            break
                                        
                                        try:
                                            lncrna_symbol = str(row.iloc[lncrna_col]).strip()
                                            target_symbol = str(row.iloc[target_col]).strip()
                                            
                                            # Skip invalid entries
                                            if lncrna_symbol in ['nan', 'NaN', ''] or target_symbol in ['nan', 'NaN', '']:
                                                continue
                                            if not lncrna_symbol or not target_symbol:
                                                continue
                                            
                                            # Determine regulation type
                                            regulation_type = 'unknown'
                                            if regulation_col is not None and regulation_col < len(row):
                                                reg_text = str(row.iloc[regulation_col]).lower()
                                                if any(word in reg_text for word in ['positive', 'activate', 'enhance', 'promote', 'up']):
                                                    regulation_type = 'positive'
                                                elif any(word in reg_text for word in ['negative', 'inhibit', 'repress', 'suppress', 'down']):
                                                    regulation_type = 'negative'
                                            
                                            regulations.append({
                                                'ncrna_symbol': lncrna_symbol,
                                                'target_symbol': target_symbol,
                                                'regulation_type': regulation_type,
                                                'mechanism': 'transcriptional',
                                                'evidence': 'experimental',
                                                'pmids': 'EVLncRNAs3.0',
                                                'confidence_score': confidence_threshold
                                            })
                                            
                                        except Exception as e:
                                            continue
                                    
                                    # If we got results from this file, break
                                    if regulations:
                                        logger.info("Successfully parsed {} regulations from {}".format(
                                            len(regulations), data_file))
                                        break
                                        
                            except Exception as e:
                                logger.warning("Error processing Excel file {}: {}".format(data_file, e))
                                continue
                        
                        else:
                            # Handle text files (existing code)
                            with zip_ref.open(data_file) as f:
                                content = f.read().decode('utf-8', errors='ignore')
                                lines = content.split('\n')
                                
                                if len(lines) < 2:
                                    continue
                                
                                # Auto-detect delimiter
                                first_line = lines[0]
                                delimiter = '\t' if '\t' in first_line else (',' if ',' in first_line else None)
                                
                                if not delimiter:
                                    continue
                                
                                # Parse header
                                headers = [h.strip().lower() for h in first_line.split(delimiter)]
                                logger.info("Headers found: {}".format(headers[:10]))
                                
                                # Find relevant columns
                                lncrna_col = None
                                target_col = None
                                regulation_col = None
                                
                                for i, header in enumerate(headers):
                                    if any(keyword in header for keyword in ['lncrna', 'longnoncoding', 'rna_name', 'regulator']):
                                        lncrna_col = i
                                    elif any(keyword in header for keyword in ['target', 'gene', 'symbol']) and 'lncrna' not in header:
                                        target_col = i
                                    elif any(keyword in header for keyword in ['regulation', 'effect', 'direction']):
                                        regulation_col = i
                                
                                if lncrna_col is None or target_col is None:
                                    logger.info("Could not find required columns in {}".format(data_file))
                                    continue
                                
                                logger.info("Using columns - lncRNA: {}, target: {}, regulation: {}".format(
                                    lncrna_col, target_col, regulation_col))
                                
                                # Process data lines
                                for line_num, line in enumerate(tqdm(lines[1:], desc="Parsing {}".format(data_file))):
                                    if not line.strip():
                                        continue
                                    
                                    if max_regs > 0 and len(regulations) >= max_regs:
                                        break
                                    
                                    try:
                                        fields = line.strip().split(delimiter)
                                        if len(fields) <= max(lncrna_col, target_col):
                                            continue
                                        
                                        lncrna_symbol = fields[lncrna_col].strip()
                                        target_symbol = fields[target_col].strip()
                                        
                                        if not lncrna_symbol or not target_symbol:
                                            continue
                                        
                                        # Determine regulation type
                                        regulation_type = 'unknown'
                                        if regulation_col is not None and regulation_col < len(fields):
                                            reg_text = fields[regulation_col].lower()
                                            if any(word in reg_text for word in ['positive', 'activate', 'enhance', 'promote', 'up']):
                                                regulation_type = 'positive'
                                            elif any(word in reg_text for word in ['negative', 'inhibit', 'repress', 'suppress', 'down']):
                                                regulation_type = 'negative'
                                        
                                        regulations.append({
                                            'ncrna_symbol': lncrna_symbol,
                                            'target_symbol': target_symbol,
                                            'regulation_type': regulation_type,
                                            'mechanism': 'transcriptional',
                                            'evidence': 'experimental',
                                            'pmids': 'EVLncRNAs3.0',
                                            'confidence_score': confidence_threshold
                                        })
                                        
                                    except Exception as e:
                                        continue
                                
                                # If we got results from this file, break
                                if regulations:
                                    logger.info("Successfully parsed {} regulations from {}".format(
                                        len(regulations), data_file))
                                    break
                                
                    except Exception as e:
                        logger.warning("Error parsing {}: {}".format(data_file, e))
                        continue
        
        except ImportError:
            logger.warning("pandas not available for Excel file reading")
        except Exception as e:
            logger.error("Error parsing EVLncRNAs: {}".format(e))
            return []
        
        logger.info("Parsed {} lncRNA-target regulations from EVLncRNAs3 (threshold: {})".format(
            len(regulations), confidence_threshold))
        return regulations
    
    def fetch_lnctard_regulations(self, confidence_threshold=None):
        """Fetch lncRNA-target regulations from lnctard2.0 data."""
        logger.info("Fetching lncRNA-target regulations from LncTard2.0...")
        
        lnctard_file = Path(self.config['regulation_sources']['lnctard']['lnctard_file'])
        
        if not lnctard_file.exists():
            logger.error("LncTard file not found: {}".format(lnctard_file))
            return []
        
        if confidence_threshold is None:
            confidence_threshold = self.thresholds.get('lnctard_confidence_threshold', 0.7)
        
        regulations = []
        max_regs = self.limits.get('max_regulations_per_source', 100000)
        
        try:
            # Handle both ZIP and text files
            if lnctard_file.suffix == '.zip':
                with zipfile.ZipFile(lnctard_file, 'r') as zip_ref:
                    file_list = zip_ref.namelist()
                    
                    # Find the data file
                    data_files = [f for f in file_list if f.endswith(('.txt', '.csv', '.tsv')) 
                                 and not f.startswith('__MACOSX') and not '._' in f]
                    
                    if not data_files:
                        logger.warning("No data files found in LncTard ZIP")
                        return []
                    
                    data_file = data_files[0]
                    logger.info("Processing LncTard file: {}".format(data_file))
                    
                    with zip_ref.open(data_file) as f:
                        content = f.read().decode('utf-8', errors='ignore')
                        lines = content.split('\n')
                        
            else:
                # Direct text file
                with open(lnctard_file, 'r', encoding='utf-8', errors='ignore') as f:
                    lines = f.readlines()
            
            if not lines:
                return []
            
            # Parse header
            header_line = lines[0].strip()
            delimiter = '\t' if '\t' in header_line else ','
            headers = [h.strip() for h in header_line.split(delimiter)]
            
            logger.info("LncTard headers: {}".format(headers[:10]))
            
            # Expected column positions for lnctard format
            # DiseaseName, Regulator, Target, RegulationDirection, ...
            if len(headers) >= 4:
                regulator_col = 1  # Regulator column
                target_col = 2     # Target column  
                direction_col = 3  # RegulationDirection column
                
                logger.info("Using LncTard columns - Regulator: {}, Target: {}, Direction: {}".format(
                    regulator_col, target_col, direction_col))
                
                # Process data lines
                for line_num, line in enumerate(tqdm(lines[1:], desc="Parsing LncTard")):
                    if not line.strip():
                        continue
                    
                    if max_regs > 0 and len(regulations) >= max_regs:
                        break
                    
                    try:
                        fields = line.strip().split(delimiter)
                        if len(fields) <= max(regulator_col, target_col, direction_col):
                            continue
                        
                        regulator = fields[regulator_col].strip()
                        target = fields[target_col].strip()
                        direction = fields[direction_col].strip() if direction_col < len(fields) else ''
                        
                        if not regulator or not target:
                            continue
                        
                        # Parse regulation direction
                        regulation_type = 'unknown'
                        direction_lower = direction.lower()
                        
                        if any(word in direction_lower for word in ['positive', 'positively', 'activate', 'enhance', 'promote', 'up']):
                            regulation_type = 'positive'
                        elif any(word in direction_lower for word in ['negative', 'negatively', 'inhibit', 'repress', 'suppress', 'down']):
                            regulation_type = 'negative'
                        
                        regulations.append({
                            'ncrna_symbol': regulator,
                            'target_symbol': target,
                            'regulation_type': regulation_type,
                            'mechanism': 'transcriptional',
                            'evidence': 'experimental',
                            'pmids': 'LncTard2.0',
                            'confidence_score': confidence_threshold
                        })
                        
                    except Exception as e:
                        continue
        
        except Exception as e:
            logger.error("Error parsing LncTard: {}".format(e))
            return []
        
        logger.info("Parsed {} lncRNA-target regulations from LncTard2.0 (threshold: {})".format(
            len(regulations), confidence_threshold))
        return regulations
    
    def fetch_all_regulations(self, mirdb_threshold=None, targetscan_threshold=None, evlncrnas_threshold=None, lnctard_threshold=None):
        """Fetch regulations from all local sources."""
        logger.info("Fetching regulations from all local sources...")
        
        all_regulations = []

        try:
            mirdb_regs = self.fetch_mirdb_regulations(score_threshold=mirdb_threshold)
            all_regulations.extend(mirdb_regs)
            logger.info("Added {} miRDB regulations".format(len(mirdb_regs)))
        except Exception as e:
            logger.warning("miRDB failed: {}".format(e))
        
        try:
            targetscan_regs = self.fetch_targetscan_regulations(context_threshold=targetscan_threshold)
            all_regulations.extend(targetscan_regs)
            logger.info("Added {} TargetScan regulations".format(len(targetscan_regs)))
        except Exception as e:
            logger.warning("TargetScan failed: {}".format(e))
        
        try:
            evlncrnas_regs = self.fetch_evlncrnas_regulations(confidence_threshold=evlncrnas_threshold)
            all_regulations.extend(evlncrnas_regs)
            logger.info("Added {} EVLncRNAs regulations".format(len(evlncrnas_regs)))
        except Exception as e:
            logger.warning("EVLncRNAs failed: {}".format(e))
        
        try:
            lnctard_regs = self.fetch_lnctard_regulations(confidence_threshold=lnctard_threshold)
            all_regulations.extend(lnctard_regs)
            logger.info("Added {} LncTard regulations".format(len(lnctard_regs)))
        except Exception as e:
            logger.warning("LncTard failed: {}".format(e))
                
        logger.info("Total regulations from all sources: {}".format(len(all_regulations)))
        return all_regulations
    
    def fetch_hpo_terms(self):
        """Fetch HPO terms from local OBO file."""
        logger.info("Fetching HPO terms from local OBO file...")
        
        obo_file = Path(self.config['data_sources']['hpo']['hpo_obo'])
        
        if not obo_file.exists():
            logger.error("HPO OBO file not found: {}".format(obo_file))
            return []
        
        hpo_terms = []
        current_term = {}
        max_hpo = self.limits.get('max_hpo_terms', 10000)
        
        try:
            with open(obo_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    if max_hpo > 0 and len(hpo_terms) >= max_hpo:
                        break
                    
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
            logger.error("Error parsing HPO OBO file: {}".format(e))
            return []
        
        logger.info("Parsed {} HPO terms from local file".format(len(hpo_terms)))
        return hpo_terms
    
    def fetch_omim_disease_associations(self):
        """Fetch disease associations from OMIM parser."""
        logger.info("Fetching disease associations from OMIM data...")
        
        try:
            import sys
            from pathlib import Path

            script_dir = Path(__file__).parent
            if str(script_dir) not in sys.path:
                sys.path.append(str(script_dir))
            
            from parse_omim import OMIMDataParser
            
            parser = OMIMDataParser(config_path=self.config_path)
            
            # Get parsed data from OMIM parser
            omim_entries, gene_omim_mapping = parser.parse_omim_file()
            
            conn = sqlite3.connect(self.db_file)
            try:
                # Insert OMIM entries 
                for entry in omim_entries:
                    conn.execute("""
                        INSERT OR REPLACE INTO omim_entries 
                        (omim_id, phenotype, title, description, location, inheritance_pattern, 
                         phenotype_MIM_number, clinical_features, molecular_genetics)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        entry['omim_id'],
                        entry['phenotype'],
                        entry['title'],
                        entry['description'],
                        entry['location'],
                        entry['inheritance_pattern'],
                        entry['phenotype_MIM_number'],
                        entry['clinical_features'],
                        entry['molecular_genetics']
                    ))
                conn.commit()
                logger.info("Inserted {} OMIM entries".format(len(omim_entries)))
            finally:
                conn.close()
            
            disease_associations = []
            
            conn = sqlite3.connect(self.db_file)
            try:
                ensembl_to_symbol = {}
                cursor = conn.execute("SELECT ensembl_gene_id, gene_symbol FROM genes")
                for row in cursor.fetchall():
                    ensembl_to_symbol[row[0]] = row[1]
                
                for mapping in gene_omim_mapping:
                    ensembl_id = mapping.get('ensembl_id')
                    gene_symbol = ensembl_to_symbol.get(ensembl_id)
                    
                    if gene_symbol:  
                        disease_associations.append({
                            'gene_symbol': gene_symbol,
                            'omim_id': mapping.get('omim_id'),
                            'association_type': 'causative',
                            'evidence_level': 'database',
                            'pubmed_ids': 'OMIM'
                        })
            finally:
                conn.close()
            
            logger.info("Converted {} disease associations from OMIM".format(len(disease_associations)))
            return disease_associations
            
        except Exception as e:
            logger.error("Error loading OMIM disease associations: {}".format(e))
            return []

    def _fast_batch_populate_regulations(self, conn, regulations):
        """ Batch processing for regulations using lookup dictionaries"""
        logger.info("Processing {} regulation relationships with fast batch matching...".format(len(regulations)))
        
        # Group regulations by ncRNA and target for batch processing
        regulation_batches = defaultdict(list)
        for reg in regulations:
            key = (reg['ncrna_symbol'], reg['target_symbol'])
            regulation_batches[key].append(reg)
        
        successful_regs = 0
        auto_created_ncrnas = 0
        auto_created_genes = 0
        failed_matches = 0
        
        # Process in batches
        batch_size = 1000
        batch_data = []
        
        for (ncrna_symbol, target_symbol), reg_list in tqdm(regulation_batches.items(), desc="Processing regulation pairs"):
            # Fast lookup for ncRNA
            ncrna_id = None
            for variation in [ncrna_symbol, ncrna_symbol.upper(), ncrna_symbol.lower()]:
                if variation in self.ncrna_lookup:
                    ncrna_id = self.ncrna_lookup[variation]
                    break
            
            # Create ncRNA if not found
            if not ncrna_id:
                try:
                    # Determine ncRNA type from symbol
                    if 'miR' in ncrna_symbol or 'mir' in ncrna_symbol:
                        ncrna_type = 'miRNA'
                    elif any(keyword in ncrna_symbol.upper() for keyword in ['LINC', 'LOC', 'ENSG']):
                        ncrna_type = 'lncRNA'
                    else:
                        # Default based on regulation source
                        reg_source = reg_list[0].get('pmids', '')
                        if 'miRDB' in reg_source or 'TargetScan' in reg_source:
                            ncrna_type = 'miRNA'
                        else:
                            ncrna_type = 'lncRNA'
                    
                    conn.execute("""
                        INSERT OR IGNORE INTO ncrna_genes 
                        (ensembl_gene_id, gene_symbol, ncrna_type, ncrna_subtype, description) 
                        VALUES (?, ?, ?, ?, ?)
                    """, (
                        "AUTO_{}".format(ncrna_symbol), ncrna_symbol, 
                        ncrna_type, ncrna_type, "Auto-added {}".format(ncrna_type)
                    ))
                    
                    result = conn.execute("SELECT ncrna_id FROM ncrna_genes WHERE gene_symbol = ?", (ncrna_symbol,)).fetchone()
                    if result:
                        ncrna_id = result[0]
                        self.ncrna_lookup[ncrna_symbol] = ncrna_id
                        auto_created_ncrnas += 1
                except:
                    failed_matches += len(reg_list)
                    continue
            
            # Fast lookup for gene
            gene_id = None
            variations = [target_symbol, target_symbol.upper(), target_symbol.lower()]
            # Add common variations
            variations.extend([
                target_symbol.replace('-', ''),
                target_symbol.replace('_', ''),
                target_symbol.replace(' ', '')
            ])
            
            for variation in variations:
                if variation in self.gene_lookup:
                    gene_id = self.gene_lookup[variation]
                    break
            
            # Create gene if not found
            if not gene_id:
                try:
                    conn.execute("""
                        INSERT OR IGNORE INTO genes 
                        (ensembl_gene_id, gene_symbol, gene_type, description) 
                        VALUES (?, ?, ?, ?)
                    """, (
                        "AUTO_{}".format(target_symbol), target_symbol,
                        'protein_coding', "Auto-added gene"
                    ))
                    
                    result = conn.execute("SELECT gene_id FROM genes WHERE gene_symbol = ?", (target_symbol,)).fetchone()
                    if result:
                        gene_id = result[0]
                        self.gene_lookup[target_symbol] = gene_id
                        auto_created_genes += 1
                except:
                    failed_matches += len(reg_list)
                    continue
            
            # Batch insert regulations
            if ncrna_id and gene_id:
                for reg in reg_list:
                    batch_data.append((
                        ncrna_id, gene_id, reg['regulation_type'],
                        reg['mechanism'], reg['evidence'], reg['pmids'], 
                        reg.get('confidence_score', 0.8)
                    ))
                    
                    if len(batch_data) >= batch_size:
                        conn.executemany("""
                            INSERT OR REPLACE INTO regulations 
                            (ncrna_id, target_gene_id, regulation_type, regulation_mechanism, 
                             evidence_level, pubmed_ids, confidence_score) 
                            VALUES (?, ?, ?, ?, ?, ?, ?)
                        """, batch_data)
                        successful_regs += len(batch_data)
                        batch_data = []
        
        # Insert remaining batch
        if batch_data:
            conn.executemany("""
                INSERT OR REPLACE INTO regulations 
                (ncrna_id, target_gene_id, regulation_type, regulation_mechanism, 
                 evidence_level, pubmed_ids, confidence_score) 
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, batch_data)
            successful_regs += len(batch_data)
        
        logger.info("Fast regulation processing completed:")
        logger.info("  Successfully inserted: {}".format(successful_regs))
        logger.info("  Auto-created ncRNAs: {}".format(auto_created_ncrnas))
        logger.info("  Auto-created genes: {}".format(auto_created_genes))
        logger.info("  Failed matches: {}".format(failed_matches))
        
        return successful_regs
    
    def populate_database(self, genes=None, ncrnas=None, regulations=None, hpo_terms=None, disease_associations=None):
        """Populate database with optimized batch processing."""
        logger.info("Populating database with optimized batch processing...")
    
        conn = sqlite3.connect(self.db_file)
        conn.execute("PRAGMA foreign_keys = ON")
        # Speed up inserts
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute("PRAGMA journal_mode = MEMORY")
    
        try:
            # Populate basic data first
            if genes:
                logger.info("Inserting {} protein-coding genes...".format(len(genes)))
                gene_columns = list(genes[0].keys())
                placeholders = ','.join(['?' for _ in gene_columns])
            
                conn.executemany(
                    "INSERT OR REPLACE INTO genes ({}) VALUES ({})".format(','.join(gene_columns), placeholders),
                    [tuple(gene[col] for col in gene_columns) for gene in genes]
                )
        
            if ncrnas:
                logger.info("Inserting {} ncRNA genes...".format(len(ncrnas)))
                ncrna_columns = list(ncrnas[0].keys())
                placeholders = ','.join(['?' for _ in ncrna_columns])
            
                conn.executemany(
                    "INSERT OR REPLACE INTO ncrna_genes ({}) VALUES ({})".format(','.join(ncrna_columns), placeholders),
                    [tuple(ncrna[col] for col in ncrna_columns) for ncrna in ncrnas]
                )
        
            if hpo_terms:
                logger.info("Inserting {} HPO terms...".format(len(hpo_terms)))
                conn.executemany(
                    "INSERT OR REPLACE INTO hpo_terms (hpo_id, hpo_name, definition, parent_terms) VALUES (?, ?, ?, ?)",
                    [(term.get('hpo_id'), term.get('hpo_name'), term.get('definition'), term.get('parent_terms'))
                     for term in hpo_terms]
                )
            
            # Commit basic data and initialize lookup
            conn.commit()
            self._initialize_gene_mapper()
    
            # Process regulations with fast batch processing
            if regulations:
                self._fast_batch_populate_regulations(conn, regulations)
            
            # Process disease associations
            if disease_associations:
                logger.info("Inserting {} disease associations...".format(len(disease_associations)))
                successful_diseases = 0
                
                for assoc in disease_associations:
                    # Fast lookup
                    gene_id = None
                    for variation in [assoc['gene_symbol'], assoc['gene_symbol'].upper(), assoc['gene_symbol'].lower()]:
                        if variation in self.gene_lookup:
                            gene_id = self.gene_lookup[variation]
                            break
                    
                    if gene_id:
                        try:
                            conn.execute("""
                                INSERT OR REPLACE INTO gene_disease_associations 
                                (gene_id, hpo_id, omim_id, association_type, evidence_level, pubmed_ids) 
                                VALUES (?, ?, ?, ?, ?, ?)
                            """, (
                                gene_id, assoc.get('hpo_id'), assoc.get('omim_id'),
                                assoc['association_type'], assoc['evidence_level'], assoc['pubmed_ids']
                            ))
                            successful_diseases += 1
                        except Exception as e:
                            continue
                
                logger.info("Successfully inserted {} disease associations".format(successful_diseases))
            
            conn.commit()
            logger.info("Database populated successfully with optimized processing!")
            
        except Exception as e:
            logger.error("Error populating database: {}".format(e))
            conn.rollback()
            raise
        finally:
            conn.close()

def main():
    parser = argparse.ArgumentParser(description='Optimized data fetcher for ncRNA regulation database')
    parser.add_argument('--source', choices=['all', 'genes', 'ncrnas', 'regulations', 'hpo', 'diseases'], 
                       default='all', help='Data source to fetch')
    parser.add_argument('--limit', type=int, help='Limit number of genes to fetch (for testing)')
    parser.add_argument('--mirdb-threshold', type=float, help='miRDB score threshold (default from config)')
    parser.add_argument('--targetscan-threshold', type=float, help='TargetScan context score threshold (default from config)')
    parser.add_argument('--evlncrnas-threshold', type=float, help='EVLncRNAs confidence threshold (default from config)')
    parser.add_argument('--lnctard-threshold', type=float, help='LncTard confidence threshold (default from config)')
    
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
            ncrnas = fetcher.fetch_gencode_lncrnas_combined()
        
        if args.source in ['all', 'regulations']:
            regulations = fetcher.fetch_all_regulations(
                mirdb_threshold=args.mirdb_threshold,
                targetscan_threshold=args.targetscan_threshold,
                evlncrnas_threshold=args.evlncrnas_threshold,
                lnctard_threshold=args.lnctard_threshold
            )
        
        if args.source in ['all', 'hpo']:
            hpo_terms = fetcher.fetch_hpo_terms()
            
        if args.source in ['all', 'diseases']:
            disease_associations = fetcher.fetch_omim_disease_associations()
        

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
        

        print("  Protein-coding genes: {}".format(gene_count))
        print("  ncRNA genes: {}".format(ncrna_count))
        print("  Regulations: {}".format(regulation_count))
        print("  HPO terms: {}".format(hpo_count))
        print("  OMIM entries: {}".format(omim_count))
        print("  Disease associations: {}".format(disease_count))
        
        conn.close()
        
    except Exception as e:
        logger.error("Data fetching failed: {}".format(e))
        raise

if __name__ == "__main__":
    main()
