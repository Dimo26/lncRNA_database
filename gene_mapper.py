"""
 Gene symbol mapper with gene2refseq integration and mappings
"""

import sqlite3
import re
import logging
import gzip
import zipfile
import json
import pickle
from pathlib import Path
import yaml
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PersistentGeneMapper:
    """
    Enhanced persistent gene symbol mapper with gene2refseq integration
    """
    
    def __init__(self, db_path, config_path="config.yaml", mapping_dir="mappings"):
        """Initialize the mapper"""
        self.db_path = Path(db_path)
        self.config_path = config_path
        self.mapping_dir = Path(mapping_dir)
        self.mapping_dir.mkdir(exist_ok=True)
        
        # Mapping files
        self.gene_mapping_file = self.mapping_dir / "gene_symbol_mappings.json"
        self.ncrna_mapping_file = self.mapping_dir / "ncrna_symbol_mappings.json"
        self.variations_file = self.mapping_dir / "symbol_variations.json"
        self.metadata_file = self.mapping_dir / "mapping_metadata.json"
        
        # Load configuration
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # Initialize mapping dictionaries
        self.gene_symbol_map = {}      
        self.ncrna_symbol_map = {}     
        self.symbol_variations = {}    
        self.refseq_mappings = {}
    
    def build_and_save_mappings(self, force_rebuild=False):
        """Build comprehensive mappings and save them to files"""
        
        if not force_rebuild and self._mappings_exist() and self._mappings_are_recent():
            logger.info("Recent mappings found. Use force_rebuild=True to rebuild.")
            return
        
        logger.info("Building comprehensive gene symbol mappings...")
        start_time = time.time()
        
        # Clear existing mappings
        self.gene_symbol_map = {}
        self.ncrna_symbol_map = {}
        self.symbol_variations = {}
        self.refseq_mappings = {}
        
        # Build mappings from all sources
        self._build_from_database()
        self._build_from_gene_info()
        self._build_from_gene2refseq()
        
        # Save mappings to files
        self._save_mappings()
        
        elapsed_time = time.time() - start_time
        logger.info("Mapping build completed in {:.2f} seconds".format(elapsed_time))
        
        # Show statistics
        stats = self.get_mapping_stats()
        logger.info("Built {:,} gene mappings and {:,} ncRNA mappings".format(
            stats['total_gene_mappings'], stats['total_ncrna_mappings']))
        logger.info("Added {:,} RefSeq mappings".format(len(self.refseq_mappings)))
    
    def load_mappings(self):
        """Load pre-built mappings from files"""
        if not self._mappings_exist():
            logger.warning("No pre-built mappings found. Building new mappings...")
            self.build_and_save_mappings()
            return
        
        logger.info("Loading pre-built gene symbol mappings...")
        start_time = time.time()
        
        try:
            # Load mappings
            with open(self.gene_mapping_file, 'r') as f:
                self.gene_symbol_map = json.load(f)
            
            with open(self.ncrna_mapping_file, 'r') as f:
                self.ncrna_symbol_map = json.load(f)
            
            with open(self.variations_file, 'r') as f:
                variations_data = json.load(f)
                # Convert back to sets
                self.symbol_variations = {k: set(v) for k, v in variations_data.items()}
            
            # Load RefSeq mappings if available
            refseq_file = self.mapping_dir / "refseq_mappings.json"
            if refseq_file.exists():
                with open(refseq_file, 'r') as f:
                    self.refseq_mappings = json.load(f)
            
            elapsed_time = time.time() - start_time
            logger.info("Mappings loaded in {:.3f} seconds".format(elapsed_time))
            
            # Show statistics
            stats = self.get_mapping_stats()
            logger.info("Loaded {:,} gene mappings and {:,} ncRNA mappings".format(
                stats['total_gene_mappings'], stats['total_ncrna_mappings']))
            
        except Exception as e:
            logger.error("Error loading mappings: {}. Rebuilding...".format(e))
            self.build_and_save_mappings()
    
    def _mappings_exist(self):
        """Check if mapping files exist"""
        return (self.gene_mapping_file.exists() and 
                self.ncrna_mapping_file.exists() and 
                self.variations_file.exists())
    
    def _mappings_are_recent(self, max_age_hours=24):
        """Check if mappings are recent"""
        if not self.metadata_file.exists():
            return False
        
        try:
            with open(self.metadata_file, 'r') as f:
                metadata = json.load(f)
            
            build_time = metadata.get('build_timestamp', 0)
            current_time = time.time()
            age_hours = (current_time - build_time) / 3600
            
            return age_hours < max_age_hours
            
        except Exception:
            return False
    
    def _save_mappings(self):
        """Save mappings to files"""
        logger.info("Saving mappings to files...")
        
        # Save gene mappings
        with open(self.gene_mapping_file, 'w') as f:
            json.dump(self.gene_symbol_map, f, indent=2)
        
        # Save ncRNA mappings  
        with open(self.ncrna_mapping_file, 'w') as f:
            json.dump(self.ncrna_symbol_map, f, indent=2)
        
        # Save variations (convert sets to lists for JSON)
        variations_data = {k: list(v) for k, v in self.symbol_variations.items()}
        with open(self.variations_file, 'w') as f:
            json.dump(variations_data, f, indent=2)
        
        # Save RefSeq mappings
        refseq_file = self.mapping_dir / "refseq_mappings.json"
        with open(refseq_file, 'w') as f:
            json.dump(self.refseq_mappings, f, indent=2)
        
        # Save metadata
        metadata = {
            'build_timestamp': time.time(),
            'total_gene_mappings': len(self.gene_symbol_map),
            'total_ncrna_mappings': len(self.ncrna_symbol_map),
            'unique_symbols': len(self.symbol_variations),
            'refseq_mappings': len(self.refseq_mappings),
            'config_file': self.config_path,
            'database_file': str(self.db_path)
        }
        
        with open(self.metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info("Mappings saved to {}".format(self.mapping_dir))
    
    def _build_from_database(self):
        """Build mappings from existing database entries"""
        if not self.db_path.exists():
            logger.warning("Database not found: {}".format(self.db_path))
            return
            
        conn = sqlite3.connect(str(self.db_path))
        
        try:
            # Get genes
            cursor = conn.execute("SELECT ensembl_gene_id, gene_symbol FROM genes")
            gene_count = 0
            for ensembl_id, gene_symbol in cursor.fetchall():
                if ensembl_id and gene_symbol:
                    self._add_gene_mapping(ensembl_id, gene_symbol)
                    gene_count += 1
            
            logger.info("Processed {} genes from database".format(gene_count))
            
            # Get ncRNAs
            cursor = conn.execute("SELECT ensembl_gene_id, gene_symbol, ncrna_type FROM ncrna_genes")
            ncrna_count = 0
            for ensembl_id, gene_symbol, ncrna_type in cursor.fetchall():
                if ensembl_id and gene_symbol:
                    self._add_ncrna_mapping(ensembl_id, gene_symbol, ncrna_type)
                    ncrna_count += 1
            
            logger.info("Processed {} ncRNAs from database".format(ncrna_count))
                    
        except Exception as e:
            logger.error("Error reading database: {}".format(e))
        finally:
            conn.close()
    
    def _build_from_gene_info(self):
        """Build comprehensive mappings from NCBI gene_info.gz file"""
        gene_info_file = Path(self.config['data_sources']['gene_info'])
        
        if not gene_info_file.exists():
            logger.warning("Gene info file not found: {}".format(gene_info_file))
            return
        
        logger.info("Processing NCBI gene_info file: {}".format(gene_info_file))
        
        try:
            gene_count = 0
            with gzip.open(str(gene_info_file), 'rt', encoding='utf-8') as f:
                header = f.readline()
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        fields = line.split('\t')
                        if len(fields) >= 15:
                            tax_id = fields[0]
                            
                            if tax_id != '9606':
                                continue
                            
                            gene_id = fields[1]
                            symbol = fields[2]
                            synonyms = fields[4]
                            dbxrefs = fields[5]
                            
                            if symbol and symbol != '-':
                                # Add primary mapping
                                self._add_gene_mapping(gene_id, symbol)
                                
                                # Add synonyms
                                if synonyms and synonyms != '-':
                                    for synonym in synonyms.split('|'):
                                        synonym = synonym.strip()
                                        if synonym:
                                            self._add_gene_mapping(synonym, symbol)
                                
                                # Add cross-references
                                if dbxrefs and dbxrefs != '-':
                                    for dbxref in dbxrefs.split('|'):
                                        dbxref = dbxref.strip()
                                        if dbxref.startswith(('RefSeq:', 'Ensembl:')):
                                            xref_id = dbxref.split(':', 1)[1]
                                            self._add_gene_mapping(xref_id, symbol)
                                            if '.' in xref_id:
                                                self._add_gene_mapping(xref_id.split('.')[0], symbol)
                                
                                gene_count += 1
                                
                    except Exception as e:
                        continue
            
            logger.info("Processed {} genes from gene_info".format(gene_count))
                        
        except Exception as e:
            logger.error("Error reading gene_info file: {}".format(e))
    
    def _build_from_gene2refseq(self):
        """Build mappings from NCBI gene2refseq.gz file"""
        gene2refseq_file = Path(self.config.get('data_sources', {}).get('gene2refseq', 
            '/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/Database_sources/gene2refseq.gz'))
        
        if not gene2refseq_file.exists():
            logger.warning("gene2refseq file not found: {}".format(gene2refseq_file))
            return
        
        logger.info("Processing NCBI gene2refseq file: {}".format(gene2refseq_file))
        
        try:
            refseq_count = 0
            with gzip.open(str(gene2refseq_file), 'rt', encoding='utf-8') as f:
                header = f.readline()
                
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    try:
                        fields = line.split('\t')
                        if len(fields) >= 16:
                            tax_id = fields[0]
                            
                            if tax_id != '9606':
                                continue
                            
                            gene_id = fields[1]
                            status = fields[2]
                            rna_acc = fields[3]
                            protein_acc = fields[5]
                            genomic_acc = fields[7]
                            symbol = fields[15] if len(fields) > 15 else ''
                            
                            # Only use validated RefSeq entries
                            if status not in ['VALIDATED', 'PROVISIONAL', 'REVIEWED']:
                                continue
                            
                            if symbol and symbol != '-':
                                # Map RefSeq accessions to gene symbol
                                if rna_acc and rna_acc != '-':
                                    self._add_gene_mapping(rna_acc, symbol)
                                    if '.' in rna_acc:
                                        self._add_gene_mapping(rna_acc.split('.')[0], symbol)
                                    self.refseq_mappings[rna_acc] = symbol
                                
                                if protein_acc and protein_acc != '-':
                                    self._add_gene_mapping(protein_acc, symbol)
                                    if '.' in protein_acc:
                                        self._add_gene_mapping(protein_acc.split('.')[0], symbol)
                                    self.refseq_mappings[protein_acc] = symbol
                                
                                if genomic_acc and genomic_acc != '-':
                                    self._add_gene_mapping(genomic_acc, symbol)
                                    if '.' in genomic_acc:
                                        self._add_gene_mapping(genomic_acc.split('.')[0], symbol)
                                    self.refseq_mappings[genomic_acc] = symbol
                                
                                # Map GeneID to symbol
                                if gene_id:
                                    self._add_gene_mapping(gene_id, symbol)
                                
                                refseq_count += 1
                                
                    except Exception as e:
                        continue
            
            logger.info("Processed {} RefSeq entries from gene2refseq".format(refseq_count))
                        
        except Exception as e:
            logger.error("Error reading gene2refseq file: {}".format(e))
    
    def _add_gene_mapping(self, identifier, gene_symbol):
        """Add gene mapping with all variations"""
        if not identifier or not gene_symbol:
            return
        
        # Standardize the gene symbol
        standard_symbol = self._standardize_gene_symbol(gene_symbol)
        
        # Add all variations including case variations
        variations = [
            identifier,
            identifier.upper(),
            identifier.lower(),
            gene_symbol,
            gene_symbol.upper(), 
            gene_symbol.lower(),
            standard_symbol,
            standard_symbol.upper(),
            standard_symbol.lower()
        ]
        
        # Clean variations (remove version numbers from IDs)
        clean_variations = []
        for var in variations:
            if var and var.startswith(('ENSG', 'ENST', 'NM_', 'XM_', 'NR_', 'XR_', 'NP_', 'XP_')):
                clean_variations.append(var.split('.')[0])
                clean_variations.append(var)
            elif var:
                clean_variations.append(var)
        
        # Add to mapping (use case-insensitive keys but preserve original symbol)
        for var in clean_variations:
            if var:
                self.gene_symbol_map[var] = standard_symbol
                self.gene_symbol_map[var.upper()] = standard_symbol
                self.gene_symbol_map[var.lower()] = standard_symbol
        
        # Track variations
        if standard_symbol not in self.symbol_variations:
            self.symbol_variations[standard_symbol] = set()
        self.symbol_variations[standard_symbol].update(clean_variations)
    
    def _add_ncrna_mapping(self, identifier, gene_symbol, ncrna_type='ncRNA'):
        """Add ncRNA mapping with all variations"""
        if not identifier or not gene_symbol:
            return
        
        # Standardize the ncRNA symbol
        standard_symbol = self._standardize_ncrna_symbol(gene_symbol, ncrna_type)
        
        # Add all variations including case variations
        variations = [
            identifier,
            identifier.upper(),
            identifier.lower(), 
            gene_symbol,
            gene_symbol.upper(),
            gene_symbol.lower(),
            standard_symbol,
            standard_symbol.upper(),
            standard_symbol.lower()
        ]
        
        # Clean variations
        clean_variations = []
        for var in variations:
            if var and var.startswith(('ENSG', 'ENST', 'NONHSAT', 'NONHSAG')):
                clean_variations.append(var.split('.')[0])
                clean_variations.append(var)
            elif var:
                clean_variations.append(var)
        
        # Add to mapping (use case-insensitive keys)
        for var in clean_variations:
            if var:
                self.ncrna_symbol_map[var] = standard_symbol
                self.ncrna_symbol_map[var.upper()] = standard_symbol
                self.ncrna_symbol_map[var.lower()] = standard_symbol
        
        # Track variations
        if standard_symbol not in self.symbol_variations:
            self.symbol_variations[standard_symbol] = set()
        self.symbol_variations[standard_symbol].update(clean_variations)
    
    def _standardize_gene_symbol(self, symbol):
        """Standardize gene symbol"""
        if not symbol:
            return symbol
        
        symbol = symbol.strip()
        
        # Remove version numbers from IDs
        if symbol.startswith(('ENSG', 'ENST', 'NM_', 'XM_', 'NR_', 'XR_', 'NP_', 'XP_')):
            symbol = symbol.split('.')[0]
        
        return symbol.upper()
    
    def _standardize_ncrna_symbol(self, symbol, ncrna_type='ncRNA'):
        """Standardize ncRNA symbol"""
        if not symbol:
            return symbol
        
        symbol = symbol.strip()
        
        # Handle miRNA names
        if ncrna_type == 'miRNA' or 'mir' in symbol.lower():
            symbol = re.sub(r'^hsa-', '', symbol, flags=re.IGNORECASE)
            symbol = re.sub(r'^mir-', 'miR-', symbol, flags=re.IGNORECASE)
        
        # Handle lncRNA names
        if symbol.startswith('lnc-'):
            symbol = symbol.replace('lnc-', '')
        
        # Remove version numbers from IDs
        if symbol.startswith(('ENSG', 'ENST', 'NONHSAT', 'NONHSAG')):
            symbol = symbol.split('.')[0]
        
        return symbol.upper()
    
    def find_gene_symbol(self, query):
        """Find standard gene symbol for a query with enhanced search"""
        if not query:
            return None
            
        query = query.strip()
        
        #  exact matches first (case-insensitive)
        for test_query in [query, query.upper(), query.lower()]:
            if test_query in self.gene_symbol_map:
                return self.gene_symbol_map[test_query]
        
        #  partial matches
        query_upper = query.upper()
        for key, value in self.gene_symbol_map.items():
            key_upper = key.upper()
            if query_upper == key_upper or (len(query) > 2 and query_upper in key_upper):
                return value
        
        return None
    
    def find_ncrna_symbol(self, query):
        """ standard ncRNA symbol for a query with enhanced search"""
        if not query:
            return None
            
        query = query.strip()
        
        # Try exact matches first (case-insensitive)
        for test_query in [query, query.upper(), query.lower()]:
            if test_query in self.ncrna_symbol_map:
                return self.ncrna_symbol_map[test_query]
        
        #  partial matches
        query_upper = query.upper()
        for key, value in self.ncrna_symbol_map.items():
            key_upper = key.upper()
            if query_upper == key_upper or (len(query) > 2 and query_upper in key_upper):
                return value
        
        return None
    
    def get_all_variations(self, symbol):
        """ known variations of a symbol"""
        standard_symbol = self.find_gene_symbol(symbol) or self.find_ncrna_symbol(symbol)
        
        if standard_symbol and standard_symbol in self.symbol_variations:
            return list(self.symbol_variations[standard_symbol])
        
        return [symbol]  # Return original if no variations found
    
    def get_mapping_stats(self):
        """ comprehensive mapping statistics"""
        stats = {
            'total_gene_mappings': len(self.gene_symbol_map),
            'total_ncrna_mappings': len(self.ncrna_symbol_map),
            'unique_symbols': len(self.symbol_variations),
            'total_variations': sum(len(variations) for variations in self.symbol_variations.values()),
            'refseq_mappings': len(self.refseq_mappings)
        }
        
        return stats
    
    def show_status(self):
        """ current mapping status"""
        print("=" * 50)
        print("GENE SYMBOL MAPPING STATUS")
        print("=" * 50)
        
        if self._mappings_exist():
            try:
                with open(self.metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                build_time = metadata.get('build_timestamp', 0)
                age_hours = (time.time() - build_time) / 3600
              
                print(" Built: {:.1f} hours ago".format(age_hours))
                print(" Gene mappings: {:,}".format(metadata.get('total_gene_mappings', 0)))
                print(" ncRNA mappings: {:,}".format(metadata.get('total_ncrna_mappings', 0)))
                print(" Unique symbols: {:,}".format(metadata.get('unique_symbols', 0)))
                print(" RefSeq mappings: {:,}".format(metadata.get('refseq_mappings', 0)))
                              
            except Exception as e:
                print(" Error reading mapping metadata: {}".format(e))
        else:
            print(" No mappings found")
        
        print("=" * 50)

def main():
    """Main function to build mappings"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Build enhanced gene symbol mappings')
    parser.add_argument('--rebuild', action='store_true', help='Force rebuild even if recent mappings exist')
    parser.add_argument('--status', action='store_true', help='Show mapping status')
    parser.add_argument('--db-path', default='database/ncrna_regulation_db.sqlite', help='Database path')
    parser.add_argument('--config-path', default='config.yaml', help='Config file path')
    
    args = parser.parse_args()
    
    try:
        mapper = PersistentGeneMapper(args.db_path, args.config_path)
        
        if args.status:
            mapper.show_status()
            return
        
        # Build and save mappings
        mapper.build_and_save_mappings(force_rebuild=args.rebuild)
        
    except Exception as e:
        print(" Error: {}".format(e))
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
