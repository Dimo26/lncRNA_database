"""
Parse OMIM data from mart_export_OMIM.txt and populate database.
"""

import sqlite3
import yaml
import logging
import re
from pathlib import Path
from typing import List, Dict, Optional

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class OMIMDataParser:
    def __init__(self, config_path="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/config.yaml"):
        """Initialize OMIM data parser."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_file = Path(self.config['database']['path']) / self.config['database']['name']
        
        
        self.omim_file = Path(self.config['data_sources']['omim_data']['mart_file'])
        
        if not self.db_file.exists():
            raise FileNotFoundError(f"Database not found: {self.db_file}. Run database_setup_final.py first.")
        
        if not self.omim_file.exists():
            raise FileNotFoundError(f"OMIM export file not found: {self.omim_file}")
    
    def parse_omim_file(self) -> List[Dict]:
        """
        Parse the mart_export_OMIM.txt file to extract OMIM entries.
        
        Expected format:
        ENSG00000XXXXX    OMIM_ID    PHENOTYPE_DESCRIPTION
        """
        logger.info(f"Parsing OMIM file: {self.omim_file}")
        
        omim_entries = []
        gene_omim_mapping = []
        processed_omim_ids = set()
        
        try:
            with open(self.omim_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    # Split by tab
                    parts = line.split('\t')
                    
                    # Expected format: ENSEMBL_ID, OMIM_ID, PHENOTYPE
                    if len(parts) >= 2:
                        ensembl_id = parts[0].strip()
                        omim_id = parts[1].strip() if len(parts) > 1 else ""
                        phenotype = parts[2].strip() if len(parts) > 2 else ""
                        
                        # Skip entries without OMIM ID
                        if not omim_id or not omim_id.isdigit():
                            continue
                        
                        # Store gene-OMIM mapping
                        if ensembl_id.startswith('ENSG'):
                            gene_omim_mapping.append({
                                'ensembl_id': ensembl_id,
                                'omim_id': omim_id,
                                'phenotype': phenotype
                            })
                        
                        # Create OMIM entry if not already processed
                        if omim_id not in processed_omim_ids:
                            processed_omim_ids.add(omim_id)
                            
                            # Parse phenotype for additional information
                            title, description = self._parse_phenotype(phenotype)
                            
                            omim_entry = {
                                'omim_id': omim_id,
                                'phenotype': phenotype,
                                'title': title,
                                'description': description,
                                'location': None, 
                                'inheritance_pattern': self._extract_inheritance(phenotype),
                                'phenotype_MIM_number': int(omim_id) if omim_id.isdigit() else None,
                                'clinical_features': self._extract_clinical_features(phenotype),
                                'molecular_genetics': None  
                            }
                            
                            omim_entries.append(omim_entry)
        
        except Exception as e:
            logger.error(f"Error parsing OMIM file at line {line_num}: {e}")
            raise
        
        logger.info(f"Parsed {len(omim_entries)} OMIM entries")
        logger.info(f"Found {len(gene_omim_mapping)} gene-OMIM mappings")
        
        return omim_entries, gene_omim_mapping
    
    def _parse_phenotype(self, phenotype: str) -> tuple:
        """Parse phenotype string to extract title and description."""
        if not phenotype:
            return "", ""
        
        # Remove common patterns
        phenotype = re.sub(r';;.*', '', phenotype)
        phenotype = re.sub(r'\s+', ' ', phenotype).strip() 
        
        # Split by semicolon or comma for title
        parts = re.split(r'[;,]', phenotype)
        title = parts[0].strip() if parts else ""
        
        # Use first 100 characters as description
        description = phenotype[:100] + "..." if len(phenotype) > 100 else phenotype
        
        return title, description
    
    def _extract_inheritance(self, phenotype: str) -> Optional[str]:
        """Extract inheritance pattern from phenotype description."""
        if not phenotype:
            return None
        
        phenotype_lower = phenotype.lower()
        
        if 'autosomal dominant' in phenotype_lower:
            return 'Autosomal dominant'
        elif 'autosomal recessive' in phenotype_lower:
            return 'Autosomal recessive'
        elif 'x-linked' in phenotype_lower or 'x linked' in phenotype_lower:
            return 'X-linked'
        elif 'mitochondrial' in phenotype_lower:
            return 'Mitochondrial'
        elif 'somatic' in phenotype_lower:
            return 'Somatic mutation'
        
        return None
    
    def _extract_clinical_features(self, phenotype: str) -> Optional[str]:
        """Extract clinical features from phenotype description."""
        if not phenotype:
            return None
        
        # Look for specific clinical terms
        clinical_terms = [
            'cancer', 'tumor', 'carcinoma', 'lymphoma', 'leukemia',
            'syndrome', 'disease', 'disorder', 'deficiency',
            'susceptibility', 'resistance', 'immunodeficiency'
        ]
        
        phenotype_lower = phenotype.lower()
        found_terms = [term for term in clinical_terms if term in phenotype_lower]
        
        if found_terms:
            return f"Associated with {', '.join(found_terms)}"
        
        return phenotype[:50] + "..." if len(phenotype) > 50 else phenotype
    
    def populate_database(self, omim_entries: List[Dict], gene_omim_mapping: List[Dict]):
        """Populate database with OMIM entries and gene mappings."""
        logger.info("Populating database with OMIM data...")
        
        conn = sqlite3.connect(self.db_file)
        conn.execute("PRAGMA foreign_keys = ON")
        
        try:
            # Insert OMIM entries
            logger.info(f"Inserting {len(omim_entries)} OMIM entries...")
            
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
            
            # Create gene-disease associations
            logger.info(f"Creating gene-disease associations...")
            
            successful_associations = 0
            for mapping in gene_omim_mapping:
                # Get gene_id from ensembl_id
                gene_result = conn.execute(
                    "SELECT gene_id FROM genes WHERE ensembl_gene_id = ?",
                    (mapping['ensembl_id'],)
                ).fetchone()
                
                if gene_result:
                    gene_id = gene_result[0]
                    
                    # Insert association
                    conn.execute("""
                        INSERT OR IGNORE INTO gene_disease_associations
                        (gene_id, omim_id, association_type, evidence_level, pubmed_ids)
                        VALUES (?, ?, 'causative', 'database', 'OMIM')
                    """, (gene_id, mapping['omim_id']))
                    
                    successful_associations += 1
            
            conn.commit()
            logger.info(f"Successfully created {successful_associations} gene-disease associations")
            
            # Print statistics
            stats = self._get_omim_statistics(conn)
            logger.info(f"Database now contains:")
            logger.info(f"  • OMIM entries: {stats['omim_count']}")
            logger.info(f"  • Disease associations: {stats['association_count']}")
            logger.info(f"  • Genes with disease associations: {stats['genes_with_diseases']}")
        
        except Exception as e:
            logger.error(f"Error populating database: {e}")
            conn.rollback()
            raise
        finally:
            conn.close()
    
    def _get_omim_statistics(self, conn) -> Dict:
        """Get OMIM-related statistics from database."""
        stats = {}
        
        stats['omim_count'] = conn.execute("SELECT COUNT(*) FROM omim_entries").fetchone()[0]
        stats['association_count'] = conn.execute("SELECT COUNT(*) FROM gene_disease_associations").fetchone()[0]
        stats['genes_with_diseases'] = conn.execute("""
            SELECT COUNT(DISTINCT gene_id) FROM gene_disease_associations WHERE gene_id IS NOT NULL
        """).fetchone()[0]
        
        return stats
    
    def run(self):
        """Run the complete OMIM data parsing and population process."""
        try:
            logger.info(" Starting OMIM data parsing...")
            
            # Parse OMIM file
            omim_entries, gene_omim_mapping = self.parse_omim_file()
            
            if not omim_entries:
                logger.warning("No OMIM entries found in the file")
                return
            
            # Populate database
            self.populate_database(omim_entries, gene_omim_mapping)
            
            logger.info(" OMIM data parsing completed successfully!")
            
        except Exception as e:
            logger.error(f" OMIM data parsing failed: {e}")
            raise

def main():
    try:
        parser = OMIMDataParser()
        parser.run()
        
        print(" OMIM data parsing completed!")

    except Exception as e:
        print(f" Error: {e}")
        raise

if __name__ == "__main__":
    main()
