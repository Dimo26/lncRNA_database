
"""
 SQLite database with schema for bidirectional seaches for genes and ncRNAs. 
"""

import sqlite3
import yaml
import logging
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DatabaseSetup:
    def __init__(self, config_path="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/config.yaml"):
        """Initialize database setup with configuration."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # Create database directory if it doesn't exist
        db_path = Path(self.config['database']['path'])
        db_path.mkdir(exist_ok=True)
        
        self.db_file = db_path / self.config['database']['name']
        self.conn = None
    
    def create_database(self):
        """Create SQLite database with complete schema."""
        try:
            self.conn = sqlite3.connect(self.db_file)
            self.conn.execute("PRAGMA foreign_keys = ON")
            
            logger.info(f"Creating database: {self.db_file}")
            
            # Create tables
            self._create_gene_tables()
            self._create_regulation_tables()
            self._create_disease_tables()
            self._create_indexes()
            
            self.conn.commit()
            logger.info("Database created successfully!")
            
        except Exception as e:
            logger.error(f"Error creating database: {e}")
            raise
        finally:
            if self.conn:
                self.conn.close()
    
    def _create_gene_tables(self):

        genes_sql = """
        CREATE TABLE IF NOT EXISTS genes (
            gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
            ensembl_gene_id TEXT UNIQUE NOT NULL,
            gene_symbol TEXT NOT NULL,
            gene_name TEXT,
            chromosome TEXT,
            start_position INTEGER,
            end_position INTEGER,
            strand TEXT CHECK(strand IN ('+', '-')),
            gene_type TEXT DEFAULT 'protein_coding',
            description TEXT
        );
        """

        ncrna_genes_sql = """
        CREATE TABLE IF NOT EXISTS ncrna_genes (
            ncrna_id INTEGER PRIMARY KEY AUTOINCREMENT,
            ensembl_gene_id TEXT UNIQUE NOT NULL,
            gene_symbol TEXT,
            gene_name TEXT,
            chromosome TEXT,
            start_position INTEGER,
            end_position INTEGER,
            strand TEXT CHECK(strand IN ('+', '-')),
            ncrna_type TEXT NOT NULL, -- lncRNA, miRNA, etc.
            ncrna_subtype TEXT, -- antisense, intronic, intergenic, etc.
            length INTEGER,
            lncipedia_id TEXT,
            description TEXT
        );
        """
        
        self.conn.execute(genes_sql)
        self.conn.execute(ncrna_genes_sql)
        logger.info("Created gene tables")
    
    def _create_regulation_tables(self):
        
        regulations_sql = """
        CREATE TABLE IF NOT EXISTS regulations (
            regulation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            ncrna_id INTEGER NOT NULL,
            target_gene_id INTEGER NOT NULL,
            regulation_type TEXT NOT NULL CHECK(regulation_type IN (
                'positive', 'negative', 'bidirectional', 'unknown'
            )),
            regulation_mechanism TEXT,
            evidence_level TEXT CHECK(evidence_level IN (
                'experimental', 'computational', 'literature', 'database'
            )),
            evidence_details TEXT,
            interaction_distance INTEGER,
            pubmed_ids TEXT,
            database_source TEXT,
            confidence_score REAL CHECK(confidence_score >= 0 AND confidence_score <= 1),
            FOREIGN KEY (ncrna_id) REFERENCES ncrna_genes (ncrna_id) ON DELETE CASCADE,
            FOREIGN KEY (target_gene_id) REFERENCES genes (gene_id) ON DELETE CASCADE,
            UNIQUE(ncrna_id, target_gene_id, regulation_type)
        );
        """
        
        self.conn.execute(regulations_sql)
        logger.info("Created regulation tables")
    
    def _create_disease_tables(self):
        
        hpo_terms_sql = """
        CREATE TABLE IF NOT EXISTS hpo_terms (
            hpo_id TEXT PRIMARY KEY,
            hpo_name TEXT NOT NULL,
            definition TEXT,
            parent_terms TEXT -- JSON array of parent HPO IDs
        );
        """
        
        omim_entries_sql = """
        CREATE TABLE IF NOT EXISTS omim_entries (
            omim_id TEXT PRIMARY KEY,
            location TEXT,
            inheritance_pattern TEXT,
            phenotype_MIM_number INTEGER, 
            phenotype TEXT,
            title TEXT,
            description TEXT,
            clinical_features TEXT,
            molecular_genetics TEXT
        );
        """

        gene_disease_sql = """
        CREATE TABLE IF NOT EXISTS gene_disease_associations (
            association_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_id INTEGER,
            ncrna_id INTEGER,
            hpo_id TEXT,
            omim_id TEXT,
            association_type TEXT CHECK(association_type IN (
                'causative', 'risk_factor', 'protective', 'modifier'
            )),
            evidence_level TEXT,
            pubmed_ids TEXT,
            FOREIGN KEY (gene_id) REFERENCES genes (gene_id) ON DELETE CASCADE,
            FOREIGN KEY (ncrna_id) REFERENCES ncrna_genes (ncrna_id) ON DELETE CASCADE,
            FOREIGN KEY (hpo_id) REFERENCES hpo_terms (hpo_id),
            FOREIGN KEY (omim_id) REFERENCES omim_entries (omim_id),
            CHECK ((gene_id IS NOT NULL AND ncrna_id IS NULL) OR 
                   (gene_id IS NULL AND ncrna_id IS NOT NULL))
        );
        """
        
        self.conn.execute(hpo_terms_sql)
        self.conn.execute(omim_entries_sql)
        self.conn.execute(gene_disease_sql)
        logger.info("Created disease association tables")
    
    def _create_indexes(self):
        """Create indexes for optimal search performance."""
        indexes = [
            # Gene indexes
            "CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(gene_symbol)",
            "CREATE INDEX IF NOT EXISTS idx_genes_ensembl ON genes(ensembl_gene_id)",
            "CREATE INDEX IF NOT EXISTS idx_genes_location ON genes(chromosome, start_position, end_position)",
            
            # ncRNA indexes
            "CREATE INDEX IF NOT EXISTS idx_ncrna_symbol ON ncrna_genes(gene_symbol)",
            "CREATE INDEX IF NOT EXISTS idx_ncrna_ensembl ON ncrna_genes(ensembl_gene_id)",
            "CREATE INDEX IF NOT EXISTS idx_ncrna_type ON ncrna_genes(ncrna_type)",
            "CREATE INDEX IF NOT EXISTS idx_ncrna_location ON ncrna_genes(chromosome, start_position, end_position)",
            
            # Regulation indexes (critical for performance)
            "CREATE INDEX IF NOT EXISTS idx_regulations_ncrna ON regulations(ncrna_id)",
            "CREATE INDEX IF NOT EXISTS idx_regulations_gene ON regulations(target_gene_id)",
            "CREATE INDEX IF NOT EXISTS idx_regulations_type ON regulations(regulation_type)",
            "CREATE INDEX IF NOT EXISTS idx_regulations_evidence ON regulations(evidence_level)",
            "CREATE INDEX IF NOT EXISTS idx_regulations_composite ON regulations(ncrna_id, target_gene_id, regulation_type)",
            
            # Disease indexes
            "CREATE INDEX IF NOT EXISTS idx_disease_gene ON gene_disease_associations(gene_id)",
            "CREATE INDEX IF NOT EXISTS idx_disease_ncrna ON gene_disease_associations(ncrna_id)",
            "CREATE INDEX IF NOT EXISTS idx_disease_hpo ON gene_disease_associations(hpo_id)",
            "CREATE INDEX IF NOT EXISTS idx_disease_omim ON gene_disease_associations(omim_id)",
        ]
        
        for index_sql in indexes:
            self.conn.execute(index_sql)
        
        logger.info("Created database indexes")
    
    def create_views(self):
        
        # Gene to ncRNA regulation view
        gene_to_ncrna_view = """
        CREATE VIEW IF NOT EXISTS gene_to_ncrna_view AS
        SELECT 
            g.gene_symbol as target_gene,
            g.ensembl_gene_id as target_ensembl_id,
            n.gene_symbol as ncrna_symbol,
            n.ensembl_gene_id as ncrna_ensembl_id,
            n.ncrna_type,
            r.regulation_type,
            r.evidence_level,
            r.confidence_score,
            ABS(n.start_position - g.start_position) as genomic_distance
        FROM regulations r
        JOIN genes g ON r.target_gene_id = g.gene_id
        JOIN ncrna_genes n ON r.ncrna_id = n.ncrna_id;
        """

        ncrna_to_gene_view = """
        CREATE VIEW IF NOT EXISTS ncrna_to_gene_view AS
        SELECT 
            n.gene_symbol as ncrna_symbol,
            n.ensembl_gene_id as ncrna_ensembl_id,
            n.ncrna_type,
            g.gene_symbol as target_gene,
            g.ensembl_gene_id as target_ensembl_id,
            r.regulation_type,
            r.evidence_level,
            r.confidence_score
        FROM regulations r
        JOIN genes g ON r.target_gene_id = g.gene_id
        JOIN ncrna_genes n ON r.ncrna_id = n.ncrna_id;
        """

        disease_view = """
        CREATE VIEW IF NOT EXISTS disease_associations_view AS
        SELECT 
            COALESCE(g.gene_symbol, n.gene_symbol) as gene_symbol,
            COALESCE(g.ensembl_gene_id, n.ensembl_gene_id) as ensembl_id,
            CASE WHEN g.gene_id IS NOT NULL THEN 'protein_coding' ELSE n.ncrna_type END as gene_type,
            h.hpo_name,
            h.hpo_id,
            o.phenotype as omim_phenotype,
            o.title as omim_title,
            o.omim_id,
            d.association_type,
            d.evidence_level,
            d.pubmed_ids
        FROM gene_disease_associations d
        LEFT JOIN genes g ON d.gene_id = g.gene_id
        LEFT JOIN ncrna_genes n ON d.ncrna_id = n.ncrna_id
        LEFT JOIN hpo_terms h ON d.hpo_id = h.hpo_id
        LEFT JOIN omim_entries o ON d.omim_id = o.omim_id;
        """
        
        try:
            self.conn = sqlite3.connect(self.db_file)
            self.conn.execute(gene_to_ncrna_view)
            self.conn.execute(ncrna_to_gene_view)
            self.conn.execute(disease_view)
            self.conn.commit()
            logger.info("Created database views")
        except Exception as e:
            logger.error(f"Error creating views: {e}")
        finally:
            if self.conn:
                self.conn.close()

def main():
    """Main function to set up database."""
    
    # Create logs directory
    Path("logs").mkdir(exist_ok=True)
    
    try:
        db_setup = DatabaseSetup()
        db_setup.create_database()
        db_setup.create_views()

    except Exception as e:
        logger.error(f"Failed to create database: {e}")
        raise

if __name__ == "__main__":
    main()
