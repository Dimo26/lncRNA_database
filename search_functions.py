"""
For search and  knockout effect predictions
"""

import sqlite3
import yaml
import pandas as pd
import argparse
import logging
from pathlib import Path
from typing import List, Dict

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NCRNADatabase:
    def __init__(self, config_path="/proj/sens2023005/nobackup/wharf/dimam/dimam-sens2023005/lncRNA_DB_codes/config.yaml"):
        """Initialize database connection and configuration."""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        self.db_file = Path(self.config['database']['path']) / self.config['database']['name']
        
        if not self.db_file.exists():
            raise FileNotFoundError(f"Database not found: {self.db_file}. Run database_setup.py first.")
    
    def _execute_query(self, query: str, params: tuple = ()) -> List[Dict]:
        """Execute SQL query, returns results as list of dictionaries."""
        try:
            conn = sqlite3.connect(self.db_file)
            conn.row_factory = sqlite3.Row  # Enable dictionary-like access
            
            cursor = conn.execute(query, params)
            results = [dict(row) for row in cursor.fetchall()]
            
            conn.close()
            return results
            
        except Exception as e:
            logger.error(f"Database query error: {e}")
            raise
    
    def search_gene_regulators(self, gene_identifier: str, 
                              id_type: str = "auto") -> pd.DataFrame:
        """
        Search for ncRNAs that regulate a specific protein-coding gene.
        """
        logger.info(f"Searching for ncRNA regulators of gene: {gene_identifier}")
        
        # Determine identifier type
        if id_type == "auto":
            if gene_identifier.startswith("ENSG"):
                id_type = "ensembl"
            else:
                id_type = "symbol"
        
        # Build query based on identifier type
        if id_type == "symbol":
            where_clause = "g.gene_symbol = ? OR g.gene_symbol LIKE ?"
            params = (gene_identifier, f"%{gene_identifier}%")
        else:
            where_clause = "g.ensembl_gene_id = ? OR g.ensembl_gene_id LIKE ?"
            params = (gene_identifier, f"{gene_identifier}%")
        
   
        query = f"""
        SELECT 
            g.gene_symbol as target_gene_symbol,
            g.ensembl_gene_id as target_ensembl_id,
            n.gene_symbol as ncrna_symbol,
            n.ensembl_gene_id as ncrna_ensembl_id,
            n.ncrna_type,
            n.ncrna_subtype,
            r.regulation_type,
            r.regulation_mechanism,
            r.evidence_level,
            r.confidence_score,
            CASE 
                WHEN r.regulation_type = 'positive' THEN 'ncRNA knockout → DOWN-regulation of target gene'
                WHEN r.regulation_type = 'negative' THEN 'ncRNA knockout → UP-regulation of target gene'
                WHEN r.regulation_type = 'bidirectional' THEN 'ncRNA knockout → Complex effects on target gene'
                ELSE 'ncRNA knockout → Unknown effect on target gene'
            END as knockout_effect,
            ABS(n.start_position - g.start_position) as genomic_distance,
            CASE 
                WHEN n.chromosome = g.chromosome THEN 'cis'
                ELSE 'trans'
            END as regulation_location,
            GROUP_CONCAT(h.hpo_name, '; ') as hpo_diseases,
            GROUP_CONCAT(h.hpo_id, '; ') as hpo_ids,
            GROUP_CONCAT(o.phenotype, '; ') as omim_diseases,
            GROUP_CONCAT(o.omim_id, '; ') as omim_ids
        FROM regulations r
        JOIN genes g ON r.target_gene_id = g.gene_id
        JOIN ncrna_genes n ON r.ncrna_id = n.ncrna_id
        LEFT JOIN gene_disease_associations gda ON g.gene_id = gda.gene_id
        LEFT JOIN hpo_terms h ON gda.hpo_id = h.hpo_id
        LEFT JOIN omim_entries o ON gda.omim_id = o.omim_id
        WHERE {where_clause}
        GROUP BY r.regulation_id, g.gene_id, n.ncrna_id
        ORDER BY r.confidence_score DESC, genomic_distance ASC
        """
        
        results = self._execute_query(query, params)
        
        if not results:
            logger.info(f"No ncRNA regulators found for gene: {gene_identifier}")
            return pd.DataFrame()
        
        df = pd.DataFrame(results)
        logger.info(f"Found {len(df)} ncRNA regulators for gene: {gene_identifier}")
        
        return df
    
    def search_ncrna_targets(self, ncrna_identifier: str, 
                            id_type: str = "auto") -> pd.DataFrame:
        
        """ Search for genes regulated by a specific ncRNA."""
        
        logger.info(f"Searching for target genes of ncRNA: {ncrna_identifier}")
        
        # Determine identifier type
        if id_type == "auto":
            if ncrna_identifier.startswith("ENSG"):
                id_type = "ensembl"
            else:
                id_type = "symbol"
        
        # Build query based on identifier type
        if id_type == "symbol":
            where_clause = "n.gene_symbol = ? OR n.gene_symbol LIKE ?"
            params = (ncrna_identifier, f"%{ncrna_identifier}%")
        else:
            where_clause = "n.ensembl_gene_id = ? OR n.ensembl_gene_id LIKE ?"
            params = (ncrna_identifier, f"{ncrna_identifier}%")

        query = f"""
        SELECT 
            n.gene_symbol as ncrna_symbol,
            n.ensembl_gene_id as ncrna_ensembl_id,
            n.ncrna_type,
            n.ncrna_subtype,
            n.length as ncrna_length,
            g.gene_symbol as target_gene_symbol,
            g.ensembl_gene_id as target_ensembl_id,
            r.regulation_type,
            r.regulation_mechanism,
            r.evidence_level,
            r.confidence_score,
            CASE 
                WHEN r.regulation_type = 'positive' THEN 'ncRNA knockout → DOWN-regulation of this target'
                WHEN r.regulation_type = 'negative' THEN 'ncRNA knockout → UP-regulation of this target'
                WHEN r.regulation_type = 'bidirectional' THEN 'ncRNA knockout → Complex effects on this target'
                ELSE 'ncRNA knockout → Unknown effect on this target'
            END as knockout_effect,
            ABS(n.start_position - g.start_position) as genomic_distance,
            CASE 
                WHEN n.chromosome = g.chromosome THEN 'cis'
                ELSE 'trans'
            END as regulation_location,
            GROUP_CONCAT(h.hpo_name, '; ') as target_hpo_diseases,
            GROUP_CONCAT(h.hpo_id, '; ') as target_hpo_ids,
            GROUP_CONCAT(o.phenotype, '; ') as target_omim_diseases,
            GROUP_CONCAT(o.omim_id, '; ') as target_omim_ids
        FROM regulations r
        JOIN genes g ON r.target_gene_id = g.gene_id
        JOIN ncrna_genes n ON r.ncrna_id = n.ncrna_id
        LEFT JOIN gene_disease_associations gda ON g.gene_id = gda.gene_id
        LEFT JOIN hpo_terms h ON gda.hpo_id = h.hpo_id
        LEFT JOIN omim_entries o ON gda.omim_id = o.omim_id
        WHERE {where_clause}
        GROUP BY r.regulation_id, g.gene_id, n.ncrna_id
        ORDER BY r.confidence_score DESC, genomic_distance ASC
        """
        
        results = self._execute_query(query, params)
        
        if not results:
            logger.info(f"No target genes found for ncRNA: {ncrna_identifier}")
            return pd.DataFrame()
        
        df = pd.DataFrame(results)
        logger.info(f"Found {len(df)} target genes for ncRNA: {ncrna_identifier}")
        
        return df
    
    def search_disease_associations(self, gene_identifier: str = None, 
                                  disease_term: str = None) -> pd.DataFrame:
        """
        Search  disease associations of genes or ncRNAs.
        """
        logger.info(f"Searching disease associations for: {gene_identifier or disease_term}")
        
        base_query = """
        SELECT 
            COALESCE(g.gene_symbol, n.gene_symbol) as gene_symbol,
            COALESCE(g.ensembl_gene_id, n.ensembl_gene_id) as ensembl_id,
            CASE WHEN g.gene_id IS NOT NULL THEN 'protein_coding' ELSE n.ncrna_type END as gene_type,
            h.hpo_name,
            h.hpo_id,
            o.phenotype as omim_phenotype,
            o.omim_id,
            d.association_type,
            d.evidence_level,
            d.pubmed_ids
        FROM gene_disease_associations d
        LEFT JOIN genes g ON d.gene_id = g.gene_id
        LEFT JOIN ncrna_genes n ON d.ncrna_id = n.ncrna_id
        LEFT JOIN hpo_terms h ON d.hpo_id = h.hpo_id
        LEFT JOIN omim_entries o ON d.omim_id = o.omim_id
        """
        
        conditions = []
        params = []
        
        if gene_identifier:
            conditions.append("(g.gene_symbol LIKE ? OR n.gene_symbol LIKE ? OR g.ensembl_gene_id LIKE ? OR n.ensembl_gene_id LIKE ?)")
            params.extend([f"%{gene_identifier}%"] * 4)
        
        if disease_term:
            conditions.append("(h.hpo_name LIKE ? OR o.phenotype LIKE ?)")
            params.extend([f"%{disease_term}%"] * 2)
        
        if conditions:
            query = base_query + " WHERE " + " AND ".join(conditions)
        else:
            query = base_query
        
        query += " ORDER BY d.association_type, gene_symbol"
        
        results = self._execute_query(query, tuple(params))
        
        if not results:
            logger.info("No disease associations found")
            return pd.DataFrame()
        
        df = pd.DataFrame(results)
        logger.info(f"Found {len(df)} disease associations")
        
        return df
    
    def get_database_statistics(self) -> Dict:
        """Get database statistics."""
        logger.info("Calculating database statistics...")
        
        stats = {}
        
        # Basic counts
        stats['total_genes'] = self._execute_query("SELECT COUNT(*) as count FROM genes")[0]['count']
        stats['total_ncrnas'] = self._execute_query("SELECT COUNT(*) as count FROM ncrna_genes")[0]['count']
        stats['total_regulations'] = self._execute_query("SELECT COUNT(*) as count FROM regulations")[0]['count']
        stats['total_hpo_terms'] = self._execute_query("SELECT COUNT(*) as count FROM hpo_terms")[0]['count']
        stats['total_omim_entries'] = self._execute_query("SELECT COUNT(*) as count FROM omim_entries")[0]['count']
        stats['total_disease_associations'] = self._execute_query("SELECT COUNT(*) as count FROM gene_disease_associations")[0]['count']
        
        # ncRNA types
        ncrna_types = self._execute_query("SELECT ncrna_type, COUNT(*) as count FROM ncrna_genes GROUP BY ncrna_type")
        stats['ncrna_types'] = {row['ncrna_type']: row['count'] for row in ncrna_types}
        
        # Regulation types
        reg_types = self._execute_query("SELECT regulation_type, COUNT(*) as count FROM regulations GROUP BY regulation_type")
        stats['regulation_types'] = {row['regulation_type']: row['count'] for row in reg_types}
        
        return stats
    
    def export_results(self, df: pd.DataFrame, filename: str, format: str = "csv"):
        """Export search results to csv file."""
        if df.empty:
            logger.warning("No data to export")
            return
        
        output_dir = Path("results")
        output_dir.mkdir(exist_ok=True)
        
        if format.lower() == "csv":
            output_file = output_dir / f"{filename}.csv"
            df.to_csv(output_file, index=False)
        elif format.lower() == "json":
            output_file = output_dir / f"{filename}.json"
            df.to_json(output_file, orient="records", indent=2)
        elif format.lower() == "excel":
            output_file = output_dir / f"{filename}.xlsx"
            df.to_excel(output_file, index=False)
        
        logger.info(f"Results exported to: {output_file}")
        return output_file

def main():
    """Command line interface for database search."""
    parser = argparse.ArgumentParser(description='Search ncRNA regulation database')
    parser.add_argument('--gene', help='Gene symbol or Ensembl ID')
    parser.add_argument('--ncrna', help='ncRNA symbol or Ensembl ID')
    parser.add_argument('--disease', help='Disease term (HPO or OMIM)')
    parser.add_argument('--search-type', 
                       choices=['gene-to-ncrna', 'ncrna-to-gene', 'disease', 'cis-trans', 'stats'],
                       default='gene-to-ncrna',
                       help='Type of search to perform')
    parser.add_argument('--export', choices=['csv', 'json', 'excel'],
                       help='Export results to file')
    parser.add_argument('--limit', type=int, help='Limit number of results')
    
    args = parser.parse_args()
    
    try:
        db = NCRNADatabase()
        
        # Perform search based on type
        if args.search_type == 'gene-to-ncrna':
            if not args.gene:
                print(" --gene required for gene-to-ncrna search")
                return
            df = db.search_gene_regulators(args.gene)
            print(f" ncRNA regulators of {args.gene}:")
            
        elif args.search_type == 'ncrna-to-gene':
            if not args.ncrna:
                print(" --ncrna required for ncrna-to-gene search")
                return
            df = db.search_ncrna_targets(args.ncrna)
            print(f" Target genes of {args.ncrna}:")
            
        elif args.search_type == 'disease':
            df = db.search_disease_associations(args.gene, args.disease)
            print(f" Disease associations:")
            
        elif args.search_type == 'stats':
            stats = db.get_database_statistics()
            print("Database Statistics: ")
            print(f" Total genes: {stats['total_genes']}")
            print(f" Total ncRNAs: {stats['total_ncrnas']}")
            print(f" Total regulations: {stats['total_regulations']}")
            print(f" Total HPO terms: {stats['total_hpo_terms']}")
            print(f" Total OMIM entries: {stats['total_omim_entries']}")
            print(f" Total disease associations: {stats['total_disease_associations']}")
            
            print(" ncRNA Types:")
            for ncrna_type, count in stats['ncrna_types'].items():
                print(f"  • {ncrna_type}: {count}")
            
            print(" Regulation Types:")
            for reg_type, count in stats['regulation_types'].items():
                print(f" • {reg_type}: {count}")
            
            return
        
        # Display results
        if not df.empty:
            if args.limit:
                df = df.head(args.limit)
            
            print(f" Found {len(df)} results:")
            print("=" * 80)
            
            # Display key columns
            for _, row in df.iterrows():
                if args.search_type == 'gene-to-ncrna':
                    print(f" ncRNA: {row['ncrna_symbol']} ({row['ncrna_type']})")
                    print(f" Regulation: {row['regulation_type']} ({row.get('regulation_mechanism', 'N/A')})")
                    print(f" Knockout effect: {row['knockout_effect']}")
                    print(f" Evidence: {row['evidence_level']} (score: {row.get('confidence_score', 'N/A')})")
                elif args.search_type == 'ncrna-to-gene':
                    print(f" Target: {row['target_gene_symbol']}")
                    print(f" Regulation: {row['regulation_type']} ({row.get('regulation_mechanism', 'N/A')})")
                    print(f" Knockout effect: {row['knockout_effect']}")
                    print(f" Evidence: {row['evidence_level']} (score: {row.get('confidence_score', 'N/A')})")
                elif args.search_type == 'disease':
                    print(f" Gene: {row['gene_symbol']} ({row['gene_type']})")
                    print(f" Disease: {row.get('hpo_name', 'N/A')} / {row.get('omim_phenotype', 'N/A')}")
                    print(f" Association: {row['association_type']} ({row['evidence_level']})")
                print("-" * 40)
            
            # Export if requested
            if args.export:
                filename = f"{args.search_type}_{args.gene or args.ncrna or 'results'}"
                db.export_results(df, filename, args.export)
        else:
            print(" No results found")
    
    except Exception as e:
        logger.error(f"Search failed: {e}")
        raise

if __name__ == "__main__":
    main()
