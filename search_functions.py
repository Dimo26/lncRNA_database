"""
Enhanced search functions with improved robustness and debugging
"""
import sqlite3
import pandas as pd
import argparse
import sys
import json
from pathlib import Path
import yaml
import logging
import re

from gene_mapper import PersistentGeneMapper

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EnhancedNCRNADatabase:
    def __init__(self, db_path="database/ncrna_regulation_db.sqlite", config_path="config.yaml"):
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError("Database not found: {}".format(db_path))
        
        self.config = {}
        if Path(config_path).exists():
            with open(config_path, 'r') as f:
                self.config = yaml.safe_load(f)
        
        self.default_thresholds = self.config.get('thresholds', {})
        
        logger.info("Loading gene symbol mappings...")
        self.mapper = PersistentGeneMapper(str(self.db_path), config_path)
        self.mapper.load_mappings()
    
    def get_connection(self):
        return sqlite3.connect(self.db_path)
    
    def search_gene_regulators(self, gene_symbol, 
                             min_confidence=0.0,
                             mirdb_threshold=None,
                             targetscan_threshold=None,
                             include_literature=True,
                             include_disease=True,
                             debug=False):
        
        logger.info("Searching regulators for gene: {}".format(gene_symbol))
        
        if debug:
            self._debug_gene_search(gene_symbol)
        
        results = self._comprehensive_gene_search(
            gene_symbol, min_confidence, mirdb_threshold, 
            targetscan_threshold, include_literature, debug
        )
        
        if not results.empty and include_disease:
            results = self._add_disease_info(results, 'target_gene_symbol')
        
        logger.info("Found {} regulators for {}".format(len(results), gene_symbol))
        return results
    
    def _debug_gene_search(self, gene_symbol):
        """Debug gene search to understand mapping issues"""
        logger.info("=== DEBUGGING GENE SEARCH FOR: {} ===".format(gene_symbol))
        
        # Check if gene exists in mappings
        standard_symbol = self.mapper.find_gene_symbol(gene_symbol)
        logger.info("Standard symbol from mapper: {}".format(standard_symbol))
        
        # Get all variations
        variations = self.mapper.get_all_variations(gene_symbol)
        logger.info("Found {} variations: {}".format(len(variations), variations[:10]))
        
        # Check database for gene existence
        conn = self.get_connection()
        try:
            # Check in genes table
            gene_count = conn.execute(
                "SELECT COUNT(*) FROM genes WHERE gene_symbol = ?", 
                (gene_symbol,)
            ).fetchone()[0]
            logger.info("Direct match in genes table: {}".format(gene_count))
            
            # Check case insensitive
            gene_count_ci = conn.execute(
                "SELECT COUNT(*) FROM genes WHERE UPPER(gene_symbol) = UPPER(?)", 
                (gene_symbol,)
            ).fetchone()[0]
            logger.info("Case-insensitive match in genes table: {}".format(gene_count_ci))
            
            # Check with variations
            if variations:
                placeholders = ','.join(['?' for _ in variations[:5]])
                variation_count = conn.execute(
                    "SELECT COUNT(DISTINCT gene_id) FROM genes WHERE gene_symbol IN ({})".format(placeholders),
                    variations[:5]
                ).fetchone()[0]
                logger.info("Variation matches in genes table: {}".format(variation_count))
            
            # Check regulations table
            reg_count = conn.execute("SELECT COUNT(*) FROM regulations").fetchone()[0]
            logger.info("Total regulations in database: {}".format(reg_count))
            
            # Sample gene symbols in database
            sample_genes = conn.execute(
                "SELECT DISTINCT gene_symbol FROM genes LIMIT 10"
            ).fetchall()
            logger.info("Sample gene symbols in database: {}".format([g[0] for g in sample_genes]))
            
        finally:
            conn.close()
        
        logger.info("=== END DEBUG ===")
    
    def _comprehensive_gene_search(self, gene_symbol, min_confidence, mirdb_threshold, 
                                 targetscan_threshold, include_literature, debug=False):
        
        all_results = []
        search_methods = [
            ("Direct exact match", self._direct_exact_search),
            ("Case insensitive match", self._case_insensitive_search), 
            ("Variation search", self._variation_search),
            ("Pattern search", self._pattern_search)
        ]
        
        for method_name, search_func in search_methods:
            logger.info("Trying search method: {}".format(method_name))
            
            results = search_func(gene_symbol, min_confidence, mirdb_threshold, 
                                targetscan_threshold, include_literature)
            
            if not results.empty:
                results['search_method'] = method_name
                all_results.append(results)
                logger.info("Found {} results using {}".format(len(results), method_name))
                
                if len(all_results) >= 2 or len(results) >= 100:
                    break
            else:
                logger.info("No results from {}".format(method_name))
        
        if all_results:
            combined_results = pd.concat(all_results, ignore_index=True)
            combined_results = combined_results.drop_duplicates(
                subset=['ncrna_symbol', 'target_gene_symbol'], keep='first'
            )
            return combined_results.sort_values('confidence_score', ascending=False)
        
        return pd.DataFrame()
    
    def _direct_exact_search(self, gene_symbol, min_confidence, mirdb_threshold, 
                           targetscan_threshold, include_literature):
        """Direct exact match search"""
        return self._execute_regulation_query(
            "g.gene_symbol = ?", [gene_symbol], min_confidence, 
            mirdb_threshold, targetscan_threshold, include_literature
        )
    
    def _case_insensitive_search(self, gene_symbol, min_confidence, mirdb_threshold, 
                               targetscan_threshold, include_literature):
        """Case insensitive search"""
        return self._execute_regulation_query(
            "UPPER(g.gene_symbol) = UPPER(?)", [gene_symbol], min_confidence, 
            mirdb_threshold, targetscan_threshold, include_literature
        )
    
    def _variation_search(self, gene_symbol, min_confidence, mirdb_threshold, 
                         targetscan_threshold, include_literature):
        """Search using mapped variations"""
        
        # Get variations from mapper
        variations = self.mapper.get_all_variations(gene_symbol)
        standard_symbol = self.mapper.find_gene_symbol(gene_symbol)
        
        # Add standard symbol if found
        if standard_symbol:
            variations.append(standard_symbol)
        
        # Add common variations
        variations.extend([
            gene_symbol.upper(),
            gene_symbol.lower(), 
            gene_symbol.replace('-', ''),
            gene_symbol.replace('_', '')
        ])
        
        # Remove duplicates and limit
        variations = list(dict.fromkeys(variations))[:20]
        
        if not variations:
            return pd.DataFrame()
        
        placeholders = ','.join(['?' for _ in variations])
        condition = "g.gene_symbol IN ({})".format(placeholders)
        
        return self._execute_regulation_query(
            condition, variations, min_confidence, 
            mirdb_threshold, targetscan_threshold, include_literature
        )
    
    def _pattern_search(self, gene_symbol, min_confidence, mirdb_threshold, 
                       targetscan_threshold, include_literature):
        """Pattern-based fuzzy search"""
        
        search_pattern = "%{}%".format(gene_symbol)
        condition = "(g.gene_symbol LIKE ? OR g.description LIKE ?)"
        params = [search_pattern, search_pattern]
        
        results = self._execute_regulation_query(
            condition, params, min_confidence, 
            mirdb_threshold, targetscan_threshold, include_literature, limit=50
        )
        
        return results
    
    def _execute_regulation_query(self, condition, params, min_confidence, 
                                mirdb_threshold, targetscan_threshold, 
                                include_literature, limit=None):
        """Execute regulation query with enhanced error handling"""
        
        conn = self.get_connection()
        
        try:
            query = """
                SELECT 
                    n.gene_symbol as ncrna_symbol,
                    n.ncrna_type,
                    n.ncrna_subtype,
                    g.gene_symbol as target_gene_symbol,
                    g.gene_type,
                    g.description as gene_description,
                    r.regulation_type,
                    r.regulation_mechanism,
                    r.evidence_level,
                    r.confidence_score,
                    r.pubmed_ids,
                    g.chromosome as gene_chromosome,
                    g.start_position as gene_start,
                    g.end_position as gene_end
                FROM genes g 
                JOIN regulations r ON g.gene_id = r.target_gene_id
                JOIN ncrna_genes n ON r.ncrna_id = n.ncrna_id  
                WHERE {}
                AND r.confidence_score >= ?
            """.format(condition)
            
            query_params = list(params) + [min_confidence]
            
            # Add additional filters
            conditions = []
            if mirdb_threshold is not None:
                conditions.append("(r.pubmed_ids != 'miRDB_v6.0' OR r.confidence_score >= ?)")
                query_params.append(mirdb_threshold / 100.0)
            
            if targetscan_threshold is not None:
                conditions.append("(r.pubmed_ids != 'TargetScan_v8.0' OR r.confidence_score >= ?)")
                query_params.append(abs(targetscan_threshold))
            
            if not include_literature:
                conditions.append("r.evidence_level != 'experimental'")
            
            if conditions:
                query += " AND " + " AND ".join(conditions)
            
            query += " ORDER BY r.confidence_score DESC"
            
            if limit:
                query += " LIMIT {}".format(limit)
            
            logger.debug("Executing query: {} with params: {}".format(query, query_params))
            
            results = pd.read_sql_query(query, conn, params=query_params)
            return results
            
        except Exception as e:
            logger.error("Error executing query: {}".format(e))
            return pd.DataFrame()
        finally:
            conn.close()
    
    def search_ncrna_targets(self, ncrna_symbol,
                           min_confidence=0.0,
                           mirdb_threshold=None,
                           targetscan_threshold=None,
                           include_disease=True,
                           debug=False):
        
        logger.info("Searching targets for ncRNA: {}".format(ncrna_symbol))
        
        if debug:
            self._debug_ncrna_search(ncrna_symbol)
        
        # Try multiple search approaches
        search_methods = [
            ("Direct exact match", lambda s: s),
            ("Case insensitive", lambda s: s),
            ("Mapper variations", self._get_ncrna_variations)
        ]
        
        for method_name, transform_func in search_methods:
            logger.info("Trying ncRNA search method: {}".format(method_name))
            
            search_terms = transform_func(ncrna_symbol)
            if not isinstance(search_terms, list):
                search_terms = [search_terms]
            
            for term in search_terms[:10]:  # Limit variations
                results = self._execute_ncrna_query(term, min_confidence, mirdb_threshold, 
                                                  targetscan_threshold, method_name)
                
                if not results.empty:
                    results['searched_symbol'] = ncrna_symbol
                    results['matched_variation'] = term
                    results['search_method'] = method_name
                    
                    if include_disease:
                        results = self._add_disease_info(results, 'target_gene_symbol')
                    
                    logger.info("Found {} targets using {}".format(len(results), method_name))
                    return results
            
            logger.info("No results from {}".format(method_name))
        
        return pd.DataFrame()
    
    def _debug_ncrna_search(self, ncrna_symbol):
        """Debug ncRNA search"""
        logger.info("=== DEBUGGING NCRNA SEARCH FOR: {} ===".format(ncrna_symbol))
        
        # Check if ncRNA exists in mappings
        standard_symbol = self.mapper.find_ncrna_symbol(ncrna_symbol)
        logger.info("Standard ncRNA symbol from mapper: {}".format(standard_symbol))
        
        # Check database for ncRNA existence
        conn = self.get_connection()
        try:
            ncrna_count = conn.execute(
                "SELECT COUNT(*) FROM ncrna_genes WHERE gene_symbol = ?", 
                (ncrna_symbol,)
            ).fetchone()[0]
            logger.info("Direct match in ncrna_genes table: {}".format(ncrna_count))
            
            # Sample ncRNA symbols
            sample_ncrnas = conn.execute(
                "SELECT DISTINCT gene_symbol FROM ncrna_genes LIMIT 10"
            ).fetchall()
            logger.info("Sample ncRNA symbols: {}".format([n[0] for n in sample_ncrnas]))
            
        finally:
            conn.close()
        
        logger.info("=== END DEBUG ===")
    
    def _get_ncrna_variations(self, ncrna_symbol):
        """Get ncRNA variations from mapper"""
        variations = self.mapper.get_all_variations(ncrna_symbol)
        standard_symbol = self.mapper.find_ncrna_symbol(ncrna_symbol)
        
        if standard_symbol:
            variations.append(standard_symbol)
        
        # Add common miRNA variations
        if 'mir' in ncrna_symbol.lower():
            variations.extend([
                ncrna_symbol.replace('mir', 'miR'),
                ncrna_symbol.replace('miR', 'mir'),
                'hsa-' + ncrna_symbol,
                ncrna_symbol.replace('hsa-', '')
            ])
        
        return list(dict.fromkeys(variations))
    
    def _execute_ncrna_query(self, ncrna_symbol, min_confidence, mirdb_threshold, 
                           targetscan_threshold, search_method):
        """Execute ncRNA target query"""
        
        conn = self.get_connection()
        
        try:
            if search_method == "Case insensitive":
                condition = "UPPER(n.gene_symbol) = UPPER(?)"
            else:
                condition = "n.gene_symbol = ?"
            
            query = """
                SELECT 
                    n.gene_symbol as ncrna_symbol,
                    n.ncrna_type,
                    n.ncrna_subtype,
                    g.gene_symbol as target_gene_symbol,
                    g.gene_type,
                    g.description as gene_description,
                    r.regulation_type,
                    r.regulation_mechanism,
                    r.evidence_level,
                    r.confidence_score,
                    r.pubmed_ids,
                    g.chromosome as gene_chromosome,
                    g.start_position as gene_start,
                    g.end_position as gene_end
                FROM ncrna_genes n
                JOIN regulations r ON n.ncrna_id = r.ncrna_id
                JOIN genes g ON r.target_gene_id = g.gene_id
                WHERE {}
                AND r.confidence_score >= ?
            """.format(condition)
            
            params = [ncrna_symbol, min_confidence]
            
            # Add filters
            conditions = []
            if mirdb_threshold is not None:
                conditions.append("(r.pubmed_ids != 'miRDB_v6.0' OR r.confidence_score >= ?)")
                params.append(mirdb_threshold / 100.0)
            
            if targetscan_threshold is not None:
                conditions.append("(r.pubmed_ids != 'TargetScan_v8.0' OR r.confidence_score >= ?)")
                params.append(abs(targetscan_threshold))
            
            if conditions:
                query += " AND " + " AND ".join(conditions)
            
            query += " ORDER BY r.confidence_score DESC LIMIT 1000"
            
            results = pd.read_sql_query(query, conn, params=params)
            return results
            
        except Exception as e:
            logger.error("Error executing ncRNA query: {}".format(e))
            return pd.DataFrame()
        finally:
            conn.close()
    
    def search_disease_associations(self, disease_name):
        conn = self.get_connection()
        
        try:
            query = """
                SELECT DISTINCT
                    g.gene_symbol,
                    g.gene_type,
                    g.description as gene_description,
                    g.chromosome,
                    g.start_position,
                    g.end_position,
                    h.hpo_id,
                    h.hpo_name as hpo_disease,
                    o.omim_id,
                    o.title as omim_disease,
                    o.inheritance_pattern,
                    o.clinical_features,
                    gda.association_type,
                    gda.evidence_level,
                    gda.pubmed_ids
                FROM genes g
                JOIN gene_disease_associations gda ON g.gene_id = gda.gene_id
                LEFT JOIN hpo_terms h ON gda.hpo_id = h.hpo_id
                LEFT JOIN omim_entries o ON gda.omim_id = o.omim_id
                WHERE (h.hpo_name LIKE ? OR o.title LIKE ? OR o.phenotype LIKE ?)
                ORDER BY g.gene_symbol
            """
            
            search_pattern = "%{}%".format(disease_name)
            params = [search_pattern, search_pattern, search_pattern]
            
            results = pd.read_sql_query(query, conn, params=params)
            return results
            
        finally:
            conn.close()
    
    def _add_disease_info(self, results, gene_column):
        """Add disease information to results"""
        if results.empty:
            return results
        
        conn = self.get_connection()
        
        try:
            gene_symbols = results[gene_column].unique()
            placeholders = ','.join(['?' for _ in gene_symbols])
            
            disease_query = """
                SELECT 
                    g.gene_symbol,
                    GROUP_CONCAT(DISTINCT h.hpo_name) as hpo_diseases,
                    GROUP_CONCAT(DISTINCT o.title) as omim_diseases
                FROM genes g
                LEFT JOIN gene_disease_associations gda ON g.gene_id = gda.gene_id
                LEFT JOIN hpo_terms h ON gda.hpo_id = h.hpo_id
                LEFT JOIN omim_entries o ON gda.omim_id = o.omim_id
                WHERE g.gene_symbol IN ({})
                GROUP BY g.gene_symbol
            """.format(placeholders)
            
            disease_info = pd.read_sql_query(disease_query, conn, params=list(gene_symbols))
            
            if not disease_info.empty:
                results = results.merge(disease_info, left_on=gene_column, right_on='gene_symbol', how='left')
        
        finally:
            conn.close()
        
        return results
    
    def export_results(self, df, filename, format="csv"):
        """Export results to file"""
        output_dir = Path("results")
        output_dir.mkdir(exist_ok=True)
        
        if format.lower() == "csv":
            output_file = output_dir / "{}.csv".format(filename)
            df.to_csv(output_file, index=False)
        elif format.lower() == "json":
            output_file = output_dir / "{}.json".format(filename)
            df.to_json(output_file, orient='records', indent=2)
        
        print("Results exported to: {}".format(output_file))
        return output_file
    
    def get_stats(self):
        """Get database statistics"""
        conn = self.get_connection()
        stats = {}
        
        try:
            stats['total_genes'] = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
            stats['total_ncrnas'] = conn.execute("SELECT COUNT(*) FROM ncrna_genes").fetchone()[0]
            stats['total_regulations'] = conn.execute("SELECT COUNT(*) FROM regulations").fetchone()[0]
            stats['total_hpo_terms'] = conn.execute("SELECT COUNT(*) FROM hpo_terms").fetchone()[0]
            stats['total_omim_entries'] = conn.execute("SELECT COUNT(*) FROM omim_entries").fetchone()[0]
            stats['total_disease_associations'] = conn.execute("SELECT COUNT(*) FROM gene_disease_associations").fetchone()[0]
            
            # Get mapping statistics
            mapping_stats = self.mapper.get_mapping_stats()
            stats.update(mapping_stats)
            
            # Get breakdowns
            ncrna_types = conn.execute("""
                SELECT ncrna_type, COUNT(*) as count 
                FROM ncrna_genes 
                GROUP BY ncrna_type
            """).fetchall()
            stats['ncrna_types'] = {row[0]: row[1] for row in ncrna_types}
            
            regulation_sources = conn.execute("""
                SELECT pubmed_ids, COUNT(*) as count
                FROM regulations
                GROUP BY pubmed_ids
            """).fetchall()
            stats['regulation_sources'] = {row[0]: row[1] for row in regulation_sources}
            
            evidence_levels = conn.execute("""
                SELECT evidence_level, COUNT(*) as count
                FROM regulations
                GROUP BY evidence_level
            """).fetchall()
            stats['evidence_levels'] = {row[0]: row[1] for row in evidence_levels}
            
            quality_ranges = conn.execute("""
                SELECT 
                    CASE 
                        WHEN confidence_score >= 0.9 THEN 'High (≥0.9)'
                        WHEN confidence_score >= 0.8 THEN 'Medium (0.8-0.9)'
                        WHEN confidence_score >= 0.7 THEN 'Low (0.7-0.8)'
                        ELSE 'Very Low (<0.7)'
                    END as quality_range,
                    COUNT(*) as count
                FROM regulations
                WHERE confidence_score IS NOT NULL
                GROUP BY quality_range
            """).fetchall()
            stats['quality_ranges'] = {row[0]: row[1] for row in quality_ranges}
            
        finally:
            conn.close()
        
        return stats

def print_stats(stats):
    print("=" * 60)
    print("DATABASE STATISTICS")
    print("=" * 60)
    print("  Total Genes: {:,}".format(stats['total_genes']))
    print("  Total ncRNAs: {:,}".format(stats['total_ncrnas']))
    print("  Total Regulations: {:,}".format(stats['total_regulations']))
    print("  HPO Terms: {:,}".format(stats['total_hpo_terms']))
    print("  OMIM Entries: {:,}".format(stats['total_omim_entries']))
    print("  Disease Associations: {:,}".format(stats['total_disease_associations']))
    print()
    print("MAPPING STATISTICS")
    print("  Gene Mappings: {:,}".format(stats['total_gene_mappings']))
    print("  ncRNA Mappings: {:,}".format(stats['total_ncrna_mappings']))
    print("  RefSeq Mappings: {:,}".format(stats.get('refseq_mappings', 0)))
    print("  Unique Symbols: {:,}".format(stats['unique_symbols']))

def print_gene_results(results, gene_symbol):
    if results.empty:
        print("No ncRNA regulators found for gene: {}".format(gene_symbol))
        return
    
    print("ncRNA Regulators for Gene: {}".format(gene_symbol))
    print("=" * 70)
    
    if 'search_method' in results.columns:
        methods = results['search_method'].unique()
        print("Search methods used: {}".format(", ".join(methods)))
        print()
    
    for ncrna_type in results['ncrna_type'].unique():
        type_results = results[results['ncrna_type'] == ncrna_type]
        print("{} Regulators ({}):".format(ncrna_type.upper(), len(type_results)))
        print("-" * 40)
        
        for _, row in type_results.head(10).iterrows():  # Show top 10
            print("  {} → {}".format(row['ncrna_symbol'], row['target_gene_symbol']))
            print("    Regulation: {} | Confidence: {:.3f}".format(
                row['regulation_type'], row['confidence_score']))
            print("    Evidence: {} | Source: {}".format(
                row['evidence_level'], row['pubmed_ids']))
            
            if 'hpo_diseases' in row and pd.notna(row['hpo_diseases']):
                print("    HPO Diseases: {}".format(row['hpo_diseases'][:100]))
            if 'omim_diseases' in row and pd.notna(row['omim_diseases']):
                print("    OMIM Diseases: {}".format(row['omim_diseases'][:100]))
            print()

def print_ncrna_results(results, ncrna_symbol):
    if results.empty:
        print("No target genes found for ncRNA: {}".format(ncrna_symbol))
        return
    
    print("Target Genes for ncRNA: {}".format(ncrna_symbol))
    print("=" * 70)
    
    if 'search_method' in results.columns:
        methods = results['search_method'].unique()
        print("Search methods used: {}".format(", ".join(methods)))
        print()
    
    for reg_type in results['regulation_type'].unique():
        type_results = results[results['regulation_type'] == reg_type]
        print("{} Regulation ({}):".format(reg_type.upper(), len(type_results)))
        print("-" * 40)
        
        for _, row in type_results.head(10).iterrows():  # Show top 10
            print("  {} → {}".format(row['ncrna_symbol'], row['target_gene_symbol']))
            print("    Confidence: {:.3f} | Mechanism: {}".format(
                row['confidence_score'], row['regulation_mechanism'] or 'Unknown'))
            print("    Evidence: {} | Gene Type: {}".format(
                row['evidence_level'], row['gene_type']))
            
            if 'hpo_diseases' in row and pd.notna(row['hpo_diseases']):
                print("    HPO Diseases: {}".format(row['hpo_diseases'][:100]))
            if 'omim_diseases' in row and pd.notna(row['omim_diseases']):
                print("    OMIM Diseases: {}".format(row['omim_diseases'][:100]))
            print()

def print_disease_results(results, disease_name):
    if results.empty:
        print("No genes found associated with disease: {}".format(disease_name))
        return
    
    print("Genes Associated with Disease: {}".format(disease_name))
    print("=" * 70)
    
    for _, row in results.iterrows():
        print("Gene: {} ({})".format(row['gene_symbol'], row['gene_type']))
        print("Location: chr{}:{}-{}".format(
            row['chromosome'], row['start_position'], row['end_position']))
        print("  Association: {} | Evidence: {}".format(
            row['association_type'], row['evidence_level']))
        
        if pd.notna(row['hpo_disease']):
            print("  HPO: {} ({})".format(row['hpo_disease'], row['hpo_id']))
        if pd.notna(row['omim_disease']):
            print("  OMIM: {} ({})".format(row['omim_disease'], row['omim_id']))
        if pd.notna(row['inheritance_pattern']):
            print("  Inheritance: {}".format(row['inheritance_pattern']))
        if pd.notna(row['clinical_features']):
            print("  Clinical Features: {}".format(row['clinical_features']))
        print()

def main():
    parser = argparse.ArgumentParser(description='Enhanced ncRNA database search')
    parser.add_argument('--search-type', required=True, 
                       choices=['stats', 'gene-to-ncrna', 'ncrna-to-gene', 'disease'],
                       help='Type of search to perform')
    parser.add_argument('--gene', help='Gene symbol to search for')
    parser.add_argument('--ncrna', help='ncRNA symbol to search for')
    parser.add_argument('--disease', help='Disease name to search for')
    parser.add_argument('--min-confidence', type=float, default=0.0,
                       help='Minimum confidence score (0.0-1.0)')
    parser.add_argument('--mirdb-threshold', type=float,
                       help='miRDB score threshold (0-100)')
    parser.add_argument('--targetscan-threshold', type=float,
                       help='TargetScan context score threshold (e.g., -0.3)')
    parser.add_argument('--no-literature', action='store_true',
                       help='Exclude literature-based regulations')
    parser.add_argument('--no-disease', action='store_true',
                       help='Exclude disease association information')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode for troubleshooting')
    parser.add_argument('--export', choices=['csv', 'json'],
                       help='Export results to file')
    parser.add_argument('--output-name', help='Custom output filename (without extension)')
    parser.add_argument('--db-path', default='database/ncrna_regulation_db.sqlite',
                       help='Path to database file')
    parser.add_argument('--config-path', default='config.yaml',
                       help='Path to configuration file')
    
    args = parser.parse_args()
    
    try:
        db = EnhancedNCRNADatabase(args.db_path, args.config_path)
        
        if args.search_type == 'stats':
            stats = db.get_stats()
            print_stats(stats)
            
            if args.export:
                stats_df = pd.DataFrame([stats])
                filename = args.output_name or "database_stats"
                db.export_results(stats_df, filename, args.export)
            
        elif args.search_type == 'gene-to-ncrna':
            if not args.gene:
                print("Gene parameter required for gene-to-ncrna search")
                sys.exit(1)
            
            results = db.search_gene_regulators(
                gene_symbol=args.gene,
                min_confidence=args.min_confidence,
                mirdb_threshold=args.mirdb_threshold,
                targetscan_threshold=args.targetscan_threshold,
                include_literature=not args.no_literature,
                include_disease=not args.no_disease,
                debug=args.debug
            )
            
            print_gene_results(results, args.gene)
            
            if args.export and not results.empty:
                filename = args.output_name or "{}_regulators".format(args.gene)
                db.export_results(results, filename, args.export)
                
        elif args.search_type == 'ncrna-to-gene':
            if not args.ncrna:
                print("ncrna parameter required for ncrna-to-gene search")
                sys.exit(1)
            
            results = db.search_ncrna_targets(
                ncrna_symbol=args.ncrna,
                min_confidence=args.min_confidence,
                mirdb_threshold=args.mirdb_threshold,
                targetscan_threshold=args.targetscan_threshold,
                include_disease=not args.no_disease,
                debug=args.debug
            )
            
            print_ncrna_results(results, args.ncrna)
            
            if args.export and not results.empty:
                filename = args.output_name or "{}_targets".format(args.ncrna)
                db.export_results(results, filename, args.export)
                
        elif args.search_type == 'disease':
            if not args.disease:
                print("Disease parameter required for disease search")
                sys.exit(1)
            
            results = db.search_disease_associations(args.disease)
            
            print_disease_results(results, args.disease)
            
            if args.export and not results.empty:
                filename = args.output_name or "{}_genes".format(args.disease.replace(' ', '_'))
                db.export_results(results, filename, args.export)
    
    except Exception as e:
        print("Error: {}".format(e))
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
