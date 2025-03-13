import pandas as pd
import gffutils
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import tempfile
import os

def process_gff_file(gff_file):
    """Extract CDS features from a GFF file"""
    # Create a temporary database
    db_path = tempfile.mktemp()
    
    try:
        # Parse GFF file
        db = gffutils.create_db(
            gff_file, 
            dbfn=db_path,
            force=True,
            keep_order=True,
            merge_strategy='merge'
        )
        
        # Extract CDS features
        cds_records = []
        for feature in db.features_of_type('CDS'):
            attrs = feature.attributes
            
            # Get product name
            product = attrs.get('product', [''])[0] if 'product' in attrs else ''
            if not product and 'locus_tag' in attrs:
                product = attrs.get('locus_tag', [''])[0]
                
            # Create record
            record = {
                'seqid': feature.seqid,
                'start': feature.start,
                'end': feature.end,
                'strand': feature.strand,
                'length': feature.end - feature.start + 1,
                'product': product
            }
            cds_records.append(record)
        
        # Convert to DataFrame
        cds_df = pd.DataFrame(cds_records)
        
    finally:
        # Clean up
        if os.path.exists(db_path):
            os.remove(db_path)
    
    return cds_df

def process_overlap_file(overlap_file):
    """Process bedtools overlap file"""
    # Define columns
    cols = ['chrom', 'start', 'end', 'read_name', 'score', 'strand', 
            'thick_start', 'thick_end', 'rgb', 'blocks', 'block_sizes', 'block_starts',
            'gff_chrom', 'gff_source', 'gff_type', 'gff_start', 'gff_end', 
            'gff_score', 'gff_strand', 'gff_phase', 'gff_attributes', 'overlap_bp']
    
    # Read file
    df = pd.read_csv(overlap_file, sep='\t', names=cols)
    
    # Filter for CDS
    cds_df = df[df['gff_type'] == 'CDS'].copy()
    
    # Extract product from attributes
    def get_product(attr_str):
        for attr in attr_str.split(';'):
            if attr.startswith('product='):
                return attr.replace('product=', '')
        return ''
    
    cds_df['product'] = cds_df['gff_attributes'].apply(get_product)
    
    # Simplify dataframe
    result = cds_df[['read_name', 'chrom', 'start', 'end', 'strand',
                     'gff_start', 'gff_end', 'gff_strand', 'product', 'overlap_bp']]
    
    return result

def plot_cds_features(cds_df, title="CDS Features"):
    """Plot CDS feature distribution"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Length distribution
    sns.histplot(cds_df['length'], bins=30, kde=True, ax=ax1)
    ax1.set_title('CDS Length Distribution')
    ax1.set_xlabel('Length (bp)')
    ax1.set_ylabel('Count')
    
    # Stats text
    stats = f"Total CDSs: {len(cds_df)}\n"
    stats += f"Mean length: {cds_df['length'].mean():.1f} bp\n"
    stats += f"Median length: {cds_df['length'].median():.1f} bp"
    ax1.text(0.95, 0.95, stats, transform=ax1.transAxes, 
             va='top', ha='right', bbox=dict(boxstyle='round', fc='white', alpha=0.7))
    
    # CDS per chromosome
    chr_counts = cds_df['seqid'].value_counts().head(10)
    sns.barplot(x=chr_counts.index, y=chr_counts.values, ax=ax2)
    ax2.set_title('CDS Count by Chromosome')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    plt.suptitle(title, fontsize=16)
    plt.subplots_adjust(top=0.85)
    
    return fig

# Main function
def analyze_gff_files(gff_dir):
    """Process all GFF files in a directory"""
    gff_dir = Path(gff_dir)
    results = {}
    
    # Find all GFF files
    gff_files = list(gff_dir.glob('*.gff'))
    print(f"Found {len(gff_files)} GFF files")
    
    for gff_file in gff_files:
        print(f"Processing {gff_file.name}...")
        cds_df = process_gff_file(gff_file)
        
        if not cds_df.empty:
            # Store results
            results[gff_file.stem] = cds_df
            
            # Plot and save
            fig = plot_cds_features(cds_df, title=f"CDS Analysis: {gff_file.stem}")
            fig.savefig(f"{gff_file.stem}_cds_analysis.png", dpi=300)
            plt.close(fig)
            
            # Save CSV
            cds_df.to_csv(f"{gff_file.stem}_cds.csv", index=False)
        else:
            print(f"No CDS features found in {gff_file.name}")
    
    return results

# Example usage in notebook:
# analyze_gff_files('refs/GCF_000196875.2')