import pandas as pd
import gffutils
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import tempfile
import os
import warnings
import argparse
warnings.filterwarnings('ignore')

# Create output directory
def create_output_dir(report_dir):
    """Create the output directory if it doesn't exist"""
    path = Path(report_dir)
    path.mkdir(exist_ok=True, parents=True)
    return path

# Process a GFF file to extract CDS features
def process_gff_file(gff_file):
    """Extract CDS features from a GFF file"""
    # Convert Path to string to avoid iteration issues
    gff_str = str(gff_file)
    db_path = tempfile.mktemp()
    
    try:
        # Parse GFF file
        db = gffutils.create_db(
            gff_str, 
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

# Process bedtools intersect output
def process_overlap_file(overlap_file):
    """Process bedtools overlap file"""
    # Define columns
    cols = ['chrom', 'start', 'end', 'read_name', 'score', 'strand', 
            'thick_start', 'thick_end', 'rgb', 'blocks', 'block_sizes', 'block_starts',
            'gff_chrom', 'gff_source', 'gff_type', 'gff_start', 'gff_end', 
            'gff_score', 'gff_strand', 'gff_phase', 'gff_attributes', 'overlap_bp']
    
    # Read file
    try:
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
    except Exception as e:
        print(f"Error processing {overlap_file}: {e}")
        return pd.DataFrame()

# Plot CDS distribution
def plot_cds_features(cds_df, report_dir, title="CDS Features", show=True):
    """Plot CDS feature distribution"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
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
    
    # Save the figure
    filename = f"{title.replace(' ', '_').replace(':', '')}.png"
    plot_path = report_dir / filename
    fig.savefig(plot_path, dpi=300)
    print(f"Saved CDS plot to {plot_path}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
        
    return plot_path

# Plot overlap summary
def plot_overlap_summary(overlap_df, report_dir, title="Read-CDS Overlap Summary", show=True):
    """Plot summary of read-CDS overlaps"""
    if overlap_df.empty:
        print("No overlap data to plot")
        return None
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Overlap length distribution
    sns.histplot(overlap_df['overlap_bp'], bins=30, kde=True, ax=ax1)
    ax1.set_title('Overlap Length Distribution')
    ax1.set_xlabel('Overlap Length (bp)')
    ax1.set_ylabel('Count')
    
    # Stats for overlap
    stats = f"Total overlaps: {len(overlap_df)}\n"
    stats += f"Unique reads: {overlap_df['read_name'].nunique()}\n"
    stats += f"Unique products: {overlap_df['product'].nunique()}\n"
    stats += f"Mean overlap: {overlap_df['overlap_bp'].mean():.1f} bp"
    ax1.text(0.95, 0.95, stats, transform=ax1.transAxes, 
             va='top', ha='right', bbox=dict(boxstyle='round', fc='white', alpha=0.7))
    
    # Top products
    product_counts = overlap_df.groupby('product')['read_name'].nunique().sort_values(ascending=False).head(10)
    sns.barplot(x=product_counts.values, y=product_counts.index, ax=ax2)
    ax2.set_title('Top 10 Products by Read Count')
    ax2.set_xlabel('Number of Unique Reads')
    
    plt.tight_layout()
    plt.suptitle(title, fontsize=16)
    plt.subplots_adjust(top=0.85)
    
    # Save the figure
    filename = f"{title.replace(' ', '_').replace(':', '')}.png"
    plot_path = report_dir / filename
    fig.savefig(plot_path, dpi=300)
    print(f"Saved overlap plot to {plot_path}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
        
    return plot_path
    
# Compare samples
def compare_samples(all_overlaps, report_dir, title="Sample Comparison", show=True):
    """Compare CDS overlaps across samples"""
    if all_overlaps['sample'].nunique() <= 1:
        print("Need multiple samples for comparison")
        return None
        
    plt.figure(figsize=(10, 6))
    sample_counts = all_overlaps.groupby('sample')['read_name'].nunique()
    sns.barplot(x=sample_counts.index, y=sample_counts.values)
    plt.title('CDS Overlaps by Sample')
    plt.xlabel('Sample')
    plt.ylabel('Number of Reads with CDS Overlaps')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    comparison_path = report_dir / 'sample_comparison.png'
    plt.savefig(comparison_path, dpi=300)
    
    if show:
        plt.show()
    else:
        plt.close()
        
    return comparison_path

# Main analysis function
def analyze_samples(gff_file, sample_files, report_dir='output/report', organism_name="Unknown", show_plots=True):
    """Analyze GFF file and sample overlap files"""
    # Create output directory
    report_dir = create_output_dir(report_dir)
    
    print(f"Processing GFF file: {gff_file}")
    cds_df = process_gff_file(gff_file)
    
    # Save CDS dataframe
    cds_csv_path = report_dir / 'genomic_cds.csv'
    cds_df.to_csv(cds_csv_path, index=False)
    print(f"Found {len(cds_df)} CDS features")
    print(f"Saved CDS data to {cds_csv_path}")
    
    # Plot CDS features
    plot_cds_features(cds_df, report_dir, title=f"{organism_name} CDS Analysis", show=show_plots)
    
    # Process sample files
    print(f"\nProcessing {len(sample_files)} sample overlap files...")
    
    # Create a summary dataframe for all overlaps
    all_overlaps = pd.DataFrame()
    
    for overlap_file in sample_files:
        overlap_path = Path(overlap_file)
        print(f"  Processing {overlap_path.name}...")
        
        overlap_df = process_overlap_file(overlap_path)
        
        if not overlap_df.empty:
            # Add sample ID
            overlap_df['sample'] = overlap_path.stem
            
            # Save individual sample results
            sample_csv_path = report_dir / f"{overlap_path.stem}_overlaps.csv"
            overlap_df.to_csv(sample_csv_path, index=False)
            print(f"    Saved to {sample_csv_path}")
            
            # Append to all overlaps
            all_overlaps = pd.concat([all_overlaps, overlap_df])
            
            # Basic stats
            print(f"    Found {len(overlap_df)} CDS overlaps")
            print(f"    Unique reads: {overlap_df['read_name'].nunique()}")
            print(f"    Unique products: {overlap_df['product'].nunique()}")
            
            # Plot individual sample summary
            plot_overlap_summary(
                overlap_df, 
                report_dir,
                title=f"{overlap_path.stem}: Read-CDS Overlap Summary",
                show=False
            )
    
    # Process combined results
    if not all_overlaps.empty:
        # Save all overlaps
        all_overlaps_csv = report_dir / 'all_overlaps.csv'
        all_overlaps.to_csv(all_overlaps_csv, index=False)
        print(f"Saved combined overlaps to {all_overlaps_csv}")
        
        # Plot overall summary
        plot_overlap_summary(
            all_overlaps, 
            report_dir,
            title="All Samples: Read-CDS Overlap Summary",
            show=show_plots
        )
        
        # Sample comparison (if multiple samples)
        if all_overlaps['sample'].nunique() > 1:
            compare_samples(all_overlaps, report_dir, show=show_plots)
            
    print(f"\nAnalysis complete! All outputs saved to {report_dir}")
    return all_overlaps
