# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from math import log, exp
warnings.filterwarnings('ignore')

# Create output directory
def create_output_dir(report_dir):
    """Create the output directory if it doesn't exist"""
    path = Path(report_dir)
    path.mkdir(exist_ok=True, parents=True)
    return path



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

# Plot overlap summary with genomic position analysis
def plot_overlap_summary(overlap_df, report_dir, title="Read-CDS Overlap Summary", show=True):
    """Plot summary of read-CDS overlaps and genomic position distribution"""
    if overlap_df.empty:
        print("No overlap data to plot")
        return None
    
    # Create a figure with 2 rows, 2 columns
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(2, 2)
    
    # Overlap length distribution (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    sns.histplot(overlap_df['overlap_bp'], bins=30, kde=True, ax=ax1)
    ax1.set_title('CDS Overlap Length Distribution')
    ax1.set_xlabel('Overlap Length (bp)')
    ax1.set_ylabel('Count')
    
    # Stats for overlap
    stats = f"Total overlaps: {len(overlap_df)}\n"
    stats += f"Unique reads: {overlap_df['read_name'].nunique()}\n"
    stats += f"Unique products: {overlap_df['product'].nunique()}\n"
    stats += f"Mean overlap: {overlap_df['overlap_bp'].mean():.1f} bp"
    ax1.text(0.95, 0.95, stats, transform=ax1.transAxes, 
             va='top', ha='right', bbox=dict(boxstyle='round', fc='white', alpha=0.7))
    
    # Top products (top right)
    ax2 = fig.add_subplot(gs[0, 1])
    product_counts = overlap_df.groupby('product')['read_name'].nunique().sort_values(ascending=False).head(10)
    sns.barplot(x=product_counts.values, y=product_counts.index, ax=ax2)
    ax2.set_title('Top 10 Products by Read Count')
    ax2.set_xlabel('Number of Unique Reads')
    
    # Analyze exact genomic positions (bottom left)
    ax3 = fig.add_subplot(gs[1, 0])
    
    # Use exact start positions for each read to identify repeated insertion sites
    position_counts = overlap_df.groupby(['chrom', 'start']).size().reset_index(name='count')
    
    # Analyze how many positions are hit multiple times
    position_distribution = position_counts['count'].value_counts().sort_index()
    
    # Convert to more readable format: X positions hit exactly Y times
    position_summary = pd.DataFrame({
        'Times_Hit': position_distribution.index,
        'Number_of_Positions': position_distribution.values
    })
    
    # Calculate total number of unique positions
    total_positions = len(position_counts)
    
    # Calculate total number of insertions
    total_insertions = position_counts['count'].sum()
    
    # Plot the distribution
    sns.barplot(x='Times_Hit', y='Number_of_Positions', data=position_summary, ax=ax3)
    ax3.set_title('Genomic Position Hit Frequency')
    ax3.set_xlabel('Number of Times Position Was Hit')
    ax3.set_ylabel('Number of Positions')
    
    # Add stats to the plot
    jackpot_stats = f"Unique genomic positions: {total_positions}\n"
    jackpot_stats += f"Total insertions: {total_insertions}\n"
    
    # Calculate some diversity metrics
    if len(position_counts) > 0:
        # Position with most hits
        max_position = position_counts.loc[position_counts['count'].idxmax()]
        jackpot_stats += f"Most hit position: {max_position['chrom']}:{max_position['start']} ({max_position['count']} hits)\n"
        
        # Calculate proportion of positions hit only once
        single_hits = position_distribution.get(1, 0)
        jackpot_stats += f"Positions hit once: {single_hits} ({single_hits/total_positions*100:.1f}%)\n"
        
        # Calculate proportion of insertions in top 10% of positions
        if len(position_counts) >= 10:
            top_positions = position_counts.nlargest(int(len(position_counts)*0.1), 'count')
            top_insertions = top_positions['count'].sum()
            jackpot_stats += f"Top 10% positions account for {top_insertions/total_insertions*100:.1f}% of insertions"
    
    ax3.text(0.95, 0.95, jackpot_stats, transform=ax3.transAxes,
             va='top', ha='right', bbox=dict(boxstyle='round', fc='white', alpha=0.7))
    
    # Genomic position heatmap by chromosome (bottom right)
    ax4 = fig.add_subplot(gs[1, 1])
    
    # Get insertion density by chromosome
    chrom_insertions = position_counts.groupby('chrom').agg({
        'count': 'sum',  # total insertions on chromosome
        'chrom': 'size'  # unique positions on chromosome
    }).rename(columns={'chrom': 'unique_positions'})
    
    # Calculate insertion density
    chrom_insertions['insertion_density'] = chrom_insertions['count'] / chrom_insertions['unique_positions']
    
    # Sort by insertion density descending
    chrom_insertions = chrom_insertions.sort_values('insertion_density', ascending=False)
    
    # Take top 15 chromosomes to plot
    top_chroms = chrom_insertions.head(15).index
    
    # Filter to just those chromosomes
    filtered_pos = position_counts[position_counts['chrom'].isin(top_chroms)]
    
    # Create a crosstab to see distribution of hit counts across chromosomes
    if not filtered_pos.empty:
        # Bin the counts to avoid too many columns
        filtered_pos['count_bin'] = pd.cut(
            filtered_pos['count'], 
            bins=[0, 1, 2, 3, 5, 10, 20, float('inf')],
            labels=['1', '2', '3', '4-5', '6-10', '11-20', '>20']
        )
        
        # Create crosstab
        chrom_hit_dist = pd.crosstab(
            filtered_pos['chrom'], 
            filtered_pos['count_bin'],
            normalize='index'
        ).fillna(0)
        
        # Plot heatmap
        sns.heatmap(chrom_hit_dist, cmap='viridis', ax=ax4, annot=True, fmt='.2f', cbar=True)
        ax4.set_title('Position Hit Distribution by Chromosome')
        ax4.set_xlabel('Times Each Position Was Hit')
        ax4.set_ylabel('Chromosome')
    else:
        ax4.text(0.5, 0.5, "Not enough data for chromosome heatmap", 
                ha='center', va='center', transform=ax4.transAxes)
    
    plt.tight_layout()
    plt.suptitle(title, fontsize=16)
    plt.subplots_adjust(top=0.95)
    
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
def analyze_samples(sample_files, report_dir='output/report', organism_name="Unknown", show_plots=True):
    """Analyze GFF file and sample overlap files"""
    # Create output directory
    report_dir = create_output_dir(report_dir)
    
    # Process sample files directly - no genomic CDS analysis
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
