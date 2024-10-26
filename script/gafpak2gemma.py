import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import subprocess
import shutil

def z_score_normalization(matrix):
    """
    Apply Z-score normalization to each column of the matrix.
    μj is the mean of the column, and σj is the standard deviation of the column.
    """
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    normalized = matrix.copy()
    
    for col in normalized.columns:
        mean_col = np.nanmean(normalized[col])
        std_col = np.nanstd(normalized[col])
        
        if std_col > 0:
            normalized[col] = (normalized[col] - mean_col) / std_col
    
    return normalized

def min_max_scale(matrix, min_value=0, max_value=2):
    """
    Scale the values of the matrix to be within the range [min_value, max_value].
    """
    scaled = matrix.copy()
    for col in scaled.columns:
        col_min = scaled[col].min()
        col_max = scaled[col].max()
        
        if col_max - col_min > 0:
            scaled[col] = (scaled[col] - col_min) / (col_max - col_min) * (max_value - min_value) + min_value
    
    return scaled

def cleanup_dir(directory):
    """Clean up a directory if it exists."""
    if directory and os.path.exists(directory):
        shutil.rmtree(directory)

def process_matrix_with_plots(input_file, output_file, plot_file_panel, dropped_columns_file='dropped_columns.txt', pheno_file=None):
    # Check if input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    if pheno_file and not os.path.exists(pheno_file):
        raise FileNotFoundError(f"Phenotype file not found: {pheno_file}")
        
    # Read the matrix from the input CSV file
    matrix = pd.read_csv(input_file, index_col=0, skiprows=1)
    print(f'Matrix size input: {matrix.shape}')

    # Ensure all data is numeric, drop non-numeric data
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    
    # Identify and save dropped columns
    nan_columns = matrix.columns[matrix.isna().all()].tolist()
    with open(dropped_columns_file, 'w') as file:
        for col in nan_columns:
            file.write(f"{col}\n")
    
    # Drop NaN rows and columns
    matrix = matrix.dropna(axis=0, how='all')
    matrix = matrix.dropna(axis=1, how='all')

    # Create plots
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    
    # Original values plot
    axes[0].hist(matrix.values.flatten(), bins=30, alpha=0.7, color='orange', edgecolor='black')
    axes[0].set_xlabel('Value')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Original Input Values')
    axes[0].set_yscale('log')
    axes[0].grid(True)
    
    # Z-score normalization
    z_normalized_matrix = z_score_normalization(matrix)
    print(f'Matrix size after Z-score normalization: {z_normalized_matrix.shape}')

    # Z-score plot
    axes[1].hist(z_normalized_matrix.values.flatten(), bins=30, alpha=0.7, color='blue', edgecolor='black')
    axes[1].set_xlabel('Value')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('After Z-Score Normalization')
    axes[1].grid(True)
    
    # Min-Max scaling
    final_normalized_matrix = min_max_scale(z_normalized_matrix, 0, 2)
    print(f'Matrix size after rescaling: {final_normalized_matrix.shape}')

    # Reset the index to make it a column
    final_normalized_matrix = final_normalized_matrix.reset_index()
    
    # Add X and Y columns after the index column
    first_col = final_normalized_matrix.columns[0]  # Get the name of the first column
    final_normalized_matrix.insert(1, 'X', 'X')
    final_normalized_matrix.insert(2, 'Y', 'Y')

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save to CSV, if you change to 2 cifres GEMMA works
    final_normalized_matrix.to_csv(output_file, index=False, float_format='%.5f')
    print(f'Output saved to: {output_file}')

    # Final distribution plot
    numeric_data = final_normalized_matrix.select_dtypes(include=[np.number])
    axes[2].hist(numeric_data.values.flatten(), bins=30, alpha=0.7, color='green', edgecolor='black')
    axes[2].set_xlabel('Value')
    axes[2].set_ylabel('Frequency')
    axes[2].set_title('After Scaling to [0, 2]')
    axes[2].grid(True)
    
    # Create plot output directory if it doesn't exist
    plot_dir = os.path.dirname(plot_file_panel)
    if plot_dir and not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
        
    # Save plot
    plt.tight_layout()
    plt.savefig(plot_file_panel)
    print(f'Plot saved to: {plot_file_panel}')

    # Run GEMMA if phenotype file is provided
    if pheno_file:
        print("\nStarting GEMMA analysis...")
        run_simple_gemma(output_file, pheno_file)
        print("GEMMA analysis completed")

def run_simple_gemma(norm_file, pheno_file):
    """Run basic GEMMA analysis with kinship calculation and association."""
    try:
        # Create output directory
        gemma_outdir = "gemma_out"
        cleanup_dir(gemma_outdir)
        os.makedirs(gemma_outdir)
        
        # Step 1: Calculate kinship matrix
        print("Running kinship calculation...")
        kinship_cmd = (
            f"gemma -g {norm_file} "
            f"-p {pheno_file} "
            f"-gk 1 "
            f"-miss 1.0 "
            f"-notsnp "
            f"-r2 1.0 "
            f"-hwe 0 "
            f"-outdir {gemma_outdir} "
            f"-o kinship"
        )
        
        print(f"Running command: {kinship_cmd}")
        result = subprocess.run(kinship_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Warning: Kinship calculation failed")
            print(f"Error output: {result.stderr}")
            return
            
        # Step 2: Run association analysis
        print("Running association analysis...")
        kinship_file = f"{gemma_outdir}/kinship.cXX.txt"
        if os.path.exists(kinship_file):
            assoc_cmd = (
                f"gemma -g {norm_file} "
                f"-p {pheno_file} "
                f"-n 1 "
                f"-k {kinship_file} "
                f"-lmm "
                f"-miss 1.0 "
                f"-notsnp "
                f"-r2 1.0 "
                f"-hwe 0 "
                f"-outdir {gemma_outdir} "
                f"-o assoc"
            )
            
            print(f"Running command: {assoc_cmd}")
            result = subprocess.run(assoc_cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Warning: Association analysis failed")
                print(f"Error output: {result.stderr}")
        else:
            print(f"Warning: Kinship matrix not found at {kinship_file}")
            
    except Exception as e:
        print(f"Error in GEMMA analysis: {str(e)}")

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Process matrix data with normalization, plotting, and optional GEMMA analysis.')
    
    # Add arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Input file path (CSV format)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path for processed data (CSV format)')
    parser.add_argument('-p', '--plot', default='value_distribution_panel.png',
                        help='Output file path for plot (default: value_distribution_panel.png)')
    parser.add_argument('-d', '--dropped', default='dropped_columns.txt',
                        help='Output file path for dropped columns (default: dropped_columns.txt)')
    parser.add_argument('--pheno', default=None,
                        help='Phenotype file for GEMMA analysis (optional)')

    # Parse arguments
    args = parser.parse_args()
    
    try:
        # Process the matrix
        process_matrix_with_plots(
            input_file=args.input,
            output_file=args.output,
            plot_file_panel=args.plot,
            dropped_columns_file=args.dropped,
            pheno_file=args.pheno
        )
    except Exception as e:
        print(f"Error: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
