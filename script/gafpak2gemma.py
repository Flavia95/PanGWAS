import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### to do : add X and Y columns to teh output file 
### to do : input file form comand line 


def z_score_normalization(matrix):
    """
    Apply Z-score normalization to each column of the matrix.
    μj is the mean of the column, and σj is the standard deviation of the column.
    
    Parameters:
    - matrix: pd.DataFrame, the input matrix to be normalized
    
    Returns:
    - normalized: pd.DataFrame, the Z-score normalized matrix
    """
    # Convert all data to numeric, forcing errors to NaN
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    
    # Create a copy of the matrix to apply normalization
    normalized = matrix.copy()
    
    # Loop through each column and apply Z-score normalization
    for col in normalized.columns:
        # Calculate the mean and standard deviation, ignoring NaN values
        mean_col = np.nanmean(normalized[col])
        std_col = np.nanstd(normalized[col])
        
        # Apply Z-score normalization only if std is not zero to avoid division by zero
        if std_col > 0:
            normalized[col] = (normalized[col] - mean_col) / std_col
    
    return normalized

def min_max_scale(matrix, min_value=0, max_value=2):
    """
    Scale the values of the matrix to be within the range [min_value, max_value].
    
    Parameters:
    - matrix: pd.DataFrame, the input matrix to be scaled
    - min_value: float, the minimum value of the scaled range (default 0)
    - max_value: float, the maximum value of the scaled range (default 2)
    
    Returns:
    - scaled: pd.DataFrame, the min-max scaled matrix
    """
    # Apply min-max scaling to each column
    scaled = matrix.copy()
    for col in scaled.columns:
        col_min = scaled[col].min()
        col_max = scaled[col].max()
        
        # Avoid division by zero in case of constant columns
        if col_max - col_min > 0:
            scaled[col] = (scaled[col] - col_min) / (col_max - col_min) * (max_value - min_value) + min_value
    
    return scaled

def process_matrix_with_plots(input_file, output_file, plot_file_panel, dropped_columns_file='dropped_columns.txt'):
    # Read the matrix from the input CSV file

    # Assume the first column is row names
    matrix = pd.read_csv(input_file, index_col=0, skiprows=1)
    print(f'Matrix size input: {matrix.shape}')


    # Ensure all data is numeric, drop non-numeric data
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    
    # Identify columns that are all NaN
    nan_columns = matrix.columns[matrix.isna().all()].tolist()
    
    # Write these column names to a file
    with open(dropped_columns_file, 'w') as file:
        for col in nan_columns:
            file.write(f"{col}\n")
    
    # Drop rows and columns that are completely NaN (non-numeric or empty)
    matrix = matrix.dropna(axis=0, how='all')  # Drop rows that are all NaN
    matrix = matrix.dropna(axis=1, how='all')  # Drop columns that are all NaN

    # Create subplots (1 row, 3 columns)
    fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    
    # Plot the distribution of original values from the input
    axes[0].hist(matrix.values.flatten(), bins=30, alpha=0.7, color='orange', edgecolor='black')
    axes[0].set_xlabel('Value')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Original Input Values')
    axes[0].set_yscale('log')  # Log scale on y-axis to handle wide range of values
    axes[0].grid(True)
    
    # Apply Z-score normalization
    z_normalized_matrix = z_score_normalization(matrix)
    print(f'Matrix size after Z-score normalization: {z_normalized_matrix.shape}')

    # Plot the distribution of values after Z-score normalization
    axes[1].hist(z_normalized_matrix.values.flatten(), bins=30, alpha=0.7, color='blue', edgecolor='black')
    axes[1].set_xlabel('Value')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('After Z-Score Normalization')
    axes[1].grid(True)
    
    # Apply Min-Max scaling to Z-score normalized data to range [0, 2]
    final_normalized_matrix = min_max_scale(z_normalized_matrix, 0, 2)
    print(f'Matrix size after rescaling: {final_normalized_matrix.shape}')

########TO DO Add two new columns 'X' and 'Y' after the first column using pd.concat for better performance
    #final_normalized_matrix(1, '', 'X')
   # matrix.insert(2, 'Y','Y')

    # Save the final normalized matrix to the output CSV file
    final_normalized_matrix.to_csv(output_file , float_format='%.5f')

    # Plot the distribution of values after scaling to [0, 2]
    axes[2].hist(final_normalized_matrix.values.flatten(), bins=30, alpha=0.7, color='green', edgecolor='black')
    axes[2].set_xlabel('Value')
    axes[2].set_ylabel('Frequency')
    axes[2].set_title('After Scaling to [0, 2]')
    axes[2].grid(True)
    
    # Adjust layout and save the combined plot
    plt.tight_layout()
    plt.savefig(plot_file_panel)
    #plt.show()

# Example usage:

input_file = "data/matrcov149_tr.geno.rpn"  # Replace with the path to your input file
output_file = "final_normalized_matrix.csv"  # Replace with the desired output file path
plot_file_panel = "value_distribution_panel.png"  # Combined plot file

process_matrix_with_plots(input_file, output_file, plot_file_panel)
