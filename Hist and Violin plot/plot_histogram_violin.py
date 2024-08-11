import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from concurrent.futures import ProcessPoolExecutor

def plot_histogram(data, column, output_file,font_size=24,tick_font_size=16):
    plt.figure(figsize=(10, 6))
    plt.hist(data[column], bins=30, color='blue', alpha=0.7, edgecolor='black')
    plt.title(f'Histogram of Negative {column}',fontsize=font_size)
    plt.xlabel(column,fontsize=font_size)
    plt.ylabel('Frequency',fontsize=font_size)
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()

    # Set tick parameters
    plt.tick_params(axis='both', which='major', labelsize=tick_font_size)

    plt.grid(True)
    plt.savefig(output_file)
    plt.close()

def plot_violin(data_chunk, column, output_file, font_size=14, tick_font_size=12):
    plt.figure(figsize=(10, 6))
    sns.set_context("notebook", rc={"axes.labelsize": font_size, "xtick.labelsize": tick_font_size, "ytick.labelsize": tick_font_size})
    sns.violinplot(y=data_chunk[column], color='skyblue')
    plt.title(f'Violin Plot of Negative {column}', fontsize=font_size)
    plt.ylabel(column, fontsize=font_size)
    
    plt.savefig(output_file)
    plt.close()

def process_column(column, data, histogram_output, violin_output):
    # Filter data for negative numbers only
    negative_data = data[data[column] < -10.0]

    if negative_data.empty:
        print(f"No negative numbers found in column '{column}'.")
        return

    # Plot histogram
    plot_histogram(negative_data, column, histogram_output)
    print(f"Histogram of negative numbers saved as {histogram_output}")

    # Plot violin plot
    plot_violin(negative_data, column, violin_output)
    print(f"Violin plot of negative numbers saved as {violin_output}")

def main(input_file, column, histogram_output, violin_output):
    # Read the data
    data = pd.read_csv(input_file, delim_whitespace=True)  # Reading space-separated file

    # Check if the column exists
    if column not in data.columns:
        print(f"Error: Column '{column}' not found in the input file.")
        print(f"Available columns are: {list(data.columns)}")
        return

    # Process the column in parallel
    with ProcessPoolExecutor() as executor:
        executor.submit(process_column, column, data, histogram_output, violin_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate histogram and violin plot from a space-separated file for negative numbers only, in parallel.")
    parser.add_argument("input_file", help="Path to the input space-separated file.")
    parser.add_argument("column", help="Column name to plot.")
    parser.add_argument("histogram_output", help="Output file name for the histogram plot.")
    parser.add_argument("violin_output", help="Output file name for the violin plot.")

    args = parser.parse_args()

    main(args.input_file, args.column, args.histogram_output, args.violin_output)

