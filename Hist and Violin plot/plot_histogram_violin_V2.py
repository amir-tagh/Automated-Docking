import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def plot_histogram(data, column, output_file, font_size=14, tick_font_size=12):
    plt.figure(figsize=(10, 6))
    plt.hist(data[column], bins=30, color='blue', alpha=0.7, edgecolor='black')
    plt.title(f'Histogram of Negative {column}', fontsize=font_size)
    plt.xlabel(column, fontsize=font_size)
    plt.ylabel('Frequency', fontsize=font_size)
    
    # Set tick parameters
    plt.tick_params(axis='both', which='major', labelsize=tick_font_size)
    
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to ensure labels fit
    plt.savefig(output_file)
    plt.close()

def plot_violin(data, column, output_file, font_size=14, tick_font_size=12):
    plt.figure(figsize=(10, 6))
    sns.set_context("notebook", rc={"axes.labelsize": font_size, "xtick.labelsize": tick_font_size, "ytick.labelsize": tick_font_size})
    sns.violinplot(y=data[column], color='skyblue')
    plt.title(f'Violin Plot of Negative {column}', fontsize=font_size)
    plt.ylabel(column, fontsize=font_size)
    
    plt.tight_layout()  # Adjust layout to ensure labels fit
    plt.savefig(output_file)
    plt.close()

def main(input_file, column, output_dir, font_size=14, tick_font_size=12):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the data
    data = pd.read_csv(input_file, sep='\s+')  # Reading space-separated file

    # Check if the column exists
    if column not in data.columns:
        print(f"Error: Column '{column}' not found in the input file.")
        print(f"Available columns are: {list(data.columns)}")
        return

    # Filter data for negative numbers only
    negative_data = data[data[column] < 0]

    if negative_data.empty:
        print(f"No negative numbers found in the data.")
        return

    # Create file names for the plots
    histogram_output = os.path.join(output_dir, f"{column}_histogram.png")
    violin_output = os.path.join(output_dir, f"{column}_violin.png")

    # Plot histogram
    plot_histogram(negative_data, column, histogram_output, font_size, tick_font_size)
    print(f"Histogram of negative numbers saved as {histogram_output}")

    # Plot violin plot
    plot_violin(negative_data, column, violin_output, font_size, tick_font_size)
    print(f"Violin plot of negative numbers saved as {violin_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and save histogram and violin plot for a single column from a space-separated file.")
    parser.add_argument("input_file", help="Path to the input space-separated file.")
    parser.add_argument("column", help="Column name to plot.")
    parser.add_argument("output_dir", help="Directory to save the output plots.")
    parser.add_argument("--font_size", type=int, default=14, help="Font size for x-label and y-label. Default is 14.")
    parser.add_argument("--tick_font_size", type=int, default=12, help="Font size for x-ticks and y-ticks. Default is 12.")

    args = parser.parse_args()

    main(args.input_file, args.column, args.output_dir, args.font_size, args.tick_font_size)

