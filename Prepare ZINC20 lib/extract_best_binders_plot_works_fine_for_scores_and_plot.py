#######################
#scripts extracts the mol2 files from the Dock 6.0 output\
#based on a threshold set by Grid Score\
#plots the histogram of score distribution\
#Amirhossein Taghavi
#UF Scripps
#07/24/2024


import re
import matplotlib.pyplot as plt

def parse_mol2_file(filename):
    """Parse the MOL2 file and extract Grid_Score values."""
    with open(filename, 'r', encoding='latin-1') as file:  # Use 'latin-1' encoding
        content = file.read()
    
    # Split content into blocks of molecule data
    mol2_blocks = content.split('@<TRIPOS>MOLECULE\n')[1:]  # Skip header part
    structure_data = []
    
    # Regex pattern to find Grid_Score and other data
    score_pattern = re.compile(r'Grid_Score:\s*(-?\d+\.\d+)')
    
    for index, block in enumerate(mol2_blocks):
        # Extract score from the block
        match = score_pattern.search(block)
        if match:
            score = float(match.group(1))
            structure_data.append((score, index))
    
    return structure_data, mol2_blocks

def filter_mol2_by_score(input_file, output_file, threshold, scores_file):
    """Filter structures based on Grid_Score threshold, save to output file, and save scores to scores_file."""
    scores, mol2_blocks = parse_mol2_file(input_file)
    
    print(f"Total structures read: {len(mol2_blocks)}")
    print(f"Total structures with scores: {len(scores)}")

    # Open the output file for writing
    with open(output_file, 'w') as file:
        # Open scores file for writing
        with open(scores_file, 'w') as score_file:
            # Write header for scores file
            score_file.write('Score Index\n')
            
            # Iterate over scores and write corresponding blocks to output file
            for score, index in scores:
                # Write scores to the scores file
                score_file.write(f"{score} {index}\n")
                
                if score <= threshold:
                    print(f"Score {score} is below or equal to threshold {threshold}, saving structure {index + 1}.")
                    file.write('@<TRIPOS>MOLECULE\n')
                    file.write(mol2_blocks[index])
                    file.write("\n")
    
    print(f"Data saved to {output_file}")
    print(f"Scores saved to {scores_file}")

def plot_histogram(scores_file):
    """Plot a histogram of the scores."""
    scores = []
    with open(scores_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            score, _ = map(float, line.split())
            scores.append(score)
    
    plt.figure(figsize=(10, 6))
    plt.hist(scores, bins=30, edgecolor='black')
    plt.title('Histogram of Grid_Scores')
    plt.xlabel('Grid_Score')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

def main():
    input_file = 'virtual.out_scored.mol2'
    output_file = 'output.mol2'
    threshold = -50.0
    scores_file = 'scores.dat'
    
    print(f"Filtering and saving structures from {input_file} to {output_file} with Grid_Score <= {threshold}")
    filter_mol2_by_score(input_file, output_file, threshold, scores_file)
    print(f"Plotting histogram of scores from {scores_file}")
    plot_histogram(scores_file)

if __name__ == "__main__":
    main()

