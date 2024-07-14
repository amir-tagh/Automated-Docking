import sys
import re
import matplotlib.pyplot as plt




input_filename = sys.argv[1]  # Replace with your file name
pattern = r'\b(Grid_Score)\b'
output_filename = 'output.txt'
data = []


# Example usage
output_filename_1 = 'scores.dat'
output_filename_2 = 'output.mol2'
start_pattern = '@<TRIPOS>MOLECULE'
end_pattern = '1 ZINC'

lines_to_save = []
between_patterns = False


#extract the mol2 file for highest docked scores
def extract_lines_between_patterns(input_filename, output_filename, start_pattern, end_pattern):
    lines_to_save = []
    between_patterns = False
    
    with open(input_filename, 'r',errors='ignore') as f_in:
        for line in f_in:
            if line.strip() == start_pattern:
                between_patterns = True
                lines_to_save.append(line)
            elif line.strip() == end_pattern:
                between_patterns = False
                lines_to_save.append(line)
            elif between_patterns:
                lines_to_save.append(line)
    
    with open(output_filename_2, 'w') as f2_out:
        for line in lines_to_save:
            f2_out.write(line)



# Open file and read line by line
with open(input_filename, 'r', errors='ignore') as f, open(output_filename_1, 'w') as f_out:
    for line in f:
        # Process each line (for demonstration, just print it)
        #print(line.strip())  # strip() removes the newline character at the end
        if re.search(pattern, line):
            #print(line.strip())  # Print the matching line (strip() removes newline)
            columns = line.split()
            if len(columns) >= 3:
                if float(columns[2]) <= -30.0:
                    third_column = columns[2]
                    print(columns[2])
                    f_out.write(third_column + '\n')
                    value = float(columns[2])
                    data.append(value)

# Plotting the histogram
plt.hist(data, bins=10, edgecolor='black')  # Adjust bins as needed
plt.xlabel('Docking score')
plt.ylabel('Frequency')
plt.gca().invert_xaxis()
plt.title('Histogram of Docking score')
plt.grid(True)
plt.savefig('docking-hist.png')
#plt.show()             

