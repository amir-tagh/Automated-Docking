import sys
import re
import matplotlib.pyplot as plt


input_filename = sys.argv[1]  # Replace with your file name
output_filename = 'best_binder.mol2'

def save_lines_between_patterns(input_filename, output_filename, start_condition, threshold, start_save, end_save):
    lines_to_save = []
    inside_block = False
    found_start_condition = False

    with open(input_filename, 'r', errors='ignore') as f_in:
        for line in f_in:
            if start_condition in line:
                found_start_condition = True
                # Extract the second column value
                try:
                    second_column_value = float(line.split()[2])
                except (IndexError, ValueError):
                    continue  # Skip lines without a second column or invalid float conversion

                if second_column_value <= threshold:
                    print(second_column_value)
                    inside_block = True
                    lines_to_save.append(line)  # Save the line with start_condition if it meets the criteria
                else:
                    inside_block = False
                continue

            if inside_block:
                lines_to_save.append(line)
                if end_save in line:
                    inside_block = False
                    lines_to_save.append('\n')  # Add a new line to separate blocks

    if found_start_condition:
        with open(output_filename, 'w') as f_out:
            for line_to_save in lines_to_save:
                f_out.write(line_to_save)

# Example usage
start_condition = 'Grid_Score'
threshold = (float(-40.0))  # Example threshold value
start_save = '@<TRIPOS>MOLECULE'
end_save = '1 ZINC'

save_lines_between_patterns(input_filename, output_filename, start_condition, threshold, start_save, end_save)

