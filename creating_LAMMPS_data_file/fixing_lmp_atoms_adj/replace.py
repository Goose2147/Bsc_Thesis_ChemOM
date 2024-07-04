def replace_coordinates(original_lines, replacement_lines):
    updated_lines = []
    for orig_line, replace_line in zip(original_lines, replacement_lines):
        # Extracting parts from original line
        parts = orig_line.split()
        # Extracting coordinates from replacement line
        coords = replace_line.split()[1:4]
        # Creating updated line with original formatting and replaced coordinates
        len_space = 7 - len(parts[0]) - 1
        ind = len_space * ' ' + str(parts[0])
        len_space2 = 5 - len(parts[2]) - 1
        ind2 = len_space2 * ' ' + str(parts[2])
        updated_line = f"{ind}   {parts[1]}{ind2} {parts[3]} {float(coords[0]):.9f} {float(coords[1]):.9f} {float(coords[2]):.9f} {' '.join(parts[7:])}\n"
        updated_lines.append(updated_line)
    return updated_lines

# Read original text file
original_file = 'original.txt'
with open(original_file, 'r') as f:
    original_lines = f.readlines()

# Read replacement text file
replacement_file = 'replacement.txt'
with open(replacement_file, 'r') as f:
    replacement_lines = f.readlines()

# Replace coordinates in the original text
updated_lines = replace_coordinates(original_lines, replacement_lines)

# Write updated text to a new file
output_file = 'updated.txt'
with open(output_file, 'w') as f:
    f.writelines(updated_lines)

print(f"Updated file '{output_file}' has been created.")
