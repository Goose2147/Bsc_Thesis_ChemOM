import re
import sys

def renumber_atom(atom_name, atom_counts, count_base):
    # Updated regex to correctly capture letters following the number
    match = re.match(r"([A-Za-z]+)(\d+)([A-Za-z]*)", atom_name)
    if match:
        prefix = match.group(1)
        num = match.group(2)
        suffix = match.group(3)
        base = prefix + num
        #print(count_base)
        if prefix not in atom_counts and base not in count_base:
            atom_counts[prefix] = 1
            count_base.append(base)
            #print(atom_name, 'finish 1')
        elif prefix in atom_counts and base not in count_base:
            atom_counts[prefix] += 1
            count_base.append(base)
            #print(atom_name, 'finish 2')
        else:
            #print(atom_name, 'finish 3')
            pass
            
        return f"{prefix}{atom_counts[prefix]}{suffix}"
    return atom_name

def process_mol2_file(input_filename):
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    atom_section = False
    atom_counts = {}
    new_lines = []
    count_base = []

    for line in lines:
        if line.startswith('@<TRIPOS>ATOM'):
            atom_section = True
            new_lines.append(line)
            continue
        if line.startswith('@<TRIPOS>BOND'):
            atom_section = False
        if atom_section:
            parts = line.split()
            if len(parts) > 1:
                atom_name = parts[1]
                new_atom_name = renumber_atom(atom_name, atom_counts, count_base)
                parts[1] = new_atom_name
                new_lines.append(" ".join(parts) + "\n")
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    return new_lines

def main():
    if len(sys.argv) != 2:
        print("Usage: python modify_mol2.py input.mol2")
        return

    input_filename = sys.argv[1]
    output_filename = 'modified_' + input_filename

    new_lines = process_mol2_file(input_filename)

    with open(output_filename, 'w') as file:
        file.writelines(new_lines)

    print(f"Modified file saved as {output_filename}")

if __name__ == "__main__":
    main()
