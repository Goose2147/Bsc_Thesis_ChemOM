import sys
from typing import List, Tuple, Dict

def read_input_file(input_filename: str) -> Tuple[
    List[str], 
    Dict[str, Tuple[str, str]], 
    Dict[str, Tuple[str, str]], 
    Dict[str, Tuple[str, str]], 
    Dict[str, List[Tuple[str, str, str]]], 
    float, 
    List[str]
]:
    """
    Reads a .par file from SwissParam and formats it into a ASE-like .par file.
    
    Args:
        input_filename (str): The name of the input file to read.

    Returns:
        Tuple: Contains header_lines (List[str]), atom_LJ (Dict[str, Tuple[str, str]]), 
               bonds (Dict[str, Tuple[str, str]]), angles (Dict[str, Tuple[str, str]]), 
               dihedrals (Dict[str, List[Tuple[str, str, str]]]), cage_cutoff (float), 
               and cage_atoms (List[str]).
    """
    # Ask the user about the composition of the organic plane
    organic_plane = input("Is the organic plane comprised of Pb and I? (y/n): ").strip().lower()
    
    if organic_plane == 'y':
        cage_atoms = ['Pb', 'I']
        cage_bond_length = 5  # No need to ask for bond length
    else:
        # Ask the user for the two atoms comprising the inorganic plane and the bond length
        cage_atom1 = input("Please specify the first atom of the inorganic plane: ").strip()
        cage_atom2 = input("Please specify the second atom of the inorganic plane: ").strip()
        cage_bond_length = float(input(f"Please specify the bond length between {cage_atom1} and {cage_atom2}: ").strip())
        cage_atoms = [cage_atom1, cage_atom2]
    
    # Calculate cage cutoff
    if cage_bond_length - int(cage_bond_length) > 0.1:
        cage_cutoff = int(cage_bond_length) + 1
    else:
        cage_cutoff = int(cage_bond_length) + 0.5
    
    header_lines = []
    atom_LJ: Dict[str, Tuple[str, str]] = {}
    bonds: Dict[str, Tuple[str, str]] = {}
    angles: Dict[str, Tuple[str, str]] = {}
    dihedrals: Dict[str, List[Tuple[str, str, str]]] = {}
    
    with open(input_filename, 'r') as infile:
        section = None
        for line in infile:
            line = line.strip()
            if line.startswith('*'):
                header_lines.append('#' + line[1:])
            elif line.startswith('ATOMS'):
                section = 'ATOMS'
            elif line.startswith('BONDS'):
                section = 'BONDS'
            elif line.startswith('ANGLES'):
                section = 'ANGLES'
            elif line.startswith('DIHEDRALS'):
                section = 'DIHEDRALS'
            elif line.startswith('IMPROPER'):
                section = 'IMPROPER'
            elif line.startswith('NONBONDED'):
                section = 'NONBONDED'
            elif section == 'NONBONDED' and line:
                if line.startswith('cutnb'):
                    continue
                else:
                    cols = line.split()
                    atom_LJ[cols[0]] = (cols[2], cols[3])
            elif section == 'BONDS' and line:
                cols = line.split()
                bond_type = cols[0] + '-' + cols[1]
                bonds[bond_type] = (cols[2], cols[3])
            elif section == 'ANGLES' and line:
                cols = line.split()
                angle_type = cols[0] + '-' + cols[1] + '-' + cols[2]
                angles[angle_type] = (cols[3], cols[4])
            elif section == 'DIHEDRALS' and line:
                cols = line.split()
                dihedral_type = cols[0] + '-' + cols[1] + '-' + cols[2] + '-' + cols[3]
                if dihedral_type not in dihedrals:
                    dihedrals[dihedral_type] = []
                dihedrals[dihedral_type].append((cols[4], cols[5], cols[6]))

    return header_lines, atom_LJ, bonds, angles, dihedrals, cage_cutoff, cage_atoms

def modify_and_write_output(
    header_lines: List[str], 
    atom_LJ: Dict[str, Tuple[str, str]], 
    bonds: Dict[str, Tuple[str, str]], 
    angles: Dict[str, Tuple[str, str]], 
    dihedrals: Dict[str, List[Tuple[str, str, str]]], 
    cage_cutoff: float, 
    cage_atoms: List[str], 
    output_filename: str
) -> None:
    """
    Modifies the extracted data and writes it to an output file.
    
    Args:
        header_lines (List[str]): Header lines from the input file.
        atom_LJ (Dict[str, Tuple[str, str]]): LJ parameters for atoms.
        bonds (Dict[str, Tuple[str, str]]): Bond parameters.
        angles (Dict[str, Tuple[str, str]]): Angle parameters.
        dihedrals (Dict[str, List[Tuple[str, str, str]]]): Dihedral parameters.
        cage_cutoff (float): Calculated cutoff for the cage atoms.
        cage_atoms (List[str]): List of cage atoms.
        output_filename (str): The name of the output file to write.
    """
    cage_atom1_type = f"{cage_atoms[0][0]}1"
    cage_atom2_type = f"{cage_atoms[1][0]}1"
    
    with open(output_filename, 'w') as outfile:
        # Write header lines
        for line in header_lines:
            outfile.write(line + '\n')

        # Write LJ parameters section
        outfile.write('# one body - LJ-parameters\n')
        outfile.write(f'{cage_atom1_type:<10} {0.000000:<10.6f} {0.000000:<10.6f}\n')
        outfile.write(f'{cage_atom2_type:<10} {0.000000:<10.6f} {0.000000:<10.6f}\n')
        for atom, params in atom_LJ.items():
            outfile.write(f'{atom:<10} {float(params[0]):<10.6f} {float(params[1]):<10.6f}\n')
            
        # Write Bonds
        outfile.write('\n# bonds\n')
        for bond_type, bond_vals in bonds.items():
            outfile.write(f"{bond_type:<15} {bond_vals[0]:<10} {bond_vals[1]:<10}\n")
        
        # Write Angles
        outfile.write('\n# angles\n')
        for angle_type, angle_vals in angles.items():
            outfile.write(f"{angle_type:<20} {angle_vals[0]:<10} {angle_vals[1]:<10}\n")
            
        # Write Dihedrals
        outfile.write('\n# dihedrals\n')
        for dihedral_type, dihedral_lst in dihedrals.items():
            for dihedral_vals in dihedral_lst:
                outfile.write(f"{dihedral_type:<25} {dihedral_vals[0]:<10} {dihedral_vals[1]:<10} {dihedral_vals[2]:<10}\n")

        # Write Cutoffs
        outfile.write('\n# cutoffs\n')
        outfile.write(f'{cage_atom1_type}-{cage_atom2_type:<22} {cage_cutoff:<10.6f}\n')
        for bond_type, bond_vals in bonds.items():
            bond_length = float(bond_vals[1])
            if bond_length - int(bond_length) > 0.1:
                bond_cutoff = int(bond_length) + 1
            else:
                bond_cutoff = int(bond_length) + 0.5
            outfile.write(f"{bond_type:<25} {bond_cutoff:<10.6f}\n")

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
    
    # Get input and output filenames from command-line arguments
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    # Read data from the input file
    header_lines, atom_LJ, bonds, angles, dihedrals, cage_cutoff, cage_atoms = read_input_file(input_filename)
    
    # Modify and write the output file
    modify_and_write_output(header_lines, atom_LJ, bonds, angles, dihedrals, cage_cutoff, cage_atoms, output_filename)
    
    print(f"File has been read from {input_filename}, modified, and written to {output_filename}.")
