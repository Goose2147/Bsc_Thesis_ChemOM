import sys
from rdkit import Chem
from collections import defaultdict
import re


def read_mol2_file(file_path):
    mol = Chem.MolFromMol2File(file_path, removeHs=False)
    if mol is None:
        raise ValueError("Failed to read the mol2 file.")
    return mol

def get_atom_names(file_path):
    atom_names = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        atom_section = False
        for line in lines:
            if line.startswith('@<TRIPOS>ATOM'):
                atom_section = True
                continue
            if line.startswith('@<TRIPOS>BOND'):
                break
            if atom_section:
                parts = line.split()
                if len(parts) > 1:
                    atom_idx = int(parts[0]) - 1  # mol2 indices are 1-based, RDKit uses 0-based
                    atom_name = parts[1]
                    atom_names[atom_idx] = atom_name
    return atom_names

def get_hydrogen_reduction_map(atom_names):
    hydrogen_map = {}
    for idx, name in atom_names.items():
        if name.startswith('H'):
            match = re.match(r'([A-Za-z]+)(\d+)([A-Za-z]+)', name)
            if match:
                prefix = match.group(1)
                number = match.group(2)
                #letter = match.group(3)
                base_name = prefix + number  # Extract the base name (e.g., "H2A" -> "H")
            else:
                base_name = name
            if base_name not in hydrogen_map:
                hydrogen_map[idx] = base_name
        else:
            hydrogen_map[idx] = name
    return hydrogen_map

def get_bonds(mol, hydrogen_map):
    bonds = []
    for bond in mol.GetBonds():
        atom1 = hydrogen_map.get(bond.GetBeginAtomIdx())
        atom2 = hydrogen_map.get(bond.GetEndAtomIdx())
        bond = (atom1, atom2)
        bonds.append(bond) if bond not in bonds else None
    return bonds

def get_angles(mol, hydrogen_map):
    angles = []
    for atom in mol.GetAtoms():
        neighbors = [hydrogen_map.get(neighbor.GetIdx(), neighbor.GetIdx()) for neighbor in atom.GetNeighbors()]
        if len(neighbors) < 2:
            continue
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                ang = (neighbors[i], hydrogen_map.get(atom.GetIdx(), atom.GetIdx()), neighbors[j])
                if ang not in angles:
                    angles.append(ang)
    return angles

def get_dihedrals(mol, hydrogen_map):
    dihedrals = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        neighbors1 = [hydrogen_map.get(neighbor.GetIdx(), neighbor.GetIdx()) for neighbor in atom1.GetNeighbors() if neighbor.GetIdx() != atom2.GetIdx()]
        neighbors2 = [hydrogen_map.get(neighbor.GetIdx(), neighbor.GetIdx()) for neighbor in atom2.GetNeighbors() if neighbor.GetIdx() != atom1.GetIdx()]
        for n1 in neighbors1:
            for n2 in neighbors2:
                dh = (n1, hydrogen_map.get(atom1.GetIdx(), atom1.GetIdx()), hydrogen_map.get(atom2.GetIdx(), atom2.GetIdx()), n2)
                if dh not in dihedrals:
                    dihedrals.append(dh)
    return dihedrals

def assign_MATCH_label(hydrogen_map):
    labeled_match_atoms = {}
    for idx, name in hydrogen_map.items():
        match_label = input(f'MATCH atom type assigned to atom {name}: ')
        labeled_match_atoms[name] = match_label
    print('Copy this dictionary to skip the assignment next time: ', labeled_match_atoms)
    return labeled_match_atoms

def read_fixed_par_file(fixed_par_file):
    with open(fixed_par_file, 'r') as file:
        lines = file.readlines()
        header = []
        atoms = {}
        bonds = {}
        angles = {}
        dihedrals = {}
        cutoffs = {}
        section = 'header'
        for line in lines:
            line = line.strip()
            if line.startswith('# one body'):
                section = 'atoms'
            elif line.startswith('# bonds'):
                section = 'bonds'
            elif line.startswith('# angles'):
                section = 'angles'
            elif line.startswith('# dihedrals'):
                section = 'dihedrals'
            elif line.startswith('# cutoffs'):
                section = 'cutoffs'            
            elif line and section == 'header':
               header.append(line)
            elif line and section == 'atoms':
               parts = line.split()
               atoms[parts[0]] = (parts[1], parts[2])
            elif line and section == 'bonds':
               parts = line.split()
               bonds[parts[0]] = (parts[1], parts[2])
            elif line and section == 'angles':
               parts = line.split()
               angles[parts[0]] = (parts[1], parts[2])
            elif line and section == 'dihedrals':
               parts = line.split()
               dihedrals[parts[0]] = (parts[1], parts[2], parts[3])
            elif line and section == 'cutoffs':
               parts = line.split()
               cutoffs[parts[0]] = parts[1]        
    return header, atoms, bonds, angles, dihedrals, cutoffs

def main(file_path, fixed_par_file):
    mol = read_mol2_file(file_path)
    atom_names = get_atom_names(file_path)
    hydrogen_map = get_hydrogen_reduction_map(atom_names)
    
    # Use assign_MATCH_label(hydrogen_map) the first time, copy the dictionary to skip the assignment.
    labeled_match_atoms = assign_MATCH_label(hydrogen_map)
    # labeled_match_atoms = {'O1': 'O301', 'N1': 'N3P3', 'H1': 'HGP2', 'C1': 'C261', 'C2': 'C261', 'H2': 'HG61', 'C3': 'C261', 'H3': 'HG61', 'C4': 'C261', 'H4': 'HG61', 'C5': 'C261', 'C6': 'C261', 'H5': 'HG61', 'C7': 'C261', 'H6': 'HG61', 'C8': 'C261', 'H7': 'HG61', 'C9': 'C261', 'H8': 'HG61', 'C10': 'C261', 'C11': 'C321', 'H9': 'HGA2', 'C12': 'C324', 'H10': 'HGA2'}
    
    bonds = get_bonds(mol, hydrogen_map)
    angles = get_angles(mol, hydrogen_map)
    dihedrals = get_dihedrals(mol, hydrogen_map)

    #print("Bonds:")
    labeled_match_bonds = {}
    for b1, b2 in bonds:
        bond_string = f'{b1}-{b2}'
        b1_match = labeled_match_atoms[b1]  
        b2_match = labeled_match_atoms[b2]
        labeled_match_bonds[bond_string] = f'{b1_match}-{b2_match}'
        #print(bond_string, ' : ', labeled_match_bonds[bond_string])
    #print(labeled_match_bonds)
    
    #print("\nAngles:")
    labeled_match_angles = {}
    for a1, a2, a3 in angles:
        angle_string = f'{a1}-{a2}-{a3}'
        a1_match = labeled_match_atoms[a1]  
        a2_match = labeled_match_atoms[a2]
        a3_match = labeled_match_atoms[a3]
        labeled_match_angles[angle_string] = f'{a1_match}-{a2_match}-{a3_match}'
        #print(angle_string, ' : ', labeled_match_angles[angle_string])
    #print(labeled_match_angles)
    
    #print("\nDihedrals:")
    labeled_match_dihedrals = {}
    for d1, d2, d3, d4 in dihedrals:
        dh_string = f'{d1}-{d2}-{d3}-{d4}'
        d1_match = labeled_match_atoms[d1]  
        d2_match = labeled_match_atoms[d2]
        d3_match = labeled_match_atoms[d3]        
        d4_match = labeled_match_atoms[d4]
        labeled_match_dihedrals[dh_string] = f'{d1_match}-{d2_match}-{d3_match}-{d4_match}'
        #print(dh_string, ' : ', labeled_match_dihedrals[dh_string])
    #print(labeled_match_dihedrals)
    
    
    # Read par file with atom types from MATCH and obtain all sections
    par_header, par_atoms, par_bonds, par_angles, par_dihedrals, par_cutoffs = read_fixed_par_file(fixed_par_file)
    
    # Write the output content to a file
    with open('forase.par', 'w') as outfile:
        # Header
        for header_line in par_header:
            outfile.write(header_line + '\n')
        
        # Write One body - LJ parameters - Charges
        outfile.write('# one body - LJ-parameters - charges\n')
        cage1 = input('Inorganic atom_name 1 (e.g. P1): ')
        cage2 = input('Inorganic atom_name 2 (e.g. I1): ')
        outfile.write(f'{cage1:<11} {0.000000:<9.6f} {0.000000:<10.6f} {0.000000:<10.6f}\n')
        outfile.write(f'{cage2:<11} {0.000000:<9.6f} {0.000000:<10.6f} {0.000000:<10.6f}\n')
        for par_atom, match_atom in labeled_match_atoms.items(): 
            param1, param2 = par_atoms[match_atom]
            outfile.write(f'{par_atom:<10} {float(param1):<10.6f} {float(param2):<10.6f} {0.000000:<10.6f}\n')
        
        # Write Bonds
        outfile.write('\n# bonds\n')
        for par_bond, match_bond in labeled_match_bonds.items():
            if match_bond in par_bonds:
                param1, param2 = par_bonds[match_bond]
            else:
                match_bond = '-'.join(match_bond.split('-')[::-1])
                param1, param2 = par_bonds[match_bond]
            outfile.write(f"{par_bond:<15} {float(param1):<10.2f} {float(param2):<10.4f}\n")
        
        # Write Angles
        outfile.write('\n# angles\n')
        for par_angle, match_angle in labeled_match_angles.items():
            if match_angle in par_angles:
                param1, param2 = par_angles[match_angle]
            else:
                match_angle = '-'.join(match_angle.split('-')[::-1])
                param1, param2 = par_angles[match_angle]
            outfile.write(f"{par_angle:<20} {float(param1):<10.2f} {float(param2):<10.2f}\n")
            
        # Write Dihedrals
        outfile.write('\n# dihedrals\n')
        for par_dh, match_dh in labeled_match_dihedrals.items():
            if match_dh in par_dihedrals:
                param1, param2, param3 = par_dihedrals[match_dh]
            else:
                match_dh = '-'.join(match_dh.split('-')[::-1])
                param1, param2, param3 = par_dihedrals[match_dh]
            outfile.write(f"{par_dh:<25} {float(param1):<10.4f} {float(param2):<10.0f} {float(param3):<10.2f}\n")

        # Write Cutoffs
        outfile.write('\n# cutoffs\n')
        cutoff_cage = cage1 + '-' + cage2
        if cutoff_cage in par_cutoffs:
            outfile.write(f'{cutoff_cage:<22} {float(par_cutoffs[cutoff_cage]):<10.6f}\n')
        else:
            cutoff_cage = '-'.join(cutoff_cage.split('-')[::-1])
            outfile.write(f'{cutoff_cage:<22} {float(par_cutoffs[cutoff_cage]):<10.6f}\n')
            
        for par_cut, match_cut in labeled_match_bonds.items():  # We use bonds for cutoffs, they are the same
            if match_cut in par_cutoffs:
                param1 = par_cutoffs[match_cut]
            else: 
                match_cut = '-'.join(match_cut.split('-')[::-1])
                param1 = par_cutoffs[match_cut]
            outfile.write(f'{par_cut:<22} {float(param1):<10.6f}\n')
            
        print('Data written to file: forase.par')
            
  
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <mol2_file_path> <fixed_match_par_file>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    fixed_par_file = sys.argv[2]
    main(file_path, fixed_par_file)
