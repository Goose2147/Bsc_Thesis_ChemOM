import sys

def make_lmp_atoms_adj(lmp_opls_file, lmp_atoms_file):
    ''' Takes the lmp_opls and lmp_atoms files and recombines them to the right lmp_atoms
    files it will be named lmp_atoms_adj '''
    
    coeffs = write_coeff(lmp_opls_file)
    
    with open(lmp_atoms_file) as f:
        next(f)
        next(f)
        for line in f:
            if 'atoms' in line:
                atoms = int(line.split()[0])
            elif 'atom types' in line:
                atom_types = int(line.split()[0])
            elif 'bonds' in line:
                bonds = int(line.split()[0])
            elif 'bond types' in line:
                bond_types = int(line.split()[0])
            elif 'angles' in line:
                angles = int(line.split()[0])
            elif 'angle types' in line:
                angle_types = int(line.split()[0])
            elif 'dihedrals' in line:
                dihedrals = int(line.split()[0])
            elif 'dihedral types' in line:
                dihedral_types = int(line.split()[0])

    with open('lmp_atoms_adj', 'w') as out:
        with open(lmp_atoms_file) as f:
            count = 0
            for line in f:
                if count != 15:
                    out.write(line)
                    count += 1
                elif 'Masses' in line:
                    out.write(line + '\n')
                    line = next(f)
                    line = next(f)
                    line = line.split()
                    line[1] = '207.200'
                    line[5] = 'Pb'
                    line = ' '.join(line)
                    out.write(' ' + line + '\n')
                    for i in range(atom_types - 1):
                        out.write(next(f))
                    with open(coeffs) as c:
                        out.write('\n')
                        for line in c:
                            out.write(line)
                    with open(lmp_atoms_file) as f:
                        out.write('\n')
                        for line in f:
                            if 'Atoms' in line:
                                out.write(line)
                                for i in range(atoms + 1):
                                    out.write(next(f))
                            elif 'Bonds' in line:
                                out.write('\n' + line)
                                for i in range(bonds + angles + dihedrals + 7):
                                    out.write(next(f))

def write_coeff(lmp_opls_file):
    ''' Writes a coeffs file which is needed by make_lmp_atoms_adj; this file may be deleted
    afterwards. Also adds a zero in each row of the dihedrals, since this is required by the opls
    dihedral_style; these rows may have to be changed later, if another dihedral_style is being used '''
    
    with open(lmp_opls_file) as f:
        with open('coeffs.txt', 'w') as out:
            out.write('Bond Coeffs' + '\n\n')
            for line in f:
                if 'bond_coeff' in line:
                    out.write(line.split(maxsplit=1)[1])
                    #out.write('\n')
    
    with open(lmp_opls_file) as f:
        with open('coeffs.txt', 'a') as out:
            out.write('\nAngle Coeffs' + '\n\n')
            for line in f:
                if 'angle_coeff' in line:
                    out.write(line.split(maxsplit=1)[1])
                    #out.write('\n')
    
    with open(lmp_opls_file) as f:
        with open('coeffs.txt', 'a') as out:
            out.write('\nDihedral Coeffs' + '\n\n')
            for line in f:
                if 'dihedral_coeff' in line:
                    line = line.split()
                    line2 = ' '.join(line[:5]) + ' 0.0 ' + ' '.join(line[5:])
                    out.write(line2.split(maxsplit=1)[1] + '\n')
    
    return 'coeffs.txt'

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merging_lmp.py <lmp_opls_file> <lmp_atoms_file>")
        sys.exit(1)
    
    lmp_opls_file = sys.argv[1]
    lmp_atoms_file = sys.argv[2]
    
    make_lmp_atoms_adj(lmp_opls_file, lmp_atoms_file)
