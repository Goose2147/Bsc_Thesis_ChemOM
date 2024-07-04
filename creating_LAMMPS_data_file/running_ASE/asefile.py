import ase
from ase.visualize import view
from opls import OPLSff, OPLSStructure

s = OPLSStructure('forase.xyz') # input xyz file
#view(s) # view with real elements
elements = { 'P1' : 'Pb' , 'I1' : 'I', 'O1' : 'O', 'N1' : 'N', 'C1' : 'C', 'C2' : 'C', 'C3' : 'C', 'C4' : 'C', 'C5' : 'C', 'C6' : 'C', 'C7' : 'C', 'C8' : 'C', 'C9' : 'C', 'C10' : 'C', 'C11' : 'C', 'C12' : 'C', 'H1' : 'H', 'H2' : 'H', 'H3' : 'H', 'H4' : 'H', 'H5' : 'H', 'H6' : 'H', 'H7' : 'H', 'H8' : 'H', 'H9' : 'H', 'H10' : 'H'}	# for NOE
#view(s.colored(elements)) # view with fake elements

opls = OPLSff('forase.par') # input par file
opls.write_lammps(s, prefix='lmp')
