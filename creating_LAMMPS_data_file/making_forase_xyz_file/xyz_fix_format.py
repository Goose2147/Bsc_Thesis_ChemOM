import numpy as np
import sys

" The file takes 5 inputs "
" Input 1: The name of the input file "
" Input 2: The box limit in x direction "
" Input 3: The box limit in y direction "
" Input 4: The box limit in z direction "
" Input 5: The name of the output file "

def main(argv):
	"------------------"
	" -----Inputs----- "
	"------------------"

	dump = open(argv[1], 'r')
	
	data_string = dump.readlines()
	nLines = len(data_string)

	dump.close()

	x_cage = float(argv[2])
	y_cage = float(argv[3])
	z_cage = float(argv[4])
	
	"-----------------------------------"
	"-----Splitting-up-the-elements-----"
	"-----------------------------------"

	" Current build assumes Pb, I, C, N, H and O in the system "

	Pb = []
	I = []
	C = []
	N = []
	H = []
	O = []
	
	for iL in range(2, nLines):
		line = data_string[iL].split()
		if line[0] == 'H':
			H.append([float(line[1]), float(line[2]), float(line[3])])

		elif line[0] == 'I':
			I.append([float(line[1]), float(line[2]), float(line[3])])

		elif line[0] == 'C':
			C.append([float(line[1]), float(line[2]), float(line[3])])
			
		elif line[0] == 'N':
			N.append([float(line[1]), float(line[2]), float(line[3])])

		elif line[0] == 'Pb':
			Pb.append([float(line[1]), float(line[2]), float(line[3])])
			
		elif line[0] == 'O':
			O.append([float(line[1]), float(line[2]), float(line[3])])
			
		else:
			print('ERROR! Element unknown.')

	Pb_pos = np.array(Pb)
	I_pos = np.array(I)
	C_pos = np.array(C)
	N_pos = np.array(N)
	H_pos = np.array(H)
	O_pos = np.array(O)

	" The amount of each species "
	Pb_num = int(len(Pb_pos[:, 0]))
	I_num = int(len(I_pos[:, 0]))
	C_num = int(len(C_pos[:, 0]))
	N_num = int(len(N_pos[:, 0]))
	O_num = int(len(O_pos[:, 0]))
	H_num = int(len(H_pos[:, 0]))


	print('==================================')
	print('\tPb \t=\t', Pb_num)
	print('\tI \t=\t', I_num)
	print('\tC \t=\t', C_num)
	print('\tN \t=\t', N_num)
	print('\tO \t=\t', O_num)
	print('\tH \t=\t', H_num)
	print('==================================')

	"--------------------------------"
	"-----Moving-values-into-box-----"
	"--------------------------------"
	"""
	" Lead "
	Pb_pos[Pb_pos[:,0] < 0, 0] = Pb_pos[Pb_pos[:,0] < 0, 0] + x_cage
	Pb_pos[Pb_pos[:,0] > x_cage, 0] = Pb_pos[Pb_pos[:,0] > x_cage, 0] - x_cage
	Pb_pos[Pb_pos[:,1] < 0, 1] = Pb_pos[Pb_pos[:,1] < 0, 1] + y_cage
	Pb_pos[Pb_pos[:,1] > y_cage, 1] = Pb_pos[Pb_pos[:,1] > y_cage, 1] - y_cage
	Pb_pos[Pb_pos[:,2] < 0, 2] = Pb_pos[Pb_pos[:,2] < 0, 2] + z_cage
	Pb_pos[Pb_pos[:,2] > z_cage, 2] = Pb_pos[Pb_pos[:,2] > z_cage, 2] - z_cage
	
	" Iodide "
	I_pos[I_pos[:,0] < 0, 0] = I_pos[I_pos[:,0] < 0, 0] + x_cage
	I_pos[I_pos[:,0] > x_cage, 0] = I_pos[I_pos[:,0] > x_cage, 0] - x_cage
	I_pos[I_pos[:,1] < 0, 1] = I_pos[I_pos[:,1] < 0, 1] + y_cage
	I_pos[I_pos[:,1] > y_cage, 1] = I_pos[I_pos[:,1] > y_cage, 1] - y_cage
	I_pos[I_pos[:,2] < 0, 2] = I_pos[I_pos[:,2] < 0, 2] + z_cage
	I_pos[I_pos[:,2] > z_cage, 2] = I_pos[I_pos[:,2] > z_cage, 2] - z_cage
	

	" Carbon "
	C_pos[C_pos[:,0] < 0, 0] = C_pos[C_pos[:,0] < 0, 0] + x_cage
	C_pos[C_pos[:,0] > x_cage, 0] = C_pos[C_pos[:,0] > x_cage, 0] - x_cage
	C_pos[C_pos[:,1] < 0, 1] = C_pos[C_pos[:,1] < 0, 1] + y_cage
	C_pos[C_pos[:,1] > y_cage, 1] = C_pos[C_pos[:,1] > y_cage, 1] - y_cage
	C_pos[C_pos[:,2] < 0, 2] = C_pos[C_pos[:,2] < 0, 2] + z_cage
	C_pos[C_pos[:,2] > z_cage, 2] = C_pos[C_pos[:,2] > z_cage, 2] - z_cage
	
	" Nitrogen "
	N_pos[N_pos[:,0] < 0, 0] = N_pos[N_pos[:,0] < 0, 0] + x_cage
	N_pos[N_pos[:,0] > x_cage, 0] = N_pos[N_pos[:,0] > x_cage, 0] - x_cage
	N_pos[N_pos[:,1] < 0, 1] = N_pos[N_pos[:,1] < 0, 1] + y_cage
	N_pos[N_pos[:,1] > y_cage, 1] = N_pos[N_pos[:,1] > y_cage, 1] - y_cage
	N_pos[N_pos[:,2] < 0, 2] = N_pos[N_pos[:,2] < 0, 2] + z_cage
	N_pos[N_pos[:,2] > z_cage, 2] = N_pos[N_pos[:,2] > z_cage, 2] - z_cage
	
	" Oxygen "
	O_pos[O_pos[:,0] < 0, 0] = O_pos[O_pos[:,0] < 0, 0] + x_cage
	O_pos[O_pos[:,0] > x_cage, 0] = O_pos[O_pos[:,0] > x_cage, 0] - x_cage
	O_pos[O_pos[:,1] < 0, 1] = O_pos[O_pos[:,1] < 0, 1] + y_cage
	O_pos[O_pos[:,1] > y_cage, 1] = O_pos[O_pos[:,1] > y_cage, 1] - y_cage
	O_pos[O_pos[:,2] < 0, 2] = O_pos[O_pos[:,2] < 0, 2] + z_cage
	O_pos[O_pos[:,2] > z_cage, 2] = O_pos[O_pos[:,2] > z_cage, 2] - z_cage
	
	
	" Hydrogen "
	H_pos[H_pos[:,0] < 0, 0] = H_pos[H_pos[:,0] < 0, 0] + x_cage
	H_pos[H_pos[:,0] > x_cage, 0] = H_pos[H_pos[:,0] > x_cage, 0] - x_cage
	H_pos[H_pos[:,1] < 0, 1] = H_pos[H_pos[:,1] < 0, 1] + y_cage
	H_pos[H_pos[:,1] > y_cage, 1] = H_pos[H_pos[:,1] > y_cage, 1] - y_cage
	H_pos[H_pos[:,2] < 0, 2] = H_pos[H_pos[:,2] < 0, 2] + z_cage
	H_pos[H_pos[:,2] > z_cage, 2] = H_pos[H_pos[:,2] > z_cage, 2] - z_cage

	"""
	"-----------------------------------"
	"-----Removing-excess-molecules-----"
	"-----------------------------------"
	" We look at all molecule types except hydrogen and check 		"
	" wether there are molecules that are too close to an identical "
	" molecule. This can happen when the xyz file is not built to 	"
	" be replicated. 												"

	" We do this first for lead and then it is repeated for the 	"
	" other molecule types in same way. 							"

	" Lead "
	" A vector that tracks which atoms should be removed "
	Pb_remove = np.zeros(Pb_num)
	for iremP in range(Pb_num - 1):
		Pbx1, Pby1, Pbz1 = Pb_pos[iremP, :]
		
		Pbrest = Pb_pos[iremP + 1:, :]
		Pbx2 = Pbrest[:, 0]
		Pby2 = Pbrest[:, 1]
		Pbz2 = Pbrest[:, 2]

		PPx = np.abs(Pbx2 - Pbx1)
		PPx[PPx > x_cage/2] = x_cage - PPx[PPx > x_cage/2]

		PPy = np.abs(Pby2 - Pby1)
		PPy[PPy > y_cage/2] = y_cage - PPy[PPy > y_cage/2]

		PPz = np.abs(Pbz2 - Pbz1)
		PPz[PPz > z_cage/2] = z_cage - PPz[PPz > z_cage/2]


		PPdist = np.sqrt(PPx**2 + PPy**2 + PPz**2)

		" We locate where the distance is too short "
		PPrem = np.argwhere(PPdist < 2) + iremP + 1
		" And change the Pb_remove value to be 1 there "
		Pb_remove[PPrem] = 1

	" We find the values that should be deleted and delet those rows "
	Pbdel = np.argwhere(Pb_remove > 0.5)
	Pb_pos = np.delete(Pb_pos, Pbdel, axis=0)
	Pb_num = int(len(Pb_pos[:, 0]))

	" Iodine "
	I_remove = np.zeros(I_num)
	for iremI in range(I_num - 1):
		Ix1, Iy1, Iz1 = I_pos[iremI, :]
		
		Irest = I_pos[iremI + 1:, :]
		Ix2 = Irest[:, 0]
		Iy2 = Irest[:, 1]
		Iz2 = Irest[:, 2]

		IIx = np.abs(Ix2 - Ix1)
		IIx[IIx > x_cage/2] = x_cage - IIx[IIx > x_cage/2]

		IIy = np.abs(Iy2 - Iy1)
		IIy[IIy > y_cage/2] = y_cage - IIy[IIy > y_cage/2]

		IIz = np.abs(Iz2 - Iz1)
		IIz[IIz > z_cage/2] = z_cage - IIz[IIz > z_cage/2]

		IIdist = np.sqrt(IIx**2 + IIy**2 + IIz*2)

		IIrem = np.argwhere(IIdist < 2) + iremI + 1
		
		I_remove[IIrem] = 1

	Idel = np.argwhere(I_remove > 0.5)
	I_pos = np.delete(I_pos, Idel, axis=0)
	I_num = int(len(I_pos[:, 0]))

	print('==================================')
	print('\tPb \t=\t', Pb_num)
	print('\tI \t=\t', I_num)
	print('\tC \t=\t', C_num)
	print('\tN \t=\t', N_num)
	print('\tO \t=\t', O_num)
	print('\tH \t=\t', H_num)
	print('==================================')

	"----------------------------"
	"-----Writing-the-output-----"
	"----------------------------"

	" We calculate the total number of atoms we will write "
	total_atoms = Pb_num + I_num + N_num + C_num + H_num +O_num
	#print np.count_nonzero(H_lab==7)
	#print H_lab
	" Finally we write our output it is structured as follows 										"
	" Line 1: The number of atoms 																	"
	" Line 2: Lattice dimensions and specifications of the lines that follow 						"
	" The following lines are specifications for all the atoms in the 								"	
	" system they have the following structure: 													"
	" Element - x coordinate - y coordinate - z coordinate - unit cell number(1) - element label 	"	
	" The element label differs between N on MA and N on BA and so forth 							"

	f = open(argv[5], 'w')
	f.write('%d\n' %total_atoms)
	f.write('Lattice="%.4f 0.0 0.0 0.0 %.4f 0.0 0.0 0.0 %.4f" Properties=species:S:1:pos:R:3:molid:I:1:type:S:1 pbc="T T T"\n' %(x_cage, y_cage, z_cage))
	for wP in range(Pb_num):
		f.write('Pb %.6f %.6f %.6f 1 P1\n' %(Pb_pos[wP, 0], Pb_pos[wP, 1], Pb_pos[wP, 2]))
	for wI in range(I_num):
		f.write('I %.6f %.6f %.6f 1 I1\n' %(I_pos[wI, 0], I_pos[wI, 1], I_pos[wI, 2]))
	for wN in range(N_num):
		f.write('N %.6f %.6f %.6f 1 N1\n' %(N_pos[wN, 0], N_pos[wN, 1], N_pos[wN, 2]))
	
	for wO in range(O_num):
		f.write('O %.6f %.6f %.6f 1 O1\n' %(O_pos[wO, 0], O_pos[wO, 1], O_pos[wO, 2]))

	for wC in range(C_num):
		f.write('C %.6f %.6f %.6f 1 Cx\n' %(C_pos[wC, 0], C_pos[wC, 1], C_pos[wC, 2]))
		
	for wH in range(H_num):
		f.write('H %.6f %.6f %.6f 1 Hx\n' %(H_pos[wH, 0], H_pos[wH, 1], H_pos[wH, 2]))
		
	f.close()


if __name__ == "__main__":
    main(sys.argv) 