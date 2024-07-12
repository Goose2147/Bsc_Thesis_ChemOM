# Bsc_Thesis_ChemOM
This repository contains all python scripts and examples I used during my Bachelor Thesis at the ChemOM department.

# LAMMPS Data File

The LAMMPS data file necessary for running simulations can be processed from a CIF file. The Crystallographic Information File (CIF) is a standard text file format for representing crystallographic information. All Python scripts and other files mentioned in this guide are available on my repository on GitHub: [https://github.com/Goose2147/Bsc_Thesis_ChemOM](https://github.com/Goose2147/Bsc_Thesis_ChemOM).

## Steps for Creating a LAMMPS Data File from the CIF File for a 2D A<sub>2</sub>M<sup>2+</sup>X<sub>4</sub> Perovskite
Where A = organic spacer, M = metal (Pb in this study), X = Halide (I in this study).

1. **Download CIF file from a database:**
   - In this study, CIF files were obtained from: [http://www.pdb.nmse-lab.ru/](http://www.pdb.nmse-lab.ru/)​[1]​.

2. **Open this CIF file in Mercury ​[2]​:**
   - Delete the molecules from inorganic plane and delete the organic spacer molecules such that only one spacer molecule remains.
   - Save this as a MOL2 file (e.g., in the case of NOE: `1_napthyl_ethyl.mol2`) and as a XYZ file (e.g., `1_napthyl_ethyl.xyz`).

3. **Check the MOL2 file:**
   - If under BOND there’s ‘un’ values, these need to be fixed. Usually this happens with aromatic bonds; then ‘ar’ can be used as a substitute.

4. **Start the process of making the PAR file for ASE (`forase.par`):**
   - Go to [http://www.swissparam.ch/](http://www.swissparam.ch/) and upload your MOL2 file.
   - Select approach MATCH and launch it ​[3], [4], [5]​.
   - Download the zip file. You only need the PAR file from this folder.

5. **Use `par_fix_format.py` to create a new PAR file with the correct format.**

6. **Use `modify_mol2.py` to create a new MOL2 file with the handy format.**
   - Open this modified MOL2 file in Mercury for visualization. 

7. **Use the new PAR and MOL2 files from step 6 and 7 with `combine_mol2_par.py` to create a `forase.par` file (without the charges):**
   - This python script will use the function assign_MATCH_label. You can copy the printed library to the python script to skip the assignment process if it's run again for the same structure.
   - Symmetrical hydrogens are not repeated, but rather contained in 1 hydrogen type.

8. **Obtain charges from Gaussian09​[6]​:**
   - Use OpenBabel to convert the XYZ file to a gjf file (Gaussian input) [​7]​.
   - Gaussian09 is available on the high-performance computer (HPC) cluster.
   - Run the Gaussian input file. It will output ’charges from ESP fit’ multiple times as it converges.
   - Take the last results from ’charges from ESP fit’ and fill them in the `forase.par` file from the previous step. Some hydrogens are symmetrical, so take those charges and use the average.

9. **The `forase.par` file is done:**
   - Make sure the sum of the charges of all atoms of the organic spacer equals the real charge (NOE: charge = +1).

## Process of Making the XYZ file for ASE (`forase.xyz`)

1. **Open the CIF file from step 2 in VESTA ​[8]​:**
   - Go to ‘Objects’, then click on ‘Boundary ...’.
   - Edit the axes so that the spacer molecules are sandwiched by the inorganic plane from the top and bottom (For NOE: axes x and y remained the same and axis z was edited to ‘-0.1’ to ‘1.0’).

2. **Clean up the CIF file:**
   - Delete all free-floating singular atoms.
   - Determine which atoms in the inorganic planes would be duplicates, if the same unit cell were added in the x, y, and z direction. If an atom in the inorganic plane is a duplicate, it can be deleted.
   - Determine which spacer molecules would be duplicates, if the same unit cell were added in the x, y, and z direction. Unlike the inorganic plane, the atoms in the organic spacer molecules cannot be deleted individually, as this would break the molecule (In the case of NOE: 8 Pb-atoms, 32 I-atoms, 16 spacer molecules).
   - Save this as a XYZ file (e.g., `NOE_final.xyz`).

3. **Use `xyz_fix_format.py` with this XYZ file to get the formatted XYZ file:**
   - This XYZ file has unlabeled Carbon and Hydrogen atom types (e.g., C3, H5).
   - Open the formatted XYZ file in VESTA and check for each atom which atom type from the `forase.par` should be assigned to each atom.
   - Check that the right amount of each atom type is represented in the XYZ file (e.g., for 16 spacer molecules there are 16 C1-atoms).
   - Save this file as `forase.xyz`.

4. **With both `forase.par` and `forase.xyz`, ASE can be run:**
   - Check `opls.py` to change ‘symlen’ when running function `read_block(name, symlen, nvalues)`.
   - Check `asefile.py` to change the element dictionary.
   - Run `asefile.py`.

5. **ASE will output:**
   - `lmp_opls`, `lmp_atoms`, `lmp_in`. The first two are needed and `lmp_in` is an example input file for LAMMPS which can be deleted.

6. **Use `merging_lmp.py` with `lmp_opls` and `lmp_atoms` to output `lmp_atoms_adj` and `coeffs`:**
   - Only `lmp_atoms_adj` is important.

7. **Fix the coordinates in `lmp_atoms_adj`:**
   - Open `lmp_atoms_adj` and go to header ‘Atoms’. The xyz-coordinates have been wrapped by ASE using PBC to remain within the cell boundaries. However, this will break molecules within the cell.
   - To fix this, the xyz-coordinates will be replaced with the ones from `forase.xyz`.
   - Copy all the lines underneath ‘Atoms’ in the `lmp_atoms_adj` file and paste them in ‘original.txt’.
   - Copy only the lines with coordinates from the `forase.xyz` file into ‘replacement.txt’ and run `replace.py`.
   - Use ‘updated.txt’ to fix the `lmp_atoms_adj`.
   - You can also run `tester.py` to check that the unit cell is charge neutral.

8. **Convert parameters under ‘Dihedral Coeffs’ in `lmp_atoms_adj`:**
   - The 3rd and 4th column need to convert from float to integer objects.

9. **Rename `lmp_atoms_adj` if desired:**
   - (NOE: `lmp_atoms_NOE***`).

10. **The LAMMPS data file is done!**

*Disclaimer: In the making of `lmp_atoms_NOE`, step 1 (process of making xyz file for ASE) was not performed. The system still works and has the right number of atoms. However, the spacers are not sandwiched on the bottom (height = z-axis). Thus, it is not the exact desired system for molecular dynamic simulations.*

## LAMMPS Pair Coefficient File

A pair coefficient file needs to be included alongside the data file to run a LAMMPS input file successfully. The creation was done with the help of the work from M. B. Fridriksson’s PhD Thesis [​9]​. On page 82-83 there are tables with Buckingham and Lennard-Jones potential coefficients for each atom pair. The Lennard-Jones coefficients are spacer-dependent. For NOE, the coefficients of POB were used as it is the most similar. For constructing your own pair coefficient file, use `paircoeff_NOE` on my repository as a reference.

## Other Files Used for MD Simulations

Examples of LAMMPS input files can also be found on the repository. These input files often create extra output files chosen by the user. These output files are further processed to visualize and quantify interesting properties. These output files and the python files used are also available.

To understand LAMMPS input files, first check the starting point: `input_NOE_nve.in`. From the saved ‘equil’ (Equilibrated) files, many other input files can continue processing, thus skipping the NVE and NTP steps.

## References
​​[1] E. I. Marchenko et al., “Database of Two-Dimensional Hybrid Perovskite Materials: Open-Access Collection of Crystal Structures, Band Gaps, and Atomic Partial Charges Predicted by Machine Learning,” Chemistry of Materials, vol. 32, no. 17, pp. 7383–7388, Sep. 2020, doi: 10.1021/acs.chemmater.0c02290. 

​[2] C. F. Macrae et al., “Mercury 4.0 : from visualization to analysis, design and prediction,” J Appl Crystallogr, vol. 53, no. 1, pp. 226–235, Feb. 2020, doi: 10.1107/S1600576719014092. 

​[3] J. D. Yesselman, D. J. Price, J. L. Knight, and C. L. Brooks, “MATCH: An atom‐typing toolset for molecular mechanics force fields,” J Comput Chem, vol. 33, no. 2, pp. 189–202, Jan. 2012, doi: 10.1002/jcc.21963. 

​[4] V. Zoete, M. A. Cuendet, A. Grosdidier, and O. Michielin, “SwissParam: A fast force field generation tool for small organic molecules,” J Comput Chem, vol. 32, no. 11, pp. 2359–2368, Aug. 2011, doi: https://doi.org/10.1002/jcc.21816. 

​[5] M. Bugnon et al., “SwissParam 2023: A Modern Web-Based Tool for Efficient Small Molecule Parametrization,” J Chem Inf Model, vol. 63, no. 21, pp. 6469–6475, Nov. 2023, doi: 10.1021/acs.jcim.3c01053. 

​[6] Inc. Gaussian, “Gaussian 09, Revision B.01.” Gaussian, Inc., Wallingford CT, 2010. Accessed: Jul. 04, 2024. [Online]. Available: https://gaussian.com 

​[7] OpenBabel, “OpenBabel (Version 3.1.0).” 2022. Accessed: Jul. 04, 2024. [Online]. Available: http://openbabel.org 

​[8] Kyoto University, “VESTA (Version 3.5.8).” Japan, 2022. Accessed: Jul. 04, 2024. [Online]. Available: https://jp-minerals.org/vesta/en/ 

​[9] M. B. Fridriksson, “Structural and Excited State Dynamics in Hybrid Halide Perovskites,” Dissertation (TU Delft), Delft University of Technology, Delft, 2020. doi: 10.4233/uuid:4d48181b-5429-4bef-976d-95b6c9535825. 

