#Perform checks to ensure data consistency

from ase.io import read #Use to collect, visualise atomic structures
import numpy as np #Mathematical calculations
from ase.data import chemical_symbols  # ASE chemical symbols
from collections import Counter

# Check full compatibility between two atomic cells
def check_cell_compatibility(atoms, reference_atoms, tolerance=1e-5):

    # Collect the cell matrices
    cell = atoms.get_cell().array
    reference_cell = reference_atoms.get_cell().array

    # Check cell dimensions and orientation
    if not np.allclose(cell, reference_cell, atol=tolerance, rtol=0):
        raise ValueError("Cell mismatch: the structures have incompatible cell dimensions or orientation.")

#Check length consistency of Bader files
def check_bader_consistency(item, ini_q, fin_q):

    # Bader vs Bader
    if len(ini_q) != len(fin_q): #i.e., Check that Bader charge lengths are consistent
        raise ValueError( #If not, raise an error
            f"Atom count mismatch in {item['transition']}: ini={len(ini_q)}, fin={len(fin_q)}")

    #Atoms in each structure
    ini_atoms = read(item["ini_structure"])
    fin_atoms = read(item["fin_structure"])

    #Extract length of atoms
    ini_len = len(ini_atoms)
    fin_len = len(fin_atoms)

    # POSCAR vs Bader atom length
    if (fin_len != len(fin_q)) or (ini_len != len(ini_q)): #i.e., Check that Bader charge and POSCAR/CONTCAR lengths are consistent
        raise ValueError( #If not, raise an error
            f"Mismatch between POSCAR/CONTCAR and Bader atoms in {item['transition']}")

    # Check cell compatibility between initial and final structures. If the images are not derived from NEB calculations, then they may still be inconsistent.
    check_cell_compatibility(ini_atoms, fin_atoms)
    
    #Check element ordering consistency
    for i, (a_ini, a_fin) in enumerate(zip(ini_atoms, fin_atoms)):
        if a_ini.symbol != a_fin.symbol:
            raise ValueError(f"Element mismatch at index {i+1} in {item['transition']}: {a_ini.symbol} (ini) != {a_fin.symbol} (fin)")

#Check that Bader coordinates and POSCAR/CONTCAR coordinates are equivalent
def check_bader_alignment(atoms, acf_coords, tol=1e-3):

    #N.B. This tolerance is to not be confused with that of Bader charge difference.
    #This difference is to check that atoms in POSCAR/CONTCAR and ACF.dat are the same.
    #Here, difference should not be too large, and hence it is defaulted.

    pos = atoms.get_positions()  # Cartesian from ASE

    #Calculate difference based on absolute distance
    diff = np.linalg.norm(pos - acf_coords, axis=1)

    #Should be within threshold
    if np.max(diff) > tol:
        raise ValueError( #Else raise an error
            f"ACF.dat coordinates do not match POSCAR/CONTCAR (max diff = {np.max(diff):.4f} Å)")

#Check that the reference structure contains the requested elements and that the number of reference atoms matches the supplied input/final structure.
def check_reference_structure(atoms, reference_atoms, reference_symbols):

    # Check symbols
    if reference_symbols is None or len(reference_symbols) == 0:
        raise ValueError("Reference symbols must contain at least one element symbol.")

    for symbol in reference_symbols:
        if symbol not in chemical_symbols[1:]:
            raise ValueError(f"Invalid chemical element symbol: '{symbol}'.")

    # Check that the full cells between the images and the reference structure are equivalent
    check_cell_compatibility(atoms, reference_atoms)

    #Check that every requested reference element exists in the reference structure.
    reference_counts = Counter(reference_atoms.get_chemical_symbols())
    #Check that the supplied initial and final structures contain the same number of reference atoms as the supplied reference structure.
    atoms_counts = Counter(atoms.get_chemical_symbols())

    for symbol in reference_symbols:

        if reference_counts[symbol] == 0: #Should contain at least one symbol
            raise ValueError(f"Reference structure does not contain the requested element '{symbol}'.")

        if atoms_counts[symbol] != reference_counts[symbol]:
            raise ValueError(f"Reference mismatch for '{symbol}': reference contains {reference_counts[symbol]} atoms, but the supplied structure contains {atoms_counts[symbol]}.")