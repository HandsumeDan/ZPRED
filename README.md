# Zeta Potential Prediction in a General Aqueous Electrolyte Solution
This software (ZPRED) has been experimentally validated and accompanies the journal
article that can be found here:
https://doi.org/10.1021/acs.langmuir.0c02031

ZPRED is a driver program in that it automates the use of various
software components to compute the electrokinetic (or zeta) potential
of a macromolecule.  Thus, successful implementation of ZPRED requires
that all software components are installed correctly.

To use ZPRED, it is necessary to specify solution conditions (solvent, temperature,
ion concentrations, etc.) along with providing a structure or ensemble of molecular
structures. Together this information is used to model the collection of solvation layers 
around the molecule (i.e. the electric double layer) and capture the average
electric potential at the molecule's predicted solvation edge.

Either .pdb and .pqr file formats can be input for structures.

