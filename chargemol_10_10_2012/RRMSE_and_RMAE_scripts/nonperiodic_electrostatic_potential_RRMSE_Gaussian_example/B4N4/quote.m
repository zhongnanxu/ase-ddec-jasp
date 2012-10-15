'Normal termination of Chargemol program.'
'Use and distribution of this program is subject to certain licensing restrictions.'
'Please see ddec.sourceforge.net for details.'
'Please cite the following journal references when using this code to compute RMSE and RRMSE values:'
'Thomas A. Manz and David S. Sholl, "Improved Atoms-in-Molecule Charge Partitioning Functional for Simultaneously Reproducing the Electrostatic Potential and Chemical States in Periodic and Non-Periodic Materials", J. Chem. Theory Comput., Vol. 8 (2012) 2844-2867.'
'Taku Watanabe, Thomas A. Manz, and David S. Sholl, "Accurate Treatment of Electrostatics during Molecular Adsorption in Nanoporous Crystals Without Assigning Point Charges to Framework Atoms", J. Phys. Chem. C, Vol. 115, 2011, pp. 4824 - 4836.'
'Exiting Chargemol'
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf ( 1, '%s\n', s );
