'Normal termination of Chargemol program.'
'Use and distribution of this program is subject to certain licensing restrictions.'
'Please see ddec.sourceforge.net for details.'
'Please cite the following journal references when using this program to compute DDEC, Hirshfeld, or ISA charges:'
'Thomas A. Manz and David S. Sholl, "Improved Atoms-in-Molecule Charge Partitioning Functional for Simultaneously Reproducing the Electrostatic Potential and Chemical States in Periodic and Non-Periodic Materials", J. Chem. Theory Comput., Vol. 8 (2012) 2844-2867.'
'Thomas A. Manz and David S. Sholl, "Chemically Meaningful Atomic Charges that Reproduce the Electrostatic Potential in Periodic and Nonperiodic Materials", J. Chem. Theory Comput., Vol. 6 (2010) 2455-2468.'
'Please cite the following journal reference when using this program to compute atomic spin moments:'
'Thomas A. Manz and David S. Sholl, "Methods for Computing Accurate Atomic Spin Moments for Collinear and Noncollinear Magnetism in Periodic and Nonperiodic Materials", J. Chem. Theory Comput., Vol. 7 (2011) 4146-4164.'
'Exiting Chargemol'
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf ( 1, '%s\n', s );
fprintf ( 1, '%s\n', version);