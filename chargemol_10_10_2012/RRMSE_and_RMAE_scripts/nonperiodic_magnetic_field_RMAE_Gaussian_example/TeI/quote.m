'Normal termination of Chargemol program.'
'Use and distribution of this program is subject to certain licensing restrictions.'
'Please see ddec.sourceforge.net for details.'
'Please cite the following journal reference when using this code to quantify how accurately a set of ASMs reproduces the magnetic field:'
'Thomas A. Manz and David S. Sholl, "Methods for Computing Accurate Atomic Spin Moments for Collinear and Noncollinear Magnetism in Periodic and Nonperiodic Materials", J. Chem. Theory Comput., Vol. 7 (2011) 4146-4164.'
'Exiting Chargemol'
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf ( 1, '%s\n', s );
