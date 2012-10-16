clear % don't delete this line

ddec = getenv('DDEC_PATH')

%%% Job options
periodicA=1 % whether the system is (1) periodic or (0) not periodic along the first direction
periodicB=1 % whether the system is (1) periodic or (0) not periodic along the second direction
periodicC=1 % whether the system is (1) periodic or (0) not periodic along the third direction
netcharge=0 % the net charge of the periodic unit cell or molecular system

input_type=2 % whether the input file is: (1) gaussian cube file format, (2) VASP format, (3) XSF format, (4) CP2K cube file format
%%% For VASP input select version 4 or version 5
%%% Select 5 for the most recent VASP versions, select 4 for earlier versions.
%%% (If one of these doesn't work, try the other one.)  
%vaspversion=4 % choose this option if there are six header lines before the number of atoms per type in the AECCAR2 file
vaspversion=5 % choose this option if there are seven header lines before the number of atoms per type in the AECCAR2 file 

% SET UP THE DIRECTORY PATHS
% It is better to use absolute paths rather than relative paths,
% because then you can run the program from any directory on your system.
%
% EXAMPLES
%
atomic_densities_directory = strcat(ddec, '/chargemol_10_10_2012/atomic_densities/')
sourcecode_directory = strcat(ddec, '/chargemol_10_10_2012/sourcecode/')
% sourcecode_directory = 'C:/mypath/chargemol_06_12_2011/sourcecode/' % an example for a fixed location on windows system
% sourcecode_directory = '~/bin/chargemol_06_12_2011/sourcecode/' % for local bin directory on linux system
%
% You get the idea ...

%%% Run the program
copyfile(strcat(sourcecode_directory,'*.*'),'.')
constants
%
%
% BEGIN CHANGE CONSTANT VALUES SECTION
% 
xsf_inputfile='valence_density.xsf' % Example of how to change input file name if you want to; otherwise program uses valence_density.xsf.
num_core(30) = 18; % Example of changing number of core electrons for an element (e.g. atomic number = 30); otherwise program uses noble gas core or reads them from input files.
%
% END CHANGE CONSTANT VALUES SECTION
%
%
chargemol
% Uncomment the following line if you want to print out the weight fractors wA(rA) for each atom.
% print_atomic_densities
file_cleanup
delete('./file_cleanup.m');

quit % whether to quit Matlab after the program finishes
