% This file contains the parameters that effect the RRMSE computation
% The inner and outer multipliers select the spatial volume
% Points closer to an atom than inner_multiplier x vdW radius are discarded
% Points farther from all atoms than outer_multiplier x vdW radius are discarded

% the following values are recommended when for nanoporous solids like zeolites and metal organic frameworks
inner_multiplier=1.3
outer_multiplier=20.0 % for nanoporous solids like metal organic frameworks and zeolites the outer multiplier is set high so that integration proceeds to center of the pores

% the following values are recommended for molecular systems, slabs, nonperiodic systems, and systems only periodic in 1 or 2 dimensions.
% inner_multiplier=1.4
% outer_multiplier=2.0


% The following parameter allows computation to be speeded by thinning the grid
% For example, skip = 2 uses half the grid points along the first, second, and third lattice vectors to select 1/8th the number of grid points
% Skip = 1 does not thin the grid at all
skip=1

%
% Select the vasp version being used
%%% Select 5 for the most recent VASP versions, select 4 for earlier versions.
%%% (If one of these doesn't work, try the other one.)  
vaspversion=4 % choose this option if there are six header lines before the number of atoms per type in the AECCAR2 file
%vaspversion=5 % choose this option if there are seven header lines before the number of atoms per type in the AECCAR2 file 
