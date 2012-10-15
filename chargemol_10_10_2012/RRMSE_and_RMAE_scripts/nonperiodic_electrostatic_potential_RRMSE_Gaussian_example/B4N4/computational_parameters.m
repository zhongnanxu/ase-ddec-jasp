% This file contains the parameters that effect the RRMSE computation
% The inner and outer multipliers select the spatial volume
% Points closer to an atom than inner_multiplier x vdW radius are discarded
% Points farther from all atoms than outer_multiplier x vdW radius are discarded
% The recommended values of 1.4 (inner multiplier) and 2.0 (outer multiplier) are 
% taken from this reference:
% Singh, U.C.; Kollman, P.A., J. Comput. Chem. 1984, Vol. 5, 129 - 145.
inner_multiplier=1.4
outer_multiplier=2.0