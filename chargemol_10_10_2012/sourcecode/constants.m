%%% Constants
% The default values should give good results for all chemical systems, so
% don't change them unless you have a compelling reason
bohrperangstrom=1.889725989;
nshells=100 % number of radial integration shells
cutoff_radius = 500 % this is in picometers
zero_tolerance=1.0e-08
integration_tolerance = 0.10
integration_tolerance_percent = 0.03
minpixelvolume = 0.0081
pixel_integration_tolerance = 0.03
charge_convergence_tolerance = 1.0e-5
spin_convergence_tolerance = 5.0e-5
distance_tolerance = 1.0e-6
reference_weighting = 3/14
rmin_cloud_penetration=200 % this is in picometers
density_decaying_exponent = 1.75
density_set_prefix = 'c2' 
spin_ref_fraction = 1/2
noble_gas_core
wA_renormalization_max=10.0

% some default input files for different types of jobs
%
% default name for xsf file
xsf_inputfile='valence_density.xsf';
%
% default direction of SAXIS in VASP
SAXIS = [1.0e-12, 0.0, 1.0]

% Parameters used only for reading gaussian basis set wfx files
preferred_grid_spacing=0.14
gaussian_overlap_tolerance=1.0e-12
analytic_alpha_cutoff = 5.0 %primitive pairs with alpha sum greater than this are treated analytically
periodic_cutoff_length = 28.0 % This is in bohr

% version information
version = ' ';



