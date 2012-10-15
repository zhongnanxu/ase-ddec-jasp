strcat('Entering ', version)
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf ( 1, '%s\n', s );
% Read the input options
% make sure the periodicity selection is valid
if periodicA ~= 1 && periodicA ~= 0
    'An invalid option has been selected for the periodicity.'
    'Program will terminate.'
    break
end
if periodicB ~= 1 && periodicB ~= 0
    'An invalid option has been selected for the periodicity.'
    'Program will terminate.'
    break
end
if periodicC ~= 1 && periodicC ~= 0
    'An invalid option has been selected for the periodicity.'
    'Program will terminate.'
    break
end
% Read the electron density, atomic coordinates, etc. and perform the
% analysis
scalefactor = 100.0*nshells/(cutoff_radius*bohrperangstrom);
run_valence_core_densities
if flag == 1
    break
end
core_iterator  
if flag == 1
    break
end
new_valence_iterator
if flag == 1
    break
end
local_multipole_moment_analysis
total_multipole_moment_analysis
cloud_penetration
generate_xyzfile
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf ( 1, '%s\n', s );
if spin_available ~= 0 
   generate_spin_lookup_tables
end   
if spin_available == 1
    collinear_spin_moments_iterator
elseif spin_available == 2
    noncollinear_spin_moments_iterator    
end      
quote
