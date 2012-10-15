if input_type == 1
    format_gaussian_cube_densities;
elseif input_type == 2
    format_vasp_densities;
elseif input_type == 3
    format_xsf_densities;
elseif input_type == 4
    format_nongaussian_cube_densities
elseif input_type == 5
    read_wfx_file
    if flag == 1
       break
    end
    generate_density_grids_from_gaussian_basis_set_coefficients
elseif input_type == 6
    format_total_cube_density
else
    'Fatal error. You selected a wrong input file type. Program will terminate.'
    break
end
if flag == 1
   break
end
ncore = sum(core_electrons)
nvalence=sum(atomic_number)- netcharge - ncore
integration_tolerance = max(integration_tolerance,integration_tolerance_percent*nvalence/100.0);
sum_valence_occupancy_correction=0;
for j=1:natoms
    sum_valence_occupancy_correction = sum_valence_occupancy_correction + occupancy_correction(j,1);
end
checkme = abs(nvalence - sum(sum(sum(valence_density(:,:,:))))*pixelvolume - sum_valence_occupancy_correction)
if (checkme > integration_tolerance) || (abs(sum_negative_density) > integration_tolerance)
    'The electrons are not properly accounted for.'
    'Either the grid in your electron density input file is too coarse, you have specified the incorrect net charge in the chargemol_job.m file,'
    'or the number of core electrons for each atom has not been setup correctly.'
    'Program will terminate'
    break
else
    'The grid spacing is adequate and all electrons are properly accounted for.'
    'Calculation will proceed.'
end
if spin_available == 2
   check_noncollinear_XC_functional    
end