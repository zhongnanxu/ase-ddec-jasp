% total multipole moment analysis
% the multipole moments are computed in two parts: an analytic part which contains the atomic nuclei and the occupancy corrections
% and a numeric part which contains the multipole moments of the valence density grid
if ((periodicA == 0) && (periodicB == 0)) && (periodicC == 0)
    compute_center_of_mass
    'Since the system is nonperiodic, total multipole moment analysis will be performed.'
    'Total multipole moments are computed using the center of mass as the origin.'
    'Dipole and quadrupole moments using the net atomic charges.'
    'This corresponds to truncatating the distributed multipole expansion at monopole order.'
    partial_charge_dipole_x = 0.0;
    partial_charge_dipole_y = 0.0;
    partial_charge_dipole_z = 0.0;
    partial_charge_quadrupole_xy = 0.0;
    partial_charge_quadrupole_xz = 0.0;
    partial_charge_quadrupole_yz = 0.0;
    partial_charge_quadrupole_x2minusy2 = 0.0;
    partial_charge_quadrupole_3z2minusr2 = 0.0;
    for j = 1:natoms
        x = coords(j,1) - center_of_mass(1);
        y = coords(j,2) - center_of_mass(2);
        z = coords(j,3) - center_of_mass(3);
        charge = final_result(j,6);
        partial_charge_dipole_x = partial_charge_dipole_x + charge*x;
        partial_charge_dipole_y = partial_charge_dipole_y + charge*y;
        partial_charge_dipole_z = partial_charge_dipole_z + charge*z;
        partial_charge_quadrupole_xy = partial_charge_quadrupole_xy + charge*x*y;
        partial_charge_quadrupole_xz = partial_charge_quadrupole_xz + charge*x*z;
        partial_charge_quadrupole_yz = partial_charge_quadrupole_yz + charge*y*z;
        partial_charge_quadrupole_x2minusy2 = partial_charge_quadrupole_x2minusy2 + charge*(x^2 - y^2);
        partial_charge_quadrupole_3z2minusr2 = partial_charge_quadrupole_3z2minusr2 + charge*(3*z*z - (x^2 + y^2 + z^2));
    end
    partial_charge_dipole_x
    partial_charge_dipole_y
    partial_charge_dipole_z
    partial_charge_dipole_magnitude = sqrt(partial_charge_dipole_x^2 + partial_charge_dipole_y^2 + partial_charge_dipole_z^2)
    partial_charge_quadrupole_xy
    partial_charge_quadrupole_xz
    partial_charge_quadrupole_yz
    partial_charge_quadrupole_x2minusy2
    partial_charge_quadrupole_3z2minusr2
    %
    'Dipole and quadrupole moments using the net atomic charges and the atomic dipoles.'
    'This corresponds to truncatating the distributed multipole expansion at dipole order.'
    dipole_x = partial_charge_dipole_x;
    dipole_y = partial_charge_dipole_y;
    dipole_z = partial_charge_dipole_z;   
    quadrupole_xy = partial_charge_quadrupole_xy;
    quadrupole_xz = partial_charge_quadrupole_xz;
    quadrupole_yz = partial_charge_quadrupole_yz;
    quadrupole_x2minusy2 = partial_charge_quadrupole_x2minusy2;
    quadrupole_3z2minusr2 = partial_charge_quadrupole_3z2minusr2;
    for j = 1:natoms
        x = coords(j,1) - center_of_mass(1);
        y = coords(j,2) - center_of_mass(2);
        z = coords(j,3) - center_of_mass(3);
        dipole_x = dipole_x + final_result(j,7);
        dipole_y = dipole_y + final_result(j,8);
        dipole_z = dipole_z + final_result(j,9);       
        quadrupole_xy = quadrupole_xy + x*final_result(j,8) + y*final_result(j,7);
        quadrupole_xz = quadrupole_xz + x*final_result(j,9) + z*final_result(j,7);
        quadrupole_yz = quadrupole_yz + y*final_result(j,9) + z*final_result(j,8);
        quadrupole_x2minusy2 = quadrupole_x2minusy2 + 2*(x*final_result(j,7) - y*final_result(j,8));
        quadrupole_3z2minusr2 = quadrupole_3z2minusr2 + 2*(2*z*final_result(j,9) - x*final_result(j,7) - y*final_result(j,8));      
    end
    dipole_x
    dipole_y
    dipole_z 
    dipole_magnitude = sqrt(dipole_x^2 + dipole_y^2 + dipole_z^2)
    quadrupole_xy
    quadrupole_xz
    quadrupole_yz
    quadrupole_x2minusy2
    quadrupole_3z2minusr2 
    %
    'Dipole and quadrupole moments using the net atomic charges, atomic dipoles, and atomic quadrupoles.'
    'This corresponds to truncatating the distributed multipole expansion at quadrupole order.'
    dipole_x
    dipole_y
    dipole_z
    dipole_magnitude = sqrt(dipole_x^2 + dipole_y^2 + dipole_z^2)
    for j = 1:natoms      
        quadrupole_xy = quadrupole_xy + final_result(j,11);
        quadrupole_xz = quadrupole_xz + final_result(j,12);
        quadrupole_yz = quadrupole_yz + final_result(j,13);
        quadrupole_x2minusy2 = quadrupole_x2minusy2 + final_result(j,14);
        quadrupole_3z2minusr2 = quadrupole_3z2minusr2 + final_result(j,15);     
    end  
    quadrupole_xy
    quadrupole_xz
    quadrupole_yz
    quadrupole_x2minusy2
    quadrupole_3z2minusr2      
end
    
    
    
    
    
    
    
  