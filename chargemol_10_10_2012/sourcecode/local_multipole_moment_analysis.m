% local multipole moment analysis
final_result=zeros(natoms,15);
count=0;
for i=1:natoms
    final_result(i,1) = i;
    final_result(i,2) = atomic_number(i);
    final_result(i,3)=coords(i,1);
    final_result(i,4)=coords(i,2);
    final_result(i,5)=coords(i,3);  
end
final_result(:,6)= atomic_number(:) - core_electrons(:) - valence_population(:);
% compute the local dipole and quadrupole moments in atomic units
temp=0.0;
sum_dipole_x = zeros(natoms,1);
sum_dipole_y = zeros(natoms,1);
sum_dipole_z = zeros(natoms,1);
sum_quadrupole_xy = zeros(natoms,1);
sum_quadrupole_xz = zeros(natoms,1);
sum_quadrupole_yz = zeros(natoms,1);
sum_quadrupole_x2minusy2 = zeros(natoms,1);
sum_quadrupole_3z2minusr2 = zeros(natoms,1);
for j=1:natoms
    Apoints=zeros(1,(2*delta_na + 1));
    for na = -delta_na:delta_na
        if periodicA
            Apoints(delta_na + na + 1) = round(mod((na + center_nabc(j,1)),totnumA) + 1);
        else
            Apoints(delta_na + na + 1) = round(na + center_nabc(j,1) + 1);
        end
    end
    Bpoints=zeros(1,(2*delta_nb + 1));
    for nb = -delta_nb:delta_nb
        if periodicB
            Bpoints(delta_nb + nb + 1) = round(mod((nb + center_nabc(j,2)),totnumB) + 1);
        else
            Bpoints(delta_nb + nb + 1) = round(nb + center_nabc(j,2) + 1);
        end
    end
    Cpoints=zeros(1,(2*delta_nc + 1));
    for nc = -delta_nc:delta_nc
        if periodicC
            Cpoints(delta_nc + nc + 1) = round(mod((nc + center_nabc(j,3)),totnumC) + 1);
        else
            Cpoints(delta_nc + nc + 1) = round(nc + center_nabc(j,3) + 1);
        end
    end
    for na = -delta_na:delta_na
        ka = Apoints(delta_na + na + 1);
        if ka < 1 || ka > totnumA
            continue
        end
        for nb = -delta_nb:delta_nb
            kb = Bpoints(delta_nb + nb + 1);
            if kb < 1 || kb > totnumB
                continue
            end
            for nc = -delta_nc:delta_nc
                kc = Cpoints(delta_nc + nc + 1);
                if kc < 1 || kc > totnumC
                    continue
                end               
                temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(j,1);
                temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(j,2);
                temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(j,3);
                distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3));
                index(1,j) = ceil(scalefactor*distance + zero_tolerance);
                if (index(1,j) <= nshells) && (total_pseudodensity(ka,kb,kc) > zero_tolerance)
                    local_density =  partial_density(j,index(1,j))*(valence_density(ka,kb,kc) + core_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc);
                    if core_pseudodensity(ka,kb,kc) >= zero_tolerance
                        local_core_density = partial_core_density(j,index(1,j))*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc);
                    else
                        local_core_density = 0.0;
                    end
                    local_valence_density = (local_density - local_core_density);
                    temp = normalization*local_valence_density*pixelvolume;
                    x = temp_vector(1);
                    y = temp_vector(2);
                    z = temp_vector(3);
                    sum_dipole_x(j) = sum_dipole_x(j) - x*temp;
                    sum_dipole_y(j) = sum_dipole_y(j) - y*temp;
                    sum_dipole_z(j) = sum_dipole_z(j) - z*temp;
                    sum_quadrupole_xy(j) = sum_quadrupole_xy(j) - x*y*temp;
                    sum_quadrupole_xz(j) = sum_quadrupole_xz(j) - x*z*temp;
                    sum_quadrupole_yz(j) = sum_quadrupole_yz(j) - y*z*temp;
                    sum_quadrupole_x2minusy2(j) = sum_quadrupole_x2minusy2(j) - (x^2 - y^2)*temp;
                    sum_quadrupole_3z2minusr2(j) = sum_quadrupole_3z2minusr2(j) - (3*z*z - (x^2 + y^2 + z^2))*temp;
                end
            end    
        end
    end
end
%
for j=1:natoms
    final_result(j,7) = sum_dipole_x(j) + occupancy_correction(j,3);
    final_result(j,8) = sum_dipole_y(j) + occupancy_correction(j,4);;
    final_result(j,9) = sum_dipole_z(j) + occupancy_correction(j,5);;
    final_result(j,10) = sqrt(final_result(j,7)^2 + final_result(j,8)^2 + final_result(j,9)^2);
    final_result(j,11) = sum_quadrupole_xy(j) + occupancy_correction(j,9);
    final_result(j,12) = sum_quadrupole_xz(j) + occupancy_correction(j,10);
    final_result(j,13) = sum_quadrupole_yz(j) + occupancy_correction(j,11);
    final_result(j,14) = sum_quadrupole_x2minusy2(j) + occupancy_correction(j,6) - occupancy_correction(j,7);
    final_result(j,15) = sum_quadrupole_3z2minusr2(j) + 2*occupancy_correction(j,8) - occupancy_correction(j,6) - occupancy_correction(j,7);
end
'Multipole analysis for each of the expansion sites.'
'XYZ coordinates, net charges, and multipoles are in atomic units. Dipoles and quadrupoles are for valence electrons.'
'center number, atomic number, x, y, z, net_charge, dipole_x, dipole_y, dipole_z, dipole_mag, Qxy, Qxz, Qyz, Q(x^2-y^2), Q(3z^2 - R^2)'
final_result

