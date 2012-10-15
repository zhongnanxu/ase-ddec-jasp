% initialize some arrays
a = 0.0;
b = 0.0;
spherical_average_atomic_spin=zeros(natoms,nshells);
local_atomic_spin = 0.0;
proportional_spin = 0.0;
L = 0.0;
Ypsilon1 = zeros(totnumA,totnumB,totnumC);
Ypsilon2 = zeros(totnumA,totnumB,totnumC);
spin_population = zeros(natoms,1);
spin_pseudodensity = zeros(totnumA,totnumB,totnumC);
tot_spin_moment = 0.0;
% iterative solution for the atomic spin distributions
'Iteratively solving for the atomic spin distributions:'
for iter = 1:500
    sum_points = zeros(natoms,nshells);
    sum_spin_density = zeros(natoms,nshells);
    sum_Ypsilon1 = zeros(totnumA,totnumB,totnumC);
    sum_spin_pseudodensity = zeros(totnumA,totnumB,totnumC);
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
                    if (index(1,j) <= nshells) && (total_pseudodensity(ka,kb,kc) > zero_tolerance) && ((valence_density(ka,kb,kc) + core_density(ka,kb,kc)) > zero_tolerance)
                        sum_points(j,index(1,j)) = sum_points(j,index(1,j)) + 1;
                        local_atomic_density =  partial_density(j,index(1,j))*(valence_density(ka,kb,kc) + core_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc);
                        proportional_spin = partial_density(j,index(1,j))*spin_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc);
                        proportional_theta_scalar = calculate_theta_scalar(local_atomic_density,proportional_spin,Xi_lookup);
                        if iter == 1
                            local_atomic_spin = proportional_spin;
                            L = proportional_theta_scalar;
                            if local_atomic_spin < 0.0
                                local_atomic_spin_sign = -1;
                            else
                                local_atomic_spin_sign = 1;
                            end
                        else
                            a = spherical_average_density(j,index(1,j));
                            b = spherical_average_atomic_spin(j,index(1,j));
                            theta_scalar = calculate_theta_scalar(a,b,Xi_lookup);
                            L = Ypsilon1(ka,kb,kc) - Ypsilon2(ka,kb,kc) + pi*(spin_density(ka,kb,kc) - spin_pseudodensity(ka,kb,kc))/(valence_density(ka,kb,kc) + core_density(ka,kb,kc)) + (1-spin_ref_fraction)*theta_scalar + spin_ref_fraction*proportional_theta_scalar;
                            mag_L_vector = abs(L);
                            if L < 0.0
                                local_atomic_spin_sign = -1;
                            else
                                local_atomic_spin_sign = 1;
                            end
                            if mag_L_vector > pi
                               mag_local_atomic_spin_vector = local_atomic_density;
                            else
                                inv_Xi = fast_calculate_inverse_Xi(2*mag_L_vector,inverse_Xi_lookup);
                                mag_local_atomic_spin_vector = local_atomic_density*inv_Xi;
                            end
                            local_atomic_spin = mag_local_atomic_spin_vector*local_atomic_spin_sign;
                        end 
                        sum_spin_density(j,index(1,j)) = sum_spin_density(j,index(1,j)) + local_atomic_spin;
                        sum_Ypsilon1(ka,kb,kc) = sum_Ypsilon1(ka,kb,kc) + L*partial_density(j,index(1,j))/total_pseudodensity(ka,kb,kc);
                        sum_spin_pseudodensity(ka,kb,kc) = sum_spin_pseudodensity(ka,kb,kc) + local_atomic_spin;
                    end
                end    
            end
        end
    end
    Ypsilon1 = sum_Ypsilon1;
    spin_pseudodensity = sum_spin_pseudodensity;
    for j = 1:natoms
        for k = 1:nshells
            if sum_points(j,k) > 0
                spherical_average_atomic_spin(j,k) = sum_spin_density(j,k)/sum_points(j,k);
            else 
                spherical_average_atomic_spin(j,k) = 0.0;
            end
        end
    end
    iter
    old=spin_population;
    spin_population = sum(sum_spin_density')*pixelvolume;
    for j=1:natoms
        spin_population(j) = spin_population(j) + occupancy_correction(j,2);
    end
    spin_population
    tot_spin_moment = sum(spin_population)
    % compute occupations
    max_change=0.0;
    for j=1:natoms
        max_change = max(max_change,abs(spin_population(j)-old(j)));
    end
    max_change
    if ((max_change < spin_convergence_tolerance) && (iter > 6)) || ((max(max(abs(spin_population)),max_change) < spin_convergence_tolerance) && (iter > 1))
        generate_spin_magnetic_moment_file
        break
    end
    sum_Ypsilon2 = zeros(totnumA,totnumB,totnumC);
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
                    if (index(1,j) <= nshells) && (total_pseudodensity(ka,kb,kc) > zero_tolerance) && ((valence_density(ka,kb,kc) + core_density(ka,kb,kc)) > zero_tolerance)
                        local_atomic_density =  partial_density(j,index(1,j))*(valence_density(ka,kb,kc) + core_density(ka,kb,kc))/total_pseudodensity(ka,kb,kc);
                        proportional_spin = partial_density(j,index(1,j))*spin_density(ka,kb,kc)/total_pseudodensity(ka,kb,kc);
                        a = spherical_average_density(j,index(1,j));
                        b = spherical_average_atomic_spin(j,index(1,j));
                        theta_scalar = calculate_theta_scalar(a,b,Xi_lookup);
                        proportional_theta_scalar = calculate_theta_scalar(local_atomic_density,proportional_spin,Xi_lookup);
                        sum_Ypsilon2(ka,kb,kc) = sum_Ypsilon2(ka,kb,kc) + ((1-spin_ref_fraction)*theta_scalar + spin_ref_fraction*proportional_theta_scalar)*partial_density(j,index(1,j))/total_pseudodensity(ka,kb,kc);
                    end
                end    
            end
        end
    end
    Ypsilon2 = sum_Ypsilon2;
end
'Final spin populations:'
spin_population
'The total spin moment of the unit cell is   '
tot_spin_moment = sum(spin_population)



