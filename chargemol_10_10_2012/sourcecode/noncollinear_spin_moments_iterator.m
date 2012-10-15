% initialize some arrays
a = 0.0;
b_vector = zeros(1,3);
spherical_average_atomic_spin_vector=zeros(natoms,nshells,3);
local_atomic_spin_vector = zeros(1,3);
proportional_spin_vector = zeros(1,3);
L_vector = zeros(1,3);
Ypsilon1_vector = zeros(totnumA,totnumB,totnumC,3);
Ypsilon2_vector = zeros(totnumA,totnumB,totnumC,3);
spin_population_vector = zeros(natoms,3);
spin_pseudodensity_vector = zeros(totnumA,totnumB,totnumC,3);
%
% iterative solution for the atomic spin distributions
'Iteratively solving for the atomic spin distributions:'
for iter = 1:500
    sum_points = zeros(natoms,nshells);
    sum_spin_density_vector = zeros(natoms,nshells,3);
    sum_Ypsilon1_vector = zeros(totnumA,totnumB,totnumC,3);
    sum_spin_pseudodensity_vector = zeros(totnumA,totnumB,totnumC,3);
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
                        for x=1:3; proportional_spin_vector(x) = partial_density(j,index(1,j))*spin_density_vector(ka,kb,kc,x)/total_pseudodensity(ka,kb,kc); end                        
                        proportional_theta_vector = calculate_theta_vector(local_atomic_density,proportional_spin_vector,Xi_lookup);
                        if iter == 1
                            local_atomic_spin_vector = proportional_spin_vector;
                            L_vector = proportional_theta_vector;
                        else
                            a = spherical_average_density(j,index(1,j));
                            for x=1:3; b_vector(x) = spherical_average_atomic_spin_vector(j,index(1,j),x); end
                            theta_vector = calculate_theta_vector(a,b_vector,Xi_lookup);
                            for x=1:3; L_vector(x) = Ypsilon1_vector(ka,kb,kc,x) - Ypsilon2_vector(ka,kb,kc,x) + pi*(spin_density_vector(ka,kb,kc,x) - spin_pseudodensity_vector(ka,kb,kc,x))/(valence_density(ka,kb,kc) + core_density(ka,kb,kc)) + (1-spin_ref_fraction)*theta_vector(x) + spin_ref_fraction*proportional_theta_vector(x); end
                            mag_L_vector = max(sqrt(L_vector*L_vector'),zero_tolerance^2);
                            if mag_L_vector > pi
                               mag_local_atomic_spin_vector = max(local_atomic_density,zero_tolerance^2);
                            else
                                inv_Xi = fast_calculate_inverse_Xi(2*mag_L_vector,inverse_Xi_lookup);
                                mag_local_atomic_spin_vector = max(local_atomic_density*inv_Xi,zero_tolerance^2);
                            end
                            local_atomic_spin_vector(:) = mag_local_atomic_spin_vector*L_vector(:)/mag_L_vector;
                        end 
                        for x=1:3; sum_spin_density_vector(j,index(1,j),x) = sum_spin_density_vector(j,index(1,j),x) + local_atomic_spin_vector(x); end
                        for x=1:3; sum_Ypsilon1_vector(ka,kb,kc,x) = sum_Ypsilon1_vector(ka,kb,kc,x) + L_vector(x)*partial_density(j,index(1,j))/total_pseudodensity(ka,kb,kc); end
                        for x=1:3; sum_spin_pseudodensity_vector(ka,kb,kc,x) = sum_spin_pseudodensity_vector(ka,kb,kc,x) + local_atomic_spin_vector(x); end
                    end
                end    
            end
        end
    end
    Ypsilon1_vector = sum_Ypsilon1_vector;
    spin_pseudodensity_vector = sum_spin_pseudodensity_vector;
    for j = 1:natoms
        for k = 1:nshells
            if sum_points(j,k) > 0
                spherical_average_atomic_spin_vector(j,k,:) = sum_spin_density_vector(j,k,:)/sum_points(j,k);
            else 
                spherical_average_atomic_spin_vector(j,k,1) = 0.0;
                spherical_average_atomic_spin_vector(j,k,2) = 0.0;
                spherical_average_atomic_spin_vector(j,k,3) = 0.0;
            end
        end
    end
    iter
    old=spin_population_vector;
    spin_population_vector = zeros(natoms,3);
    for j = 1:natoms
        for k = 1:nshells
            for x=1:3; spin_population_vector(j,x) = spin_population_vector(j,x) + sum_spin_density_vector(j,k,x)*pixelvolume; end
        end
    end 
    tot_spin_moment_vector = zeros(1,3);
    for j = 1:natoms
       for x=1:3; tot_spin_moment_vector(x) = tot_spin_moment_vector(x) + spin_population_vector(j,x); end
    end
    spin_population_vector = spin_population_vector
    tot_spin_moment_vector = tot_spin_moment_vector
    % compute occupations
    max_change=0.0;
    for j=1:natoms
        for i = 1:3
            max_change = max(max_change,abs(spin_population_vector(j,i)-old(j,i)));
        end
    end
    max_change
    if (max_change < spin_convergence_tolerance) && (iter > 6)
        generate_spin_magnetic_moment_file
        break
    end
    if iter == 1
       'Initial spin populations based on proportional partitioning of the spin density are:'
       spin_population_vector
       'The total spin moment of the unit cell is   '
       tot_spin_moment_vector
    end
    sum_Ypsilon2_vector = zeros(totnumA,totnumB,totnumC,3);
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
                        for x=1:3; proportional_spin_vector(x) = partial_density(j,index(1,j))*spin_density_vector(ka,kb,kc,x)/total_pseudodensity(ka,kb,kc); end                       
                        a = spherical_average_density(j,index(1,j));
                        for x=1:3; b_vector(x) = spherical_average_atomic_spin_vector(j,index(1,j),x); end
                        theta_vector = calculate_theta_vector(a,b_vector,Xi_lookup);
                        proportional_theta_vector = calculate_theta_vector(local_atomic_density,proportional_spin_vector,Xi_lookup);
                        for x=1:3; sum_Ypsilon2_vector(ka,kb,kc,x) = sum_Ypsilon2_vector(ka,kb,kc,x) + ((1-spin_ref_fraction)*theta_vector(x) + spin_ref_fraction*proportional_theta_vector(x))*partial_density(j,index(1,j))/total_pseudodensity(ka,kb,kc); end
                    end
                end    
            end
        end
    end
    Ypsilon2_vector = sum_Ypsilon2_vector;
end
'Final spin populations:'
spin_population_vector
'The total spin moment of the unit cell is   '
tot_spin_moment_vector = tot_spin_moment_vector


