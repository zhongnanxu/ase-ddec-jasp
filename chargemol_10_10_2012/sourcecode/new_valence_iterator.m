'Iteratively solving for the atomic partial charge distributions:'
wA_renormalization=ones(natoms,1);
old_old_old_density_change = 0.1;
old_old_density_change = 0.1; 
old_density_change = 0.1;
max_density_change = 0.1;
nonoverlapping_atom_tolerance = charge_convergence_tolerance/(10.0*wA_renormalization_max);
for iter = 1:2000
    total_pseudodensity = zeros(totnumA,totnumB,totnumC);
    ref_pseudodensity = zeros(totnumA,totnumB,totnumC);
    normalized_ref_pseudodensity = zeros(totnumA,totnumB,totnumC); 
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
                    if (index(1,j) <= nshells)
                        ref_pseudodensity(ka,kb,kc) = ref_pseudodensity(ka,kb,kc) + combined_oxidation_density(j,index(1,j));
                        normalized_ref_pseudodensity(ka,kb,kc) = normalized_ref_pseudodensity(ka,kb,kc) + normalized_oxidation_density(j,index(1,j));
                        total_pseudodensity(ka,kb,kc) = total_pseudodensity(ka,kb,kc) + partial_density(j,index(1,j));
                    end
                end    
            end
        end
    end
    sum_points = zeros(natoms,nshells);
    sum_valence_density = zeros(natoms,nshells);
    sum_density = zeros(natoms,nshells);
    sum_PtoWref = zeros(natoms,nshells);
    sum_sqrt_Wref = zeros(natoms,nshells);
    sum_wA_over_sqrt_Wref = zeros(natoms,nshells);
    wA_renormalization_slope=zeros(natoms,1); 
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
                        total_density = valence_density(ka,kb,kc) + core_density(ka,kb,kc);
                        local_density =  partial_density(j,index(1,j))*total_density/total_pseudodensity(ka,kb,kc);
                        sum_density(j,index(1,j)) = sum_density(j,index(1,j)) + local_density;
                        sum_points(j,index(1,j)) = sum_points(j,index(1,j)) + 1;
                        sum_PtoWref(j,index(1,j)) = sum_PtoWref(j,index(1,j)) + total_density/ref_pseudodensity(ka,kb,kc);
                        sum_sqrt_Wref(j,index(1,j)) = sum_sqrt_Wref(j,index(1,j)) + sqrt(normalized_ref_pseudodensity(ka,kb,kc));
                        sum_wA_over_sqrt_Wref(j,index(1,j)) = sum_wA_over_sqrt_Wref(j,index(1,j)) + normalized_oxidation_density(j,index(1,j))/sqrt(normalized_ref_pseudodensity(ka,kb,kc));
                        wA_renormalization_slope(j) = wA_renormalization_slope(j) + pixelvolume*(1.0 - partial_density(j,index(1,j))/total_pseudodensity(ka,kb,kc))*local_density;
                        if core_pseudodensity(ka,kb,kc) > zero_tolerance
                            local_core_density = partial_core_density(j,index(1,j))*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc);
                        else
                            local_core_density = 0.0;
                        end
                        local_valence_density = (local_density - local_core_density);
                        sum_valence_density(j,index(1,j)) = sum_valence_density(j,index(1,j)) + local_valence_density;
                    end
                end    
            end
        end
    end
    old_spherical_average_density = spherical_average_density; 
    for j = 1:natoms
        for k = 1:nshells
            if sum_points(j,k) > 0
                spherical_average_density(j,k) = sum_density(j,k)/sum_points(j,k);
                avg_PtoWref(j,k) = sum_PtoWref(j,k)/sum_points(j,k);
                avg_sqrt_Wref(j,k) = sum_sqrt_Wref(j,k)/sum_points(j,k);
                avg_wA_over_sqrt_Wref(j,k) = sum_wA_over_sqrt_Wref(j,k)/sum_points(j,k);
            else 
                spherical_average_density(j,k) = 0.0;
                avg_PtoWref(j,k) = 1.0;
                avg_sqrt_Wref(j,k) = sqrt(zero_tolerance);
                avg_wA_over_sqrt_Wref(j,k) = sqrt(zero_tolerance);
            end
        end
    end
    iter
    % compute occupations and population changes
    old=valence_population;
    old_old_old_change = old_old_change;
    old_old_change = old_change;
    old_change = max_change;
    old_old_old_density_change = old_old_density_change;
    old_old_density_change = old_density_change;
    old_density_change = max_density_change;
    valence_population = sum(sum_valence_density')*pixelvolume;
    for j = 1:natoms
        valence_population(j) = valence_population(j) + occupancy_correction(j,1);
    end
    normalization=nvalence/sum(valence_population)
    valence_population=normalization*valence_population
    change = valence_population - old;
    max_change = max(abs(change))
    max_density_change = max(max(abs(old_spherical_average_density - spherical_average_density)))
    temp_vector = [max_change,old_change,old_old_change,old_old_old_change,max_density_change,old_density_change,old_old_density_change,old_old_old_density_change];
    if  (max(temp_vector) < charge_convergence_tolerance) && (iter > 9)
        break
    end
    if iter == 1
       'Information for noniterative Hirshfeld method will be printed now.'
       local_multipole_moment_analysis
       'Hirshfeld analysis finished, calculation of iterative AIM will proceed.'
       Hirshfeld_population = valence_population
    end
    %update atomic densities
    if reference_weighting == 0
        partial_density = spherical_average_density;
    else
        first_integration_sum = zeros(natoms,1);
        weighted_points=zeros(natoms,1);
        sigma = zeros(natoms,nshells);
        for j = 1:natoms
            for k = 1:nshells
                ISA_part = spherical_average_density(j,k)^(1-reference_weighting);
                normalized_oxidation_density(j,k) = combined_oxidation_density(j,k)*avg_PtoWref(j,k);
                IH_part = normalized_oxidation_density(j,k)^reference_weighting;
                sigma(j,k) = ISA_part*IH_part;
                first_integration_sum(j) = first_integration_sum(j) + sigma(j,k)*sum_points(j,k);
                weighted_points(j) = weighted_points(j) + sum_points(j,k)*sqrt(sigma(j,k));
            end
        end
        % initialize the partial density array with an estimate
        partial_density = sigma;
        new_update_atomic_densities
        if flag == 1
          break
        end
        if iter > 3
          add_coefficient = 0.0;
          continue_flag = 1;
          for trial_num = 1:30
            if continue_flag == 0
               'Iterations to converge reshaping: '
               trial_num - 1
               break
            end
            continue_flag = 0;
            second_integration_sum = zeros(natoms,1);
            for j = 1:natoms
              partial_density(j,1) = max(partial_density(j,1),0.0);
              temp = partial_density(j,1);
              for k = 2:nshells
                constraint_term = 1 - (avg_wA_over_sqrt_Wref(j,k)/avg_sqrt_Wref(j,k))^2;
                % The decay exponent, convert to per radial shell
                exp_const = exp(-density_decaying_exponent*constraint_term*cutoff_radius*bohrperangstrom/(100*nshells));
                partial_density(j,k) = max(partial_density(j,k),0.0);
                if (sum_points(j,k) > 0) && (temp > 0)
                    partial_density(j,k) = min(partial_density(j,k),temp*exp_const);
                end
                temp = partial_density(j,k);
              end
              for k = 1:nshells
                second_integration_sum(j) = second_integration_sum(j) + partial_density(j,k)*sum_points(j,k);
              end
              % compute the constant to add to each point
              if weighted_points(j) > zero_tolerance
                 add_coefficient = (first_integration_sum(j) - second_integration_sum(j))/weighted_points(j);
              else
                 add_coefficient = 0.0
              end
              if abs(add_coefficient) > zero_tolerance
                continue_flag = 1;
                for k = 1:nshells
                  partial_density(j,k) = partial_density(j,k) + add_coefficient*sqrt(sigma(j,k));
                end
              end
            end
          end
          for j = 1:natoms  
            %    
            % the following constraint makes sure the number of valence electrons assigned to each atom is non-negative
            if wA_renormalization_slope(j) > nonoverlapping_atom_tolerance
               wA_renormalization(j) = max(1.0, (wA_renormalization(j) - valence_population(j)/wA_renormalization_slope(j)));
            else
               wA_renormalization(j) = 1.0;
            end
            for k = 1:nshells
               partial_density(j,k) = partial_density(j,k)*wA_renormalization(j);
            end
          end
        end
    end
    if flag == 1
       break
    end
end
if iter == 2000
    'Sorry, the calculation of atomic charges did not converge within 2000 iterations.'
end
