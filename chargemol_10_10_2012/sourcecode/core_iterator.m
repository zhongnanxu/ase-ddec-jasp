spherical_avg_core_density= zeros(natoms,nshells);
% The decay exponent is 2 per bohr, convert to per radial shell
exp_const = exp(-2.0*cutoff_radius*bohrperangstrom/(100*nshells));
'Iteratively solving for the core charge distributions:'
for iter = 1:200
    core_pseudodensity(:,:,:) = zeros(totnumA,totnumB,totnumC);
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
                    if index(1,j) <= nshells
                        core_pseudodensity(ka,kb,kc) = core_pseudodensity(ka,kb,kc) + partial_core_density(j,index(1,j));
                    end
                end    
            end
        end
    end 
    sum_core_density = zeros(natoms,nshells);
    sum_points = zeros(natoms,nshells);
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
                    if (index(1,j) <= nshells) && (core_pseudodensity(ka,kb,kc) > zero_tolerance)
                        sum_core_density(j,index(1,j)) = sum_core_density(j,index(1,j)) + partial_core_density(j,index(1,j))*core_density(ka,kb,kc)/core_pseudodensity(ka,kb,kc);
                        sum_points(j,index(1,j)) = sum_points(j,index(1,j)) + 1;
                    end
                end    
            end
        end
    end 
    for j = 1:natoms
        for k = 1:nshells
            if sum_points(j,k) > 0
                spherical_avg_core_density(j,k) = sum_core_density(j,k)/sum_points(j,k);
            else 
                spherical_avg_core_density(j,k) = 0.0;
            end
        end
    end
    iter
    old=core_population;
    core_population=pixelvolume*(sum(sum_core_density'));
    % compute core occupations
    change=0.0;
    for j=1:natoms
        change = max(change,abs(core_population(j)-old(j)));
    end
    change
    if (change < charge_convergence_tolerance) && (iter > 5)
        break
    end
    %update atomic core density factors
    for j = 1:natoms
        partial_core_density(j,1) = spherical_avg_core_density(j,1);
        temp = partial_core_density(j,1);
        for k = 2:nshells
            if (sum_points(j,k) > 0) && (temp > 0)
                partial_core_density(j,k) = min(spherical_avg_core_density(j,k),temp*exp_const);
            else
                partial_core_density(j,k) = spherical_avg_core_density(j,k);
            end    
            temp = partial_core_density(j,k);   
        end
    end  
end

