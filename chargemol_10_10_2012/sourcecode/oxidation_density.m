function [output_oxidation_density,flag] = oxidation_density(atomic_densities_directory,density_set_prefix,Atomic_Num,atom_charge,cutoff_radius,nshells)
temp_oxidation_density = zeros(1,nshells);
flag = 0;
if atom_charge < 0
    sign_charge = -1;
else
    sign_charge = 1;
end    
for q = 0:sign_charge:atom_charge
    % construct the file name to read
    atom_no_string=num2str(Atomic_Num,'%03i');
    nelectron=Atomic_Num - q;
    nelectron_string=num2str(nelectron,'%03i');    
    cutoff_radius_string=num2str(cutoff_radius,'%03i');
    nshells_string=num2str(nshells,'%03i');
    combinedstring=strcat(atomic_densities_directory,density_set_prefix,'_',atom_no_string,'_',atom_no_string,'_',nelectron_string,'_',cutoff_radius_string,'_',nshells_string,'.txt');
    fid = fopen(combinedstring,'r');
    if fid < 0
        combinedstring
        'Could not find a suitable reference density. Program will terminate.'
        flag = 1;
        break
    end
    % read in the density
    data = textscan(fid, '%f',nshells,'headerlines',12);
    fclose(fid);
    temp = cell2mat(data);
    % compute the normalization factor for the reference density
    summed_oxidation_density = 0.0;
    summed_delta_oxidation_density = 0.0;
    for k = 1:nshells
        temp(k) = max(temp(k),1.0e-16);
        summed_oxidation_density = summed_oxidation_density + (k-1)^2*temp(k);
    end 
    %make the oxidation density monotonically decreasing
    temp1 = 1.0e-16;
    for k=nshells:-1:1
        if temp(k) > temp1
            temp1 = temp(k);
        else
            summed_delta_oxidation_density = summed_delta_oxidation_density + (temp1 - temp(k))*(k-1)^2;
            temp(k) = temp1;
        end
    end
    for k = 1:nshells
        temp(k) = temp(k)*summed_oxidation_density/(summed_oxidation_density + summed_delta_oxidation_density);
    end  
    % make the density difference monotonically increasing
    if q ~= 0
        for refine_step = 1:100
            summed_delta_oxidation_density = 0.0;
            delta_rho_min = 0.0;
            for k = nshells:-1:1
                delta_rho_min = max(delta_rho_min,sign_charge*(temp_oxidation_density(k) - temp(k)));
                density_shift = temp_oxidation_density(k) - sign_charge*delta_rho_min - temp(k);
                temp(k) = temp(k) + density_shift;
                summed_delta_oxidation_density = summed_delta_oxidation_density + density_shift*(k-1)^2;
            end
            for k = 1:nshells
                temp(k) = temp(k)*summed_oxidation_density/(summed_oxidation_density + summed_delta_oxidation_density);
            end
            if summed_delta_oxidation_density == 0
                break
            end    
        end    
    end
    temp_oxidation_density(:) = temp(:);
    clear data temp fid;
end
output_oxidation_density = temp_oxidation_density;