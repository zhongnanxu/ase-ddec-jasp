% create the oxidation state reference densities
lower_oxidation_density = zeros(natoms,nshells);
upper_oxidation_density = zeros(natoms,nshells);
combined_oxidation_density=zeros(natoms,nshells);
% flag = 1 if a reference atomic density file is not found
flag = 0;
% read the lower and upper state reference densities
for j = 1:natoms
    atom_charge = ceil(atomic_number(j) - core_electrons(j) - valence_population(j));
    [temp_oxidation_density,flag] = oxidation_density(atomic_densities_directory,density_set_prefix,atomic_number(j),atom_charge,cutoff_radius,nshells);
    if flag == 1
       break
    end
    lower_oxidation_density(j,:) = temp_oxidation_density(:);
    atom_charge = atom_charge - 1;
    [temp_oxidation_density,flag] = oxidation_density(atomic_densities_directory,density_set_prefix,atomic_number(j),atom_charge,cutoff_radius,nshells);
    if flag == 1
       break
    end
    upper_oxidation_density(j,:) = temp_oxidation_density(:);
end 

if flag == 1
    break
end

% combine the spherical average densities with the reference oxidation states
for j = 1:natoms
    for k = 1:nshells
        f =valence_population(j) - floor(valence_population(j));
        combined_oxidation_density(j,k) = upper_oxidation_density(j,k)*f + lower_oxidation_density(j,k)*(1.0 - f);
    end
end
