index_min_penetration = ceil(rmin_cloud_penetration*nshells/cutoff_radius);
% calculation of electron cloud penetration terms
'Calculation of the electron cloud penetration terms'
slope=zeros(natoms,1);
intercept=zeros(natoms,1);
Rsquared=zeros(natoms,1);
for j = 1:natoms
    Sx = 0.0;
    Sxx = 0.0;
    Sxy = 0.0;
    Sy = 0.0;
    Syy = 0.0;
    n = 0;
    for k = index_min_penetration:nshells
        x = (k-0.5)*(cutoff_radius*bohrperangstrom/(100*nshells));
        if (spherical_average_density(j,k) + spherical_avg_core_density(j,k)) < zero_tolerance^1.5
           continue
        end
        y = log(spherical_average_density(j,k) + spherical_avg_core_density(j,k));
        Sx = Sx + x;
        Sxx = Sxx + x*x;
        Sxy = Sxy + x*y;
        Sy = Sy + y;
        Syy = Syy + y*y;
        n = n + 1;
    end
    atom_number = j;
    slope(j) = (n*Sxy - Sx*Sy)/(n*Sxx - Sx*Sx);
    intercept(j) = (Sy - slope(j)*Sx)/n;
    Rsquared(j) = (n*Sxy - Sx*Sy)^2/((n*Sxx - Sx*Sx)*(n*Syy - Sy*Sy));
end    
