% check the minimum pixel volume
if pixelvolume > minpixelvolume
   'The grid spacing is too coarse.' 
   'Please correct your input files and re-submit. Program will terminate.'
   flag = 1;
   break
end
% Initialize the partial density array with a guess
partial_density=zeros(natoms,nshells);
for i = 1:natoms
    for j = 1:nshells
        partial_density(i,j) = 0.05*((1 - (j-1)/(nshells - 1)))^3;
    end
end
% Check to see whether the grid spacing is adequate
index=zeros(1,natoms);
population=zeros(1,natoms);
sum_density = zeros(natoms,nshells);
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
                if (index(1,j) <= nshells)
                    sum_density(j,index(1,j)) = sum_density(j,index(1,j)) + partial_density(j,index(1,j));
                    sum_points(j,index(1,j)) = sum_points(j,index(1,j)) + 1;
                end
            end    
        end
    end
end 
flag = 0;
population=sum(sum_density')*pixelvolume;
for j = 1:natoms
    if population(j) - sum(population)/natoms > integration_tolerance
        'Integration volumes are not sufficiently accurate.'
        'Either your electron density input file(s) are too coarsely grained' 
        '    or they do not include large enough distance(s) along the nonperiodic direction(s) (if any).' 
        'Please correct your input files and re-submit. Program will terminate.'
        sum_density
        sum_points
        flag = 1;
        break
    end
end
if flag == 0
    'The grid spacing in your electron density input file is adequate.'
else
    break
end