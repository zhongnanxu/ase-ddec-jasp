% first determine the atomic positions and check the grid spacing
charge_center_positions
parallelpiped
check_grid_spacing
'Checking to see that all core electrons are accounted for:'
%
ref_core_density = zeros(natoms,nshells);
for j = 1:natoms
    if missing_core(j) == 0
        continue
    end
    % construct the file name to read
    string0='_';
    string1=num2str(atomic_number(j),'%03i');
    string2=num2str(missing_core(j),'%03i'); 
    string4=num2str(cutoff_radius, '%03i');
    string5=num2str(nshells, '%03i');
    string6='.txt';
    combinedstring=strcat(atomic_densities_directory,'core_',string1,string0,string1,string0,string2,string0,string4,string0,string5,string6);
    fid = fopen(combinedstring,'r');
    if fid == -1
         combinedstring
         'The atomic density file listed above does not exist. Program will terminate.'
         flag = 1
         break
    end
    % read in the density
    data = textscan(fid, '%f',nshells,'headerlines',12);
    temp = cell2mat(data);
    for k = 1:nshells
        ref_core_density(j,k) = temp(k);
    end
    fclose(fid);
    clear data temp fid;
end
if flag == 1
   break
end
for j=1:natoms
    if missing_core(j) == 0
        continue
    end
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
                   core_density(ka,kb,kc) = core_density(ka,kb,kc) + ref_core_density(j,index(1,j));
                end
            end    
        end
    end
end
if flag == 1
   break
end
'Finished the check for missing core electrons.' 

