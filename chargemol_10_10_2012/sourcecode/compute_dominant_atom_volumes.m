% This file assigns each point in the unit cell to the atom with the
% highest density (as defined by dominant_atom_weight)
%
dominant_atom_points = zeros(totnumA,totnumB,totnumC);
maximum_dominant_atom_weight = zeros(totnumA,totnumB,totnumC);
active_points = 0;
num_changes = 0;
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
                    if dominant_atom_weight(j,index(1,j)) > maximum_dominant_atom_weight(ka,kb,kc)
                       dominant_atom_points(ka,kb,kc) = j;
                    end
                    maximum_dominant_atom_weight(ka,kb,kc) = max(maximum_dominant_atom_weight(ka,kb,kc),dominant_atom_weight(j,index(1,j)));
                end
            end    
        end
    end
end
clear temp_variable maximum_dominant_atom_weight dominant_atom_weight