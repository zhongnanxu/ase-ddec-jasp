clear
% this file is for Gaussian input
% load constants
zero_tolerance = 0.01
truncated_icosahedron_vertices
vdW_radii
bohrperangstrom=1.889725989;
kcalmolhartree = 627.5095
% load spin data
atomic_spin_moments_data
partial_charge_position = zeros(1,3);
natoms=length(atomic_spin_moment)
coords=zeros(natoms,3);
inputfile='spin_density.cube'
fid1 = fopen(inputfile,'r')
data = textscan(fid1, '%f %f %f %f',4,'headerlines',2);
parameters = cell2mat(data)
clear data;
fclose(fid1);
fid1 = fopen(inputfile,'r');
data = textscan(fid1, '%f %f %f %f %f', natoms,'headerlines',6);
temparray=cell2mat(data);
atomic_number(:) = temparray(:,1);
coords(:,1) = temparray(:,3);
coords(:,2) = temparray(:,4);
coords(:,3) = temparray(:,5);
clear data;
clear temparray;
fclose(fid1);
fid1 = fopen(inputfile,'r'); 
data = textscan(fid1,'%f','headerlines',(natoms+6));
unformatted_spin_density= cell2mat(data);
fclose(fid1);
clear fid1;
% construct the potential array:
totnumA = parameters(2,1);
totnumB = parameters(3,1);
totnumC = parameters(4,1);
boundary(1,1)=parameters(2,2);
boundary(1,2)=parameters(2,3);
boundary(1,3)=parameters(2,4);
boundary(2,1)=parameters(3,2);
boundary(2,2)=parameters(3,3);
boundary(2,3)=parameters(3,4);
boundary(3,1)=parameters(4,2);
boundary(3,2)=parameters(4,3);
boundary(3,3)=parameters(4,4);
pixelvolume=det(boundary);
vector1=parameters(2,1)*[parameters(2,2),parameters(2,3),parameters(2,4)];
vector2=parameters(3,1)*[parameters(3,2),parameters(3,3),parameters(3,4)];
vector3=parameters(4,1)*[parameters(4,2),parameters(4,3),parameters(4,4)];
x_origin=parameters(1,2);
y_origin=parameters(1,3);
z_origin=parameters(1,4);
M = totnumA*totnumB*totnumC;
spin_density = zeros(totnumA,totnumB,totnumC);
for j = 1:M
    Q = j - 1;
    numC = mod(Q,totnumC);
    temp = (Q - numC)/totnumC;
    numB = mod(temp,totnumB);
    numA = (temp - numB)/totnumB;
    spin_density(numA+1,numB+1,numC+1) = unformatted_spin_density(round(j));
end
clear unformatted_spin_density
% positions of the charge centers in terms of na, nb, nc
b=inv(boundary');
center_nabc=zeros(natoms,3);
xyz = [0.00,0.00,0.00];
origin = [x_origin,y_origin,z_origin];
for k = 1:natoms
    xyz(:)=coords(k,:);
    temp_vector = b*(xyz' - origin');
    center_nabc(k,:) = round(temp_vector(:));
    center_shift(k,:) = boundary'*(temp_vector(:) - round(temp_vector(:)));
end
center_nabc
% compute the total spin moment for the valid integration points
total_moment=0.0;
for ka_int = 1:totnumA
    for kb_int = 1:totnumB
        for kc_int = 1:totnumC
            x_int = ((ka_int-1)*boundary(1,1) + (kb_int-1)*boundary(2,1) + (kc_int-1)*boundary(3,1) + x_origin);
            y_int = ((ka_int-1)*boundary(1,2) + (kb_int-1)*boundary(2,2) + (kc_int-1)*boundary(3,2) + y_origin);
            z_int = ((ka_int-1)*boundary(1,3) + (kb_int-1)*boundary(2,3) + (kc_int-1)*boundary(3,3) + z_origin);
            for j_int = 1:natoms
                partial_charge_position_int(1) = coords(j_int,1);
                partial_charge_position_int(2) = coords(j_int,2);
                partial_charge_position_int(3) = coords(j_int,3);
                distance_int = sqrt((partial_charge_position_int(1) - x_int)^2 + (partial_charge_position_int(2) - y_int)^2 + (partial_charge_position_int(3) - z_int)^2);
                if distance_int < 2.4*vdW_radius(atomic_number(j_int))
                    total_moment = total_moment + spin_density(ka_int,kb_int,kc_int)*pixelvolume;
                    break
                end
            end
        end
    end
end 
'The integration of the total magnetic moment was'
total_moment
% compute the magnetic field at each of the valid grid points
sum_MAE = 0.0;
sum_points=0;
sum_sq_field = 0.0;
count = 0;
sum_test_diff = 0.0;
for j = 1:natoms
    for shell = 1:3
        if shell == 1
           radial_multiplier = 3.0 + 1/6;
        elseif shell == 2
           radial_multiplier = 3.0 + 1/2;
        else
           radial_multiplier = 3.0 + 5/6;
        end
        for vert_num = 1:60
            flag = 0;
            % construct the grid point
            x = coords(j,1) + buckeyball_vertices(vert_num,1)*vdW_radius(atomic_number(j))*radial_multiplier;
            y = coords(j,2) + buckeyball_vertices(vert_num,2)*vdW_radius(atomic_number(j))*radial_multiplier;
            z = coords(j,3) + buckeyball_vertices(vert_num,3)*vdW_radius(atomic_number(j))*radial_multiplier;
            % determine if the grid point is valid (i.e. lies outside 3x van der Waals radii of all atoms)
            for k = 1:natoms
                partial_charge_position(1) = coords(k,1);
                partial_charge_position(2) = coords(k,2);
                partial_charge_position(3) = coords(k,3);
                distance = sqrt((partial_charge_position(1) - x)^2 + (partial_charge_position(2) - y)^2 + (partial_charge_position(3) - z)^2);
                if (distance + zero_tolerance) < radial_multiplier*vdW_radius(atomic_number(k))
                    flag = 1;
                    break
                end
            end
            % if flag == 0 then the integration point is valid
            if flag == 0
                count = count + 1;
                % compute the magnetic field at the valid grid point using the ab initio spin density
                magnetic_field_x = 0.0;
                magnetic_field_y = 0.0;
                magnetic_field_z = 0.0;
                test1 = 0.0;
                for ka_int = 1:totnumA
                    for kb_int = 1:totnumB
                        for kc_int = 1:totnumC
                            x_int = ((ka_int-1)*boundary(1,1) + (kb_int-1)*boundary(2,1) + (kc_int-1)*boundary(3,1) + x_origin);
                            y_int = ((ka_int-1)*boundary(1,2) + (kb_int-1)*boundary(2,2) + (kc_int-1)*boundary(3,2) + y_origin);
                            z_int = ((ka_int-1)*boundary(1,3) + (kb_int-1)*boundary(2,3) + (kc_int-1)*boundary(3,3) + z_origin);
                            for j_int = 1:natoms
                                partial_charge_position_int(1) = coords(j_int,1);
                                partial_charge_position_int(2) = coords(j_int,2);
                                partial_charge_position_int(3) = coords(j_int,3);
                                distance_int = sqrt((partial_charge_position_int(1) - x_int)^2 + (partial_charge_position_int(2) - y_int)^2 + (partial_charge_position_int(3) - z_int)^2);
                                if distance_int < 2.4*vdW_radius(atomic_number(j_int))
                                    delta_X = x - x_int;
                                    delta_Y = y - y_int;
                                    delta_Z = z - z_int;
                                    delta_distance = sqrt(delta_X^2 + delta_Y^2 + delta_Z^2);
                                    magnetic_field_x = magnetic_field_x + 3*delta_X*delta_Z*spin_density(ka_int,kb_int,kc_int)*pixelvolume/delta_distance^5;
                                    magnetic_field_y = magnetic_field_y + 3*delta_Y*delta_Z*spin_density(ka_int,kb_int,kc_int)*pixelvolume/delta_distance^5;
                                    magnetic_field_z = magnetic_field_z + 3*delta_Z*delta_Z*spin_density(ka_int,kb_int,kc_int)*pixelvolume/delta_distance^5 - spin_density(ka_int,kb_int,kc_int)*pixelvolume/delta_distance^3;
                                    test1 = test1 + spin_density(ka_int,kb_int,kc_int)*pixelvolume;
                                    break
                                end
                            end
                        end
                    end
                end
                % compute the magnetic field from the atomic spin moments 
                pred_field_x = 0.0;
                pred_field_y = 0.0;
                pred_field_z = 0.0;
                pred_scalar_potential = 0.0;
                test2 = 0.0;
                for k = 1:natoms
                    delta_X = x - coords(k,1);
                    delta_Y = y - coords(k,2);
                    delta_Z = z - coords(k,3);                
                    distance = sqrt(delta_X^2 + delta_Y^2 + delta_Z^2);
                    pred_field_x = pred_field_x + 3*delta_X*delta_Z*atomic_spin_moment(k)/distance^5;
                    pred_field_y = pred_field_y + 3*delta_Y*delta_Z*atomic_spin_moment(k)/distance^5;
                    pred_field_z = pred_field_z + 3*delta_Z*delta_Z*atomic_spin_moment(k)/distance^5 - atomic_spin_moment(k)/distance^3;
                    test2 = test2 + atomic_spin_moment(k);
                end
                sum_test_diff = sum_test_diff + abs(test1 - test2);
                % integrate the abs difference between the predicted and ab initio magnetic fields over the volume between 3x and 4x vdW surfaces
                % the square of the radial_multiplier is proportional to the integration voluem for the grid point
                sum_points = sum_points + (radial_multiplier*vdW_radius(atomic_number(j)))^2;
                sum_MAE = sum_MAE + (radial_multiplier*vdW_radius(atomic_number(j)))^2*sqrt((pred_field_x -magnetic_field_x)^2 +  (pred_field_y -magnetic_field_y)^2 + (pred_field_z -magnetic_field_z)^2);
                sum_sq_field = sum_sq_field + (radial_multiplier*vdW_radius(atomic_number(j)))^2*sqrt(magnetic_field_x^2 + magnetic_field_y^2 + magnetic_field_z^2); 
            end
        end
    end  
end
%
'The total number of valid grid points were'
count
'The mean absolute error without atomic spin moments was'
sum_sq_field/sum_points
'The mean absolute error with atomic spin moments was'
MAE = sum_MAE/sum_points
if sum_sq_field > 0.0
    'The ratio of the two was'
    MAE*sum_points/sum_sq_field
end  
'The integration of the total magnetic moment was'
total_moment
'The average difference in test conditions was'
sum_test_diff/count
quote
