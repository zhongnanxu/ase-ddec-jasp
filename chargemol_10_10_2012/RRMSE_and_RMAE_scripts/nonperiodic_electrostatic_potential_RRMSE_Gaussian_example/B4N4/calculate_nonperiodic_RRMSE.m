clear
computational_parameters
% this file is for Gaussian input
partial_charges_data
partial_charge_position = zeros(1,3);
natoms=length(partial_charge)
sq_pot_error= 0.0;
sq_pot_error_dipole= 0.0;
sq_pot = 0.0;
sum_pot = 0.0;
sum_pred_pot = 0.0;
sum_pred_pot_dipole = 0.0
sum_points=0;
coords=zeros(natoms,3);
sum_pot_error = 0.0;
sum_pot_error_dipole = 0.0;
% convert distances to atomic units
bohrperangstrom=1.889725989;
multiply_by = 1.0;
kcalmolhartree = 627.5095
vdW_radii
% open and read the cube file
inputfile='potential.cube'
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
unformatted_potential= cell2mat(data);
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
abinit_potential = zeros(totnumA,totnumB,totnumC);
for j = 1:M
    Q = j - 1;
    numC = mod(Q,totnumC);
    temp = (Q - numC)/totnumC;
    numB = mod(temp,totnumB);
    numA = (temp - numB)/totnumB;
    abinit_potential(numA+1,numB+1,numC+1) = unformatted_potential(round(j));
end
clear unformatted_potential
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
% accumulate the error sums
temp=0.0;
for ka = 1:totnumA
    ka
    for kb = 1:totnumB    
        for kc = 1:totnumC
            flag = 0;
            count = 0;
            x = ((ka-1)*boundary(1,1) + (kb-1)*boundary(2,1) + (kc-1)*boundary(3,1) + x_origin);
            y = ((ka-1)*boundary(1,2) + (kb-1)*boundary(2,2) + (kc-1)*boundary(3,2) + y_origin);
            z = ((ka-1)*boundary(1,3) + (kb-1)*boundary(2,3) + (kc-1)*boundary(3,3) + z_origin);
            pred_pot = 0.0;
            pred_pot_dipole = 0.0;
            for j = 1:natoms
                partial_charge_position(1) = coords(j,1);
                partial_charge_position(2) = coords(j,2);
                partial_charge_position(3) = coords(j,3);
                distance = sqrt((partial_charge_position(1) - x)^2 + (partial_charge_position(2) - y)^2 + (partial_charge_position(3) - z)^2);
                if distance < inner_multiplier*vdW_radius(atomic_number(j))
                    flag = 1;
                elseif distance < outer_multiplier*vdW_radius(atomic_number(j))
                    count = count + 1;
                end
                dipole_potential = (atomic_dipole(j,1)*(x - partial_charge_position(1)) + atomic_dipole(j,2)*(y - partial_charge_position(2)) + atomic_dipole(j,3)*(z - partial_charge_position(3)))/(distance^3*multiply_by^2);
                pred_pot = pred_pot + partial_charge(j)/(distance*multiply_by);
                pred_pot_dipole = pred_pot_dipole + partial_charge(j)/(distance*multiply_by) + dipole_potential;
            end
            if flag == 0 && count > 0.5
                sum_points = sum_points + 1;
                sq_pot_error= sq_pot_error + (abinit_potential(ka,kb,kc)-pred_pot)^2;
                sq_pot_error_dipole= sq_pot_error_dipole + (abinit_potential(ka,kb,kc)-pred_pot_dipole)^2;                
                sq_pot = sq_pot + (abinit_potential(ka,kb,kc))^2; 
                sum_pred_pot = sum_pred_pot + pred_pot;
                sum_pred_pot_dipole = sum_pred_pot_dipole + pred_pot_dipole;
                sum_pot = sum_pot + abinit_potential(ka,kb,kc);
                sum_pot_error = sum_pot_error + abinit_potential(ka,kb,kc)-pred_pot;
                sum_pot_error_dipole = sum_pot_error_dipole + abinit_potential(ka,kb,kc) - pred_pot_dipole;
            end    
        end    
    end
end
%
'The inner vdW radius multiplier was'
inner_multiplier
'The outer vdW radius multiplier was'
outer_multiplier
'The total number of valid grid points was'
sum_points
mean_pot = sum_pot/sum_points;
rms_sq_pot = sqrt(sq_pot/sum_points);
ave_pot_error = (sum_pot_error/sum_points);
ave_pot_error_dipole = (sum_pot_error_dipole/sum_points);
'The rms potential error without charges in kcal/mol is'
no_charges_rms_error = kcalmolhartree*sqrt(rms_sq_pot^2 - mean_pot^2)
'The rms potential error with partial charges in kcal/mol is'
RMSE_monopole = kcalmolhartree*sqrt(sq_pot_error/sum_points - ave_pot_error^2)
'The RRMSE value at monopole order'
RRMSE_monopole = RMSE_monopole/no_charges_rms_error
'The rms potential error with partial charges and atomic dipoles in kcal/mol is'
RMSE_dipole = kcalmolhartree*sqrt(sq_pot_error_dipole/sum_points - ave_pot_error_dipole^2)
'The RRMSE value at dipole order'
RRMSE_dipole = RMSE_dipole/no_charges_rms_error
quote
