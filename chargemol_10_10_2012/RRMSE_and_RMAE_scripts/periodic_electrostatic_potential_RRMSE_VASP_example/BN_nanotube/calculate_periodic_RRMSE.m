clear
computational_parameters
if vaspversion==4
   nheader=6
elseif vaspversion==5
   nheader=7
end
% convert distances to atomic units
bohrperangstrom=1.889725989;
eV_per_au = 27.2113961;
kcalmolhartree = 627.5095
vdW_radii
% Ewald summation parameters
Rdamp = 10.0*bohrperangstrom;
alpha = sqrt(pi)/Rdamp
cutoff_length = 3.0/alpha
% this file is for VASP input
computational_parameters
temp = 0.0;
% specify the input file with partial atomic charges
partial_charges_data
partial_charges_position = zeros(1,3);
natoms=length(partial_charge)
sq_pot_error= 0.0;
sq_pot_error_dipole = 0.0;
sq_pot = 0.0;
sum_pot = 0.0;
sum_pred_pot = 0.0;
sum_points=0;
sum_pot_error = 0.0;
sum_pot_error_dipole = 0.0;

% constant for erfc approximation
a = 8*(pi - 3)/(3*pi*(4-pi))
% open the potential file
inputfile='LOCPOT'
fid1 = fopen(inputfile,'r')
% read the data from the LOCPOT file
data = textscan(fid1, '%f',1,'headerlines',1);
latticevectorfactor = cell2mat(data);
clear data
data = textscan(fid1, '%f %f %f',3);
vectors = cell2mat(data)*latticevectorfactor*bohrperangstrom
clear data;
fclose(fid1);
fid1 = fopen(inputfile,'r');
for i = 1:nheader
    data = fgetl(fid1);
end    
atomspertype = str2num(data);
natoms=sum(atomspertype);
natomtypes=length(atomspertype);
clear data
textscan(fid1, '%s', 1);
data = textscan(fid1, '%f %f %f', natoms);
direct_coords=cell2mat(data)
clear data;
fclose(fid1);
fid1 = fopen(inputfile,'r');
data = textscan(fid1, '%f %f %f', 1,'headerlines',(2+nheader+natoms));
totnum = cell2mat(data)
clear data;
data = textscan(fid1,'%f');
raw_pot= cell2mat(data);
fclose(fid1);
clear fid1;
% read the vector and grid information
totnumA = totnum(1);
totnumB = totnum(2);
totnumC = totnum(3);
boundary(1,1)=vectors(1,1)/totnumA;
boundary(1,2)=vectors(1,2)/totnumA;
boundary(1,3)=vectors(1,3)/totnumA;
boundary(2,1)=vectors(2,1)/totnumB;
boundary(2,2)=vectors(2,2)/totnumB;
boundary(2,3)=vectors(2,3)/totnumB;
boundary(3,1)=vectors(3,1)/totnumC;
boundary(3,2)=vectors(3,2)/totnumC;
boundary(3,3)=vectors(3,3)/totnumC;
pixelvolume=det(boundary);
vector1=[vectors(1,1),vectors(1,2),vectors(1,3)];
vector2=[vectors(2,1),vectors(2,2),vectors(2,3)];
vector3=[vectors(3,1),vectors(3,2),vectors(3,3)];
% cutoffs
realsumlimitA = ceil(cutoff_length/sqrt(vectors(1,1)^2 + vectors(1,2)^2 + vectors(1,3)^2))
realsumlimitB = ceil(cutoff_length/sqrt(vectors(2,1)^2 + vectors(2,2)^2 + vectors(2,3)^2))
realsumlimitC = ceil(cutoff_length/sqrt(vectors(3,1)^2 + vectors(3,2)^2 + vectors(3,3)^2))
kmaxA = ceil(2*sqrt(vector1*vector1')/Rdamp)
kmaxB = ceil(2*sqrt(vector2*vector2')/Rdamp)
kmaxC = ceil(2*sqrt(vector3*vector3')/Rdamp)
for k=1:natoms
    coords(k,1)=direct_coords(k,1)*vector1(1) + direct_coords(k,2)*vector2(1)+ direct_coords(k,3)*vector3(1);
    coords(k,2)=direct_coords(k,1)*vector1(2) + direct_coords(k,2)*vector2(2)+ direct_coords(k,3)*vector3(2);
    coords(k,3)=direct_coords(k,1)*vector1(3) + direct_coords(k,2)*vector2(3)+ direct_coords(k,3)*vector3(3);
end
M = totnumA*totnumB*totnumC;
abinit_potential = zeros(totnumA,totnumB,totnumC);
for j = 1:M
    Q = j - 1;
    numA = mod(Q,totnumA);
    temp = (Q - numA)/totnumA;
    numB = mod(temp,totnumB);
    numC = (temp - numB)/totnumB;
    abinit_potential(numA+1,numB+1,numC+1) = (0.0 - raw_pot(round(j)))/eV_per_au;
end
clear raw_pot
% read the atomic numbers from the POTCAR file
inputfile = 'POTCAR'
fid1 = fopen(inputfile,'r');
nvalence_type=zeros(1,natomtypes);
atomic_number_type=zeros(1,natomtypes);
for i = 1:natomtypes
    tline = fgetl(fid1);
    data = textscan(fid1, '%f',1);
    nvalence_type(1,i)=cell2mat(data);
    clear data
    tline = fgetl(fid1);
    tline = fgetl(fid1);
    tline = fgetl(fid1);
    start_position=strfind(tline, '=')+1;
    end_position=strfind(tline, ':')-1;
    atomic_symbol = strtrim(tline(start_position:end_position));
    atomic_symbol_to_number;
    atomic_number_type(1,i) = z;
    for j = 1:100000
        tline = fgetl(fid1);
        if strncmpi(strtrim(tline),'End of Dataset',14)
            break
        end
    end
end  
fclose(fid1);
clear fid1
count = 0;
atomic_number = zeros(natoms,1);
for i = 1:natomtypes
    for j = 1:atomspertype(i)
        count = count + 1;
        atomic_number(count,1) = atomic_number_type(1,i);
    end
end
% Calculation of the reciprocal space vectors
reciprocal_vectors = 2*pi*((vectors')^-1)
%
% accumulate the error sums
for ka = 1:skip:totnumA
    ka
    for kb = 1:skip:totnumB    
        for kc = 1:skip:totnumC
            flag = 0;
            count = 0;
            x = ((ka-1)*boundary(1,1) + (kb-1)*boundary(2,1) + (kc-1)*boundary(3,1));
            y = ((ka-1)*boundary(1,2) + (kb-1)*boundary(2,2) + (kc-1)*boundary(3,2));
            z = ((ka-1)*boundary(1,3) + (kb-1)*boundary(2,3) + (kc-1)*boundary(3,3));
            pred_pot = 0.0;
            pred_pot_dipole = 0.0;
            dipole_pot = 0.0;
            B0 = 0.0;
            B1 = 0.0;
            for j = 1:natoms
                % The real summation part
                if flag == 1
                   break
                end 
                for repeata = -realsumlimitA:realsumlimitA
                    if flag == 1
                       break
                    end 
                    for repeatb = -realsumlimitB:realsumlimitB
                        if flag == 1
                           break
                        end
                        for repeatc = -realsumlimitC:realsumlimitC
                            if flag == 1
                               break
                            end
                            partial_charges_position(1) = coords(j,1) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1);
                            partial_charges_position(2) = coords(j,2) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2);
                            partial_charges_position(3) = coords(j,3) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3);
                            distance = sqrt((partial_charges_position(1) - x)^2 + (partial_charges_position(2) - y)^2 + (partial_charges_position(3) - z)^2);
                            if distance < inner_multiplier*vdW_radius(atomic_number(j))
                                flag = 1;
                                continue
                            elseif distance < outer_multiplier*vdW_radius(atomic_number(j))
                                count = count + 1;    
                            end
                            % erfc algorithm
                            s = distance*alpha; % argument for the erfc function
                            value_erfc = 1 - sqrt(1 - exp(-(s^2)*(4/pi + a*s^2)/(1+a*s^2))); % use this approximation of erfc
                            % value_erfc = erfc(s); % uncomment to use matlab builtin function
                            B0 = value_erfc/distance;
                            B1 = (B0 + 2*alpha*exp(-((alpha*distance)^2))/sqrt(pi))/(distance^2);
                            factor1 = atomic_dipole(j,1)*(x - partial_charges_position(1)) + atomic_dipole(j,2)*(y - partial_charges_position(2)) + atomic_dipole(j,3)*(z - partial_charges_position(3));
                            dipole_pot = factor1*B1;
                            pred_pot = pred_pot + B0*partial_charge(j);
                            pred_pot_dipole = pred_pot_dipole + B0*partial_charge(j) + dipole_pot; 
                        end
                    end
                end
            end
            if flag == 1 || count < 0.5
               continue
            end
            % The reciprocal space part
            temp1=0.0;
            temp2=0.0;
            for j = 1:natoms
                % The imagination summation part
                for repeata = -kmaxA:kmaxA
                    for repeatb = -kmaxB:kmaxB
                        for repeatc = -kmaxC:kmaxC
                            if (abs(repeata) + abs(repeatb) + abs(repeatc)) == 0
                               continue
                            end
                            delta_x = x - coords(j,1);
                            delta_y = y - coords(j,2);
                            delta_z = z - coords(j,3);
                            kx = repeata*reciprocal_vectors(1,1) + repeatb*reciprocal_vectors(2,1) + repeatc*reciprocal_vectors(3,1);
                            ky = repeata*reciprocal_vectors(1,2) + repeatb*reciprocal_vectors(2,2) + repeatc*reciprocal_vectors(3,2);
                            kz = repeata*reciprocal_vectors(1,3) + repeatb*reciprocal_vectors(2,3) + repeatc*reciprocal_vectors(3,3);
                            ksquared = kx*kx + ky*ky + kz*kz;
                            temp1 = cos(kx*delta_x + ky*delta_y + kz*delta_z);
                            temp2 = sin(kx*delta_x + ky*delta_y + kz*delta_z);
                            Ak = (4*pi*exp(-ksquared/(4*alpha*alpha)))/((pixelvolume*totnumA*totnumB*totnumC)*ksquared);
                            factor1 = (atomic_dipole(j,1)*kx + atomic_dipole(j,2)*ky + atomic_dipole(j,3)*kz);
                            pred_pot = pred_pot + temp1*partial_charge(j)*Ak;
                            pred_pot_dipole = pred_pot_dipole + temp1*partial_charge(j)*Ak + Ak*temp2*factor1;
                        end
                    end 
                end
            end            
            if flag == 0 && count > 0.5
                sum_points = sum_points + 1;
                sq_pot_error= sq_pot_error + (abinit_potential(ka,kb,kc)-pred_pot)^2;
                sum_pot_error = sum_pot_error + abinit_potential(ka,kb,kc)-pred_pot;
                sq_pot = sq_pot + (abinit_potential(ka,kb,kc))^2; 
                sum_pot = sum_pot + abinit_potential(ka,kb,kc);
                sum_pred_pot = sum_pred_pot + pred_pot;
                sum_pot_error_dipole = sum_pot_error_dipole + abinit_potential(ka,kb,kc) - pred_pot_dipole;
                sq_pot_error_dipole = sq_pot_error_dipole + (abinit_potential(ka,kb,kc)-pred_pot_dipole)^2;
            end
        end    
    end
end
%
'The number of valid grid points were'
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
'The Ewald summation parameters were:'
'alpha = '
alpha
'realsumlimitA = '
realsumlimitA
'realsumlimitB = '
realsumlimitB
'realsumlimitC = '
realsumlimitC
'kmaxA = '
kmaxA
'kmaxB = '
kmaxB
'kmaxC = '
kmaxC
'skip = '
skip
quote
