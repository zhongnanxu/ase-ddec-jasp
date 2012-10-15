% read in the file containing the valence density
inputfile='valence_density.cube'
fid1 = fopen(inputfile,'r')
data = textscan(fid1, '%f %f %f %f',4,'headerlines',2);
parameters = cell2mat(data)
clear data;
fclose(fid1);
fid1 = fopen(inputfile,'r');
natoms=parameters(1,1);
atomic_number=zeros(natoms,1);
effective_nuclear_charge=zeros(natoms,1);
data = textscan(fid1, '%f %f %f %f %f', natoms,'headerlines',6);
temparray=cell2mat(data);
atomic_number(:,1) = temparray(:,1);
effective_nuclear_charge(:,1) = temparray(:,2);
coords(:,1) = temparray(:,3);
coords(:,2) = temparray(:,4);
coords(:,3) = temparray(:,5);
clear data;
fclose(fid1);
fid1 = fopen(inputfile,'r'); 
data = textscan(fid1,'%f','headerlines',(natoms+6));
raw_valence = cell2mat(data);
fclose(fid1);
clear fid1;
% read in the file containing the total density
flag = 0;
inputfile='total_density.cube'
fid1 = fopen(inputfile,'r')
if fid1 < 0
   core_available = 0;
else
   core_available = 1;
end
if (core_available == 1)
  data = textscan(fid1, '%f %f %f %f',4,'headerlines',2);
  if sum(sum(abs(parameters - cell2mat(data)))) > zero_tolerance
    flag = 1;
  end    
  clear data;
  fclose(fid1);
  fid1 = fopen(inputfile,'r');
  data = textscan(fid1, '%f %f %f %f %f', natoms,'headerlines',6);
  new_temparray=cell2mat(data);
  if sum(sum(abs(new_temparray - temparray))) > zero_tolerance
    flag = 1;
  end  
  clear data new_temparray;
  fclose(fid1);   
  if flag == 1
    'The core and valence density files do not contain the same lattice vectors, grid spacing, or atom positions.'
    'Please correct them and try again.'
    'Program will terminate.'
    break
  end 
  fid1 = fopen(inputfile,'r'); 
  data = textscan(fid1,'%f','headerlines',(natoms+6));
  raw_total = cell2mat(data);
  fclose(fid1);
end
clear fid1;
% construct the atomic density array:
unformatted_valence_density=raw_valence';
clear raw_valence
if core_available == 1
  unformatted_total_density=raw_total';
  clear raw_total
end
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
origin = [parameters(1,2), parameters(1,3), parameters(1,4)];
sum_negative_density = 0.0;
M = totnumA*totnumB*totnumC;
valence_density = zeros(totnumA,totnumB,totnumC);
core_density = zeros(totnumA,totnumB,totnumC);
for j = 1:M
    Q = j - 1;
    numC = mod(Q,totnumC);
    temp = (Q - numC)/totnumC;
    numB = mod(temp,totnumB);
    numA = (temp - numB)/totnumB;
    v = max(unformatted_valence_density(round(j)),0.0);
    if v*pixelvolume > pixel_integration_tolerance
        'Pixel with large density found. Program will correct density of this pixel.'
        correction_size = v*pixelvolume - pixel_integration_tolerance
        v = pixel_integration_tolerance/pixelvolume;
    end
    if unformatted_valence_density(round(j)) < 0.0
        sum_negative_density = sum_negative_density - pixelvolume*(v - unformatted_valence_density(round(j)));
    end
    valence_density(numA+1,numB+1,numC+1) = v;
    if core_available == 1
      core_density(numA+1,numB+1,numC+1) = max((unformatted_total_density(round(j)) - v),0.0);
    end
end
clear unformatted_valence_density
sum_negative_density
if core_available == 1
  missing_core = atomic_number - effective_nuclear_charge;
  clear unformatted_total_density
  % estimate the number of core electrons for each element
  core_electrons = zeros(natoms,1);
  if abs(sum(sum(sum(core_density))))*pixelvolume <= integration_tolerance
    core_electrons = atomic_number - effective_nuclear_charge;
  else
    for j=1:natoms
        if effective_nuclear_charge(j) == atomic_number(j)
           core_electrons(j) = num_core(atomic_number(j));
        else
           core_electrons(j) = atomic_number(j) - effective_nuclear_charge(j);
        end
    end
  end
else
  for j = 1:natoms
    core_electrons(j) = num_core(atomic_number(j));
  end
  missing_core = core_electrons; 
end
add_missing_core_density
%%%%%%%%%%%%%%
% open and read the spin density information
% read in the file containing the spin density
spin_flag = 0;
inputfile='spin_density.cube'
fid1 = fopen(inputfile,'r')
if fid1 > 0
    spin_available = 1;
else
    fid2 = fopen('spin_density_x.cube','r');
    fid3 = fopen('spin_density_y.cube','r');
    fid4 = fopen('spin_density_z.cube','r');
    if (fid2 > 0) && (fid3 > 0) && (fid4 > 0)
        spin_available = 2;
    else
        spin_available = 0;
    end
end
if spin_available == 1
    data = textscan(fid1, '%f %f %f %f',4,'headerlines',2);
    if sum(sum(abs(parameters - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end    
    clear data;
    fclose(fid1);
    fid1 = fopen(inputfile,'r');
    data = textscan(fid1, '%f %f %f %f %f', natoms,'headerlines',6);
    new_temparray=cell2mat(data);
    if sum(sum(abs(new_temparray - temparray))) > zero_tolerance
        spin_flag = 1;
        'The spin density file does not contain the same lattice vectors, grid spacing, or atom positions as the density files.'
        'Program will skip the spin moment analysis.'
    end  
    clear data new_temparray 
    fclose(fid1);
    fid1 = fopen(inputfile,'r'); 
    data = textscan(fid1,'%f','headerlines',(natoms+6));
    raw_spin = cell2mat(data);
    clear data
    fclose(fid1);
    clear fid1; 
    if spin_flag == 1
        spin_available = 0;
        break
    end
    unformatted_spin_density=raw_spin';
    clear raw_spin
    spin_density = zeros(totnumA,totnumB,totnumC);
    for j = 1:M
        Q = j - 1;
        numC = mod(Q,totnumC);
        temp = (Q - numC)/totnumC;
        numB = mod(temp,totnumB);
        numA = (temp - numB)/totnumB;
        p = valence_density(numA+1,numB+1,numC+1) + core_density(numA+1,numB+1,numC+1);
        if (p == 0.0) || (unformatted_spin_density(round(j)) == 0.0)
            spin_density(numA+1,numB+1,numC+1) = 0.0;
        else
            s = unformatted_spin_density(round(j));            
            spin_density(numA+1,numB+1,numC+1) = s*min(p/abs(s),1.0);
        end     
    end
    clear unformatted_spin_density
elseif spin_available == 2
    % load the spin_density_x.cube file
    inputfile = 'spin_density_x.cube'
    data = textscan(fid2, '%f %f %f %f',4,'headerlines',2);
    if sum(sum(abs(parameters - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end    
    clear data;
    fclose(fid2);
    fid2 = fopen(inputfile,'r');
    data = textscan(fid2, '%f %f %f %f %f', natoms,'headerlines',6);
    new_temparray=cell2mat(data);
    if sum(sum(abs(new_temparray - temparray))) > zero_tolerance
        spin_flag = 1;
        'The spin_density_x.cube file does not contain the same lattice vectors, grid spacing, or atom positions as the density files.'
        'Program will skip the spin moment analysis.'
    end  
    clear data new_temparray 
    fclose(fid2);
    fid2 = fopen(inputfile,'r'); 
    data = textscan(fid2,'%f','headerlines',(natoms+6));
    raw_spin_x = cell2mat(data);
    clear data
    fclose(fid2);
    clear fid2; 
    if spin_flag == 1
        spin_available = 0;
        break
    end
    unformatted_spin_density_x = raw_spin_x';
    clear raw_spin_x
    % load the spin_density_y.cube file
    inputfile = 'spin_density_y.cube'
    data = textscan(fid3, '%f %f %f %f',4,'headerlines',2);
    if sum(sum(abs(parameters - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end    
    clear data;
    fclose(fid3);
    fid3 = fopen(inputfile,'r');
    data = textscan(fid3, '%f %f %f %f %f', natoms,'headerlines',6);
    new_temparray=cell2mat(data);
    if sum(sum(abs(new_temparray - temparray))) > zero_tolerance
        spin_flag = 1;
        'The spin_density_y.cube file does not contain the same lattice vectors, grid spacing, or atom positions as the density files.'
        'Program will skip the spin moment analysis.'
    end  
    clear data new_temparray 
    fclose(fid3);
    fid3 = fopen(inputfile,'r'); 
    data = textscan(fid3,'%f','headerlines',(natoms+6));
    raw_spin_y = cell2mat(data);
    clear data
    fclose(fid3);
    clear fid3; 
    if spin_flag == 1
        spin_available = 0;
        break
    end
    unformatted_spin_density_xy= raw_spin_y';
    clear raw_spin_y    
    % load the spin_density_z.cube file
    inputfile = 'spin_density_z.cube'
    data = textscan(fid4, '%f %f %f %f',4,'headerlines',2);
    if sum(sum(abs(parameters - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end    
    clear data;
    fclose(fid4);
    fid4 = fopen(inputfile,'r');
    data = textscan(fid4, '%f %f %f %f %f', natoms,'headerlines',6);
    new_temparray=cell2mat(data);
    if sum(sum(abs(new_temparray - temparray))) > zero_tolerance
        spin_flag = 1;
        'The spin_density_z.cube file does not contain the same lattice vectors, grid spacing, or atom positions as the density files.'
        'Program will skip the spin moment analysis.'
    end  
    clear data new_temparray 
    fclose(fid4);
    fid4 = fopen(inputfile,'r'); 
    data = textscan(fid4,'%f','headerlines',(natoms+6));
    raw_spin_z = cell2mat(data);
    clear data
    fclose(fid4);
    clear fid4; 
    if spin_flag == 1
        spin_available = 0;
        break
    end
    unformatted_spin_density_z = raw_spin_z';
    clear raw_spin_z    
    % construct the spin density vector arrays
    spin_density_vector = zeros(totnumA,totnumB,totnumC,3);
    for j = 1:M
        Q = j - 1;
        numC = mod(Q,totnumC);
        temp = (Q - numC)/totnumC;
        numB = mod(temp,totnumB);
        numA = (temp - numB)/totnumB;
        p = valence_density(numA+1,numB+1,numC+1) + core_density(numA+1,numB+1,numC+1);
        s = sqrt(unformatted_spin_density_x(round(j))^2 + unformatted_spin_density_y(round(j))^2 + unformatted_spin_density_z(round(j))^2);
        if (p == 0.0) || (s == 0.0)
            spin_density_vector(numA+1,numB+1,numC+1,1) = 0.0;
            spin_density_vector(numA+1,numB+1,numC+1,2) = 0.0;
            spin_density_vector(numA+1,numB+1,numC+1,3) = 0.0;            
        else           
            spin_density_vector(numA+1,numB+1,numC+1,1) = unformatted_spin_density_x(round(j))*min(p/s,1.0);
            spin_density_vector(numA+1,numB+1,numC+1,2) = unformatted_spin_density_y(round(j))*min(p/s,1.0);
            spin_density_vector(numA+1,numB+1,numC+1,3) = unformatted_spin_density_z(round(j))*min(p/s,1.0);
        end     
    end
    clear unformatted_spin_density_x unformatted_spin_density_y unformatted_spin_density_z   
end
if flag == 1
    break
end
initialize_atomic_densities
if flag == 1
    break
end
% valence occupancy corrections can be added here if desired
occupancy_correction = zeros(natoms,11);    
