% initialize some things here
num_periodic_directions = periodicA + periodicB + periodicC;
% open and read the xsf file containing the valence density
fid1 = fopen(xsf_inputfile,'r');
if fid1 <= 0
    'Could not find input xsf file. Program will terminate.'
    flag=1
    break
else
    fileInfo = dir(xsf_inputfile);
    fileSize = fileInfo.bytes;
    if fileSize < 10
       'Input xsf file appears to be empty. Program will terminate.'
       flag=1
       break
    end 
end
spins_loaded=0;
while ~feof(fid1)
   data = fgetl(fid1);
   if (length(lower(strtrim(data))) == length('primvec'))
       if (lower(strtrim(data)) == 'primvec')
           data = textscan(fid1, '%f %f %f',num_periodic_directions);
           periodic_vectors = cell2mat(data)*bohrperangstrom
           found_primvec = 1;
           clear data;
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('primcoord'))
       if (lower(strtrim(data)) == 'primcoord')
           clear data;
           natoms = cell2mat(textscan(fid1, '%u',1))
           junk = textscan(fid1, '%u',1);
           atomic_number=zeros(natoms,1);
           data = textscan(fid1, '%f %f %f %f',natoms);
           temparray=cell2mat(data);
           atomic_number(:,1) = int16(temparray(:,1))
           coords(:,1) = temparray(:,2)*bohrperangstrom;
           coords(:,2) = temparray(:,3)*bohrperangstrom;
           coords(:,3) = temparray(:,4)*bohrperangstrom
           clear data;
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('atoms'))
       if (lower(strtrim(data)) == 'atoms')
           clear data
           if num_periodic_directions == 3
              'You indicated a periodic system, but the grid file does not contain the unit cell information.'
              'Program will terminate. Please prepare a new input file and try again.'
              break
           end    
           for i = 1:10000
              data = fgetl(fid1);
              temp = str2num(data);
              if isempty(temp)
                 clear temp;
                 break
              end
              atomic_number(i,1) = temp(1);
              coords(i,1) = temp(2)*bohrperangstrom;
              coords(i,2) = temp(3)*bohrperangstrom;
              coords(i,3) = temp(4)*bohrperangstrom;          
              clear data
              continue
           end
       end        
   end
   if (length(lower(strtrim(data))) == length('begin_block_datagrid_3d'))
       if (lower(strtrim(data)) == 'begin_block_datagrid_3d')
           clear data
           tline = fgetl(fid1);
           while ischar(tline)              
              if (length(lower(strtrim(tline))) == length('begin_datagrid_3d_rho:spin_1'))
                  if (lower(strtrim(tline)) == 'begin_datagrid_3d_rho:spin_1')
                      data = str2num(fgetl(fid1));
                      if isempty(data)
                         clear data;
                         'Error reading xsf file. Program will terminate.'
                         flag = 1;
                         break
                      end
                      totnumA = data(1) - 1;
                      totnumB = data(2) - 1;
                      totnumC = data(3) - 1;
                      clear data;
                      data = textscan(fid1, '%f %f %f',1);
                      origin = cell2mat(data)*bohrperangstrom;
                      clear data;
                      data = textscan(fid1, '%f %f %f',3);
                      vectors = cell2mat(data)*bohrperangstrom;
                      clear data;
                      if (num_periodic_directions == 3) && found_primvec
                          if sum(sum(abs(periodic_vectors - vectors))) > zero_tolerance
                              flag = 1;
                              'The periodic vectors and the grid vectors must match for a periodic system.'
                              'Program will terminate. Please prepare a new input file and try again.'
                              break
                          end
                      end
                      npoints = (totnumA+1)*(totnumB+1)*(totnumC+1);
                      data = textscan(fid1,'%f', npoints);
                      raw_spin_1= cell2mat(data);
                      spins_loaded = spins_loaded + 1
                      clear data
                  elseif (lower(strtrim(tline)) == 'begin_datagrid_3d_rho:spin_2')
                      data = str2num(fgetl(fid1));
                      if isempty(data)
                         clear data;
                         'Error reading xsf file. Program will terminate.'
                         flag = 1;
                         break
                      end
                      totnumA = data(1) - 1;
                      totnumB = data(2) - 1;
                      totnumC = data(3) - 1;
                      clear data;
                      data = textscan(fid1, '%f %f %f',1);
                      origin = cell2mat(data)*bohrperangstrom;
                      clear data;
                      data = textscan(fid1, '%f %f %f',3);
                      vectors = cell2mat(data)*bohrperangstrom;
                      clear data;
                      if (num_periodic_directions == 3) && found_primvec
                          if sum(sum(abs(periodic_vectors - vectors))) > zero_tolerance
                              flag = 1;
                              'The periodic vectors and the grid vectors must match for a periodic system.'
                              'Program will terminate. Please prepare a new input file and try again.'
                              break
                          end
                      end
                      npoints = (totnumA+1)*(totnumB+1)*(totnumC+1);
                      data = textscan(fid1,'%f', npoints);
                      raw_spin_2= cell2mat(data);
                      spins_loaded = spins_loaded + 1
                      clear data    
                  elseif (lower(strtrim(tline)) == 'begin_datagrid_3d_rho:spin_3')
                      data = str2num(fgetl(fid1));
                      if isempty(data)
                         clear data;
                         'Error reading xsf file. Program will terminate.'
                         flag = 1;
                         break
                      end
                      totnumA = data(1) - 1;
                      totnumB = data(2) - 1;
                      totnumC = data(3) - 1;
                      clear data;
                      data = textscan(fid1, '%f %f %f',1);
                      origin = cell2mat(data)*bohrperangstrom;
                      clear data;
                      data = textscan(fid1, '%f %f %f',3);
                      vectors = cell2mat(data)*bohrperangstrom;
                      clear data;
                      if (num_periodic_directions == 3) && found_primvec
                          if sum(sum(abs(periodic_vectors - vectors))) > zero_tolerance
                              flag = 1;
                              'The periodic vectors and the grid vectors must match for a periodic system.'
                              'Program will terminate. Please prepare a new input file and try again.'
                              break
                          end
                      end
                      npoints = (totnumA+1)*(totnumB+1)*(totnumC+1);
                      data = textscan(fid1,'%f', npoints);
                      raw_spin_3 = cell2mat(data);
                      spins_loaded = spins_loaded + 1
                      clear data
                  elseif  (lower(strtrim(tline)) == 'begin_datagrid_3d_rho:spin_4')
                      data = str2num(fgetl(fid1));
                      if isempty(data)
                         clear data;
                         'Error reading xsf file. Program will terminate.'
                         flag = 1;
                         break
                      end
                      totnumA = data(1) - 1;
                      totnumB = data(2) - 1;
                      totnumC = data(3) - 1;
                      clear data;
                      data = textscan(fid1, '%f %f %f',1);
                      origin = cell2mat(data)*bohrperangstrom;
                      clear data;
                      data = textscan(fid1, '%f %f %f',3);
                      vectors = cell2mat(data)*bohrperangstrom;
                      clear data;
                      if (num_periodic_directions == 3) && found_primvec
                          if sum(sum(abs(periodic_vectors - vectors))) > zero_tolerance
                              flag = 1;
                              'The periodic vectors and the grid vectors must match for a periodic system.'
                              'Program will terminate. Please prepare a new input file and try again.'
                              break
                          end
                      end
                      npoints = (totnumA+1)*(totnumB+1)*(totnumC+1);
                      data = textscan(fid1,'%f', npoints);
                      raw_spin_4 = cell2mat(data);
                      spins_loaded = spins_loaded + 1
                      clear data
                  end
              end
              tline = fgetl(fid1);
           end
           fclose(fid1);
           clear fid1;
           break
       end
   end    
end
clear tline
totnumA
totnumB
totnumC
origin
vectors
temp = size(coords);
natoms = temp(1);
clear temp
if flag == 1
    break
end
% construct the spin density matrices from the individual spin components
if spins_loaded == 1
    spin_available = 0;
    'No spin moment analysis will be performed.'
    raw_valence = raw_spin_1;
    clear raw_spin_1;
elseif spins_loaded == 2
    spin_available = 1;
    'Collinear spin moment analysis will be performed.'
    raw_valence = raw_spin_1 + raw_spin_2;
    raw_spin = raw_spin_1 - raw_spin_2;
    clear raw_spin_1 raw_spin_2;    
elseif spins_loaded == 4
    'Warning: Non-collinear spin feature is still being tested. Results may not be accurate!'
    spin_available = 2;
    'Noncollinear spin moment analysis will be performed.'
    raw_valence = raw_spin_1 + raw_spin_2;
    raw_spin_x = 2.0*raw_spin_3;
    raw_spin_y = -2.0*raw_spin_4;
    raw_spin_z = raw_spin_1 - raw_spin_2;
    clear raw_spin_1 raw_spin_2 raw_spin_3 raw_spin_4; 
else
    'Error reading the density and spin information in the xsf file.'
    'Program will terminate.'
    flag = 1;
    break
end      
% initialize the valence density arrays
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
sum_negative_density = 0.0;
valence_density = zeros(totnumA,totnumB,totnumC);
for j = 1:npoints
    Q = j - 1;
    numA = mod(Q,(totnumA+1));
    temp = (Q - numA)/(totnumA + 1);
    numB = mod(temp,(totnumB+1));
    numC = (temp - numB)/(totnumB + 1);
    v=max(raw_valence(round(j)),0.0);
    if v*pixelvolume > pixel_integration_tolerance
        'Pixel with large density found. Program will correct density of this pixel.'
        correction_size = v*pixelvolume - pixel_integration_tolerance
        v = pixel_integration_tolerance/pixelvolume;
    end
    if raw_valence(round(j)) < 0.0
        sum_negative_density = sum_negative_density - (pixelvolume*v - raw_valence(round(j))/(totnumA*totnumB*totnumC));
    end
    numA = mod(numA,totnumA);
    numB = mod(numB,totnumB);
    numC = mod(numC,totnumC);
    valence_density(numA+1,numB+1,numC+1) = v; % valence density   
end
clear raw_core raw_valence
sum_negative_density
for j = 1:natoms
  core_electrons(j) = num_core(atomic_number(j));
end
core_available = 0;
core_density = zeros(totnumA,totnumB,totnumC);
% assign the number of core electrons and effective nuclear charge if core fitting is not performed
if core_available == 0
  for j = 1:natoms
    core_electrons(j) = num_core(atomic_number(j));
  end
  missing_core=core_electrons;
end
add_missing_core_density
% process the spin densities
if spin_available == 1
    spin_density = zeros(totnumA,totnumB,totnumC);
    for j = 1:npoints
        Q = j - 1;
        numA = mod(Q,(totnumA+1));
        temp = (Q - numA)/(totnumA + 1);
        numB = mod(temp,(totnumB+1));
        numC = (temp - numB)/(totnumB + 1);
        numA = mod(numA,totnumA);
        numB = mod(numB,totnumB);
        numC = mod(numC,totnumC);
        p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1);
        if (abs(p) == 0.0) || (raw_spin(round(j)) == 0.0)
            spin_density(numA+1,numB+1,numC+1) = 0.0;
        else
            s = raw_spin(round(j));
            spin_density(numA+1,numB+1,numC+1) = s*min(abs(p)/abs(s),1.0);
        end    
    end
    checkspin=sum(sum(sum(spin_density)))*pixelvolume
    clear raw_spin   
elseif spin_available == 2
    spin_density_vector = zeros(totnumA,totnumB,totnumC,3);
    for j = 1:npoints
        Q = j - 1;
        numA = mod(Q,(totnumA+1));
        temp = (Q - numA)/(totnumA + 1);
        numB = mod(temp,(totnumB+1));
        numC = (temp - numB)/(totnumB + 1);
        numA = mod(numA,totnumA);
        numB = mod(numB,totnumB);
        numC = mod(numC,totnumC);
        p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1);
        s = sqrt(raw_spin_x(round(j))^2 + raw_spin_y(round(j))^2 + raw_spin_z(round(j))^2);
        if (p == 0.0) || (s == 0.0)
            spin_density_vector(numA+1,numB+1,numC+1,1) = 0.0;
            spin_density_vector(numA+1,numB+1,numC+1,2) = 0.0;
            spin_density_vector(numA+1,numB+1,numC+1,3) = 0.0;
        else           
            spin_density_vector(numA+1,numB+1,numC+1,1) = raw_spin_x(round(j))*min(p/s,1.0);
            spin_density_vector(numA+1,numB+1,numC+1,2) = raw_spin_y(round(j))*min(p/s,1.0);
            spin_density_vector(numA+1,numB+1,numC+1,3) = raw_spin_z(round(j))*min(p/s,1.0);
        end    
    end
    clear raw_spin_x raw_spin_y raw_spin_z    
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
