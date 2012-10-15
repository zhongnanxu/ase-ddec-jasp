% set the number of headerlines before the atom types line
if vaspversion==4
   nheader=6;
elseif vaspversion==5
   nheader=7;
end 
% open and read the valence density file
inputfile='AECCAR2'
fid1 = fopen(inputfile,'r');
test = fid1;
if test < 0
   inputfile='CHGCAR'
   clear fid1
   fid1 = fopen(inputfile,'r');
   if fid1 < 0
     'Could not find input valence density file. Program will terminate.'
     flag = 1
     break
   else
     'Could not find the AECCAR2 file so valence (pseudo)density will be read from the CHGCAR file.'
     'For best results, please use the PAW potentials together with the INCAR keyword LAECHG = .TRUE. to generate the AECCAR2 file.'
     'Calculation will proceed.'
   end
end
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
data = textscan(fid1, '%f %f %f', 1);
totnum = cell2mat(data);
totnumA = totnum(1);
totnumB = totnum(2);
totnumC = totnum(3);
clear data;
npoints = totnumA*totnumB*totnumC;
data = textscan(fid1,'%f',npoints);
raw_valence= cell2mat(data);
fclose(fid1);
clear fid1;

% open and read the core density file
flag = 0;
inputfile='AECCAR0'
fid1 = fopen(inputfile,'r');
if fid1 > 0
   core_available = 1;
else
   core_available = 0;
end
core_flag = 0;
if (core_available == 1)
  data = textscan(fid1, '%f',1,'headerlines',1);
  if sum(sum(abs(latticevectorfactor - cell2mat(data)))) > zero_tolerance
     core_flag = 1;
  end
  clear data
  data = textscan(fid1, '%f %f %f',3);
  if sum(sum(abs(vectors - cell2mat(data)*latticevectorfactor*bohrperangstrom))) > zero_tolerance
    core_flag = 1;
  end
  clear data;
  fclose(fid1);
  fid1 = fopen(inputfile,'r');
  for i = 1:nheader
     data = fgetl(fid1);
  end   
  if sum(sum(abs(atomspertype - str2num(data)))) > zero_tolerance
    core_flag = 1;
  end
  clear data
  textscan(fid1, '%s', 1);
  data = textscan(fid1, '%f %f %f', natoms);
  if sum(sum(abs(direct_coords - cell2mat(data)))) > zero_tolerance
    core_flag = 1;
  end
  clear data;
  data = textscan(fid1, '%f %f %f', 1);
  if sum(sum(abs(totnum - cell2mat(data)))) > zero_tolerance
    core_flag = 1;
  end
  clear data;  
  data = textscan(fid1,'%f');
  raw_core= cell2mat(data);
  fclose(fid1);
  clear fid1;
  clear data;
end
if core_flag == 1
    'The core density (AECCAR0) file does not contain the same lattice vectors, grid spacing, or atom positions as the valence density (AECCAR2) file.'
    'This can occur if the AECCAR0 and AECCAR2 were produced from a geometry optimization rather than from a single-point calculation.'
    'The program will generate its own core densities.'
    core_available = 0;
end  
% open and read the valence pseudodensity
inputfile='CHG'
fid1 = fopen(inputfile,'r');
if fid1 <= 0
   inputfile='CHG_noncollinear'
end
clear fid1
fid1 = fopen(inputfile,'r');
if fid1 > 0
    fileInfo = dir(inputfile);
    fileSize = fileInfo.bytes;
    if fileSize > 10
        valence_grid_correct = 1;
    else
        valence_grid_correct = 0;
    end 
else
   valence_grid_correct = 0;
end
valence_grid_correct_flag = 0;
if (valence_grid_correct == 1)
  data = textscan(fid1, '%f',1,'headerlines',1);
  if sum(sum(abs(latticevectorfactor - cell2mat(data)))) > zero_tolerance
     valence_grid_correct_flag = 1;
  end
  clear data
  data = textscan(fid1, '%f %f %f',3);
  if sum(sum(abs(vectors - cell2mat(data)*latticevectorfactor*bohrperangstrom))) > zero_tolerance
    valence_grid_correct_flag = 1;
  end
  clear data;
  fclose(fid1);
  fid1 = fopen(inputfile,'r');
  for i = 1:nheader
     data = fgetl(fid1);
  end   
  if sum(sum(abs(atomspertype - str2num(data)))) > zero_tolerance
    valence_grid_correct_flag = 1;
  end
  clear data
  textscan(fid1, '%s', 1);
  data = textscan(fid1, '%f %f %f', natoms);
  if sum(sum(abs(direct_coords - cell2mat(data)))) > zero_tolerance
    valence_grid_correct_flag = 1;
  end
  clear data;
  data = textscan(fid1, '%f %f %f', 1);
  if sum(sum(abs(totnum - cell2mat(data)))) > zero_tolerance
    valence_grid_correct_flag = 1;
  end
  clear data;  
  data = textscan(fid1,'%f');
  raw_valence_pseudodensity = cell2mat(data);
  fclose(fid1);
  clear fid1;
  clear data;
end
if valence_grid_correct_flag == 1
    'The CHG file does not contain the same lattice vectors, grid spacing, or atom positions as the valence density (AECCAR2) file.'
    'This can occur if the files were produced from a geometry optimization rather than from a fixed geometry calculation.'
    'For best results, please regenerate the CHG and AECCAR2 files using a fixed geometry calculation. Charge analysis will proceed.'
    valence_grid_correct = 0;
end  
% construct the density array:
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
coords=zeros(natoms,3);
for k=1:natoms
    coords(k,1)=direct_coords(k,1)*vector1(1) + direct_coords(k,2)*vector2(1)+ direct_coords(k,3)*vector3(1);
    coords(k,2)=direct_coords(k,1)*vector1(2) + direct_coords(k,2)*vector2(2)+ direct_coords(k,3)*vector3(2);
    coords(k,3)=direct_coords(k,1)*vector1(3) + direct_coords(k,2)*vector2(3)+ direct_coords(k,3)*vector3(3);
end
origin = [0.0, 0.0, 0.0];
sum_negative_density = 0.0;
M = totnumA*totnumB*totnumC;
core_density = zeros(totnumA,totnumB,totnumC);
valence_density = zeros(totnumA,totnumB,totnumC);
if (valence_grid_correct == 1)
   valence_pseudodensity = zeros(totnumA,totnumB,totnumC);
end
for j = 1:M
    Q = j - 1;
    numA = mod(Q,totnumA);
    temp = (Q - numA)/totnumA;
    numB = mod(temp,totnumB);
    numC = (temp - numB)/totnumB;
    v=max(raw_valence(round(j)),0.0)/(pixelvolume*totnumA*totnumB*totnumC);
    if v*pixelvolume > pixel_integration_tolerance
        'Pixel with large density found. Program will correct density of this pixel.'
        correction_size = v*pixelvolume - pixel_integration_tolerance
        v = pixel_integration_tolerance/pixelvolume;
    end
    valence_density(numA+1,numB+1,numC+1) = v; % valence density
    if (valence_grid_correct == 1)
        valence_pseudodensity(numA+1,numB+1,numC+1) = raw_valence_pseudodensity(round(j))/(pixelvolume*totnumA*totnumB*totnumC);      
    end    
    if (core_available == 1)
      core_density(numA+1,numB+1,numC+1) = max(raw_core(round(j)),0.0)/(pixelvolume*totnumA*totnumB*totnumC);
    end
    if raw_valence(round(j)) < 0.0
        sum_negative_density = sum_negative_density - (pixelvolume*v - raw_valence(round(j))/(totnumA*totnumB*totnumC));
    end
end
clear raw_core raw_valence raw_valence_pseudodensity
sum_negative_density
% load the atomic numbers and numbers of valence and core electrons
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
core_electrons = zeros(natoms,1);
for i = 1:natomtypes
    for j = 1:atomspertype(i)
        count = count + 1;
        atomic_number(count,1) = atomic_number_type(1,i);
        core_electrons(count,1) = atomic_number(count,1) - nvalence_type(1,i);
    end
end
if core_available == 1
   missing_core=zeros(1,natoms);
   effective_nuclear_charge = atomic_number;
else
   missing_core = core_electrons;
end
% assign the number of core electrons and effective nuclear charge if core fitting is not performed
%% add the missing core densities
add_missing_core_density
%%%%%%%%%%%%%%
% open and read the charge and spin information from the CHG file
spin_flag = 0;
inputfile='CHG';
fid1 = fopen(inputfile,'r');
if fid1 > 0
    fileInfo = dir('CHG');
    fileSize = fileInfo.bytes;
    if fileSize > 10
        spin_available = 1;
    else
        spin_available = 0;
    end    
else
    fid2 = fopen('CHG_noncollinear','r');
    if (fid2 > 0)
        spin_available = 2;
    else
        spin_available = 0;
    end
end
while spin_available == 1
    data = textscan(fid1, '%f',1,'headerlines',1);
    if sum(sum(abs(latticevectorfactor - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data
    data = textscan(fid1, '%f %f %f',3);
    if sum(sum(abs(vectors - cell2mat(data)*latticevectorfactor*bohrperangstrom))) > zero_tolerance
        spin_flag = 1;
    end
    clear data;
    fclose(fid1);
    fid1 = fopen(inputfile,'r');
    for i = 1:nheader
       data = fgetl(fid1);
    end    
    if sum(sum(abs(atomspertype - str2num(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data
    textscan(fid1, '%s', 1);
    data = textscan(fid1, '%f %f %f', natoms);
    if sum(sum(abs(direct_coords - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data;
    data = textscan(fid1, '%f %f %f', 1);
    if sum(sum(abs(totnum - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data;
    if spin_flag == 1
        'The CHG and AECCAR0 files does not contain the same lattice vectors, grid spacing, or atom positions.'
        'This error can occur if the densities were produced from a geometry optimization rather than from a single-point calculation'
        'Please recalculate the input density files using a fixed geometry.'
        'Program will skip the computation of atomic spin magnetic moments.'
    end 
    fclose(fid1);
    fid1 = fopen(inputfile,'r');
    nheaderlines = nheader + 4 + natoms + ceil(totnumA*totnumB*totnumC/10);
    data = textscan(fid1,'%f','headerlines',nheaderlines);
    if sum(abs(length(cell2mat(data)) - totnumA*totnumB*totnumC)) > zero_tolerance
        % try this because sometimes a blank line occurs between the total and spin density blocks
        clear data
        fclose(fid1);
        fid1 = fopen(inputfile,'r');
        data = textscan(fid1,'%f','headerlines',(nheaderlines+1));      
    end 
    if length(cell2mat(data)) == totnumA*totnumB*totnumC
        raw_spin = cell2mat(data);
    else 
        spin_flag = 1;   
    end  
    fclose(fid1);
    clear fid1;
    if spin_flag == 1
        spin_available = 0;
        break
    end
    spin_density = zeros(totnumA,totnumB,totnumC);
    checkspin1_positive = 0.0;
    checkspin2_positive = 0.0;
    checkspin1_negative = 0.0;
    checkspin2_negative = 0.0;
    for j = 1:M
        Q = j - 1;
        numA = mod(Q,totnumA);
        temp = (Q - numA)/totnumA;
        numB = mod(temp,totnumB);
        numC = (temp - numB)/totnumB;
        p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1);
        if (p <= 0.0) || (raw_spin(round(j)) == 0.0) || (valence_pseudodensity(numA+1,numB+1,numC+1) <= 0.0)
            spin_density(numA+1,numB+1,numC+1) = 0.0;
        else
            s = raw_spin(round(j))/(pixelvolume*totnumA*totnumB*totnumC);
            s = s*min(abs(p/valence_pseudodensity(numA+1,numB+1,numC+1)),10.0);
            spin_density(numA+1,numB+1,numC+1) = s*min(abs(p)/abs(s),1.0);
        end
        if raw_spin(round(j)) > 0
           checkspin1_positive = checkspin1_positive + raw_spin(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_positive = checkspin2_positive + spin_density(numA+1,numB+1,numC+1)*pixelvolume;
        else
           checkspin1_negative = checkspin1_negative - raw_spin(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_negative = checkspin2_negative - spin_density(numA+1,numB+1,numC+1)*pixelvolume;   
        end
    end
    for i=1:3
       checkspin3_positive=0.0;
       checkspin3_negative=0.0;
       for ka=1:totnumA
           for kb=1:totnumB
               for kc=1:totnumC
                   p = core_density(ka,kb,kc) + valence_density(ka,kb,kc);
                   if (spin_density(ka,kb,kc) > 0) && (checkspin2_positive > zero_tolerance)
                      s = spin_density(ka,kb,kc)*checkspin1_positive/checkspin2_positive;
                      spin_density(ka,kb,kc) = s*min(abs(p)/abs(s),1.0);
                      checkspin3_positive = checkspin3_positive + spin_density(ka,kb,kc)*pixelvolume;
                   elseif (spin_density(ka,kb,kc) < 0) && (checkspin2_negative > zero_tolerance)
                      s = spin_density(ka,kb,kc)*checkspin1_negative/checkspin2_negative;
                      spin_density(ka,kb,kc) = s*min(abs(p)/abs(s),1.0);
                      checkspin3_negative = checkspin3_negative - spin_density(ka,kb,kc)*pixelvolume;
                   end
               end
           end
       end
       checkspin2_positive=checkspin3_positive;
       checkspin2_negative=checkspin3_negative;
    end
    checkspin3=checkspin3_positive - checkspin3_negative
    clear raw_spin
    break
end  
while spin_available == 2
    % read the 'CHG_noncollinear' file
    data = textscan(fid2, '%f',1,'headerlines',1);
    if sum(sum(abs(latticevectorfactor - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data
    data = textscan(fid2, '%f %f %f',3);
    if sum(sum(abs(vectors - cell2mat(data)*latticevectorfactor*bohrperangstrom))) > zero_tolerance
        spin_flag = 1;
    end
    fclose(fid2);
    clear data fid2;
    fid2 = fopen('CHG_noncollinear','r');
    for i = 1:nheader
        data = fgetl(fid2);
    end   
    if sum(sum(abs(atomspertype - str2num(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data
    textscan(fid2, '%s', 1);
    data = textscan(fid2, '%f %f %f', natoms);
    if sum(sum(abs(direct_coords - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data;
    data = textscan(fid2, '%f %f %f', 1);
    if sum(sum(abs(totnum - cell2mat(data)))) > zero_tolerance
        spin_flag = 1;
    end
    clear data;
    if spin_flag == 1
        'The CHG_noncollinear and AECCAR0 files does not contain the same lattice vectors, grid spacing, or atom positions.'
        'This error can occur if the densities were produced from a geometry optimization rather than from a single-point calculation'
        'Please recalculate the input density files using a fixed geometry.'
        'Program will skip the computation of atomic spin magnetic moments.'
    end
    % read the first spin component
    fclose(fid2);
    fid2 = fopen('CHG_noncollinear','r');
    nheaderlines = nheader + 4 + natoms + ceil(totnumA*totnumB*totnumC/10);
    data = textscan(fid2,'%f',1,'headerlines',nheaderlines);
    if abs(cell2mat(data) - totnumA) < zero_tolerance
        extra_lines=1;
    else
        extra_lines=0;
    end
    clear data
    nheaderlines = nheaderlines + extra_lines;
    fclose(fid2);
    fid2 = fopen('CHG_noncollinear','r');
    data = textscan(fid2,'%f',totnumA*totnumB*totnumC,'headerlines',nheaderlines);
    if length(cell2mat(data)) == totnumA*totnumB*totnumC
        raw_spin_1 = cell2mat(data);
    else 
        spin_flag = 1;   
    end  
    fclose(fid2);
    clear fid2;
    if spin_flag == 1
        spin_available = 0;
        break
    end
    % read the second spin component
    fid2 = fopen('CHG_noncollinear','r');
    nheaderlines = nheaderlines + 1 + ceil(totnumA*totnumB*totnumC/10);
    data = textscan(fid2,'%f',1,'headerlines',nheaderlines);
    if abs(cell2mat(data) - totnumA) < zero_tolerance
        extra_lines=1;
    else
        extra_lines=0;
    end
    clear data
    nheaderlines = nheaderlines + extra_lines;
    fclose(fid2);
    fid2 = fopen('CHG_noncollinear','r');
    data = textscan(fid2,'%f',totnumA*totnumB*totnumC,'headerlines',nheaderlines);
    if length(cell2mat(data)) == totnumA*totnumB*totnumC
        raw_spin_2 = cell2mat(data);
    else 
        spin_flag = 1;   
    end  
    fclose(fid2);
    clear fid2;
    if spin_flag == 1
        spin_available = 0;
        break
    end 
    % read the third spin component
    fid2 = fopen('CHG_noncollinear','r');
    nheaderlines = nheaderlines + 1 + ceil(totnumA*totnumB*totnumC/10);
    data = textscan(fid2,'%f',1,'headerlines',nheaderlines);
    if abs(cell2mat(data) - totnumA) < zero_tolerance
        extra_lines=1;
    else
        extra_lines=0;
    end
    clear data
    nheaderlines = nheaderlines + extra_lines;
    fclose(fid2);
    fid2 = fopen('CHG_noncollinear','r');
    data = textscan(fid2,'%f','headerlines',nheaderlines);
    if length(cell2mat(data)) == totnumA*totnumB*totnumC
        raw_spin_3 = cell2mat(data);
    else 
        spin_flag = 1;   
    end  
    fclose(fid2);
    clear fid2;
    if spin_flag == 1
        spin_available = 0;
        break
    end 
    % compute the x, y, and z spin components using the SAXIS value
    alpha = atan2(SAXIS(2),SAXIS(1));
    beta = atan2(sqrt(SAXIS(1)^2 + SAXIS(2)^2),SAXIS(3));
    if (abs(alpha) > 1.58) || (abs(beta) > 1.58)
        'Error reading SAXIS value. Program will terminate.'
        flag = 1;
        break
    end
    raw_spin_x = cos(beta)*cos(alpha)*raw_spin_1 - sin(alpha)*raw_spin_2 + sin(beta)*cos(alpha)*raw_spin_3;
    raw_spin_y = cos(beta)*sin(alpha)*raw_spin_1 + cos(alpha)*raw_spin_2 + sin(beta)*sin(alpha)*raw_spin_3;
    raw_spin_z = -sin(beta)*raw_spin_1 + cos(beta)*raw_spin_3;
    clear raw_spin_1 raw_spin_2 raw_spin_3
    % construct the spin density matrices
    spin_density_vector = zeros(totnumA,totnumB,totnumC,3);
    checkspin1_positive = zeros(1,3);
    checkspin2_positive = zeros(1,3);
    checkspin1_negative = zeros(1,3);
    checkspin2_negative = zeros(1,3);
    for j = 1:M
        Q = j - 1;
        numA = mod(Q,totnumA);
        temp = (Q - numA)/totnumA;
        numB = mod(temp,totnumB);
        numC = (temp - numB)/totnumB;
        p = core_density(numA+1,numB+1,numC+1) + valence_density(numA+1,numB+1,numC+1);
        s = sqrt(raw_spin_x(round(j))^2 + raw_spin_y(round(j))^2 + raw_spin_z(round(j))^2)/(pixelvolume*totnumA*totnumB*totnumC);
        if (p <= 0.0) || (s == 0.0) || (valence_pseudodensity(numA+1,numB+1,numC+1) <= 0.0)
            spin_density_vector(numA+1,numB+1,numC+1,1) = 0.0;
            spin_density_vector(numA+1,numB+1,numC+1,2) = 0.0;
            spin_density_vector(numA+1,numB+1,numC+1,3) = 0.0;
        else
            enhancement_factor = min(abs(p/valence_pseudodensity(numA+1,numB+1,numC+1)),10.0);
            spin_density_vector(numA+1,numB+1,numC+1,1) = enhancement_factor*raw_spin_x(round(j))*min(p/(enhancement_factor*s),1.0)/(pixelvolume*totnumA*totnumB*totnumC);
            spin_density_vector(numA+1,numB+1,numC+1,2) = enhancement_factor*raw_spin_y(round(j))*min(p/(enhancement_factor*s),1.0)/(pixelvolume*totnumA*totnumB*totnumC);
            spin_density_vector(numA+1,numB+1,numC+1,3) = enhancement_factor*raw_spin_z(round(j))*min(p/(enhancement_factor*s),1.0)/(pixelvolume*totnumA*totnumB*totnumC);
        end
        if raw_spin_x(round(j)) > 0
           checkspin1_positive(1) = checkspin1_positive(1) + raw_spin_x(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_positive(1) = checkspin2_positive(1) + spin_density_vector(numA+1,numB+1,numC+1,1)*pixelvolume;
        else
           checkspin1_negative(1) = checkspin1_negative(1) - raw_spin_x(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_negative(1) = checkspin2_negative(1) - spin_density_vector(numA+1,numB+1,numC+1,1)*pixelvolume;   
        end
        if raw_spin_y(round(j)) > 0
           checkspin1_positive(2) = checkspin1_positive(2) + raw_spin_y(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_positive(2) = checkspin2_positive(2) + spin_density_vector(numA+1,numB+1,numC+1,2)*pixelvolume;
        else
           checkspin1_negative(2) = checkspin1_negative(2) - raw_spin_y(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_negative(2) = checkspin2_negative(2) - spin_density_vector(numA+1,numB+1,numC+1,2)*pixelvolume;   
        end
        if raw_spin_z(round(j)) > 0
           checkspin1_positive(3) = checkspin1_positive(3) + raw_spin_z(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_positive(3) = checkspin2_positive(3) + spin_density_vector(numA+1,numB+1,numC+1,3)*pixelvolume;
        else
           checkspin1_negative(3) = checkspin1_negative(3) - raw_spin_z(round(j))/(totnumA*totnumB*totnumC);
           checkspin2_negative(3) = checkspin2_negative(3) - spin_density_vector(numA+1,numB+1,numC+1,3)*pixelvolume;   
        end
    end
    for i=1:3
       checkspin3_positive=zeros(1,3);
       checkspin3_negative=zeros(1,3);
       for ka=1:totnumA
           for kb=1:totnumB
               for kc=1:totnumC
                  p = core_density(ka,kb,kc) + valence_density(ka,kb,kc);
                  s = sqrt(spin_density_vector(ka,kb,kc,1)^2 + spin_density_vector(ka,kb,kc,2)^2 + spin_density_vector(ka,kb,kc,3)^2);
                  s_trial(1)=spin_density_vector(ka,kb,kc,1);
                  s_trial(2)=spin_density_vector(ka,kb,kc,2);
                  s_trial(3)=spin_density_vector(ka,kb,kc,3);
                  for component=1:3
                      if (spin_density_vector(ka,kb,kc,component) > 0) && (checkspin2_positive(component) > zero_tolerance)
                         s_trial(component) = spin_density_vector(ka,kb,kc,component)*checkspin1_positive(component)/checkspin2_positive(component);
                      elseif (spin_density_vector(ka,kb,kc,component) < 0) && (checkspin2_negative(component) > zero_tolerance)
                         s_trial(component) = spin_density_vector(ka,kb,kc,component)*checkspin1_negative(component)/checkspin2_negative(component);
                      end
                  end
                  mag_s_trial = sqrt(s_trial(1)^2 + s_trial(2)^2 + s_trial(3)^2);
                  if mag_s_trial > 0
                     spin_density_vector(ka,kb,kc,1) = s_trial(1)*min(p/mag_s_trial,1.0);
                     spin_density_vector(ka,kb,kc,2) = s_trial(2)*min(p/mag_s_trial,1.0);
                     spin_density_vector(ka,kb,kc,3) = s_trial(3)*min(p/mag_s_trial,1.0);
                  end
                  for component=1:3
                      if spin_density_vector(ka,kb,kc,component) > 0 
                         checkspin3_positive(component) = checkspin3_positive(component) + spin_density_vector(ka,kb,kc,component)*pixelvolume;
                      else
                         checkspin3_negative(component) = checkspin3_negative(component) - spin_density_vector(ka,kb,kc,component)*pixelvolume;
                      end
                  end
               end
           end
       end
       checkspin2_positive=checkspin3_positive;
       checkspin2_negative=checkspin3_negative;
    end
    checkspin3=checkspin3_positive - checkspin3_negative
    checkspin3_mag = sqrt(checkspin3(1)^2 + checkspin3(2)^2 + checkspin3(3)^2)
    clear raw_spin_x raw_spin_y raw_spin_z
    break
end
if flag == 1
    break
end
initialize_atomic_densities
if flag == 1
    break
end
% compute the valence occupancy corrections
occupancy_correction = zeros(natoms,11);
if (valence_grid_correct == 1)
    dominant_atom_weight = neutral_density;
    compute_dominant_atom_volumes
    for ka = 1:totnumA
        for kb = 1:totnumB
            for kc = 1:totnumC
                if dominant_atom_points(ka,kb,kc) > 0
                   occupancy_correction(round(dominant_atom_points(ka,kb,kc)),1) = occupancy_correction(round(dominant_atom_points(ka,kb,kc)),1) + pixelvolume*(valence_pseudodensity(ka,kb,kc) - valence_density(ka,kb,kc));
                end   
            end
        end
    end
    'The largest occupancy correction was'
    max(max(abs(occupancy_correction)))
end    


