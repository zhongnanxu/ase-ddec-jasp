valence_population=zeros(1,natoms);
core_population=zeros(1,natoms);
spherical_average_density=zeros(natoms,nshells);
spherical_avg_ref_pseudodensity=zeros(natoms,nshells);
avg_PtoWref = ones(natoms,nshells);
change=zeros(1,natoms);
max_change = 1.0;
old_change = 1.0;
old_old_change = 1.0;
neutral_density = zeros(natoms,nshells);
flag = 0;
test = 0;
for j = 1:natoms
    % construct the file name to read
    string0='_';
    string1=num2str(atomic_number(j),'%03i');   
    string4=num2str(cutoff_radius,'%03i');
    string5=num2str(nshells,'%03i');
    string6='.txt';
    combinedstring=strcat(atomic_densities_directory,density_set_prefix,string0,string1,string0,string1,string0,string1,string0,string4,string0,string5,string6);
    fid = fopen(combinedstring,'r');
    test = fid;
    if test == -1
        clear fid
        string0='_';
        string1=num2str(atomic_number(j));   
        string4=num2str(cutoff_radius);
        string5=num2str(nshells);
        string6='.txt';
        combinedstring=strcat(atomic_densities_directory,string1,string0,string1,string0,string1,string0,string4,string0,string5,string6);
        fid = fopen(combinedstring,'r'); 
        if fid == -1
            combinedstring
            'Could not find a suitable reference density. Program will terminate.'
            flag = 1;
            break
        end
    end
    % read in the density
    data = textscan(fid, '%f',nshells,'headerlines',12);
    temp = cell2mat(data);
    for k = 1:nshells
        neutral_density(j,k) = temp(k) + 1.0e-16;
    end
    fclose(fid);
    clear data temp fid;
end
if flag == 1
    break
end
% initialize some atomic density distributions
combined_oxidation_density = neutral_density;
partial_density = neutral_density;
normalized_oxidation_density = neutral_density;
partial_core_density = neutral_density;