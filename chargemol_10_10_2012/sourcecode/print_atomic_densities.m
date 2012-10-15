cutoff_radius_string=num2str(cutoff_radius,'%03i');
nshells_string=num2str(nshells,'%03i'); 
for j=1:natoms
    atom_no_string=num2str(j,'%03i');
    combinedstring=strcat('atom_',atom_no_string,'_partial_density_',cutoff_radius_string,'_',nshells_string,'.txt');
    fid = fopen(combinedstring,'w');
    for k = 1:nshells
       temp = fprintf(fid,'%14.8E\n',partial_density(j,k)); 
    end
    fclose(fid);
    clear fid;
end