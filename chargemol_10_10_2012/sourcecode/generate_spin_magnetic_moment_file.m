V1 = vector1*periodicA/bohrperangstrom;
V2 = vector2*periodicB/bohrperangstrom;
V3 = vector3*periodicC/bohrperangstrom;
filename = 'DDEC_atomic_spin_moments.xyz';
fid = fopen(filename,'wt');
count=fprintf(fid, '%3.0f\n', natoms);
if periodicA || periodicB || periodicC
    count=fprintf(fid, 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell [{ %f %f %f }, { %f %f %f }, { %f %f %f }]\n',V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),V3(1),V3(2),V3(3));    
else    
    count=fprintf(fid, '%s\n', 'nonperiodic system');
end

%
for i = 1:natoms
   z = final_result(i,2);
   atomic_number_to_symbol
   count=fprintf(fid, '%s', atomic_symbol);
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,3)/bohrperangstrom);
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,4)/bohrperangstrom);
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,5)/bohrperangstrom);
   count=fprintf(fid, '%s', '   ');
   if spin_available == 2
       count=fprintf(fid, '%f', sqrt(spin_population_vector(i,1)^2 + spin_population_vector(i,2)^2 + spin_population_vector(i,3)^2));
       count=fprintf(fid, '%s', '   ');
       count=fprintf(fid, '%f', spin_population_vector(i,1));
       count=fprintf(fid, '%s', '   ');
       count=fprintf(fid, '%f', spin_population_vector(i,2));
       count=fprintf(fid, '%s', '   ');
       count=fprintf(fid, '%f\n', spin_population_vector(i,3));
   elseif spin_available == 1
       count=fprintf(fid, '%f\n', spin_population(1,i));   
   end    
end
% Print information at the bottom of the file
count=fprintf(fid, '%s\n', '   ');
if spin_available == 2
    count=fprintf(fid, '%s\n', 'Noncollinear spin population analysis was performed ');
    count=fprintf(fid, '%s', 'The total spin magnetic moment vector of the unit cell is   ');
    count=fprintf(fid, '%f', tot_spin_moment_vector(1));
    count=fprintf(fid, '%s', '   ');
    count=fprintf(fid, '%f', tot_spin_moment_vector(2));
    count=fprintf(fid, '%s', '   ');
    count=fprintf(fid, '%f\n', tot_spin_moment_vector(3));    
elseif spin_available == 1 
    count=fprintf(fid, '%s\n', 'Collinear spin population analysis was performed ');
    count=fprintf(fid, '%s', 'The total spin magnetic moment of the unit cell is   ');
    count=fprintf(fid, '%f\n', tot_spin_moment);
end
% Print information at the bottom of the file
count=fprintf(fid, '%s\n', '   ');
%fprintf (fid, '%s\n', version );
fprintf (fid, '%s\n', 'See ddec.sourceforge.net for latest version.' );
count=fprintf(fid, '%s\n', '   ');
fprintf (fid, '%s\n', 'Computational parameters:' );
count=fprintf(fid, '%s', 'spin_ref_fraction = ');
count=fprintf(fid, '%3.2f\n', spin_ref_fraction);
count=fprintf(fid, '%s', 'Number of radial integration shells = ');
count=fprintf(fid, '%i\n', nshells);
count=fprintf(fid, '%s', 'Cutoff radius (pm) = ');
count=fprintf(fid, '%i\n', cutoff_radius);
count=fprintf(fid, '%s', 'spin convergence tolerance = ');
count=fprintf(fid, '%f\n', spin_convergence_tolerance);
count=fprintf(fid, '%s', 'Number of iterations to convergence = ');
count=fprintf(fid, '%i\n', iter);
count=fprintf(fid, '%s\n', '   ');
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf (fid, '%s\n', s );
fclose(fid);
clear fid;


