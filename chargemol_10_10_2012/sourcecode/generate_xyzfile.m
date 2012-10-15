t = now;
c = datevec ( t );
V1 = vector1*periodicA/bohrperangstrom;
V2 = vector2*periodicB/bohrperangstrom;
V3 = vector3*periodicC/bohrperangstrom;
filename = 'DDEC_net_atomic_charges.xyz'; 
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
   count=fprintf(fid, '%f\n', final_result(i,6));
end
% Print information at the bottom of the file
count=fprintf(fid, '%s\n', '   ');
%fprintf (fid, '%s\n', version );
fprintf (fid, '%s\n', 'See ddec.sourceforge.net for latest version.' );
count=fprintf(fid, '%s\n', '   ');
fprintf (fid, '%s\n', 'Computational parameters:' );
count=fprintf(fid, '%s', 'reference_weighting = ');
count=fprintf(fid, '%3.2f\n', reference_weighting);
count=fprintf(fid, '%s', 'Number of radial integration shells = ');
count=fprintf(fid, '%i\n', nshells);
count=fprintf(fid, '%s', 'Cutoff radius (pm) = ');
count=fprintf(fid, '%i\n', cutoff_radius);
count=fprintf(fid, '%s', 'Error in the integrated total number of electrons before renormalization (e) = ');
count=fprintf(fid, '%i\n', checkme);
count=fprintf(fid, '%s', 'charge convergence tolerance = ');
count=fprintf(fid, '%f\n', charge_convergence_tolerance);
count=fprintf(fid, '%s', 'minimum radius for electron cloud penetration fitting (pm) = ');
count=fprintf(fid, '%f\n', rmin_cloud_penetration);
count=fprintf(fid, '%s', 'Decay exponent for electron density of buried atom tails = ');
count=fprintf(fid, '%f\n', density_decaying_exponent);
count=fprintf(fid, '%s', 'Number of iterations to convergence = ');
count=fprintf(fid, '%i\n', iter);
count=fprintf(fid, '%s\n', '   ');
count=fprintf(fid, '%s\n', 'The following XYZ coordinates are in angstroms. The atomic dipoles and quadrupoles are in atomic units.');
count=fprintf(fid, '%s\n', 'atomic symbol, x, y, z, net_charge, dipole_x, dipole_y, dipole_z, dipole_mag, Qxy, Qxz, Qyz, Q(x^2-y^2), Q(3z^2 - R^2)');
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
   count=fprintf(fid, '%f', final_result(i,6));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,7));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,8));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,9));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,10));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,11));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,12));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,13));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', final_result(i,14));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f\n', final_result(i,15));
end
count=fprintf(fid, '%s\n', '   ');
count=fprintf(fid, '%s\n', 'The sperically averaged electron density of each atom fit to a function of the form exp(a - br) for r >= rmin_cloud_penetration');
count=fprintf(fid, '%s\n', 'atomic symbol, x, y, z, a, b, Rsquared where a and b are in atomic units and Rsquared is the squared correlation coefficient');
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
   count=fprintf(fid, '%f', intercept(i));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f', -slope(i));
   count=fprintf(fid, '%s', '   ');
   count=fprintf(fid, '%f\n', Rsquared(i));
end
count=fprintf(fid, '%s\n', '   ');
s = datestr ( c, 0 );
fprintf (fid, '%s\n', s );
fclose(fid);
clear fid;


