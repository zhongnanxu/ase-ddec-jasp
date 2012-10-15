% Generates valence, spin, and core densities from natural orbital
% coefficients with gaussian primitive basis set
%
% Step 1: Choose the Grid Points
%
% Choose the pixel boundary used to generate the grid
boundary=zeros(3,3);
temp_vector1=zeros(1,3);
temp_vector2=zeros(1,3);
temp_vector3=zeros(1,3);
if num_periodic_directions == 0
    boundary = [preferred_grid_spacing, 0.0, 0.0; 0.0, preferred_grid_spacing, 0.0; 0.0, 0.0, preferred_grid_spacing];
else
    vector_length=zeros(num_periodic_directions,1);
    for i=1:num_periodic_directions
        vector_length(i)= sqrt(periodic_vectors(i,1)^2 + periodic_vectors(i,2)^2 + periodic_vectors(i,3)^2);
        boundary(i,:) = periodic_vectors(i,:)/ceil(vector_length(i)/preferred_grid_spacing);        
    end
    if num_periodic_directions == 1
       if dot(periodic_vectors,[ 1, 0, 0]) > dot(periodic_vectors,[ 0, 1, 0])    
           temp_vector2=cross(periodic_vectors,[ 0, 1, 0]);
       else
           temp_vector2=cross(periodic_vectors,[ 1, 0, 0]);
       end
       boundary(2,:) = temp_vector2(:)*preferred_grid_spacing/norm(temp_vector2);
    end       
    if num_periodic_directions < 3
        temp_vector1(:)=boundary(1,:);
        temp_vector2(:)=boundary(2,:);
        temp_vector3(:)=cross(temp_vector1,temp_vector2);
        boundary(3,:)= temp_vector3(:)*preferred_grid_spacing/norm(temp_vector3);
    end
end
periodicA=0;
periodicB=0;
periodicC=0;
if num_periodic_directions > 0.5
    periodicA=1;
end
if num_periodic_directions > 1.5
    periodicB=1;
end
if num_periodic_directions > 2.5
    periodicC=1;
end
temp_vector1(:)=boundary(1,:);
temp_vector2(:)=boundary(2,:);
temp_vector3(:)=boundary(3,:);
pixelvolume=det(boundary);
'Initial setup of atoms using temporary origin (0,0,0)'
origin = [0.0, 0.0, 0.0];
charge_center_positions
min_na=min(center_nabc(:,1));
min_nb=min(center_nabc(:,2));
min_nc=min(center_nabc(:,3));
max_na=max(center_nabc(:,1));
max_nb=max(center_nabc(:,2));
max_nc=max(center_nabc(:,3));
parallelpiped
'Compute a better origin and number of grid points along each grid direction'
if periodicA == 0
    totnumA = max_na - min_na + 2*delta_na
    origin = origin + (min_na - delta_na)*temp_vector1;
else
    totnumA = ceil(vector_length(1)/preferred_grid_spacing)
end
if periodicB == 0
    totnumB = max_nb - min_nb + 2*delta_nb
    origin = origin + (min_nb - delta_nb)*temp_vector2;
else
    totnumB = ceil(vector_length(2)/preferred_grid_spacing)
end
if periodicC == 0
    totnumC = max_nc - min_nc + 2*delta_nc
    origin = origin + (min_nc - delta_nc)*temp_vector3;
else
    totnumC = ceil(vector_length(3)/preferred_grid_spacing)
end
'The new origin is'
origin
vector1=temp_vector1*totnumA;
vector2=temp_vector2*totnumB;
vector3=temp_vector3*totnumC;
'Compute atomic positions using the better origin'
charge_center_positions
% The
if num_periodic_directions > 0.5
    periodicsumlimitA = ceil(periodic_cutoff_length/sqrt(vector1(1)^2 + vector1(2)^2 + vector1(3)^2))
else
    periodicsumlimitA = 0;
end
if num_periodic_directions > 1.5
    periodicsumlimitB = ceil(periodic_cutoff_length/sqrt(vector2(1)^2 + vector2(2)^2 + vector2(3)^2))
else
    periodicsumlimitB = 0;
end
if num_periodic_directions > 2.5
    periodicsumlimitC = ceil(periodic_cutoff_length/sqrt(vector3(1)^2 + vector3(2)^2 + vector3(3)^2))
else
    periodicsumlimitC = 0;
end
nimages=(2*periodicsumlimitA + 1)*(2*periodicsumlimitB + 1)*(2*periodicsumlimitC + 1)
% Step 2: Compute the density matrices for alpha and beta electrons
alpha_1PDM = zeros(nprimitives,nprimitives);
beta_1PDM = zeros(nprimitives,nprimitives);
for j=1:nprimitives
    for k = j:nprimitives
        for i=1:norbitals
            alpha_1PDM(j,k) = alpha_1PDM(j,k) + alpha_beta_occupation(i,1)*orbital_coefficients(i,j)*orbital_coefficients(i,k);
            beta_1PDM(j,k) = beta_1PDM(j,k) + alpha_beta_occupation(i,2)*orbital_coefficients(i,j)*orbital_coefficients(i,k);
        end
    end    
end 
'Setup of grid points was successful.'
% Step 3: Use cutoff condition to determine which primitive pairs (i,j) to include
included_primitive_pairs=zeros((nprimitives*(nprimitives+1)*nimages/2),1);
analytic_pair_count=0;
numeric_pair_count=0;
total_count=0;
for i = 1:nprimitives
    for repeata=-periodicsumlimitA:periodicsumlimitA
        for repeatb=-periodicsumlimitB:periodicsumlimitB
            for repeatc=-periodicsumlimitC:periodicsumlimitC
                for j = i:nprimitives
                    total_count = total_count + 1;
                    alpha1 = primitive_exponents(i,2);
                    alpha2 = primitive_exponents(j,2);
                    X1=basis_set_centers(round(primitive_exponents(i,1)),3);
                    Y1=basis_set_centers(round(primitive_exponents(i,1)),4);
                    Z1=basis_set_centers(round(primitive_exponents(i,1)),5);
                    X2=basis_set_centers(round(primitive_exponents(j,1)),3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1);
                    Y2=basis_set_centers(round(primitive_exponents(j,1)),4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2);
                    Z2=basis_set_centers(round(primitive_exponents(j,1)),5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3);   
                    % Round powers to even number so overlap test is positive definite
                    Lx1=primitive_exponents(i,3) - mod(primitive_exponents(i,3),2);
                    Ly1=primitive_exponents(i,4) - mod(primitive_exponents(i,4),2);
                    Lz1=primitive_exponents(i,5) - mod(primitive_exponents(i,5),2);
                    Lx2=primitive_exponents(j,3) - mod(primitive_exponents(j,3),2);
                    Ly2=primitive_exponents(j,4) - mod(primitive_exponents(j,4),2);
                    Lz2=primitive_exponents(j,5) - mod(primitive_exponents(j,5),2);
                    max_coeff = abs(alpha_1PDM(i,j)) + abs(beta_1PDM(i,j));
                    test = max_coeff*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2);
                    if test < gaussian_overlap_tolerance
                        continue
                    else 
                        if (alpha1*min(basis_set_centers(primitive_exponents(i,1),2),1) + alpha2*min(basis_set_centers(primitive_exponents(j,1),2),1)) > analytic_alpha_cutoff
                            analytic_pair_count = analytic_pair_count + 1;
                            included_primitive_pairs(total_count)= 1.0;                
                        else
                            numeric_pair_count = numeric_pair_count + 1;
                            included_primitive_pairs(total_count)= 2.0;
                        end     
                    end
                end
            end
        end   
    end
end
'The number of primitive products to treat analytically is'
num_analytic_primitive_pairs = analytic_pair_count
'The number of primitive products to treat numerically is'
num_numeric_primitive_pairs = numeric_pair_count
'The total number of included primitives is'
num_included_primitive_pairs= num_analytic_primitive_pairs + num_numeric_primitive_pairs
% Step 4: Check whether the primitive pairs give the correct total number of electons
check_included_electrons=0;
pair_overlap=zeros((nprimitives*(nprimitives+1)*nimages/2),1);
total_count = 0;
for i = 1:nprimitives
    for repeata=-periodicsumlimitA:periodicsumlimitA
        for repeatb=-periodicsumlimitB:periodicsumlimitB
            for repeatc=-periodicsumlimitC:periodicsumlimitC
                for j = i:nprimitives
                    total_count = total_count + 1;
                    if included_primitive_pairs(total_count)== 0
                        continue
                    end    
                    alpha1=primitive_exponents(i,2);
                    alpha2=primitive_exponents(j,2);
                    X1=basis_set_centers(round(primitive_exponents(i,1)),3);
                    Y1=basis_set_centers(round(primitive_exponents(i,1)),4);
                    Z1=basis_set_centers(round(primitive_exponents(i,1)),5);
                    X2=basis_set_centers(round(primitive_exponents(j,1)),3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1);
                    Y2=basis_set_centers(round(primitive_exponents(j,1)),4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2);
                    Z2=basis_set_centers(round(primitive_exponents(j,1)),5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3);          
                    Lx1 = primitive_exponents(i,3);
                    Ly1 = primitive_exponents(i,4);
                    Lz1 = primitive_exponents(i,5);
                    Lx2 = primitive_exponents(j,3);
                    Ly2 = primitive_exponents(j,4);
                    Lz2 = primitive_exponents(j,5); 
                    pair_overlap(total_count) = gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2);
                    if j > i 
                        check_included_electrons = check_included_electrons + 2*(alpha_1PDM(i,j) + beta_1PDM(i,j))*pair_overlap(total_count);
                    else
                        check_included_electrons = check_included_electrons + (alpha_1PDM(i,j) + beta_1PDM(i,j))*pair_overlap(total_count);
                    end
                end
            end    
        end
    end
end
total_count=0;
if abs(check_included_electrons - included_electrons) > 0.001
   'The one-particle density matrix is not properly normalized.'
    check_included_electrons - included_electrons
    'Computing the self-overlap (should equal one) for each natural orbital.'
    MO_normalization_check=zeros(norbitals,1);
    total_count = 0;
    for i = 1:nprimitives
        for repeata=-periodicsumlimitA:periodicsumlimitA
            for repeatb=-periodicsumlimitB:periodicsumlimitB
                for repeatc=-periodicsumlimitC:periodicsumlimitC
                    for j = i:nprimitives
                        total_count = total_count + 1;
                        if included_primitive_pairs(total_count)== 0
                            continue
                        end 
                        for k=1:norbitals
                            MO_normalization_check(k) = MO_normalization_check(k) + orbital_coefficients(k,i)*orbital_coefficients(k,j)*pair_overlap(total_count);
                        end
                    end
                end
            end    
        end    
    end
    MO_normalization_check 
    'Program will terminate.'  
    flag=1
    break
else
    'The one-particle density matrix is properly normalized.'
end
if flag == 1
    break
end    
% Collect key data for each primitive pair
included_primitive_data=zeros(num_numeric_primitive_pairs,33);
count = 0;
analytic_pair_count=0;
total_count = 0;
temp_nabc=zeros(1,3);
temp_center_shift=zeros(1,3);
occupancy_correction = zeros(natoms,11);
for i = 1:nprimitives
    for repeata=-periodicsumlimitA:periodicsumlimitA
        for repeatb=-periodicsumlimitB:periodicsumlimitB
            for repeatc=-periodicsumlimitC:periodicsumlimitC
                for j = i:nprimitives
                    total_count = total_count + 1;
                    if included_primitive_pairs(total_count)== 0
                        continue
                    elseif included_primitive_pairs(total_count) == 1 % analytic pairs
                        analytic_pair_count = analytic_pair_count + 1;
                        alpha1=primitive_exponents(i,2);
                        alpha2=primitive_exponents(j,2);
                        X1=basis_set_centers(round(primitive_exponents(i,1)),3);
                        Y1=basis_set_centers(round(primitive_exponents(i,1)),4);
                        Z1=basis_set_centers(round(primitive_exponents(i,1)),5);
                        X2=basis_set_centers(round(primitive_exponents(j,1)),3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1);
                        Y2=basis_set_centers(round(primitive_exponents(j,1)),4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2);
                        Z2=basis_set_centers(round(primitive_exponents(j,1)),5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3);        
                        Lx1 = primitive_exponents(i,3);
                        Ly1 = primitive_exponents(i,4);
                        Lz1 = primitive_exponents(i,5);
                        Lx2 = primitive_exponents(j,3);
                        Ly2 = primitive_exponents(j,4);
                        Lz2 = primitive_exponents(j,5); 
                        if alpha1*min(basis_set_centers(primitive_exponents(i,1),2),1) > alpha2*min(basis_set_centers(primitive_exponents(j,1),2),1)
                            partial_pair_weight_i = 1;
                            partial_pair_weight_j = 0;
                        elseif alpha1*min(basis_set_centers(primitive_exponents(i,1),2),1) < alpha2*min(basis_set_centers(primitive_exponents(j,1),2),1)
                            partial_pair_weight_i = 0;
                            partial_pair_weight_j = 1;
                        else
                            partial_pair_weight_i = 0.5;
                            partial_pair_weight_j = 0.5;                
                        end   
                        if j > i 
                            % factor of 2 to account for symmetry of 1-PDM and index summation over i <= j
                            partial_pair_weight_i = partial_pair_weight_i*2;
                            partial_pair_weight_j = partial_pair_weight_j*2;
                        end
                        full_1PDM = alpha_1PDM(i,j) + beta_1PDM(i,j);
                        atom_index_for_i = basis_set_centers(round(primitive_exponents(i,1)),6);
                        atom_index_for_j = basis_set_centers(round(primitive_exponents(j,1)),6);
                        if atom_index_for_i > 0.5
                            % number of valence electron correction
                            occupancy_correction(atom_index_for_i,1) = occupancy_correction(atom_index_for_i,1) + full_1PDM*pair_overlap(total_count)*partial_pair_weight_i;
                            % spin magnetic moment correction
                            occupancy_correction(atom_index_for_i,2) = occupancy_correction(atom_index_for_i,2) + (alpha_1PDM(i,j) - beta_1PDM(i,j))*pair_overlap(total_count)*partial_pair_weight_i;
                            % dipole corrections
                            occupancy_correction(atom_index_for_i,3) = occupancy_correction(atom_index_for_i,3) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2);
                            occupancy_correction(atom_index_for_i,4) = occupancy_correction(atom_index_for_i,4) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1+1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2);
                            occupancy_correction(atom_index_for_i,5) = occupancy_correction(atom_index_for_i,5) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1+1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2);
                            % quadrupole corrections
                            occupancy_correction(atom_index_for_i,6) = occupancy_correction(atom_index_for_i,6) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+2,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2); % x^2
                            occupancy_correction(atom_index_for_i,7) = occupancy_correction(atom_index_for_i,7) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1+2,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2); % y^2
                            occupancy_correction(atom_index_for_i,8) = occupancy_correction(atom_index_for_i,8) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1+2,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2); % z^2
                            occupancy_correction(atom_index_for_i,9) = occupancy_correction(atom_index_for_i,9) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+1,Ly1+1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2); % xy
                            occupancy_correction(atom_index_for_i,10) = occupancy_correction(atom_index_for_i,10) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1+1,Ly1,Lz1+1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2); % xz
                            occupancy_correction(atom_index_for_i,11) = occupancy_correction(atom_index_for_i,11) - full_1PDM*partial_pair_weight_i*gaussian_overlap_integral(alpha1,Lx1,Ly1+1,Lz1+1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2); % yz
                        end
                        if atom_index_for_j > 0.5
                            % number of valence electron correction
                            occupancy_correction(atom_index_for_j,1) = occupancy_correction(atom_index_for_j,1) + full_1PDM*pair_overlap(total_count)*partial_pair_weight_j;
                            % spin magnetic moment correction
                            occupancy_correction(atom_index_for_j,2) = occupancy_correction(atom_index_for_j,2) + (alpha_1PDM(i,j) - beta_1PDM(i,j))*pair_overlap(total_count)*partial_pair_weight_j;
                            % dipole corrections
                            occupancy_correction(atom_index_for_j,3) = occupancy_correction(atom_index_for_j,3) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+1,Ly2,Lz2,X2,Y2,Z2);
                            occupancy_correction(atom_index_for_j,4) = occupancy_correction(atom_index_for_j,4) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2+1,Lz2,X2,Y2,Z2);
                            occupancy_correction(atom_index_for_j,5) = occupancy_correction(atom_index_for_j,5) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2+1,X2,Y2,Z2);
                            % quadrupole corrections
                            occupancy_correction(atom_index_for_j,6) = occupancy_correction(atom_index_for_j,6) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+2,Ly2,Lz2,X2,Y2,Z2); % x^2
                            occupancy_correction(atom_index_for_j,7) = occupancy_correction(atom_index_for_j,7) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2+2,Lz2,X2,Y2,Z2); % y^2
                            occupancy_correction(atom_index_for_j,8) = occupancy_correction(atom_index_for_j,8) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2+2,X2,Y2,Z2); % z^2
                            occupancy_correction(atom_index_for_j,9) = occupancy_correction(atom_index_for_j,9) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+1,Ly2+1,Lz2,X2,Y2,Z2); % xy
                            occupancy_correction(atom_index_for_j,10) = occupancy_correction(atom_index_for_j,10) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2+1,Ly2,Lz2+1,X2,Y2,Z2); % xz
                            occupancy_correction(atom_index_for_j,11) = occupancy_correction(atom_index_for_j,11) - full_1PDM*partial_pair_weight_j*gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2+1,Lz2+1,X2,Y2,Z2); % yz 
                        end                               
                    end
                    count = count + 1;
                    included_primitive_data(count,1) = i; % original basis function number
                    included_primitive_data(count,2) = j; % original basis function number
                    if i == j
                        included_primitive_data(count,3) = alpha_1PDM(i,j);
                        included_primitive_data(count,4) = beta_1PDM(i,j);
                    else
                        included_primitive_data(count,3) = 2*alpha_1PDM(i,j);
                        included_primitive_data(count,4) = 2*beta_1PDM(i,j);
                    end
                    alpha1=primitive_exponents(i,2);
                    alpha2=primitive_exponents(j,2);
                    X1=basis_set_centers(round(primitive_exponents(i,1)),3);
                    Y1=basis_set_centers(round(primitive_exponents(i,1)),4);
                    Z1=basis_set_centers(round(primitive_exponents(i,1)),5);
                    X2=basis_set_centers(round(primitive_exponents(j,1)),3) + repeata*vector1(1) + repeatb*vector2(1) + repeatc*vector3(1);
                    Y2=basis_set_centers(round(primitive_exponents(j,1)),4) + repeata*vector1(2) + repeatb*vector2(2) + repeatc*vector3(2);
                    Z2=basis_set_centers(round(primitive_exponents(j,1)),5) + repeata*vector1(3) + repeatb*vector2(3) + repeatc*vector3(3);         
                    Lx1 = primitive_exponents(i,3);
                    Ly1 = primitive_exponents(i,4);
                    Lz1 = primitive_exponents(i,5);
                    Lx2 = primitive_exponents(j,3);
                    Ly2 = primitive_exponents(j,4);
                    Lz2 = primitive_exponents(j,5); 
                    included_primitive_data(count,5) = pair_overlap(total_count);
                    included_primitive_data(count,6)=alpha1;
                    included_primitive_data(count,7)= Lx1;
                    included_primitive_data(count,8)= Ly1;
                    included_primitive_data(count,9)= Lz1;
                    included_primitive_data(count,10)=X1;
                    included_primitive_data(count,11)=Y1;
                    included_primitive_data(count,12)=Z1;
                    included_primitive_data(count,13)=alpha2;
                    included_primitive_data(count,14)=Lx2;
                    included_primitive_data(count,15)=Ly2;
                    included_primitive_data(count,16)=Lz2;
                    included_primitive_data(count,17)=X2;
                    included_primitive_data(count,18)=Y2;
                    included_primitive_data(count,19)=Z2;
                    max_coeff = abs(alpha_1PDM(i,j)) + abs(beta_1PDM(i,j));
                    r12 = sqrt((X1-X2)^2 + (Y1-Y2)^2 + (Z1-Z2)^2);
                    radius = sqrt(max((-log(gaussian_overlap_tolerance)/(alpha1+alpha2) - alpha1*alpha2*(r12^2)/((alpha1+alpha2)^2)),0.25));
                    included_primitive_data(count,20)=ceil(radius*100*delta_na/(cutoff_radius*bohrperangstrom)) + 1; % delta_na value for pair
                    included_primitive_data(count,21)=ceil(radius*100*delta_nb/(cutoff_radius*bohrperangstrom)) + 1; % delta_nb value for pair
                    included_primitive_data(count,22)=ceil(radius*100*delta_nb/(cutoff_radius*bohrperangstrom)) + 1; % delta_nc value for pair
                    % compute the center (na,nb,nc) position and the center shifts
                    XD=(alpha1*X1 + alpha2*X2)/(alpha1 + alpha2); % coordinates of pair product center
                    YD=(alpha1*Y1 + alpha2*Y2)/(alpha1 + alpha2);
                    ZD=(alpha1*Z1 + alpha2*Z2)/(alpha1 + alpha2);
                    included_primitive_data(count,23)=XD;
                    included_primitive_data(count,24)=YD;
                    included_primitive_data(count,25)=ZD;        
                    xyz=[XD,YD,ZD];
                    temp_vector = inv_boundary*(xyz' - origin');
                    included_primitive_data(count,26)=round(temp_vector(1)); % pair product center Na, Nb, Nc values
                    included_primitive_data(count,27)=round(temp_vector(2));
                    included_primitive_data(count,28)=round(temp_vector(3));
                    temp_center_shift(:) = boundary'*(temp_vector(:) - round(temp_vector(:))); % pair product xyz center shifts
                    included_primitive_data(count,29)=temp_center_shift(1); 
                    included_primitive_data(count,30)=temp_center_shift(2);
                    included_primitive_data(count,31)=temp_center_shift(3);
                    % Pairwise renormalization factor goes in column 32. Set
                    % initial estimate to 1.
                    included_primitive_data(count,32)=1.0;
                    % Whether the pair is treated numerically or analytically
                    if included_primitive_pairs(total_count) == 1
                       included_primitive_data(count,33)=1.0;
                    end
                end
            end
        end        
    end
end    
'The transformation of gaussian primitives was succesful.'
% Compute the valence_density and spin_density values at each grid point
% The first iteration computes the pairwise normalization factor, which is
% used to renormalize the pair contributions in the second pass, but only
% if the adjustment would lead to a change greater than the gaussian
% overlap tolerance.
valence_density=zeros(totnumA,totnumB,totnumC);
spin_density=zeros(totnumA,totnumB,totnumC);
core_density=zeros(totnumA,totnumB,totnumC);   
temp = 0.0;
alpha_smooth=1/(preferred_grid_spacing^2);
temp_pair_overlap=0.0;
for renormalization_step = 1:2
    for count=1:num_included_primitive_pairs
        if mod(count,100) == 0
            count
        end  
        analytic_pair_overlap=included_primitive_data(count,5);
        max_coeff=abs(included_primitive_data(count,3)) + abs(included_primitive_data(count,4));
        if (renormalization_step == 2) && (abs(max_coeff*analytic_pair_overlap)*included_primitive_data(count,32) < gaussian_overlap_tolerance)
            continue
        end  
        overlap_sum=0.0;
        abs_overlap_sum = 0.0;
        alpha1=included_primitive_data(count,6);
        Lx1=included_primitive_data(count,7);
        Ly1=included_primitive_data(count,8);
        Lz1=included_primitive_data(count,9);
        X1=included_primitive_data(count,10);
        Y1=included_primitive_data(count,11);
        Z1=included_primitive_data(count,12);
        alpha2=included_primitive_data(count,13);
        Lx2=included_primitive_data(count,14);
        Ly2=included_primitive_data(count,15);
        Lz2=included_primitive_data(count,16);
        X2=included_primitive_data(count,17);
        Y2=included_primitive_data(count,18);
        Z2=included_primitive_data(count,19); 
        delta_na_pair = included_primitive_data(count,20);
        delta_nb_pair = included_primitive_data(count,21);
        delta_nc_pair = included_primitive_data(count,22);
        Apoints=zeros(1,(2*delta_na_pair + 1));
        for na = -delta_na_pair:delta_na_pair
            if periodicA
                Apoints(delta_na_pair + na + 1) = round(mod((na + included_primitive_data(count,26)),totnumA) + 1);
            else
                Apoints(delta_na_pair + na + 1) = round(na + included_primitive_data(count,26) + 1);
            end
        end
        Bpoints=zeros(1,(2*delta_nb_pair + 1));
        for nb = -delta_nb_pair:delta_nb_pair
            if periodicB
                Bpoints(delta_nb_pair + nb + 1) = round(mod((nb + included_primitive_data(count,27)),totnumB) + 1);
            else
                Bpoints(delta_nb_pair + nb + 1) = round(nb + included_primitive_data(count,27) + 1);
            end
        end
        Cpoints=zeros(1,(2*delta_nc_pair + 1));
        for nc = -delta_nc_pair:delta_nc_pair
            if periodicC
                Cpoints(delta_nc_pair + nc + 1) = round(mod((nc + included_primitive_data(count,28)),totnumC) + 1);
            else
                Cpoints(delta_nc_pair + nc + 1) = round(nc + included_primitive_data(count,28) + 1);
            end
        end
        for na = -delta_na_pair:delta_na_pair
            ka = Apoints(delta_na_pair + na + 1);
            if ka < 1 || ka > totnumA
                continue
            end
            for nb = -delta_nb_pair:delta_nb_pair
                kb = Bpoints(delta_nb_pair + nb + 1);
                if kb < 1 || kb > totnumB
                    continue
                end
                for nc = -delta_nc_pair:delta_nc_pair
                    kc = Cpoints(delta_nc_pair + nc + 1);
                    if kc < 1 || kc > totnumC
                        continue
                    end          
                    X = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) + included_primitive_data(count,23) - included_primitive_data(count,29);
                    Y = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) + included_primitive_data(count,24) - included_primitive_data(count,30);
                    Z = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) + included_primitive_data(count,25) - included_primitive_data(count,31);
                    temp=(X-X1)^Lx1*(Y-Y1)^Ly1*(Z-Z1)^Lz1*(X-X2)^Lx2*(Y-Y2)^Ly2*(Z-Z2)^Lz2*exp(-alpha1*((X-X1)*(X-X1) + (Y-Y1)*(Y-Y1) + (Z-Z1)*(Z-Z1)) - alpha2*((X-X2)*(X-X2) + (Y-Y2)*(Y-Y2) + (Z-Z2)*(Z-Z2))); 
                    overlap_sum = overlap_sum + temp*pixelvolume;
                    abs_overlap_sum = abs_overlap_sum + abs(temp)*pixelvolume;
                    if renormalization_step == 1
                       temp_scalar=1.0;
                    else
                       temp_scalar=included_primitive_data(count,32);
                       temp = abs(temp);
                    end
                    spin_contribution = temp_scalar*(included_primitive_data(count,3) - included_primitive_data(count,4))*temp;
                    spin_density(ka,kb,kc) = spin_density(ka,kb,kc) + spin_contribution;
                    if included_primitive_data(count,33) == 0
                       valence_density(ka,kb,kc) = valence_density(ka,kb,kc) + temp_scalar*(included_primitive_data(count,3) + included_primitive_data(count,4))*temp;
                    elseif (included_primitive_data(count,33) == 1)
                       core_density(ka,kb,kc) = core_density(ka,kb,kc) + temp_scalar*(included_primitive_data(count,3) + included_primitive_data(count,4))*temp;
                       basis_index_for_i =included_primitive_data(count,1);
                       basis_index_for_j =included_primitive_data(count,2);
                       center_index_for_i=primitive_exponents(basis_index_for_i,1);
                       center_index_for_j=primitive_exponents(basis_index_for_j,1);
                       atom_index_for_i = basis_set_centers(center_index_for_i,6);
                       atom_index_for_j = basis_set_centers(center_index_for_j,6); 
                       if alpha1*min(basis_set_centers(center_index_for_i,2),1) > alpha2*min(basis_set_centers(center_index_for_j,2),1)
                            partial_pair_weight_i = 1;
                            partial_pair_weight_j = 0;
                       elseif alpha1*min(basis_set_centers(center_index_for_i,2),1) < alpha2*min(basis_set_centers(center_index_for_j,2),1)
                            partial_pair_weight_i = 0;
                            partial_pair_weight_j = 1;
                       else
                            partial_pair_weight_i = 0.5;
                            partial_pair_weight_j = 0.5;                
                       end   
                       if atom_index_for_i > 0.5
                            % spin magnetic moment correction
                            occupancy_correction(atom_index_for_i,2) = occupancy_correction(atom_index_for_i,2) - spin_contribution*pixelvolume*partial_pair_weight_i;
                       end
                       if atom_index_for_j > 0.5
                            % spin magnetic moment correction
                            occupancy_correction(atom_index_for_j,2) = occupancy_correction(atom_index_for_j,2) - spin_contribution*pixelvolume*partial_pair_weight_j;
                       end                              
                    end            
                end       
            end
        end
        % Update the renormalization factor
        if ((included_primitive_data(count,33) == 0) && (abs_overlap_sum > gaussian_overlap_tolerance))
           % the following line makes sure the correction is not more than 5 percent of the overlap sum
           included_primitive_data(count,32)=max(min(((analytic_pair_overlap - overlap_sum)/abs_overlap_sum),0.05),-0.05); 
        else
            included_primitive_data(count,32) = 0;
        end    
    end
end   
'Generation of the valence and spin density for each grid point was successful.'
% Generate the core density grid
if core_available == 1
    for j=1:n_edf_primitives
        edf_primitives(j,7)= basis_set_centers(edf_primitives(j,1),3); %edf center X
        edf_primitives(j,8)= basis_set_centers(edf_primitives(j,1),4); %edf center Y
        edf_primitives(j,9)= basis_set_centers(edf_primitives(j,1),5); %edf center Z
        xyz=[edf_primitives(j,7),edf_primitives(j,8),edf_primitives(j,9)];
        temp_vector = inv_boundary*(xyz' - origin');
        edf_primitives(j,10)=round(temp_vector(1)); %edf center na
        edf_primitives(j,11)=round(temp_vector(2)); %edf center nb
        edf_primitives(j,12)=round(temp_vector(3)); %edf center nc
        temp_center_shift(:) = boundary'*(temp_vector(:) - round(temp_vector(:)));
        edf_primitives(j,13)= temp_center_shift(1); %edf center shift x
        edf_primitives(j,14)= temp_center_shift(2); %edf center shift y
        edf_primitives(j,15)= temp_center_shift(3); %edf center shift z
    end
    for j=1:n_edf_primitives
        Apoints=zeros(1,(2*delta_na + 1));
        for na = -delta_na:delta_na
            if periodicA
                Apoints(delta_na + na + 1) = round(mod((na + edf_primitives(j,10)),totnumA) + 1);
            else
                Apoints(delta_na + na + 1) = round(na + edf_primitives(j,10) + 1);
            end
        end
        Bpoints=zeros(1,(2*delta_nb + 1));
        for nb = -delta_nb:delta_nb
            if periodicB
                Bpoints(delta_nb + nb + 1) = round(mod((nb + edf_primitives(j,11)),totnumB) + 1);
            else
                Bpoints(delta_nb + nb + 1) = round(nb + edf_primitives(j,11) + 1);
            end
        end
        Cpoints=zeros(1,(2*delta_nc + 1));
        for nc = -delta_nc:delta_nc
            if periodicC
                Cpoints(delta_nc + nc + 1) = round(mod((nc + edf_primitives(j,12)),totnumC) + 1);
            else
                Cpoints(delta_nc + nc + 1) = round(nc + edf_primitives(j,12) + 1);
            end
        end
        for na = -delta_na:delta_na
            ka = Apoints(delta_na + na + 1);
            if ka < 1 || ka > totnumA
                continue
            end
            for nb = -delta_nb:delta_nb
                kb = Bpoints(delta_nb + nb + 1);
                if kb < 1 || kb > totnumB
                    continue
                end
                for nc = -delta_nc:delta_nc
                    kc = Cpoints(delta_nc + nc + 1);
                    if kc < 1 || kc > totnumC
                        continue
                    end          
                    temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - edf_primitives(j,13);
                    temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - edf_primitives(j,14);
                    temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - edf_primitives(j,15);
                    distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3));
                    temp_index = ceil(scalefactor*distance + zero_tolerance);
                    if temp_index <= nshells
                        alpha=edf_primitives(j,2);
                        Lx=edf_primitives(j,3);
                        Ly=edf_primitives(j,4);
                        Lz=edf_primitives(j,5);
                        Ax=edf_primitives(j,7);
                        Ay=edf_primitives(j,8);
                        Az=edf_primitives(j,9);
                        X = temp_vector(1) + edf_primitives(j,7);
                        Y = temp_vector(2) + edf_primitives(j,8);
                        Z = temp_vector(3) + edf_primitives(j,9);
                        temp_scalar = gaussian_value(edf_primitives(j,6),alpha,Lx,Ly,Lz,Ax,Ay,Az,X,Y,Z);
                        core_density(ka,kb,kc) = core_density(ka,kb,kc) + temp_scalar;
                    end
                end    
            end
        end
    end     
end    
if flag == 1
    break
end
% Make sure the core and total electron densities are positive definite
for ka=1:totnumA
    for kb=1:totnumB
        for kc=1:totnumC
            core_density(ka,kb,kc) = max(core_density(ka,kb,kc),0.0);
            valence_density(ka,kb,kc) = max(valence_density(ka,kb,kc),-core_density(ka,kb,kc));
        end
    end
end
add_missing_core_density
initialize_atomic_densities
if flag == 1
    break
end
sum_negative_density=0.0;