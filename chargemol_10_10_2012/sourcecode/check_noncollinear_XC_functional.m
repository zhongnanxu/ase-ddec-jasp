% This program separately computes the integrated changes of the directional and magnitude
% components to the spin magnetization density. These quantities can be
% used to determine whether a two or four-component non-collinear DFT
% functional is needed.
%
sum_sq_gradient_mag = 0.0;
sum_sq_gradient_tot = 0.0;
for ka = 1:(totnumA - 1)
    for kb = 1:(totnumB - 1)
        for kc = 1:(totnumC - 1)
            % spin density magnitude changes
            spin_density_vector_magnitude = sqrt(spin_density_vector(ka,kb,kc,1)^2 + spin_density_vector(ka,kb,kc,2)^2 + spin_density_vector(ka,kb,kc,3)^2);
            spin_density_vector_magnitude_a = sqrt(spin_density_vector((ka+1),kb,kc,1)^2 + spin_density_vector((ka+1),kb,kc,2)^2 + spin_density_vector((ka+1),kb,kc,3)^2);
            spin_density_vector_magnitude_b = sqrt(spin_density_vector(ka,(kb+1),kc,1)^2 + spin_density_vector(ka,(kb+1),kc,2)^2 + spin_density_vector(ka,(kb+1),kc,3)^2);
            spin_density_vector_magnitude_c = sqrt(spin_density_vector(ka,kb,(kc+1),1)^2 + spin_density_vector(ka,kb,(kc+1),2)^2 + spin_density_vector(ka,kb,(kc+1),3)^2);
            delta_magnitude_x = (spin_density_vector_magnitude_a - spin_density_vector_magnitude)*inv_boundary(1,1) + (spin_density_vector_magnitude_b - spin_density_vector_magnitude)*inv_boundary(2,1) + (spin_density_vector_magnitude_c - spin_density_vector_magnitude)*inv_boundary(3,1);
            delta_magnitude_y = (spin_density_vector_magnitude_a - spin_density_vector_magnitude)*inv_boundary(1,2) + (spin_density_vector_magnitude_b - spin_density_vector_magnitude)*inv_boundary(2,2) + (spin_density_vector_magnitude_c - spin_density_vector_magnitude)*inv_boundary(3,2);
            delta_magnitude_z = (spin_density_vector_magnitude_a - spin_density_vector_magnitude)*inv_boundary(1,3) + (spin_density_vector_magnitude_b - spin_density_vector_magnitude)*inv_boundary(2,3) + (spin_density_vector_magnitude_c - spin_density_vector_magnitude)*inv_boundary(3,3);
            sq_gradient_mag = delta_magnitude_x^2 + delta_magnitude_y^2 + delta_magnitude_z^2;
            % mx changes 
            delta_mx_x = (spin_density_vector((ka+1),kb,kc,1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(1,1) + (spin_density_vector(ka,(kb+1),kc,1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(2,1) + (spin_density_vector(ka,kb,(kc+1),1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(3,1);
            delta_mx_y = (spin_density_vector((ka+1),kb,kc,1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(1,2) + (spin_density_vector(ka,(kb+1),kc,1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(2,2) + (spin_density_vector(ka,kb,(kc+1),1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(3,2);
            delta_mx_z = (spin_density_vector((ka+1),kb,kc,1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(1,3) + (spin_density_vector(ka,(kb+1),kc,1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(2,3) + (spin_density_vector(ka,kb,(kc+1),1) - spin_density_vector(ka,kb,kc,1))*inv_boundary(3,3);             
            sq_gradient_mx = delta_mx_x^2 + delta_mx_y^2 + delta_mx_z^2;            
            % my changes 
            delta_my_x = (spin_density_vector((ka+1),kb,kc,2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(1,1) + (spin_density_vector(ka,(kb+1),kc,2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(2,1) + (spin_density_vector(ka,kb,(kc+1),2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(3,1);
            delta_my_y = (spin_density_vector((ka+1),kb,kc,2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(1,2) + (spin_density_vector(ka,(kb+1),kc,2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(2,2) + (spin_density_vector(ka,kb,(kc+1),2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(3,2);
            delta_my_z = (spin_density_vector((ka+1),kb,kc,2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(1,3) + (spin_density_vector(ka,(kb+1),kc,2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(2,3) + (spin_density_vector(ka,kb,(kc+1),2) - spin_density_vector(ka,kb,kc,2))*inv_boundary(3,3);             
            sq_gradient_my = delta_my_x^2 + delta_my_y^2 + delta_my_z^2; 
            % mz changes 
            delta_mz_x = (spin_density_vector((ka+1),kb,kc,3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(1,1) + (spin_density_vector(ka,(kb+1),kc,3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(2,1) + (spin_density_vector(ka,kb,(kc+1),3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(3,1);
            delta_mz_y = (spin_density_vector((ka+1),kb,kc,3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(1,2) + (spin_density_vector(ka,(kb+1),kc,3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(2,2) + (spin_density_vector(ka,kb,(kc+1),3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(3,2);
            delta_mz_z = (spin_density_vector((ka+1),kb,kc,3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(1,3) + (spin_density_vector(ka,(kb+1),kc,3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(2,3) + (spin_density_vector(ka,kb,(kc+1),3) - spin_density_vector(ka,kb,kc,3))*inv_boundary(3,3);             
            sq_gradient_mz = delta_mz_x^2 + delta_mz_y^2 + delta_mz_z^2;
            if spin_density_vector_magnitude > zero_tolerance
               sum_sq_gradient_mag = sum_sq_gradient_mag + sq_gradient_mag*pixelvolume;
               sum_sq_gradient_tot = sum_sq_gradient_tot + (sq_gradient_mx + sq_gradient_my + sq_gradient_mz)*pixelvolume;
            end
        end
    end     
end               
% 
'Changes in the magnitude of the spin magnetization density contribute an integrated weight of '
sum_sq_gradient_mag
'The total changes in the spin magnetization density contribute an integrated weight of'
sum_sq_gradient_tot
'Based on these calculations, directional changes of the spin magnetization density contribe an integrated weight of'
sum_sq_gradient_dir = max((sum_sq_gradient_tot - sum_sq_gradient_mag), 0.0)
