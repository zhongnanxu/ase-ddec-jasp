% compute the three side lengths of a parallelepiped that enclosed each
% atom, where each side of the parallelpiped lies along the chosen axes
temp_vector1=zeros(1,3);
temp_vector2=zeros(1,3);
temp_vector3=zeros(1,3);
temp_vector1(:) = boundary(1,:);
temp_vector2(:) = boundary(2,:);
temp_vector3(:) = boundary(3,:);
a_dot_a = temp_vector1*temp_vector1';
a_dot_b = temp_vector1*temp_vector2';
a_dot_c = temp_vector1*temp_vector3';
b_dot_b = temp_vector2*temp_vector2';
b_dot_c = temp_vector2*temp_vector3';
c_dot_c = temp_vector3*temp_vector3';
R = nshells/scalefactor;
cos_a_b = (a_dot_b)/sqrt(a_dot_a*b_dot_b);
cos_a_c = (a_dot_c)/sqrt(a_dot_a*c_dot_c);
cos_b_c = (b_dot_c)/sqrt(b_dot_b*c_dot_c);
sin_a_b = sqrt(1 - cos_a_b*cos_a_b);
sin_a_c = sqrt(1 - cos_a_c*cos_a_c);
sin_b_c = sqrt(1 - cos_b_c*cos_b_c);
delta_na = ceil(max((R/sin_a_b), (R/sin_a_c))/sqrt(a_dot_a)) + 1
delta_nb = ceil(max((R/sin_a_b), (R/sin_b_c))/sqrt(b_dot_b)) + 1
delta_nc = ceil(max((R/sin_a_c), (R/sin_b_c))/sqrt(c_dot_c)) + 1
% 