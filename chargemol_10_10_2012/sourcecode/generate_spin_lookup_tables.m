num_lookup_points = 10000;
Xi_lookup = zeros(num_lookup_points,1);
inverse_Xi_lookup = zeros(num_lookup_points,1);
for i=1:num_lookup_points
    Xi_lookup(i) = calculate_Xi(i/num_lookup_points);
    inverse_Xi_lookup(i) = calculate_inverse_Xi(i*2*pi/num_lookup_points);
end    