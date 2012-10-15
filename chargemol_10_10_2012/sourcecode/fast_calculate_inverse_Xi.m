function inv_Xi = fast_calculate_inverse_Xi(Xi_value,inverse_Xi_lookup)

num_lookup_points = length(inverse_Xi_lookup);
threshold = 1/num_lookup_points;
scaled_Xi_value = Xi_value/(2*pi);
if scaled_Xi_value <= 0.0
    inv_Xi = 0.0;
elseif scaled_Xi_value >= 1
    inv_Xi = 1.0;
elseif scaled_Xi_value < threshold
    inv_Xi = Xi_value*0.75/pi;
else  
    i = floor(scaled_Xi_value*num_lookup_points);
    f = scaled_Xi_value*num_lookup_points - i;
    inv_Xi = (1-f)*inverse_Xi_lookup(i) + f*inverse_Xi_lookup(i+1);
end    