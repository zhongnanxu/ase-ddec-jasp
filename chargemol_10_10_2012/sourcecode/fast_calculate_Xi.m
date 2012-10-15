function Xi = fast_calculate_Xi(tau,Xi_lookup)

num_lookup_points = length(Xi_lookup);
threshold = 1/num_lookup_points;
if tau <= 0.0
    Xi = 0.0;
elseif tau >= 1
    Xi = 2*pi;
elseif tau < threshold
    Xi = 4*pi*tau/3;
else  
    i = floor(tau*num_lookup_points);
    f = tau*num_lookup_points - i;
    Xi = (1-f)*Xi_lookup(i) + f*Xi_lookup(i+1);
end    