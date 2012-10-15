function Xi_derivative = calculate_Xi_derivative(tau)

threshold = 10^-3;

if tau < threshold
    Xi_derivative = 4*pi/3;
elseif tau > (1 - threshold)
    Xi_derivative = 50.0;
else    
    Xi_derivative = 2*pi*(log((1+tau)/(1-tau)) - 2*tau)/(tau^3);
end    