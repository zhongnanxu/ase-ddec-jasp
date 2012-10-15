function inv_Xi = calculate_inverse_Xi(Xi_value)

zero_tolerance = 10^-6;

if Xi_value < 0.0
    inv_Xi = 0.0;
elseif Xi_value < zero_tolerance
    inv_Xi = 3*Xi_value/(4*pi);
elseif Xi_value > (2*pi - zero_tolerance)
    inv_Xi = 1.0;
else    
    tau_estimate = (Xi_value - 0.036835*Xi_value*Xi_value)/(4.178319 - 0.136129*Xi_value + 0.038148*Xi_value*Xi_value);
    Xi_derivative = calculate_Xi_derivative(tau_estimate);
    for i = 1:4
        Xi_estimate = calculate_Xi(tau_estimate);
        tau_estimate = tau_estimate + (Xi_value - Xi_estimate)/Xi_derivative;
    end 
    inv_Xi = tau_estimate;
end    