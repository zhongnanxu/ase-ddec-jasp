function Xi = calculate_Xi(tau)

threshold = 10^-3;

if tau < 0
    Xi = 0.0;
elseif tau > 1
    Xi = 2*pi;
elseif tau < threshold
    Xi = 4*pi*tau/3;
elseif tau > (1 - threshold)
    Xi = 2*pi + (tau - 1)*49.956 + (tau - 1)^2*8370.6;
else    
    Xi = pi*((tau^(-2) - 1)*log((1 - tau)/(1 + tau)) + 2/tau);
end      