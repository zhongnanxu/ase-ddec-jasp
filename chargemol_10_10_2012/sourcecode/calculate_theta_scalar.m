function theta_scalar = calculate_theta_scalar(a,b,Xi_lookup)

zero_tolerance = 10^-6;
if (abs(b) < zero_tolerance) || (a < zero_tolerance)
    theta_scalar = 0.0;
else
    tau = (abs(b)/a);
    Xi = fast_calculate_Xi(tau,Xi_lookup);
    theta_scalar = b*Xi/(2*abs(b));
end    