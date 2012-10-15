function theta_vector = calculate_theta_vector(a,b_vector,Xi_lookup)

zero_tolerance = 10^-6;
mag_b_vector = sqrt(b_vector*b_vector');
if (mag_b_vector > zero_tolerance) && (a > zero_tolerance)
    tau = (mag_b_vector/a);
    Xi = fast_calculate_Xi(tau,Xi_lookup);
    theta_vector = b_vector*Xi/(2*mag_b_vector);
else
    theta_vector = zeros(1,3);
end    