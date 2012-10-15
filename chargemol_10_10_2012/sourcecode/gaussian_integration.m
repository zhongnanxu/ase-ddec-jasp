function value = gaussian_integration(alpha,Lx)

temp=1;
if mod(Lx,2) == 1 
    value = 0;
else
    for i=1:2:Lx
        temp=temp*i;
    end
    value =sqrt(pi/alpha)*temp/((2*alpha)^((Lx)/2));
end    
