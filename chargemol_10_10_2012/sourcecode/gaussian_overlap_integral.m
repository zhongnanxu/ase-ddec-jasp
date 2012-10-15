function value = gaussian_overlap_integral(alpha1,Lx1,Ly1,Lz1,X1,Y1,Z1,alpha2,Lx2,Ly2,Lz2,X2,Y2,Z2)

% Computes the spatial overlap integral between two gaussians f1 and f2
% f1 =(X-X1)^Lx1*(Y-Y1)^Ly1*(Z-Z1)^Lz1*exp(-alpha1*((X-X1)^2 +(Y-Y1)^2 +(Z-Z1)^2))
% f2 =(X-X2)^Lx2*(Y-Y2)^Ly2*(Z-Z2)^Lz2*exp(-alpha2*((X-X2)^2 +(Y-Y2)^2 +(Z-Z2)^2)) 
factor1=exp(-alpha1*alpha2/(alpha1 + alpha2)*((X2-X1)^2 + (Y2-Y1)^2 + (Z2-Z1)^2));
XD=(alpha1*X1 + alpha2*X2)/(alpha1 + alpha2);
YD=(alpha1*Y1 + alpha2*Y2)/(alpha1 + alpha2);
ZD=(alpha1*Z1 + alpha2*Z2)/(alpha1 + alpha2);
factor2=0;
for i=0:Lx1
    for j=0:Lx2
        if mod((i+j),2) == 1
           continue
        else
           temp=(factorial(Lx1)/(factorial(Lx1-i)*factorial(i)));
           temp=temp*(factorial(Lx2)/(factorial(Lx2-j)*factorial(j)));
           temp=temp*((XD - X1)^(Lx1-i))*((XD - X2)^(Lx2-j));
           SART = gaussian_integration((alpha1+alpha2),(i+j));
        end
        factor2 = factor2 + temp*SART;
    end
end
factor3=0;
for i=0:Ly1
    for j=0:Ly2
        if mod((i+j),2) == 1
           continue
        else
           temp=(factorial(Ly1)/(factorial(Ly1-i)*factorial(i)));
           temp=temp*(factorial(Ly2)/(factorial(Ly2-j)*factorial(j)));
           temp=temp*((YD - Y1)^(Ly1-i))*((YD - Y2)^(Ly2-j));
           SART = gaussian_integration((alpha1+alpha2),(i+j));
        end
        factor3 = factor3 + temp*SART;
    end
end
factor4=0;
for i=0:Lz1
    for j=0:Lz2
        if mod((i+j),2) == 1
           continue
        else
           temp=(factorial(Lz1)/(factorial(Lz1-i)*factorial(i)));
           temp=temp*(factorial(Lz2)/(factorial(Lz2-j)*factorial(j)));
           temp=temp*((ZD - Z1)^(Lz1-i))*((ZD - Z2)^(Lz2-j));
           SART = gaussian_integration((alpha1+alpha2),(i+j));
        end
        factor4 = factor4 + temp*SART;
    end
end
value=factor1*factor2*factor3*factor4;
