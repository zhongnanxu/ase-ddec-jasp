function value = gaussian_value(constant,alpha,Lx,Ly,Lz,Ax,Ay,Az,x,y,z)

value = constant*((x-Ax)^Lx)*((y-Ay)^Ly)*((z-Az)^Lz)*exp(-alpha*((x-Ax)^2 + (y-Ay)^2 + (z-Az)^2));
   
