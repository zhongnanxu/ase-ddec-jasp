standard_atomic_weights;
sumvector=zeros(1,3);
Mtot=0.0;
for j=1:natoms
   sumvector(1,1) = sumvector(1,1) + atomic_weight(atomic_number(j))*coords(j,1); 
   sumvector(1,2) = sumvector(1,2) + atomic_weight(atomic_number(j))*coords(j,2);
   sumvector(1,3) = sumvector(1,3) + atomic_weight(atomic_number(j))*coords(j,3);
   Mtot = Mtot + atomic_weight(atomic_number(j));
end
center_of_mass=sumvector/Mtot