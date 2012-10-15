% Compute the positions of the charge centers in terms of na, nb, nc
inv_boundary=inv(boundary');
center_nabc=zeros(natoms,3);
xyz = [0.00,0.00,0.00];
for i = 1:natoms
    xyz(:)=coords(i,:);
    temp_vector = inv_boundary*(xyz' - origin');
    center_nabc(i,:) = round(temp_vector(:));
    center_shift(i,:) = boundary'*(temp_vector(:) - round(temp_vector(:)));
end
center_nabc
 