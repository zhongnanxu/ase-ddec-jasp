% open and read the xsf file containing the valence density
spin_available=0;
num_periodic_directions=0;
core_available=0;
fid1 = fopen(wfx_inputfile,'r');
if fid1 <= 0
    'Could not find input wfx file. Program will terminate.'
    flag=1
    break
else
    fileInfo = dir(wfx_inputfile);
    fileSize = fileInfo.bytes;
    if fileSize < 10
       'Input wfx file appears to be empty. Program will terminate.'
       flag=1
       break
    end 
end
while ~feof(fid1)
   if flag == 1
       break
   end
   data = fgetl(fid1);
   if (length(lower(strtrim(data))) == length('<keywords>'))
       if (lower(strtrim(data)) == '<keywords>')
           basis_set_type = fgetl(fid1);
           if (length(lower(strtrim(basis_set_type))) == length('gto'))
              if (lower(strtrim(basis_set_type)) == 'gto')
                  'GTO basis set type is supported. Calculation will proceed.'
              else
                  'Basis set type is not supported. Calculation will terminate.'
                  flag=1;
                  break
              end
           end   
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Net Charge>'))
       if (lower(strtrim(data)) == '<net charge>')
           temp_string = fgetl(fid1);
           netcharge=str2num(temp_string)
           clear temp_string
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Number of Translation Vectors>'))
       if (lower(strtrim(data)) == '<number of translation vectors>')
           temp_string = fgetl(fid1);
           num_periodic_directions=str2num(temp_string)
           clear temp_string
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Translation Vectors>'))
       if (lower(strtrim(data)) == '<translation vectors>')
           periodic_vectors=zeros(num_periodic_directions,3);
           for i=1:num_periodic_directions
               data = textscan(fid1, '%f %f %f',1);
               temp_array=cell2mat(data);
               periodic_vectors(i,1)=temp_array(1);
               periodic_vectors(i,2)=temp_array(2);
               periodic_vectors(i,3)=temp_array(3);
           end    
           clear data temp_array;
           periodic_vectors
           continue
       end
   end   
   if (length(lower(strtrim(data))) == length('<Number of Nuclei>'))
       if (lower(strtrim(data)) == '<number of nuclei>')
           temp_string = fgetl(fid1);
           ncenters=str2num(temp_string)
           basis_set_centers=zeros(ncenters,6);
           clear temp_string
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Atomic Numbers>'))
       if (lower(strtrim(data)) == '<atomic numbers>')          
           for i=1:ncenters
               temp_string = fgetl(fid1);
               basis_set_centers(i,1)=int16(str2num(temp_string));
               clear temp_string
           end    
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Nuclear Charges>'))
       if (lower(strtrim(data)) == '<nuclear charges>')          
           for i=1:ncenters
               temp_string = fgetl(fid1);
               basis_set_centers(i,2)=str2num(temp_string);
               clear temp_string
           end    
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Nuclear Cartesian Coordinates>'))
       if (lower(strtrim(data)) == '<nuclear cartesian coordinates>')          
           data = textscan(fid1, '%f %f %f',ncenters);
           temp_array=cell2mat(data);
           basis_set_centers(:,3) = temp_array(:,1);
           basis_set_centers(:,4) = temp_array(:,2);
           basis_set_centers(:,5) = temp_array(:,3);    
           clear data temp_array;
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Number of Primitives>'))
       if (lower(strtrim(data)) == '<number of primitives>')          
           temp_string = fgetl(fid1);
           nprimitives=str2num(temp_string)
           primitive_exponents=zeros(nprimitives,5);
           clear temp_string
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<Primitive Centers>'))
       if (lower(strtrim(data)) == '<primitive centers>')          
           data = textscan(fid1, '%d',nprimitives);
           temp_array=cell2mat(data);
           for i=1:nprimitives
               primitive_exponents(i,1)=temp_array(i);
           end    
           clear data temp_array
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<Primitive Exponents>'))
       if (lower(strtrim(data)) == '<primitive exponents>')          
           data = textscan(fid1, '%f',nprimitives);
           temp_array=cell2mat(data);
           for i=1:nprimitives
               primitive_exponents(i,2)=temp_array(i);
           end    
           clear data temp_array
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<Primitive Types>'))
       if (lower(strtrim(data)) == '<primitive types>')          
           data = textscan(fid1, '%d',nprimitives);
           temp_array=cell2mat(data);
           for i=1:nprimitives
               if temp_array(i) == 1
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 2
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 3
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 0;                
               elseif temp_array(i) == 4
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 5
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 6
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 0;                  
               elseif temp_array(i) == 7
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 2;
               elseif temp_array(i) == 8
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 0;                
               elseif temp_array(i) == 9
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 10
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 11
                  primitive_exponents(i,3) = 3;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 0;                 
               elseif temp_array(i) == 12
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 3;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 13
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 3;                
               elseif temp_array(i) == 14
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 15
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 16
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 1;                  
               elseif temp_array(i) == 17
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 18
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 2;                
               elseif temp_array(i) == 19
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 2;
               elseif temp_array(i) == 20
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 21
                  primitive_exponents(i,3) = 4;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 0;       
               elseif temp_array(i) == 22
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 4;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 23
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 4;                
               elseif temp_array(i) == 24
                  primitive_exponents(i,3) = 3;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 25
                  primitive_exponents(i,3) = 3;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 26
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 3;
                  primitive_exponents(i,5) = 0;                  
               elseif temp_array(i) == 27
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 3;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 28
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 3;                
               elseif temp_array(i) == 29
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 3;
               elseif temp_array(i) == 30
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 31
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 2;         
               elseif temp_array(i) == 32
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 2;
               elseif temp_array(i) == 33
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 1;                
               elseif temp_array(i) == 34
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 35
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 2;
               elseif temp_array(i) == 36
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 5;                  
               elseif temp_array(i) == 37
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 4;
               elseif temp_array(i) == 38
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 3;                
               elseif temp_array(i) == 39
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 3;
                  primitive_exponents(i,5) = 2;
               elseif temp_array(i) == 40
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 4;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 41
                  primitive_exponents(i,3) = 0;
                  primitive_exponents(i,4) = 5;
                  primitive_exponents(i,5) = 0;                 
               elseif temp_array(i) == 42
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 4;
               elseif temp_array(i) == 43
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 3;                
               elseif temp_array(i) == 44
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 2;
               elseif temp_array(i) == 45
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 3;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 46
                  primitive_exponents(i,3) = 1;
                  primitive_exponents(i,4) = 4;
                  primitive_exponents(i,5) = 0;                  
               elseif temp_array(i) == 47
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 3;
               elseif temp_array(i) == 48
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 2;                
               elseif temp_array(i) == 49
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 50
                  primitive_exponents(i,3) = 2;
                  primitive_exponents(i,4) = 3;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 51
                  primitive_exponents(i,3) = 3;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 2;       
               elseif temp_array(i) == 52
                  primitive_exponents(i,3) = 3;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 53
                  primitive_exponents(i,3) = 3;
                  primitive_exponents(i,4) = 2;
                  primitive_exponents(i,5) = 0;                
               elseif temp_array(i) == 54
                  primitive_exponents(i,3) = 4;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 1;
               elseif temp_array(i) == 55
                  primitive_exponents(i,3) = 4;
                  primitive_exponents(i,4) = 1;
                  primitive_exponents(i,5) = 0;
               elseif temp_array(i) == 56
                  primitive_exponents(i,3) = 5;
                  primitive_exponents(i,4) = 0;
                  primitive_exponents(i,5) = 0;
               else
                   'Only basis set primitives up to H shells (L=5) are supported. Program will terminate.'
                   flag=1;
                   break
               end
           end    
           clear data temp_array
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<Number of EDF Primitives>'))
       if (lower(strtrim(data)) == '<number of edf primitives>')          
           temp_string = fgetl(fid1);
           n_edf_primitives=str2num(temp_string)
           edf_primitives=zeros(n_edf_primitives,15);
           core_available=1;
           clear temp_string
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<EDF Primitive Centers>'))
       if (lower(strtrim(data)) == '<edf primitive centers>')          
           data = textscan(fid1, '%d',n_edf_primitives);
           temp_array=cell2mat(data);
           for i=1:n_edf_primitives
               edf_primitives(i,1)=temp_array(i);
           end    
           clear data temp_array
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<EDF Primitive Exponents>'))
       if (lower(strtrim(data)) == '<edf primitive exponents>')          
           data = textscan(fid1, '%f',n_edf_primitives);
           temp_array=cell2mat(data);
           for i=1:n_edf_primitives
               edf_primitives(i,2)=temp_array(i);
           end    
           clear data temp_array
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<EDF Primitive Coefficients>'))
       if (lower(strtrim(data)) == '<edf primitive coefficients>')          
           data = textscan(fid1, '%f',n_edf_primitives);
           temp_array=cell2mat(data);
           for i=1:n_edf_primitives
               edf_primitives(i,6)=temp_array(i);
           end    
           clear data temp_array
           continue
       end   
   end   
   if (length(lower(strtrim(data))) == length('<EDF Primitive Types>'))
       if (lower(strtrim(data)) == '<edf primitive types>')          
           data = textscan(fid1, '%d',n_edf_primitives);
           temp_array=cell2mat(data);
           for i=1:n_edf_primitives
               if temp_array(i) == 1
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 2
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 3
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 0;                
               elseif temp_array(i) == 4
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 5
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 6
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 0;                  
               elseif temp_array(i) == 7
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 2;
               elseif temp_array(i) == 8
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 0;                
               elseif temp_array(i) == 9
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 10
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 11
                  edf_primitives(i,3) = 3;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 0;                 
               elseif temp_array(i) == 12
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 3;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 13
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 3;                
               elseif temp_array(i) == 14
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 15
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 16
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 1;                  
               elseif temp_array(i) == 17
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 18
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 2;                
               elseif temp_array(i) == 19
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 2;
               elseif temp_array(i) == 20
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 21
                  edf_primitives(i,3) = 4;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 0;       
               elseif temp_array(i) == 22
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 4;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 23
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 4;                
               elseif temp_array(i) == 24
                  edf_primitives(i,3) = 3;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 25
                  edf_primitives(i,3) = 3;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 26
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 3;
                  edf_primitives(i,5) = 0;                  
               elseif temp_array(i) == 27
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 3;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 28
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 3;                
               elseif temp_array(i) == 29
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 3;
               elseif temp_array(i) == 30
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 31
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 2;         
               elseif temp_array(i) == 32
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 2;
               elseif temp_array(i) == 33
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 1;                

               elseif temp_array(i) == 34
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 35
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 2;
               elseif temp_array(i) == 36
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 5;                  
               elseif temp_array(i) == 37
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 4;
               elseif temp_array(i) == 38
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 3;                
               elseif temp_array(i) == 39
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 3;
                  edf_primitives(i,5) = 2;
               elseif temp_array(i) == 40
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 4;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 41
                  edf_primitives(i,3) = 0;
                  edf_primitives(i,4) = 5;
                  edf_primitives(i,5) = 0;                 
               elseif temp_array(i) == 42
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 4;
               elseif temp_array(i) == 43
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 3;                
               elseif temp_array(i) == 44
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 2;
               elseif temp_array(i) == 45
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 3;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 46
                  edf_primitives(i,3) = 1;
                  edf_primitives(i,4) = 4;
                  edf_primitives(i,5) = 0;                  
               elseif temp_array(i) == 47
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 3;
               elseif temp_array(i) == 48
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 2;                
               elseif temp_array(i) == 49
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 50
                  edf_primitives(i,3) = 2;
                  edf_primitives(i,4) = 3;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 51
                  edf_primitives(i,3) = 3;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 2;       
               elseif temp_array(i) == 52
                  edf_primitives(i,3) = 3;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 53
                  edf_primitives(i,3) = 3;
                  edf_primitives(i,4) = 2;
                  edf_primitives(i,5) = 0;                
               elseif temp_array(i) == 54
                  edf_primitives(i,3) = 4;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 1;
               elseif temp_array(i) == 55
                  edf_primitives(i,3) = 4;
                  edf_primitives(i,4) = 1;
                  edf_primitives(i,5) = 0;
               elseif temp_array(i) == 56
                  edf_primitives(i,3) = 5;
                  edf_primitives(i,4) = 0;
                  edf_primitives(i,5) = 0;
               else
                   'Only basis set primitives up to H shells (L=5) are supported. Program will terminate.'
                   flag=1;
                   break
               end
           end    
           clear data temp_array
           continue
       end   
   end
   if (length(lower(strtrim(data))) == length('<Number of Occupied Molecular Orbitals>'))
       if (lower(strtrim(data)) == '<number of occupied molecular orbitals>')
           temp_string = fgetl(fid1);
           norbitals=str2num(temp_string)
           clear temp_string
           continue
       end
   end  
   if (length(lower(strtrim(data))) == length('<Molecular Orbital Occupation Numbers>'))
       if (lower(strtrim(data)) == '<molecular orbital occupation numbers>')
           temp_occupation=zeros(norbitals,1);
           alpha_beta_occupation=zeros(norbitals,2);
           for i=1:norbitals
               temp_string = fgetl(fid1);
               temp_occupation(i)=str2num(temp_string);
               clear temp_string
           end    
           continue
       end
   end  
   if (length(lower(strtrim(data))) == length('<Molecular Orbital Spin Types>'))
       if (lower(strtrim(data)) == '<molecular orbital spin types>')
           for i=1:norbitals
               temp_string = fgetl(fid1);
               if (length(lower(strtrim(temp_string))) == 5)
                  if temp_occupation(i) < 1.01
                      alpha_beta_occupation(i,1)=temp_occupation(i);
                      alpha_beta_occupation(i,2)=0.0;
                      spin_available=1;
                  else
                      'Alpha occupation numbers greater than one are unphysical.'
                      'The program you used to generate the wfx file contains a bug.'
                      'Please regenerate the wfx file using the latest version of your quantum chemistry program.'
                      'If the error persists, please inform your quantum chemistry program vendor of the bug.'
                      'Program will terminate.'
                      flag=1;
                      break
                  end
               elseif (length(lower(strtrim(temp_string))) == 4)
                  alpha_beta_occupation(i,1)=0.0;
                  alpha_beta_occupation(i,2)=temp_occupation(i);
               else
                  alpha_beta_occupation(i,1)=temp_occupation(i)/2.0;
                  alpha_beta_occupation(i,2)=temp_occupation(i)/2.0;
               end 
           end 
           'Natural orbital occupation numbers for alpha (1st column) and beta (2nd column) electrons.'
           alpha_beta_occupation 
           'The total number of electrons in natural orbitals is'
           included_electrons=sum(sum(alpha_beta_occupation))
           continue
       end
   end
   if (length(lower(strtrim(data))) == length('<Molecular Orbital Primitive Coefficients>'))
       if (lower(strtrim(data)) == '<molecular orbital primitive coefficients>')
           orbital_coefficients=zeros(norbitals,nprimitives);
           for i=1:norbitals
               data = fgetl(fid1);
               if length(data) == 0
                   fgetl(fid1);
               end    
               clear data
               if str2num(fgetl(fid1)) == i
                   data = fgetl(fid1);
                   clear data
                   data = textscan(fid1, '%f',nprimitives);
                   temp_array=cell2mat(data);
                   for j=1:nprimitives
                       orbital_coefficients(i,j)=temp_array(j);
                   end    
                   clear data temp_array                   
               else
                   'Fatal error reading molecular orbital coefficients. Program will terminate.'
                   flag=1;
                   break
               end
           end
           continue
       end
   end   
               
                  
end
clear tline
if flag == 1
    break
end
% Check to make sure things were read in
if exist('basis_set_type') == 0
    'Basis set type could not be read from wfx file. Program will terminate.'
    flag=1;
    break
elseif exist('netcharge') == 0
    'Net charge could not be read from wfx file. Program will terminate.'
    flag=1;
    break
elseif exist('ncenters') == 0
    'Number of basis set centers could not be read from wfx file. Program will terminate.'
    flag=1;
    break
elseif exist('nprimitives') == 0
    'Number of basis set primitives could not be read from wfx file. Program will terminate.'
    flag=1;
    break 
elseif exist('norbitals') == 0
    'Number of natural orbitals could not be read from wfx file. Program will terminate.'
    flag=1;
    break  
elseif exist('alpha_beta_occupation') == 0
    'Natural orbital occupation numbers could not be read from wfx file. Program will terminate.'
    flag=1;
    break
elseif exist('orbital_coefficients') == 0
    'Natural orbital coefficients could not be read from wfx file. Program will terminate.'
    flag=1;
    break
end    
% Determine the number of atoms
natoms=0;
for i=1:ncenters
    if basis_set_centers(i,2) > 0.5
       natoms = natoms + 1;
    end
end
natoms
atomic_number=zeros(natoms,1);
core_electrons = zeros(natoms,1);
missing_core=zeros(1,natoms);
coords=zeros(natoms,3);
% Set the atomic coordinates, core charges,etc.
count=0;
for i=1:ncenters
    if basis_set_centers(i,2) > 0.5
       count = count + 1;
       atomic_number(count) = basis_set_centers(i,1);
       core_electrons(count) = basis_set_centers(i,1) - basis_set_centers(i,2);
       coords(count,1) = basis_set_centers(i,3);
       coords(count,2) = basis_set_centers(i,4);
       coords(count,3) = basis_set_centers(i,5);
       basis_set_centers(i,6) = count; % set the assigned atom in the last column of basis set centers
       if core_available == 1
          missing_core(count) = 0;
       else
          missing_core(count) = core_electrons(count);
       end
    end    
end
% Run consistency checks for the number of valence electrons
if abs(sum(atomic_number)- netcharge - sum(core_electrons) - sum(temp_occupation)) > 0.001
    'The number of valence electrons differs by more than 0.001e from the sum of orbital occupations.'
    'The quantum chemistry program you used to generate the wfx file contains a bug. Please contact the vendor of that program.'
    'Program will terminate.'
    flag == 1
    break
else
    'The number of valence electrons is consistent with the sum of orbital occupations.'
end
