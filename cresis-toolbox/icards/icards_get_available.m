function [A] = get_available(filename)
% function [A] = get_available(filename)
%
% this function returns a matrix that describes the data in 'filename'
% returns an empty matrix if the file doesn't contain any data. (I can't think of 
% an instance where this would happen)

A = zeros(3,0);

fid=fopen (filename,'r');
if fid==-1
   disp('Could not open file');
   return
end

fseek(fid,64,'cof');

while 1
   [temp num] = fread(fid,3,'int32');
   if num ~= 3  % reached the end of the file
      fclose(fid);
      A = A';
      return
   end
   A=[A temp];
   fseek(fid,temp(2)*temp(3),'cof');
end
