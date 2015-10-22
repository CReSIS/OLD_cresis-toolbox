function [A] = get_data(filename,datatype,offset,number)
% function [A] = get_data(filename,datatype,offset,number)
%
% given the filename and datatype, this function will return the all the data 
% associated with that datatype.  If you provide an offset,  it will return all 
% the data between the offset and the end.  If you provide the offset and the number, it 
% will start at the offset and return "number" records

if (nargin < 2)
   disp('get_data: you must specify a filename and datatype');
elseif (nargin == 2)
   offset = 0;
   number = inf;
elseif (nargin == 3)
   number = inf;
end

A = [];

if (round(offset) ~= offset) | (round(number) ~= number)
   disp('offset and number must be an integer');
   return
end 

avail = icards_get_available(filename);

if isempty(avail)
   return
end   

fid1=fopen (filename,'r');
if fid1==-1
   disp('Could not open file');
   return
end
   
if datatype==0  % return the header
   A=fread(fid1,2,'float32');
   A=[A;fread(fid1,14,'uint32')];
   fclose(fid1);
   return
end

% check to see if the data is available
datarow=find(avail(:,1)==datatype);
if isempty(datarow)
   disp('Data is not available');
   fclose(fid1);
   return;
end

%jump past the header
fseek(fid1,64,'cof');    % 16 - 4 byte 

[datasize datastring] = icards_trans_type(datatype); % find the datasize and datastring corresponding to datatype

for x=1:datarow-1,
   fseek(fid1,12,'cof'); %move past the 3 long secondary header
   fseek(fid1,avail(x,2)*avail(x,3),'cof'); %move past data set
end

fseek(fid1,12,'cof'); %move past the 3 long secondary header
num_of_records = avail(datarow,3);
if offset > num_of_records
   disp('get_data: offset is larger than number of records');
   return
end

fseek(fid1,avail(datarow,2)*offset,'cof'); %move to the offset point
num_of_records = num_of_records-offset;
if number > num_of_records
   if (isfinite(number))  % if user specified the number and it goes beyond the number of record contained in the data
      disp('get_data: there were not that many records; returning what I have');
   end
   number = num_of_records;
end

A = fread(fid1,[avail(datarow,2)/datasize number],datastring);

fclose(fid1);