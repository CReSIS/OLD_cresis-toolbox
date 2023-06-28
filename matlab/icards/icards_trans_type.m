function [datasize,datastring] = trans_type(datatype)
% function [datasize,datastring] = trans_type(datatype)
%
% Given the datatype, this function return datasize (number of bytes in one piece of data)
% and datastring (this is used by fread and fwrite)
 
switch datatype
case 1,	% Incoherent Data
   datasize = 2;
   datastring = 'uint16';
case 2,	% I channel Coherent Data
   datasize = 2;
   datastring = 'int16';
case 3,	% Q channel Coherent Data
   datasize = 2;
   datastring = 'int16';
case 4,	% GPS string
   datasize = 1;
   datastring = 'uint8';
case 5,	% Computer Time
   datasize = 1;
   datastring = 'uint8';
case 10, % Corrected Time
   datasize = 4;
   datastring = 'float32';
case 11, % Corrected Lat
   datasize = 4;
   datastring = 'float32';
case 12, % Corrected Lon
   datasize = 4;
   datastring = 'float32';
case 13, % Corrected Ellipsoidal Height
   datasize = 4;
   datastring = 'float32';
case 14, % Pitch
   datasize = 4;
   datastring = 'float32';
case 15, % Roll
   datasize = 4;
   datastring = 'float32';
case 20,	% Top curve
   datasize = 4;
   datastring = 'float32';
case 21,	% Bottom curve
   datasize = 4;
   datastring = 'float32';
case 22, % Laser Altimeter Height Measurement
   datasize = 4;          
   datastring = 'float32';
case 23, % Slope in the X direction
   datasize = 4;
   datastring = 'float32';
case 24, % Slope in the Y direction
   datasize = 4;
   datastring = 'float32';
end