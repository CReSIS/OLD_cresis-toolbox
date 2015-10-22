% 07/06/2011
% Author: David Mattione
% 
% Function for reading data sets and paramters from paramter files to csv
% file
% 
% Column 1 - Date (YYYYMMDD)
% Column 2 - Data Set
% Column 3 - Processing Frames
% Column 4 - PRF
% Column 5 - Sampling Frequency
% Column 6 - Start Frequency
% Column 7 - Stop Frequency
% Column 8 - Fmult
% Column 9 - Flo
% Column 10 - Pulse Length
% Column 11 - ADC Bits
% Column 12 - VPP Scale
% Column 13 - Time Delay
% Column 14 - Time Correction
% Column 15 - Bandwidth
% Column 16 - Start Number
% Column 17 - GPS Filename
% Column 18 - Base Directory
% Column 19 - Folder Directory

param_snow_all_greenland_images_2010;




% header = {};
% header{1} = 'Date(YYYYMMDD)';
% header{2} = 'Segment';
% header{3} = 'Processing Frames';
% header{4} = 'PRF';
% header{5} = 'Fs';
% header{6} = 'F0';
% header{7} = 'F1';
% header{8} = 'FMult';
% header{9} = 'FLo';
% header{10} = 'TPD';
% header{11} = 'ADC Bits';
% header{12} = 'VPP Scale';
% header{13} = 'Time Delay';
% header{14} = 'Time_Corr';
% header{15} = 'Bandwidth';
% 
% 
% for header_idx = 1:length(header)
%   if header_idx == length(header)
%     fprintf(fid,'%s\n',header{header_idx});
%   elseif header_idx < length(header)
%     fprintf(fid,'%s,',header{header_idx});
%   end
% end

fid = fopen('/cresis/scratch2/mdce/dmattione/cresis-toolbox/fmcw/test_files/radar.csv','w');
fprintf(fid,'%s\n','Command');
for i1 = 1:length(param_files)  
  param = feval(param_files{i1}); 
  for i2 = 1:length(param.data_set)
    fprintf(fid,'%s,%i,,,,%s\n',datestr([param.year,param.month,param.day,0,0,0],'yyyymmdd'),param.data_set(i2)+1,param.flight_name);
  end
end

fprintf(fid,'\n%s\n','Vectors');
for i1 = 1:length(param_files)  
  param = feval(param_files{i1}); 
  for i2 = 1:length(param.data_set)
    fprintf(fid,'%s,%i,%s,%s,%s,,%s,%f\n',datestr([param.year,param.month,param.day,0,0,0],'yyyymmdd'),param.data_set(i2)+1,param.output_dir,...
      datestr([param.year,param.month,param.day,0,0,0],'yyyymmdd'),'data0',param.gps_file,param.time_corr(i2));
  end
end

fprintf(fid,'\n%s\n','Qlook');
for i1 = 1:length(param_files)  
  param = feval(param_files{i1}); 
  for i2 = 1:length(param.data_set)
    fprintf(fid,'%s,%i,%i,%i,,,,%i,,%i,%i,%i,,%i,,%i,%i,,\n',datestr([param.year,param.month,param.day,0,0,0],'yyyymmdd'),param.data_set(i2)+1,1,1,4,5,3,1,...
      param.nyquist_zone,2,0);
  end
end

fprintf(fid,'\n%s\n','Post');
for i1 = 1:length(param_files)  
  param = feval(param_files{i1}); 
  for i2 = 1:length(param.data_set)
    fprintf(fid,'%s,%i,%i,%i,,,,%i,,,,%s,%i\n',datestr([param.year,param.month,param.day,0,0,0],'yyyymmdd'),param.data_set(i2)+1,1,1,1,'jpg',90);
  end
end

fprintf(fid,'\n%s\n','Radar');
for i1 = 1:length(param_files)  
  param = feval(param_files{i1}); 
  for i2 = 1:length(param.data_set)
    fprintf(fid,'%s,%i,,%i,%i,,,%i,%i,%f,%i,%i,%f\n',datestr([param.year,param.month,param.day,0,0,0],'yyyymmdd'),param.data_set(i2)+1,param.prf(i2),param.sampling_freq(i2),...
      32,0,param.pulse_len(i2),14,2,4.5e-8);
  end   
end


 fclose(fid);
 
 
 
 
 
 
 