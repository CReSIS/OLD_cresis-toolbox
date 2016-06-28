
function [valid_data_file]=icards_data_ignore_list(fns,full_dir);
% This function is for the convenience to ignore some wrong file of certain
% days-----Qi Shi
switch full_dir
  case 'Z:\ICARDS\2002\may24\'%to ignore "fiberdly" of 20020524
    valid_data_file=[2:205]; 
  case 'Z:\ICARDS\2002\may20\'%to ignore "seaice"   of 20020520
    valid_data_file=[1:43];
  case 'Z:\ICARDS\1997\may21\';%there is a "may21_98.043" in this folder
    valid_data_file=[1:212];
  case 'Z:\ICARDS\1997\may17\'%there is a "may17_98.011" in this folder
    valid_data_file=[1:103];
    
  case '/cresis/snfs1/data/ICARDS/2002/may24/'%to ignore "fiberdly" of 20020524(Linux)
    valid_data_file=[2:205]; 
  case '/cresis/snfs1/data/ICARDS/2002/may20/'%to ignore "seaice"   of 20020520(Linux)
    valid_data_file=[1:43];
  case '/cresis/snfs1/data/ICARDS/1997/may21/';%there is a "may21_98.043" in this folder(Linux)
    valid_data_file=[1:212];
  case '/cresis/snfs1/data/ICARDS/1997/may17/'%there is a "may17_98.011" in this folder(Linux)
    valid_data_file=[1:103];
    
  otherwise
    valid_data_file=[1:length(fns)];
end

end
  
