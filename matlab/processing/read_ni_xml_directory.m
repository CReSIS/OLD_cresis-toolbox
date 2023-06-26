function [settings,settings_enc] = read_ni_xml_directory(path_dir,file_prefix,changes_only)
% [settings,settings_enc] = read_ni_xml_directory(path_dir,file_prefix,changes_only)
%
% Function for reading in a directory full of National Instruments
% (NI) XML files.
%
% Input:
%  path_dir = directory path
%  file_prefix = 'DDS', 'radar', or 'mcords4'
%  changes_only = logical which returns only settings that represent a
%    change from the last settings file read
% Output:
%  settings = struct array of settings including time stamp of when the
%    settings were first applied
%  settings_enc = struct array of settings using encoded filenames so that
%    write_ni_xml_object function can use it
%
% Example:
%   path_dir = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/';
%   settings = read_ni_xml_directory(path_dir,'mcords4');
%
%   path_dir = '/N/dc2/projects/cresis/2013_Greenland_P3/20130410/mcords/';
%   settings = read_ni_xml_directory(path_dir,'DDS');
%
% Author: John Paden
%
% See also: print_ni_xml_object.m, ni_xml_name_decode.m,
%   ni_xml_name_encode.m, read_ni_xml_object.m, write_ni_xml_object.m

xml_fns = get_filenames(path_dir,file_prefix,'','.xml');

first_time = true;
settings = [];
settings_enc = [];
for fn_idx = 1:length(xml_fns)
  xml_fn = xml_fns{fn_idx};
  try
    [new_settings,new_settings_enc] = read_cresis_xml(xml_fn);
    if first_time
      settings = new_settings;
      settings_enc = new_settings_enc;
      first_time = false;
    elseif ~changes_only || compare_structs(rmfield(new_settings,{'datenum','fn'}), rmfield(settings(end),{'datenum','fn'}))
      settings(end+1) = new_settings;
      settings_enc(end+1) = new_settings_enc;
    end
  catch ME
    warning('%s failed, skipping', xml_fn);
    ME.getReport
  end
  
%   xDoc = xmlread(xml_fn);
%   new_settings = read_ni_xml_object(xDoc);
%   new_settings.fn = xml_fn;
%   [tmp xml_fn_name] = fileparts(xml_fn);
%   time_stamp_idx = find(xml_fn_name == '_',1) + 1;
%   year = str2double(xml_fn_name(time_stamp_idx + (0:3)));
%   month = str2double(xml_fn_name(time_stamp_idx + 4 + (0:1)));
%   day = str2double(xml_fn_name(time_stamp_idx + 6 + (0:1)));
%   hour = str2double(xml_fn_name(time_stamp_idx + 9 + (0:1)));
%   min = str2double(xml_fn_name(time_stamp_idx + 11 + (0:1)));
%   sec = str2double(xml_fn_name(time_stamp_idx + 13 + (0:1)));
%   %fprintf('  %s: %d\n', new_settings.fn, length(new_settings.Configuration.Waveforms));
%   new_settings.datenum = datenum(year,month,day,hour,min,sec);
end

return;
