function [settings,settings_enc] = read_cresis_xml(xml_fn)
% [settings,settings_enc] = read_cresis_xml(xml_fn)
%
% Function for reading in a cresis XML file (from NI system)
%
% See also read_ni_xml_directory

%% Read in XML file
xDoc = xmlread(xml_fn);
[settings,settings_enc] = read_ni_xml_object(xDoc);

settings.fn = xml_fn;
settings_enc.fn = xml_fn;

if isfield(settings,'sys')
  settings_enc.DDCZ20Ctrl = settings_enc.sys.DDCZ20Ctrl;
  settings_enc.DDSZ5FSetup = settings_enc.sys.DDSZ5FSetup;
  settings_enc.XMLZ20FileZ20Path = settings_enc.sys.XMLZ20FileZ20Path;
  settings_enc.xmlversion = settings_enc.sys.xmlversion;
  settings_enc = rmfield(settings_enc,'sys');
  settings.DDC_Ctrl = settings.sys.DDC_Ctrl;
  settings.DDS_Setup = settings.sys.DDS_Setup;
  settings.XML_File_Path = settings.sys.XML_File_Path;
  settings.xmlversion = settings.sys.xmlversion;
  settings= rmfield(settings,'sys');
end

%% Convert filename time stamp to datenum that is included with settings
[tmp xml_fn_name] = fileparts(xml_fn);
try
  time_stamp_idx = find(xml_fn_name == '_',1) + 1;
  year = str2double(xml_fn_name(time_stamp_idx + (0:3)));
  month = str2double(xml_fn_name(time_stamp_idx + 4 + (0:1)));
  day = str2double(xml_fn_name(time_stamp_idx + 6 + (0:1)));
  hour = str2double(xml_fn_name(time_stamp_idx + 9 + (0:1)));
  min = str2double(xml_fn_name(time_stamp_idx + 11 + (0:1)));
  sec = str2double(xml_fn_name(time_stamp_idx + 13 + (0:1)));
  %fprintf('  %s: %d\n', new_settings.fn, length(new_settings.Configuration.Waveforms));
  settings.datenum = datenum(year,month,day,hour,min,sec);
  settings_enc.datenum = datenum(year,month,day,hour,min,sec);
catch
  settings.datenum = NaN;
  settings_enc.datenum = NaN;
end

return;
