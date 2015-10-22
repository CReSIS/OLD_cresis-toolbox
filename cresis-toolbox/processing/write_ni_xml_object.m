function write_ni_xml_object(settings, fid, encoded, param)
% write_ni_xml_object(settings, fid, encoded, param)
%
% Function for writing a National Instruments (NI) XML object.
%
% Input:
%  settings = struct from read_ni_xml_object.m
%  fid = file identifier (e.g. from fopen)
%  encoded = boolean for whether or not to use ni_xml_name_decode on the
%    struct "settings" field names
%  param = optional structure
%   .array_list = cell array of field names that should be treated like
%     arrays even if they are of length 1. For example:
%     param.array_list = {'Waveforms'}
%   .enum_list = cell array of field names that should be treated like
%     enumerations. For example:
%     param.enum_list = {'DDCZ20sel'}
%
% Example:
%   xml_fn = 'C:\tmp\2013_Greenland_P3\mcords_fast time gain meas_20130403\DDS_20130403_174348.xml'
%   xDoc = xmlread(xml_fn);
%   [settings,settings_enc] = read_ni_xml_object(xDoc);
%   fid = fopen('C:\tmp\test.xml','w');
%   fprintf(fid,'<?xml version=''1.0'' standalone=''yes'' ?>\n');
%   fprintf(fid,'<LVData xmlns="http://www.ni.com/LVData">\n');
%   write_ni_xml_object(settings_enc,fid,true);
%   fprintf(fid,'</LVData>');
%   fclose(fid);
%
% Author: John Paden
%
% See also: print_ni_xml_object.m, ni_xml_name_decode.m,
%   ni_xml_name_encode.m, read_ni_xml_object.m, write_ni_xml_object.m

% Input checking
if ~exist('param','var')
  param = [];
end

if ~isfield(param,'array_list')
  param.array_list = {};
end

names = fieldnames(settings);

for name_idx = 1:length(names)
  if encoded
    name_dec = ni_xml_name_decode(names{name_idx});
  else
    name_dec = names{name_idx};
  end
  
  %% Handle '>' and '<' in the field names by encoding as '&gt;' and '&lt;'
  % (i.e. escape these characters)
  name_esc = '';
  out_idx = 0;
  for idx = 1:length(name_dec)
    if name_dec(idx) == '>'
      out_idx = out_idx + 1;
      name_esc(out_idx:out_idx+3) = '&gt;';
      out_idx = out_idx + 3;
    elseif name_dec(idx) == '<'
      out_idx = out_idx + 1;
      name_esc(out_idx:out_idx+3) = '&lt;';
      out_idx = out_idx + 3;
    else
      out_idx = out_idx + 1;
      name_esc(out_idx) = name_dec(idx);
    end
  end
  
  if isstruct(settings.(names{name_idx}))
    fn = fieldnames(settings.(names{name_idx}));
    if ~isempty(strmatch(names{name_idx},param.enum_list,'exact'))
      %% Ring and Enumerated Type Controls (EW) type
        fprintf(fid,'<EW>\n');
        fprintf(fid,'<Name>%s</Name>\n', name_esc);
        for choice_idx = 1:length(settings.(names{name_idx}).Choice)
          fprintf(fid,'<Choice>%s</Choice>\n', settings.(names{name_idx}).Choice{choice_idx});
        end
        if ischar(settings.(names{name_idx}).Val)
          fprintf(fid,'<Val>%s</Val>\n', settings.(names{name_idx}).Val);
        else
          fprintf(fid,'<Val>%d</Val>\n', settings.(names{name_idx}).Val);
        end
        fprintf(fid,'</EW>\n');
      
    else
      %% Structure: Recurse
      if (~isempty(strmatch(names{name_idx},param.array_list,'exact')) && ~isempty(settings.(names{name_idx}))) ...
          || length(settings.(names{name_idx})) > 1
        fprintf(fid,'<Array>\n');
        fprintf(fid,'<Name>%s</Name>\n', name_esc);
        fprintf(fid,'<Dimsize>%d</Dimsize>\n', length(settings.(names{name_idx})));
        is_array = true;
      else
        is_array = false;
      end
      for array_idx = 1:length(settings.(names{name_idx}))
        fprintf(fid,'<Cluster>\n');
        if ~is_array
          fprintf(fid,'<Name>%s</Name>\n', name_esc);
        else
          fprintf(fid,'<Name>%s</Name>\n', '');
        end
        fprintf(fid,'<NumElts>%d</NumElts>\n', length(fieldnames(settings.(names{name_idx}))));
        write_ni_xml_object(settings.(names{name_idx})(array_idx), fid, encoded,param);
        fprintf(fid,'</Cluster>\n');
      end
      if is_array
        fprintf(fid,'</Array>\n');
      end
    end
    
  elseif ischar(settings.(names{name_idx}))
    %% String
    fprintf(fid,'<%s>%s</%s>\n', name_esc, settings.(names{name_idx}), name_esc);
    
  elseif isnumeric(settings.(names{name_idx})) || islogical(settings.(names{name_idx}))
    %% Numeric or Logical
    int_type = true;
    if strcmpi(class(settings.(names{name_idx})), 'double')
      type_str = 'DBL';
      int_type = false;
    elseif strcmpi(class(settings.(names{name_idx})), 'single')
      type_str = 'SGL';
      int_type = false;
    elseif strcmpi(class(settings.(names{name_idx})), 'uint8')
      type_str = 'U8';
    elseif strcmpi(class(settings.(names{name_idx})), 'int8')
      type_str = 'I8';
    elseif strcmpi(class(settings.(names{name_idx})), 'uint16')
      type_str = 'U16';
    elseif strcmpi(class(settings.(names{name_idx})), 'int16')
      type_str = 'I16';
    elseif strcmpi(class(settings.(names{name_idx})), 'uint32')
      type_str = 'U32';
    elseif strcmpi(class(settings.(names{name_idx})), 'int32')
      type_str = 'I32';
    elseif strcmpi(class(settings.(names{name_idx})), 'logical')
      type_str = 'Boolean';
    else
      keyboard
    end
    if length(settings.(names{name_idx})) > 1
      fprintf(fid,'<Array>\n');
      fprintf(fid,'<Name>%s</Name>\n', name_esc);
      fprintf(fid,'<Dimsize>%d</Dimsize>\n', length(settings.(names{name_idx})));
      is_array = true;
    else
      is_array = false;
    end
    for array_idx = 1:length(settings.(names{name_idx}))
      fprintf(fid,'<%s>\n',type_str);
      if ~is_array
        fprintf(fid,'<Name>%s</Name>\n', name_esc);
      else
        fprintf(fid,'<Name>%s</Name>\n', '');
      end
      if int_type
        fprintf(fid,'<Val>%.0f</Val>\n', settings.(names{name_idx})(array_idx));
      else
        fprintf(fid,'<Val>%.16g</Val>\n', settings.(names{name_idx})(array_idx));
      end
      fprintf(fid,'</%s>\n',type_str);
    end
    if length(settings.(names{name_idx})) > 1
      fprintf(fid,'</Array>\n');
    end
  else
    %% Unhandled case
    keyboard
  end
end

% fprintf(fid,'\n')



return;
