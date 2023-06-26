function [struct_val,struct_val_enc] = read_ni_xml_object(obj, parent_string)
% [struct_val,struct_val_enc] = read_ni_xml_object(obj, parent_string)
%
% Function for reading in National Instruments (NI) XML object.
%
% Input:
%  obj = xDoc.item(0) from xmlread (see example below)
%  parent_string = just for debugging recursion (do not use this input)
% Output:
%  struct_val: struct containing parameters described in XML file. A simple
%    field name mapping is applied to remove spaces and other characters
%    that are not allowed in struct field names.
%  struct_val_enc: struct with the same fields as struct_val, but using an
%    encoding for the field names so that the original field names (encoded)
%    can be used and written back to the file with write_ni_xml_object
%    The encoding is similar to URL encoding, but the special character is
%    "Z" since a-z and A-Z are the only fully unrestricted characters in
%    Matlab struct field names.  See ni_xml_name_encode/decode M-files.
%
% Example:
%   xml_fn = get_filename('/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/','mcords4','','xml');
%   xDoc = xmlread(xml_fn);
%   settings = read_ni_xml_object(xDoc);
%
% Author: John Paden
%
% See also: print_ni_xml_object.m, ni_xml_name_decode.m,
%   ni_xml_name_encode.m, read_ni_xml_object.m, write_ni_xml_object.m
%
% 1. xmlread returns a tree of nodes.
% 2. We parse the tree based on the type of each node
% POSSIBLE TYPES (Cluster, Array, DBL, U8, U32, I8, Boolean, etc with special type "Version")
%   Version
%    Special case where we just get the text contents
%   Cluster->Recursively creates a struct
%    item(1): Name
%    item(3): NumElts
%    item(5) and onward: Vals (values can be any type)
%   Array->Recursively creates an array (could be an array of structs or scalars)
%    item(1): Name
%    item(3): Dimsize
%    item(5) and onward: Vals (values can only be structs or scalars)
%   DBL, U16, etc->Creates a scalar
%    item(1): Name
%    item(3): Val (apply str2double)
% To traverse nodes:
%  obj.getNextSibling: Get next element
%  Only process nodes of type 3 (skip all other types)
% To pull fields out you need the name and contents... we have hard coded
% the position of the tag names, so we don't always pull the tag name
%  obj.getTagName: Cluster, Version, Array, DBL, U8, etc.
%    We don't grab the tag name when it is: Name, Val, NumElts, DimSize
%  obj.getTextContent: '12.0.1f3', 'Phase Offset', 0.000, 255, etc

if obj.getNodeType == 9
  %% This is the main document, skip to the first useful element
  obj = obj.item(0).item(0);
  parent_string = 'xDoc';
end

%% For debugging only
if ~exist('parent_string','var')
  parent_string = 'base';
end

%% Iterate through each element until we come to an empty object
first_time = 1;
while length(obj)
  if obj.getNodeType == 1 % This is the node type we want to read in
    
    %% Read in tag name (holds the "field type")
    tag_name = char(obj.getTagName);
    
    %% Parse each field type
    % Note: We have hard coded the position of "Name", "Val", "NumElts", "DimSize"
    %   which technically should be found by looking at the names of each node
    %   rather than assuming they occur at fixed positions.
    if strcmp(tag_name,'Version')
      %% Handle special "Version" field type
      struct_val.Version = char(obj.getTextContent);
      struct_val_enc.Version = char(obj.getTextContent);
      obj = obj.getNextSibling;
      % fprintf('%s\n', sprintf('%s.%s', parent_string, tag_name));  % DEBUG
      continue;
    else
      %% All other field types: read in the name
      name = char(obj.item(1).getTextContent);
      name_enc = ni_xml_name_encode(name);
      name(name==' ') = '_';
      name(name=='#') = '';
      name(name=='(') = '';
      name(name==')') = '';
      name(name=='-') = '';
      name(name=='>') = '';
    end
    
    % cur_string supports debugging
    if isempty(name)
      cur_string = sprintf('%s.%s', parent_string, tag_name);
    else
      cur_string = sprintf('%s.%s|%s', parent_string, name, tag_name);
    end
    %fprintf('%s\n', cur_string); % DEBUG
    
    if strcmp(tag_name,'Cluster')
      num_el = str2double(char(obj.item(3).getTextContent));
      [val,val_enc] = read_ni_xml_object(obj.item(5),cur_string);
    elseif strcmp(tag_name,'DBL')
      val = str2double(char(obj.item(3).getTextContent));
      val_enc = val;
    elseif strcmp(tag_name,'SGL')
      val = single(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'U8')
      val = uint8(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'U16')
      val = uint16(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'U32')
      val = uint32(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'I8')
      val = int8(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'I16')
      val = int16(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'I32')
      val = int32(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif strcmp(tag_name,'Boolean')
      val = logical(str2double(char(obj.item(3).getTextContent)));
      val_enc = val;
    elseif any(strcmp(tag_name,{'String','Path'}))
      val = {struct('type',tag_name,'values',[])};
      val{1}.values = {obj.item(3).getTextContent.char};
      val_enc = val;
    elseif strcmp(tag_name,'Array')
      dim_size = str2double(char(obj.item(3).getTextContent));
      [val,val_enc] = read_ni_xml_object(obj.item(5),cur_string);
    elseif strcmp(tag_name,'EW')
      % Ring and Enumerated Type Controls
      contents = textscan(get(obj,'TextContent'),'%s','Delimiter','\n');
      val = struct('Val', uint32(str2double(contents{1}{end})));
      val.Choice = contents{1}(3:end-1);
      val_enc = val;
    else
      warning('Unsupported type %s', tag_name);
      keyboard
    end
    
    %% Assign the value to this field
    if length(name) > 0
      % Struct assignment
      struct_val.(name) = val;
      struct_val_enc.(name_enc) = val_enc;
    else
      % Array assignment
      if first_time
        struct_val = val;
        struct_val_enc = val_enc;
        first_time = 0;
      else
        if iscell(val)
          struct_val{1}.values{end+1} = val{1}.values{1};
          struct_val_enc{1}.values{end+1} = val_enc{1}.values{1};
        else
          struct_val(end+1) = val;
          struct_val_enc(end+1) = val_enc;
        end
      end
    end
  end
  
  %% Iterate to next element
  obj = obj.getNextSibling;
end

return
