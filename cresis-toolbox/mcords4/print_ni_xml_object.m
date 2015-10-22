function print_ni_xml_object(obj, parent_string, idx)
% struct_val = print_ni_xml_object(obj, parent_string, idx)
%
% Function for printing a National Instruments (NI) XML object.
% Primarily used for debugging the XML object.
%
% Input:
%  obj = xDoc from xmlread (see example below)
%  parent_string = whatever you want to show up as the first part of the
%    node name (default is "root")
%  idx = index to start on (default is 1)
% Output:
%  struct containing parameters described in XML file
%
% Example:
%   xml_fn = get_filename('/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/','mcords4','','xml');
%   xDoc = xmlread(xml_fn);
%   print_ni_xml_object(xDoc,'xDoc',1);
%
% Author: John Paden
%
% See also: read_ni_xml_object.m, read_ni_xml_directory.m

if ~exist('parent_string','var')
  parent_string = 'root';
end

if ~exist('idx','var')
  idx = 1;
end

while length(obj)
  node_type = obj.getNodeType;
  try
    tag_name = char(obj.getTagName);
  catch
    tag_name = '';
  end
  try
    text_content = char(obj.getTextContent);
    text_content(text_content == 10) = ' ';
  catch
    text_content = '';
  end
  
  cur_name = sprintf('%s(%d)', parent_string, idx);
  
  % Print current node
  fprintf('%s: NodeType{%d} TagName{%s} TextContent{%s}\n', cur_name, node_type, tag_name, text_content);
  
  % Print any children nodes
  print_ni_xml_object(obj.item(0), [cur_name,'.item'], 0);
  
  % Iterate to next sibling
  obj = obj.getNextSibling;
  idx = idx + 1;
end

return
