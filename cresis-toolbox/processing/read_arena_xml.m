function out = read_arena_xml(xml_fn,print_flag)
% out = read_arena_xml(xml_fn,print_flag)
%
% out = read_arena_xml('c:\temp\config_GHSR.xml',0);
%
% Author: John Paden

if ~exist('print_flag','var')
  print_flag = 0;
end

obj = xmlread(xml_fn);

  out = read_arena_xml_recurse(obj,'DOC',print_flag);

end

function out = read_arena_xml_recurse(obj,xml_path,print_flag)

out = [];
name_list = {};
for idx = 1:obj.getLength()
  if obj.item(idx-1).getNodeType == 1
    name = obj.item(idx-1).getNodeName.toCharArray.';
    attrib = read_arena_xml_attrib(obj.item(idx-1));
    name_idx = sum(strcmp(name,name_list))+1;
    name_list{end+1} = name;
    if obj.item(idx-1).getLength > 1
      out.(name){name_idx} = read_arena_xml_recurse(obj.item(idx-1), ...
        sprintf('%s.%s{%d}',xml_path,name,name_idx),print_flag);
      out.(name){name_idx} ...
        = read_arena_xml_attrib_print(attrib,out.(name){name_idx}, ...
        sprintf('%s.%s{%d}',xml_path,name,name_idx),print_flag);
    else
      val = obj.item(idx-1).getTextContent.toCharArray.';
      if isempty(attrib)
        if print_flag
          fprintf('%s.%s = {''%s''};\n', xml_path, name, val);
        end
        out.(name){name_idx} = val;
      else
        if print_flag
          fprintf('%s.%s.name = ''%s'';\n', xml_path, name, val);
        end
        out.(name){name_idx}.name = val;
        out.(name){name_idx} = read_arena_xml_attrib_print(attrib,out.(name){name_idx}, ...
          [xml_path '.' name],print_flag);
      end
    end
  end
end

end

function attrib = read_arena_xml_attrib(obj)

attrib = [];
for attrib_idx = 1:obj.getAttributes.getLength
  attrib(attrib_idx).name = obj.getAttributes.item(attrib_idx-1).getName.toCharArray.';
  attrib(attrib_idx).value = obj.getAttributes.item(attrib_idx-1).getValue.toCharArray.';
end

end

function out = read_arena_xml_attrib_print(attrib,out,xml_path,print_flag)

for attrib_idx = 1:length(attrib)
  if print_flag
    fprintf('%s.%s = ''%s'';\n', xml_path, attrib(attrib_idx).name, attrib(attrib_idx).value);
  end
  out.(attrib(attrib_idx).name) = attrib(attrib_idx).value;
end

end
