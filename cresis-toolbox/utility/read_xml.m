function document = read_xml(xml_obj,plot_flag,depth,parent_name)
% document = read_xml(xml_obj,depth,parent_name,plot_flag)
%
% Parses a structure returned from Matlab's built in xmlread function.
% Can be used to read KML files for example.
%
% xml_obj = object returned from xmlread (DOMNODE)
% plot_flag = boolean flag indicating whether or not to plot node names
%   and values (default is 0)
% depth,parent_name = used internally by this command when it calls itself
%   recursively
%
% document = structure representing the xml_obj
%
% Examples:
%
% fn = 'C:\Users\radar\Documents\Travel\Greenland11\flightplans\cryoland.kml';
% xDoc = xmlread(fn);
% document = read_xml(xDoc)
% pos = textscan(document.kml{1}.Document{1}.Placemark{1}.LineString{1}.coordinates{1}.text{1}.node_val, ...
%   '%f%f%f','Delimiter',',');
% lat = pos{1};
% lon = pos{2};
%
% Authors: John Paden
%
% See also xmlread.m

if ~exist('plot_flag','var') || isempty(plot_flag)
  plot_flag = 0;
end
if ~exist('depth','var') || isempty(depth)
  depth = 0;
end
if ~exist('parent_name','var') || isempty(parent_name)
  parent_name = '';
end

% Print out current node
if ~isempty(xml_obj.getNodeName)
  node_name = xml_obj.getNodeName.toCharArray.';
else
  node_name = '';
end
node_name(node_name=='#') = '';
if ~isempty(xml_obj.getNodeValue)
  node_val = xml_obj.getNodeValue.toCharArray.';
else
  node_val = '';
end
if plot_flag
  if ~all(isspace(node_name)) && ~all(isspace(node_val))
    fprintf('%s%s: %s:\n', parent_name, node_name, node_val);
  elseif ~all(isspace(node_name))
    fprintf('%s%s: --:\n', parent_name, node_name);
  elseif ~all(isspace(node_val))
    fprintf('%s--: %s:\n', parent_name, node_val);
  end
end

document.node_name = node_name;
document.node_val = node_val;

% Traverse children
for idx = 0:xml_obj.getLength-1
  if ~isempty(xml_obj.item(idx))
    if depth == 0
      new_parent_name = [node_name];
    elseif isempty(parent_name)
      new_parent_name = [parent_name '.' '--'];
    else
      new_parent_name = [parent_name '.' node_name];
    end
    tmp = read_xml(xml_obj.item(idx), plot_flag, depth+1, new_parent_name);
    if ~isempty(tmp.node_name)
      if isfield(document,tmp.node_name)
        document.(tmp.node_name){end+1} = tmp;
      else
        document.(tmp.node_name){1} = tmp;
      end
    end
  end
end

return;


