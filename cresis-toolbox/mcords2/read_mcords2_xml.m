% xml_fn = '/mnt/samba-ro2/mcords/dds_20110317_062937.xml';
% xml_fn = '/mnt/samba-ro2/mcords/dds_20110317_070122.xml';
% xml_fn = '/mnt/samba-ro2/mcords/dds_20110317_070609.xml';
% xml_fn = '/mnt/samba-ro2/mcords/dds_20110317_071050.xml';
% xml_fn = '/mnt/samba-ro2/mcords/dds_20110317_071641.xml';
xml_fn = '/mnt/samba-ro2/DDS_20110326_123745.xml';
% xml_fn = '/mnt/samba-ro2/DDS_20110326_125430.xml';

xDoc = xmlread(xml_fn);

% xDoc.item(0).item(3)
% xDoc.item(0).item(3).item(7)
% xDoc.item(0).item(3).item(7).item(3)
% xDoc.item(0).item(3).item(7).item(3).getFirstChild
% xDoc.item(0).item(3).item(7).item(3).getFirstChild.getData

fprintf('=======================================================\n');

config = struct();
for idx = 0:xDoc.item(0).item(3).getLength-1
  node_name = xDoc.item(0).item(3).item(idx).getNodeName.toCharArray.';
  if node_name(1) ~= '#' && xDoc.item(0).item(3).item(idx).getLength >= 5
    fprintf('%d: %s\n', idx, node_name)
    if strcmpi(node_name,'U32') ...
      || strcmpi(node_name,'I32') ...
      || strcmpi(node_name,'U16') ...
      || strcmpi(node_name,'I16') ...
      || strcmpi(node_name,'U8') ...
      || strcmpi(node_name,'I8') ...
      || strcmpi(node_name,'DBL')
      field_name = xDoc.item(0).item(3).item(idx).item(1).item(0).getNodeValue.toCharArray.';
      field_val = str2double(xDoc.item(0).item(3).item(idx).item(3).item(0).getNodeValue.toCharArray.');
      field_name(field_name == ' ') = '_';
      fprintf('  %s = %d\n', field_name, field_val);
      config = setfield(config,field_name,field_val);
      
    elseif strcmpi(node_name,'ARRAY')
      if isempty(xDoc.item(0).item(3).item(idx).item(1).item(0))
        field_name = 'temp';
      else
        field_name = xDoc.item(0).item(3).item(idx).item(1).item(0).getNodeValue.toCharArray.';
      end
      field_name(field_name == ' ') = '_';
      first_parenthesis = strfind(field_name,'(');
      if ~isempty(first_parenthesis)
        field_name = field_name(1:first_parenthesis-2);
      end
      
      % Find all DimSize variables
      dim = [];
      for idx2 = 1:xDoc.item(0).item(3).item(idx).getLength
        if ~isempty(xDoc.item(0).item(3).item(idx).item(idx2))
          if strcmpi('DimSize', ...
              xDoc.item(0).item(3).item(idx).item(idx2).getNodeName.toCharArray.')
            dim(end+1) = str2double(xDoc.item(0).item(3).item(idx).item(idx2).item(0).getNodeValue.toCharArray.');
          end
        end
      end
      if length(dim) == 1
        dim(end+1) = 1;
      end
      fprintf('  Array %d by %d\n', dim(1), dim(2));
      
      % Collect all U? variables and assign to array
      A = zeros(dim(1), dim(2));
      A_idx = 0;
      for idx2 = 1:xDoc.item(0).item(3).item(idx).getLength
        if ~isempty(xDoc.item(0).item(3).item(idx).item(idx2))
          node_name = xDoc.item(0).item(3).item(idx).item(idx2).getNodeName.toCharArray.';
          if strcmpi(node_name,'U32') ...
            || strcmpi(node_name,'I32') ...
            || strcmpi(node_name,'U16') ...
            || strcmpi(node_name,'I16') ...
            || strcmpi(node_name,'U8') ...
            || strcmpi(node_name,'I8') ...
            || strcmpi(node_name,'DBL')
            A_idx = A_idx + 1;
            A(A_idx) = str2double(xDoc.item(0).item(3).item(idx).item(idx2).item(3).item(0).getNodeValue.toCharArray.');
          end
        end
      end
      fprintf('  %s = \n', field_name);
      A
      config = setfield(config,field_name,A);
    end
    
  end
end


return;

