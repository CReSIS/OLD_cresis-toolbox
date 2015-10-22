function kml_write_cresis(in_fn,out_fn,doc_name,type,sample_spacing)
% kml_write_cresis(in_fn,out_fn,doc_name,type,sample_spacing)
%
% A very limited KML writer for the specific purpose of reading
% in a CReSIS CSV file and writing a KML file with the format below.
%
% in_fn: input filename (cresis CSV file)
% out_fn: output filename (KML file)
% doc_name: document name to be placed inside KML file
% type: string indicating the type of KML file to create:
%   'vectors': creates points
%   'segment': creates lines
% sample_spacing: only used for type "segment"
%   2 element vector: [max_step_size num_points_in_line]
%     At least one must be finite.
%   max_step_size: can be inf in which case it is ignored, this sets the
%     maximum number of points that can be skipped (overrides
%     num_points_in_line).
%   num_points_in_line: set the step size to make this many points in the
%     line (useful for really long CSV lines where you don't need a lot of
%     detail). Set to inf to ignore.
%
% File Format For Type "Segment"
%
% <?xml version="1.0" encoding="utf-8"?>
% <kml xmlns="http://earth.google.com/kml/2.1">
%    <Document>
%       <name>20091016_01</name>
%       <Placemark>
%          <Snippet maxLines="0"> </Snippet>
%          <name>'20091016_01_001'</name>
%          <LineString>
%             <coordinates>-71.0869964418275,42.3500973174683,0.0
%             -71.0969964418275,42.3400973174683,0.0
%             -71.1069964418275,42.3200973174683,0.0
%             -71.1169964418275,42.3000973174683,0.0</coordinates>
%          </LineString>
%       </Placemark>
%       <Placemark>
%          <Snippet maxLines="0"> </Snippet>
%          <name>'20091016_01_002'</name>
%          <LineString>
%             <coordinates>-71.5869964418275,42.3500973174683,0.0
%             -71.5969964418275,42.3400973174683,0.0
%             -71.6069964418275,42.3200973174683,0.0
%             -71.6169964418275,42.3000973174683,0.0</coordinates>
%          </LineString>
%       </Placemark>
%    </Document>
% </kml>
%
% Author: John Paden

if strcmpi(type,'vectors')
  fid = fopen(in_fn);
  C = textscan(fid, '%f%f%f%f%s','delimiter',',','headerlines',1,'EndOfLine','\r\n');
  [LAT,LON,ELEVATION,TIME,FRAME] = C{:};
  fclose(fid);
else
  fid = fopen(in_fn);
  C = textscan(fid, '%f%f%f%f%f%s%f%f%f','delimiter',',','headerlines',1,'EndOfLine','\r\n');
  [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = C{:};
  fclose(fid);
end

out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
[fid,msg] = fopen(out_fn,'w');
if fid < 1
  fprintf('Could not open file %s\n', out_fn);
  error(msg);
end

fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid,'<kml xmlns="http://earth.google.com/kml/2.1">\n');
fprintf(fid,'   <Document>\n');
fprintf(fid,'      <name>%s</name>\n', doc_name);

fprintf(fid,'      <Style id="s_frame_normal">\n');
fprintf(fid,'        <LineStyle>\n');
fprintf(fid,'          <color>ff0000ff</color>\n');
fprintf(fid,'          <width>2</width>\n');
fprintf(fid,'        </LineStyle>\n');
fprintf(fid,'      </Style>\n');

fprintf(fid,'      <Style id="s_frame_highlight">\n');
fprintf(fid,'        <LineStyle>\n');
fprintf(fid,'          <color>ffff0000</color>\n');
fprintf(fid,'          <width>3</width>\n');
fprintf(fid,'        </LineStyle>\n');
fprintf(fid,'      </Style>\n');

fprintf(fid,'      <StyleMap id="m_frame">\n');
fprintf(fid,'        <Pair>\n');
fprintf(fid,'          <key>normal</key>\n');
fprintf(fid,'          <styleUrl>#s_frame_normal</styleUrl>\n');
fprintf(fid,'        </Pair>\n');
fprintf(fid,'        <Pair>\n');
fprintf(fid,'          <key>highlight</key>\n');
fprintf(fid,'          <styleUrl>#s_frame_highlight</styleUrl>\n');
fprintf(fid,'        </Pair>\n');
fprintf(fid,'      </StyleMap>\n');
  
if strcmpi(type,'mission')
  % Truncate to just day_seg so that placemarks are group by segment
  for idx = 1:length(FRAME)
    FRAME{idx} = FRAME{idx}(1:10);
  end
else
  % Do nothing, placemarks are group by frame
end

% Find all the unique placemarks
frm_ids = unique(FRAME);

% Create each placemark
for frm_idx = 1:length(frm_ids)
  
  fprintf(fid,'      <Placemark>\n');
  fprintf(fid,'         <Snippet maxLines="0"> </Snippet>\n');
	fprintf(fid,'         <styleUrl>#m_frame</styleUrl>\n');
  
  if strcmpi(type,'mission')
    fprintf(fid,'         <name>%s_%s</name>\n', ...
      frm_ids{frm_idx}(1:8), frm_ids{frm_idx}(9:10));
  else
    fprintf(fid,'         <name>%s_%s_%s</name>\n', ...
      frm_ids{frm_idx}(1:8), frm_ids{frm_idx}(9:10), frm_ids{frm_idx}(11:end));
  end

  if strcmpi(type,'vectors')
    fprintf(fid,'         <Point>\n');
    fprintf(fid,'            <coordinates>');
    
    pnt_idxs = strmatch(frm_ids{frm_idx},FRAME,'exact');
    
    fprintf(fid,'            %.13f,%.13f,%.3f\n', [LON(pnt_idxs(1)),LAT(pnt_idxs(1)),ELEVATION(pnt_idxs(1))].');
    
    fprintf(fid,'            </coordinates>\n');
    
    fprintf(fid,'            <description>');
    utc_time_datenum = datenum(str2double(FRAME{pnt_idxs(1)}(1:4)), ...
      str2double(FRAME{pnt_idxs(1)}(5:6)), ...
      str2double(FRAME{pnt_idxs(1)}(7:8)), ...
      0, 0, TIME(pnt_idxs(1)));
    fprintf(fid,'UTC Time: %s\n', datestr(utc_time_datenum, 'yyyymmdd HH:MM:SS'));
    fprintf(fid,'            </description>\n');
    
    fprintf(fid,'         </Point>\n');
    
  else
    fprintf(fid,'         <LineString>\n');
    fprintf(fid,'            <coordinates>');
    
    pnt_idxs = strmatch(frm_ids{frm_idx},FRAME);
    sample_spacing_actual = max(1,min(sample_spacing(1),round(length(pnt_idxs)/sample_spacing(2))));
    keep_idxs = [1:sample_spacing_actual:length(pnt_idxs)];
    if keep_idxs(end) ~= length(pnt_idxs)
      keep_idxs(end+1) = length(pnt_idxs);
    end
    pnt_idxs = pnt_idxs(keep_idxs);
    
    fprintf(fid,'            %.13f,%.13f,%.3f\n', [LON(pnt_idxs),LAT(pnt_idxs),THICK(pnt_idxs)].');
    
    fprintf(fid,'            </coordinates>\n');
    
    fprintf(fid,'         </LineString>\n');
  end
  fprintf(fid,'      </Placemark>\n');
end

fprintf(fid,'   </Document>\n');
fprintf(fid,'</kml>\n');

fclose(fid);

return;
