function [pts_struc] = kml_read_shapefile(filename,outname,varargin)
% [pts_struc] = kml_read_shapefile(filename,outname,varargin)
%
% Matlab script to read KML file into a mapstruct structure array and 
% offers the option of writing to a Shapefile. 
% KML files must be of a uniform type: Point, LineString, LinearRing or Polygon. 
% The output, called 'pts_struc', contains 5-6 fields, Geometry, X, Y, altitude,
% 'Id' and 'name'. The final field is created and populated only if there is 
% exactly one <name> tag for each feature. This mapstruct array can be easily 
% displayed in Matlab using the 'geoshow' command. This output format allows 
% exporting to a Shapefile using the Matlab command, 'shapewrite', which is 
% initiated by entering an output filename into the function. 
%
% INPUTS:
%   map_structure = kml_read_shapefile('myKML.kml');
%   map_structure = kml_read_shapefile('myKML.kml','mySHP.shp');
%
% OUTPUTS:
%   The code outputs a mapstruct structure array if only the filename is 
%   provided by the user. The mapstruct is put into memory. A shapefile is 
%   also created if an output shapefile filename is provided. This file is  
%   not brought into memory, however.
%
% by Michael Toomey, University of California Santa Barbara, 
% Department of Geography, Santa Barbara, CA 93106
% mtoomey@geog.ucsb.edu
% last modified: June 1, 2010
%
% See also run_kml_read_shapefile.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%ARGUMENTS IN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% determine number of arguments and whether they are the appropriate format
if nargin == 1
    fnm =ischar(filename);
    if fnm == 0
        disp('Filename must be a string')
    end    
elseif nargin == 2
    fnm=ischar(filename);
    onm=ischar(outname);
    if fnm == 0
        disp('Filename must be a string')
    end    
    if onm == 0 
        disp('Output Shapefile name must be a string');
    end
elseif nargin >= 3
    disp('Too many arguments. Arguments should include input KML filename and optionally, output Shapefile name')
end

%%%%%%%%%%%%%%%%%OPEN FILE AND EXTRACT CONTENTS AS CELL ARRAY%%%%%%%%%%%%%%
% file input name and open file, scan text, then close file
fid =fopen(filename);
if fid == -1
    disp('Cannot find file')
end
text = textscan(fid, '%s');
fclose(fid);
% determine size of scanned KML file - each line is one subelement in the
% cell array
sizes = numel(text{1});
body = zeros(sizes,1);

%%%%%%%%%%%%%%%%%%%FINDING TAGS FOR ALL REPORTED FEATURES%%%%%%%%%%%%%%%%%%
% start FOR loop to go through all text statements in the cell array and
% find instances of reported coordinates as well as description of object
% geometry (e.g. point, line, polygon)
for i=1:sizes
    % find indices of cell array, 'text', containing the string,
    % '<coordinates>', assign to 'body1st', and lines containing the string,
    % '</coordinates>', assign to 'body2nd'
    ans0 = strfind(char(text{1}(i)),'<coordinates>');
    body1st(i) = isempty(ans0);
    ans1 = strfind(char(text{1}(i)),'</coordinates>');
    body2nd(i) = isempty(ans1);
    % determine geometry of the objects by looking for one of the three
    % common Geometry tags
    ans2 = strfind(char(text{1}(i)),'<Point>');
    ispoint(i) = isempty(ans2);
    ans3 = strfind(char(text{1}(i)),'<Polygon>');
    ispolygon(i) = isempty(ans3);
    ans4 = strfind(char(text{1}(i)),'<LineString>');
    isline(i) = isempty(ans4);
    ans5 = strfind(char(text{1}(i)),'<LinearRing>');
    islinering(i) = isempty(ans5);
    % here we search for the <Placemark> and <name> tags, if present, to
    % create a name field for the output Shapefile.
    ans6 = strfind(char(text{1}(i)),'<Placemark>');
    isplacemk(i) = isempty(ans6);
    ans7 = strfind(char(text{1}(i)),'<name>');
    isname(i) = isempty(ans7);
    ans8 = strfind(char(text{1}(i)),'</name>');
    isendname(i) = isempty(ans8);
end

% determine which geometry types were reported and assign proper geometry
% to variable, 'geom' which will be later used to assign to the mapstruct
pt_inst = find(ispoint==0);
poly_inst = find(ispolygon==0);
line_inst = find(isline==0);
linering_inst = find(islinering==0);
if numel(pt_inst) >= 1
    geom = 'Point';
elseif numel(poly_inst) >= 1
    geom = 'Polygon';
elseif numel(line_inst) >= 1
    geom = 'Line';    
elseif numel(linering_inst) >= 1
    geom = 'Line';    
end
% determine locations within cell array of reported coordinates
reported_coords1st = find(body1st == 0);
reported_coords2nd = find(body2nd == 0);

% determine locations of 'Placemark' and 'name'
% tags. Instances where name comes right after Placemark, is an indication
% of a reported title/name for each feature. 'intersect' is used to
% determine whether and where the name tag immediately follows the 
% placemark tag, indicating it is an object name and not just a folder name
% or something similar
placemk_inst = find(isplacemk==0);
name_inst = find(isname==0);
endname_inst = find(isendname==0);
[intseca,intsecb,intsecc]=intersect(placemk_inst,name_inst-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%CREATE MAPSTRUCT ARRAY%%%%%%%%%%%%%%%%%%%%%%%%%
% determine whether there are as many names available as there are features 
% to map out. if affirmative, 'names'=1 and we include 'name' field
% in the mapstruct structure array
if numel(reported_coords1st)==numel(intsecc)
    pts_struc = struct('Geometry',{},'X',{},'Y',{},'altitude',{},'Id',{},'name',{}); 
    names = 1;
else
    pts_struc = struct('Geometry',{},'X',{},'Y',{},'altitude',{},'Id',{}); 
    names = 0;
end

%%%%%%%%%%%%%%%%%EXTRACTING POINTS/NODES FROM KML TAGS%%%%%%%%%%%%%%%%%%%%%
% if all instances of '<coordinate>' and '</coordinate>' are on same line,
% is easier to deal with than if on separate lines. compare
% 'reported_coords1st' and 'reported_coords2nd' to see if equivalent  
equivtest=isequalwithequalnans(reported_coords2nd,reported_coords1st);
% if vectors are the same, assign these lines of text to untrim_coords
% because all information is found on the same line. if not, each point is 
% reported on a separate line and will be dealt with later 
if equivtest == 1
    untrim_coords = text{1}(reported_coords1st);
    % convert to character array
    untrim_coords = char(untrim_coords);
    % remove all instances of 'coordinates' from each string and then convert
    % to number for insertion into double array. 
    for i=1:numel(reported_coords1st)
        temp1 = strrep(untrim_coords(i,:),'<coordinates>','');
        temp2 = strrep(temp1,'</coordinates>','');
        % the first token extracted is the longitude - write into the
        % structure array
        [temp3,temp4] = strtok(temp2,','); 
        pts_struc(i).X = str2num(temp3);
        % the second token extracted is the latitude
        [temp5,temp6] = strtok(temp4,','); 
        pts_struc(i).Y = str2num(temp5);
        % the third token extracted is the altitude
        [temp7,temp8] = strtok(temp6,','); 
        pts_struc(i).altitude = str2num(temp7);
        % also insert the correct geometry into the structure
        pts_struc(i).Geometry = geom;
        % insert appropriate Id number
        pts_struc(i).Id = i;
        % if names=1 (i.e. there is a name tag for each object, then strip 
        % tags and insert appropriate text into 'name' field of pts_struc,
        % with some fancy character array flippage
        if names ==1
            nametxt=text{1}(name_inst(intsecc(i)):endname_inst(intsecc(i)));
            test=char(nametxt);ans=test';ans=rot90(ans(:));
            % strip tags
            ans=strrep(ans,'<name>','');
            pts_struc(i).name = strrep(ans,'</name>','');
        end
    end
% here we deal with the case if coordinates, with their tags are not all
% reported on the same line
elseif equivtest == 0
    % determine how many points in each coordinates cluster
    numpts = (reported_coords2nd-reported_coords1st)-1;
    for i=1:numel(numpts)
        temp1=zeros(numpts(i),3);
        temp2=char(text{1}(reported_coords1st(i)+1:reported_coords2nd(i)-1));
        for j=1:numpts(i)
            % the first token extracted is the longitude 
            [temp3,temp4] = strtok(temp2(j,:),','); 
            temp1(j,1) = str2num(temp3);
            % the second token extracted is the latitude
            [temp5,temp6] = strtok(temp4,','); 
            temp1(j,2)= str2num(temp5);
            % the third token extracted is the altitude
            [temp7,temp8] = strtok(temp6,','); 
            temp1(j,3) = str2num(temp7);
        end
        % load node coordinates for each object, stored in 'temp1' into the
        % corresponding place in the mapstruct
        pts_struc(i).X = temp1(:,1);
        pts_struc(i).Y = temp1(:,2);
        % you can uncomment this next line if you want to have an altitude
        % assigned to each node in a multipoint line or polygon feature,
        % but that does not aid the mapshow presentation and crashes the
        % shapewrite.m program
        %pts_struc(i).altitude = temp1(:,3);
        pts_struc(i).altitude = temp1(1,3);
        % also insert the correct geometry into the structure
        pts_struc(i).Geometry = geom;
        % insert appropriate Id number
        pts_struc(i).Id = i;
        % if names=1 (i.e. there is a name tag for each object, then strip 
        % tags and insert appropriate text into 'name' field of pts_struc,
        % with some fancy character array flippage
        if names ==1
            nametxt=text{1}(name_inst(intsecc(i)):endname_inst(intsecc(i)));
            test=char(nametxt);ans=test';ans=rot90(ans(:));
            % strip tags
            ans=strrep(ans,'<name>','');
            pts_struc(i).name = strrep(ans,'</name>','');
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%WRITING SHAPEFILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write Shapefile to disk if output name has been provided
if nargin ==2
    shapewrite(pts_struc,outname);
    comment = cat(2,'Shapefile ',outname,' written to disk');
    disp(comment)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

