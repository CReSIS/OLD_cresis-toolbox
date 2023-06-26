% script csv_good_loader
%
% Loads and plots the CSV files from:
%   ftp://data.cresis.ku.edu/data/rds/csv_good/
%
% Identifies a few bad segments/frames and removes these

data_dir = 'C:\tmp\csv_good';
fns = get_filenames(data_dir,'Browse','Greenland','csv');

Latitude = [];
Longitude = [];
Thickness = [];
Elevation = [];
Frame = {};
Surface = [];
Bottom = [];
Quality = [];
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  fid = fopen(fn);
  C = textscan(fid, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
  fclose(fid);
  [LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = deal(C{:});
  Latitude = cat(2,Latitude,LAT.');
  Longitude = cat(2,Longitude,LON.');
  Thickness = cat(2,Thickness,THICK.');
  Elevation = cat(2,Elevation,ELEVATION.');
  Frame = cat(2,Frame,FRAME.');
  Surface = cat(2,Surface,SURFACE.');
  Bottom = cat(2,Bottom,BOTTOM.');
  Quality = cat(2,Quality,QUALITY.');
end

Latitude = Latitude(1:2:end);
Longitude = Longitude(1:2:end);
Elevation = Elevation(1:2:end);
Thickness = Thickness(1:2:end);
Frame = Frame(1:2:end);

% fprintf('%d points\n', numel(Latitude))
% figure(1); clf;
% scatter(Longitude,Latitude,[],Thickness,'.');
% hbar = colorbar;

bad_mask = zeros(size(Latitude));
idxs = find(Latitude > 76.94 & Latitude < 77.04 & Longitude > -38.4 & Longitude < -37)
Frame(idxs)
idxs = find(Latitude > 75.11 & Latitude < 75.16 & Longitude > -30.96 & Longitude < -30.8)
Frame(idxs)
idxs = find(Latitude > 69.111 & Latitude < 69.112 & Longitude > -49.348 & Longitude < -49.344)
Frame(idxs)
idxs = find(Latitude > 60.85 & Latitude < 61 & Longitude > -46.8 & Longitude < -46.6)
unique(Frame(idxs))

idxs = find(Latitude > 61 & Latitude < 61.2 & Longitude > -47.2 & Longitude < -46.8)
unique(Frame(idxs))


% '2006061001003' '2006061001006' <-- Frames with data from the Qassimiut Lobe

bad_mask = bad_mask | strncmp(Frame,'2003050901018',length('2003050901018'));
bad_mask = bad_mask | strncmp(Frame,'1997052101',length('1997052101'));
% bad_mask = bad_mask | strncmp(Frame,'1999051901',length('1999051901')); % Does not appear to be bad
bad_mask = bad_mask | strncmp(Frame,'2007091001009',length('2007091001009'));

% Bad Jakobshavn:
% 2005051402015

fprintf('%d points\n', numel(Latitude(~bad_mask)))
figure(1); clf;
scatter(Longitude(~bad_mask),Latitude(~bad_mask),[],Thickness(~bad_mask),'.');
hbar = colorbar;

