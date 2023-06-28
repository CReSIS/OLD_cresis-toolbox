% script oneudaq_gps_track

fn_dir = '/mnt/scratch2/output/accum/2011_Greenland_P3/CSARP_qlook/20110412_01/';
day_seg = '20110412_01';

file_nums = [290:298];

load_data_en = true;
redo_surface_en = true;

physical_constants;

for file_idx = 1:length(file_nums)
  file_num = file_nums(file_idx);
  fn_name = sprintf('Data_%s_%03d.mat',day_seg,file_num);
  fn = fullfile(fn_dir,fn_name);
  if file_idx == 1
    load(fn);
  else
    tmp = load(fn);
    if load_data_en
      Data = cat(2,Data,tmp.Data);
    end
    Surface = cat(2,Surface,tmp.Surface);
    Latitude = cat(2,Latitude,tmp.Latitude);
    Longitude = cat(2,Longitude,tmp.Longitude);
    Elevation = cat(2,Elevation,tmp.Elevation);
  end
end

if load_data_en && redo_surface_en
  if 0
    [tmp surf_bins] = max(Data);
  else
    surf_bins = zeros(1,size(Data,2));
    for rline = 1:size(Data,2)
      tmp = find(lp(Data(:,rline)) > -35,1);
      if isempty(tmp)
        surf_bins(rline) = NaN;
      else
        surf_bins(rline) = tmp;
      end
    end
  end
  Surface = interp1(1:length(Time), Time, surf_bins);
end

SurfaceD = detrend(Surface);
ElevationD = detrend(Elevation);

along_track = geodetic_to_along_track(Latitude,Longitude,Elevation);

figure(1); clf;
plot(along_track/1e3,SurfaceD*c/2);
hold on;
plot(along_track/1e3,ElevationD,'r');
hold off;
xlabel('Along-track (km)');
ylabel('Relative height (m)');
legend('Accum Altimeter','GPS Elevation');

if load_data_en
  figure(2); clf;
  imagesc([],Time*1e6,lp(Data));
  colormap(1-gray(256));
  hold on;
  plot(Surface*1e6);
  hold off;
end

return;
