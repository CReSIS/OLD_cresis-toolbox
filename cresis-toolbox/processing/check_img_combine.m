function mdata = check_img_combine(data_fn, rlines)
% mdata = check_img_combine(data_fn, rlines)
%
% Function for checking the combining of images.
%
% mdata =
% check_img_combine('~/Scratch/rds/2016_Greenland_Polar6/CSARP_qlook/20160413_17/Data_20160413_17_005');
% mdata = check_img_combine('~/Scratch/rds/2016_Greenland_Polar6/CSARP_standard/20160426_03/Data_20160426_03_001.mat');
%
% mdata = check_img_combine('/cresis/snfs1/dataproducts/ct_data/accum/2017_Greenland_P3/CSARP_standard/20170511_01/Data_20170511_01_072.mat');
%
% Author: John Paden

if ~exist('rlines','var')
  rlines = 1;
end

[data_fn_dir,data_fn_name] = fileparts(data_fn);

clear mdata;

mdata{1} = load(data_fn);

h_fig = figure; clf(h_fig);
h_axes = axes('parent',h_fig);
imagesc([], mdata{1}.Time*1e6, lp(mdata{1}.Data), 'Parent', h_axes);
xlabel(h_axes, 'Range line');
ylabel(h_axes, 'Time (us)');
grid(h_axes,'on');

done = false; img = 1;
while ~done
  if data_fn_name(6) == 'i'
    data_img_fn = data_fn_name;
    data_img_fn(10:11) = sprintf('%02d', img);
    data_img_fn = fullfile(data_fn_dir,[data_img_fn '.mat']);
  else
    data_img_fn = fullfile(data_fn_dir,[data_fn_name(1:5) sprintf('img_%02d_', img) data_fn_name(6:end) '.mat']);
  end
  if exist(data_img_fn,'file')
    mdata{img+1} = load(data_img_fn);
    img = img + 1;
  else
    done = true;
  end
end

h_fig = figure; clf(h_fig);
h_axes = axes('parent',h_fig);
for rline = rlines(:).'
  for img = 1:length(mdata)
    h_plot = plot(h_axes,mdata{img}.Time*1e6, lp(mdata{img}.Data(:,rline)));
    hold(h_axes,'on');
    if img == 1
      legend_str{img} = 'Combined';
      set(h_plot,'LineWidth',2);
    else
      legend_str{img} = sprintf('Img %d', img-1);
    end
  end
end
xlabel(h_axes,'Time (us)');
ylabel(h_axes,'Relative power (dB)');
grid(h_axes,'on');
legend(h_axes,legend_str);

return
