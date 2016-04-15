function mdata = check_img_combine(data_fn, rlines)
% mdata = check_img_combine(data_fn, rlines)
%
% Function for checking the combining of images.
%
% mdata = check_img_combine('~/Scratch/rds/2016_Greenland_Polar6/CSARP_qlook/20160413_17/Data_20160413_17_005');
%
% Author: John Paden

if ~exist('rlines','var')
  rlines = 1;
end

[data_fn_dir,data_fn_name] = fileparts(data_fn);

clear mdata;

mdata{1} = load(data_fn);

figure(1); clf;
imagesc([], mdata{1}.Time, lp(mdata{1}.Data));
xlabel('Range line');
ylabel('Time (us)');
grid on;

done = false; img = 1;
while ~done
  data_img_fn = fullfile(data_fn_dir,[data_fn_name(1:5) sprintf('img_%02d_', img) data_fn_name(6:end) '.mat']);
  if exist(data_img_fn,'file')
    mdata{img+1} = load(data_img_fn);
    img = img + 1;
  else
    done = true;
  end
end

figure(2); clf;
for rline = rlines(:).'
  for img = 1:length(mdata)
    plot(mdata{img}.Time*1e6, lp(mdata{img}.Data(:,rline)));
    hold on;
    if img == 1
      legend_str{img} = 'Combined';
    else
      legend_str{img} = sprintf('Img %d', img-1);
    end
  end
end
xlabel('Time (us)');
ylabel('Relative power (dB)');
grid on;
legend(legend_str);

return
