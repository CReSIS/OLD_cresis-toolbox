% script read_DMS_photogrammetry
%
% Example of how to read DMS photogrammetry data

combine = false;
if combine == true
  final_RGB = [];
  final_xaxis = [];
  final_yaxis = [];
end
std_vals = [];
num_vals = [];
hist_centers = -6:0.01:10;
n = zeros(size(hist_centers));
for line = 1:15
  geotiff_fns = get_filenames(fullfile('C:\tmp\DMS_Photogrammetry\03-23-2011\',sprintf('Line%d',line)),'','DEM','.tif');
  %geotiff_fns = get_filenames('C:\tmp\DMS_Photogrammetry\03-23-2011\Line5\','','DEM','.tif');
  
  for fn_idx = 1:length(geotiff_fns)
    %geotiff_fn = 'C:\tmp\DMS_Photogrammetry\03-23-2011\Line1\Line1_DEM_A_01.tif';
    geotiff_fn = geotiff_fns{fn_idx};
    fprintf('Line %d of %d File %d of %d: %s\n', line, 15, fn_idx, length(geotiff_fns), geotiff_fn);
    
    % Get the projection information
    proj = geotiffinfo(geotiff_fn);
    
    % Read the image
    [RGB, R, tmp] = geotiffread(geotiff_fn);
    RGB(RGB == -32767) = NaN;
    
    xaxis = R(3,1) + R(2,1)*(0:size(RGB,2)-1);
    yaxis = R(3,2) + R(1,2)*(0:size(RGB,2)-1);
    
    if combine == true
      if isempty(final_RGB)
        final_RGB = RGB;
        final_xaxis = xaxis;
        final_yaxis = yaxis;
        xaxis_step = R(2,1);
        yaxis_step = R(1,2);
      else
        return
        
        new_xaxis = [fliplr(final_xaxis(1):-xaxis_step:min([xaxis final_xaxis])) ...
          final_xaxis(1)+xaxis_step:xaxis_step:max([xaxis final_xaxis])];
        figure(1); clf;
        plot(new_xaxis,'r');
        hold on;
        plot(final_xaxis,'b--');
        plot(xaxis,'b--');
        hold off;
        
        if yaxis_step < 0
          new_yaxis = [fliplr(final_yaxis(1):-yaxis_step:max([yaxis final_yaxis])) ...
            final_yaxis(1)+yaxis_step:yaxis_step:min([yaxis final_yaxis])];
        else
          new_yaxis = [fliplr(final_yaxis(1):-yaxis_step:min([yaxis final_yaxis])) ...
            final_yaxis(1)+yaxis_step:yaxis_step:max([yaxis final_yaxis])];
        end
        
        min(xaxis,final_xaxis):xaxis_step:max(xaxis,final_xaxis)
        new_RGB = interp1
      end
      
      imagesc(final_xaxis,final_yaxis,final_RGB);
      colorbar;
    else
%       imagesc(xaxis,yaxis,RGB);
%       colorbar;
    end
    
    std_vals{line}(fn_idx) = nanstd(RGB(:));
    num_vals{line}(fn_idx) = nanstd(RGB(:));
    
    n_new = hist(RGB(~isnan(RGB)) - nanmean(RGB(:)),hist_centers);
    n = n + n_new;
    
  end
  plot(std_vals{line})
  drawnow;
end

std_vals_all = cell2mat(std_vals);
num_vals_all = cell2mat(num_vals);

return
pdf_mean = sum(hist_centers.*n) / sum(n)
pdf_std = sqrt(sum(abs(hist_centers).^2.*n) / sum(n))
threshold = 0.4;
pdf_mean = sum(hist_centers(abs(hist_centers)<threshold).*n(abs(hist_centers)<threshold)) / sum(n(abs(hist_centers)<threshold))
pdf_std = sqrt(sum(abs(hist_centers(abs(hist_centers)<threshold)).^2.*n(abs(hist_centers)<threshold)) / sum(n(abs(hist_centers)<threshold)))

% pdf_mean = 0;
% pdf_std = 0.15;

x_vals = linspace(-4,8,1001);
y_vals = 1 / (pdf_std*sqrt(2*pi)) * exp( -(x_vals - pdf_mean).^2/(2*pdf_std.^2) );

figure(1); clf;
h = bar(hist_centers, n/sum(n) / (hist_centers(2)-hist_centers(1)), 'hist')
hold on;
plot(x_vals,y_vals,'r')
hold off;
xlabel('RMS height (m)');
ylabel('Probability');
grid on;
xlim([-0.75 0.75]);

sqrt(mean(std_vals_all(~isnan(std_vals_all)).^2))
nanmean(std_vals_all)

% find(~isnan(mean(RGB,1)),1)
% find(~isnan(mean(RGB,1)),1,'last')
% find(~isnan(mean(RGB,2)),1)
% find(~isnan(mean(RGB,2)),1,'last')

% From first loaded image
RGB = RGB(214:950,743:end);

RGB_corr = xcorr2(RGB - nanmean(RGB(:)));
figure(2); clf;
imagesc(RGB_corr);
grid on;

corr_vals = RGB_corr(737,:);
corr_vals = corr_vals / max(corr_vals);
figure(3); clf;
plot(((1:1019)-510)*0.2, corr_vals);
x_vals = linspace(-10,10,1001);
hold on
y_vals = exp(-abs(x_vals)/1.2);
plot(x_vals,y_vals,'r')
hold off;
grid on;
xlim([-4 4]);
xlabel('Offset (m)');
ylabel('Correlation');



% set(obj.map_panel.h_image,'XData', xaxis, ...
%   'YData', yaxis, ...
%   'CData', A, ...
%   'Visible', 'on');
% set(obj.map_panel.h_axes, 'Xlim', sort(xaxis([1 end])), ...
%   'Ylim', sort(yaxis([1 end])), ...
%   'YDir', 'normal', ...
%   'Visible', 'on');

return

