function clip_and_resample_image(h_img,h_axes,max_size_MB)
% downsize_image(h_img,h_axes,max_size_MB)
%
% Pass in an image and axes handle and the maximum memory size for the
% final product.
%
% Function which takes an image and resamples it using the current axis and
% specified resolution. Useful when saving figures with very large images
% where the resolution is much higher than needed or much of the image is
% clipped.
%
% Interpolation is carefully done so that no large memory copies are done.
% The only copy that is done is to copy two columns out of the original
% image and these are linearly interpolated.
%
% imresize would have also worked, but requires the image processing
% toolbox
%
% Author: John Paden

CData = get(h_img,'CData');
XData = get(h_img,'XData');
YData = get(h_img,'YData');
xlims = xlim(h_axes);
ylims = ylim(h_axes);

x_extent = max(xlims) - min(xlims);
y_extent = max(ylims) - min(ylims);

if strcmp(class(CData),'double')
  sample_size = 16;
else
  sample_size = 8;
end
resolution = sqrt(x_extent*y_extent / (max_size_MB*2^20/sample_size/size(CData,3)));

new_x_axis = min(xlims) : resolution : max(xlims)+resolution;
new_y_axis = (min(ylims) : resolution : max(ylims)+resolution);
XData = linspace(min(XData),max(XData),size(CData,2));
YData = linspace(min(YData),max(YData),size(CData,1)).';
if strcmpi(get(h_axes,'YDir'),'normal')
  YData = flipud(YData);
end

x_mask = XData > min(xlims)-x_extent*0.1 & XData < max(xlims)+x_extent*0.1;
y_mask = YData > min(ylims)-y_extent*0.1 & YData < max(ylims)+y_extent*0.1;
XData = XData(x_mask);
YData = YData(y_mask);
CData = CData(y_mask,x_mask,:);

class_fh = eval(sprintf('@%s',class(CData)));
CData_new = zeros(length(new_y_axis),length(new_x_axis),size(CData,3),class(CData));
for idx=1:size(CData_new,3)
  for col=1:size(CData_new,2)
    x_idx = find(XData < new_x_axis(col),1,'last');
    if isempty(x_idx)
      CData_new(:,col,idx) = NaN;
    else
      x_idxs = max(1,x_idx-1) : min(length(XData),x_idx+1);
      [ygrid,xgrid] = meshgrid(YData,XData(x_idxs));
      if sample_size == 16
        CData_new(:,col,idx) = class_fh(interp2(ygrid,xgrid,CData(:,x_idxs,idx).',new_y_axis,new_x_axis(col)));
      else
        CData_new(:,col,idx) = class_fh(interp2(ygrid,xgrid,single(CData(:,x_idxs,idx).'),new_y_axis,new_x_axis(col)));
      end
    end
  end
end

set(h_img,'XData',new_x_axis,'YData',new_y_axis,'CData',CData_new);

end
