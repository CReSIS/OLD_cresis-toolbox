function clip_and_resample_image(h_img,h_axes,max_size_MB)
% clip_and_resample_image(h_img,h_axes,max_size_MB)
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
% Example of running function when you do not already have the handles of
% interest for the active figure:
%   h_children = get(gca,'Children');
%   h_img = h_children(end);
%   clip_and_resample_image(h_img,gca,10);
%
% Author: John Paden

CData = get(h_img,'CData');
XData = get(h_img,'XData');
YData = get(h_img,'YData');

size_x = size(CData,2);
size_y = size(CData,1);

dx = (XData(end)-XData(1)) / (size_x-1);
dy = (YData(end)-YData(1)) / (size_y-1);

xlims_all = dx * size(CData,2);
ylims_all = dy * size(CData,1);

xlims = xlim;
ylims = ylim;

size_x_clip = diff(xlims);
size_y_clip = diff(ylims);

Nx = size_x_clip / dx;
Ny = size_y_clip / dy;

if strcmp(class(CData),'double')
  sample_size = 16;
else
  sample_size = 8;
end
scale_factor = sqrt(Nx*Ny/ (max_size_MB*2^20/sample_size/size(CData,3)));

if scale_factor > 1
  Nx = Nx / scale_factor;
  Ny = Ny / scale_factor;
end

dx_new = diff(xlims)/Nx;
dy_new = diff(ylims)/Ny;
new_x_axis = linspace(xlims(1)+0.5*dx_new,xlims(end)-0.5*dx_new,Nx);
new_y_axis = linspace(ylims(1)+0.5*dy_new,ylims(end)-0.5*dy_new,Ny);

XData = XData(1):dx:XData(end);
YData = YData(1):dy:YData(end);
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
