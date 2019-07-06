function music_quick_loader

% 1. Specify filename in fn
% 2. Specify start range bin to plot (rbin) and how many range bins(rbins) of
% leave blank.

fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_music.mat';
% fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_20140429_01_005_wf2_music.mat';
% fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_20140429_01_067_wf2_music.mat';
rbin = [];
rbins = [];

% fn = 'Data_img_03_20190416_01_052.mat';
% fn = 'Data_img_03_20190416_01_053.mat';
% mdata.rbin = 700;
% mdata.rbins = 500;

mdata = load(fn);

mdata.Nt = size(mdata.Tomo.img,1);
mdata.Nx = size(mdata.Tomo.img,3);
mdata.rline = 1;
if ~exist('rbin','var') || isempty(rbin)
  rbin = 1;
end
mdata.rbin = rbin;
if ~exist('rbins','var') || isempty(rbins)
  rbins = mdata.Nt;
end
mdata.rbins = rbins;

mdata.h_fig = figure;
mdata.h_axes = axes('parent',mdata.h_fig);
mdata.h_image = imagesc([],mdata.rbin+(0:mdata.rbins-1),10*log10(mdata.Tomo.img(mdata.rbin+(0:mdata.rbins-1),:,mdata.rline)), 'parent', mdata.h_axes);
xlabel(mdata.h_axes, 'Elevation angle bin');
ylabel(mdata.h_axes, 'Range bin');
set(mdata.h_fig,'WindowKeyPressFcn',@music_quick_loader_keyboard);
title(mdata.h_axes, sprintf('rline %d',mdata.rline));
set(mdata.h_fig,'UserData',mdata);

% Cycle through all range bins
fprintf('Run "clear global mdata" to remove this large variable from memory when finished.\n');
fprintf('Press "->" to move forward one range line.\n');
fprintf('Press "<-" to move backward one range line.\n');
fprintf('Press "up" to move up one range bin.\n');
fprintf('Press "down" to move down one range bin.\n');
fprintf('Press shift when moving to move ten bins or lines.\n');

end

%% music_quick_loader_keyboard
function music_quick_loader_keyboard(h_object,event)

if any(strcmpi(event.Modifier,'shift'))
  step = 10;
else
  step = 1;
end

if ~isempty(event.Key)
  mdata = h_object.UserData;
  
  % see event.Modifier for modifiers
  switch event.Key
    
    case 'uparrow' % Up arrow
      if mdata.rbin > step
        mdata.rbin = mdata.rbin - step;
      end
      
    case 'downarrow' % Down arrow
      if mdata.rbin < mdata.Nt-mdata.rbins+2-step
        mdata.rbin = mdata.rbin + step;
      end
      
    case 'rightarrow' % Right arrow
      if mdata.rline < mdata.Nx+1-step
        mdata.rline = mdata.rline + step;
      end
      
    case 'leftarrow' % Left arrow
      if mdata.rline > step
        mdata.rline = mdata.rline - step;
      end
      
    otherwise
      return;
  end
  set(mdata.h_image,'CData',10*log10(mdata.Tomo.img(mdata.rbin+(0:mdata.rbins-1),:,mdata.rline)))
  set(mdata.h_image,'YData',mdata.rbin+(0:mdata.rbins-1));
  ylim(mdata.rbin+[0 mdata.rbins-1]);
  title(mdata.h_axes, sprintf('rline %d',mdata.rline));
  set(h_object,'UserData',mdata);
end

end
