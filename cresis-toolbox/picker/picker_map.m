function picker_map(cmd,param)
% picker_map(cmd,param)
%
% Runs commands related to the picker map.
%
% Called from picker.m and support functions.

global gCtrl;
global hui;

normal_color = [0 0 1];
selection_color = [1 0 0];
pick_color = [0 1 0];
view_color = [1 0 1];

% =================================================================
% =================================================================
% Update map selection (also called at startup)
% =================================================================
% =================================================================
if cmd == 1
  
  figure(hui.mapfig.handle); clf;
  cur_sel = get(hui.fig.ctrl_panel.mapsLB,'Value');
  
  if cur_sel == 1
    % --------------------------------------------------------------
    % No map is selected
    hui.mapfig.axes.handle = axes;
    hold on;

    % Create X,Y variables for each frame
    gCtrl.source.X = [];
    gCtrl.source.Y = [];
    gCtrl.source.F = [];
    for idx = 1:length(gCtrl.source.geo)
      gCtrl.source.geo{idx}.X = gCtrl.source.geo{idx}.Longitude;
      gCtrl.source.geo{idx}.Y = gCtrl.source.geo{idx}.Latitude;
    end
  else
    % --------------------------------------------------------------
    % Geotiff is selected
    
    if ~exist(gCtrl.geotiff.fns{cur_sel},'file')
      fprintf('%s does not exist\n', gCtrl.geotiff.fns{cur_sel});
      set(hui.fig.ctrl_panel.mapsLB,'Value',gCtrl.geotiff.cur);
      return;
    end
    [path name] = fileparts(gCtrl.geotiff.fns{cur_sel});
    gCtrl.tic = tic;
    fprintf('Loading and plotting the geotiff %s\n', name);
    % Obtain the projection structure.
    gCtrl.geotiff.proj = geotiffinfo(gCtrl.geotiff.fns{cur_sel});
    % Read the image
    [RGB, R, tmp] = geotiffread(gCtrl.geotiff.fns{cur_sel});
    if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
      RGB = uint8(RGB);
    end
    R = R/1e3;
    mapshow(RGB, R);
    hui.mapfig.axes.handle = gca;
    fprintf('  Done (%.1f sec)\n', toc(gCtrl.tic));
    
    % Create X,Y variables for each frame
    gCtrl.source.X = [];
    gCtrl.source.Y = [];
    gCtrl.source.F = [];
    for idx = 1:length(gCtrl.source.geo)
      [gCtrl.source.geo{idx}.X,gCtrl.source.geo{idx}.Y] ...
        = projfwd(gCtrl.geotiff.proj,gCtrl.source.geo{idx}.Latitude, ...
        gCtrl.source.geo{idx}.Longitude);
      gCtrl.source.geo{idx}.X = gCtrl.source.geo{idx}.X/1e3;
      gCtrl.source.geo{idx}.Y = gCtrl.source.geo{idx}.Y/1e3;
    end
  end
  
  % --------------------------------------------------------------
  % Plot flight lines and create quick search database
  hold on;
  for idx = 1:length(gCtrl.source.geo)
    gCtrl.source.geo{idx}.h = plot(gCtrl.source.geo{idx}.X,gCtrl.source.geo{idx}.Y);
    set(gCtrl.source.geo{idx}.h,'Color',normal_color)
    gCtrl.source.geo{idx}.h0 = plot(gCtrl.source.geo{idx}.X(1),gCtrl.source.geo{idx}.Y(1),'g.');
    set(gCtrl.source.geo{idx}.h0,'Color',normal_color)
    gCtrl.source.X = cat(2,gCtrl.source.X,gCtrl.source.geo{idx}.X);
    gCtrl.source.Y = cat(2,gCtrl.source.Y,gCtrl.source.geo{idx}.Y);
    gCtrl.source.F = cat(2,gCtrl.source.F, ...
      idx * ones(size(gCtrl.source.geo{idx}.X)) );
  end
    
  % Highlight current marker position
  hui.mapfig.pick_cursor.h = plot(gCtrl.source.geo{idx}.X(1),gCtrl.source.geo{idx}.Y(1),'kx','LineWidth',2,'MarkerSize',12);
  hui.mapfig.view_cursor.h = plot(gCtrl.source.geo{idx}.X(1),gCtrl.source.geo{idx}.Y(1),'kx','LineWidth',2,'MarkerSize',12);
  
  hold off;

  % Highlight current frame
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'Color',view_color)
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'Color',pick_color)
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'Color',selection_color)
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'LineWidth',4);
  
  % Set callback properties
  zoom on; zoom off;
  set(hui.mapfig.handle,'Interruptible','off');
  set(hui.mapfig.handle,'WindowButtonUpFcn',@picker_map_button)
  set(hui.mapfig.handle,'WindowButtonMotionFcn',[])
  set(hui.mapfig.handle,'WindowKeyPressFcn',@picker_map_key);

  gCtrl.geotiff.cur = get(hui.fig.ctrl_panel.mapsLB,'Value');

  if cur_sel == 1
    xlabel('Longitude (deg,E)');
    ylabel('Latitude (deg,N)');
  else
    xlabel('X (km)');
    ylabel('Y (km)');
  end

% =================================================================
% =================================================================
% Update current frame selection
% =================================================================
% =================================================================
elseif cmd == 2
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'Color',normal_color)
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'LineWidth',1);
  gCtrl.source.cur_sel = param;
  set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel);
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'Color',view_color)
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'Color',pick_color)
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'Color',selection_color)
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'LineWidth',4);

  menuString = gCtrl.source.src_disp{gCtrl.source.cur_sel};
  if gCtrl.source.cur_src > length(menuString)
    gCtrl.source.cur_src = 1;
    set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src);
  end
  set(hui.fig.ctrl_panel.sourceLB,'String',menuString);

% =================================================================
% =================================================================
% Update current frame pick
% =================================================================
% =================================================================
elseif cmd == 3
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'Color',normal_color)
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'LineWidth',1);
  gCtrl.source.cur_pick = param;
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'Color',view_color)
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'Color',pick_color)
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'Color',selection_color)
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'LineWidth',4);

% =================================================================
% =================================================================
% Update current frame view
% =================================================================
% =================================================================
elseif cmd == 4
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'Color',normal_color)
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'LineWidth',1);
  gCtrl.source.cur_view = param;
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'Color',view_color)
  set(gCtrl.source.geo{gCtrl.source.cur_view}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'Color',pick_color)
  set(gCtrl.source.geo{gCtrl.source.cur_pick}.h,'LineWidth',4);
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'Color',selection_color)
  set(gCtrl.source.geo{gCtrl.source.cur_sel}.h,'LineWidth',4);
  
end

return;
