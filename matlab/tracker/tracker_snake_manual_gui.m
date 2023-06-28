function [layer,dataPnts] = tracker_snake_manual_gui(A,dataPnts,startAxis)
%
% Interactive GUI wrapper for tracker_snake
%
% With mouse hovering over figure, you can use the following key commands:
% d = add layer point
% k = keyboard (debug)
% q = quit
% backspace = delete last layer point
% z = zoom mode (press any key to leave zoom mode)
% 1 = mouse button click (read off current pixel)
%
% Example:
% A = medfilt2(IMAGE,[3 3]);
% [layer,pnts] = testTrack(A);

global h1 h1_axis;
global pntbase;

pntbase.A = A;
pntbase.offset_abs = 0;

% Plot geo data
h1 = figure(1); clf;
% Ensure that special tools are not enabled before setting callbacks
zoom on;
zoom off;
% Unset callback properties
set(h1,'Interruptible','on');
set(h1,'WindowButtonUpFcn',[]);
set(h1,'WindowButtonMotionFcn',[]);
set(h1,'WindowKeyPressFcn',[]);

imagesc(pntbase.A);
h1_axis = gca;
colormap(jet(256));

% ---------------------------------------------------------------------
% Parse dataPnts input
if ~exist('dataPnts','var') || isempty(dataPnts)
  pntbase.dataPnts = [];
  layer = zeros(1,size(pntbase.A,2));
else
  pntbase.dataPnts = dataPnts;
  layer = tracker_snake(pntbase.A,pntbase.dataPnts,pntbase.offset_abs);
end
if ~exist('startAxis','var') || isempty(startAxis)
  pntbase.startAxis = [1 size(A,2) 1 size(A,1)];
else
  pntbase.startAxis = startAxis;
end
axis(pntbase.startAxis);

% ---------------------------------------------------------------------
% Create point handles

hold on;
pntbase.layerHandle = plot(layer,'k');
for ind = 1:length(pntbase.dataPnts)
  pntbase.dataPnts(ind).handle = plot(pntbase.dataPnts(ind).col, ...
    pntbase.dataPnts(ind).row,'kx');
end
hold off;

% ---------------------------------------------------------------------
% Set mouse pointer
set(h1,'Pointer','crosshair');

% Set callback properties
set(h1,'Interruptible','off');
set(h1,'WindowButtonUpFcn',@gui_windowbuttonupfcn)
set(h1,'WindowKeyPressFcn',@gui_windowkeypressfcn);

% Wait for quit signal
set(h1,'UserData',0);
waitfor(h1,'UserData');

layer = tracker_snake(pntbase.A,pntbase.dataPnts,pntbase.offset_abs);

return;

% =====================================================================
% WindowButtonUpFcn call back function
% =====================================================================
function gui_windowbuttonupfcn(src,event)

global h1 h1_axis;
global pntbase;

[x,y,but] = get_mouse_info(h1,h1_axis);
if but ~= 1
  return;
end

title(sprintf('1: x = %.2f, y = %.2f, but = %d\n', x, y, but));

return;

% =====================================================================
% WindowKeyPressFcn call back function
% =====================================================================
function gui_windowkeypressfcn(src,event)

if strcmpi(event.Key,'f1')
  fprintf('Manual tracker: \n');
  fprintf('  left mouse button click then CMD: add point;\n');
  fprintf('    s: snake point\n');
  fprintf('    t: threshold point\n');
  fprintf('    c: constant point\n');
  fprintf('  backspace: delete last point;\n');
  fprintf('  q: quit and continue;\n');
  fprintf('  k: debug tracker code;\n');
  return;
end

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~isempty(event.Character)

  global h1 h1_axis;
  global pntbase;
  dataPntsInd = length(pntbase.dataPnts);

  [x,y] = get_mouse_info(h1,h1_axis);

  if 0 % This code for debugging
    if ischar(event.Key)
      fprintf('x = %.2f, y = %.2f, key = %s\n', x, y, event.Key);
    else
      fprintf('x = %.2f, y = %.2f, key = %d\n', x, y, event.Key);
    end
    if ~isempty(event.Modifier)
      fprintf('  Modifiers ');
      for ind = 1:length(event.Modifier)
        fprintf('%s ', event.Modifier{ind});
      end
      fprintf('\n');
    end
  end

  % see event.Modifier for modifiers
  switch event.Character
    case 8 % Backspace
      if length(pntbase.dataPnts) > 0
        fprintf('Deleted last point\n');
        delete(pntbase.dataPnts(dataPntsInd).handle);
        dataPntsInd = dataPntsInd - 1;
        pntbase.dataPnts = pntbase.dataPnts(1:dataPntsInd);
      end
    case {'s','t','c'}
      % Finish chip command
      dataPntsInd = dataPntsInd + 1;
      pntbase.dataPnts(dataPntsInd).col = round(x);
      pntbase.dataPnts(dataPntsInd).row = round(y);
      pntbase.dataPnts(dataPntsInd).method = event.Character;
      fprintf('Added point\n');
      fprintf('  Row = %d Col = %d\n', pntbase.dataPnts(dataPntsInd).row(1), ...
        pntbase.dataPnts(dataPntsInd).col(1));
      hold on;
      pntbase.dataPnts(dataPntsInd).handle ...
        = plot(pntbase.dataPnts(dataPntsInd).col, ...
          pntbase.dataPnts(dataPntsInd).row,'kx');
      hold off;
    case 'k'
      % Code for testing simple tracker
      global A;
      global A_pnt;
      global A_axis;
      A = pntbase.A;
      A_pnt = pntbase.dataPnts;
      A_axis = pntbase.startAxis;
      keyboard
      error('Testing simple_tracker');
      % Copy and paste:
      global A;
      global A_pnt;
      global A_axis;
      surfBins = simple_tracker(A,A_pnt,A_axis);
    case 'z'
      zoom;
    case 'q'
      if length(pntbase.dataPnts) < 1
        fprintf('Need to specify at least one data point\n');
      else
        fprintf('Quitting\n');
        figure(h1);
        title('');
        % Set mouse pointer
        set(h1,'Pointer','arrow');

        % Set callback properties
        set(h1,'Interruptible','on');
        set(h1,'WindowButtonUpFcn',[]);
        set(h1,'WindowKeyPressFcn',[]);
        set(h1,'UserData',2);
      end
  end
  if event.Character == 't' || event.Character == 'c' || event.Character == 's' || event.Character == 8
    layer = tracker_snake(pntbase.A,pntbase.dataPnts,pntbase.offset_abs);
    delete(pntbase.layerHandle);
    hold on;
    pntbase.layerHandle = plot(layer,'k');
    hold off;
  end
end

return;


