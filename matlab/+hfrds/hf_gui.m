% Class hf_gui
%
% Class which talks to the HF sounder
%
% obj = hf_gui(param_fn);

classdef (HandleCompatible = true) hf_gui < handle
  properties
    h_fig % Figure handle
    h_axes % Axes handle
    h_image % Image handle
    h_gui % Structure of graphics handles
    
    serial_dev % matlab serial device
    param % Parameters (used for file loading and saving)
    test_mode
    
    param_fn % Filename to save lines to
  end
  
  methods
    function obj = hf_gui(param_fn)
      if ~exist('param_fn','var')
        param_fn = '';
      end
      
      obj.param_fn = param_fn;
      obj.test_mode = 0;
      
      if exist(param_fn,'file')
        obj.param = load(obj.param_fn);
      else
        obj.param.serial = 'COM8';
        obj.param.baud_rate = 115200;
        obj.param.input_buffer_size = 16777216;
        obj.param.wf = 1;
        obj.param.wf_delay = 1850;
        obj.param.pri = 65535;
        obj.param.record_start = 512;
        obj.param.record_stop = 2304;
        obj.param.presums = 512;
        obj.param.bit_shifts = 6;
        obj.param.decimation = 0;
        for row = 1:4
          obj.param.digital(row).inv = 1;
          obj.param.digital(row).start = 725;
          obj.param.digital(row).stop = 1032;
        end
        obj.param.vga = [0 48 56 56];
        obj.param.dac = 65280;
        obj.param.duration = 16384;
        obj.param.config = 2;
      end
      
      obj.h_fig = figure('Name','HF Sounder');
      h_fig_pos = get(obj.h_fig,'Position');
      set(obj.h_fig,'Position',h_fig_pos + [0 -185 250 185]);
     
      % h_gui.h_lpanel.{handle,...}
      % h_gui.h_rpanel.{handle,...}
      % h_gui.h_axes
      % h_gui.h_table
      % h_gui.h_lpanel_table
      % h_gui.h_rpanel_table
      %% Create widgets of main table
      obj.h_gui.h_lpanel.handle = uipanel('Parent',obj.h_fig);
      set(obj.h_gui.h_lpanel.handle,'HighlightColor',[0.8 0.8 0.8]);
      set(obj.h_gui.h_lpanel.handle,'ShadowColor',[0.6 0.6 0.6]);
      obj.h_gui.h_rpanel.handle = uipanel('Parent',obj.h_fig);
      set(obj.h_gui.h_rpanel.handle,'HighlightColor',[0.8 0.8 0.8]);
      set(obj.h_gui.h_rpanel.handle,'ShadowColor',[0.6 0.6 0.6]);
      
      %% Setup main table
      obj.h_gui.h_table.ui = obj.h_fig;
      obj.h_gui.h_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_table.false_height = NaN*zeros(30,30);
      row = 1; col = 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.h_lpanel.handle;
      obj.h_gui.h_table.width(row,col)     = 240;
      obj.h_gui.h_table.height(row,col)    = inf;
      obj.h_gui.h_table.false_height(row,col) = 0;
      col = col + 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.h_rpanel.handle;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = inf;
      obj.h_gui.h_table.false_height(row,col) = 0;
      obj.h_gui.h_table.width_margin(row,col) = 0;
      obj.h_gui.h_table.height_margin(row,col) = 0;
      obj.h_gui.h_table.width_margin ...
        = obj.h_gui.h_table.width_margin(1:row,1:col);
      obj.h_gui.h_table.height_margin ...
        = obj.h_gui.h_table.height_margin(1:row,1:col);
      obj.h_gui.h_table.false_height ...
        = obj.h_gui.h_table.false_height(1:row,1:col);
      clear row col
      table_draw(obj.h_gui.h_table);
      
      %% Setup right panel GUI objects
      obj.h_axes = axes('Parent',obj.h_gui.h_rpanel.handle);
      hold(obj.h_axes,'on');
      
      obj.h_gui.r_panel.h_plot = plot(1,1,'b-','Parent',obj.h_axes);
      set(obj.h_gui.r_panel.h_plot,'XData',[]);
      set(obj.h_gui.r_panel.h_plot,'YData',[]);

      obj.h_gui.h_rpanel.statusText = uicontrol('parent',obj.h_gui.h_rpanel.handle);
      set(obj.h_gui.h_rpanel.statusText,'Style','text');
      set(obj.h_gui.h_rpanel.statusText,'HorizontalAlignment','left');
      set(obj.h_gui.h_rpanel.statusText,'String','Ready.');

      %% Setup right panel table
      obj.h_gui.h_rpanel_table.ui = obj.h_gui.h_rpanel.handle;
      obj.h_gui.h_rpanel_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_rpanel_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_rpanel_table.false_height = NaN*zeros(30,30);
      row = 1; col = 1;
      obj.h_gui.h_rpanel_table.handles{row,col}   = obj.h_axes;
      obj.h_gui.h_rpanel_table.width(row,col)     = inf;
      obj.h_gui.h_rpanel_table.height(row,col)    = inf;
      obj.h_gui.h_rpanel_table.width_margin(row,col) = 30;
      obj.h_gui.h_rpanel_table.height_margin(row,col) = 30;
      obj.h_gui.h_rpanel_table.false_height(row,col) = 0;
      row = row + 1; col = 1;
      obj.h_gui.h_rpanel_table.handles{row,col}   = obj.h_gui.h_rpanel.statusText;
      obj.h_gui.h_rpanel_table.width(row,col)     = inf;
      obj.h_gui.h_rpanel_table.height(row,col)    = 20;
      obj.h_gui.h_rpanel_table.false_height(row,col) = 0;
      obj.h_gui.h_rpanel_table.width_margin(row,col) = 0;
      obj.h_gui.h_rpanel_table.height_margin(row,col) = 0;
      obj.h_gui.h_rpanel_table.width_margin ...
        = obj.h_gui.h_rpanel_table.width_margin(1:row,1:col);
      obj.h_gui.h_rpanel_table.height_margin ...
        = obj.h_gui.h_rpanel_table.height_margin(1:row,1:col);
      obj.h_gui.h_rpanel_table.false_height ...
        = obj.h_gui.h_rpanel_table.false_height(1:row,1:col);
      clear row col
      table_draw(obj.h_gui.h_rpanel_table);
      
      %% Setup left subpanel GUI objects
      % 	Style: [ {pushbutton} | togglebutton | radiobutton | 
      % checkbox | edit | text | slider | frame | listbox | popupmenu ]

      obj.h_gui.h_lpanel.serialText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.serialText,'Style','Text');
      set(obj.h_gui.h_lpanel.serialText,'String','Serial Port');

      obj.h_gui.h_lpanel.serialEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.serialEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.serialEdit,'String','');

      obj.h_gui.h_lpanel.wfText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.wfText,'Style','Text');
      set(obj.h_gui.h_lpanel.wfText,'String','Waveform');

      obj.h_gui.h_lpanel.wfEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.wfEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.wfEdit,'String','');

      obj.h_gui.h_lpanel.wf_delayText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.wf_delayText,'Style','Text');
      set(obj.h_gui.h_lpanel.wf_delayText,'String','Waveform Delay');

      obj.h_gui.h_lpanel.wf_delayEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.wf_delayEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.wf_delayEdit,'String','');

      obj.h_gui.h_lpanel.priText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.priText,'Style','Text');
      set(obj.h_gui.h_lpanel.priText,'String','PRI');

      obj.h_gui.h_lpanel.priEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.priEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.priEdit,'String','');

      obj.h_gui.h_lpanel.record_startText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.record_startText,'Style','Text');
      set(obj.h_gui.h_lpanel.record_startText,'String','Record Start/Stop');

      obj.h_gui.h_lpanel.record_startEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.record_startEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.record_startEdit,'String','');

      obj.h_gui.h_lpanel.record_stopEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.record_stopEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.record_stopEdit,'String','');

      obj.h_gui.h_lpanel.presumsText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.presumsText,'Style','Text');
      set(obj.h_gui.h_lpanel.presumsText,'String','Presums/Bit-shifts');

      obj.h_gui.h_lpanel.presumsEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.presumsEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.presumsEdit,'String','');

      obj.h_gui.h_lpanel.bit_shiftsEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.bit_shiftsEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.bit_shiftsEdit,'String','');

      obj.h_gui.h_lpanel.digitalText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.digitalText,'Style','Text');
      set(obj.h_gui.h_lpanel.digitalText,'String','Digital (Invert/Start/Stop)');

      for row = 1:4
        obj.h_gui.h_lpanel.digital(row,1) = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
        set(obj.h_gui.h_lpanel.digital(row,1),'Style','RadioButton');
        set(obj.h_gui.h_lpanel.digital(row,1),'String','Invert');
        set(obj.h_gui.h_lpanel.digital(row,1),'Value',0);
        for col = 2:3
          obj.h_gui.h_lpanel.digital(row,col) = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
          set(obj.h_gui.h_lpanel.digital(row,col),'Style','Edit');
          set(obj.h_gui.h_lpanel.digital(row,col),'String','');
        end
      end

      for row = 1:2
        obj.h_gui.h_lpanel.vga(row,1) = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
        set(obj.h_gui.h_lpanel.vga(row,1),'Style','Text');
        set(obj.h_gui.h_lpanel.vga(row,1),'String',sprintf('VGA [%d %d]',(row-1)*2+[1 2]));
        obj.h_gui.h_lpanel.vga(row,2) = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
        set(obj.h_gui.h_lpanel.vga(row,2),'Style','Edit');
        set(obj.h_gui.h_lpanel.vga(row,2),'String','');
        obj.h_gui.h_lpanel.vga(row,3) = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
        set(obj.h_gui.h_lpanel.vga(row,3),'Style','Edit');
        set(obj.h_gui.h_lpanel.vga(row,3),'String','');
      end

      obj.h_gui.h_lpanel.decimationText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.decimationText,'Style','Text');
      set(obj.h_gui.h_lpanel.decimationText,'String','Decimation');

      obj.h_gui.h_lpanel.decimationEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.decimationEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.decimationEdit,'String','');

      obj.h_gui.h_lpanel.dacText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.dacText,'Style','Text');
      set(obj.h_gui.h_lpanel.dacText,'String','DAC Gain');

      obj.h_gui.h_lpanel.dacEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.dacEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.dacEdit,'String','');

      obj.h_gui.h_lpanel.durationText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.durationText,'Style','Text');
      set(obj.h_gui.h_lpanel.durationText,'String','Record Duration (sec)');

      obj.h_gui.h_lpanel.durationEdit = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.durationEdit,'Style','Edit');
      set(obj.h_gui.h_lpanel.durationEdit,'String','');
      
      obj.h_gui.h_lpanel.configText = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.configText,'Style','Text');
      set(obj.h_gui.h_lpanel.configText,'String','Config Profile');

      obj.h_gui.h_lpanel.configPM = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.configPM,'Style','popupmenu');
      set(obj.h_gui.h_lpanel.configPM,'String',{'1','2','3','4','5','6','7','8'});
      set(obj.h_gui.h_lpanel.configPM,'Value',1);
      set(obj.h_gui.h_lpanel.configPM,'Min',1);
      
      obj.h_gui.h_lpanel.programPB = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.programPB,'Style','PushButton');
      set(obj.h_gui.h_lpanel.programPB,'String','Program');
      set(obj.h_gui.h_lpanel.programPB,'Callback',@obj.programPB_callback);
      
      obj.h_gui.h_lpanel.quickPB = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.quickPB,'Style','PushButton');
      set(obj.h_gui.h_lpanel.quickPB,'String','Quick Look');
      set(obj.h_gui.h_lpanel.quickPB,'Callback',@obj.quickPB_callback);
      
      obj.h_gui.h_lpanel.startPB = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.startPB,'Style','PushButton');
      set(obj.h_gui.h_lpanel.startPB,'String','Start Radar');
      set(obj.h_gui.h_lpanel.startPB,'Callback',@obj.startPB_callback);
      
      obj.h_gui.h_lpanel.stopPB = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.stopPB,'Style','PushButton');
      set(obj.h_gui.h_lpanel.stopPB,'String','Stop Radar');
      set(obj.h_gui.h_lpanel.stopPB,'Callback',@obj.stopPB_callback);
      
      obj.h_gui.h_lpanel.openPB = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.openPB,'Style','PushButton');
      set(obj.h_gui.h_lpanel.openPB,'String','Load Settings');
      set(obj.h_gui.h_lpanel.openPB,'Callback',@obj.openPB_callback);
      
      obj.h_gui.h_lpanel.savePB = uicontrol('Parent',obj.h_gui.h_lpanel.handle);
      set(obj.h_gui.h_lpanel.savePB,'Style','PushButton');
      set(obj.h_gui.h_lpanel.savePB,'String','Save Settings');
      set(obj.h_gui.h_lpanel.savePB,'Callback',@obj.savePB_callback);
      
      %% Setup left subpanel table
      obj.h_gui.h_lpanel_table.ui = obj.h_gui.h_lpanel.handle;
      obj.h_gui.h_lpanel_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_lpanel_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_lpanel_table.false_height = NaN*zeros(30,30);
      
      text_width = 80;
      
      row = 1; col = 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.serialText;
      obj.h_gui.h_lpanel_table.width(row,col)     = text_width;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.serialEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.programPB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.wfText;
      obj.h_gui.h_lpanel_table.width(row,col)     = text_width;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.wfEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.quickPB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.wf_delayText;
      obj.h_gui.h_lpanel_table.width(row,col)     = text_width;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.wf_delayEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.startPB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.priText;
      obj.h_gui.h_lpanel_table.width(row,col)     = text_width;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.priEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.stopPB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.record_startText;
      obj.h_gui.h_lpanel_table.width(row,col)     = text_width;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.record_startEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.record_stopEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.presumsText;
      obj.h_gui.h_lpanel_table.width(row,col)     = text_width;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.presumsEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.bit_shiftsEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.digitalText;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      for sub_row = 1:4
        for col = 1:3
          obj.h_gui.h_lpanel_table.handles{row+sub_row,col}   = obj.h_gui.h_lpanel.digital(sub_row,col);
          obj.h_gui.h_lpanel_table.width(row+sub_row,col)     = inf;
          obj.h_gui.h_lpanel_table.height(row+sub_row,col)    = 25;
        end
      end
      row = row + sub_row;
      
      for sub_row = 1:2
        for col = 1:3
          obj.h_gui.h_lpanel_table.handles{row+sub_row,col}   = obj.h_gui.h_lpanel.vga(sub_row,col);
          obj.h_gui.h_lpanel_table.height(row+sub_row,col)    = 25;
          if col == 1
            obj.h_gui.h_lpanel_table.width(row+sub_row,col)     = text_width;
            obj.h_gui.h_lpanel_table.false_height(row+sub_row,col) = 5;
          else
            obj.h_gui.h_lpanel_table.width(row+sub_row,col)     = inf;
          end
        end
      end
      row = row + sub_row;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.decimationText;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.decimationEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.dacText;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.dacEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.durationText;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.durationEdit;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.configText;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      obj.h_gui.h_lpanel_table.false_height(row,col) = 5;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.configPM;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;

      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.openPB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_lpanel.savePB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 25;

      obj.h_gui.h_lpanel_table.width_margin ...
        = obj.h_gui.h_lpanel_table.width_margin(1:row,1:3);
      obj.h_gui.h_lpanel_table.height_margin ...
        = obj.h_gui.h_lpanel_table.height_margin(1:row,1:3);
      obj.h_gui.h_lpanel_table.false_height ...
        = obj.h_gui.h_lpanel_table.false_height(1:row,1:3);
      clear row col
      table_draw(obj.h_gui.h_lpanel_table);
      
      %% Set up general handles
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      
      %% Load GUI with parameters (from param_fn or the defaults)
      obj.set_gui_parameters();
    end
    
    function delete(obj)
      % Delete the map figure handle
      try
        delete(obj.h_fig);
      end
      % Close the serial port/device
      try
        fclose(obj.serial_dev);
      end
    end
    
    function close_win(obj,h_obj,event)
      delete(obj);
    end
    
    function button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes)
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      % Check to make sure mouse clicked inside of obj.h_axes.handle
      %   Since extends the full y-length, just check to the right of minimum x
      set(obj.h_gui.h_rpanel.handle,'Units','normalized');
      uipanel_pos = get(obj.h_gui.h_rpanel.handle,'Position');
      set(obj.h_gui.h_rpanel.handle,'Units','Points');
      if mouse_pos(1) <= uipanel_pos(1)
        return
      end
      
      if but == 1
      elseif but == 3
      end
    end
    
    function button_down(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes)
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      % Check to make sure mouse clicked inside of obj.h_axes.handle
      %   Since extends the full y-length, just check to the right of minimum x
      set(obj.h_gui.h_rpanel.handle,'Units','normalized');
      uipanel_pos = get(obj.h_gui.h_rpanel.handle,'Position');
      set(obj.h_gui.h_rpanel.handle,'Units','Points');
      if mouse_pos(1) <= uipanel_pos(1)
        return
      end
      
      if but == 1
        %rbbox;
      elseif but == 3
      end
    end
    
    function button_scroll(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes)
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      % Check to make sure mouse clicked inside of obj.h_axes.handle
      %   Since extends the full y-length, just check to the right of minimum x
      set(obj.h_gui.h_rpanel.handle,'Units','normalized');
      uipanel_pos = get(obj.h_gui.h_rpanel.handle,'Position');
      set(obj.h_gui.h_rpanel.handle,'Units','Points');
      if mouse_pos(1) <= uipanel_pos(1)
        return
      end
      
      zooms = -1 + (event.VerticalScrollCount/2);
      
      cur_axis = [get(obj.h_axes,'Xlim') ...
        get(obj.h_axes,'YLim')];
      y_extent = cur_axis(4) - cur_axis(3);
      x_extent = cur_axis(2) - cur_axis(1);
      
      % Zoom so that the mouse pointer's position in the echogram does not change
      x_percent = (x-cur_axis(1))/x_extent;
      y_percent = (y-cur_axis(3))/y_extent;
      xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
      ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];

      xlim(xlims);
      ylim(ylims);
      
    end
    
    function programPB_callback(obj,h_obj,event)
      set(obj.h_gui.h_rpanel.statusText,'String','Programming...');
      obj.get_gui_parameters();
      obj.set_gui_parameters();
      obj.open_serial();
      obj.clear_serial();
      
      combined = obj.param.wf*2^16 + obj.param.wf_delay;
      obj.write(['P',0],'uint8');
      obj.write(combined,'uint32'); pause(0.25);
      
      obj.write(['P',1],'uint8');
      obj.write(obj.param.pri,'uint32'); pause(0.25);
      obj.write(['P',2],'uint8');
      obj.write(obj.param.record_start,'uint32'); pause(0.25);
      obj.write(['P',3],'uint8');
      obj.write(obj.param.record_stop,'uint32'); pause(0.25);
      obj.write(['P',4],'uint8');
      obj.write(obj.param.presums,'uint32'); pause(0.25);
      
      combined = obj.param.bit_shifts*2^4 + obj.param.decimation;
      obj.write(['P',5],'uint8');
      obj.write(combined,'uint32'); pause(0.25);

      for row = 1:4
        combined = obj.param.digital(1).inv*2^31 + obj.param.digital(row).start*2^16 + obj.param.digital(row).stop;
        obj.write(['P',5+row],'uint8');
        obj.write(combined,'uint32'); pause(0.25);
      end
      
      combined = sum(obj.param.vga .* 2.^[24 16 8 0]);
      obj.write(['P',10],'uint8');
      obj.write(combined,'uint32'); pause(0.25);
      
      obj.write(['P',11],'uint8');
      obj.write(obj.param.dac,'uint32'); pause(0.25);      
      obj.write(['P',12],'uint8');
      obj.write(obj.param.duration,'uint32'); pause(0.25);
      
      set(obj.h_gui.h_rpanel.statusText,'String','Programming... Done.');
    end
    
    function quickPB_callback(obj,h_obj,event)
      set(obj.h_gui.h_rpanel.statusText,'String','Getting quick look...');
      obj.get_gui_parameters();
      obj.open_serial();
      obj.clear_serial();
      
      obj.write(['Q'],'char'); pause(1);
        
        response = fread(obj.serial_dev,1,'uint8');
        if (response == 'Q')
          if 0
            start_read_tic = tic;
            time_out_sec = 10;
            samples_read = 0;
            while samples_read < obj.param.record_stop && toc(start_read_tic) < time_out_sec
              pause(1);
              samples_to_read = floor(obj.serial_dev.BytesAvailable/2);
              if samples_to_read > 0
                data = fread(obj.serial_dev,samples_to_read,'uint16');
                samples_read = samples_read + numel(data);
                set(obj.h_gui.h_rpanel.statusText,'String',sprintf('Getting quick look... %d samples read.', samples_read));
                drawnow;
              end
            end
            if samples_read < obj.param.record_stop
              set(obj.h_gui.h_rpanel.statusText,'String','Getting quick look... timed out.');
              return;
            end
          else
            pause(3);
            samples_to_read = floor(obj.serial_dev.BytesAvailable/2);
            data = fread(obj.serial_dev,samples_to_read,'uint16');
          end
        end
        hdr.fifo  = data(4);
        hdr.epri  = 2^16*data(5) + data(6);
        hdr.delay = data(7);
        hdr.gps   = char([floor(data(08)/256),mod(data(08),256), ...
          floor(data(09)/256),mod(data(09),256), ...
          floor(data(10)/256),mod(data(10),256)]);
        hdr.frac  = 2^16*data(11) + data(12);
        hdr.pri   = data(13);
        hdr.start = data(14);
        hdr.stop  = data(15);
        hdr.bshft = floor(data(16)/4096);
        hdr.dec   = mod(floor(data(16)/1024),4);
        hdr.pre   = mod(data(16),1024);
        wave      = mod(data(17:end)-mod((hdr.pre+1)*8192/(2^hdr.bshft),65536)+32768,65536)-32768;

        set(obj.h_gui.r_panel.h_plot,'XData',1:length(wave),'YData',wave);
      
      set(obj.h_gui.h_rpanel.statusText,'String','Getting quick look... Done.');
    end
    
    function startPB_callback(obj,h_obj,event)
      set(obj.h_gui.h_rpanel.statusText,'String','Starting radar...');
      obj.get_gui_parameters();
      obj.open_serial();
      obj.clear_serial();
      obj.write(['G','0'+obj.param.config],'char'); pause(3);
      set(obj.h_gui.h_rpanel.statusText,'String','Starting radar... Done.');
    end
    
    function stopPB_callback(obj,h_obj,event)
      set(obj.h_gui.h_rpanel.statusText,'String','Stopping radar...');
      obj.get_gui_parameters();
      obj.open_serial();
      obj.clear_serial();
      obj.write('X','char'); pause(3);
      set(obj.h_gui.h_rpanel.statusText,'String','Stopping radar... Done.');
    end
    
    function clear_serial(obj)
      bytes_to_read = get(obj.serial_dev,'BytesAvailable');
      if bytes_to_read > 0
        fread(obj.serial_dev,bytes_to_read,'uint8');
      end
    end
    
    function open_serial(obj)
      if obj.test_mode < 2
        if isempty(obj.serial_dev) || ~strcmpi(get(obj.serial_dev,'Port'),obj.param.serial)
          fprintf('Getting serial device %s handle\n', obj.param.serial);
          obj.serial_dev = serial(obj.param.serial,'BaudRate',obj.param.baud_rate, ...
            'InputBufferSize',obj.param.input_buffer_size);
        end
        if ~strcmpi(get(obj.serial_dev,'Status'),'open')
          fprintf('Opening serial device %s\n', obj.param.serial);
          fopen(obj.serial_dev);
        end
      end
    end
    
    function write(obj,write_buf,type)
      if obj.test_mode
        fprintf('Writing %s\n', type);
        if strcmpi(class(write_buf),'char')
          for idx = 1:numel(write_buf)
            fprintf('  %d: %s\n', idx, write_buf(idx));
          end
        else
          for idx = 1:numel(write_buf)
            fprintf('  %d: %g\n', idx, write_buf(idx));
          end
        end
        write_buf = double(write_buf);
        hex_str = dec2hex(write_buf);
        for row = 1:size(hex_str,1)
          fprintf('  %d: 0x%s (%d)\n', row, hex_str(row,:), write_buf(row));
        end
      else
        fwrite(obj.serial_dev,write_buf,type);
        if strcmpi(type,'uint8')
          read_buf = fread(obj.serial_dev,length(write_buf),type);
          if any(reshape(read_buf,[1 numel(read_buf)]) ~= write_buf)
            set(obj.h_gui.h_rpanel.statusText,'String','Error writing');
          end
        elseif strcmpi(type,'uint32')
          read_buf = fread(obj.serial_dev,length(write_buf),type);
          read_buf = dec2hex(read_buf,8);
          write_buf = dec2hex(write_buf,8);
          for row = 1:size(read_buf,1)
            read_buf(row,:) = read_buf(row,[7 8 5 6 3 4 1 2]);
            if any(read_buf(row,:) ~= write_buf(row,:))
              set(obj.h_gui.h_rpanel.statusText,'String','Error writing');
            end
          end
        end
      end
    end
    
    function savePB_callback(obj,h_obj,event)
      set(obj.h_gui.h_rpanel.statusText,'String','Saving...');
      [filename, pathname] = uiputfile( ...
        {'*.mat', 'All MATLAB Files (*.mat)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as', obj.param_fn);
      if isempty(filename) || filename(1) == 0
        return;
      end
      
      obj.param_fn = fullfile(pathname, filename);
      
      obj.get_gui_parameters();
      param = obj.param;
      save(obj.param_fn,'-struct','param');
      set(obj.h_gui.h_rpanel.statusText,'String',sprintf('Saving... %s', obj.param_fn));
    end
    
    function openPB_callback(obj,h_obj,event)
      set(obj.h_gui.h_rpanel.statusText,'String','Opening...');
      [filename, pathname] = uigetfile( ...
        {'*.mat','MAT-files (*.mat)'}, ...
        'Pick a file', ...
        obj.param_fn, 'MultiSelect', 'off');
      
      if isempty(filename) || filename(1) == 0
        return;
      end
      
      obj.param_fn = fullfile(pathname, filename);
              
      obj.param = load(obj.param_fn);
      
      obj.set_gui_parameters();
      set(obj.h_gui.h_rpanel.statusText,'String',sprintf('Opening... %s', obj.param_fn));
    end
    
    function set_gui_parameters(obj)
      set(obj.h_gui.h_lpanel.serialEdit,'String',sprintf('%s',obj.param.serial));
      set(obj.h_gui.h_lpanel.wfEdit,'String',sprintf('%d',obj.param.wf));
      set(obj.h_gui.h_lpanel.wf_delayEdit,'String',sprintf('%d',obj.param.wf_delay));
      set(obj.h_gui.h_lpanel.priEdit,'String',sprintf('%d',obj.param.pri));
      set(obj.h_gui.h_lpanel.record_startEdit,'String',sprintf('%d',obj.param.record_start));
      set(obj.h_gui.h_lpanel.record_stopEdit,'String',sprintf('%d',obj.param.record_stop));
      set(obj.h_gui.h_lpanel.presumsEdit,'String',sprintf('%d',obj.param.presums));
      set(obj.h_gui.h_lpanel.bit_shiftsEdit,'String',sprintf('%d',obj.param.bit_shifts));
      set(obj.h_gui.h_lpanel.decimationEdit,'String',sprintf('%d',obj.param.decimation));
      for row = 1:4
        set(obj.h_gui.h_lpanel.digital(row,1),'Value',obj.param.digital(row).inv);
        set(obj.h_gui.h_lpanel.digital(row,2),'String',sprintf('%d',obj.param.digital(row).start));
        set(obj.h_gui.h_lpanel.digital(row,3),'String',sprintf('%d',obj.param.digital(row).stop));
      end
      for row = 1:2
        for col = 2:3
          set(obj.h_gui.h_lpanel.vga(row,col),'String',sprintf('%d',obj.param.vga((row-1)*2 + col-1)));
        end
      end
      set(obj.h_gui.h_lpanel.dacEdit,'String',sprintf('%d',obj.param.dac));
      set(obj.h_gui.h_lpanel.durationEdit,'String',sprintf('%d',obj.param.duration));
      set(obj.h_gui.h_lpanel.configPM,'Value',obj.param.config);
    end
    
    function get_gui_parameters(obj)
      obj.param.serial = get(obj.h_gui.h_lpanel.serialEdit,'String');
      obj.param.wf = str2double(get(obj.h_gui.h_lpanel.wfEdit,'String'));
      obj.param.wf_delay = str2double(get(obj.h_gui.h_lpanel.wf_delayEdit,'String'));
      obj.param.pri = str2double(get(obj.h_gui.h_lpanel.priEdit,'String'));
      obj.param.record_start = str2double(get(obj.h_gui.h_lpanel.record_startEdit,'String'));
      obj.param.record_stop = str2double(get(obj.h_gui.h_lpanel.record_stopEdit,'String'));
      obj.param.presums = str2double(get(obj.h_gui.h_lpanel.presumsEdit,'String'));
      obj.param.bit_shifts = str2double(get(obj.h_gui.h_lpanel.bit_shiftsEdit,'String'));
      obj.param.decimation = str2double(get(obj.h_gui.h_lpanel.decimationEdit,'String'));
      for row = 1:4
        obj.param.digital(row).inv = get(obj.h_gui.h_lpanel.digital(row,1),'Value');
        obj.param.digital(row).start = str2double(get(obj.h_gui.h_lpanel.digital(row,2),'String'));
        obj.param.digital(row).stop = str2double(get(obj.h_gui.h_lpanel.digital(row,3),'String'));
      end
      for row = 1:2
        for col = 2:3
          tmp = str2double(get(obj.h_gui.h_lpanel.vga(row,col),'string'));
          if isempty(tmp)
            tmp = 0;
          end
          obj.param.vga((row-1)*2 + col-1) = tmp;
          if (obj.param.vga((row-1)*2 + col-1)) < 0
            obj.param.vga((row-1)*2 + col-1) = 0;
          elseif (obj.param.vga((row-1)*2 + col-1)) > 255
            obj.param.vga((row-1)*2 + col-1) = 255;
          end
        end
      end
      obj.param.dac = str2double(get(obj.h_gui.h_lpanel.dacEdit,'String'));
      obj.param.duration = str2double(get(obj.h_gui.h_lpanel.durationEdit,'String'));
      obj.param.config = get(obj.h_gui.h_lpanel.configPM,'Value');
    end
    
  end
end



