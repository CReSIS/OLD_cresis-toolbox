function create_ui_crandall(obj,xpos,ypos)

% create_ui_crandall(obj,xpos,ypos)
%
% Creates components for the snake param window's UI when the crandall pick
% tool is selected. Plots the window at (xpos,ypos).
%

set(obj.h_fig,'visible','off');

figure(obj.h_fig);
clf;
obj.table = [];

%==========================================================================
% top panel
obj.top_panel.handle = uipanel('Parent',obj.h_fig);
set(obj.top_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.top_panel.handle,'ShadowColor',[0.6 0.6 0.6]);
%set(obj.top_panel.handle,'visible','off');
% bottom panel
obj.bottom_panel.handle= uipanel('Parent',obj.h_fig);
set(obj.bottom_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(obj.bottom_panel.handle,'ShadowColor',[0.6 0.6 0.6]);
%set(obj.bottom_panel.handle,'visible','off');

%--------------------------------------
% table
% set default position (changed when window accessed)
set(obj.h_fig,'Units','Pixels');
set(obj.h_fig,'Position',[xpos ypos obj.w obj.crandall_h]);

obj.table.ui=obj.h_fig;

row = 1; col = 1;
obj.table.handles{row,col}   = obj.top_panel.handle;
obj.table.width(row,col)     = inf;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col) = 0;
obj.table.height_margin(row,col) = 0;

row = 2; col = 1;
obj.table.handles{row,col}   = obj.bottom_panel.handle;
obj.table.width(row,col)     = inf;
obj.table.height(row,col)    = inf;
obj.table.width_margin(row,col) = 0;
obj.table.height_margin(row,col) = 0;

clear row col
table_draw(obj.table);

%============================================================================================
% top panel table contents

%----snake tool list box
obj.top_panel.tool_PM = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.tool_PM,'Style','popupmenu');
set(obj.top_panel.tool_PM,'String',{'basic','crandall','panton'});
set(obj.top_panel.tool_PM,'Value',2)
set(obj.top_panel.tool_PM,'Callback',@obj.toolPM_callback);

%-----mode label
obj.top_panel.mode_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.mode_label,'Style','text');
set(obj.top_panel.mode_label,'String','Mode');

%----insert range
obj.top_panel.insert_range_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.insert_range_label,'Style','text');
set(obj.top_panel.insert_range_label,'String','Max point range:');
%----param1 box
obj.top_panel.insert_range_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.insert_range_TE,'Style','edit');
set(obj.top_panel.insert_range_TE,'String',obj.in_rng_sv);
%----param2 name
obj.top_panel.snake_range_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.snake_range_label,'Style','text');
set(obj.top_panel.snake_range_label,'String','Snake range:');
%----param2 box
obj.top_panel.snake_range_TE = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.snake_range_TE,'Style','edit');
set(obj.top_panel.snake_range_TE,'String',obj.sn_rng_sv);
%----reinterp mode enable label
obj.top_panel.reinterp_mode_label = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.reinterp_mode_label,'Style','text');
set(obj.top_panel.reinterp_mode_label,'String','Reinterpolation (beta):');
%----reinterp mode enable cbox
obj.top_panel.reinterp_mode_cbox = uicontrol('Parent',obj.top_panel.handle);
set(obj.top_panel.reinterp_mode_cbox,'Style','checkbox');


%---------------------------------------------------------------------------------------------
% set up top panel table
obj.top_panel.table.ui=obj.top_panel.handle;
obj.top_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.top_panel.table.height_margin = NaN*zeros(30,30);
obj.top_panel.table.false_width = NaN*zeros(30,30);
obj.top_panel.table.false_height = NaN*zeros(30,30);
obj.top_panel.table.offset = [0 0];

row = 1; col = 1; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.mode_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 1; col = 2;
obj.top_panel.table.handles{row,col}   = obj.top_panel.tool_PM;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 2; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.insert_range_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 2; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.insert_range_TE;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 3; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.snake_range_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 3; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.snake_range_TE;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 4; col = 1;
obj.top_panel.table.handles{row,col}   = obj.top_panel.reinterp_mode_label;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

row = 4; col = 2; 
obj.top_panel.table.handles{row,col}   = obj.top_panel.reinterp_mode_cbox;
obj.top_panel.table.width(row,col)     = inf;
obj.top_panel.table.height(row,col)    = inf;
obj.top_panel.table.width_margin(row,col)= 1.5;
obj.top_panel.table.height_margin(row,col)=1.5;

clear row col

% Draw table
table_draw(obj.top_panel.table);

%=============================================================================================
% bottom_panel table contents
%---------------------------------------------------------------------------------------------
%---- top smooth 
obj.bottom_panel.topSmooth = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.topSmooth,'Style','text');
set(obj.bottom_panel.topSmooth,'String','Top smooth:');
%---- bottom smooth 
obj.bottom_panel.bottomSmooth = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.bottomSmooth,'Style','text');
set(obj.bottom_panel.bottomSmooth,'String','Bottom smooth:');

%---- top smooth slider
obj.bottom_panel.topSmoothSlider = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.topSmoothSlider,'Style','slider','Min',0,'Max',1,...
          'Value',obj.top_sm_sv,'SliderStep',[0.05 0.2]);
%---- bottom smooth slider
obj.bottom_panel.bottomSmoothSlider = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.bottomSmoothSlider,'Style','slider','Min',0,'Max',1,...
          'Value',obj.bot_sm_sv,'SliderStep',[0.05 0.2]);      
      
%---- top peakRatio 
obj.bottom_panel.topPeakRatio = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.topPeakRatio,'Style','text');
set(obj.bottom_panel.topPeakRatio,'String','Top peak ratio:');
%---- bottom PeakRatio 
obj.bottom_panel.bottomPeakRatio = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.bottomPeakRatio,'Style','text');
set(obj.bottom_panel.bottomPeakRatio,'String','Bottom peak ratio:');

%---- top PeakRatio slider
obj.bottom_panel.topPeakRatioSlider = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.topPeakRatioSlider,'Style','slider','Min',0,'Max',1,...
          'Value',obj.top_pk_sv,'SliderStep',[0.05 0.2]);
%---- bottom PeakRatio slider
obj.bottom_panel.bottomPeakRatioSlider = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.bottomPeakRatioSlider,'Style','slider','Min',0,'Max',1,...
          'Value',obj.bot_pk_sv,'SliderStep',[0.05 0.2]);   

%==============================
%---- Repulse name label
obj.bottom_panel.repulseLabel = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.repulseLabel,'Style','text');
set(obj.bottom_panel.repulseLabel,'String','Repulse:');
%---- Repulse slider
obj.bottom_panel.repulseSlider = uicontrol('Parent',obj.bottom_panel.handle);
set(obj.bottom_panel.repulseSlider,'Style','slider','Min',0,'Max',1,...
          'Value',obj.rep_sv,'SliderStep',[0.05 0.2]);   
        
%----------------------------------------------------------------------------------------------
% set up bottom panel table

obj.bottom_panel.table.ui=obj.bottom_panel.handle;
obj.bottom_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.bottom_panel.table.height_margin = NaN*zeros(30,30);
obj.bottom_panel.table.false_width = NaN*zeros(30,30);
obj.bottom_panel.table.false_height = NaN*zeros(30,30);
obj.bottom_panel.table.offset = [0 0];


row = 1; col = 1;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.topSmooth;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 1; col = 2;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.topSmoothSlider;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 2; col = 1;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.bottomSmooth;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 2; col=2;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.bottomSmoothSlider;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 3; col = 1;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.topPeakRatio;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 3; col = 2;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.topPeakRatioSlider;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 4; col = 1;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.bottomPeakRatio;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 4; col = 2;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.bottomPeakRatioSlider;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 5; col = 1;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.repulseLabel;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;

row = 5; col=2;
obj.bottom_panel.table.handles{row,col}   = obj.bottom_panel.repulseSlider;
obj.bottom_panel.table.width(row,col)     = inf;
obj.bottom_panel.table.height(row,col)    = inf;
obj.bottom_panel.table.width_margin(row,col)= 1.5;
obj.bottom_panel.table.height_margin(row,col)=1.5;



clear row col

%----------------------------------------------------------------------------------------------
% Draw table
table_draw(obj.bottom_panel.table);
set(obj.h_fig,'visible','on');

return