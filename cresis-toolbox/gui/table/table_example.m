function container_gui

global hui;

hui.panel.figure = figure(1); clf;
set(hui.panel.figure,'Pointer','arrow');
set(hui.panel.figure,'Visible','off');
set(hui.panel.figure,'Tag','echo');
set(hui.panel.figure,'NumberTitle','on');
set(hui.panel.figure,'Name','CrATE : Control Panel');
set(hui.panel.figure,'Color',[0.4 0.4 0.4]);
set(hui.panel.figure,'Visible','on');

hui.panel.mapPanel = uipanel('Parent',hui.panel.figure);
set(hui.panel.mapPanel,'Position',[0.01 0.01 0.48 0.98]);
set(hui.panel.mapPanel,'Title','Map Controls');
set(hui.panel.mapPanel,'TitlePosition','CenterTop');
set(hui.panel.mapPanel,'HighlightColor',[0.8 0.8 0.8]);
set(hui.panel.mapPanel,'ShadowColor',[0.6 0.6 0.6]);

hui.mapPanel.baseMapText = uicontrol('Parent',hui.panel.mapPanel);
set(hui.mapPanel.baseMapText,'Style','Text');
set(hui.mapPanel.baseMapText,'String','Current Basemap');
set(hui.mapPanel.baseMapText,'HorizontalAlignment','Left');

hui.mapPanel.baseMapMenu = uicontrol('Parent',hui.panel.mapPanel);
set(hui.mapPanel.baseMapMenu,'Style','PopupMenu');
set(hui.mapPanel.baseMapMenu,'Callback',[]);
set(hui.mapPanel.baseMapMenu,'BackgroundColor',[1 1 1]);
set(hui.mapPanel.baseMapMenu,'ForegroundColor',[0 0 0]);
menuString{1} = 'ZERO';
menuString{2} = 'ONE';
set(hui.mapPanel.baseMapMenu,'String',menuString);
clear menuString

hui.table.ui = hui.panel.mapPanel;
hui.table.handles(1,1) = hui.mapPanel.baseMapText;
hui.table.width(1) = 220;
hui.table.height(1) = 20;
hui.table.handles(2,1) = hui.mapPanel.baseMapMenu;
hui.table.height(2) = 20;
table_draw(hui.table);

return;

