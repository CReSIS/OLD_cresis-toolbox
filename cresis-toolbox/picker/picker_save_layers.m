function picker_save_layers(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if ishandle(hui.pickfig.handle)
  fn = gCtrl.source.fns{gCtrl.source.cur_pick};
  [path name] = fileparts(fn);
  fprintf('Saving layer %s\n', name);
  gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{1}.value{1}.data ...
    = get(hui.pickfig.layer_h(1),'YData') / 1e6;
  gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{1}.value{2}.data ...
    = get(hui.pickfig.layer_h(2),'YData') / 1e6;
  gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{2}.value{1}.data ...
    = get(hui.pickfig.layer_h(3),'YData') / 1e6;
  gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{2}.value{2}.data ...
    = get(hui.pickfig.layer_h(4),'YData') / 1e6;
  for cur_layer = 1:2
    good = get(hui.pickfig.quality_h(3*(cur_layer-1)+1),'YData');
    moderate = get(hui.pickfig.quality_h(3*(cur_layer-1)+2),'YData');
    derived = get(hui.pickfig.quality_h(3*(cur_layer-1)+3),'YData');
    qual = 1*isfinite(good) + 2*isfinite(moderate) + 3*isfinite(derived);
    qual(qual < 1 | qual > 3) = 1;
    gCtrl.source.layers{gCtrl.source.cur_pick}.layerData{cur_layer}.quality = qual;
  end

  layerData = gCtrl.source.layers{gCtrl.source.cur_pick}.layerData;

  if 0
    % Debug code for comparing to currently stored data
    tmp = load(fn);
    figure(11); clf;
    subplot(4,1,1);
    plot(tmp.layerData{1}.quality);
    subplot(4,1,2);
    plot(layerData{1}.quality,'kx');
    subplot(4,1,3);
    plot(tmp.layerData{2}.quality,'r');
    subplot(4,1,4);
    plot(layerData{2}.quality,'gx');
    figure(12); clf;
    plot(tmp.layerData{1}.value{1}.data,'bo');
    hold on;
    plot(tmp.layerData{1}.value{2}.data,'b--');
    plot(layerData{1}.value{1}.data+3e-6,'kx');
    plot(layerData{1}.value{2}.data+3e-6,'k--');
    figure(13); clf;
    plot(tmp.layerData{2}.value{1}.data,'bo');
    hold on;
    plot(tmp.layerData{2}.value{2}.data,'b--');
    plot(layerData{2}.value{1}.data+3e-6,'kx');
    plot(layerData{2}.value{2}.data+3e-6,'k--');
  else
    save(fn,'-APPEND','layerData');
    % Update frames list and picker title since this frame has no
    % longer been modified since last saved.
    gCtrl.source.modified(gCtrl.source.cur_pick) = ' ';
    menuString = get(hui.fig.ctrl_panel.framesLB,'String');
    menuString{gCtrl.source.cur_pick}(end) = gCtrl.source.modified(gCtrl.source.cur_pick);
    set(hui.fig.ctrl_panel.framesLB,'String',menuString);
    title_str = get(get(hui.pickfig.axes.handle,'Title'),'String');
    title_str(end) = gCtrl.source.modified(gCtrl.source.cur_pick);
    set(get(hui.pickfig.axes.handle,'Title'),'String',title_str);
  end
  
  fprintf('  Done\n');
end

end