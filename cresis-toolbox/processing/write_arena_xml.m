function doc = write_arena_xml(doc,node,arena)
% doc = write_arena_xml(doc,node,arena)
%
% obj = xmlread(xml_fn);

% Author: John Paden

if ~strcmpi(arena.version,'1')
  error('Unsupported version %s', arena.version);
end

%% Initialization:
if strcmpi(node,'init')
  doc = com.mathworks.xml.XMLUtils.createDocument('system');
  system = doc.getDocumentElement;
  
  configs = doc.createElement('configs');
  system.appendChild(configs);
end

%% Central Timing Unit:
if strcmpi(node,'ctu_0013')
  system = doc.getFirstChild;
  configs = system.getFirstChild;
  
  config = doc.createElement('config');
  configs.appendChild(config);
  config.setAttribute('type','ctu_0013');
  
  child = doc.createElement('name'); config.appendChild(child);
  child.appendChild(doc.createTextNode('ctuConfig0'));
  
  child = doc.createElement('description'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  
  output = 0;
  for PA = 1:8
    child = doc.createElement('signalAlias'); config.appendChild(child);
    child.appendChild(doc.createTextNode(''));
    grandchild = doc.createElement('signal'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('output %d',output))); output = output + 1;
    grandchild = doc.createElement('name'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('PA ENA %d',PA)));
  end
  
  for signal_name = {'T/R','ISO','EPRI','PRI','EPRI','PRI','EPRI','PRI'}
    child = doc.createElement('signalAlias'); config.appendChild(child);
    child.appendChild(doc.createTextNode(''));
    grandchild = doc.createElement('signal'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('output %d',output))); output = output + 1;
    grandchild = doc.createElement('name'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(signal_name{1}));
  end
  
  child = doc.createElement('numSegments'); config.appendChild(child);
  child.appendChild(doc.createTextNode('4'));
  
  segment_states{1} = [
    0 1 1 0 % PA ENA 1
    0 1 1 0 % PA ENA 2
    0 1 1 0 % PA ENA 3
    0 1 1 0 % PA ENA 4
    0 1 1 0 % PA ENA 5
    0 1 1 0 % PA ENA 6
    0 1 1 0 % PA ENA 7
    0 1 1 0 % PA ENA 8
    0 1 1 0 % T/R
    0 1 1 0 % ISO
    1 0 0 0 % EPRI
    0 1 0 0 % PRI
    1 0 0 0 % EPRI
    0 1 0 0 % PRI
    1 0 0 0 % EPRI
    0 1 0 0 % PRI
    ];
  segment_states{2} = [
    0 1 1 0 % PA ENA 1
    0 1 1 0 % PA ENA 2
    0 1 1 0 % PA ENA 3
    0 1 1 0 % PA ENA 4
    0 1 1 0 % PA ENA 5
    0 1 1 0 % PA ENA 6
    0 1 1 0 % PA ENA 7
    0 1 1 0 % PA ENA 8
    0 1 1 0 % T/R
    0 1 1 0 % ISO
    0 0 0 0 % EPRI
    0 1 0 0 % PRI
    0 0 0 0 % EPRI
    0 1 0 0 % PRI
    0 0 0 0 % EPRI
    0 1 0 0 % PRI
    ];
  
  num_modes = 0;
  for wf = 1:length(arena.wfs)
    if arena.wfs(wf).Tpd > 0e-6
      segment_times = [0.1 0.2 arena.wfs(wf).Tpd*1e6+2.2 arena.PRI*1e6];
    else
      segment_times = [0.1 0.2 arena.wfs(wf).Tpd*1e6+1 arena.PRI*1e6];
    end
    if wf == 1
      wf_modes = 3;
      segment_states_idx = [1 2 2];
    else
      wf_modes = 2;
      segment_states_idx = [2 2 2];
    end
    
    for mode = 0:wf_modes-1
      child = doc.createElement('mode'); config.appendChild(child);
      child.appendChild(doc.createTextNode(''));
      grandchild = doc.createElement('id'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%d',num_modes+mode)));
      grandchild = doc.createElement('segmentTimes'); child.appendChild(grandchild);
      segment_time_str = [sprintf('%g', segment_times(1)), sprintf(' %g', segment_times(2:end))];
      grandchild.appendChild(doc.createTextNode(segment_time_str));
      grandchild = doc.createElement('segmentStates'); child.appendChild(grandchild);
      idx = segment_states_idx(mode+1);
      segment_state_str = lower(dec2hex(bin2dec(char(segment_states{idx}(end:-1:1,1).'+48)),8));
      for state_idx = 2:size(segment_states{idx},2)
        segment_state_str = cat(2,segment_state_str, ' ', ...
          lower(dec2hex(bin2dec(char(segment_states{idx}(end:-1:1,state_idx).'+48)),8)));
      end
      grandchild.appendChild(doc.createTextNode(segment_state_str));
    end
    num_modes = num_modes + wf_modes;
  end
  
  child = doc.createElement('pps'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  grandchild = doc.createElement('input'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  grandchild = doc.createElement('polarity'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  
  child = doc.createElement('nmea'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  grandchild = doc.createElement('input'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  grandchild = doc.createElement('polarity'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  
  child = doc.createElement('pscIntr'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  grandchild = doc.createElement('input'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  grandchild = doc.createElement('polarity'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  
  child = doc.createElement('encoder'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  grandchild = doc.createElement('input'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  grandchild = doc.createElement('polarity'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('2'));
  
  child = doc.createElement('marker'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  grandchild = doc.createElement('input'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  grandchild = doc.createElement('polarity'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('1'));
  
  child = doc.createElement('direction'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  grandchild = doc.createElement('input'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0'));
  grandchild = doc.createElement('polarity'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('1'));
  
end

%% DAC:
if strcmpi(node,'dac-ad9129_0014')
  system = doc.getFirstChild;
  configs = system.getFirstChild;
  
  for dac = 0:7
    config = doc.createElement('config'); configs.appendChild(config);
    config.setAttribute('type','dac-ad9129_0014');
    
    child = doc.createElement('name'); config.appendChild(child);
    child.appendChild(doc.createTextNode(sprintf('dacConfig%d',dac)));
    
    child = doc.createElement('description'); config.appendChild(child);
    child.appendChild(doc.createTextNode(''));
    
    
    num_modes = 0;
    for wf = 1:length(arena.wfs)
      Tpd = round(arena.wfs(wf).Tpd*1e6);
      if wf == 1
        wf_modes = 3;
      else
        wf_modes = 2;
      end
      
      for mode = 0:wf_modes-1
        child = doc.createElement('mode'); config.appendChild(child);
        child.appendChild(doc.createTextNode(''));
        grandchild = doc.createElement('id'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%d',num_modes+mode)));
        grandchild = doc.createElement('enabled'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('1'));
        grandchild = doc.createElement('delay'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('0.000000'));
        grandchild = doc.createElement('vDelayEnabled'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('0'));
        grandchild = doc.createElement('vDelay'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('0.000000'));
        grandchild = doc.createElement('config'); child.appendChild(grandchild);
        if arena.wfs(wf).enabled(dac+1)
          if mod(mode,2)
            grandchild.appendChild(doc.createTextNode(sprintf('waveformCh%d_%dus_180',dac,Tpd)));
          else
            grandchild.appendChild(doc.createTextNode(sprintf('waveformCh%d_%dus',dac,Tpd)));
          end
        else
          grandchild.appendChild(doc.createTextNode('No_Tx'));
        end
        grandchild.setAttribute('type','dac-ad9129_0014_waveform');
      end
      num_modes = num_modes + wf_modes;
    end
  end
end

%% Waveforms:
if strcmpi(node,'dac-ad9129_0014_waveform')
  system = doc.getFirstChild;
  configs = system.getFirstChild;

  wf = 1; dac = 1; zeropimod = 0;
  fs = arena.fs*1e6;
  fc = arena.wfs(wf).fc;
  BW = arena.wfs(wf).BW;
  Tpd = arena.wfs(wf).Tpd;
  alpha = arena.wfs(wf).tukey;
  equal.delay = arena.wfs(wf).delay;
  equal.phase = arena.wfs(wf).phase;
  equal.scale = arena.wfs(wf).scale;
  Nt = round((arena.wfs(wf).Tpd+equal.delay/1e9) * fs);
  
  config = doc.createElement('config'); configs.appendChild(config);
  config.setAttribute('type','dac-ad9129_0014_waveform');
  
  child = doc.createElement('name'); config.appendChild(child);
  child.appendChild(doc.createTextNode('No_Tx'));
  
  child = doc.createElement('description'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  
  child = doc.createElement('sampFreq'); config.appendChild(child);
  child.appendChild(doc.createTextNode(sprintf('%f',fs/1e6)));
  
  child = doc.createElement('pulse'); config.appendChild(child);
  grandchild = doc.createElement('name'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('Pulse'));
  grandchild = doc.createElement('centerFreq'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',fc(dac+1)/1e6)));
  grandchild = doc.createElement('bandwidth'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',BW(dac+1)/1e6)));
  grandchild = doc.createElement('initialDelay'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.delay(dac+1)*1e-3)));
  grandchild = doc.createElement('initialPhase'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.phase(dac+1)+zeropimod)));
  grandchild = doc.createElement('afterPulseDelay'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('1.000000'));
  grandchild = doc.createElement('taper'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('Tukey'));
  grandchild = doc.createElement('alpha'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',alpha)));
  grandchild = doc.createElement('scale'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0.000000'));
  grandchild = doc.createElement('numPoints'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%d',Nt(dac+1))));
  grandchild = doc.createElement('Filename'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(''));  

  waveform_names = {};
  for wf = 1:length(arena.wfs)
    fs = arena.fs*1e6;
    fc = arena.wfs(wf).fc;
    BW = arena.wfs(wf).BW;
    Tpd = arena.wfs(wf).Tpd;
    alpha = arena.wfs(wf).tukey;
    equal.delay = arena.wfs(wf).delay;
    equal.phase = arena.wfs(wf).phase;
    equal.scale = arena.wfs(wf).scale;
    Nt = round((arena.wfs(wf).Tpd+equal.delay/1e9) * fs);
    
    %[0.63 ]
    %  0.1652 0.326800 0.511500 0.63 0.6300 0.511500 0.3268 0.1652
    %    chebwin(8,30)
    
    for dac=0:7
      for zeropimod = [0 180]
        if zeropimod == 0
          new_waveform_name = sprintf('waveformCh%d_%.0fus',dac,Tpd*1e6);
        else
          new_waveform_name = sprintf('waveformCh%d_%.0fus_180',dac,Tpd*1e6);
        end
        if any(strcmpi(new_waveform_name,waveform_names))
          continue;
        end
        waveform_names{end+1} = new_waveform_name;
        
        config = doc.createElement('config'); configs.appendChild(config);
        config.setAttribute('type','dac-ad9129_0014_waveform');
        
        child = doc.createElement('name'); config.appendChild(child);
        child.appendChild(doc.createTextNode(waveform_names{end}));
        
        child = doc.createElement('description'); config.appendChild(child);
        child.appendChild(doc.createTextNode(''));
        
        child = doc.createElement('sampFreq'); config.appendChild(child);
        child.appendChild(doc.createTextNode(sprintf('%f',fs/1e6)));
        
        child = doc.createElement('pulse'); config.appendChild(child);
        grandchild = doc.createElement('name'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('Pulse'));
        grandchild = doc.createElement('centerFreq'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',fc(dac+1)/1e6)));
        grandchild = doc.createElement('bandwidth'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',BW(dac+1)/1e6)));
        grandchild = doc.createElement('initialDelay'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.delay(dac+1)*1e-3)));
        grandchild = doc.createElement('initialPhase'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.phase(dac+1)+zeropimod)));
        grandchild = doc.createElement('afterPulseDelay'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('1.000000'));
        grandchild = doc.createElement('taper'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('Tukey'));
        grandchild = doc.createElement('alpha'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',alpha)));
        grandchild = doc.createElement('scale'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.scale(dac+1))));
        grandchild = doc.createElement('numPoints'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%d',Nt(dac+1))));
        grandchild = doc.createElement('Filename'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(''));
      end
    end
  end
  
end

%% Primary sequence controller:
if strcmpi(node,'psc_0001')
  system = doc.getFirstChild;
  configs = system.getFirstChild;
  
  config = doc.createElement('config'); configs.appendChild(config);
  config.setAttribute('type','psc_0001');
  
  child = doc.createElement('name'); config.appendChild(child);
  child.appendChild(doc.createTextNode('pscConfig0'));
  
  child = doc.createElement('description'); config.appendChild(child);
  child.appendChild(doc.createTextNode(''));
  
  num_psc = 0;
  for wf = 1:length(arena.wfs)
    if wf == 1
      if arena.wfs(wf).presums == 1
        psc_repeat = [0 0 0];
        psc_name = {sprintf('%.0fus, EPRI, Zero',arena.wfs(wf).Tpd*1e6)};
      elseif arena.wfs(wf).presums == 2
        psc_repeat = [0 0 0;
          1 0 0];
        psc_name = {sprintf('%.0fus, EPRI, Zero',arena.wfs(wf).Tpd*1e6)
          sprintf('%.0fus, PRI, PI',arena.wfs(wf).Tpd*1e6)};
      elseif ~mod(arena.wfs(wf).presums,2)
        psc_repeat = [0 0 0;
          1 0 0;
          2 1 (arena.wfs(wf).presums-4)/2;
          1 0 0];
        psc_name = {sprintf('%.0fus, EPRI, Zero',arena.wfs(wf).Tpd*1e6)
          sprintf('%.0fus, PRI, PI',arena.wfs(wf).Tpd*1e6)
          sprintf('%.0fus, PRI, Zero',arena.wfs(wf).Tpd*1e6)
          sprintf('%.0fus, PRI, PI',arena.wfs(wf).Tpd*1e6)};
      else
        error('Odd number of presums, %d, greater than 1 not supported.', ...
          arena.wfs(wf).presums);
      end
    else
      if arena.wfs(wf).presums == 1
        psc_repeat = [2*wf-1 0 0];
        psc_name = {sprintf('%.0fus, PRI, Zero',arena.wfs(wf).Tpd*1e6)};
      elseif ~mod(arena.wfs(wf).presums,2)
        psc_repeat = [2*wf-1 0 0;
          2*wf num_psc (arena.wfs(wf).presums-2)/2];
        psc_name = {sprintf('%.0fus, PRI, Zero',arena.wfs(wf).Tpd*1e6)
          sprintf('%.0fus, PRI, PI',arena.wfs(wf).Tpd*1e6)};
      else
        error('Odd number of presums, %d, greater than 1 not supported.', ...
          arena.wfs(wf).presums);
      end
    end
    num_psc = num_psc + size(psc_repeat,1);
    for psc_idx = 1:size(psc_repeat,1)
      mode = psc_repeat(psc_idx,1);
      repeat_mode = psc_repeat(psc_idx,2);
      repeat_count = psc_repeat(psc_idx,3);
      child = doc.createElement('sequence'); config.appendChild(child);
      child.setAttribute('type','primary');
      grandchild = doc.createElement('mode'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%d',mode)));
      grandchild = doc.createElement('name'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(psc_name{psc_idx}));
      grandchild = doc.createElement('period'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%g',arena.PRI*1e6)));
      grandchild = doc.createElement('repeatTo'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%d',repeat_mode)));
      grandchild = doc.createElement('repeatCount'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%d',repeat_count)));
      grandchild = doc.createElement('interruptible'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode('0'));
    end
  end
  child = doc.createElement('interruptEna'); config.appendChild(child);
  child.appendChild(doc.createTextNode('0'));
end

%% Subsystems
if strcmpi(node,'subsystems')
  system = doc.getFirstChild;
  configs = system.getFirstChild;
  
  % AWG
  desiredAlignMin = [0 10 -5 10 0 10 0 10];
  desiredAlignMax = desiredAlignMin+20;
  for awg = 0:3
    config = doc.createElement('subSystem'); system.appendChild(config);
    
    child = doc.createElement('name'); config.appendChild(child);
    child.appendChild(doc.createTextNode(sprintf('arena-awg%d',awg)));
    
    child = doc.createElement('interface'); config.appendChild(child);
    child.appendChild(doc.createTextNode('eth0'));
    child.setAttribute('type','nic');
    
    child = doc.createElement('port'); config.appendChild(child);
    child.appendChild(doc.createTextNode('10000'));
    
    child = doc.createElement('config'); config.appendChild(child);
    child.appendChild(doc.createTextNode(''));
    child.setAttribute('type','psc_0001');
    
    child = doc.createElement('subSystem'); config.appendChild(child);
    grandchild = doc.createElement('name'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('awg%d',awg*2)));
    grandchild = doc.createElement('disableSync'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode('0'));
    grandchild = doc.createElement('dacClk'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode('1600.000000'));
    grandchild = doc.createElement('mixMode'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode('0'));
    grandchild = doc.createElement('desiredAlignMin'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('%d',desiredAlignMin(awg*2+1))));
    grandchild = doc.createElement('desiredAlignMax'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('%d',desiredAlignMax(awg*2+1))));
    grandchild = doc.createElement('config'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('dacConfig%d',awg*2)));
    grandchild.setAttribute('type','dac-ad9129_0014');
    
    child = doc.createElement('subSystem'); config.appendChild(child);
    grandchild = doc.createElement('name'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('awg%d',awg*2+1)));
    grandchild = doc.createElement('disableSync'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode('0'));
    grandchild = doc.createElement('dacClk'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode('1600.000000'));
    grandchild = doc.createElement('mixMode'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode('0'));
    grandchild = doc.createElement('desiredAlignMin'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('%d',desiredAlignMin(awg*2+2))));
    grandchild = doc.createElement('desiredAlignMax'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('%d',desiredAlignMax(awg*2+2))));
    grandchild = doc.createElement('config'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('dacConfig%d',awg*2+1)));
    grandchild.setAttribute('type','dac-ad9129_0014');
  end
  
  % CTU
  config = doc.createElement('subSystem'); system.appendChild(config);
  
  child = doc.createElement('name'); config.appendChild(child);
  child.appendChild(doc.createTextNode('arena-ctu0'));
  
  child = doc.createElement('interface'); config.appendChild(child);
  child.appendChild(doc.createTextNode('eth0'));
  child.setAttribute('type','nic');
  
  child = doc.createElement('port'); config.appendChild(child);
  child.appendChild(doc.createTextNode('10000'));
  
  child = doc.createElement('config'); config.appendChild(child);
  child.appendChild(doc.createTextNode('pscConfig0'));
  child.setAttribute('type','psc_0001');
  
  child = doc.createElement('subSystem'); config.appendChild(child);
  grandchild = doc.createElement('name'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('ctu0'));
  grandchild = doc.createElement('config'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('ctuConfig0'));
  grandchild.setAttribute('type','ctu_0013');
end

return;
