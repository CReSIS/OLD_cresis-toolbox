function doc = write_arena_xml(doc,node,arena)
% doc = write_arena_xml(doc,node,arena)
%
% obj = xmlread(xml_fn);
%
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
  for signal_name = arena.TTL_names
    child = doc.createElement('signalAlias'); config.appendChild(child);
    child.appendChild(doc.createTextNode(''));
    grandchild = doc.createElement('signal'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(sprintf('output %d',output))); output = output + 1;
    grandchild = doc.createElement('name'); child.appendChild(grandchild);
    grandchild.appendChild(doc.createTextNode(signal_name{1}));
  end
  
  child = doc.createElement('numSegments'); config.appendChild(child);
  child.appendChild(doc.createTextNode(sprintf('%d',size(arena.TTL_states{1},2))));
  
  num_modes = 0;
  for wf = 1:length(arena.wfs)
    segment_times = [arena.TTL_time(1:2) arena.wfs(wf).Tpd*1e6+arena.TTL_time(3) arena.PRI*1e6];
    if wf == 1
      wf_modes = length(arena.zeropimods)+1;
      segment_states_idx = [1 2*ones(size(arena.zeropimods))];
    else
      wf_modes = length(arena.zeropimods);
      segment_states_idx = [2*ones(size(arena.zeropimods))];
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
      segment_state_str = lower(dec2hex(bin2dec(char(arena.TTL_states{idx}(end:-1:1,1).'+48)),8));
      for state_idx = 2:size(arena.TTL_states{idx},2)
        segment_state_str = cat(2,segment_state_str, ' ', ...
          lower(dec2hex(bin2dec(char(arena.TTL_states{idx}(end:-1:1,state_idx).'+48)),8)));
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

%% DAC ACCUM:
if strcmpi(node,'dac-ad9129_0014')
  system = doc.getFirstChild;
  configs = system.getFirstChild;
  
  for dac = arena.dacs
    config = doc.createElement('config'); configs.appendChild(config);
    config.setAttribute('type','dac-ad9129_0014');
    
    child = doc.createElement('name'); config.appendChild(child);
    child.appendChild(doc.createTextNode(sprintf('dacConfig%d',dac)));
    
    child = doc.createElement('description'); config.appendChild(child);
    child.appendChild(doc.createTextNode(''));
    
    num_modes = 0;
    for wf = 1:length(arena.wfs)
      zeropimods = arena.wfs(wf).zeropimods(:).';
      Tpd = round(arena.wfs(wf).Tpd*1e6);
      if wf == 1
        wf_modes = 1+length(zeropimods);
      else
        wf_modes = length(zeropimods);
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
          grandchild.appendChild(doc.createTextNode(sprintf('waveformCh%d_%dus_%.0f',dac,Tpd,zeropimods(1+mod(mode,length(zeropimods))) )));
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

  wf = 1; dac_idx = 1; dac = arena.dacs(dac_idx); zeropimod = 0;
  fs = arena.dacs_sampFreq(dac_idx);
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
  grandchild.appendChild(doc.createTextNode(sprintf('%f',fc(dac_idx)/1e6)));
  grandchild = doc.createElement('bandwidth'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',BW(dac_idx)/1e6)));
  grandchild = doc.createElement('initialDelay'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.delay(dac_idx)*1e-3)));
  grandchild = doc.createElement('initialPhase'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.phase(dac_idx)+zeropimod)));
  grandchild = doc.createElement('afterPulseDelay'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('1.000000'));
  grandchild = doc.createElement('taper'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('Tukey'));
  grandchild = doc.createElement('alpha'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%f',alpha)));
  grandchild = doc.createElement('scale'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode('0.000000'));
  grandchild = doc.createElement('numPoints'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(sprintf('%d',Nt(dac_idx))));
  grandchild = doc.createElement('Filename'); child.appendChild(grandchild);
  grandchild.appendChild(doc.createTextNode(''));  

  waveform_names = {};
  for wf = 1:length(arena.wfs)
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
    
    for dac_idx=1:length(arena.dacs)
      dac = arena.dacs(dac_idx);
      fs = arena.dacs_sampFreq(dac_idx);
      for zeropimod = arena.wfs(wf).zeropimods(:).'
        new_waveform_name = sprintf('waveformCh%d_%.0fus_%.0f',dac,Tpd*1e6,zeropimod);
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
        grandchild.appendChild(doc.createTextNode(sprintf('%f',fc(dac_idx)/1e6)));
        grandchild = doc.createElement('bandwidth'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',BW(dac_idx)/1e6)));
        grandchild = doc.createElement('initialDelay'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.delay(dac_idx)*1e-3)));
        grandchild = doc.createElement('initialPhase'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.phase(dac_idx)+zeropimod)));
        grandchild = doc.createElement('afterPulseDelay'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('1.000000'));
        grandchild = doc.createElement('taper'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode('Tukey'));
        grandchild = doc.createElement('alpha'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',alpha)));
        grandchild = doc.createElement('scale'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%f',equal.scale(dac_idx))));
        grandchild = doc.createElement('numPoints'); child.appendChild(grandchild);
        grandchild.appendChild(doc.createTextNode(sprintf('%d',Nt(dac_idx))));
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
  num_wf = 0;
  for wf = 1:length(arena.wfs)
    zeropimods = arena.wfs(wf).zeropimods(:).';
    if wf == 1
      if arena.wfs(wf).presums == 1
        psc_repeat = [0 0 0];
        zeropimod = zeropimods(1);
        psc_name = {sprintf('%.0fus, EPRI, %d',arena.wfs(wf).Tpd*1e6, zeropimod)};
        num_pulse_states = num_pulse_states + 1;
      elseif arena.wfs(wf).presums == length(zeropimods)
        psc_repeat = zeros(length(zeropimods),3);
        psc_repeat(:,1) = 0:length(zeropimods)-1;
        psc_name{1} = sprintf('%.0fus, EPRI, %d',arena.wfs(wf).Tpd*1e6, zeropimods(1));
        for zeropimod = zeropimods(2:end)
          psc_name{end+1} = sprintf('%.0fus, PRI, %d',arena.wfs(wf).Tpd*1e6, zeropimod);
        end
      elseif mod(arena.wfs(wf).presums,length(zeropimods)) == 0
        psc_repeat = zeros(2*length(zeropimods),3);
        psc_repeat(:,1) = [0:length(zeropimods), 1:length(zeropimods)-1];
        psc_repeat(end,2:3) = [length(zeropimods), (arena.wfs(wf).presums-2*length(zeropimods))/length(zeropimods)];
        psc_name{1} = sprintf('%.0fus, EPRI, %d',arena.wfs(wf).Tpd*1e6, zeropimods(1));
        for zeropimod = zeropimods(2:end)
          psc_name{end+1} = sprintf('%.0fus, PRI, %d',arena.wfs(wf).Tpd*1e6, zeropimod);
        end
        for zeropimod = zeropimods
          psc_name{end+1} = sprintf('%.0fus, PRI, %d',arena.wfs(wf).Tpd*1e6, zeropimod);
        end
      else
        error('Presums %d>1 so presums must be a multiple of zeropimods length=%d.', ...
          arena.wfs(wf).presums, length(zeropimods));
      end
    else
      if arena.wfs(wf).presums == 1
        psc_repeat = [length(zeropimods)*(wf-1)+1 0 0];
        psc_name = {sprintf('%.0fus, PRI, Zero',arena.wfs(wf).Tpd*1e6)};
      elseif ~mod(arena.wfs(wf).presums,2)
        psc_repeat = zeros(length(zeropimods),3);
        psc_repeat(:,1) = length(zeropimods)*(wf-1) + (1:length(zeropimods));
        psc_repeat(end,2:3) = [num_psc, (arena.wfs(wf).presums-length(zeropimods))/length(zeropimods)];
        for zeropimod = zeropimods
          psc_name{end+1} = sprintf('%.0fus, PRI, %d',arena.wfs(wf).Tpd*1e6, zeropimod);
        end
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
  for awg_idx = 1:length(arena.awg)
    awg = arena.awg(awg_idx).awg;
    
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

    for dac_idx = 1:length(arena.awg(awg_idx).dacs)
      dac = arena.awg(awg_idx).dacs(dac_idx);
      child = doc.createElement('subSystem'); config.appendChild(child);
      grandchild = doc.createElement('name'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('awg%d',dac)));
      grandchild = doc.createElement('disableSync'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode('0'));
      grandchild = doc.createElement('dacClk'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%g',arena.awg(awg_idx).dacClk(dac_idx)/1e6)));
      grandchild = doc.createElement('mixMode'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode('0'));
      grandchild = doc.createElement('desiredAlignMin'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%d',arena.awg(awg_idx).desiredAlignMin(dac_idx))));
      grandchild = doc.createElement('desiredAlignMax'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('%d',arena.awg(awg_idx).desiredAlignMax(dac_idx))));
      grandchild = doc.createElement('config'); child.appendChild(grandchild);
      grandchild.appendChild(doc.createTextNode(sprintf('dacConfig%d',dac)));
      grandchild.setAttribute('type','dac-ad9129_0014');
    end
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
