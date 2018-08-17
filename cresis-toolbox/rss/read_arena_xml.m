function settings = read_arena_xml(xml_fn)
% settings = read_arena_xml(xml_fn)
%
% settings = read_arena_xml('c:\temp\config_GHSR.xml');
%
% Author: John Paden

[xml_fn_dir,xml_fn_name] = fileparts(xml_fn);
xml_type_idx = find(xml_fn_name=='_',1,'last');
config_xml_fn = fullfile(xml_fn_dir,[xml_fn_name(1:xml_type_idx), 'config.xml']);

settings.xml_fn = xml_fn;
settings.xml_fname = fname_info_arena(xml_fn);
settings.config_xml_fn = config_xml_fn;
settings.max_num_bins = 0;

doc = xmlread(xml_fn);
%   sys = doc.getElementsByTagName('system'); sys = sys.item(0);

doc_cfg = xmlread(config_xml_fn);
%   sys = doc.getElementsByTagName('system'); sys = sys.item(0);

%   for idx = 1:sys.getLength
%     it = sys.item(idx-1);
%     if it.getNodeType==1
%       name = it.getNodeName;
%       type = it.getAttribute('type');
%       fprintf('%s type %s\n', name.toCharArray, type.toCharArray);
%     end
%   end

% Import the XPath classes
import javax.xml.xpath.*

% Create an XPath expression.
factory = XPathFactory.newInstance;
xpath = factory.newXPath;

% Get the system type
expression = xpath.compile('//system/@type');
nodeList = expression.evaluate(doc,XPathConstants.NODESET);
radar_name = nodeList.item(0).getTextContent.toCharArray;
radar_name = radar_name(:).';
settings.radar_name = radar_name;

% Get all the ADCs
expression = xpath.compile('//subSystem[starts-with(@type,"adc")]');
% expression = xpath.compile('//subSystem[@type="daq"]');
adcList = expression.evaluate(doc,XPathConstants.NODESET);

for adc_idx = 1:adcList.getLength
  
  % 1. Get the name of the ADC
  expression = xpath.compile('//subSystem[starts-with(@type,"adc")]/name');
  % expression = xpath.compile('//subSystem[@type="daq"]');
  nodeList = expression.evaluate(adcList.item(adc_idx-1),XPathConstants.NODESET);
  name = nodeList.item(0).getTextContent.toCharArray;
  name = name(:).';
  
  % 2. Search for the subSystem in the config XML
  expression = xpath.compile(sprintf('//subSystem/subSystem[name="%s"]',name));
  nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
  match = nodeList.item(0);
  
  % 3. Get the config name and type for this ADC
  expression = xpath.compile('config/@type');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_type = nodeList.item(0).getTextContent.toCharArray;
  config_type = config_type(:).';
  expression = xpath.compile('config');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_name = nodeList.item(0).getTextContent.toCharArray;
  config_name = config_name(:).';
  
  % Get the config associated with this ADC
  expression = xpath.compile(sprintf('//configs/config[(@type="%s" and name="%s")]',config_type,config_name));
  nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
  adc_cfg = nodeList.item(0);
  
  % Find the longest possible record size
  if strcmpi(config_type,'adc-ads42lb69_0010')
    % TOHFSounder
    
    % Get each subchannel
    expression = xpath.compile('subChannels/subChannel');
    subchannel_nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
    for subchannel_idx = 1:subchannel_nodeList.getLength
      subchannel_cfg = subchannel_nodeList.item(subchannel_idx-1);
      if isempty(subchannel_cfg)
        continue;
      end
      
      expression = xpath.compile('id');
      nodeList = expression.evaluate(subchannel_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      subchannel = nodeList.item(0);
      subchannel = subchannel.getTextContent.toCharArray;
      subchannel = str2double(subchannel(:).');
      
      expression = xpath.compile('mode');
      mode_nodeList = expression.evaluate(subchannel_cfg,XPathConstants.NODESET);
      for mode_idx = 1:mode_nodeList.getLength
        mode_cfg = mode_nodeList.item(mode_idx-1);
        if isempty(mode_cfg)
          continue;
        end
        
        expression = xpath.compile('id');
        nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
        if nodeList.getLength < 1 || isempty(nodeList.item(0))
          continue;
        end
        mode_latch = nodeList.item(0);
        mode_latch = mode_latch.getTextContent.toCharArray;
        mode_latch = str2double(mode_latch(:).');
        
        expression = xpath.compile('ncoFreq');
        nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
        if nodeList.getLength < 1 || isempty(nodeList.item(0))
          continue;
        end
        ncoFreq = nodeList.item(0);
        ncoFreq = ncoFreq.getTextContent.toCharArray;
        ncoFreq = str2double(ncoFreq(:).');
        settings.adc{mode_latch+1,subchannel+1}.ncoFreq = ncoFreq;
        
        expression = xpath.compile('cicDecimation');
        nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
        if nodeList.getLength < 1 || isempty(nodeList.item(0))
          continue;
        end
        cicDecimation = nodeList.item(0);
        cicDecimation = cicDecimation.getTextContent.toCharArray;
        cicDecimation = str2double(cicDecimation(:).');
        settings.adc{mode_latch+1,subchannel+1}.cicDecimation = cicDecimation;
      end
    end
    
    % Get each subchannel
    expression = xpath.compile('processing/subChannel');
    subchannel_nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
    for subchannel_idx = 1:subchannel_nodeList.getLength
      subchannel_cfg = subchannel_nodeList.item(subchannel_idx-1);
      if isempty(subchannel_cfg)
        continue;
      end
      
      expression = xpath.compile('id');
      nodeList = expression.evaluate(subchannel_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      subchannel = nodeList.item(0);
      subchannel = subchannel.getTextContent.toCharArray;
      subchannel = str2double(subchannel(:).');
      
      expression = xpath.compile('mode');
      mode_nodeList = expression.evaluate(subchannel_cfg,XPathConstants.NODESET);
      for mode_idx = 1:mode_nodeList.getLength
        mode_cfg = mode_nodeList.item(mode_idx-1);
        if isempty(mode_cfg)
          continue;
        end
        
        expression = xpath.compile('id');
        nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
        if nodeList.getLength < 1 || isempty(nodeList.item(0))
          continue;
        end
        mode_latch = nodeList.item(0);
        mode_latch = mode_latch.getTextContent.toCharArray;
        mode_latch = str2double(mode_latch(:).');
        
        expression = xpath.compile('digRx_RG');
        nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
        if nodeList.getLength < 1 || isempty(nodeList.item(0))
          continue;
        end
        digRx_RG = nodeList.item(0);
        digRx_RG = digRx_RG.getTextContent.toCharArray;
        digRx_RG = digRx_RG(:).';
        % Assumes simple range gate format "start:stop"
        [start,stop] = strtok(digRx_RG,':'); stop=stop(2:end);
        num_bins = str2double(stop) - str2double(start) + 1;
        if num_bins > settings.max_num_bins
          settings.max_num_bins = num_bins;
        end
        
        settings.adc{mode_latch+1,subchannel+1}.digRx_RG = digRx_RG;
      end
    end
    
  elseif strcmpi(config_type,'adc-isla214p50_0005')
    % KUSnow
    expression = xpath.compile('//subChannels/subChannel/mode/rg');
    nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
    for mode_idx = 1:nodeList.getLength
      modes = nodeList.item(mode_idx-1);
      range_gates = modes.getTextContent.toCharArray;
      range_gates = range_gates(:).';
      % Assumes simple range gate format "start:stop"
      [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
      num_bins = str2double(stop) - str2double(start) + 1;
      if num_bins > settings.max_num_bins
        settings.max_num_bins = num_bins;
      end
    end
    
  elseif strcmpi(config_type,'adc-ads42lb69_0010')
    % DopplerScat
    expression = xpath.compile('//processing/subChannel/mode/digRx_RG');
    nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
    for mode_idx = 1:nodeList.getLength
      modes = nodeList.item(mode_idx-1);
      range_gates = modes.getTextContent.toCharArray;
      range_gates = range_gates(:).';
      % Assumes simple range gate format "start:stop"
      [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
      num_bins = str2double(stop) - str2double(start) + 1;
      if num_bins > settings.max_num_bins
        settings.max_num_bins = num_bins;
      end
    end
    
  elseif strcmpi(config_type,'adc-ad9680_0017')
    % BAS Accumulation Radar, Dome Fuji RDS
    expression = xpath.compile('//subChannels/subChannel/integrator/rg');
    nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
    for mode_idx = 1:nodeList.getLength
      modes = nodeList.item(mode_idx-1);
      range_gates = modes.getTextContent.toCharArray;
      range_gates = range_gates(:).';
      % Assumes simple range gate format "start:stop"
      [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
      num_bins = str2double(stop) - str2double(start) + 1;
      if num_bins > settings.max_num_bins
        settings.max_num_bins = num_bins;
      end
    end
    
  else
    error('ADC type %s not supported.', config_type);
  end
end


% Get all the CTUs
expression = xpath.compile('//subSystem[starts-with(@type,"ctu")]');
ctuList = expression.evaluate(doc,XPathConstants.NODESET);

for ctu_idx = 1:ctuList.getLength
  
  % 1. Get the name of the CTU
  expression = xpath.compile('name');
  nodeList = expression.evaluate(ctuList.item(ctu_idx-1),XPathConstants.NODESET);
  name = nodeList.item(0).getTextContent.toCharArray;
  name = name(:).';
  
  % 2. Search for the subSystem in the config XML
  expression = xpath.compile(sprintf('//subSystem/subSystem[name="%s"]',name));
  nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
  match = nodeList.item(0);
  
  % 3. Get the config name and type for this ADC
  expression = xpath.compile('config/@type');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_type = nodeList.item(0).getTextContent.toCharArray;
  config_type = config_type(:).';
  expression = xpath.compile('config');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_name = nodeList.item(0).getTextContent.toCharArray;
  config_name = config_name(:).';
  
  % Get the config associated with this ADC
  expression = xpath.compile(sprintf('//configs/config[(@type="%s" and name="%s")]',config_type,config_name));
  nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
  ctu_cfg = nodeList.item(0);
  
  % Find the longest possible record size
  if strcmpi(config_type,'ctu_0013')
    % TOHFSounder
    
    % Get each subchannel
    expression = xpath.compile('mode');
    mode_nodeList = expression.evaluate(ctu_cfg,XPathConstants.NODESET);
    for mode_idx = 1:mode_nodeList.getLength
      mode_cfg = mode_nodeList.item(mode_idx-1);
      if isempty(mode_cfg)
        continue;
      end
      
      expression = xpath.compile('id');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      mode_latch = nodeList.item(0);
      mode_latch = mode_latch.getTextContent.toCharArray;
      mode_latch = str2double(mode_latch(:).');
      
      expression = xpath.compile('segmentTimes');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      segmentTimes = nodeList.item(0);
      segmentTimes = segmentTimes.getTextContent.toCharArray;
      segmentTimes = segmentTimes(:).';
      settings.ctu{mode_latch+1}.segmentTimes = segmentTimes;
      
      expression = xpath.compile('segmentStates');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      segmentStates = nodeList.item(0);
      segmentStates = segmentStates.getTextContent.toCharArray;
      segmentStates = segmentStates(:).';
      settings.ctu{mode_latch+1}.segmentStates = segmentStates;
    end
  end
  
end


% Get all the DACs
expression = xpath.compile('//subSystem[starts-with(@type,"dac")]');
dacList = expression.evaluate(doc,XPathConstants.NODESET);

for dac_idx = 1:dacList.getLength
  
  % 1. Get the name of the DAC
  expression = xpath.compile('name');
  nodeList = expression.evaluate(dacList.item(dac_idx-1),XPathConstants.NODESET);
  name = nodeList.item(0).getTextContent.toCharArray;
  name = name(:).';
  
  % 2. Search for the subSystem in the config XML
  expression = xpath.compile(sprintf('//subSystem/subSystem[name="%s"]',name));
  nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
  match = nodeList.item(0);
  if isempty(match)
    warning('DAC %s has no config assigned. Skipping.', name);
    continue;
  end
  
  % 3. Get the config name and type for this ADC
  expression = xpath.compile('config/@type');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_type = nodeList.item(0).getTextContent.toCharArray;
  config_type = config_type(:).';
  expression = xpath.compile('config');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_name = nodeList.item(0).getTextContent.toCharArray;
  config_name = config_name(:).';
  
  % Get the config associated with this ADC
  expression = xpath.compile(sprintf('//configs/config[(@type="%s" and name="%s")]',config_type,config_name));
  nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
  dac_cfg = nodeList.item(0);
  
  % Find the longest possible record size
  if strcmpi(config_type,'dac-ad9129_0012')
    % TOHFSounder
    
    % Get each subchannel
    expression = xpath.compile('mode');
    mode_nodeList = expression.evaluate(dac_cfg,XPathConstants.NODESET);
    for mode_idx = 1:mode_nodeList.getLength
      mode_cfg = mode_nodeList.item(mode_idx-1);
      if isempty(mode_cfg)
        continue;
      end
      
      expression = xpath.compile('id');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      mode_latch = nodeList.item(0);
      mode_latch = mode_latch.getTextContent.toCharArray;
      mode_latch = str2double(mode_latch(:).');
      
      expression = xpath.compile('delay');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      if nodeList.getLength < 1 || isempty(nodeList.item(0))
        continue;
      end
      delay = nodeList.item(0);
      delay = delay.getTextContent.toCharArray;
      delay = str2double(delay(:).');
      settings.dac{mode_latch+1}.delay = delay;
      
      % 3. Get the config name and type for this ADC
      expression = xpath.compile('config/@type');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      config_type = nodeList.item(0).getTextContent.toCharArray;
      config_type = config_type(:).';
      expression = xpath.compile('config');
      nodeList = expression.evaluate(mode_cfg,XPathConstants.NODESET);
      config_name = nodeList.item(0).getTextContent.toCharArray;
      config_name = config_name(:).';
      
      % Get the config associated with this ADC
      expression = xpath.compile(sprintf('//configs/config[(@type="%s" and name="%s")]',config_type,config_name));
      nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
      dac_wf_cfg = nodeList.item(0);
      
      % Find the longest possible record size
      if strcmpi(config_type,'dac-ad9129_0012_waveform')
        % TOHFSounder
        
        expression = xpath.compile('name');
        nodeList = expression.evaluate(dac_wf_cfg,XPathConstants.NODESET);
        name = nodeList.item(0);
        name = name.getTextContent.toCharArray;
        name = name(:).';
        settings.dac{mode_latch+1}.name = name;
        
        expression = xpath.compile('sampFreq');
        nodeList = expression.evaluate(dac_wf_cfg,XPathConstants.NODESET);
        sampFreq = nodeList.item(0);
        sampFreq = sampFreq.getTextContent.toCharArray;
        sampFreq = str2double(sampFreq(:).');
        settings.dac{mode_latch+1}.sampFreq = sampFreq;
        
        % Get each subchannel
        expression = xpath.compile('pulse');
        pulse_nodeList = expression.evaluate(dac_wf_cfg,XPathConstants.NODESET);
        for pulse_idx = 1:pulse_nodeList.getLength
          pulse_cfg = pulse_nodeList.item(pulse_idx-1);
          if isempty(pulse_cfg)
            continue;
          end
          
          expression = xpath.compile('name');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          name = nodeList.item(0);
          name = name.getTextContent.toCharArray;
          name = name(:).';
          settings.dac{mode_latch+1}.wfs{pulse_idx}.name = name;
          
          expression = xpath.compile('centerFreq');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          centerFreq = nodeList.item(0);
          centerFreq = centerFreq.getTextContent.toCharArray;
          centerFreq = str2double(centerFreq(:).');
          settings.dac{mode_latch+1}.wfs{pulse_idx}.centerFreq = centerFreq;
          
          expression = xpath.compile('bandwidth');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          bandwidth = nodeList.item(0);
          bandwidth = bandwidth.getTextContent.toCharArray;
          bandwidth = str2double(bandwidth(:).');
          settings.dac{mode_latch+1}.wfs{pulse_idx}.bandwidth = bandwidth;
          
          expression = xpath.compile('taper');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          taper = nodeList.item(0);
          taper = taper.getTextContent.toCharArray;
          taper = taper(:).';
          settings.dac{mode_latch+1}.wfs{pulse_idx}.taper = taper;
          
          expression = xpath.compile('alpha');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          alpha = nodeList.item(0);
          alpha = alpha.getTextContent.toCharArray;
          alpha = str2double(alpha(:).');
          settings.dac{mode_latch+1}.wfs{pulse_idx}.alpha = alpha;
          
          expression = xpath.compile('scale');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          scale = nodeList.item(0);
          scale = scale.getTextContent.toCharArray;
          scale = str2double(scale(:).');
          settings.dac{mode_latch+1}.wfs{pulse_idx}.scale = scale;
          
          expression = xpath.compile('numPoints');
          nodeList = expression.evaluate(pulse_cfg,XPathConstants.NODESET);
          numPoints = nodeList.item(0);
          numPoints = numPoints.getTextContent.toCharArray;
          numPoints = str2double(numPoints(:).');
          settings.dac{mode_latch+1}.wfs{pulse_idx}.numPoints = numPoints;
        end
      end
      
      
    end
  end
  
end


