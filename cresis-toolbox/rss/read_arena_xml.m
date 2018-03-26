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
  expression = xpath.compile('//config/@type');
  nodeList = expression.evaluate(match,XPathConstants.NODESET);
  config_type = nodeList.item(0).getTextContent.toCharArray;
  config_type = config_type(:).';
  expression = xpath.compile('//config/name');
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
        range_gates = nodeList.item(0);
        range_gates = range_gates.getTextContent.toCharArray;
        range_gates = range_gates(:).';
        % Assumes simple range gate format "start:stop"
        [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
        num_bins = str2double(stop) - str2double(start) + 1;
        if num_bins > settings.max_num_bins
          settings.max_num_bins = num_bins;
        end
        
        settings.hdrs{mode_latch+1,subchannel+1} = range_gates;
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
    
  else
    error('ADC type %s not supported.', config_type);
  end
end
