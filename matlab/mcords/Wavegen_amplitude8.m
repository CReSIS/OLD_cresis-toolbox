% This routine is for the 8 channel waveform generator
% designed for UAV2008. It is used by all MCoRDS platforms
% to generate the custom waveshape files.
%
% Carl Leuschen did the scaling up, 
% John Ledford did the rest.

% Clear up everything:
fprintf('=================================================\n');
fprintf('Wavegen\n');
clear;
close all;    %close all open figures

% -------------------------------------------------------
% User Parameters
% -------------------------------------------------------

% Choose one pulse duration setting
if 0
  % 10 us setting
  addLossDDS(1) = 11.7; % Additional loss from maximum output, in dB
  addLossDDS(2) = 5.82; % 1e6 for value will drive output to zero
  addLossDDS(3) = 3; % These should never be negative!
  addLossDDS(4) = 9.605;
  addLossDDS(5) = 17.72;
  addLossDDS(6) = 30;
  addLossDDS(7) = 30;
  addLossDDS(8) = 30;
  % Set up the waveform parameters (1d0 = 1.0, 2d5 = 2.5, etc)
  outfile='Paden_1d0us_20percentTukey_10usPulse.wf';
  numPoints = 222; % Get from the waveform generator feedback
  % for DC8 setup:
  %  222 = 1us   Basepulse
  %  555 = 2.5us Basepulse
  %  666 = 3us   Basepulse
  
  window_type = 3;    % 1 = rectangular
  % 2 = hanning
  % 3 = tukey 20%
elseif 0
  % 20 us setting
  addLossDDS(1) = 8.55; % Additional loss from maximum output, in dB
  addLossDDS(2) = 6.36; % 1e6 for value will drive output to zero
  addLossDDS(3) = 3; % These should never be negative!
  addLossDDS(4) = 9.6425;
  addLossDDS(5) = 10.8575;
  addLossDDS(6) = 30;
  addLossDDS(7) = 30;
  addLossDDS(8) = 30;
  % Set up the waveform parameters (1d0 = 1.0, 2d5 = 2.5, etc)
  outfile='Paden_1d0us_10percentTukey_20usPulse.wf';
  numPoints = 222; % Get from the waveform generator feedback
  % for DC8 setup:
  %  222 = 1us   Basepulse
  %  555 = 2.5us Basepulse
  %  666 = 3us   Basepulse
  
  window_type = 4;    % 1 = rectangular
  % 2 = hanning
  % 3 = tukey 20%
  % 4 = tukey 10%
elseif 0
  % 3 us setting
  addLossDDS(1) = 11.4; % Additional loss from maximum output, in dB
  addLossDDS(2) = 5.64; % 1e6 for value will drive output to zero
  addLossDDS(3) = 3; % These should never be negative!
  addLossDDS(4) = 9.98;
  addLossDDS(5) = 18.19;
  addLossDDS(6) = 30;
  addLossDDS(7) = 30;
  addLossDDS(8) = 30;
  % Set up the waveform parameters (1d0 = 1.0, 2d5 = 2.5, etc)
  outfile='Paden_1d0us_20percentTukey_3usPulse.wf';
  numPoints = 222; % Get from the waveform generator feedback
  % for DC8 setup:
  %  222 = 1us   Basepulse
  %  555 = 2.5us Basepulse
  %  666 = 3us   Basepulse
  
  window_type = 3;    % 1 = rectangular
  % 2 = hanning
  % 3 = tukey 20%
elseif 0
  % 30 us setting
  addLossDDS(1) = 8.55; % Additional loss from maximum output, in dB
  addLossDDS(2) = 6.36; % 1e6 for value will drive output to zero
  addLossDDS(3) = 3; % These should never be negative!
  addLossDDS(4) = 9.6425;
  addLossDDS(5) = 10.8575;
  addLossDDS(6) = 30;
  addLossDDS(7) = 30;
  addLossDDS(8) = 30;
  % Set up the waveform parameters (1d0 = 1.0, 2d5 = 2.5, etc)
  outfile='Paden_2d5us_10percentTukey_30usPulse.wf';
  numPoints = 555; % Get from the waveform generator feedback
  % for DC8 setup:
  %  222 = 1us   Basepulse
  %  555 = 2.5us Basepulse
  %  666 = 3us   Basepulse
  
  window_type = 4;    % 1 = rectangular
  % 2 = hanning
  % 3 = tukey 20%
  % 4 = tukey 10%
elseif 1
  % Rx Lab Gain Tests
  addLossDDS(1) = 11.7; % Additional loss from maximum output, in dB
  addLossDDS(2) = 5.82; % 1e6 for value will drive output to zero
  addLossDDS(3) = 3; % These should never be negative!
  addLossDDS(4) = 9.605;
  addLossDDS(5) = 17.72;
  addLossDDS(6) = 30;
  addLossDDS(7) = 30;
  addLossDDS(8) = 30;
  % Set up the waveform parameters (1d0 = 1.0, 2d5 = 2.5, etc)
  outfile='Paden_1d0us_boxcar_10usPulse.wf';
  numPoints = 222; % Get from the waveform generator feedback
  % for DC8 setup:
  %  222 = 1us   Basepulse
  %  555 = 2.5us Basepulse
  %  666 = 3us   Basepulse
  
  window_type = 1;    % 1 = rectangular
  % 2 = hanning
  % 3 = tukey 20%
end

% -------------------------------------------------------
% Don't touch the stuff below unless you really know what your doing
% -------------------------------------------------------

if any(addLossDDS < 0)
  error('addLossDDS must all be positive')
end

numDDS = 8;  % Just in case I want to make it expandable
    
    
if  (window_type == 1) % rectangular
    template = ones(numPoints,1);
end

if  (window_type == 2) % Hanning
    template = hanning(numPoints);
end

if  (window_type == 3) % Tukey 20%
    template = tukeywin(numPoints,0.2);
end

if  (window_type == 4) % Tukey 10%
    template = tukeywin(numPoints,0.1);
end

% Set last point = 0, to guarantee that the signal is off
% during the off period. The problem is that the DDS chip flips the 
% signal, which is corrected for in the output exporting phase, but 
% for now it means we actually set the first sample to zero.
template(1)=0;

% Scale up - By Carl Leuschen
template = template*((2^32)-1);
template = (2^18)*floor(template/(2^18));

% plot(template);


% Now make the individual amplitudes

i=1;
while (i <= numDDS)
    fprintf('Calculating DDS[%d] Shape\n',i)
    wfAmp{i} = template/(10^(addLossDDS(i)/20));
    i=i+1;
end


figure(1); clf;
plot(wfAmp{1},'b'); hold on;
plot(wfAmp{2},'g');
plot(wfAmp{3},'r');
plot(wfAmp{4},'c');
plot(wfAmp{5},'m');
plot(wfAmp{6},'y');
plot(wfAmp{7},'k');
plot(wfAmp{8},'b');




% Open file (force little-endian formating == PC byte ordering)
fid = fopen(outfile,'w','ieee-be');

% Write waveshape points, reversing the order since the DDS plays them
% "backwards" (in John's opinion).
i=1;
while (i <= numDDS)
    for index = 1:length(wfAmp{i})
        fwrite(fid,wfAmp{i}((length(wfAmp{i})+1)-index),'uint32');
    end
    i=i+1;
end



% Close file
fclose(fid);

return;
