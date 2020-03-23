% script run_load_sar_data
%
% Load SAR processed data (e.g. output of csarp stage)
%
% Author: John Paden, Logan Smith

warning('This is an example file, copy to personal directory and remove this warning/return to use');
%return;

% =======================================================================
% User Settings
% =======================================================================

param = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'20191231_04');

% param.load_sar_data.fn = ''; % Leave empty for default 'sar'

% Start and stop chunk to load (inf for second element loads to the end)
% param.load_sar_data.chunk = []; % Leave empty for all chunks

% param.load_sar_data.sar_type = ''; % Leave empty for default 'fk'

param.load_sar_data.frm = 1; % Specify data frame to load

% param.load_sar_data.subap = []; % Leave empty for default (all subapertures)

% {Images} with [wf,adc] pairs to load
% param.load_sar_data.imgs = {[1 1; 1 2; 1 3; 1 4; 1 5; 1 6]};
% param.load_sar_data.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14]};
param.load_sar_data.imgs = {[1*ones(8,1) (1:8)'],[2*ones(8,1) (1:8)'],[3*ones(8,1) (1:8)']};

% Debug level (1 = default)
param.load_sar_data.debug_level = 2;

% Combine receive channels
% param.load_sar_data.combine_channels = 0;

% Combine waveforms parameters
% param.load_sar_data.wf_comb = 10e-6;

% Take abs()^2 of the data (only runs if combine_channels runs)
% param.load_sar_data.incoherent = 0;

% Combine waveforms (only runs if incoherent runs)
% param.load_sar_data.combine_imgs = 0;

% Parameters for local_detrend (cmd == 5 disables, only runs if incoherent runs)
% param.load_sar_data.detrend.cmd = 3;
% param.load_sar_data.detrend.B_noise = [100 200];
% param.load_sar_data.detrend.B_sig = [1 10];
% param.load_sar_data.detrend.minVal = -inf;

[data,metadata] = load_sar_data(param);

return


%% Examples of manipulating multilook data
good_mask = metadata.wfs(1).time > 0e-6 & metadata.wfs(1).time < 30e-6;

% Step through multilooks
subap = 0;
figure(1); clf;
while 1
  subap = mod(subap,size(data{1},4))+1;
  imagesc([],metadata.wfs(1).time(good_mask),lp(mean(data{1}(good_mask,:,:,subap),3)),[-150 -23])
  title(sprintf('%.0f',subap));
  if subap == 1
    beep
  end
  pause; % pause for keystroke advance or pause(0.1) for animation
end


figure(2); clf;
imagesc([],metadata.wfs(1).time(good_mask),lp(local_detrend(mean(abs(mean(data{1}(good_mask,:,:,[6:16]),3)).^2,4),[40 200],[1 3])))

figure(2); clf;
imagesc([],metadata.wfs(1).time(good_mask),lp(mean(abs(mean(data{1}(good_mask,:,:,[6:16]),3)).^2,4)),[-150 -23])

figure(2); clf;
imagesc([],metadata.wfs(1).time(good_mask),lp(mean(abs(mean(data{1}(good_mask,:,:,[1:6]),3)).^2,4)),[-150 -23])

figure(3); clf;
imagesc([],metadata.wfs(1).time(good_mask),lp(mean(abs(mean(data{1}(good_mask,:,:,[end-5:end]),3)).^2,4)),[-150 -23])







