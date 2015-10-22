% script basic_check_mcords_data
%
% Checks EPRI and time variables from the header of each ADC to see if
% the data is okay. Certain MCoRDS settings cause the digital system
% to be unstable and this code checks for that. You should have no
% EPRI or timing jumps in a "good" state (occasionally a "good" state
% does have 1 or 2 errors, but the majority of files should be clean).
%
% Author: John Paden

fprintf('Checking mcords data\n');

adcs = [3];
base_path = '/cresis/data3/MCoRDS/2010_Chile/20101106/chan1/seg2/';
base_path = 'D:\data\chan3\seg40\';
radar_num = 1;
file_num = 0;
for adc = adcs
  fn_midfix = sprintf('r%d-%d',radar_num,adc);
  fn_suffix = sprintf('.%04d.dat',file_num);
  fprintf('  Path: %s\n', base_path);
  fprintf('  Match: mcords*%s*%s\n', fn_midfix, fn_suffix);
  fn = get_filename(base_path,'mcords',fn_midfix,fn_suffix);
  if isempty(fn)
    fprintf('  Could not find any files which match\n');
    return;
  end
  fprintf('  Loading file %s\n', fn);
  if file_num == 0
    [hdr,data] = basic_load_mcords(fn,struct('clk',1e9/9,'first_byte',2^26));
  else
    [hdr,data] = basic_load_mcords(fn,struct('clk',1e9/9));
  end
  finfo = frame_sync_info(fn);
    
  figure(adc); clf;
  plot(diff(finfo.syncs));
  title(sprintf('Record size should be constant at %d',hdr.finfo.rec_size));
  xlabel('Record');
  ylabel('Records Size');
  figure(10+adc); clf;
  plot(diff(hdr.epri));
  ylim([0 4]);
  title('EPRI spacing should be constant at 1');
  xlabel('Record');
  ylabel('Diff EPRI');
  figure(20+adc); clf;
  plot(diff(hdr.time_sod),'.');
  ylim([0 0.02]);
  title('Time seconds of day spacing should be almost constant');
end

