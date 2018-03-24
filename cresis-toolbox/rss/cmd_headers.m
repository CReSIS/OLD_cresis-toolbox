%% Print cmd headers
fprintf('\ncmd WorkSheet\n');
fprintf('======Copy lines below and paste in cmd sheet===========\n\n')
fprintf('Version\t''%s\n',version_num)
fprintf('Radar\t%s\n',radar_name)
fprintf('Season\t%s\n',season_name)
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
'Date','','frms','create_records','create_frames','get_heights',...
'csarp','combine_wf_chan','generic','mission_names','notes');
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
'YYYYMMDD','Segment','r','b','b','b','b','b','r','t','t');