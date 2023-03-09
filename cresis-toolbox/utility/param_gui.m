function params = param_gui

% an utility gui to select params
% Author: Hara Madhav Talasila
%

global gRadar;
N_irritation  = 3;
fns_param     = dir(fullfile(gRadar.param_path, '*.xls'));
fns_param     = {fns_param.name};
fn_param      = [];

try
  user_pref = load(fullfile(gRadar.tmp_path,'user_pref.mat'));
  param_init_value = user_pref.param_init_value;
  if ~isempty(param_init_value)
    param_init_value = (param_init_value<=length(fns_param)) * param_init_value;
  end
catch
  param_init_value = 1;
end

for patience = 1:N_irritation
  
  [idx,okay] = listdlg('ListString', fns_param,...%     'SelectionMode', 'single', ...
    'ListSize', [420,902], ...
    'InitialValue', param_init_value, ...
    'Name', sprintf('%s param_gui', mfilename), ...
    'PromptString', {'Which radar_param_YYYY_Campaign_Platform.xls ???'}, ...
    'OKString', 'Enter', ...
    'CancelString', sprintf('You will be back! %d',N_irritation-patience));
  
  if okay
    fn_param = fullfile(gRadar.param_path,fns_param{idx});
    try
      user_pref.param_init_value = idx;
      ct_save(fullfile(gRadar.tmp_path,'user_pref.mat'), '-struct', 'user_pref');
    end
    sheets    = sheet_gui(fn_param);
    day_segs  = day_segs_gui(fn_param);
    params = read_param_xls(fn_param,day_segs,sheets);
    break;
  else
    warning('Only way forward is through!');
  end
  
end

end


%%

function sheets = sheet_gui(fn_param)

% an utility gui to select (generic) sheets
% Author: Hara Madhav Talasila
%

global gRadar;
N_irritation  = 3;

sheets = setdiff(sheetnames(fn_param), ...
  ["cmd", "records", "qlook", "sar", "array", "radar", "post"] );

try
  user_pref = load(fullfile(gRadar.tmp_path,'user_pref.mat'));
  sheet_init_value = user_pref.sheet_init_value;
  if ~isempty(sheet_init_value)
    sheet_init_value = sheet_init_value(find(sheet_init_value<=length(sheets)));
  end
catch
  sheet_init_value = 1;
end

for patience = 1:N_irritation
  
  [idx,okay] = listdlg('ListString', sheets,...
    'SelectionMode', 'multiple', ...
    'ListSize', [420,369], ...
    'InitialValue', sheet_init_value, ...
    'Name', sprintf('%s sheet_gui', mfilename), ...
    'PromptString', {'Extra cheese ???'}, ...
    'OKString', 'Enter', ...
    'CancelString', sprintf('No, thanKU! %d',N_irritation-patience));
  
  if okay || patience == N_irritation
    sheets = sheets(idx);
    if size(sheets,1)>1
      if 1
        for she = 1:size(sheets,1)
          sheesh{she,1} = sheets(she);
        end
        sheets = sheesh;
      else
        sheets = {sheets};
      end
    end
    try
      user_pref.sheet_init_value = idx;
      ct_save(fullfile(gRadar.tmp_path,'user_pref.mat'), '-struct', 'user_pref');
    end
    break;
  else
    warning('Only way forward is through!');
  end
  
end

end

%%

function day_segs  = day_segs_gui(fn_param)

% an utility gui to select day_segs
% Author: Hara Madhav Talasila
%

global gRadar;
N_irritation  = 3;

[numbers, txtt] = xlsread(fn_param,"cmd");
mission_names_column = find(~cellfun(@isempty,regexp(txtt(4,:), 'mission_names')));
if isempty(mission_names_column)
  mission_names_en = 0;
else
  mission_names_en = 1;
end

N_day_segs = size(numbers,1) - 5; % day_segs start at 6
for ds_idx = 1:N_day_segs
  day_segs{ds_idx} = sprintf('%d_%02d', numbers(ds_idx+5,1), numbers(ds_idx+5,2) );
  if mission_names_en
    if ~isempty(regexp(txtt{ds_idx+5, mission_names_column+1}, 'do not process'))
      day_segs_display{ds_idx} = sprintf('%d_%02d  %s DO NOT PROCESS', numbers(ds_idx+5,1), numbers(ds_idx+5,2), txtt{ds_idx+5, mission_names_column} );
    else
      day_segs_display{ds_idx} = sprintf('%d_%02d %s', numbers(ds_idx+5,1), numbers(ds_idx+5,2), txtt{ds_idx+5, mission_names_column} );
    end
  end
end

if ~mission_names_en
  day_segs_display = day_segs;
end
  
try
  user_pref = load(fullfile(gRadar.tmp_path,'user_pref.mat'));
  day_segs_init_value = user_pref.day_segs_init_value;
  if ~isempty(day_segs_init_value)
    day_segs_init_value = day_segs_init_value(find(day_segs_init_value<=length(day_segs)));
  end
catch
  day_segs_init_value = 1;
end

for patience = 1:N_irritation
  
  [idx,okay] = listdlg('ListString', day_segs_display,...
    'SelectionMode', 'multiple', ...
    'ListSize', [420,902], ...
    'InitialValue', day_segs_init_value, ...
    'Name', sprintf('%s day_segs_gui', mfilename), ...
    'PromptString', {'Which day_segs ???'}, ...
    'OKString', 'Enter', ...
    'CancelString', sprintf('You will be back! %d',N_irritation-patience));
  
  if okay
    day_segs = day_segs(idx);
    try
      user_pref.day_segs_init_value = idx;
      ct_save(fullfile(gRadar.tmp_path,'user_pref.mat'), '-struct', 'user_pref');
    end
    break;
  else
    warning('Only way forward is through!');
  end
  
end

end

