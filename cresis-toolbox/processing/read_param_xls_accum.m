function [params] = read_param_xls_accum(param_fn, day_seg_filter)
% [params] = read_param_xls_accum(param_fn, day_seg_filter)
%
% Support function for read_param_xls (for accumulation radar).
%
% Author: Brady Maasen, John Paden
%
% See also: read_param_xls

cell_boolean = @read_param_xls_boolean;
cell_text = @read_param_xls_text;
cell_read = @read_param_xls_general;

% ======================================================================%
% CREATING THE PARAM STRUCTURE ARRAY FROM PARAM_STARTER.XLS
% ======================================================================%
warning('off','MATLAB:xlsread:Mode');

% =======================================================================
% Create Command Parameters
% =======================================================================
sheet_name = 'command';
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);
[num txt] = xlsread(param_fn,sheet_name,'','basic');

param_file_version        = cell_text(1,2,num,txt);
radar_name                = cell_text(2,2,num,txt);
season_name              = cell_text(3,2,num,txt);
sw_version                = current_software_version;

num_header_rows = 5;
rows = max(size(num,1), size(txt,1)) - num_header_rows;

for idx = 1:rows
  params(idx).radar_name                        = radar_name;
  params(idx).season_name                      = season_name;
  row = idx + num_header_rows;
  params(idx).day_seg                           = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
  col = 3;
  params(idx).cmd.frms                          = cell_read(row,col,num,txt); col = col + 1;
  params(idx).cmd.create_vectors                = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).cmd.create_records                = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).cmd.qlook                         = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).cmd.post                          = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).cmd.mission_name                  = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.sw_version                   = sw_version;
end


% =======================================================================
% Create Vectors Parameters
% =======================================================================
sheet_name = 'vectors';
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);
[num txt] = xlsread(param_fn,sheet_name,'','basic');

num_header_rows = 2;
rows = max(size(num,1), size(txt,1)) - num_header_rows;

for idx = 1:rows
  row = idx + num_header_rows;
  day_seg = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
  if ~strcmpi(params(idx).day_seg,day_seg)
    error('The order of sheets cmd and %s do not match at row %d in the excel spreadsheet. Each sheet must have the same date and segment list.',sheet_name,row);
  end
  col = 3;
  params(idx).vectors.file.start_idx            = cell_read(row,col,num,txt); col = col + 1;
  params(idx).vectors.file.stop_idx             = cell_read(row,col,num,txt); col = col + 1;
  params(idx).vectors.file.base_dir             = cell_text(row,col,num,txt); col = col + 1;
  params(idx).vectors.file.adc_folder_name      = cell_text(row,col,num,txt); col = col + 1;
  params(idx).vectors.file.file_prefix          = cell_text(row,col,num,txt); col = col + 1;
  params(idx).vectors.out_fn                    = cell_text(row,col,num,txt); col = col + 1;
  params(idx).vectors.gps.fn                    = cell_text(row,col,num,txt); col = col + 1;
  params(idx).vectors.gps.time_offset           = cell_read(row,col,num,txt); col = col + 1;
  params(idx).vectors.gps.utc_time_halved       = cell_read(row,col,num,txt); col = col + 1;
  params(idx).vectors.verification              = cell_text(row,col,num,txt); col = col + 1;
end

% =======================================================================
% QLook Parameters
% =======================================================================
sheet_name = 'qlook';
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);
[num txt] = xlsread(param_fn,sheet_name,'','basic');

num_header_rows = 2;
rows = max(size(num,1), size(txt,1)) - num_header_rows;

for idx = 1:rows
  row = idx + num_header_rows;
  day_seg = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
  if ~strcmpi(params(idx).day_seg,day_seg)
    error('The order of sheets cmd and %s do not match at row %d in the excel spreadsheet. Each sheet must have the same date and segment list.',sheet_name,row);
  end
  col = 3;
  
  params(idx).qlook.out.dir                   = cell_text(row,col,num,txt); col = col + 1;
  params(idx).qlook.out.en                    = cell_boolean(row,col,num,txt); col = col + 1;
  
  params(idx).qlook.gps.en                    = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).qlook.records_en                = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).qlook.gps.fn                    = params(idx).vectors.gps.fn;
  params(idx).qlook.gps.time_offset           = params(idx).vectors.gps.time_offset;
  
  params(idx).qlook.file.start_idx            = params(idx).vectors.file.start_idx;
  params(idx).qlook.file.stop_idx             = params(idx).vectors.file.stop_idx;
  params(idx).qlook.file.base_dir             = params(idx).vectors.file.base_dir;
  params(idx).qlook.file.adc_folder_name      = params(idx).vectors.file.adc_folder_name;
  params(idx).qlook.file.file_prefix          = params(idx).vectors.file.file_prefix;
  
  params(idx).qlook.tukey                     = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.band_window_func          = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.td_window_func            = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.window_func               = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.sw_presums                = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.incoh_ave                 = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.decimate                  = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.plot.en                   = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.elev_comp.en              = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.detrend_poly_order        = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.surf.min_bin              = cell_read(row,col,num,txt); col = col + 1;
  params(idx).qlook.surf.search_rng           = cell_read(row,col,num,txt); col = col + 1;
  
  params(idx).qlook.lever_arm_fh              = cell_read(row,col,num,txt); col = col + 1;

  params(idx).qlook.sw_version                = current_software_version;
end

% =======================================================================
% Posting Parameters
% =======================================================================
sheet_name = 'post';
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);
[num txt] = xlsread(param_fn,sheet_name,'','basic');

num_header_rows = 2;
rows = max(size(num,1), size(txt,1)) - num_header_rows;

for idx = 1:rows
  row = idx + num_header_rows;
  day_seg = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
  if ~strcmpi(params(idx).day_seg,day_seg)
    error('The order of sheets cmd and %s do not match at row %d in the excel spreadsheet. Each sheet must have the same date and segment list.',sheet_name,row);
  end
  col = 3;
  
  params(idx).post.in_dir                     = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.out_dir                    = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.gps.en                     = cell_boolean(row,col,num,txt); col = col + 1;
  params(idx).post.gps.fn                     = params(idx).vectors.gps.fn;
  params(idx).post.gps.time_offset            = params(idx).vectors.gps.time_offset;

  params(idx).post.num_frm_combine            = cell_read(row,col,num,txt); col = col + 1;
  params(idx).post.depth_rng                   = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.map.en                     = cell_read(row,col,num,txt); col = col + 1;
    
  params(idx).post.map.type                   = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.map.vectors_fn             = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.map.location               = cell_text(row,col,num,txt); col = col + 1;
  
  params(idx).post.er_ice                     = cell_read(row,col,num,txt); col = col + 1;
  params(idx).post.mission_name               = params(idx).cmd.mission_name;
  params(idx).post.img_type                   = cell_text(row,col,num,txt); col = col + 1;
  params(idx).post.img_dpi                    = cell_read(row,col,num,txt); col = col + 1;
    
  params(idx).post.plot_params                = cell_read(row,col,num,txt); col = col + 1;

  params(idx).post.sw_version                 = current_software_version;
end

% =======================================================================
% Radar and general processing parameters
% =======================================================================
sheet_name = 'radar';
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);
[num txt] = xlsread(param_fn,sheet_name,'','basic');

num_header_rows = 2;
rows = max(size(num,1), size(txt,1)) - num_header_rows;

for idx = 1:rows
  row = idx + num_header_rows;
  day_seg = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
  if ~strcmpi(params(idx).day_seg,day_seg)
    error('The order of sheets cmd and %s do not match at row %d in the excel spreadsheet. Each sheet must have the same date and segment list.',sheet_name,row);
  end
  col = 3;
  
  params(idx).radar.wfs         = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.prf         = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.fs          = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.f0          = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.fLO         = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.f_step      = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.BW          = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.Tpd         = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.adc_bits    = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.Vpp_scale   = cell_read(row,col,num,txt); col = col + 1;
  params(idx).radar.td          = cell_read(row,col,num,txt); col = col + 1;
  
  num_rx_wfs                  = 16;
  
  % Channel equalization, each channel data DIVIDED by this
  chan_equal_mag = [];
  for wf = 1:num_rx_wfs
    % Power
    chan_equal_mag(wf) = cell_read(row,col,num,txt); col = col + 1;
    chan_equal_ang(wf) = 0;
  end
  for wf = 1:num_rx_wfs
    params(idx).radar.chan_equal(wf) ...
      = 10^(chan_equal_mag(wf)/20) * exp(1j*chan_equal_ang(wf)/180*pi);
  end
end

return

