function viterbi_segs (params, options)
% function viterbi_segs (params, options)
%
% Builds a large matrix containing multiple segments of 2D data and applies
%  the Viterbi layer-tracking program on it. 
% Loads crossovers from OPS to use as ground truth for the tracking.
% Loads labelled bottom layer into the OPS.
%
% See also: run_viterbi_segs.m
%
% Authors: John Paden, Victor Berger


%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar;

physical_constants;
clear('param_override');

% Input checking
if ~exist('params','var')
    error('Use run_master: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
else
    param_override = gRadar;
end

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
    param = params(param_idx);
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
        continue;
    end
    
    param = merge_structs(param,param_override);
    
    % Automated section
    big_matrix          = {};
    big_matrix.Data     = [];
    big_matrix.Surface  = [];
    big_matrix.GPS_time = [];
    
    bottom_bin    = -1;
    egt_weight    = -1;
    mu_size       = 31;
    mu            = log10(exp(-(-(mu_size-1)/2 : (mu_size-1)/2).^4/1));
    mu(mu<-30)    = -30;
    mu            = mu - mean(mu);
    sigma         = sum(abs(mu))/10*ones(1,mu_size);
    smooth_weight = 1;
    smooth_var    = inf;
    repulsion     = 150000;
    ice_bin_thr   = 10;
    
    % Load frames file
    load(ct_filename_support(param,param.records.frames_fn,'frames'));
    
    if isempty(param.cmd.frms)
        param.cmd.frms = 1:length(frames.frame_idxs);
    end
    
    % Remove frames that do not exist from param.cmd.frms list
    [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
    if length(valid_frms) ~= length(param.cmd.frms)
        bad_mask = ones(size(param.cmd.frms));
        bad_mask(keep_idxs) = 0;
        warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
            param.cmd.frms(find(bad_mask,1)));
        param.cmd.frms = valid_frms;
    end
    
    data_fn_dir = ct_filename_out(param, options.name, '');
    for frm = param.cmd.frms
        data_fn_name        = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
        data_fn             = fullfile(data_fn_dir, data_fn_name );
        cur_matrix          = load(data_fn);
        big_matrix.Data     = horzcat(big_matrix.Data, cur_matrix.Data);
        big_matrix.Surface  = horzcat(big_matrix.Surface, cur_matrix.Surface);
        big_matrix.GPS_time = horzcat(big_matrix.GPS_time, cur_matrix.GPS_time);
        big_matrix.Time     = cur_matrix.Time;
    end
    
    Nx             = size(big_matrix.Data, 2);
    mask           = ones([1 Nx]);
    slope          = round(diff(big_matrix.Surface));
    viterbi_weight = ones([1 Nx]);
    
    opsAuthenticate(param,false);
    layer_name                   = 'bottom';
    sys                          = ct_output_dir(param.radar_name);
    ops_param                    = struct('properties',[]);
    ops_param.properties.season  = param.season_name;
    ops_param.properties.segment = param.day_seg;
    [status,ops_frames]          = opsGetSegmentInfo(sys,ops_param);
    
    query = sprintf('SELECT rds_segments.id FROM rds_seasons,rds_segments where rds_seasons.name=''%s'' and rds_seasons.id=rds_segments.season_id and rds_segments.name=''%s''',param.season_name,param.day_seg);
    [status,tables] = opsQuery(query);
    segment_id = tables{1};
    
    ops_param                       = struct('properties',[]);
    ops_param.properties.location   = param.post.ops.location;
    ops_param.properties.lyr_name   = layer_name;
    ops_param.properties.frame      = ops_frames.properties.frame;
    ops_param.properties.segment_id = ones(size(ops_param.properties.frame)) ...
        *double(segment_id);
    
    [status,data] = opsGetCrossovers(sys,ops_param);
    
    ops_param = struct('properties',[]);
    ops_param.properties.location = param.post.ops.location;
    ops_param.properties.season = param.season_name;
    ops_param.properties.start_gps_time = ops_frames.properties.start_gps_time(1);
    ops_param.properties.stop_gps_time = ops_frames.properties.stop_gps_time(end);
    ops_param.properties.nativeGeom = true;
    [~,ops_data] = opsGetPath(sys,ops_param);
    
    %%
    rline = [];
    rows  = [];
    cols  = [];
    gps_time = [];
    season_name = {};
    for i = 1 : length(data.properties.source_point_path_id)
        if ~isnan(data.properties.twtt(i))
            new_rline = find(ops_data.properties.id == data.properties.source_point_path_id(i));
            new_gps_time = ops_data.properties.gps_time(new_rline);
            new_season_name = data.properties.season_name{i};
            if big_matrix.GPS_time(1) <= new_gps_time ...
                    && big_matrix.GPS_time(end) >= new_gps_time ...
                    && str2double(new_season_name(1:4)) >= 2006
                rline(end+1) = new_rline;
                gps_time(end+1) = new_gps_time;
                season_name{end+1} = new_season_name;
                [~, cols(end+1)] = min(abs(big_matrix.GPS_time - ops_data.properties.gps_time(rline(end))));
                twtt = data.properties.twtt(i);
                twtt = twtt + (data.properties.source_elev(i)-data.properties.cross_elev(i))/(c/2);
                rows(end+1) = round(interp1(big_matrix.Time, 1:length(big_matrix.Time), twtt));
            end
        end
    end
    
    surf_bins = round(interp1(big_matrix.Time, 1:length(big_matrix.Time), big_matrix.Surface));
    
    
    %%
    figure(1);
    clf;
    imagesc(lp(big_matrix.Data));
    hold on;
    plot(cols, rows,'kx','LineWidth', 2);
    plot(surf_bins);
    hold off;
    
    % Set variable echogram tracking parameters
    gt     = [cols(:).'; rows(:).'];
    bounds = [0 Nx];
    big_matrix.Data = lp(big_matrix.Data);

    if options.compile
        clc
        tmp = pwd;
        cd ~/scripts/cresis-toolbox/cresis-toolbox/+tomo/
        mex -largeArrayDims viterbi.cpp
        cd(tmp);
    end

    % Call viterbi.cpp
    tic
    labels = tomo.viterbi(double(big_matrix.Data), ...
        double(surf_bins), double(bottom_bin), ...
        double(gt), double(mask), ...
        double(mu), double(sigma), double(egt_weight), ...
        double(smooth_weight), double(smooth_var), double(slope), ...
        int64(bounds), double(viterbi_weight), ...
        double(repulsion), double(ice_bin_thr));
    toc
    
    figure(2);
    clf;
    imagesc(big_matrix.Data);
    hold on;
    plot(1:Nx, labels, 'Color', 'red'); drawnow;
end