function MCMC_HMM_automated_frames(MCMC_HMM_param)
   
% MCMC_HMM_automated_frames(MCMC_HMM_param) 
%
%   Script to run a sequence of frames from the data .mat files and
%       call the MCMC (or HMM) algorithms.
%
%   Called from run_MCMC_HMM_automated_frames.m
%
%   See also: run_MCMC_HMM_automated_frames.m
%
%   Author: Victor Berger                              

% Algorithm parameters:
top_smooth=1000;
bottom_smooth=1000;
top_peak = 0.5;
bottom_peak = 0.5;
repulse = 10;
pts1=[];
pts2=[];

input_params = xlsread(MCMC_HMM_param.param_dir);
[input_rows, ~] = size(input_params);

mcmc = strcmp(MCMC_HMM_param.alg, 'MCMC');

for idx = 1:input_rows
    
    % Check which frames have a generic setting of '1'
    
    if (isnan(input_params(idx,1)) || isnan(input_params(idx,9)) ...
            || input_params(idx,9) ~= 1)
        continue
    end
    
    frame_name_num = input_params(idx,1);
    frame_name_str = num2str(frame_name_num);
    segment_num = input_params(idx,2);
    segment_num_str = num2str(segment_num);
    
    if MCMC_HMM_param.debug
        fprintf('\nMatch found: loop idx %i, frame %s, segment %s\n', idx, ...
            frame_name_str, segment_num_str);
    end
    
    if segment_num < 10
        segment_num_str = strcat('0', segment_num_str);
    end
    
    full_dir = strcat(MCMC_HMM_param.data_dir, '/', frame_name_str, '_', segment_num_str, '/');
    
    if ~exist(full_dir, 'dir')
        fprintf('Data folder not found.\nPath: %s\n', full_dir);
        continue;
    end
    
    dir_list_total = dir(full_dir);
    dir_list_match = {};
    name_prototype = strcat('Data_', frame_name_str, '_', segment_num_str);
    
    for i = 1:size(dir_list_total)
        if strncmp(dir_list_total(i).name, name_prototype, length(name_prototype))
            dir_list_match{end+1} = dir_list_total(i).name;
        end
    end
    
    if MCMC_HMM_param.png_output
        if(~isempty(dir_list_match))
            img_frame_dir_name = strcat(frame_name_str, '_', segment_num_str);
            mkdir (MCMC_HMM_param.img_output_dir, img_frame_dir_name);
        else continue
        end
    end
    
    for j = 1:length(dir_list_match)        
        
        fprintf('Running %s algorithm on frame %s\n', MCMC_HMM_param.alg, dir_list_match{j});
        
        fname = strcat(full_dir, dir_list_match{j});
        
        match_cat = dir_list_match{j};
        match_cat = match_cat(6:20);
        
        if MCMC_HMM_param.png_output
            outpath = strcat(MCMC_HMM_param.img_output_dir, '/', img_frame_dir_name, '/', match_cat, '.png');
        else outpath = '';
        end
        
        % Call run_ice to run desired algorithm
        
        [path, ~, ~] = run_ice(fname, top_smooth, bottom_smooth, ...
            top_peak, bottom_peak, repulse, pts1, pts2, mcmc, MCMC_HMM_param.display, MCMC_HMM_param.png_output, outpath);
        
        if MCMC_HMM_param.save_layer_data
            
            % Save surface and bottom paths to layerData file
            
            layerData_file_in = strcat(MCMC_HMM_param.original_layer_dir, '/', frame_name_str, '_', ...
                segment_num_str, '/Data_', match_cat);
            layerData_struct = load(layerData_file_in);
            
            [~, data_size_alg] = size(path(1,:));
            [~, data_size_orig] = size(layerData_struct.layerData{1}.value{1}.data);
            
            if MCMC_HMM_param.debug
                fprintf('Size of original layer: %i, size of vector returned by algorithm: %i\n', ...
                    data_size_orig, data_size_alg);
            end
            
            data_dir_full = strcat(MCMC_HMM_param.data_dir, '/', frame_name_str, '_', segment_num_str, ...
                '/', dir_list_match{j});
            
            echo_data = load(data_dir_full, 'GPS_time', 'Time');
            
            if(data_size_alg ~= data_size_orig)
                
                % Interpolate to match GPS time of layer data
                % with GPS time of data file
                
                fprintf('Running interpolation on vector resulting from %s algorithm...\n', MCMC_HMM_param.alg);
                
                layerData_struct.layerData{1}.value{2}.data = ...
                    interp1(echo_data.GPS_time, path(1,:), layerData_struct.GPS_time);
                
                layerData_struct.layerData{2}.value{2}.data = ...
                    interp1(echo_data.GPS_time, path(2,:), layerData_struct.GPS_time);                
                
            else
                layerData_struct.layerData{1}.value{2}.data = path(1,:);
                layerData_struct.layerData{2}.value{2}.data = path(2,:);
            end
            
            % Interpolate back from row number to TWTT
            
            layerData_struct.layerData{1}.value{2}.data = ...
                interp1(1:length(echo_data.Time), echo_data.Time, layerData_struct.layerData{1}.value{2}.data);
            
            layerData_struct.layerData{2}.value{2}.data = ...
                interp1(1:length(echo_data.Time), echo_data.Time, layerData_struct.layerData{2}.value{2}.data);
                       
            % Generate output file name and save layer data
            
            layerData_file_out = strcat(MCMC_HMM_param.output_layer_dir, '/', frame_name_str, '_', ...
                segment_num_str, '/Data_', match_cat, '.mat');
            save(layerData_file_out, '-struct', 'layerData_struct');
            
            clear path;
        end
    end
end
end
