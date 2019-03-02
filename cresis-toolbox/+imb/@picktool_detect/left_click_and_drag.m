function cmds = left_click_and_drag(obj,param)
% cmds = left_click_and_drag(obj,param)
%
% Detect tool
%
% Compile with
%   mex -largeArrayDims viterbi.cpp

image_x = param.image_x;
image_y = param.image_y;
image_c = param.image_c;
cur_layers = param.cur_layers;
x = param.x;
y = param.y;
cmds = [];

fprintf('Detecting (HMM) points %f to %f, %f to %f\n', x, y);

param.x_bounds = 3;
param.y_bounds = 1;

tool_idx = get(obj.top_panel.tool_PM,'Value');
if tool_idx == 1
    
    %=========================================================================
    
    for layer_idx = 1:length(cur_layers)
        cur_layer = cur_layers(layer_idx);
        
        [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer);
        
        if length(manual_idxs) < 1
            warning('Insufficient points to detect');
            continue;
        elseif ~isempty(auto_idxs)
            
            % Nx: number of along track records/range lines
            Nx = length(image_x);
            custom_data.mu = [11.2575 11.3748 11.4393 11.4555 11.4323   11.3666   11.2668   11.1332   10.9900 10.8484   10.6916];
            custom_data.sigma = [5.4171    5.2945    5.2187    5.1939    5.2174    5.3247    5.4643    5.6571    5.8428 6.0477    6.2935];
            
            % Interpolate surface layer to match image x-axis coordinates
            surf_bins = interp1(param.layer.x,param.layer.y{1},image_x);
            % Interpolate surface layer y-axis units to image pixels
            surf_bins = interp1(image_y, 1:length(image_y),surf_bins);
            % Interpolate all non-finite values using surrounding data
            surf_bins = interp_finite(surf_bins);
            
            % Match GT points with axis coordinates
            gt = [param.layer.x(manual_idxs); interp1(image_y, 1:length(image_y),param.layer.y{cur_layer}(manual_idxs))];
            
            % Echogram Parameters
            detect_data    = image_c;
            bottom_bin     = -1;
            mask           = ones([1 Nx]);
            egt_weight     = -1;
            slope          = round(diff(surf_bins));
            bounds         = [];
            viterbi_weight = ones([1 Nx]);
            mu_size        = 31;
            mu             = log10(exp(-(-(mu_size-1)/2 : (mu_size-1)/2).^4/1));
            mu(mu<-30)     = -30;
            mu             = mu - mean(mu);
            sigma          = sum(abs(mu))/10*ones(1,mu_size);
            smooth_weight  = 1;
            smooth_var     = inf;
            repulsion      = 150000;
            ice_bin_thr    = 10;
            
            %%%% TO COMPILE
            if 0
                tmp = pwd;
                cd ~/scripts/cresis-toolbox/cresis-toolbox/+tomo/
                mex -largeArrayDims viterbi.cpp
                cd(tmp);
            end
            %%%%
            
            offset = gt(1,1) - x(1);
            if(gt(1,2) + offset > length(image_x))
                error('Index exceeds matrix dimensions, select smaller left-click-and-drag window.');
            end

            % Call viterbi.cpp
            tic
            labels = tomo.viterbi(double(detect_data), ...
                double(surf_bins), double(bottom_bin), ...
                double(gt), double(mask), ...
                double(mu), double(sigma), double(egt_weight), ...
                double(smooth_weight), double(smooth_var), double(slope), ...
                int64(bounds), double(viterbi_weight), ...
                double(repulsion), double(ice_bin_thr));
            toc
            keyboard
            % Correct layer offset
            y_new  = labels(round(x(1)+offset) : round(x(2)+offset));
            
            % Interpolate layer to match image y-axis
            y_new  = interp1(1:length(image_y), image_y, y_new);
            cmds(end+1).undo_cmd = 'insert';
            cmds(end).undo_args = {cur_layer, auto_idxs, ...
                param.layer.y{cur_layer}(auto_idxs), ...
                param.layer.type{cur_layer}(auto_idxs), ...
                param.layer.qual{cur_layer}(auto_idxs)};
            cmds(end).redo_cmd = 'insert';
            cmds(end).redo_args = {cur_layer, auto_idxs, ...
                y_new, ...
                2*ones(size(auto_idxs)), param.cur_quality*ones(size(auto_idxs))};
        end
    end
else
    %%% Do nothing for now
end

return

