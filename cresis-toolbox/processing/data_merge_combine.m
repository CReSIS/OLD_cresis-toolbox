function [hdr,data] = data_merge_combine(param,hdr,data)

%% Combine wf-adc pairs into a single wf-adc pair
% =========================================================================

% Apply motion compensation






  %% Roll compensation
%   if param.qlook.roll_correction
%     % Apply roll-only motion compensation
%     for wf_adc = 1:size(param.load.imgs{img},1)
%       wf = abs(param.load.imgs{img}(wf_adc,1));
%       adc = abs(param.load.imgs{img}(wf_adc,2));
%       rx = wfs(wf).rx_paths(adc);
%       for rline = 1:size(g_data,2)
%         drange = radar_lever_arm(2,wf_adc) * -tan(out_records.roll(rline));
%         dphase = drange * 2 * 2*pi/lambda_fc;
%         g_data(:,rline,wf_adc) = g_data(:,rline,wf_adc) * exp(1j*dphase);
%       end
%     end
%     g_data = mean(g_data,3);
%   end
%   
%   %% Elevation compensation
%   if param.qlook.elev_correction && ~simple_firdec
%     % Remove elevation variations (just a phase shift, not a time shift)
%     %  - With simple_firdec there is no point in elevation compensation
%     %    because the along-track averaging has already been done
%     drange = out_records.elev-mean(out_records.elev);
%     dphase = drange * 2 * 2*pi/lambda_fc;
%     for rline = 1:size(g_data,2)
%       g_data(:,rline) = g_data(:,rline) .* exp(1j*dphase(rline));
%     end
%   end




% Create new trajectory for combined phase centers

% Combined
%data{state.img(img_comb_idx,accum_idx)}(Nt,Nx,COMBINE_TO_CREATE_NEW)

% Combine chain is a cell which contains lists for which wf-adc idx to
% include in the current cell index, also includes a set of Tsys,
% chan_equal style coefficients to apply to each wf-adc idx before
% combining.

%% Combine images into a single image
% =========================================================================
%data{state.img(COMBINE_TO_CREATE_SINGLE,accum_idx)}(Nt,Nx,Nc)



end