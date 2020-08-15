function collate_snapshots(param,param_override)
% collate_snapshots(param,param_override)
%
% Function that takes snapshots and pulls out calibration pixels and stores
% with corresponding surface ARRIVAL ANGLES (not surface intersection
% angles, i.e. roll-corrected)
%
% TO DO:  ADD FUNCTIONALITY TO STORE MASK IN A MEANINGFULL WAY FOR BROWSING
% IN THE PICKER.
%
% param: struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override: parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_FUNCTION.m for how to run this function directly.
%  This function may be called from the run_master.m script using the
%  param spreadsheet and the cmd.generic column.
%
% Authors: Theresa Moore
%
% See also: collate_snapshots.m

%% General Setup
% =====================================================================

param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =====================================================================
if ~isfield(param,mfilename) || isempty(param.(mfilename))
  param.(mfilename) = [];
end

if ~isfield(param.collate_snapshots,'scr_threshold') || isempty(param.collate_snapshots.scr_threshold)
  param.collate_snapshots.scr_threshold = 15;
end

if ~isfield(param.collate_snapshots,'out_path') || isempty(param.collate_snapshots.out_path)
  param.collate_snapshots.out_path = 'collate_snapshots';
end

if ~isfield(param.collate_snapshots,'in_path') || isempty(param.collate_snapshots.in_path)
  param.collate_snapshots.in_path = 'snapshot';
end

if ~isfield(param.collate_snapshots,'layer_params') || isempty(param.collate_snapshots.layer_params)
  param.collate_snapshots.layer_params.name = 'surface';
  param.collate_snapshots.layer_params.source = 'layerdata';
  param.collate_snapshots.layer_params.layer_path = 'layer';
end

if ~isfield(param.collate_snapshots.layer_params,'layer_path') || isempty(param.collate_snapshots.layer_params.layer_path)
  param.collate_snapshots.layer_params.layer_path = 'layer';
end

if ~isfield(param.collate_snapshots,'remove_multiple') || isempty(param.collate_snapshots.remove_multiple)
  param.collate_snapshots.remove_multiple = true;
  param.collate_snapshots.multiple_guard_bins = 2;
end

if ~isfield(param.collate_snapshots,'multiple_guard_bins') || isempty(param.collate_snapshots.multiple_guard_bins)
  param.collate_snapshots.multiple_guard_bins = 2;
end

if ~isfield(param.collate_snapshots,'img_list')
  param.collate_snapshots.img_list = [];
end

if ~isfield(param.collate_snapshots,'plot_enable') || isempty(param.collate_snapshots.plot_enable)
  param.collate_snapshots.dtheta = 1;
end

if ~isfield(param.collate_snapshots,'debug') || isempty(param.collate_snapshots.debug)
  param.collate_snapshots.debug = 0;
end

if ~isfield(param.collate_snapshots,'same_doa_guard') || isempty(param.collate_snapshots.same_doa_guard)
  param.collate_snapshots.same_doa_guard = 5; % degrees
end

% 2014-05-06 16:13:17.84
% 2014-05-06 16:13:36.82

%%
% Load frames file
frames = frames_load(param);

param.cmd.frms = frames_param_cmd_frms(param,frames);

for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  %   fprintf('Collate snapshots %s_%03d \n', param.day_seg,frm)
  fprintf('\n\n');
  fprintf('collate snapshots:  %s_%03d (%s)\n', param.day_seg,frm, datestr(now));
  fprintf('---------------------------------------------------------------------\n');
  
  % Building input and output directories for the frame
  snp_dir = ct_filename_out(param,param.collate_snapshots.in_path);
  sv_offset_dir = ct_filename_out(param,param.collate_snapshots.out_path);
  
  % Check for existence of layerdata
  layer_fn = fullfile(ct_filename_out(param,param.collate_snapshots.layer_params.layer_path,''),sprintf('Data_%s_%03d.mat', param.day_seg,frm));
  
  % If surface layer doesn't exist, disable the flag to remove the multiple
  if ~(exist(layer_fn,'file') == 2)
    param.collate_snapshots.remove_multiple = false;
    fprintf('Layer file does not exist: %s,\n', layer_fn);
    fprintf('Surface multiple not removed\n');
    layer = [];
  else
    % Load in surface
    surf_param = param.collate_snapshots.layer_params;
    fprintf('Loading surface:  %s (%s)\n', layer_fn, datestr(now));
    layer = opsLoadLayers(param,surf_param);
  end
  
  
  % Loop over images in the list
  for img_idx = 1:length(param.collate_snapshots.img_list)
    % Initialize variables
    Tomo = [];
    Time = [];
    Elevation = [];
    GPS_time = [];
    Heading = [];
    Latitude = [];
    Longitude = [];
    Roll = [];
    Heading =[];
    Latitude =[];
    Longitude =[];
    Pitch =[];
    Roll=[];
    param_array =[];
    param_records=[];
    param_sar = [];
    file_version=[];
    
    
    % Building the input filename
    img= param.collate_snapshots.img_list(img_idx);
    file = sprintf('Data_img_%02d_%s_%03d.mat',img,param.day_seg,frm);
    fn = fullfile(snp_dir,file);
    
    % Ensure file exists
    if exist(fn,'file') == 2
      fprintf('Loading snapshot file:  %s (%s)\n', fn, datestr(now));
      load(fn,'Tomo','Time','Elevation','GPS_time','Heading','Latitude','Longitude','Pitch','Roll','param_array','param_records','param_sar','file_version');
      [Nt, Nc, Nx]  = size(Tomo.img);
      
      % Reformat data
      Tomo.img = permute(Tomo.img,[2 1 3]);
      Tomo.power = permute(Tomo.power,[2 1 3]);
      Tomo.surf_theta = permute(Tomo.surf_theta,[2 1 3]);
      Tomo.surf_ice_mask = permute(Tomo.surf_ice_mask,[2 1 3]);
      Tomo.two_src_mask = squeeze(sum(~isnan(Tomo.surf_theta),1))==2;
      Tomo.surf_ice_mask = squeeze(any(Tomo.surf_ice_mask,1));
      Tomo.surf_ice_mask = cumsum(Tomo.surf_ice_mask,1) >= 1;
      Tomo.surf_multiple_mask = true(size(Tomo.surf_ice_mask));
      
      % Interpolate surface onto snapshot slow time axis
      if ~isempty(layer)
        surface = interp_finite(interp1(layer.gps_time,layer.twtt,GPS_time));
        for rline_idx = 1:Nx
          surf_bin = round(interp_finite(interp1(Time,1:Nt,surface(rline_idx))));
          surf_bins = surf_bin + [-param.collate_snapshots.multiple_guard_bins:param.collate_snapshots.multiple_guard_bins];
          if min(surf_bins) > 1 & max(surf_bins) < Nt
            Tomo.surf_multiple_mask(surf_bins,rline_idx) = false;
          end
        end
      end
      
      % Assign nans to -inf
      Tomo.power(isnan(Tomo.power)) = -inf;
      
      % Sort powers in descending order
      [Tomo.power,sort_idxs] = sort(Tomo.power,1,'descend');
      Tomo.power = 10*log10(abs(Tomo.power));
      
      
      % Convert surface angle to DOA
      Tomo.surf_doa = bsxfun(@minus,Tomo.surf_theta,permute(Roll*180/pi,[1 3 2]));
      
      % Maximum number of sources
      Ndoa_max = size(Tomo.surf_doa,1);
      doa_sort_idxs = bsxfun(@plus,bsxfun(@plus,sort_idxs,Ndoa_max*(0:Nt-1)),permute(Ndoa_max*Nt*(0:Nx-1),[1 3 2]));
      
      if Ndoa_max > 1
        % Threshold
        Tomo.power_mask = squeeze(Tomo.power(1,:,:) > Tomo.power(2,:,:) + param.collate_snapshots.scr_threshold);
      else
        Tomo.power_mask = true(size(squeeze(Tomo.power)));
      end
      
      % Number of sources per pixel
      %       nsrc_per_pixel = squeeze(sum(~isnan(Tomo.surf_theta),1));
      
      %       Tomo.surf_theta = Tomo.surf_theta(theta_sort_idxs);
      Tomo.surf_theta = Tomo.surf_theta(doa_sort_idxs);
      
      Tomo.surf_doa   = Tomo.surf_doa(doa_sort_idxs);
      Tomo.same_doa_mask = squeeze(abs(Tomo.surf_doa(1,:,:)-Tomo.surf_doa(2,:,:)) < param.collate_snapshots.same_doa_guard);
      
      Tomo.mask = Tomo.two_src_mask & ~Tomo.surf_ice_mask & (Tomo.power_mask | Tomo.same_doa_mask) & Tomo.surf_multiple_mask;
      
      
      %% Create simplified output
      
      % Pick out the snapshots of the good pixels
      snapshots.sv_list = Tomo.img(:,Tomo.mask);
      
      % Pick out the good surface intersections and doas
      snapshots.surf_doa = Tomo.surf_doa(1,Tomo.mask);
      snapshots.surf_theta = Tomo.surf_theta(1,Tomo.mask);
      
      % Pick out corresponding power
      snapshots.power = Tomo.power(1,Tomo.mask);
      
      % Pick out gps time and twtt of the culled pixels
      twtt_matrix     = single(repmat(Time,1,Nx));
      gps_time_matrix = single(repmat(GPS_time,Nt,1));
      elev_matrix     = single(repmat(Elevation,Nt,1));
      pitch_matrix    = single(repmat(Pitch,Nt,1));
      heading_matrix  = single(repmat(Heading,Nt,1));
      roll_matrix     = single(repmat(Roll,Nt,1));
      lat_matrix      = single(repmat(Latitude,Nt,1));
      lon_matrix      = single(repmat(Longitude,Nt,1));
      
      snapshots.twtt = twtt_matrix(Tomo.mask);
      snapshots.gps_time = gps_time_matrix(Tomo.mask);
      snapshots.elev = elev_matrix(Tomo.mask);
      snapshots.lat = lat_matrix(Tomo.mask);
      snapshots.lon = lon_matrix(Tomo.mask);
      snapshots.pitch = pitch_matrix(Tomo.mask);
      snapshots.heading = heading_matrix(Tomo.mask);
      snapshots.roll = roll_matrix(Tomo.mask);
      
      snapshots.twtt = snapshots.twtt(:).';
      snapshots.gps_time = snapshots.gps_time(:).';
      snapshots.elev = snapshots.elev(:).';
      snapshots.lat = snapshots.lat(:).';
      snapshots.lon = snapshots.lon(:).';
      snapshots.pitch = snapshots.pitch(:).';
      snapshots.heading = snapshots.heading(:).';
      snapshots.roll = snapshots.roll(:).';
      
      snapshots.param_collate_snapshots = param.collate_snapshots;
      snapshots.param_array = param_array;
      snapshots.param_sar = param_sar;
      snapshots.param_records = param.records;
      snapshots.file_version = file_version;
      snapshots.file_type = 'snapshots';
      
      %% Save the result
      % =====================================================================
      out_fn_dir = ct_filename_out(param,param.collate_snapshots.out_path,'');
      if ~isdir(out_fn_dir)
        mkdir(out_fn_dir)
      end
      out_fn = fullfile(out_fn_dir,sprintf('snapshots_simp_img_%02d_%s_%03d.mat', img, param.day_seg, frm));
      fprintf('Saving %s (%s)\n', out_fn, datestr(now));
      ct_save(out_fn,'-struct','snapshots');
    end
  end
end
% Phase centers
%       param_array.gps_source = param_records.gps_source;
%       for wf_adc = 1:size(param_array.array.imgs{img}{1},1)
%         wf = param_array.array.imgs{img}{1}(wf_adc,1);
%         adc = param_array.array.imgs{img}{1}(wf_adc,2);
%         rxchannel = param_array.radar.wfs(wf).rx_paths(adc);
%         phase_centers(:,wf_adc) = lever_arm(param_array, param_array.radar.wfs(wf).tx_weights, rxchannel);
%       end
%
%       % Create angle bins
%       bins = [min(param_array.array_proc.surface_theta):param.collate_snapshots.dtheta:max(param_array.array_proc.surface_theta)];
%
%       % Compute nominal steering vector forthe bin
%       Nsv = {'theta'};
%       Nsv{2} = bins/180*pi;
%       [~,sv] = array_proc_sv(Nsv, param_array.radar.wfs(1).fc, phase_centers(2,:).', phase_centers(3,:).');
%
%       % Preallocate memory to story mean, median, variance and number of
%       % samples for each bin.
%       sv_mean{img_idx,frm_idx} = zeros(size(sv_list,1),length(bins));
%       sv_var{img_idx,frm_idx} = zeros(size(sv_list,1),length(bins));
%       sv_median{img_idx,frm_idx} = zeros(size(sv_list,1),length(bins));
%       num_sam{img_idx,frm_idx} = zeros(1,length(bins));
%
%       % Loop over bins
%       for bin_idx = 1:length(bins)
%         bin_angle = bins(bin_idx);
%
%         bin_mask = abs(theta-bins(bin_idx)) < param.collate_snapshots.dtheta/2;
%
%         sv_list_bin{bin_idx} = sv_list(:,bin_mask);
%
%         for sv_idx = 1:size
%           % Apply correction
%
%         end
%
%
%         % Average relative amplitude and phase offset for the bin
%         rel_amplitude = sqrt(mean(abs(sv_list(:,bin_mask)).^2,2));
%         rel_angle = angle(mean(sv_list(:,bin_mask),2));
%
%         sv_mean{img_idx,frm_idx}(:,bin_idx) = rel_amplitude .* exp(1i*rel_angle);
%         num_sam{img_idx,frm_idx}(bin_idx) = sum(bin_mask);
%
%         % Debug plotting
%         if param.collate_snapshots.debug && param.collate_snapshots.plot_enable
%           if param.collate_snapshots.plot_enable && 0
%
%             figure(1);clf
%             hold on
%             grid on
%             grid minor
%
%             color_vec = hsv(Nc);
%             for id = 1:Nc
%               plot(real(sv_list_bin{bin_idx}(id,:)),imag(sv_list_bin{bin_idx}(id,:)),'Marker','.','Color',hsv(id));
%             end
%
%             figure(5); clf;
%             figure(4); clf;
%
%             % Debug print outs and plots
%             fprintf('%.2f deg %d\n', mean(theta(bin_mask)), sum(bin_mask));
%             fprintf('  %.1f - %.1f dB %.1f\n', min(power(bin_mask)), max(power(bin_mask)), 10*log10(mean(10.^(power(bin_mask)/10))));
%
%             ff=unwrap(angle( SVmean{img,frm}(:,bin_idx) .* sv(:,bin_idx) ));
%             %       ff=unwrap(angle( SVmean(:,bin_idx) ));
%             %       ff=unwrap(angle( sv(:,bin_idx) ));
%             gg=SVmean{img,frm}(:,bin_idx);
%
%             figure(5);
%             plot((ff-ff(ref_channel))*180/pi);
%             grid on;
%             hold on;
%             title(sprintf('%.2f deg\n', bin_angle));
%             xlabel('Rx path');
%             ylabel('Angle of steering vector (deg)');
%
%             figure(4);
%             plot(lp(gg) - lp(gg(ref_channel)));
%             grid on;
%             hold on;
%             title(sprintf('%.2f deg\n', bin_angle));
%             xlabel('Rx path');
%             ylabel('Power steering vector (dB)');
%             %               pause
%           end
%
%
%           figure(img*100+frm); clf;
%           subplot(3,1,1:2);
%           aa = gca;
%           imagesc(bins,[],angle(sv_mean{img_idx,frm_idx})*180/pi);
%           colorbar
%           subplot(3,1,3);
%           aa(2) = gca;
%           plot(bins, num_sam{img,frm});
%           linkaxes(aa,'x');
%         end
%
%       end
%
%       % Save results
%       wf = param_array.array_proc.imgs{img}(1,1);
%       adc_list = param_array.array_proc.imgs{img}(:,2);
%       rx_list
%       file = sprintf('sv_offsets_img_%02d_%s_%03d.mat',param_array.array_proc.imgs{img}(1,1),param.day_seg,frm);
%       out_dir = ct_filename_out(param,param.collate_snapshots.out_path)
%       if ~isdir(out_dir)
%         mkdir(out_dir);
%       end
%
%       out_fn = fullfile(param.collate_snapshots.out_path,file);
%       fprintf('  Saving cal data %s (%s)\n', out_fn, datestr(now));
%       %             ct_save('-v7.3',out_fn,'SVmatrix','Nsrcs_frm','doa_max_source','power_max_source','ice_mask_max','doa_source_two','power_source_two','ice_mask_source_two','Latitude','Longitude','Elevation','GPS_time', ...
%       %               'Time','param_array','param_records', ...
%       %               'param_sar', 'Roll', 'Pitch', 'Heading','file_version');
%
% %     else
% %       fprintf('File does not exist: %s (%s)\n\n', fn, datestr(now));
% %
% %     end
% %
%       ref_channel = param.collate_snapshots.ref_chan;
%
%       % Relative offset from the reference channel
%       sv_list = bsxfun(@times, sv_list, exp(-1i*angle(sv_list(ref_channel,:))) );

%
%     keyboard
%
%   end
%
%% Save the result
%     % =====================================================================
%     out_fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.out_path, ''));
%     out_fn = fullfile(out_fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
%     fprintf('Saving %s (%s)\n', out_fn, datestr(now));
%     ct_save(out_fn,'-struct','noise_simp');



