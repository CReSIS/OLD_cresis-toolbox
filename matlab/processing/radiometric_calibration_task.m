function [success] = radiometric_calibration_task(param)
% [success] = radiometric_calibration_task(param)
%
% Cluster task for radiometric_calibration. Responsible for develoment of
% radar equation variables and application of gains to data.
%
%
%
% Author: Anthony Hoch
% note: variable structure/definition created for qLook and SAR products
% is created to match structure in source script. (i.e. qlook data
% variables are copied from get_heights_task
%
% See also radiometric_calibration, get_heights_task

physical_constants;

% load data to apply correction
load(param.load.data_file_in)
% get_heights/qlook files include:      Surface, Latitude, Longitude, Elevation, GPS_time, Data, Time, Depth, param_records, param_get_heights
% CSARP/SAR files include:              param_csarp, param_combine_wf_chan, Latitude, Longitude, Elevation, GPS_time, Data, Time, Depth, array_param, param_records, Surface, Bottom


% bring in records for GPS/INS and surface information
%         records_fn = ct_filename_support(param,param.load.records_fn,'records');
%         records_ver = load(records_fn,'ver');
%         if isfield(records_ver,'ver')
%             records = load(records_fn);
%         else
%             load(records_fn, 'records');
%         end
% THIS IS NOT NEEDED, CONVERT EVERYTHING TO USE GPS/INS AND SURFACE FROM'fcs'


% if not a get_heights file, need to create time and surface variables
dataType = 'get_heights';
if ~exist('Time')
    dataType = 'CSARP';
    frame_str_loc = findstr(param.load.data_file_in,'wf');
    wf                  = str2num(param.load.data_file_in(frame_str_loc+[3:1:4]));
    Time                = wfs(wf).time;
    Surface             = fcs.surface;
    Data                = fk_data;
end



% create/pre-allocate processing variables
totalRange              = 0.*Time;          % range (not using 'range' to avoid function overloading)
totalDepth              = 0.*Time;          % range from surface, for ice volume attenuation calculation
permProfile             = 0.*Time;          % permitivity at each range sample ... at range, not at depth
lossProfile             = 0.*Time;          % attenuation/loss at each range sample, linear
lossCummulative         = 0.*Time;          % attenuation/loss accumulated from total path, two way, linear

permAir                 = 1;                % default assumptions
lossAir                 = 0;
permIce                 = 3.15;
lossIce                 = 0;                % currently unused, listed for completeness


% create standard variables for radar equations



% assign variables from spreadsheet
f                       = mean([param.radar.wfs(1).f0,param.radar.wfs(1).f1]);  % center frequency
Pt                      = 1;                % take value from radar_config.tx_weights( given recievers )
G                       = 1;                % fixed to 1, for isotropic radiator assumption??


                                            % this next section either returns a locc constant or a loss profile, need
                                            % to clean up ... dramatically


% assign/update ice perm and loss and attenuation vector
use_constant = false;
if ~isempty(param.radiometric.ice_loss_constant)            % if constant is available, then use it
    use_constant = true;
    if param.radiometric.ice_loss_constant >= 0             % convert dielectric to linear loss
        permIce = param.radiometric.ice_loss_constant;
        lossIce = 1-10^((-20/10)/(1000*2));                 % uses fixed -20 dB/km, current does not work, needs more parameters
    else                                                    % convert dB loss to linear loss
        % permIce stays ate default value
        lossIce = 1-10^((param.radiometric.ice_loss_constant/10)/(1000*2));
    end
else % profile create from outside function
    %[dielectric.depth,dielectric.er] = summitPerm(wfs(wf).fc);
    [filePath,fileName, fileExt] = fileparts(param.radiometric.ice_loss_fh);
    eval(['fh = @' fileName ';']);
    freq        = wfs(wf).fc;
    depth = linspace(0,c*Time(end)/2,10001).';  % this should cover all values of depth in this process
    [perm.depth,perm.er] = fh(freq,depth);
    perm.depthInt = (depth(1:end-1) + depth(2:end))/2;
    
    [TWtime,gain] = genPropProfileFromPerm(perm.depth,perm.er,freq);
end




for trace_idx = 1:size(Data,2)
    
    % determine lossCummulative  and totalRange for this trace
    if use_constant
        % use permConstant and lossConstant instead of vector profiles
        % this block is not loop sensitive (depth sensitive), but in else statement because other option is depth sensitive
        for range_idx = 1:length(Time)
            if Time(range_idx) < 0      % pre-transmit, use air values
                totalRange(range_idx)   = 0;
                lossProfile(range_idx)  = lossAir;
            elseif Time(range_idx) < Surface(trace_idx) % between transmit and surface
                totalRange(range_idx)   = c.*Time(range_idx)/(sqrt(real(permAir))*2);
                lossProfile(range_idx)  = lossAir;
            else % at or after surface
                totalRange(range_idx)   = c.*Surface(trace_idx)/(sqrt(real(permAir))*2) + c.*diff([Surface(trace_idx),Time(range_idx)])/(sqrt(real(permIce))*2);
                lossProfile(range_idx)  = lossIce;
            end
            lossCummulative(range_idx) = prod(ones(size(lossProfile(1:range_idx))) - lossProfile(1:range_idx).*diff([totalRange(1); totalRange(1:range_idx)]));
        end        
        
    else
        % this is slightly 'backward' because the permittivity profile is in depth and this script needs depth to determine permittivity
        % so, assuming continuity, will calculate one at a time and use previous depth to determine current permittiivity
        
                                        % should furst turn permittivity profile into time using velocity,
                                        % then project ... will have the same result as this does, but will
                                        % remove the interlaced time and epth domain processing
        
        
        for range_idx = 1:length(Time)
            if Time(range_idx) < 0      % pre-transmit, use air values
                totalRange(range_idx)   = 0;
                totalDepth = 0;
                permProfile(range_idx) = permAir;
                lossProfile(range_idx)  = lossAir;
                lossCummulative(range_idx) = prod(ones(size(lossProfile(1:range_idx))) - lossProfile(1:range_idx).*diff([totalRange(1); totalRange(1:range_idx)]));
            elseif Time(range_idx) <= Surface(trace_idx) % between transmit and surface
                totalRange(range_idx)   = c.*Time(range_idx)/(sqrt(real(permAir))*2);
                totalDepth = 0;
                permProfile(range_idx) = permAir;
                lossProfile(range_idx)  = lossAir;
                lossCummulative(range_idx) = prod(ones(size(lossProfile(1:range_idx))) - lossProfile(1:range_idx).*diff([totalRange(1); totalRange(1:range_idx)]));
            else % at or after surface
                % note: this does not work if range_idx < 2 or Time(range_idx) < 0
                totalRange(range_idx)   = totalRange(range_idx-1) + c.*diff(Time([range_idx-1 range_idx]))/(sqrt(real(permProfile(range_idx-1)))*2);
                totalDepth = totalRange(range_idx) - c.*Surface(trace_idx)/(sqrt(real(permAir))*2);
                permProfile(range_idx)  = interp1(perm.depthInt,perm.er,totalDepth,'spline','extrap');
                lossProfile(range_idx)  = 0; % not used past the surface in this implementation (may calculate for debug)
                lossCummulative(range_idx) = interp1(perm.depthInt,gain,totalDepth,'spline','extrap'); % without 'extrap' the value at zero depth returns NaN
            end
        end
    end
    
    
    
    % calculate total correction to trace PRIMARY EQUATION DEPENDENCE
    % BREAK THIS APART A LITTLE MORE SO THAT IT IS FAR EASIER TO READ
    if strcmpi(param.radiometric.type,'specular')
        totalCorrection = ( (64*pi^2*f^2*sqrt(real(permAir))^2) / (Pt * G^2 * c^2) ) .* totalRange.^2 ./ lossCummulative;
    elseif strcmpi(param.radiometric.type,'distributed')
        totalCorrection = ( (64*pi^3*f^2*sqrt(real(permAir))^2) / (Pt * G^2 * c^2) ) .* totalRange.^4 ./ lossCummulative;
    else
        totalCorrection = ( (64*pi^2*f^2*sqrt(real(permAir))^2) / (Pt * G^2 * c^2) ) .* totalRange.^2 ./ lossCummulative;
        totalCorrection = ones(size(totalCorrection));
        % this is just the specular correction to get the correct size, then converted to ones due to the failure
    end
    % note: permAir use here implies that the antenna is recieving in air, would not be accurate for ground missions
    
    
    % apply total correction to trace
    Data(:,trace_idx) = totalCorrection .* Data(:,trace_idx);
    % note: data is already convereted to power and incoherently integrated
    % !!! this script needs to account for integration, but fir_dec does not alter signal power (should not change result)
end







% save information
% THIS PART WILL BE REMOVED WHEN SCRIPT BLOCK IS MOVED, THIS
% CREATES SEPERATE COPIES OF THE RADIOMETRIC CALIBRATION DATA, LEAVING THE STANDARD PRODUCTS ALONE

[file_path file_name file_ext] = fileparts(param.load.data_file_out);
if ~exist(file_path,'dir')
    mkdir(file_path);
end

if      strcmpi(dataType,'get_heights')
    save(param.load.data_file_out, 'Data', 'Time', 'Surface', 'GPS_time', 'Latitude', 'Longitude', 'Elevation');
elseif  strcmpi(dataType,'CSARP')
    fk_data = Data;
    param_radiometric = param.radiometric;
    save(param.load.data_file_out, 'wfs', 'fcs', 'fk_data', 'param_records', 'param_csarp', 'param_radiometric');
end
    




success = true;

return;
