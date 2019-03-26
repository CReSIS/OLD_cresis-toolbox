function [block] = tdms2block(tdms)
    % Michael Christofferen
    % May 2017
    % Converts TDMS files from radar to block files for picking in pickgui
    
    %% 1D data record -> 2D
    [ch0_chirp,ch0_bark,ch1_chirp,ch1_bark] = utua_rds.slice_tdms(tdms);
    block.ch0 = ch0_chirp; %*tdms.radar.ch0.Props.NI_Scale_0__Linear_Slope + tdms.radar.ch0.Props.NI_Scale_0__Linear_Y_Intercept;
    block.clutter = zeros(size(block.ch0));
    block.ch1 = ch1_chirp; %*tdms.radar.ch1.Props.NI_Scale_0__Linear_Slope + tdms.radar.ch1.Props.NI_Scale_0__Linear_Y_Intercept;
    %block.ch0_bark = ch0_bark;%*tdms.radar.ch0.Props.NI_Scale_0__Linear_Slope + tdms.radar.ch0.Props.NI_Scale_0__Linear_Y_Intercept;;
    %block.ch1_bark = ch1_bark;%*tdms.radar.ch1.Props.NI_Scale_0__Linear_Slope + tdms.radar.ch1.Props.NI_Scale_0__Linear_Y_Intercept;;
    
    %% Metadata
    block.dt = tdms.Props.dt;
    block.name = tdms.Props.name;
    
    %if(isfield(tdms.Props,'pulse_len'))
    %    block.bark.len = tdms.Props.pulse_len;
    %    block.bark.delay = tdms.Props.pulse_delay;
    %else
    %    block.bark.len = tdms.Props.bark_len;
    %    block.bark.delay = tdms.Props.bark_delay;
    %end
    %if(isfield(tdms.Props,'pulse'))
    %    block.bark = tdms.Props.pulse;
    %else
    %    block.bark = tdms.Props.bark;
    %end
    
    % Chirp info
    block.chirp.len = tdms.Props.chirp_len;
    block.chirp.cf = tdms.Props.chirp_cf;
    block.chirp.bw = tdms.Props.chirp_bw;
    block.chirp.amp = tdms.Props.chirp_amp;
    block.prf = tdms.Props.prf;
    
    block.stack = tdms.Props.stacking;
    block.num_sample = tdms.Props.record_len;
    block.num_trace = length(tdms.radar.ch0.data)/block.num_sample;
    
    % Preamp
    block.ritec.gain = tdms.Props.ritec_gain;
    block.ritec.lo_filter = tdms.Props.ritec_lo_filter;
    block.ritec.hi_filter = tdms.Props.ritec_hi_filter;
    block.ritec.impedance = tdms.Props.ritec_impedance;
    block.comment = tdms.Props.comment;
    
    % Location, the first if clause is for an odd bug in Aug 2017 data
    if(length(tdms.meta.lat.data) > 2*block.num_trace)
        block.lat = tdms.meta.lat.data(1:2:end-1);
        block.lon = tdms.meta.lon.data(1:2:end-1);
        block.elev_air = tdms.meta.elev.data(1:2:end-1);
        tdms.meta.time.data = tdms.meta.time.data(1:2:end-1);
    else
        % Otherwise just strip the one extra location/time off the end
        block.lat = tdms.meta.lat.data(1:end-1);
        block.lon = tdms.meta.lon.data(1:end-1);
        block.elev_air = tdms.meta.elev.data(1:end-1);
        tdms.meta.time.data = tdms.meta.time.data(1:end-1);
    end
    
    % Time
    for i=1:length(tdms.meta.time.data)-1
        if(tdms.meta.time.data(i) == 0)
            tdms.meta.time.data(i) = (tdms.meta.time.data(i-1)+tdms.meta.time.data(i+1))/2;
        end
    end
    tdms.meta.time.data(end) = tdms.meta.time.data(end-1) + (tdms.meta.time.data(end-1)-tdms.meta.time.data(end-2));
    block.time = zeros(1,length(tdms.meta.time.data));
    st = datenum(tdms.Props.start_time,'yyyy-mm-ddTHH:MM:SS.FFF');
    st = datetime(st,'ConvertFrom','datenum');
    
    %Display start time and number of traces
    disp(st)
    disp(length(tdms.meta.time.data));
    
    for i=1:length(tdms.meta.time.data)
       t = st + ((tdms.meta.time.data(i) - tdms.meta.time.data(1))/86400);
       %block.time(i) = datetime(t,'Format','yyyy-MM-dd''T''HH:mm:ss.SSSSSS');
       block.time(i) = datenum(t);
    end
    
    if(isfield(tdms.radar.ch0.Props,'range'))
        block.ch0_meta.range = tdms.radar.ch0.Props.range;
        block.ch1_meta.range = tdms.radar.ch1.Props.range;
    end
    
    if(isfield(tdms.radar.ch0.Props,'NI_Scale_0__Linear_Slope'))
        block.ch0_meta.scale_slope = tdms.radar.ch0.Props.NI_Scale_0__Linear_Slope;
        block.ch0_meta.scale_intercept = tdms.radar.ch0.Props.NI_Scale_0__Linear_Y_Intercept;
        block.ch1_meta.scale_slope = tdms.radar.ch1.Props.NI_Scale_0__Linear_Slope;
        block.ch1_meta.scale_intercept = tdms.radar.ch1.Props.NI_Scale_0__Linear_Y_Intercept;
    end
    
    block.ch0_meta.atten = tdms.radar.ch0.Props.atten;
    block.ch1_meta.atten = tdms.radar.ch1.Props.atten;
    block.amp = utua_rds.pcsynthetic(block);
    
    % Chop off non unique time values, again weird Aug 2017 bug
    if(length(block.time) ~= length(unique(block.time)))
        [block.time,ia,ic] = unique(block.time);
        %block.lat = block.lat(ia);
        %block.lon = block.lon(ia);
        %block.elev_air = block.elev_air(ia);
        %block.x = block.x(ia);
        %block.y = block.y(ia);
        block.amp = block.amp(ia);
        block.ch0 = block.ch0(:,ia);
        block.ch1 = block.ch1(:,ia);
        %block.ch0_bark = block.ch0_bark(:,ia);
        %block.ch1_bark = block.ch1_bark(:,ia);
        block.num_trace = length(ia);
    end
    
    block.gps_corrected = 0;
    
    % Convert lat/lon to utm
    block.x = block.lat;
    block.y = block.lat;
    [block.x,block.y] = utua_rds.wgs2utm(block.lat,block.lon,6,'N');
    block.x = block.x/1000; % km for pickgui
    block.y = block.y/1000;
    
    block.dist = block.x;
    block.dist(1) = 0;
    for i=2:length(block.dist)
        block.dist(i) = block.dist(i-1) + sqrt((block.x(i)-block.x(i-1))^2+(block.y(i)-block.y(i-1))^2);
    end
    block.dist = block.dist;
    block.dist_lin = [0:length(block.dist)-1]/(length(block.dist)-1);
    block.dist_lin = block.dist_lin*block.dist(end);
    
    %block.lat_init = tdms.meta.lat.data;
    %block.lon_init = tdms.meta.lon.data;
    %block.elev_air_init = tdms.meta.elev.data;
    
    % Resample to 20ns x 2048 samples for pickgui, resampling to double for resample func
    factor = ((20e-9)/block.dt); %calculate the new sampling ratio
    block.ch0 = double(block.ch0);
    block.ch1 = double(block.ch1);
    block.amp = double(block.amp);
    %block.clutter = double(block.clutter);

    block.ch0 = resample(block.ch0,1,factor);
    block.ch1 = resample(block.ch1,1,factor);
    block.amp = resample(block.amp,1,factor);
    %block.clutter = resample(block.clutter,1,factor);
    
    % Back to single precision to save space
    block.dt = 20e-9;
    block.ch0 = single(block.ch0);
    block.ch1 = single(block.ch1);
    block.amp = single(block.amp);
    block.clutter = single(block.clutter);
   
    if size(block.ch0,1) >= 2048 %if new ch0 is more than 2048 traces long, truncate. Otherwise, pad out to 2048 with NaN
        block.ch0 = block.ch0(1:2048,:);
    else
        samples0 = size(block.ch0,1);
        block.ch0(samples0:2048,:) = NaN;
    end
   
   if size(block.ch1,1) >= 2048 %if new ch1 is more than 2048 traces long, truncate. Otherwise, pad out to 2048 with NaN
       block.ch1 = block.ch1(1:2048,:);
   else
       samples0 = size(block.ch1,1);
       block.ch1(samples0:2048,:) = NaN;
   end
   
   if size(block.amp,1) >= 2048 %if new amp is more than 2048 traces long, truncate. Otherwise, pad out to 2048 with NaN
       block.amp = block.amp(1:2048,:);
   else
       samples0 = size(block.amp,1);
       block.amp(samples0:2048,:) = NaN;
   end
   
   if size(block.clutter,1) >= 2048 %if new clutter is more than 2048 traces long, truncate. Otherwise, pad out to 2048 with NaN
       block.clutter = block.clutter(1:2048,:);
   else
       samples0 = size(block.clutter,1);
       block.clutter(samples0:2048,:) = NaN;
   end
   
   block.num_sample = 2048;
   block.twtt = [0:block.dt:(size(block.amp,1)-1)*block.dt]';
   block.twtt_surf = zeros(1,size(block.amp,2));
   block.ind_overlap = [NaN, NaN];
end
