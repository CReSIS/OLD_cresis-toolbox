function [ch0_chirp2d,ch0_bark2d,ch1_chirp2d,ch1_bark2d] = slice_tdms(tdms)
    % Michael Christoffersen
    % May 2017
    % Function to slice 1D records from the radar in to 2D arrays
    
    reclen = tdms.Props.record_len;
    ntrace = length(tdms.radar.ch0.data)/reclen;
    
    if(isfield(tdms.Props,'bark'))
        bark = tdms.Props.bark;
    else
        bark = tdms.Props.pulse;
    end
    
    if(bark)
        if(isfield(tdms.Props,'bark_len'))
            barklen = ceil((tdms.Props.bark_len+tdms.Props.bark_delay)/tdms.Props.dt); 
        else
            barklen = ceil((tdms.Props.pulse_len+tdms.Props.pulse_delay)/tdms.Props.dt); 
        end
        ch0_chirp2d = zeros(reclen-barklen, ntrace);
        ch1_chirp2d = zeros(reclen-barklen, ntrace);
        ch0_bark2d = zeros(barklen,ntrace);
        ch1_bark2d = zeros(barklen,ntrace);
        for i=1:ntrace
            ch0_bark2d(:,i) = tdms.radar.ch0.data(1+reclen*(i-1):1+reclen*(i-1)+barklen-1);
            ch1_bark2d(:,i) = tdms.radar.ch1.data(1+reclen*(i-1):1+reclen*(i-1)+barklen-1);
            ch0_chirp2d(:,i) = tdms.radar.ch0.data(1+reclen*(i-1)+barklen:reclen*(i));  
            ch1_chirp2d(:,i) = tdms.radar.ch1.data(1+reclen*(i-1)+barklen:reclen*(i));   
        end
    else
        ch0_chirp2d = zeros(reclen, ntrace);
        ch1_chirp2d = zeros(reclen, ntrace);
        ch0_bark2d = 0;
        ch1_bark2d = 0;
        for i=1:ntrace
            ch0_chirp2d(:,i) = tdms.radar.ch0.data(1+reclen*(i-1):reclen*(i));
            ch1_chirp2d(:,i) = tdms.radar.ch1.data(1+reclen*(i-1):reclen*(i));
    end
end