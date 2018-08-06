function [hdr,data] = data_resample(hdr,data,pq)
% [hdr,data] = data_resample(hdr,data,pq)
%
% https://ops.cresis.ku.edu/wiki/index.php/Data_load#data_resample.m
%
% Author: John Paden

% Fast time resampling
if ~isequal(pq(1,:),[1 1])
  for img = 1:length(data)
    % Matlab R2015a resample has bug when working with 3+ dimensions, we must do
    % reshape to work around bug.
    s = size(data{img});
    data{img} = single(resample(double(data{img}), ...
      pq(1,1), pq(1,2)));
    s(1) = size(data{img},1);
    data{img} = reshape(data{img},s);
    
    % Resample time axis
    dt = hdr.time{img}(2) - hdr.time{img}(1);
    Nt = size(data{img},1);
    hdr.time{img} = hdr.time{img}(1) + pq(1,2)/pq(1,1)*dt * (0:Nt-1).';
    
    % Resample frequency axis
    df = 1/(Nt*dt);
    hdr.freq{img} = hdr.freq{img}(1) + df * ifftshift(floor(-Nt/2) : floor((Nt-1)/2)).';
  end
end

% Along track resampling
if ~isequal(pq(2,:),[1 1])
  for img = 1:length(data)
    % Matlab R2015a resample has bug when working with 3+ dimensions, we must do
    % reshape to work around bug.
    s = size(data{img});
    s([1 2]) = s([2 1]); % because of permute
    data{img} = single(resample(double(permute(data{img},[2 1 3:ndims(data{img})])), ...
      pq(1,1), pq(1,2)));
    s(1) = size(data{img},1);
    data{img} = reshape(data{img},s);
    data{img} = permute(data{img},[2 1 3:ndims(data{img})]);
    
    % Resample along-track position vectors
    hdr.records{img}.gps_time = resample(hdr.records{img}.gps_time, pq(2,1), pq(2,2));
    hdr.records{img}.lat = resample(hdr.records{img}.lat, pq(2,1), pq(2,2));
    hdr.records{img}.lon = resample(hdr.records{img}.lon, pq(2,1), pq(2,2));
    hdr.records{img}.elev = resample(hdr.records{img}.elev, pq(2,1), pq(2,2));
    hdr.records{img}.roll = resample(hdr.records{img}.roll, pq(2,1), pq(2,2));
    hdr.records{img}.pitch = resample(hdr.records{img}.pitch, pq(2,1), pq(2,2));
    hdr.records{img}.heading = resample(hdr.records{img}.heading, pq(2,1), pq(2,2));
    hdr.records{img}.surface = resample(hdr.records{img}.surface, pq(2,1), pq(2,2));
    
  end
end
