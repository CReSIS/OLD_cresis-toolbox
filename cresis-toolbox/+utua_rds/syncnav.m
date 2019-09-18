function [ lat,lon,elev ] = syncnav( block,fname )
    % Michael Christoffersen
    % May 2017
    % syncnav Uses the time in a block file to extract better location info
    times = block.time;
    s = zeros(1,length(block.time));

    for i=1:length(block.time)
        t = times(i);
        t = split(t,'T');
        t = t(2);
        t = split(t,':');
        % 37 second to change from GPS to UTC
        s(i) = str2double(t(1))*60*60 + str2double(t(2))*60 + str2double(t(3))-37;
    end

    tdata = dlmread(fname,' ',1,0);
    lat = interp1(tdata(:,1),tdata(:,2),s);
    lon = interp1(tdata(:,1),tdata(:,3),s);
    elev = interp1(tdata(:,1),tdata(:,4),s);

end

