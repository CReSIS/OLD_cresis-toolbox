% Michael Christoffersen
% May 2017
% Wrapper for syncnav

datadir = '/media/mchristo/blind-swamp-eel/paden/radar';
gpsdir = '/media/mchristo/blind-swamp-eel/paden/gps';
datafiles = dir(datadir);
nfiles = length(datafiles);

for i=1:nfiles
    dfname = datafiles(i).name;
    disp(dfname)
    if(~endsWith(dfname,'.mat'))
        continue
    end
    dfname = [datadir,'/',dfname];
    load(dfname)

    ref = split(block.time(1),'T');
    ref = ref(1);
    ref = split(ref,'-');

    % Relies on specific directory structure here
    year = ref(1);
    month = ref(2);
    day = ref(3);
    gfname = [gpsdir,'/','18',char(month),char(day),'.pos'];

    [block.lat,block.lon,block.elev_air] = syncnav(block,gfname);
    block.gps_corrected = 1;
    save(dfname,'block');
end

clear block day dfname gfname i month year ref datadir gpsdir datafiles nfiles
