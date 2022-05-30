% Melt_CSV_shortening: Produces shortened resampled scripts with data
% points resampled every other points, every 4th and every 8th point.
% Script outputs multiple spreadsheets that are then ingested into ARCMAP
% for plotting CReSIS radar melt rates.

cd 'C:\Users\c262b531\Documents\Thesis related documents\Dynamic Thinning Correction\Finished Melt CSV\P2\'
my_files = pwd;
% For loop for loading
melt_CSV = dir(fullfile(my_files,'P2M*.csv'));
for i = 1:numel(melt_CSV)
    import_data = fullfile(my_files, melt_CSV(i).name);
    melt_data(i).data = readmatrix(import_data);
    melt_data(i).AT = melt_data(i).data(:,1).';
    melt_data(i).Lons = melt_data(i).data(:,2).';
    melt_data(i).Lats = melt_data(i).data(:,3).';
    melt_data(i).X = melt_data(i).data(:,4).';
    melt_data(i).Y = melt_data(i).data(:,5).';
    melt_data(i).total_melt = melt_data(i).data(:,6).';
    melt_data(i).annual_melt = melt_data(i).data(:,7).';
    melt_data(i).annual_melt_s150 = melt_data(i).data(:,8).';
    melt_data(i).annual_melt_s100 = melt_data(i).data(:,9).';
    melt_data(i).surface = melt_data(i).data(:,10).';
    melt_data(i).bed = melt_data(i).data(:,11).';
    melt_data(i).thickness = melt_data(i).data(:,12).';
    melt_data(i).SMB = melt_data(i).data(:,13).';
    melt_data(i).Dynamic_thinning = melt_data(i).data(:,14).';
end
   
% produce smoothed versions of melt signals for annual and total
for i = 1:numel(melt_CSV)
    melt_data(i).total_melt_s10 = smooth(melt_data(i).total_melt,10).';
    melt_data(i).annual_melt_s10 = smooth(melt_data(i).annual_melt,10).';
    melt_data(i).annual_melt_s150_s10 = smooth(melt_data(i).annual_melt_s150,10).';
    melt_data(i).annual_melt_s100_s10 = smooth(melt_data(i).annual_melt_s100,10).';
    
    % 20 data point smoothing
    melt_data(i).total_melt_s20 = smooth(melt_data(i).total_melt,20).';
    melt_data(i).annual_melt_s20 = smooth(melt_data(i).annual_melt,20).';
    melt_data(i).annual_melt_s150_s20 = smooth(melt_data(i).annual_melt_s150,20).';
    melt_data(i).annual_melt_s100_s20 = smooth(melt_data(i).annual_melt_s100,20).';
    
    % 50 data point smoothing
    melt_data(i).total_melt_s50 = smooth(melt_data(i).total_melt,50).';
    melt_data(i).annual_melt_s50 = smooth(melt_data(i).annual_melt,50).';
    melt_data(i).annual_melt_s150_s50 = smooth(melt_data(i).annual_melt_s150,50).';
    melt_data(i).annual_melt_s100_s50 = smooth(melt_data(i).annual_melt_s100,50).';
end

% Produce Export arrays with different 
for i = 1:numel(melt_CSV)
    export_melt_s10(i).data_frame = cat(2, melt_data(i).AT.', melt_data(i).Lons.', melt_data(i).Lats.',...
        melt_data(i).X.', melt_data(i).Y.', melt_data(i).total_melt_s10.', melt_data(i).annual_melt_s10.',...
        melt_data(i).annual_melt_s150_s10.', melt_data(i).annual_melt_s100_s10.', melt_data(i).surface.', ...
        melt_data(i).bed.', melt_data(i).thickness.', melt_data(i).SMB.', melt_data(i).Dynamic_thinning.');
    
    export_melt_s20(i).data_frame = cat(2, melt_data(i).AT.', melt_data(i).Lons.', melt_data(i).Lats.',...
        melt_data(i).X.', melt_data(i).Y.', melt_data(i).total_melt_s20.', melt_data(i).annual_melt_s20.',...
        melt_data(i).annual_melt_s150_s20.', melt_data(i).annual_melt_s100_s20.', melt_data(i).surface.', ...
        melt_data(i).bed.', melt_data(i).thickness.', melt_data(i).SMB.', melt_data(i).Dynamic_thinning.');
    
    export_melt_s50(i).data_frame = cat(2, melt_data(i).AT.', melt_data(i).Lons.', melt_data(i).Lats.',...
        melt_data(i).X.', melt_data(i).Y.', melt_data(i).total_melt_s50.', melt_data(i).annual_melt_s50.',...
        melt_data(i).annual_melt_s150_s50.', melt_data(i).annual_melt_s100_s50.',melt_data(i).surface.', ...
        melt_data(i).bed.', melt_data(i).thickness.', melt_data(i).SMB.', melt_data(i).Dynamic_thinning.');
end

% Resample each exportfile to every 10, 20, and 50 points 
for i = 1:numel(melt_CSV)
    export_file(i).matrix_10 =  export_melt_s10(i).data_frame(1:10:end,:);
    export_file(i).matrix_20 =  export_melt_s20(i).data_frame(1:20:end,:);  
    export_file(i).matrix_50 =  export_melt_s50(i).data_frame(1:50:end,:);    
end

%     melt_data(i).short_2 = melt_data(i).data(1:2:end,:);
%     melt_data(i).short_4 = melt_data(i).data(1:4:end,:);
%     melt_data(i).short_8 = melt_data(i).data(1:8:end,:);
%     melt_data(i).name = melt_CSV(i).name;


%% Define header and fields, export each resampled file into new shortened file
cheader = {'AT Distance','Lons', 'Lats', 'X', 'Y', 'smoothed total melt', 'smoothed annual melt', 'smoothed annual 150 melt', 'smoothed annual 100 melt', 'surface','Bed', 'Thickness', 'SMB', 'Dynamic thinning'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

% change folder
cd 'C:\Users\c262b531\Documents\Thesis related documents\Dynamic Thinning Correction\Finished Melt CSV\P2'

% Change file name per file prior to exporting, comment out excess file
% names if there are fewer years than 5 being exported to CSV

% 2-point Shortened file
for j = 1:numel(melt_CSV)
    name(j).file = append(string(melt_CSV(j).name),'_10s_point','.csv');
    fid = fopen(name(j).file,'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    dlmwrite(name(j).file, export_file(j).matrix_10, '-append');
end

% 4-point Shortened file
for j = 1:numel(melt_CSV)
    name(j).file = append(string(melt_CSV(j).name),'_20s_point','.csv');
    fid = fopen(name(j).file,'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    dlmwrite(name(j).file, export_file(j).matrix_20, '-append');
end

% 8-point Shortened file
for j = 1:numel(melt_CSV)
    name(j).file = append(string(melt_CSV(j).name),'_50s_point','.csv');
    fid = fopen(name(j).file,'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    dlmwrite(name(j).file, export_file(j).matrix_50, '-append');
end
