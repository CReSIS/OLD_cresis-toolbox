% function Basal_crevasse_paths(input_string)

% This function loads in point data corresponding to each flight line
% grouping for calculation of annual rate of crevasse formation. Each
% if/end statement loads individual spreadsheets for the chosen flight line
% for further use in multipass.Basal_crevasse_Formation_rates. 

%--------------------------------------------------------------------------
if strcmpi(input_string,'Petermann_line1')
    file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Extended_Passes\');
    myfiles = pwd; 
    Velocity_CSV_1 = dir(fullfile(myfiles,'P1M_div*.csv'));
        for i = 1:numel(Velocity_CSV_1)
            import_profile = fullfile(myfiles,Velocity_CSV_1(i).name);
            XY_datas(i).data = readmatrix(import_profile);
            XY_datas(i).fid = XY_datas(i).data(1:4436,1).';
            XY_datas(i).X = XY_datas(i).data(1:4436,2).';
            XY_datas(i).Y = XY_datas(i).data(1:4436,3).';
            XY_datas(i).Lons = XY_datas(i).data(1:4436,4).';
            XY_datas(i).Lats = XY_datas(i).data(1:4436,5).';
            XY_datas(i).Surf = XY_datas(i).data(1:4436,6).';
            XY_datas(i).Bed = XY_datas(i).data(1:4436,7).';
            XY_datas(i).Thickness = XY_datas(i).data(1:4436,8).';
            XY_datas(i).AT_dist = XY_datas(i).data(1:4436,9).';
            XY_datas(i).Vel_mag = XY_datas(i).data(1:4436,10).';
            XY_datas(i).div_09_10_2000m = XY_datas(i).data(1:4436,11).';
            XY_datas(i).div_12_13_2000m = XY_datas(i).data(1:4436,12).';
            XY_datas(i).div_14_15_2000m = XY_datas(i).data(1:4436,13).';
            XY_datas(i).div_15_16_2000m = XY_datas(i).data(1:4436,14).';
            XY_datas(i).div_16_17_2000m = XY_datas(i).data(1:4436,15).';
            XY_datas(i).div_17_18_2000m = XY_datas(i).data(1:4436,16).';
            XY_datas(i).SMB_10 = XY_datas(i).data(1:4436,17).';
            XY_datas(i).SMB_11 = XY_datas(i).data(1:4436,18).';
            XY_datas(i).SMB_12 = XY_datas(i).data(1:4436,19).';
            XY_datas(i).SMB_13 = XY_datas(i).data(1:4436,20).';
            XY_datas(i).SMB_14 = XY_datas(i).data(1:4436,21).';
            XY_datas(i).SMB_15 = XY_datas(i).data(1:4436,22).';        
            XY_datas(i).SMB_16 = XY_datas(i).data(1:4436,23).';
            XY_datas(i).SMB_17 = XY_datas(i).data(1:4436,24).';
            XY_datas(i).SMB_18 = XY_datas(i).data(1:4436,25).';
            XY_datas(i).SMB_19 = XY_datas(i).data(1:4436,26).';       
            XY_datas(i).AT_vel_corrected = XY_datas(1).AT_dist + XY_datas(i).Vel_mag;
            XY_datas(i).name_parts = strsplit(string(Velocity_CSV_1(i).name), '_');
            XY_datas(i).name_parts2 = strsplit(string(XY_datas(i).name_parts(4)),'.');
            XY_datas(i).year = str2double(XY_datas(i).name_parts2(1));
        end
end
%     file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Melt_CSV\');
%     myfiles = pwd; 
%     Velocity_CSV_1 = dir(fullfile(myfiles,'P1M_melt_E_new_*.csv'));
%         for i = 1:numel(Velocity_CSV_1)
%             import_profile = fullfile(myfiles,Velocity_CSV_1(i).name);
%             XY_datas(i).data = readmatrix(import_profile);
%             XY_datas(i).AT_dist = XY_datas(i).data(1:end,1).';
%             XY_datas(i).Lons = XY_datas(i).data(1:end,2).';
%             XY_datas(i).Lats = XY_datas(i).data(1:end,3).';
%             XY_datas(i).X = XY_datas(i).data(1:end,4).';
%             XY_datas(i).Y = XY_datas(i).data(1:end,5).';
%             XY_datas(i).Vel_mag = XY_datas(i).data(1:end,6).';
%             XY_datas(i).SMB_10 = XY_datas(i).data(1:end,7).';
%             XY_datas(i).SMB_11 = XY_datas(i).data(1:end,8).';
%             XY_datas(i).SMB_12 = XY_datas(i).data(1:end,9).';
%             XY_datas(i).SMB_13 = XY_datas(i).data(1:end,10).';
%             XY_datas(i).SMB_14 = XY_datas(i).data(1:end,11).';
%             XY_datas(i).SMB_15 = XY_datas(i).data(1:end,12).';        
%             XY_datas(i).SMB_16 = XY_datas(i).data(1:end,13).';
%             XY_datas(i).SMB_17 = XY_datas(i).data(1:end,14).';
%             XY_datas(i).SMB_18 = XY_datas(i).data(1:end,15).';
%             XY_datas(i).SMB_19 = XY_datas(i).data(1:end,16).';  
%             XY_datas(i).total_melt = XY_datas(i).data(1:end,17).';
%             XY_datas(i).annual_melt = XY_datas(i).data(1:end,18).'; 
%             XY_datas(i).annual_melt_150s = XY_datas(i).data(1:end,19).';
%             XY_datas(i).annual_melt_100s = XY_datas(i).data(1:end,20).';  
%             XY_datas(i).Surf = XY_datas(i).data(1:end,21).';
%             XY_datas(i).Bed = XY_datas(i).data(1:end,22).';
%             XY_datas(i).Thickness = XY_datas(i).data(1:end,23).';
%             XY_datas(i).SMB_component = XY_datas(i).data(1:end,24).';  
%             XY_datas(i).dynamic_thinning_component = XY_datas(i).data(1:end,25).';
%             XY_datas(i).AT_vel_corrected = XY_datas(1).AT_dist + XY_datas(i).Vel_mag;
%             XY_datas(i).name_parts = strsplit(string(Velocity_CSV_1(i).name), '_');
%             XY_datas(i).name_parts2 = strsplit(string(XY_datas(i).name_parts(4)),'.');
%             XY_datas(i).year = str2double(XY_datas(i).name_parts(5));
%         end        

