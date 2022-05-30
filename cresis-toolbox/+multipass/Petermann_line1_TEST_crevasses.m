%% Clip all profiles to grounding line and x-axis query
% Run this line once or uncomment to run again
% Bed
AT_data.interp_data.P11 = AT_data.interp_data.P11(479:end);
AT_data.interp_data.P14 = AT_data.interp_data.P14(479:end);
AT_data.interp_data.P18 = AT_data.interp_data.P18(479:end);
AT_data.query_array.P11 = AT_data.query_array.P11(479:end);

% Surface
AT_data.interp_data.P11_Surf = AT_data.interp_data.P11_Surf(479:end);
AT_data.interp_data.P14_Surf = AT_data.interp_data.P14_Surf(479:end);
AT_data.interp_data.P18_Surf = AT_data.interp_data.P18_Surf(479:end);

% Ice Thickness
AT_data.ice_thickness_interp.P11 = AT_data.ice_thickness_interp.P11(479:end);
AT_data.ice_thickness_interp.P14 = AT_data.ice_thickness_interp.P14(479:end);
AT_data.ice_thickness_interp.P18 = AT_data.ice_thickness_interp.P18(479:end);

% Velocity data
AT_data.interp_data.P11_VEL = AT_data.interp_data.P11_VEL(479:end);
AT_data.interp_data.P14_VEL = AT_data.interp_data.P14_VEL(479:end);
AT_data.interp_data.P18_VEL = AT_data.interp_data.P18_VEL(479:end);

% Latitude
AT_data.interp_data.P11_LAT = AT_data.interp_data.P11_LAT(479:end);
AT_data.interp_data.P14_LAT = AT_data.interp_data.P14_LAT(479:end);
AT_data.interp_data.P18_LAT = AT_data.interp_data.P18_LAT(479:end);
 
% Longitude
AT_data.interp_data.P11_LON = AT_data.interp_data.P11_LON(479:end);
AT_data.interp_data.P14_LON = AT_data.interp_data.P14_LON(479:end);
AT_data.interp_data.P18_LON = AT_data.interp_data.P18_LON(479:end);

% Distance
AT_data.interp_data.P11_DIST = AT_data.interp_data.P11_DIST(479:end);
AT_data.interp_data.P14_DIST = AT_data.interp_data.P14_DIST(479:end);
AT_data.interp_data.P18_DIST = AT_data.interp_data.P18_DIST(479:end);
%% Concatenate Lat/Lon into one variable and export to csv
% run if needing to export, otherwise skip step
% make into one table and export to csv
AT_data.lat_lon_cat.p11 = cat(2,AT_data.interp_data.P11_LON.', AT_data.interp_data.P11_LAT.');
AT_data.lat_lon_cat.p14 = cat(2,AT_data.interp_data.P14_LON.', AT_data.interp_data.P14_LAT.');
AT_data.lat_lon_cat.p18 = cat(2,AT_data.interp_data.P18_LON.', AT_data.interp_data.P18_LAT.');
writematrix(AT_data.lat_lon_cat.p11, 'clipped_latlon_2011.csv','Delimiter','comma')
writematrix(AT_data.lat_lon_cat.p14, 'clipped_latlon_2014.csv','Delimiter','comma')
writematrix(AT_data.lat_lon_cat.p18, 'clipped_latlon_2018.csv','Delimiter','comma')
%% Export of difference in ice thickness data to CSV
% data is clipped so that the grounding line is not in profile, melt rates
% include column 1) 2011-2014, 2014-2018, and 2011-2018, also the average
% anual melt rate rate between 2011-2018 as the last column
% export of lat/lon data is from the master profile, thus 2014.
% ice thickness data and annual change between, clipped to idx 917
AT_data.thickness_change.P11_14 = AT_data.interp_data.P14 - AT_data.interp_data.P11;
AT_data.thickness_change.P14_18 = AT_data.interp_data.P18 - AT_data.interp_data.P14;
AT_data.thickness_change.P11_18 = AT_data.interp_data.P18 - AT_data.interp_data.P11;
AT_data.thickness_change.P11_18_annual = AT_data.thickness_change.P11_18./7;
% index clipping
AT_data.thickness_change.P11_14 = AT_data.thickness_change.P11_14(917:end);
AT_data.thickness_change.P14_18 = AT_data.thickness_change.P14_18(917:end);
AT_data.thickness_change.P11_18 = AT_data.thickness_change.P11_18(917:end);
AT_data.thickness_change.P11_18_annual = AT_data.thickness_change.P11_18_annual(917:end);
AT_data.interp_data.P14_LAT = AT_data.interp_data.P14_LAT(917:end);
AT_data.interp_data.P14_LON = AT_data.interp_data.P14_LON(917:end);
export_file = cat(2, AT_data.interp_data.P14_LON.', AT_data.interp_data.P14_LAT.', AT_data.thickness_change.P11_14.', ...
    AT_data.thickness_change.P14_18.', AT_data.thickness_change.P11_18.', AT_data.thickness_change.P11_18_annual.');
writematrix(export_file, 'melt_data_for_ARCMAP.csv','Delimiter','comma');
%% read in velocity arrays and correct length
% read excel in and then take transpose
AT_data.extracted_vel.p11 = xlsread('2011_vel.xlsx');
AT_data.extracted_vel.p14 = xlsread('2014_vel.xlsx');
AT_data.extracted_vel.p18 = xlsread('2018_vel.xlsx');
AT_data.extracted_vel.p11 = AT_data.extracted_vel.p11.';
AT_data.extracted_vel.p14 = AT_data.extracted_vel.p14.';
AT_data.extracted_vel.p18 = AT_data.extracted_vel.p18.';
% find difference in array size and make nan_pad array for padding
size_dif = length(AT_data.ice_thickness_interp.P11)- length(AT_data.extracted_vel.p11);
AT_data.nan_pad = NaN(1,size_dif);
% concatenate pad array and original arrays together
AT_data.interp_data.p11_vel_corrected = cat(2,AT_data.extracted_vel.p11,AT_data.nan_pad);
AT_data.interp_data.p14_vel_corrected = cat(2,AT_data.extracted_vel.p14,AT_data.nan_pad);
AT_data.interp_data.p18_vel_corrected = cat(2,AT_data.extracted_vel.p18,AT_data.nan_pad);
%% Basal Crevasse Apex Picking section:
% Crevasse Apex calculation and binary file formation
AT_data.crevasse.apex_bin.P11 = islocalmax(AT_data.interp_data.P11);
AT_data.crevasse.apex_bin.P14 = islocalmax(AT_data.interp_data.P14);
AT_data.crevasse.apex_bin.P18 = islocalmax(AT_data.interp_data.P18);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
AT_data.crevasse.apex_pic_P11 = AT_data.interp_data.P11.*AT_data.crevasse.apex_bin.P11;
AT_data.crevasse.apex_pic_P11(AT_data.crevasse.apex_pic_P11 == 0) = NaN;
AT_data.crevasse.apex_pic_P14 = AT_data.interp_data.P14.*AT_data.crevasse.apex_bin.P14;
AT_data.crevasse.apex_pic_P14(AT_data.crevasse.apex_pic_P14 == 0) = NaN;
AT_data.crevasse.apex_pic_P18 = AT_data.interp_data.P18.*AT_data.crevasse.apex_bin.P18;
AT_data.crevasse.apex_pic_P18(AT_data.crevasse.apex_pic_P18 == 0) = NaN;

%% Basal Crevasse Wall Base Picking:
% Crevasse wall base calculation and binary file formation
AT_data.crevasse.base_bin.P11 = islocalmin(AT_data.interp_data.P11);
AT_data.crevasse.base_bin.P14 = islocalmin(AT_data.interp_data.P14);
AT_data.crevasse.base_bin.P18 = islocalmin(AT_data.interp_data.P18);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
AT_data.crevasse.base_pic_P11 = AT_data.interp_data.P11.*AT_data.crevasse.base_bin.P11;
AT_data.crevasse.base_pic_P11(AT_data.crevasse.base_pic_P11 == 0) = NaN;
AT_data.crevasse.base_pic_P14 = AT_data.interp_data.P14.*AT_data.crevasse.base_bin.P14;
AT_data.crevasse.base_pic_P14(AT_data.crevasse.base_pic_P14 == 0) = NaN;
AT_data.crevasse.base_pic_P18 = AT_data.interp_data.P18.*AT_data.crevasse.base_bin.P18;
AT_data.crevasse.base_pic_P18(AT_data.crevasse.base_pic_P18 == 0) = NaN;

%% AFTER MANUAL EDITTING EXPORT AND SAVE AS CSV
% after mannually editting of each profile, save as new variable and export
% variable as a csv 
% DO NO DO IF CSV ALREADY EXISTS
% COMMENT SECTION OUT AFTER FINALS PICKS AS TO NO OVERWRITE EXISTING FILES
crevasse_2011 = AT_data.crevasse.base_pic_P11;
crevasse_2014 = AT_data.crevasse.base_pic_P14;
crevasse_2018 = AT_data.crevasse.base_pic_P18;
writematrix(crevasse_2011, 'P1_Crevasse_2011.csv','Delimiter','comma')
writematrix(crevasse_2014, 'P1_Crevasse_2014.csv','Delimiter','comma')
writematrix(crevasse_2018, 'P1_Crevasse_2018.csv','Delimiter','comma')
%% Reload new csv in after editting
AT_data.crevasse.base_pic_P11 = csvread('P1_Crevasse_2011.csv');
AT_data.crevasse.base_pic_P14 = csvread('P1_Crevasse_2014.csv');
AT_data.crevasse.base_pic_P18 = csvread('P1_Crevasse_2018.csv');
%% 2011 Basal Crevasse picks
AT_data.crevasse_pick.P11.crevasse1 = AT_data.interp_data.P11(27:56);
AT_data.crevasse_pick.P11.crevasse2 = AT_data.interp_data.P11(56:70);
AT_data.crevasse_pick.P11.crevasse3_1 = AT_data.interp_data.P11(70:104);
AT_data.crevasse_pick.P11.crevasse3_2 = AT_data.interp_data.P11(104:117); %117 or 124 end
AT_data.crevasse_pick.P11.crevasse3_3 = AT_data.interp_data.P11(128:151); %124 or 128 start
AT_data.crevasse_pick.P11.crevasse4 = AT_data.interp_data.P11(151:172);
AT_data.crevasse_pick.P11.crevasse5_1 = AT_data.interp_data.P11(182:205); %was 216
AT_data.crevasse_pick.P11.crevasse5_2 = AT_data.interp_data.P11(216:232);
AT_data.crevasse_pick.P11.crevasse5_3 = AT_data.interp_data.P11(240:250);
AT_data.crevasse_pick.P11.crevasse6_1 = AT_data.interp_data.P11(258:274);
AT_data.crevasse_pick.P11.crevasse6_2 = AT_data.interp_data.P11(286:318);
AT_data.crevasse_pick.P11.crevasse7 = AT_data.interp_data.P11(324:339);
AT_data.crevasse_pick.P11.crevasse8 = AT_data.interp_data.P11(339:359);
AT_data.crevasse_pick.P11.crevasse9 = AT_data.interp_data.P11(374:396);
AT_data.crevasse_pick.P11.crevasse10 = AT_data.interp_data.P11(396:408);
%% 2014 Basal Crevasse picks
AT_data.crevasse_pick.P14.crevasse1 = AT_data.interp_data.P14(21:48);
AT_data.crevasse_pick.P14.crevasse2 = AT_data.interp_data.P14(53:86);
AT_data.crevasse_pick.P14.crevasse3 = AT_data.interp_data.P14(86:148); %121 where combined hump ends
AT_data.crevasse_pick.P14.crevasse4 = AT_data.interp_data.P14(154:186);
AT_data.crevasse_pick.P14.crevasse5_1 = AT_data.interp_data.P14(198:215);
AT_data.crevasse_pick.P14.crevasse5_2 = AT_data.interp_data.P14(215:249);
AT_data.crevasse_pick.P14.crevasse6_1 = AT_data.interp_data.P14(249:273);
AT_data.crevasse_pick.P14.crevasse6_2 = AT_data.interp_data.P14(280:297);
AT_data.crevasse_pick.P14.crevasse7 = AT_data.interp_data.P14(304:318);
AT_data.crevasse_pick.P14.crevasse8 = AT_data.interp_data.P14(325:360);
AT_data.crevasse_pick.P14.crevasse9 = AT_data.interp_data.P14(378:403); %crevasse 9 none-existant in 2018 data
AT_data.crevasse_pick.P14.crevasse10 = AT_data.interp_data.P14(:);

%% 2018 Basal Crevasse Picks
AT_data.crevasse_pick.P18.crevasse1 = AT_data.interp_data.P18(39:67);
AT_data.crevasse_pick.P18.crevasse2 = AT_data.interp_data.P18(67:105);
AT_data.crevasse_pick.P18.crevasse3 = AT_data.interp_data.P18(118:168);
AT_data.crevasse_pick.P18.crevasse4 = AT_data.interp_data.P18(176:214);
AT_data.crevasse_pick.P18.crevasse5 = AT_data.interp_data.P18(214:251);
AT_data.crevasse_pick.P18.crevasse6 = AT_data.interp_data.P18(257:292); 
AT_data.crevasse_pick.P18.crevasse7 = AT_data.interp_data.P18(319:338);
AT_data.crevasse_pick.P18.crevasse8 = AT_data.interp_data.P18(338:382);
%AT_data.crevasse_pick.P18.crevasse9 = AT_data.interp_data.P18(:); % does not exist 
AT_data.crevasse_pick.P18.crevasse10 = AT_data.interp_data.P18(:);


%% Find max values 2011
AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse1(AT_data.crevasse_pick.P11.crevasse1...
    == max(AT_data.crevasse_pick.P11.crevasse1));
AT_data.crevasse_max_idx_C1_11 = find(AT_data.crevasse_pick.P11.crevasse1 == AT_data.crevasse_max_val)+27;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse2(AT_data.crevasse_pick.P11.crevasse2...
    == max(AT_data.crevasse_pick.P11.crevasse2));
AT_data.crevasse_max_idx_C2_11 = find(AT_data.crevasse_pick.P11.crevasse2 == AT_data.crevasse_max_val)+56;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse3_1(AT_data.crevasse_pick.P11.crevasse3_1...
    == max(AT_data.crevasse_pick.P11.crevasse3_1));
AT_data.crevasse_max_idx_C3_1_11 = find(AT_data.crevasse_pick.P11.crevasse3_1 == AT_data.crevasse_max_val)+70;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse3_2(AT_data.crevasse_pick.P11.crevasse3_2...
    == max(AT_data.crevasse_pick.P11.crevasse3_2));
AT_data.crevasse_max_idx_C3_2_11 = find(AT_data.crevasse_pick.P11.crevasse3_2 == AT_data.crevasse_max_val)+104;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse3_3(AT_data.crevasse_pick.P11.crevasse3_3...
    == max(AT_data.crevasse_pick.P11.crevasse3_3));
AT_data.crevasse_max_idx_C3_3_11 = find(AT_data.crevasse_pick.P11.crevasse3_3 == AT_data.crevasse_max_val)+128;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse4(AT_data.crevasse_pick.P11.crevasse4...
    == max(AT_data.crevasse_pick.P11.crevasse4));
AT_data.crevasse_max_idx_C4_11 = find(AT_data.crevasse_pick.P11.crevasse4 == AT_data.crevasse_max_val)+151;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse5_1(AT_data.crevasse_pick.P11.crevasse5_1...
    == max(AT_data.crevasse_pick.P11.crevasse5_1));
AT_data.crevasse_max_idx_C5_1_11 = find(AT_data.crevasse_pick.P11.crevasse5_1 == AT_data.crevasse_max_val)+182;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse5_2(AT_data.crevasse_pick.P11.crevasse5_2...
    == max(AT_data.crevasse_pick.P11.crevasse5_2));
AT_data.crevasse_max_idx_C5_2_11 = find(AT_data.crevasse_pick.P11.crevasse5_2 == AT_data.crevasse_max_val)+216;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse5_3(AT_data.crevasse_pick.P11.crevasse5_3...
    == max(AT_data.crevasse_pick.P11.crevasse5_3));
AT_data.crevasse_max_idx_C5_3_11 = find(AT_data.crevasse_pick.P11.crevasse5_3 == AT_data.crevasse_max_val)+240;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse6_1(AT_data.crevasse_pick.P11.crevasse6_1...
    == max(AT_data.crevasse_pick.P11.crevasse6_1));
AT_data.crevasse_max_idx_C6_1_11 = find(AT_data.crevasse_pick.P11.crevasse6_1 == AT_data.crevasse_max_val)+258;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse6_2(AT_data.crevasse_pick.P11.crevasse6_2...
    == max(AT_data.crevasse_pick.P11.crevasse6_2));
AT_data.crevasse_max_idx_C6_2_11 = find(AT_data.crevasse_pick.P11.crevasse6_2 == AT_data.crevasse_max_val)+286;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse7(AT_data.crevasse_pick.P11.crevasse7...
    == max(AT_data.crevasse_pick.P11.crevasse7));
AT_data.crevasse_max_idx_C7_11 = find(AT_data.crevasse_pick.P11.crevasse7 == AT_data.crevasse_max_val)+324;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse8(AT_data.crevasse_pick.P11.crevasse8...
    == max(AT_data.crevasse_pick.P11.crevasse8));
AT_data.crevasse_max_idx_C8_11 = find(AT_data.crevasse_pick.P11.crevasse8 == AT_data.crevasse_max_val)+339;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P11.crevasse9(AT_data.crevasse_pick.P11.crevasse9...
    == max(AT_data.crevasse_pick.P11.crevasse9));
AT_data.crevasse_max_idx_C9_11 = find(AT_data.crevasse_pick.P11.crevasse9 == AT_data.crevasse_max_val)+374;

%% Find max values 2014
AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse1(AT_data.crevasse_pick.P14.crevasse1...
    == max(AT_data.crevasse_pick.P14.crevasse1));
AT_data.crevasse_max_idx_C1_14 = find(AT_data.crevasse_pick.P14.crevasse1 == AT_data.crevasse_max_val)+21;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse2(AT_data.crevasse_pick.P14.crevasse2...
    == max(AT_data.crevasse_pick.P14.crevasse2));
AT_data.crevasse_max_idx_C2_14 = find(AT_data.crevasse_pick.P14.crevasse2 == AT_data.crevasse_max_val)+53;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse3(AT_data.crevasse_pick.P14.crevasse3...
    == max(AT_data.crevasse_pick.P14.crevasse3));
AT_data.crevasse_max_idx_C3_14 = find(AT_data.crevasse_pick.P14.crevasse3 == AT_data.crevasse_max_val)+86;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse4(AT_data.crevasse_pick.P14.crevasse4...
    == max(AT_data.crevasse_pick.P14.crevasse4));
AT_data.crevasse_max_idx_C4_14 = find(AT_data.crevasse_pick.P14.crevasse4 == AT_data.crevasse_max_val)+148;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse5_1(AT_data.crevasse_pick.P14.crevasse5_1...
    == max(AT_data.crevasse_pick.P14.crevasse5_1));
AT_data.crevasse_max_idx_C5_1_14 = find(AT_data.crevasse_pick.P14.crevasse5_1 == AT_data.crevasse_max_val)+186;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse5_2(AT_data.crevasse_pick.P14.crevasse5_2...
    == max(AT_data.crevasse_pick.P14.crevasse5_2));
AT_data.crevasse_max_idx_C5_2_14 = find(AT_data.crevasse_pick.P14.crevasse5_2 == AT_data.crevasse_max_val)+215;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse6_1(AT_data.crevasse_pick.P14.crevasse6_1...
    == max(AT_data.crevasse_pick.P14.crevasse6_1));
AT_data.crevasse_max_idx_C6_1_14 = find(AT_data.crevasse_pick.P14.crevasse6_1 == AT_data.crevasse_max_val)+249;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse6_2(AT_data.crevasse_pick.P14.crevasse6_2...
    == max(AT_data.crevasse_pick.P14.crevasse6_2));
AT_data.crevasse_max_idx_C6_2_14 = find(AT_data.crevasse_pick.P14.crevasse6_2 == AT_data.crevasse_max_val)+280;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse7(AT_data.crevasse_pick.P14.crevasse7...
    == max(AT_data.crevasse_pick.P14.crevasse7));
AT_data.crevasse_max_idx_C7_14 = find(AT_data.crevasse_pick.P14.crevasse7 == AT_data.crevasse_max_val)+304;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse8(AT_data.crevasse_pick.P14.crevasse8...
    == max(AT_data.crevasse_pick.P14.crevasse8));
AT_data.crevasse_max_idx_C8_14 = find(AT_data.crevasse_pick.P14.crevasse8 == AT_data.crevasse_max_val)+325;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P14.crevasse9(AT_data.crevasse_pick.P14.crevasse9...
    == max(AT_data.crevasse_pick.P14.crevasse9));
AT_data.crevasse_max_idx_C9_14 = find(AT_data.crevasse_pick.P14.crevasse9 == AT_data.crevasse_max_val)+378;

%% Find max values 2018
AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse1(AT_data.crevasse_pick.P18.crevasse1...
    == max(AT_data.crevasse_pick.P18.crevasse1));
AT_data.crevasse_max_idx_C1_18 = find(AT_data.crevasse_pick.P18.crevasse1 == AT_data.crevasse_max_val)+39;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse2(AT_data.crevasse_pick.P18.crevasse2...
    == max(AT_data.crevasse_pick.P18.crevasse2));
AT_data.crevasse_max_idx_C2_18 = find(AT_data.crevasse_pick.P18.crevasse2 == AT_data.crevasse_max_val)+67;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse3(AT_data.crevasse_pick.P18.crevasse3...
    == max(AT_data.crevasse_pick.P18.crevasse3));
AT_data.crevasse_max_idx_C3_18 = find(AT_data.crevasse_pick.P18.crevasse3 == AT_data.crevasse_max_val)+118;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse4(AT_data.crevasse_pick.P18.crevasse4...
    == max(AT_data.crevasse_pick.P18.crevasse4));
AT_data.crevasse_max_idx_C4_18 = find(AT_data.crevasse_pick.P18.crevasse4 == AT_data.crevasse_max_val)+176;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse5(AT_data.crevasse_pick.P18.crevasse5...
    == max(AT_data.crevasse_pick.P18.crevasse5));
AT_data.crevasse_max_idx_C5_18 = find(AT_data.crevasse_pick.P18.crevasse5 == AT_data.crevasse_max_val)+214;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse6(AT_data.crevasse_pick.P18.crevasse6...
    == max(AT_data.crevasse_pick.P18.crevasse6));
AT_data.crevasse_max_idx_C6_18 = find(AT_data.crevasse_pick.P18.crevasse6 == AT_data.crevasse_max_val)+257;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse7(AT_data.crevasse_pick.P18.crevasse7...
    == max(AT_data.crevasse_pick.P18.crevasse7));
AT_data.crevasse_max_idx_C7_18 = find(AT_data.crevasse_pick.P18.crevasse7 == AT_data.crevasse_max_val)+319;

AT_data.crevasse_max_val = AT_data.crevasse_pick.P18.crevasse8(AT_data.crevasse_pick.P18.crevasse8...
    == max(AT_data.crevasse_pick.P18.crevasse8));
AT_data.crevasse_max_idx_C8_18 = find(AT_data.crevasse_pick.P18.crevasse8 == AT_data.crevasse_max_val)+338;

%% 2011 Melt crevasse apex values from total melt
AT_data.apex_melt.P11.C1 = [total_melt1(AT_data.crevasse_max_idx_C1_11)];
AT_data.apex_melt.P11.C2 = [total_melt1(AT_data.crevasse_max_idx_C2_11)];
AT_data.apex_melt.P11.C3 = [total_melt1(AT_data.crevasse_max_idx_C3_1_11), total_melt1(AT_data.crevasse_max_idx_C3_2_11),...
    total_melt1(AT_data.crevasse_max_idx_C3_3_11)];
AT_data.apex_melt.P11.C4 = [total_melt1(AT_data.crevasse_max_idx_C4_11)];
AT_data.apex_melt.P11.C5 = [total_melt1(AT_data.crevasse_max_idx_C5_1_11), total_melt1(AT_data.crevasse_max_idx_C5_2_11),...
    total_melt1(AT_data.crevasse_max_idx_C5_3_11)];
AT_data.apex_melt.P11.C6 = [total_melt1(AT_data.crevasse_max_idx_C6_1_11), total_melt1(AT_data.crevasse_max_idx_C6_2_11)];
AT_data.apex_melt.P11.C7 = [total_melt1(AT_data.crevasse_max_idx_C7_11)];
AT_data.apex_melt.P11.C8 = [total_melt1(AT_data.crevasse_max_idx_C8_11)];
AT_data.apex_melt.P11.C9 = [total_melt1(AT_data.crevasse_max_idx_C9_11)];

%% 2014 Melt crevasse apex values from total melt
AT_data.apex_melt.P14.C1 = [total_melt1(AT_data.crevasse_max_idx_C1_14)];
AT_data.apex_melt.P14.C2 = [total_melt1(AT_data.crevasse_max_idx_C2_14)];
AT_data.apex_melt.P14.C3 = [total_melt1(AT_data.crevasse_max_idx_C3_14)];
AT_data.apex_melt.P14.C4 = [total_melt1(AT_data.crevasse_max_idx_C4_14)];
AT_data.apex_melt.P14.C5 = [total_melt1(AT_data.crevasse_max_idx_C5_1_14), total_melt1(AT_data.crevasse_max_idx_C5_2_14)];
AT_data.apex_melt.P14.C6 = [total_melt1(AT_data.crevasse_max_idx_C6_1_14), total_melt1(AT_data.crevasse_max_idx_C6_2_14)];
AT_data.apex_melt.P14.C7 = [total_melt1(AT_data.crevasse_max_idx_C7_14)];
AT_data.apex_melt.P14.C8 = [total_melt1(AT_data.crevasse_max_idx_C8_14)];
AT_data.apex_melt.P14.C9 = [total_melt1(AT_data.crevasse_max_idx_C9_14)];

%% 2018 Melt crevasse apex values from total melt
AT_data.apex_melt.P18.C1 = [total_melt1(AT_data.crevasse_max_idx_C1_18)];
AT_data.apex_melt.P18.C2 = [total_melt1(AT_data.crevasse_max_idx_C2_18)];
AT_data.apex_melt.P18.C3 = [total_melt1(AT_data.crevasse_max_idx_C3_18)];
AT_data.apex_melt.P18.C4 = [total_melt1(AT_data.crevasse_max_idx_C4_18)];
AT_data.apex_melt.P18.C5 = [total_melt1(AT_data.crevasse_max_idx_C5_18)];
AT_data.apex_melt.P18.C6 = [total_melt1(AT_data.crevasse_max_idx_C6_18)];
AT_data.apex_melt.P18.C7 = [total_melt1(AT_data.crevasse_max_idx_C7_18)];
AT_data.apex_melt.P18.C8 = [total_melt1(AT_data.crevasse_max_idx_C8_18)];

%% 
AT_data.apex_melt.C1 = [total_melt1(AT_data.crevasse_max_idx_C1_11), total_melt1(AT_data.crevasse_max_idx_C1_14),...
    total_melt1(AT_data.crevasse_max_idx_C1_18)];
AT_data.apex_melt.C2 = [total_melt1(AT_data.crevasse_max_idx_C2_11), total_melt1(AT_data.crevasse_max_idx_C2_14),...
    total_melt1(AT_data.crevasse_max_idx_C2_18)];
AT_data.apex_melt.C3 = [total_melt1(AT_data.crevasse_max_idx_C3_1_11), total_melt1(AT_data.crevasse_max_idx_C3_2_11),...
    total_melt1(AT_data.crevasse_max_idx_C3_3_11), total_melt1(AT_data.crevasse_max_idx_C1_14),...
    total_melt1(AT_data.crevasse_max_idx_C1_18)];
AT_data.apex_melt.C4 = [total_melt1(AT_data.crevasse_max_idx_C4_11), total_melt1(AT_data.crevasse_max_idx_C4_14),...
    total_melt1(AT_data.crevasse_max_idx_C4_18)];
AT_data.apex_melt.C5 = [total_melt1(AT_data.crevasse_max_idx_C5_1_11), total_melt1(AT_data.crevasse_max_idx_C5_2_11),...
    total_melt1(AT_data.crevasse_max_idx_C5_3_11), total_melt1(AT_data.crevasse_max_idx_C5_1_14), ...
    total_melt1(AT_data.crevasse_max_idx_C5_2_14), total_melt1(AT_data.crevasse_max_idx_C5_18)];
AT_data.apex_melt.C6 = [total_melt1(AT_data.crevasse_max_idx_C6_1_11), total_melt1(AT_data.crevasse_max_idx_C6_2_11), ...
    total_melt1(AT_data.crevasse_max_idx_C6_1_14), total_melt1(AT_data.crevasse_max_idx_C6_2_14),...
    total_melt1(AT_data.crevasse_max_idx_C4_18)];
AT_data.apex_melt.C7 = [total_melt1(AT_data.crevasse_max_idx_C7_11), total_melt1(AT_data.crevasse_max_idx_C7_14),...
    total_melt1(AT_data.crevasse_max_idx_C7_18)];
AT_data.apex_melt.C8 = [total_melt1(AT_data.crevasse_max_idx_C8_11), total_melt1(AT_data.crevasse_max_idx_C8_14),...
    total_melt1(AT_data.crevasse_max_idx_C8_18)];
AT_data.apex_melt.C9 = [total_melt1(AT_data.crevasse_max_idx_C9_11), total_melt1(AT_data.crevasse_max_idx_C9_14)];

%% Along track axis for individual crevasses and Normalization of each to start point 0
% crevasse 1
AT_data.along_track_crevasse.C1_11 = AT_data.query_array.P11(27:56)/1e3;
AT_data.along_track_crevasse.C1_14 = AT_data.query_array.P11(21:48)/1e3;
AT_data.along_track_crevasse.C1_18 = AT_data.query_array.P11(39:67)/1e3;
AT_data.along_track_crevasse_normalized.C1_11_n = (AT_data.along_track_crevasse.C1_11 - min(AT_data.along_track_crevasse.C1_11))*1000;
AT_data.along_track_crevasse_normalized.C1_14_n = (AT_data.along_track_crevasse.C1_14 - min(AT_data.along_track_crevasse.C1_14))*1000;
AT_data.along_track_crevasse_normalized.C1_18_n = (AT_data.along_track_crevasse.C1_18 - min(AT_data.along_track_crevasse.C1_18))*1000;

% crevasse 2
AT_data.along_track_crevasse.C2_11 = AT_data.query_array.P11(56:70)/1e3;
AT_data.along_track_crevasse.C2_14 = AT_data.query_array.P11(53:86)/1e3;
AT_data.along_track_crevasse.C2_18 = AT_data.query_array.P11(67:105)/1e3;
AT_data.along_track_crevasse_normalized.C2_11_n = (AT_data.along_track_crevasse.C2_11 - min(AT_data.along_track_crevasse.C2_11))*1000;
AT_data.along_track_crevasse_normalized.C2_14_n = (AT_data.along_track_crevasse.C2_14 - min(AT_data.along_track_crevasse.C2_14))*1000;
AT_data.along_track_crevasse_normalized.C2_18_n = (AT_data.along_track_crevasse.C2_18 - min(AT_data.along_track_crevasse.C2_18))*1000;

% crevasse 3
AT_data.along_track_crevasse.C3_1_11 = AT_data.query_array.P11(70:104)/1e3;
AT_data.along_track_crevasse.C3_2_11 = AT_data.query_array.P11(104:117)/1e3;
AT_data.along_track_crevasse.C3_3_11 = AT_data.query_array.P11(128:151)/1e3;
AT_data.along_track_crevasse.C3_14 = AT_data.query_array.P11(86:148)/1e3;
AT_data.along_track_crevasse.C3_18 = AT_data.query_array.P11(118:168)/1e3;
AT_data.along_track_crevasse_normalized.C3_1_11_n = (AT_data.along_track_crevasse.C3_1_11 - min(AT_data.along_track_crevasse.C3_1_11))*1000;
AT_data.along_track_crevasse_normalized.C3_2_11_n = ((AT_data.along_track_crevasse.C3_2_11 - min(AT_data.along_track_crevasse.C3_2_11)))*1000 + max(AT_data.along_track_crevasse_normalized.C3_1_11_n);
AT_data.along_track_crevasse_normalized.C3_3_11_n = ((AT_data.along_track_crevasse.C3_3_11 - min(AT_data.along_track_crevasse.C3_3_11)))*1000 + max(AT_data.along_track_crevasse_normalized.C3_2_11_n);
AT_data.along_track_crevasse_normalized.C3_14_n = (AT_data.along_track_crevasse.C3_14 - min(AT_data.along_track_crevasse.C3_14))*1000;
AT_data.along_track_crevasse_normalized.C3_18_n = (AT_data.along_track_crevasse.C3_18 - min(AT_data.along_track_crevasse.C3_18))*1000;


% figure(33)
% h1 = plot(AT_data.along_track_crevasse_normalized.C3_1_11_n, AT_data.crevasse_pick.P11.crevasse3_1);
% hold on
% h2 = plot(AT_data.along_track_crevasse_normalized.C3_2_11_n, AT_data.crevasse_pick.P11.crevasse3_2);
% h3 = plot(AT_data.along_track_crevasse_normalized.C3_3_11_n, AT_data.crevasse_pick.P11.crevasse3_3);
% h4 = plot(AT_data.along_track_crevasse_normalized.C3_14_n, AT_data.crevasse_pick.P14.crevasse3);
% h5 = plot(AT_data.along_track_crevasse_normalized.C3_18_n, AT_data.crevasse_pick.P18.crevasse3);
% cross_area3_1_11 = patch(AT_data.along_track_crevasse_normalized.C3_1_11_n, AT_data.crevasse_pick.P11.crevasse3_1,'r','LineWidth',0.5);
% cross_area3_2_11 = patch(AT_data.along_track_crevasse_normalized.C3_2_11_n, AT_data.crevasse_pick.P11.crevasse3_2,'r','LineWidth',0.5);
% cross_area3_3_11 = patch(AT_data.along_track_crevasse_normalized.C3_3_11_n, AT_data.crevasse_pick.P11.crevasse3_3,'r','LineWidth',0.5);
% cross_area3_14 = patch(AT_data.along_track_crevasse_normalized.C3_14_n, AT_data.crevasse_pick.P14.crevasse3,'b','LineWidth',0.5);
% cross_area3_18 = patch(AT_data.along_track_crevasse_normalized.C3_18_n, AT_data.crevasse_pick.P18.crevasse3,'g','LineWidth',0.5);
% title('Basal Crevasse 3');
% xlabel('Along Track (km)');
% ylabel('Elevation (m)');
% %legend('basal crevasse 3.1 P11', 'basal crevasse 3.2 P11', 'basal crevasse 3.3 P11', ...
% %    'basal crevasse 3 P14', 'basal crevasse 3 P18', 'Location', 'northwest');
% alpha(0.3)
% hold off


% crevasse 4
AT_data.along_track_crevasse.C4_11 = AT_data.query_array.P11(151:172)/1e3;
AT_data.along_track_crevasse.C4_14 = AT_data.query_array.P11(154:186)/1e3;
AT_data.along_track_crevasse.C4_18 = AT_data.query_array.P11(176:214)/1e3;
AT_data.along_track_crevasse_normalized.C4_11_n = (AT_data.along_track_crevasse.C4_11 - min(AT_data.along_track_crevasse.C4_11))*1000;
AT_data.along_track_crevasse_normalized.C4_14_n = (AT_data.along_track_crevasse.C4_14 - min(AT_data.along_track_crevasse.C4_14))*1000;
AT_data.along_track_crevasse_normalized.C4_18_n = (AT_data.along_track_crevasse.C4_18 - min(AT_data.along_track_crevasse.C4_18))*1000;

% crevasse 5
AT_data.along_track_crevasse.C5_1_11 = AT_data.query_array.P11(182:205)/1e3;
AT_data.along_track_crevasse.C5_2_11 = AT_data.query_array.P11(216:232)/1e3;
AT_data.along_track_crevasse.C5_3_11 = AT_data.query_array.P11(240:250)/1e3;
AT_data.along_track_crevasse.C5_1_14 = AT_data.query_array.P11(198:215)/1e3;
AT_data.along_track_crevasse.C5_2_14 = AT_data.query_array.P11(215:249)/1e3;
AT_data.along_track_crevasse.C5_18 = AT_data.query_array.P11(214:251)/1e3;
AT_data.along_track_crevasse_normalized.C5_1_11_n = (AT_data.along_track_crevasse.C5_1_11 - min(AT_data.along_track_crevasse.C5_1_11))*1000;
AT_data.along_track_crevasse_normalized.C5_2_11_n = ((AT_data.along_track_crevasse.C5_2_11 - min(AT_data.along_track_crevasse.C5_2_11)))*1000 + max(AT_data.along_track_crevasse_normalized.C5_1_11_n);
AT_data.along_track_crevasse_normalized.C5_3_11_n = ((AT_data.along_track_crevasse.C5_3_11 - min(AT_data.along_track_crevasse.C5_3_11)))*1000 + max(AT_data.along_track_crevasse_normalized.C5_2_11_n);
AT_data.along_track_crevasse_normalized.C5_1_14_n = (AT_data.along_track_crevasse.C5_1_14 - min(AT_data.along_track_crevasse.C5_1_14))*1000;
AT_data.along_track_crevasse_normalized.C5_2_14_n = ((AT_data.along_track_crevasse.C5_2_14 - min(AT_data.along_track_crevasse.C5_2_14)))*1000 + max(AT_data.along_track_crevasse_normalized.C5_1_14_n);
AT_data.along_track_crevasse_normalized.C5_18_n = (AT_data.along_track_crevasse.C5_18 - min(AT_data.along_track_crevasse.C5_18))*1000;

% crevasse 6
AT_data.along_track_crevasse.C6_1_11 = AT_data.query_array.P11(258:274)/1e3;
AT_data.along_track_crevasse.C6_2_11 = AT_data.query_array.P11(286:318)/1e3;
AT_data.along_track_crevasse.C6_1_14 = AT_data.query_array.P11(249:273)/1e3;
AT_data.along_track_crevasse.C6_2_14 = AT_data.query_array.P11(280:297)/1e3;
AT_data.along_track_crevasse.C6_18 = AT_data.query_array.P11(257:292)/1e3;
AT_data.along_track_crevasse_normalized.C6_1_11_n = (AT_data.along_track_crevasse.C6_1_11 - min(AT_data.along_track_crevasse.C6_1_11))*1000;
AT_data.along_track_crevasse_normalized.C6_2_11_n = ((AT_data.along_track_crevasse.C6_2_11 - min(AT_data.along_track_crevasse.C6_2_11)))*1000 + max(AT_data.along_track_crevasse_normalized.C6_1_11_n);
AT_data.along_track_crevasse_normalized.C6_1_14_n = (AT_data.along_track_crevasse.C6_1_14 - min(AT_data.along_track_crevasse.C6_1_14))*1000;
AT_data.along_track_crevasse_normalized.C6_2_14_n = ((AT_data.along_track_crevasse.C6_2_14 - min(AT_data.along_track_crevasse.C6_2_14)))*1000 + max(AT_data.along_track_crevasse_normalized.C6_1_14_n);
AT_data.along_track_crevasse_normalized.C6_18_n = (AT_data.along_track_crevasse.C6_18 - min(AT_data.along_track_crevasse.C6_18))*1000;

% crevasse 7
AT_data.along_track_crevasse.C7_11 = AT_data.query_array.P11(324:339)/1e3;
AT_data.along_track_crevasse.C7_14 = AT_data.query_array.P11(304:318)/1e3;
AT_data.along_track_crevasse.C7_18 = AT_data.query_array.P11(319:338)/1e3;
AT_data.along_track_crevasse_normalized.C7_11_n = (AT_data.along_track_crevasse.C7_11 - min(AT_data.along_track_crevasse.C7_11))*1000;
AT_data.along_track_crevasse_normalized.C7_14_n = (AT_data.along_track_crevasse.C7_14 - min(AT_data.along_track_crevasse.C7_14))*1000;
AT_data.along_track_crevasse_normalized.C7_18_n = (AT_data.along_track_crevasse.C7_18 - min(AT_data.along_track_crevasse.C7_18))*1000;

% crevasse 8
AT_data.along_track_crevasse.C8_11 = AT_data.query_array.P11(339:359)/1e3;
AT_data.along_track_crevasse.C8_14 = AT_data.query_array.P11(325:360)/1e3;
AT_data.along_track_crevasse.C8_18 = AT_data.query_array.P11(338:382)/1e3;
AT_data.along_track_crevasse_normalized.C8_11_n = (AT_data.along_track_crevasse.C8_11 - min(AT_data.along_track_crevasse.C8_11))*1000;
AT_data.along_track_crevasse_normalized.C8_14_n = (AT_data.along_track_crevasse.C8_14 - min(AT_data.along_track_crevasse.C8_14))*1000;
AT_data.along_track_crevasse_normalized.C8_18_n = (AT_data.along_track_crevasse.C8_18 - min(AT_data.along_track_crevasse.C8_18))*1000;

% crevasse 9
AT_data.along_track_crevasse.C9_11 = AT_data.query_array.P11(374:396)/1e3;
AT_data.along_track_crevasse.C9_14 = AT_data.query_array.P11(378:403)/1e3;
%AT_data.along_track_crevasse.C9_18 = AT_data.query_array.P11(338:382)/1e3;
AT_data.along_track_crevasse_normalized.C9_11_n = (AT_data.along_track_crevasse.C9_11 - min(AT_data.along_track_crevasse.C9_11))*1000;
AT_data.along_track_crevasse_normalized.C9_14_n = (AT_data.along_track_crevasse.C9_14 - min(AT_data.along_track_crevasse.C9_14))*1000;
%AT_data.along_track_crevasse_normalized.C9_18_n = (AT_data.along_track_crevasse.C9_18 - min(AT_data.along_track_crevasse.C9_11))*1000;

%% Subplots 3x3
figure(100)
subplot(3,3,1)
h1 = plot(AT_data.along_track_crevasse_normalized.C1_11_n, AT_data.crevasse_pick.P11.crevasse1);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C1_14_n, AT_data.crevasse_pick.P14.crevasse1);
h3 = plot(AT_data.along_track_crevasse_normalized.C1_18_n, AT_data.crevasse_pick.P18.crevasse1);
cross_area1_11 = patch(AT_data.along_track_crevasse_normalized.C1_11_n, AT_data.crevasse_pick.P11.crevasse1,'r','LineWidth',0.5);
cross_area1_14 = patch(AT_data.along_track_crevasse_normalized.C1_14_n, AT_data.crevasse_pick.P14.crevasse1,'b','LineWidth',0.5);
cross_area1_18 = patch(AT_data.along_track_crevasse_normalized.C1_18_n, AT_data.crevasse_pick.P18.crevasse1,'g','LineWidth',0.5);
title('Basal Crevasse 1');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 1 P11', 'basal crevasse 1 P14', 'basal crevasse 1 P18', 'Location', 'southoutside');
alpha(0.3)
hold off

subplot(3,3,2)
h1 = plot(AT_data.along_track_crevasse_normalized.C2_11_n, AT_data.crevasse_pick.P11.crevasse2);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C2_14_n, AT_data.crevasse_pick.P14.crevasse2);
h3 = plot(AT_data.along_track_crevasse_normalized.C2_18_n, AT_data.crevasse_pick.P18.crevasse2);
cross_area2_11 = patch(AT_data.along_track_crevasse_normalized.C2_11_n, AT_data.crevasse_pick.P11.crevasse2,'r','LineWidth',0.5);
cross_area2_14 = patch(AT_data.along_track_crevasse_normalized.C2_14_n, AT_data.crevasse_pick.P14.crevasse2,'b','LineWidth',0.5);
cross_area2_18 = patch(AT_data.along_track_crevasse_normalized.C2_18_n, AT_data.crevasse_pick.P18.crevasse2,'g','LineWidth',0.5);
title('Basal Crevasse 2');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 2 P11', 'basal crevasse 2 P14', 'basal crevasse 2 P18', 'Location', 'northwest');
alpha(0.3)
hold off

subplot(3,3,3)
h1 = plot(AT_data.along_track_crevasse_normalized.C3_1_11_n, AT_data.crevasse_pick.P11.crevasse3_1);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C3_2_11_n, AT_data.crevasse_pick.P11.crevasse3_2);
h3 = plot(AT_data.along_track_crevasse_normalized.C3_3_11_n, AT_data.crevasse_pick.P11.crevasse3_3);
h4 = plot(AT_data.along_track_crevasse_normalized.C3_14_n, AT_data.crevasse_pick.P14.crevasse3);
h5 = plot(AT_data.along_track_crevasse_normalized.C3_18_n, AT_data.crevasse_pick.P18.crevasse3);
cross_area3_1_11 = patch(AT_data.along_track_crevasse_normalized.C3_1_11_n, AT_data.crevasse_pick.P11.crevasse3_1,'r','LineWidth',0.5);
cross_area3_2_11 = patch(AT_data.along_track_crevasse_normalized.C3_2_11_n, AT_data.crevasse_pick.P11.crevasse3_2,'r','LineWidth',0.5);
cross_area3_3_11 = patch(AT_data.along_track_crevasse_normalized.C3_3_11_n, AT_data.crevasse_pick.P11.crevasse3_3,'r','LineWidth',0.5);
cross_area3_14 = patch(AT_data.along_track_crevasse_normalized.C3_14_n, AT_data.crevasse_pick.P14.crevasse3,'b','LineWidth',0.5);
cross_area3_18 = patch(AT_data.along_track_crevasse_normalized.C3_18_n, AT_data.crevasse_pick.P18.crevasse3,'g','LineWidth',0.5);
title('Basal Crevasse 3');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 3.1 P11', 'basal crevasse 3.2 P11', 'basal crevasse 3.3 P11', ...
%    'basal crevasse 3 P14', 'basal crevasse 3 P18', 'Location', 'northwest');
alpha(0.3)
hold off

subplot(3,3,4)
h1 = plot(AT_data.along_track_crevasse_normalized.C4_11_n, AT_data.crevasse_pick.P11.crevasse4);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C4_14_n, AT_data.crevasse_pick.P14.crevasse4);
h3 = plot(AT_data.along_track_crevasse_normalized.C4_18_n, AT_data.crevasse_pick.P18.crevasse4);
cross_area4_11 = patch(AT_data.along_track_crevasse_normalized.C4_11_n, AT_data.crevasse_pick.P11.crevasse4,'r','LineWidth',0.5);
cross_area4_14 = patch(AT_data.along_track_crevasse_normalized.C4_14_n, AT_data.crevasse_pick.P14.crevasse4,'b','LineWidth',0.5);
cross_area4_18 = patch(AT_data.along_track_crevasse_normalized.C4_18_n, AT_data.crevasse_pick.P18.crevasse4,'g','LineWidth',0.5);
title('Basal Crevasse 4');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 4 P11', 'basal crevasse 4 P14', 'basal crevasse 4 P18', 'Location', 'northwest');
alpha(0.3)
hold off

subplot(3,3,5)
h1 = plot(AT_data.along_track_crevasse_normalized.C5_1_11_n, AT_data.crevasse_pick.P11.crevasse5_1);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C5_2_11_n, AT_data.crevasse_pick.P11.crevasse5_2);
h3 = plot(AT_data.along_track_crevasse_normalized.C5_3_11_n, AT_data.crevasse_pick.P11.crevasse5_3);
h4 = plot(AT_data.along_track_crevasse_normalized.C5_1_14_n, AT_data.crevasse_pick.P14.crevasse5_1);
h5 = plot(AT_data.along_track_crevasse_normalized.C5_2_14_n, AT_data.crevasse_pick.P14.crevasse5_2);
h6 = plot(AT_data.along_track_crevasse_normalized.C5_18_n, AT_data.crevasse_pick.P18.crevasse5);
cross_area5_1_11 = patch(AT_data.along_track_crevasse_normalized.C5_1_11_n, AT_data.crevasse_pick.P11.crevasse5_1,'r','LineWidth',0.5);
cross_area5_2_11 = patch(AT_data.along_track_crevasse_normalized.C5_2_11_n, AT_data.crevasse_pick.P11.crevasse5_2,'r','LineWidth',0.5);
cross_area5_3_11 = patch(AT_data.along_track_crevasse_normalized.C5_3_11_n, AT_data.crevasse_pick.P11.crevasse5_3,'r','LineWidth',0.5);
cross_area5_1_14 = patch(AT_data.along_track_crevasse_normalized.C5_1_14_n, AT_data.crevasse_pick.P14.crevasse5_1,'b','LineWidth',0.5);
cross_area5_2_14 = patch(AT_data.along_track_crevasse_normalized.C5_2_14_n, AT_data.crevasse_pick.P14.crevasse5_2,'b','LineWidth',0.5);
cross_area5_18 = patch(AT_data.along_track_crevasse_normalized.C5_18_n, AT_data.crevasse_pick.P18.crevasse5,'g','LineWidth',0.5);
title('Basal Crevasse 5');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 5.1 P11', 'basal crevasse 5.2 P11', 'basal crevasse 5.3 P11',...
%    'basal crevasse 5.1 P14', 'basal crevasse 5.2 P14', 'basal crevasse 5 P18', ...
%    'Location', 'northwest');
alpha(0.3)
hold off

subplot(3,3,6)
h1 = plot(AT_data.along_track_crevasse_normalized.C6_1_11_n, AT_data.crevasse_pick.P11.crevasse6_1);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C6_2_11_n, AT_data.crevasse_pick.P11.crevasse6_2);
h3 = plot(AT_data.along_track_crevasse_normalized.C6_1_14_n, AT_data.crevasse_pick.P14.crevasse6_1);
h4 = plot(AT_data.along_track_crevasse_normalized.C6_2_14_n, AT_data.crevasse_pick.P14.crevasse6_2);
h5 = plot(AT_data.along_track_crevasse_normalized.C6_18_n, AT_data.crevasse_pick.P18.crevasse6);
cross_area6_1_11 = patch(AT_data.along_track_crevasse_normalized.C6_1_11_n, AT_data.crevasse_pick.P11.crevasse6_1,'r','LineWidth',0.5);
cross_area6_2_11 = patch(AT_data.along_track_crevasse_normalized.C6_2_11_n, AT_data.crevasse_pick.P11.crevasse6_2,'r','LineWidth',0.5);
cross_area6_1_14 = patch(AT_data.along_track_crevasse_normalized.C6_1_14_n, AT_data.crevasse_pick.P14.crevasse6_1,'b','LineWidth',0.5);
cross_area6_2_14 = patch(AT_data.along_track_crevasse_normalized.C6_2_14_n, AT_data.crevasse_pick.P14.crevasse6_2,'b','LineWidth',0.5);
cross_area6_18 = patch(AT_data.along_track_crevasse_normalized.C6_18_n, AT_data.crevasse_pick.P18.crevasse6,'g','LineWidth',0.5);
title('Basal Crevasse 6');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 6.1 P11', 'basal crevasse 6.2 P11', 'basal crevasse 6.1 P14', ...
%    'basal crevasse 6.2 P14', 'basal crevasse 6 P18', 'Location', 'northwest');
alpha(0.3)
hold off

subplot(3,3,7)
h1 = plot(AT_data.along_track_crevasse_normalized.C7_11_n, AT_data.crevasse_pick.P11.crevasse7);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C7_14_n, AT_data.crevasse_pick.P14.crevasse7);
h3 = plot(AT_data.along_track_crevasse_normalized.C7_18_n, AT_data.crevasse_pick.P18.crevasse7);
cross_area7_11 = patch(AT_data.along_track_crevasse_normalized.C7_11_n, AT_data.crevasse_pick.P11.crevasse7,'r','LineWidth',0.5);
cross_area7_14 = patch(AT_data.along_track_crevasse_normalized.C7_14_n, AT_data.crevasse_pick.P14.crevasse7,'b','LineWidth',0.5);
cross_area7_18 = patch(AT_data.along_track_crevasse_normalized.C7_18_n, AT_data.crevasse_pick.P18.crevasse7,'g','LineWidth',0.5);
title('Basal Crevasse 7');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 7 P11', 'basal crevasse 7 P14', 'basal crevasse 7 P18', 'Location', 'northwest');
alpha(0.3)
hold off

subplot(3,3,8)
h1 = plot(AT_data.along_track_crevasse_normalized.C8_11_n, AT_data.crevasse_pick.P11.crevasse8,'r','LineWidth',0.5);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C8_14_n, AT_data.crevasse_pick.P14.crevasse8, 'b','LineWidth',0.5);
h3 = plot(AT_data.along_track_crevasse_normalized.C8_18_n, AT_data.crevasse_pick.P18.crevasse8, 'g','LineWidth',0.5);
cross_area8_11 = patch(AT_data.along_track_crevasse_normalized.C8_11_n, AT_data.crevasse_pick.P11.crevasse8,'r','LineWidth',0.5);
cross_area8_14 = patch(AT_data.along_track_crevasse_normalized.C8_14_n, AT_data.crevasse_pick.P14.crevasse8,'b','LineWidth',0.5);
cross_area8_18 = patch(AT_data.along_track_crevasse_normalized.C8_18_n, AT_data.crevasse_pick.P18.crevasse8,'g','LineWidth',0.5);
title('Basal Crevasse 8');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
legend('2011 Season Basal Crevasses ', '2014 Season Basal Crevasses', '2018 Season Basal Crevasses', 'Location', 'southoutside');
alpha(0.3)
hold off

subplot(3,3,9)
h1 = plot(AT_data.along_track_crevasse_normalized.C9_11_n, AT_data.crevasse_pick.P11.crevasse9);
hold on
h2 = plot(AT_data.along_track_crevasse_normalized.C9_14_n, AT_data.crevasse_pick.P14.crevasse9);
%h3 = plot(AT_data.along_track_crevasse_normalized.C9_18_n, AT_data.crevasse_pick.P18.crevasse9);
cross_area9_11 = patch(AT_data.along_track_crevasse_normalized.C9_11_n, AT_data.crevasse_pick.P11.crevasse9,'r','LineWidth',0.5);
cross_area9_14 = patch(AT_data.along_track_crevasse_normalized.C9_14_n, AT_data.crevasse_pick.P14.crevasse9,'b','LineWidth',0.5);
%cross_area9_18 = patch(AT_data.along_track_crevasse_normalized.C9_18_n, AT_data.crevasse_pick.P18.crevasse9,'g','LineWidth',0.5);
title('Basal Crevasse 9');
xlabel('Feature Width (m)');
ylabel('Elevation (m)');
xlim([0 700])
ylim([-450 -100])
%legend('basal crevasse 9 P11', 'basal crevasse 9 P14', 'Location', 'northwest');
alpha(0.3)
hold off

suptitle('Basal Crevasse Morphology for Petermann Flight Line 1')

%% Integrated Cross-Sectional Area
% crevasse 1
AT_data.cross_area.C1_11 = trapz(AT_data.along_track_crevasse.C1_11, -1.*AT_data.crevasse_pick.P11.crevasse1);
AT_data.cross_area.C1_14 = trapz(AT_data.along_track_crevasse.C1_14, -1.*AT_data.crevasse_pick.P14.crevasse1);
AT_data.cross_area.C1_18 = trapz(AT_data.along_track_crevasse.C1_18, -1.*AT_data.crevasse_pick.P18.crevasse1);
% crevasse 2
AT_data.cross_area.C2_11 = trapz(AT_data.along_track_crevasse.C2_11, -1.*AT_data.crevasse_pick.P11.crevasse2);
AT_data.cross_area.C2_14 = trapz(AT_data.along_track_crevasse.C2_14, -1.*AT_data.crevasse_pick.P14.crevasse2);
AT_data.cross_area.C2_18 = trapz(AT_data.along_track_crevasse.C2_18, -1.*AT_data.crevasse_pick.P18.crevasse2);
% crevasse 3
AT_data.cross_area.C3_1_11 = trapz(AT_data.along_track_crevasse.C3_1_11, -1.*AT_data.crevasse_pick.P11.crevasse3_1);
AT_data.cross_area.C3_2_11 = trapz(AT_data.along_track_crevasse.C3_2_11, -1.*AT_data.crevasse_pick.P11.crevasse3_2);
AT_data.cross_area.C3_3_11 = trapz(AT_data.along_track_crevasse.C3_3_11, -1.*AT_data.crevasse_pick.P11.crevasse3_3);
AT_data.cross_area.C3_14 = trapz(AT_data.along_track_crevasse.C3_14, -1.*AT_data.crevasse_pick.P14.crevasse3);
AT_data.cross_area.C3_18 = trapz(AT_data.along_track_crevasse.C3_18, -1.*AT_data.crevasse_pick.P18.crevasse3);
% crevasse 4
AT_data.cross_area.C4_11 = trapz(AT_data.along_track_crevasse.C4_11, -1.*AT_data.crevasse_pick.P11.crevasse4);
AT_data.cross_area.C4_14 = trapz(AT_data.along_track_crevasse.C4_14, -1.*AT_data.crevasse_pick.P14.crevasse4);
AT_data.cross_area.C4_18 = trapz(AT_data.along_track_crevasse.C4_18, -1.*AT_data.crevasse_pick.P18.crevasse4);
% crevasse 5
AT_data.cross_area.C5_1_11 = trapz(AT_data.along_track_crevasse.C5_1_11, -1.*AT_data.crevasse_pick.P11.crevasse5_1);
AT_data.cross_area.C5_2_11 = trapz(AT_data.along_track_crevasse.C5_2_11, -1.*AT_data.crevasse_pick.P11.crevasse5_2);
AT_data.cross_area.C5_3_11 = trapz(AT_data.along_track_crevasse.C5_3_11, -1.*AT_data.crevasse_pick.P11.crevasse5_3);
AT_data.cross_area.C5_1_14 = trapz(AT_data.along_track_crevasse.C5_1_14, -1.*AT_data.crevasse_pick.P14.crevasse5_1);
AT_data.cross_area.C5_2_14 = trapz(AT_data.along_track_crevasse.C5_2_14, -1.*AT_data.crevasse_pick.P14.crevasse5_2);
AT_data.cross_area.C5_18 = trapz(AT_data.along_track_crevasse.C5_18, -1.*AT_data.crevasse_pick.P18.crevasse5);
% crevasse 6
AT_data.cross_area.C6_1_11 = trapz(AT_data.along_track_crevasse.C6_1_11, -1.*AT_data.crevasse_pick.P11.crevasse6_1);
AT_data.cross_area.C6_2_11 = trapz(AT_data.along_track_crevasse.C6_2_11, -1.*AT_data.crevasse_pick.P11.crevasse6_2);
AT_data.cross_area.C6_1_14 = trapz(AT_data.along_track_crevasse.C6_1_14, -1.*AT_data.crevasse_pick.P14.crevasse6_1);
AT_data.cross_area.C6_2_14 = trapz(AT_data.along_track_crevasse.C6_2_14, -1.*AT_data.crevasse_pick.P14.crevasse6_2);
AT_data.cross_area.C6_18 = trapz(AT_data.along_track_crevasse.C6_18, -1.*AT_data.crevasse_pick.P18.crevasse6);
% crevasse 7
AT_data.cross_area.C7_11 = trapz(AT_data.along_track_crevasse.C7_11, -1.*AT_data.crevasse_pick.P11.crevasse7);
AT_data.cross_area.C7_14 = trapz(AT_data.along_track_crevasse.C7_14, -1.*AT_data.crevasse_pick.P14.crevasse7);
AT_data.cross_area.C7_18 = trapz(AT_data.along_track_crevasse.C7_18, -1.*AT_data.crevasse_pick.P18.crevasse7);
% crevasse 8
AT_data.cross_area.C8_11 = trapz(AT_data.along_track_crevasse.C8_11, -1.*AT_data.crevasse_pick.P11.crevasse8);
AT_data.cross_area.C8_14 = trapz(AT_data.along_track_crevasse.C8_14, -1.*AT_data.crevasse_pick.P14.crevasse8);
AT_data.cross_area.C8_18 = trapz(AT_data.along_track_crevasse.C8_18, -1.*AT_data.crevasse_pick.P18.crevasse8);
% crevasse 9
AT_data.cross_area.C9_11 = trapz(AT_data.along_track_crevasse.C9_11, -1.*AT_data.crevasse_pick.P11.crevasse9);
AT_data.cross_area.C9_14 = trapz(AT_data.along_track_crevasse.C9_14, -1.*AT_data.crevasse_pick.P14.crevasse9);
%AT_data.cross_area.C8_18 = trapz(AT_data.along_track_crevasse.C8_18, -1.*AT_data.crevasse_pick.P18.crevasse8);

%% Widths and heights of each crevasse (in meters)
% Crevasse 1 - Width
AT_data.crevasse_width.C1_11 = (max(AT_data.along_track_crevasse.C1_11) - min(AT_data.along_track_crevasse.C1_11))*1000;
AT_data.crevasse_width.C1_14 = (max(AT_data.along_track_crevasse.C1_14) - min(AT_data.along_track_crevasse.C1_14))*1000;
AT_data.crevasse_width.C1_18 = (max(AT_data.along_track_crevasse.C1_18) - min(AT_data.along_track_crevasse.C1_18))*1000;
% Crevasse 2 - Width
AT_data.crevasse_width.C2_11 = (max(AT_data.along_track_crevasse.C2_11) - min(AT_data.along_track_crevasse.C2_11))*1000;
AT_data.crevasse_width.C2_14 = (max(AT_data.along_track_crevasse.C2_14) - min(AT_data.along_track_crevasse.C2_14))*1000;
AT_data.crevasse_width.C2_18 = (max(AT_data.along_track_crevasse.C2_18) - min(AT_data.along_track_crevasse.C2_18))*1000;
% Crevasse 3 - Width
AT_data.crevasse_width.C3_1_11 = (max(AT_data.along_track_crevasse.C3_1_11) - min(AT_data.along_track_crevasse.C3_1_11))*1000;
AT_data.crevasse_width.C3_2_11 = (max(AT_data.along_track_crevasse.C3_2_11) - min(AT_data.along_track_crevasse.C3_2_11))*1000;
AT_data.crevasse_width.C3_3_11 = (max(AT_data.along_track_crevasse.C3_3_11) - min(AT_data.along_track_crevasse.C3_3_11))*1000;
AT_data.crevasse_width.C3_14 = (max(AT_data.along_track_crevasse.C3_14) - min(AT_data.along_track_crevasse.C3_14))*1000;
AT_data.crevasse_width.C3_18 = (max(AT_data.along_track_crevasse.C3_18) - min(AT_data.along_track_crevasse.C3_18))*1000;
% Crevasse 4 - Width
AT_data.crevasse_width.C4_11 = (max(AT_data.along_track_crevasse.C4_11) - min(AT_data.along_track_crevasse.C4_11))*1000;
AT_data.crevasse_width.C4_14 = (max(AT_data.along_track_crevasse.C4_14) - min(AT_data.along_track_crevasse.C4_14))*1000;
AT_data.crevasse_width.C4_18 = (max(AT_data.along_track_crevasse.C4_18) - min(AT_data.along_track_crevasse.C4_18))*1000;
% Crevasse 5 - Width
AT_data.crevasse_width.C5_1_11 = (max(AT_data.along_track_crevasse.C5_1_11) - min(AT_data.along_track_crevasse.C5_1_11))*1000;
AT_data.crevasse_width.C5_2_11 = (max(AT_data.along_track_crevasse.C5_2_11) - min(AT_data.along_track_crevasse.C5_2_11))*1000;
AT_data.crevasse_width.C5_3_11 = (max(AT_data.along_track_crevasse.C5_3_11) - min(AT_data.along_track_crevasse.C5_3_11))*1000;
AT_data.crevasse_width.C5_1_14 = (max(AT_data.along_track_crevasse.C5_1_14) - min(AT_data.along_track_crevasse.C5_1_14))*1000;
AT_data.crevasse_width.C5_2_14 = (max(AT_data.along_track_crevasse.C5_2_14) - min(AT_data.along_track_crevasse.C5_2_14))*1000;
AT_data.crevasse_width.C5_18 = (max(AT_data.along_track_crevasse.C5_18) - min(AT_data.along_track_crevasse.C3_18))*1000;
% Crevasse 6 - Width
AT_data.crevasse_width.C6_1_11 = (max(AT_data.along_track_crevasse.C6_1_11) - min(AT_data.along_track_crevasse.C6_1_11))*1000;
AT_data.crevasse_width.C6_2_11 = (max(AT_data.along_track_crevasse.C6_2_11) - min(AT_data.along_track_crevasse.C6_2_11))*1000;
AT_data.crevasse_width.C6_1_14 = (max(AT_data.along_track_crevasse.C6_1_14) - min(AT_data.along_track_crevasse.C6_1_14))*1000;
AT_data.crevasse_width.C6_2_14 = (max(AT_data.along_track_crevasse.C6_2_14) - min(AT_data.along_track_crevasse.C6_2_14))*1000;
AT_data.crevasse_width.C6_18 = (max(AT_data.along_track_crevasse.C6_18) - min(AT_data.along_track_crevasse.C6_18))*1000;
% Crevasse 7 - Width
AT_data.crevasse_width.C7_11 = (max(AT_data.along_track_crevasse.C7_11) - min(AT_data.along_track_crevasse.C7_11))*1000;
AT_data.crevasse_width.C7_14 = (max(AT_data.along_track_crevasse.C7_14) - min(AT_data.along_track_crevasse.C7_14))*1000;
AT_data.crevasse_width.C7_18 = (max(AT_data.along_track_crevasse.C7_18) - min(AT_data.along_track_crevasse.C7_18))*1000;
% Crevasse 8 - Width
AT_data.crevasse_width.C8_11 = (max(AT_data.along_track_crevasse.C8_11) - min(AT_data.along_track_crevasse.C8_11))*1000;
AT_data.crevasse_width.C8_14 = (max(AT_data.along_track_crevasse.C8_14) - min(AT_data.along_track_crevasse.C8_14))*1000;
AT_data.crevasse_width.C8_18 = (max(AT_data.along_track_crevasse.C8_18) - min(AT_data.along_track_crevasse.C8_18))*1000;
% Crevasse 9 - Width
AT_data.crevasse_width.C9_11 = (max(AT_data.along_track_crevasse.C9_11) - min(AT_data.along_track_crevasse.C9_11))*1000;
AT_data.crevasse_width.C9_14 = (max(AT_data.along_track_crevasse.C9_14) - min(AT_data.along_track_crevasse.C9_14))*1000;
%AT_data.crevasse_width.C9_18 = (max(AT_data.along_track_crevasse.C9_18) - min(AT_data.along_track_crevasse.C9_18))*1000;

% Crevasse 1 - Height
AT_data.crevasse_height.C1_11 = (max(AT_data.crevasse_pick.P11.crevasse1) - min(AT_data.crevasse_pick.P11.crevasse1));
AT_data.crevasse_height.C1_14 = (max(AT_data.crevasse_pick.P14.crevasse1) - min(AT_data.crevasse_pick.P14.crevasse1));
AT_data.crevasse_height.C1_18 = (max(AT_data.crevasse_pick.P18.crevasse1) - min(AT_data.crevasse_pick.P18.crevasse1));
% Crevasse 2 - Height
AT_data.crevasse_height.C2_11 = (max(AT_data.crevasse_pick.P11.crevasse2) - min(AT_data.crevasse_pick.P11.crevasse2));
AT_data.crevasse_height.C2_14 = (max(AT_data.crevasse_pick.P14.crevasse2) - min(AT_data.crevasse_pick.P14.crevasse2));
AT_data.crevasse_height.C2_18 = (max(AT_data.crevasse_pick.P18.crevasse2) - min(AT_data.crevasse_pick.P18.crevasse2));
% Crevasse 3 - Height
AT_data.crevasse_height.C3_1_11 = (max(AT_data.crevasse_pick.P11.crevasse3_1) - min(AT_data.crevasse_pick.P11.crevasse3_1));
AT_data.crevasse_height.C3_2_11 = (max(AT_data.crevasse_pick.P11.crevasse3_2) - min(AT_data.crevasse_pick.P11.crevasse3_2));
AT_data.crevasse_height.C3_3_11 = (max(AT_data.crevasse_pick.P11.crevasse3_3) - min(AT_data.crevasse_pick.P11.crevasse3_3));
AT_data.crevasse_height.C3_14 = (max(AT_data.crevasse_pick.P14.crevasse3) - min(AT_data.crevasse_pick.P14.crevasse3));
AT_data.crevasse_height.C3_18 = (max(AT_data.crevasse_pick.P18.crevasse3) - min(AT_data.crevasse_pick.P18.crevasse3));
% Crevasse 4 - Height
AT_data.crevasse_height.C4_11 = (max(AT_data.crevasse_pick.P11.crevasse4) - min(AT_data.crevasse_pick.P11.crevasse4));
AT_data.crevasse_height.C4_14 = (max(AT_data.crevasse_pick.P14.crevasse4) - min(AT_data.crevasse_pick.P14.crevasse4));
AT_data.crevasse_height.C4_18 = (max(AT_data.crevasse_pick.P18.crevasse4) - min(AT_data.crevasse_pick.P18.crevasse4));
% Crevasse 5 - Height
AT_data.crevasse_height.C5_1_11 = (max(AT_data.crevasse_pick.P11.crevasse5_1) - min(AT_data.crevasse_pick.P11.crevasse5_1));
AT_data.crevasse_height.C5_2_11 = (max(AT_data.crevasse_pick.P11.crevasse5_2) - min(AT_data.crevasse_pick.P11.crevasse5_2));
AT_data.crevasse_height.C5_3_11 = (max(AT_data.crevasse_pick.P11.crevasse5_3) - min(AT_data.crevasse_pick.P11.crevasse5_3));
AT_data.crevasse_height.C5_1_14 = (max(AT_data.crevasse_pick.P14.crevasse5_1) - min(AT_data.crevasse_pick.P14.crevasse5_1));
AT_data.crevasse_height.C5_2_14 = (max(AT_data.crevasse_pick.P14.crevasse5_2) - min(AT_data.crevasse_pick.P14.crevasse5_2));
AT_data.crevasse_height.C5_18 = (max(AT_data.crevasse_pick.P18.crevasse5) - min(AT_data.crevasse_pick.P18.crevasse5));
% Crevasse 6 - Height
AT_data.crevasse_height.C6_1_11 = (max(AT_data.crevasse_pick.P11.crevasse6_1) - min(AT_data.crevasse_pick.P11.crevasse6_1));
AT_data.crevasse_height.C6_2_11 = (max(AT_data.crevasse_pick.P11.crevasse6_2) - min(AT_data.crevasse_pick.P11.crevasse6_2));
AT_data.crevasse_height.C6_1_14 = (max(AT_data.crevasse_pick.P14.crevasse6_1) - min(AT_data.crevasse_pick.P14.crevasse6_1));
AT_data.crevasse_height.C6_2_14 = (max(AT_data.crevasse_pick.P14.crevasse6_2) - min(AT_data.crevasse_pick.P14.crevasse6_2));
AT_data.crevasse_height.C6_18 = (max(AT_data.crevasse_pick.P18.crevasse6) - min(AT_data.crevasse_pick.P18.crevasse6));
% Crevasse 7 - Height
AT_data.crevasse_height.C7_11 = (max(AT_data.crevasse_pick.P11.crevasse7) - min(AT_data.crevasse_pick.P11.crevasse7));
AT_data.crevasse_height.C7_14 = (max(AT_data.crevasse_pick.P14.crevasse7) - min(AT_data.crevasse_pick.P14.crevasse7));
AT_data.crevasse_height.C7_18 = (max(AT_data.crevasse_pick.P18.crevasse7) - min(AT_data.crevasse_pick.P18.crevasse7));
% Crevasse 8 - Height
AT_data.crevasse_height.C8_11 = (max(AT_data.crevasse_pick.P11.crevasse8) - min(AT_data.crevasse_pick.P11.crevasse8));
AT_data.crevasse_height.C8_14 = (max(AT_data.crevasse_pick.P14.crevasse8) - min(AT_data.crevasse_pick.P14.crevasse8));
AT_data.crevasse_height.C8_18 = (max(AT_data.crevasse_pick.P18.crevasse8) - min(AT_data.crevasse_pick.P18.crevasse8));
% Crevasse 9 - Height
AT_data.crevasse_height.C9_11 = (max(AT_data.crevasse_pick.P11.crevasse9) - min(AT_data.crevasse_pick.P11.crevasse9));
AT_data.crevasse_height.C9_14 = (max(AT_data.crevasse_pick.P14.crevasse9) - min(AT_data.crevasse_pick.P14.crevasse9));
%AT_data.crevasse_height.C9_18 = (max(AT_data.crevasse_pick.P18.crevasse9) - min(AT_data.crevasse_pick.P18.crevasse9));

%% Total Lists 
% Area lists
AT_data.all_season_area.C1 = [AT_data.cross_area.C1_11, AT_data.cross_area.C1_14, AT_data.cross_area.C1_18];
AT_data.all_season_area.C2 = [AT_data.cross_area.C2_11, AT_data.cross_area.C2_14, AT_data.cross_area.C2_18];
AT_data.all_season_area.C9 = [AT_data.cross_area.C3_1_11, AT_data.cross_area.C3_2_11, AT_data.cross_area.C3_3_11,...
    AT_data.cross_area.C3_14, AT_data.cross_area.C3_18];
AT_data.all_season_area.C4 = [AT_data.cross_area.C4_11, AT_data.cross_area.C4_14, AT_data.cross_area.C4_18];
AT_data.all_season_area.C5 = [AT_data.cross_area.C5_1_11, AT_data.cross_area.C5_2_11, AT_data.cross_area.C5_3_11,...
    AT_data.cross_area.C5_1_14, AT_data.cross_area.C5_2_14, AT_data.cross_area.C5_18];
AT_data.all_season_area.C6 = [AT_data.cross_area.C6_1_11, AT_data.cross_area.C6_2_11, AT_data.cross_area.C6_1_14, ...
    AT_data.cross_area.C6_2_14, AT_data.cross_area.C6_18];
AT_data.all_season_area.C7 = [AT_data.cross_area.C7_11, AT_data.cross_area.C7_14, AT_data.cross_area.C7_18];
AT_data.all_season_area.C8 = [AT_data.cross_area.C8_11, AT_data.cross_area.C8_14, AT_data.cross_area.C8_18];
AT_data.all_season_area.C9 = [AT_data.cross_area.C9_11, AT_data.cross_area.C9_14];

% Year areas
AT_data.years.y11 = 2011;
AT_data.years.y14 = 2014;
AT_data.years.y18 = 2018;

AT_data.all_dates.array.C1 = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C2 = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C3 = [AT_data.years.y11, AT_data.years.y11, AT_data.years.y11, ...
    AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C3_sub = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C4 = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C5 = [AT_data.years.y11, AT_data.years.y11, AT_data.years.y11, ...
    AT_data.years.y14, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C5_sub = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C6 = [AT_data.years.y11, AT_data.years.y11, AT_data.years.y14, ...
    AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C6_sub = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C7 = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C8 = [AT_data.years.y11, AT_data.years.y14, AT_data.years.y18];
AT_data.all_dates.array.C9 = [AT_data.years.y11, AT_data.years.y14];

% Total widths for each feature
AT_data.all_season_widths.C1 = [AT_data.crevasse_width.C1_11, AT_data.crevasse_width.C1_14, AT_data.crevasse_width.C1_18];
AT_data.all_season_widths.C2 = [AT_data.crevasse_width.C2_11, AT_data.crevasse_width.C2_14, AT_data.crevasse_width.C2_18];
AT_data.all_season_widths.C3 = [AT_data.crevasse_width.C3_1_11, AT_data.crevasse_width.C3_2_11, AT_data.crevasse_width.C3_3_11,...
    AT_data.crevasse_width.C1_14, AT_data.crevasse_width.C1_18];
AT_data.all_season_widths.C4 = [AT_data.crevasse_width.C4_11, AT_data.crevasse_width.C4_14, AT_data.crevasse_width.C4_18];
AT_data.all_season_widths.C5 = [AT_data.crevasse_width.C5_1_11, AT_data.crevasse_width.C5_2_11, AT_data.crevasse_width.C5_3_11, ...
    AT_data.crevasse_width.C5_1_14, AT_data.crevasse_width.C5_2_14, AT_data.crevasse_width.C1_18];
AT_data.all_season_widths.C6 = [AT_data.crevasse_width.C6_1_11, AT_data.crevasse_width.C6_2_11, AT_data.crevasse_width.C6_1_14, ...
    AT_data.crevasse_width.C6_2_14, AT_data.crevasse_width.C6_18];
AT_data.all_season_widths.C7 = [AT_data.crevasse_width.C7_11, AT_data.crevasse_width.C7_14, AT_data.crevasse_width.C7_18];
AT_data.all_season_widths.C8 = [AT_data.crevasse_width.C8_11, AT_data.crevasse_width.C8_14, AT_data.crevasse_width.C8_18];
AT_data.all_season_widths.C9 = [AT_data.crevasse_width.C9_11, AT_data.crevasse_width.C9_14];

% Total height differences for each feature annually
AT_data.all_season_heights.C1 = [AT_data.crevasse_height.C1_11, AT_data.crevasse_height.C1_14, AT_data.crevasse_height.C1_18];
AT_data.all_season_heights.C2 = [AT_data.crevasse_height.C2_11, AT_data.crevasse_height.C2_14, AT_data.crevasse_height.C2_18];
AT_data.all_season_heights.C3 = [AT_data.crevasse_height.C3_1_11, AT_data.crevasse_height.C3_2_11, AT_data.crevasse_height.C3_3_11,...
    AT_data.crevasse_height.C3_14, AT_data.crevasse_height.C3_18];
AT_data.all_season_heights.C4 = [AT_data.crevasse_height.C4_11, AT_data.crevasse_height.C4_14, AT_data.crevasse_height.C4_18];
AT_data.all_season_heights.C5 = [AT_data.crevasse_height.C5_1_11, AT_data.crevasse_height.C5_2_11, AT_data.crevasse_height.C5_3_11,...
    AT_data.crevasse_height.C5_1_14, AT_data.crevasse_height.C5_2_14, AT_data.crevasse_height.C5_18];
AT_data.all_season_heights.C6 = [AT_data.crevasse_height.C6_1_11, AT_data.crevasse_height.C6_2_11, AT_data.crevasse_height.C6_1_14, ...
    AT_data.crevasse_height.C6_1_14, AT_data.crevasse_height.C6_18];
AT_data.all_season_heights.C7 = [AT_data.crevasse_height.C7_11, AT_data.crevasse_height.C7_14, AT_data.crevasse_height.C7_18];
AT_data.all_season_heights.C8 = [AT_data.crevasse_height.C8_11, AT_data.crevasse_height.C8_14, AT_data.crevasse_height.C8_18];
AT_data.all_season_heights.C9 = [AT_data.crevasse_height.C9_11, AT_data.crevasse_height.C9_14];

% Total apex heights for each feature annually
AT_data.season_apex.C1 = [max(AT_data.crevasse_pick.P11.crevasse1), max(AT_data.crevasse_pick.P14.crevasse1),...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C2 = [max(AT_data.crevasse_pick.P11.crevasse2), max(AT_data.crevasse_pick.P14.crevasse2),...
    max(AT_data.crevasse_pick.P18.crevasse2)];
AT_data.season_apex.C3 = [max(AT_data.crevasse_pick.P11.crevasse3_1), max(AT_data.crevasse_pick.P11.crevasse3_2), ...
    max(AT_data.crevasse_pick.P11.crevasse3_3), max(AT_data.crevasse_pick.P14.crevasse1),max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C3_1 = [max(AT_data.crevasse_pick.P11.crevasse3_1), max(AT_data.crevasse_pick.P14.crevasse1), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C3_2 = [max(AT_data.crevasse_pick.P11.crevasse3_2), max(AT_data.crevasse_pick.P14.crevasse1), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C3_3 = [max(AT_data.crevasse_pick.P11.crevasse3_3), max(AT_data.crevasse_pick.P14.crevasse1), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C4 = [max(AT_data.crevasse_pick.P11.crevasse4), max(AT_data.crevasse_pick.P14.crevasse4),...
    max(AT_data.crevasse_pick.P18.crevasse4)];
AT_data.season_apex.C5 = [max(AT_data.crevasse_pick.P11.crevasse5_1), max(AT_data.crevasse_pick.P11.crevasse5_2),...
    max(AT_data.crevasse_pick.P11.crevasse5_3), max(AT_data.crevasse_pick.P14.crevasse5_1), ...
    max(AT_data.crevasse_pick.P14.crevasse5_2), max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C5_1 = [max(AT_data.crevasse_pick.P11.crevasse5_1), max(AT_data.crevasse_pick.P14.crevasse5_1), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C5_2 = [max(AT_data.crevasse_pick.P11.crevasse5_2), max(AT_data.crevasse_pick.P14.crevasse5_2), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C5_3 = [max(AT_data.crevasse_pick.P11.crevasse5_3), max(AT_data.crevasse_pick.P14.crevasse5_2), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C6 = [max(AT_data.crevasse_pick.P11.crevasse6_1), max(AT_data.crevasse_pick.P11.crevasse6_2),...
    max(AT_data.crevasse_pick.P14.crevasse6_1), max(AT_data.crevasse_pick.P14.crevasse6_2), max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C6_1 = [max(AT_data.crevasse_pick.P11.crevasse6_1), max(AT_data.crevasse_pick.P14.crevasse6_1), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C6_2 = [max(AT_data.crevasse_pick.P11.crevasse6_2), max(AT_data.crevasse_pick.P14.crevasse6_2), ...
    max(AT_data.crevasse_pick.P18.crevasse1)];
AT_data.season_apex.C7 = [max(AT_data.crevasse_pick.P11.crevasse7), max(AT_data.crevasse_pick.P14.crevasse7),...
    max(AT_data.crevasse_pick.P18.crevasse7)];
AT_data.season_apex.C8 = [max(AT_data.crevasse_pick.P11.crevasse8), max(AT_data.crevasse_pick.P14.crevasse8),...
    max(AT_data.crevasse_pick.P18.crevasse8)];
AT_data.season_apex.C9 = [max(AT_data.crevasse_pick.P11.crevasse9), max(AT_data.crevasse_pick.P14.crevasse9)];
%%
figure(45)
h1 = plot(AT_data.all_dates.array.C1, AT_data.all_season_widths.C1, '-o');
hold on
h2 = plot(AT_data.all_dates.array.C2, AT_data.all_season_widths.C2, '-o');
h3 = plot(AT_data.all_dates.array.C3, AT_data.all_season_widths.C3, '-o');
h5 = plot(AT_data.all_dates.array.C4, AT_data.all_season_widths.C4, '-o');
h6 = plot(AT_data.all_dates.array.C5, AT_data.all_season_widths.C5, '-o');
h7 = plot(AT_data.all_dates.array.C6, AT_data.all_season_widths.C6, '-o');
h8 = plot(AT_data.all_dates.array.C7, AT_data.all_season_widths.C7, '-o');
h9 = plot(AT_data.all_dates.array.C8, AT_data.all_season_widths.C8, '-o');
h10 = plot(AT_data.all_dates.array.C9, AT_data.all_season_widths.C9, '-o');
title('Annual Change in Crevasse Width');
xlabel('Survey Year (yr)');
ylabel('crevasse width (m)');
xlim([2010 2020]);
ylim([0 500]);
hold off

figure(46)
h1 = plot(AT_data.all_dates.array.C1, AT_data.all_season_heights.C1, '-o');
hold on
h2 = plot(AT_data.all_dates.array.C2, AT_data.all_season_heights.C2, '-o');
h3 = plot(AT_data.all_dates.array.C3, AT_data.all_season_heights.C3, '-o');
h5 = plot(AT_data.all_dates.array.C4, AT_data.all_season_heights.C4, '-o');
h6 = plot(AT_data.all_dates.array.C5, AT_data.all_season_heights.C5, '-o');
h7 = plot(AT_data.all_dates.array.C6, AT_data.all_season_heights.C6, '-o');
h8 = plot(AT_data.all_dates.array.C7, AT_data.all_season_heights.C7, '-o');
h9 = plot(AT_data.all_dates.array.C8, AT_data.all_season_heights.C8, '-o');
h10 = plot(AT_data.all_dates.array.C9, AT_data.all_season_heights.C9, '-o');
title('Annual Change in Crevasse Height');
xlabel('Survey Year (yr)');
ylabel('crevasse height between base and apex (m)');
xlim([2010 2020]);
%ylim([0 500]);
hold off

figure(47)
h1 = plot(AT_data.all_dates.array.C1, AT_data.season_apex.C1, '-o');
hold on
h2 = plot(AT_data.all_dates.array.C2, AT_data.season_apex.C2, '-o');
h3 = plot(AT_data.all_dates.array.C3, AT_data.season_apex.C3, '-o');
h5 = plot(AT_data.all_dates.array.C4, AT_data.season_apex.C4, '-o');
h6 = plot(AT_data.all_dates.array.C5, AT_data.season_apex.C5, '-o');
h7 = plot(AT_data.all_dates.array.C6, AT_data.season_apex.C6, '-o');
h8 = plot(AT_data.all_dates.array.C7, AT_data.season_apex.C7, '-o');
h9 = plot(AT_data.all_dates.array.C8, AT_data.season_apex.C8, '-o');
h10 = plot(AT_data.all_dates.array.C9, AT_data.season_apex.C9, '-o');
h11 = plot(AT_data.all_dates.array.C1, AT_data.apex_melt.C1, '-d');
h12 = plot(AT_data.all_dates.array.C2, AT_data.apex_melt.C2, '-d');
h13 = plot(AT_data.all_dates.array.C3, AT_data.apex_melt.C3, '-d');
h14 = plot(AT_data.all_dates.array.C4, AT_data.apex_melt.C4, '-d');
h15 = plot(AT_data.all_dates.array.C5, AT_data.apex_melt.C5, '-d');
h16 = plot(AT_data.all_dates.array.C6, AT_data.apex_melt.C6, '-d');
h17 = plot(AT_data.all_dates.array.C7, AT_data.apex_melt.C7, '-d');
h18 = plot(AT_data.all_dates.array.C8, AT_data.apex_melt.C8, '-d');
h19 = plot(AT_data.all_dates.array.C9, AT_data.apex_melt.C9, '-d');
title('Annual Change in Crevasse Apex elevation');
xlabel('Survey Year (yr)');
ylabel('crevasse apex elevation (m)');
xlim([2010 2020]);
%ylim([0 500]);
hold off


figure(48)
h1 = plot(AT_data.all_dates.array.C1, AT_data.season_apex.C1, '-o','Color','r');
hold on
h2 = plot(AT_data.all_dates.array.C2, AT_data.season_apex.C2, '-o','Color','#A2142F');
h3 = plot(AT_data.all_dates.array.C3_sub, AT_data.season_apex.C3_1, '-o','Color','#D95319');
h4 = plot(AT_data.all_dates.array.C4, AT_data.season_apex.C4, '-o','Color','#EDB120');
h5 = plot(AT_data.all_dates.array.C5_sub, AT_data.season_apex.C5_1, '-o','Color','g');
h8 = plot(AT_data.all_dates.array.C6_sub, AT_data.season_apex.C6_1, '-o','Color','#77AC30');
h7 = plot(AT_data.all_dates.array.C7, AT_data.season_apex.C7, '-o','Color','b');
h8 = plot(AT_data.all_dates.array.C8, AT_data.season_apex.C8, '-o','Color',	'#0072BD');
h9 = plot(AT_data.all_dates.array.C9, AT_data.season_apex.C9, '-o','Color','#7E2F8E');
h4 = plot(AT_data.all_dates.array.C3_sub, AT_data.season_apex.C3_2, '-o','Color','#D95319'); 
h5 = plot(AT_data.all_dates.array.C3_sub, AT_data.season_apex.C3_3, '-o','Color','#D95319');
h6 = plot(AT_data.all_dates.array.C5_sub, AT_data.season_apex.C5_2, '-o','Color','g');
h7 = plot(AT_data.all_dates.array.C5_sub, AT_data.season_apex.C5_3, '-o','Color','g');
h8 = plot(AT_data.all_dates.array.C6_sub, AT_data.season_apex.C6_2, '-o','Color','#77AC30');
title('Annual Change in Crevasse Apex elevation');
xlabel('Survey Year (yr)');
ylabel('crevasse apex elevation (m)');
legend('Crevasse 1', 'Crevasse 2', 'Crevasse 3', 'Crevasse 4', 'Crevasse 5', ...
    'Crevasse 6', 'Crevasse 7', 'Crevasse 8', 'Crevasse 9', 'Location','best');
xlim([2010 2020]);
%ylim([0 500]);
hold off
%% Basal Melt rate Calculation using Lebrocq et al, 2013 forumula  
% Divided into 3 sections, 2011-2014, 2014-2018, and 2011-2018
% pull thickness arrays from beginning and plug into formula below.
% accumulation is believed to be 0, thus accum variable = 0, otherwise
% specificy a value
accum = 0; % switch to value if known
% thickness profiles and Velocties
H1 = AT_data.ice_thickness_interp.P11;
H2 = AT_data.ice_thickness_interp.P14;
H3 = AT_data.ice_thickness_interp.P18;
VEL1 = AT_data.interp_data.p11_vel_corrected;
VEL2 = AT_data.interp_data.p14_vel_corrected;
VEL3 = AT_data.interp_data.p18_vel_corrected;
dist1 = AT_data.interp_data.P11_DIST;
dist2 = AT_data.interp_data.P14_DIST;
dist3 = AT_data.interp_data.P18_DIST;

% Annual Melt rates
tmp_melt1 = accum - ((H2.*VEL2)-(H1.*VEL1));
tmp_melt2 = accum - ((H3.*VEL3)-(H2.*VEL2));
tmp_melt3 = accum - ((H3.*VEL3)-(H1.*VEL1));
melt_rate1 = tmp_melt1./dist1; 
melt_rate2 = tmp_melt2./dist2;
melt_rate3 = tmp_melt3./dist3;

% total melt
total_melt1 = (2014-2011)*melt_rate1;
total_melt2 = (2018-2014)*melt_rate2;
total_melt3 = (2018-2011)*melt_rate3;

figure(1)
subplot(3,1,1)
h1 = plot(AT_data.query_array.P11./1000, (AT_data.interp_data.P14 - AT_data.interp_data.P11),'r');
hold on 
h2 = plot(AT_data.query_array.P11./1000, total_melt1,'k');
legend('Observed Elevation change','Calculated melt');
title('2011-2014 Observed elevation change vs Calculated melt');
xlabel('along track distance (km)');
ylabel('elevation change (m)')
ylim([-100 300]);
grid on;
hold off

subplot(3,1,2)
h1 = plot(AT_data.query_array.P11./1000, (AT_data.interp_data.P18 - AT_data.interp_data.P14),'g');
hold on 
h2 = plot(AT_data.query_array.P11./1000, total_melt2,'k');
legend('Observed Elevation change','Calculated melt');
title('2014-2018 Observed elevation change vs Calculated melt');
xlabel('along track distance (km)');
ylabel('elevation change (m)')
ylim([-100 300]);
grid on;
hold off

subplot(3,1,3)
h1 = plot(AT_data.query_array.P11./1000, (AT_data.interp_data.P18 - AT_data.interp_data.P11),'b');
hold on 
h2 = plot(AT_data.query_array.P11./1000, total_melt3,'k');
legend('Observed Elevation change','Calculated melt');
title('2014-2018 Observed elevation change vs Calculated melt');
xlabel('along track distance (km)');
ylabel('elevation change (m)')
ylim([-100 300]);
grid on;
hold off
%%
figure(8)
subplot(3,1,1)
h1 = plot(AT_data.crevasse_width.C1_11, AT_data.crevasse_height.C1_11,'*','MarkerSize',10);
hold on;
h2 = plot(AT_data.crevasse_width.C2_11, AT_data.crevasse_height.C2_11,'*','MarkerSize',10);
h3 = plot(AT_data.crevasse_width.C3_1_11, AT_data.crevasse_height.C3_1_11,'*','MarkerSize',10);
h4 = plot(AT_data.crevasse_width.C3_2_11, AT_data.crevasse_height.C3_2_11,'*','MarkerSize',10);
h5 = plot(AT_data.crevasse_width.C3_3_11, AT_data.crevasse_height.C3_3_11,'*','MarkerSize',10);
h6 = plot(AT_data.crevasse_width.C4_11, AT_data.crevasse_height.C4_11,'*','MarkerSize',10);
h7 = plot(AT_data.crevasse_width.C5_1_11, AT_data.crevasse_height.C5_1_11,'*','MarkerSize',10);
h8 = plot(AT_data.crevasse_width.C5_2_11, AT_data.crevasse_height.C5_2_11,'*','MarkerSize',10);
h9 = plot(AT_data.crevasse_width.C5_3_11, AT_data.crevasse_height.C5_3_11,'*','MarkerSize',10);
h10 = plot(AT_data.crevasse_width.C6_1_11, AT_data.crevasse_height.C6_1_11,'*','MarkerSize',10);
h11 = plot(AT_data.crevasse_width.C6_2_11, AT_data.crevasse_height.C6_2_11,'*','MarkerSize',10);
h12 = plot(AT_data.crevasse_width.C7_11, AT_data.crevasse_height.C7_11,'*','MarkerSize',10);
h13 = plot(AT_data.crevasse_width.C8_11, AT_data.crevasse_height.C8_11,'*','MarkerSize',10);
h14 = plot(AT_data.crevasse_width.C9_11, AT_data.crevasse_height.C9_11,'*','MarkerSize',10);
title('2011 Basal crevasse height vs. width');
xlabel('crevasse width (m)');
ylabel('crevasse heigh (m)');
%xlim([0 200]);
%ylim([0 0.7]);
%legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
%  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
%  'Location', 'southeast');
hold off

subplot(3,1,2)
h1 = plot(AT_data.crevasse_width.C1_14, AT_data.crevasse_height.C1_14,'*','MarkerSize',10);
hold on;
h2 = plot(AT_data.crevasse_width.C2_14, AT_data.crevasse_height.C2_14,'*','MarkerSize',10);
h3 = plot(AT_data.crevasse_width.C3_14, AT_data.crevasse_height.C3_14,'*','MarkerSize',10);
h4 = plot(AT_data.crevasse_width.C4_14, AT_data.crevasse_height.C4_14,'*','MarkerSize',10);
h5 = plot(AT_data.crevasse_width.C5_1_14, AT_data.crevasse_height.C5_1_14,'*','MarkerSize',10);
h6 = plot(AT_data.crevasse_width.C5_2_14, AT_data.crevasse_height.C5_2_14,'*','MarkerSize',10);
h7 = plot(AT_data.crevasse_width.C6_1_14, AT_data.crevasse_height.C6_1_14,'*','MarkerSize',10);
h8 = plot(AT_data.crevasse_width.C6_2_14, AT_data.crevasse_height.C6_2_14,'*','MarkerSize',10);
h9 = plot(AT_data.crevasse_width.C7_14, AT_data.crevasse_height.C7_14,'*','MarkerSize',10);
h10 = plot(AT_data.crevasse_width.C8_14, AT_data.crevasse_height.C8_14,'*','MarkerSize',10);
h11 = plot(AT_data.crevasse_width.C9_14, AT_data.crevasse_height.C9_14,'*','MarkerSize',10);
title('2014 Basal crevasse height vs. width');
xlabel('crevasse width (m)');
ylabel('crevasse heigh (m)');
%xlim([0 200]);
%ylim([0 0.7]);
%legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
%  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
%  'Location', 'southeast');
hold off

subplot(3,1,3)
h1 = plot(AT_data.crevasse_width.C1_18, AT_data.crevasse_height.C1_18,'*','MarkerSize',10);
hold on;
h2 = plot(AT_data.crevasse_width.C2_18, AT_data.crevasse_height.C2_18,'*','MarkerSize',10);
h3 = plot(AT_data.crevasse_width.C3_18, AT_data.crevasse_height.C3_18,'*','MarkerSize',10);
h4 = plot(AT_data.crevasse_width.C4_18, AT_data.crevasse_height.C4_18,'*','MarkerSize',10);
h5 = plot(AT_data.crevasse_width.C5_18, AT_data.crevasse_height.C5_18,'*','MarkerSize',10);
h6 = plot(AT_data.crevasse_width.C6_18, AT_data.crevasse_height.C6_18,'*','MarkerSize',10);
h7 = plot(AT_data.crevasse_width.C7_18, AT_data.crevasse_height.C7_18,'*','MarkerSize',10);
h8 = plot(AT_data.crevasse_width.C8_18, AT_data.crevasse_height.C8_18,'*','MarkerSize',10);
%h9 = plot(AT_data.crevasse_width.C9_18, AT_data.crevasse_height.C9_18,'*','MarkerSize',10);
title('2018 Basal crevasse height vs. width');
xlabel('crevasse width (m)');
ylabel('crevasse heigh (m)');
%xlim([0 200]);
%ylim([0 0.7]);
%legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
%  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
%  'Location', 'southeast');
hold off
%% Height to width ratio (H/W)
%Crevasse 1
AT_data.H2W.C1_11 = (AT_data.crevasse_height.C1_11)/(AT_data.crevasse_width.C1_11);
AT_data.H2W.C1_14 = (AT_data.crevasse_height.C1_14)/(AT_data.crevasse_width.C1_14);
AT_data.H2W.C1_18 = (AT_data.crevasse_height.C1_18)/(AT_data.crevasse_width.C1_18);
%Crevasse 2
AT_data.H2W.C2_11 = (AT_data.crevasse_height.C2_11)/(AT_data.crevasse_width.C2_11);
AT_data.H2W.C2_14 = (AT_data.crevasse_height.C2_14)/(AT_data.crevasse_width.C2_14);
AT_data.H2W.C2_18 = (AT_data.crevasse_height.C2_18)/(AT_data.crevasse_width.C2_18);
%Crevasse 3
AT_data.H2W.C3_1_11 = (AT_data.crevasse_height.C3_1_11)/(AT_data.crevasse_width.C3_1_11);
AT_data.H2W.C3_2_11 = (AT_data.crevasse_height.C3_2_11)/(AT_data.crevasse_width.C3_2_11);
AT_data.H2W.C3_3_11 = (AT_data.crevasse_height.C3_3_11)/(AT_data.crevasse_width.C3_3_11);
AT_data.H2W.C3_14 = (AT_data.crevasse_height.C3_14)/(AT_data.crevasse_width.C3_14);
AT_data.H2W.C3_18 = (AT_data.crevasse_height.C3_18)/(AT_data.crevasse_width.C3_18);
%Crevasse 4
AT_data.H2W.C4_11 = (AT_data.crevasse_height.C4_11)/(AT_data.crevasse_width.C4_11);
AT_data.H2W.C4_14 = (AT_data.crevasse_height.C4_14)/(AT_data.crevasse_width.C4_14);
AT_data.H2W.C4_18 = (AT_data.crevasse_height.C4_18)/(AT_data.crevasse_width.C4_18);
%Crevasse 5
AT_data.H2W.C5_1_11 = (AT_data.crevasse_height.C5_1_11)/(AT_data.crevasse_width.C5_1_11);
AT_data.H2W.C5_2_11 = (AT_data.crevasse_height.C5_2_11)/(AT_data.crevasse_width.C5_2_11);
AT_data.H2W.C5_3_11 = (AT_data.crevasse_height.C5_3_11)/(AT_data.crevasse_width.C5_3_11);
AT_data.H2W.C5_1_14 = (AT_data.crevasse_height.C5_1_14)/(AT_data.crevasse_width.C5_1_14);
AT_data.H2W.C5_2_14 = (AT_data.crevasse_height.C5_2_14)/(AT_data.crevasse_width.C5_2_14);
AT_data.H2W.C5_18 = (AT_data.crevasse_height.C5_18)/(AT_data.crevasse_width.C5_18);
%Crevasse 6
AT_data.H2W.C6_1_11 = (AT_data.crevasse_height.C6_1_11)/(AT_data.crevasse_width.C6_1_11);
AT_data.H2W.C6_2_11 = (AT_data.crevasse_height.C6_2_11)/(AT_data.crevasse_width.C6_2_11);
AT_data.H2W.C6_1_14 = (AT_data.crevasse_height.C6_1_14)/(AT_data.crevasse_width.C6_1_14);
AT_data.H2W.C6_2_14 = (AT_data.crevasse_height.C6_2_14)/(AT_data.crevasse_width.C6_2_14);
AT_data.H2W.C6_18 = (AT_data.crevasse_height.C6_18)/(AT_data.crevasse_width.C6_18);
%Crevasse 7
AT_data.H2W.C7_11 = (AT_data.crevasse_height.C7_11)/(AT_data.crevasse_width.C7_11);
AT_data.H2W.C7_14 = (AT_data.crevasse_height.C7_14)/(AT_data.crevasse_width.C7_14);
AT_data.H2W.C7_18 = (AT_data.crevasse_height.C7_18)/(AT_data.crevasse_width.C7_18);
%Crevasse 8
AT_data.H2W.C8_11 = (AT_data.crevasse_height.C8_11)/(AT_data.crevasse_width.C8_11);
AT_data.H2W.C8_14 = (AT_data.crevasse_height.C8_14)/(AT_data.crevasse_width.C8_14);
AT_data.H2W.C8_18 = (AT_data.crevasse_height.C8_18)/(AT_data.crevasse_width.C8_18);
%Crevasse 9
AT_data.H2W.C9_11 = (AT_data.crevasse_height.C9_11)/(AT_data.crevasse_width.C9_11);
AT_data.H2W.C9_14 = (AT_data.crevasse_height.C9_14)/(AT_data.crevasse_width.C9_14);
%AT_data.H2W.C9_18 = (AT_data.crevasse_height.C9_18)/(AT_data.crevasse_width.C9_18);

% Years
AT_data.years.p11 = 2011;
AT_data.years.p14 = 2014;
AT_data.years.p18 = 2018;


%% Cross sectional area vs H2W
figure(9)
subplot(3,1,1)
h1 = plot(AT_data.cross_area.C1_11, AT_data.H2W.C1_11,'*','MarkerSize',10);
hold on;
h2 = plot(AT_data.cross_area.C2_11, AT_data.H2W.C2_11,'*','MarkerSize',10);
h3 = plot(AT_data.cross_area.C3_1_11, AT_data.H2W.C3_1_11,'*','MarkerSize',10);
h4 = plot(AT_data.cross_area.C3_2_11, AT_data.H2W.C3_2_11,'*','MarkerSize',10);
h5 = plot(AT_data.cross_area.C3_3_11, AT_data.H2W.C3_3_11,'*','MarkerSize',10);
h6 = plot(AT_data.cross_area.C4_11, AT_data.H2W.C4_11,'*','MarkerSize',10);
h7 = plot(AT_data.cross_area.C5_1_11, AT_data.H2W.C5_1_11,'*','MarkerSize',10);
h8 = plot(AT_data.cross_area.C5_2_11, AT_data.H2W.C5_2_11,'*','MarkerSize',10);
h9 = plot(AT_data.cross_area.C5_3_11, AT_data.H2W.C5_3_11,'*','MarkerSize',10);
h10 = plot(AT_data.cross_area.C6_1_11, AT_data.H2W.C6_1_11,'*','MarkerSize',10);
h11 = plot(AT_data.cross_area.C6_2_11, AT_data.H2W.C6_2_11,'*','MarkerSize',10);
h12 = plot(AT_data.cross_area.C7_11, AT_data.H2W.C7_11,'*','MarkerSize',10);
h13 = plot(AT_data.cross_area.C8_11, AT_data.H2W.C8_11,'*','MarkerSize',10);
h14 = plot(AT_data.cross_area.C9_11, AT_data.H2W.C9_11,'*','MarkerSize',10);
title('2011 Basal crevasse cross-sectional area vs. height-width value');
xlabel('cross sectional area (m^2)');
ylabel('Height-Width Ratio');
xlim([0 200]);
ylim([0 0.7]);
%legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
%  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
%  'Location', 'southeast');
hold off

subplot(3,1,2)
h1 = plot(AT_data.cross_area.C1_14, AT_data.H2W.C1_14,'*','MarkerSize',10);
text(84, 0.572,'C1');
hold on;
h2 = plot(AT_data.cross_area.C2_14, AT_data.H2W.C2_14,'*','MarkerSize',10);
h3 = plot(AT_data.cross_area.C3_14, AT_data.H2W.C3_14,'*','MarkerSize',10);
h4 = plot(AT_data.cross_area.C4_14, AT_data.H2W.C4_14,'*','MarkerSize',10);
h5 = plot(AT_data.cross_area.C5_1_14, AT_data.H2W.C5_1_14,'*','MarkerSize',10);
h6 = plot(AT_data.cross_area.C5_2_14, AT_data.H2W.C5_2_14,'*','MarkerSize',10);
h7 = plot(AT_data.cross_area.C6_1_14, AT_data.H2W.C6_1_14,'*','MarkerSize',10);
h8 = plot(AT_data.cross_area.C6_2_14, AT_data.H2W.C6_2_14,'*','MarkerSize',10);
h9 = plot(AT_data.cross_area.C7_14, AT_data.H2W.C7_14,'*','MarkerSize',10);
h10 = plot(AT_data.cross_area.C8_14, AT_data.H2W.C8_14,'*','MarkerSize',10);
h11 = plot(AT_data.cross_area.C9_14, AT_data.H2W.C9_14,'*','MarkerSize',10);
title('2014 Basal crevasse cross-sectional area vs. height-width value');
xlabel('cross sectional area (m^2)');
ylabel('Height-Width Ratio');
xlim([0 200]);
ylim([0 0.7]);
%legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
%  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
%  'Location', 'southeast');
hold off

subplot(3,1,3)
h1 = plot(AT_data.cross_area.C1_18, AT_data.H2W.C1_18,'*','MarkerSize',10);
hold on;
h2 = plot(AT_data.cross_area.C2_18, AT_data.H2W.C2_18,'*','MarkerSize',10);
h3 = plot(AT_data.cross_area.C3_18, AT_data.H2W.C3_18,'*','MarkerSize',10);
h4 = plot(AT_data.cross_area.C4_18, AT_data.H2W.C4_18,'*','MarkerSize',10);
h5 = plot(AT_data.cross_area.C5_18, AT_data.H2W.C5_18,'*','MarkerSize',10);
h6 = plot(AT_data.cross_area.C6_18, AT_data.H2W.C6_18,'*','MarkerSize',10);
h7 = plot(AT_data.cross_area.C7_18, AT_data.H2W.C7_18,'*','MarkerSize',10);
h8 = plot(AT_data.cross_area.C8_18, AT_data.H2W.C8_18,'*','MarkerSize',10);
%h15 = plot(AT_data.cross_area.C9_18, AT_data.H2W.C9_18,'*','MarkerSize',10);
title('2018 Basal crevasse cross-sectional area vs. height-width value');
xlabel('cross sectional area (m^2)');
ylabel('Height-Width Ratio');
xlim([0 200]);
ylim([0 0.7]);
%legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
%  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
%  'Location', 'southeast');
hold off

%% Full Profiles with surface
figure(10)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on;
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11_Surf,'r');
h5 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14_Surf,'b');
h6 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18_Surf,'g');
title('Surface and Bed Elevation Profiles - Petermann Petermann Grouping 1');
xlabel('Along Track distance (km)');
ylabel('Elevation (m)');
ylim([-450 100])
legend('2011', '2014', '2018', ...
  'Location', 'southeast');
hold off

%% vertical elevation change profiles
figure(11)
subplot(2,1,1)
h1 = plot(AT_data.query_array.P11/1e3, (AT_data.interp_data.P14 - AT_data.interp_data.P11),'k');
hold on;
title('Vertical Elevation Change between 2011 and 2014');
xlabel('Along Track distance (km)');
ylabel('Elevation change (m)');
grid on
hold off

subplot(2,1,2)
h1 = plot(AT_data.query_array.P11/1e3, (AT_data.interp_data.P18 - AT_data.interp_data.P14),'k');
hold on
title('Vertical Elevation Change between 2014 and 2018');
xlabel('Along Track distance (km)');
ylabel('Elevation change (m)');
grid on
%% Basal Crevasse Plotting
figure(10)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on;
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.apex_pic_P11,'k*');
h5 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.apex_pic_P14,'k*');
h6 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.apex_pic_P18,'k*');
title('Basal Crevasse Apex Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
  'Apex Pick 2011', 'Apex Pick 2014', 'Apex Pick 2018', ...
  'Location', 'southeast');
hold off
%% all profiles apex picks
figure(11)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on;
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P11,'k*');
h5 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P14,'k*');
h6 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P18,'k*');
title('Basal Crevasse Base Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
  'Base Pick 2011', 'Base Pick 2014', 'Base Pick 2018', ...
  'Location', 'southeast');
hold off;

%% All profiles base picks
figure(25)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on;
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
title('Basal Crevasse Base Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'corrected 2014', 'corrected 2018', ...
  'Location', 'southeast');
hold off;
%% Annual Profiles with local min picks
figure(12)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11,'r');
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P11,'k*');
title('Basal Crevasse Base Picks (2011) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2011', 'Base Pick 2011', 'Location', 'southeast');

figure(13)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14,'b');
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P14,'k*');
title('Basal Crevasse Base Picks (2014) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2014', 'Base Pick 2014', 'Location', 'southeast');

figure(14)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18,'g');
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.crevasse.base_pic_P18,'k*');
title('Basal Crevasse Base Picks (2018) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2018', 'Base Pick 2018', 'Location', 'southeast');
%% Crevassse 1

figure(15)
h1 = plot(AT_data.along_track_crevasse.C1_11, AT_data.crevasse_pick.P11.crevasse1);
hold on
h2 = plot(AT_data.along_track_crevasse.C1_14, AT_data.crevasse_pick.P14.crevasse1);
h3 = plot(AT_data.along_track_crevasse.C1_18, AT_data.crevasse_pick.P18.crevasse1);
cross_area1_11 = patch(AT_data.along_track_crevasse.C1_11, AT_data.crevasse_pick.P11.crevasse1,'r','LineWidth',0.5);
cross_area1_14 = patch(AT_data.along_track_crevasse.C1_14, AT_data.crevasse_pick.P14.crevasse1,'b','LineWidth',0.5);
cross_area1_18 = patch(AT_data.along_track_crevasse.C1_18, AT_data.crevasse_pick.P18.crevasse1,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 1 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 1 P11', 'basal crevasse 1 P14', 'basal crevasse 1 P18', 'Location', 'northwest');
alpha(0.3)
hold off

%Crevasse 2
figure(16)
h1 = plot(AT_data.along_track_crevasse.C2_11, AT_data.crevasse_pick.P11.crevasse2);
hold on
h2 = plot(AT_data.along_track_crevasse.C2_14, AT_data.crevasse_pick.P14.crevasse2);
h3 = plot(AT_data.along_track_crevasse.C2_18, AT_data.crevasse_pick.P18.crevasse2);
cross_area2_11 = patch(AT_data.along_track_crevasse.C2_11, AT_data.crevasse_pick.P11.crevasse2,'r','LineWidth',0.5);
cross_area2_14 = patch(AT_data.along_track_crevasse.C2_14, AT_data.crevasse_pick.P14.crevasse2,'b','LineWidth',0.5);
cross_area2_18 = patch(AT_data.along_track_crevasse.C2_18, AT_data.crevasse_pick.P18.crevasse2,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 2 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 2 P11', 'basal crevasse 2 P14', 'basal crevasse 2 P18', 'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 3
figure(17)
h1 = plot(AT_data.along_track_crevasse.C3_1_11, AT_data.crevasse_pick.P11.crevasse3_1);
hold on
h2 = plot(AT_data.along_track_crevasse.C3_2_11, AT_data.crevasse_pick.P11.crevasse3_2);
h3 = plot(AT_data.along_track_crevasse.C3_3_11, AT_data.crevasse_pick.P11.crevasse3_3);
h4 = plot(AT_data.along_track_crevasse.C3_14, AT_data.crevasse_pick.P14.crevasse3);
h5 = plot(AT_data.along_track_crevasse.C3_18, AT_data.crevasse_pick.P18.crevasse3);
cross_area3_1_11 = patch(AT_data.along_track_crevasse.C3_1_11, AT_data.crevasse_pick.P11.crevasse3_1,'r','LineWidth',0.5);
cross_area3_2_11 = patch(AT_data.along_track_crevasse.C3_2_11, AT_data.crevasse_pick.P11.crevasse3_2,'r','LineWidth',0.5);
cross_area3_3_11 = patch(AT_data.along_track_crevasse.C3_3_11, AT_data.crevasse_pick.P11.crevasse3_3,'r','LineWidth',0.5);
cross_area3_14 = patch(AT_data.along_track_crevasse.C3_14, AT_data.crevasse_pick.P14.crevasse3,'b','LineWidth',0.5);
cross_area3_18 = patch(AT_data.along_track_crevasse.C3_18, AT_data.crevasse_pick.P18.crevasse3,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 3 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 3.1 P11', 'basal crevasse 3.2 P11', 'basal crevasse 3.3 P11', ...
    'basal crevasse 3 P14', 'basal crevasse 3 P18', 'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 4
figure(18)
h1 = plot(AT_data.along_track_crevasse.C4_11, AT_data.crevasse_pick.P11.crevasse4);
hold on
h2 = plot(AT_data.along_track_crevasse.C4_14, AT_data.crevasse_pick.P14.crevasse4);
h3 = plot(AT_data.along_track_crevasse.C4_18, AT_data.crevasse_pick.P18.crevasse4);
cross_area4_11 = patch(AT_data.along_track_crevasse.C4_11, AT_data.crevasse_pick.P11.crevasse4,'r','LineWidth',0.5);
cross_area4_14 = patch(AT_data.along_track_crevasse.C4_14, AT_data.crevasse_pick.P14.crevasse4,'b','LineWidth',0.5);
cross_area4_18 = patch(AT_data.along_track_crevasse.C4_18, AT_data.crevasse_pick.P18.crevasse4,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 2 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 4 P11', 'basal crevasse 4 P14', 'basal crevasse 4 P18', 'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 5
figure(19)
h1 = plot(AT_data.along_track_crevasse.C5_1_11, AT_data.crevasse_pick.P11.crevasse5_1);
hold on
h2 = plot(AT_data.along_track_crevasse.C5_2_11, AT_data.crevasse_pick.P11.crevasse5_2);
h3 = plot(AT_data.along_track_crevasse.C5_3_11, AT_data.crevasse_pick.P11.crevasse5_3);
h4 = plot(AT_data.along_track_crevasse.C5_1_14, AT_data.crevasse_pick.P14.crevasse5_1);
h5 = plot(AT_data.along_track_crevasse.C5_2_14, AT_data.crevasse_pick.P14.crevasse5_2);
h6 = plot(AT_data.along_track_crevasse.C5_18, AT_data.crevasse_pick.P18.crevasse5);
cross_area5_1_11 = patch(AT_data.along_track_crevasse.C5_1_11, AT_data.crevasse_pick.P11.crevasse5_1,'r','LineWidth',0.5);
cross_area5_2_11 = patch(AT_data.along_track_crevasse.C5_2_11, AT_data.crevasse_pick.P11.crevasse5_2,'r','LineWidth',0.5);
cross_area5_3_11 = patch(AT_data.along_track_crevasse.C5_3_11, AT_data.crevasse_pick.P11.crevasse5_3,'r','LineWidth',0.5);
cross_area5_1_14 = patch(AT_data.along_track_crevasse.C5_1_14, AT_data.crevasse_pick.P14.crevasse5_1,'b','LineWidth',0.5);
cross_area5_2_14 = patch(AT_data.along_track_crevasse.C5_2_14, AT_data.crevasse_pick.P14.crevasse5_2,'b','LineWidth',0.5);
cross_area5_18 = patch(AT_data.along_track_crevasse.C5_18, AT_data.crevasse_pick.P18.crevasse5,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 5 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 5.1 P11', 'basal crevasse 5.2 P11', 'basal crevasse 5.3 P11',...
    'basal crevasse 5.1 P14', 'basal crevasse 5.2 P14', 'basal crevasse 5 P18', ...
    'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 6
figure(20)
h1 = plot(AT_data.along_track_crevasse.C6_1_11, AT_data.crevasse_pick.P11.crevasse6_1);
hold on
h2 = plot(AT_data.along_track_crevasse.C6_2_11, AT_data.crevasse_pick.P11.crevasse6_2);
h3 = plot(AT_data.along_track_crevasse.C6_1_14, AT_data.crevasse_pick.P14.crevasse6_1);
h4 = plot(AT_data.along_track_crevasse.C6_2_14, AT_data.crevasse_pick.P14.crevasse6_2);
h5 = plot(AT_data.along_track_crevasse.C6_18, AT_data.crevasse_pick.P18.crevasse6);
cross_area6_1_11 = patch(AT_data.along_track_crevasse.C6_1_11, AT_data.crevasse_pick.P11.crevasse6_1,'r','LineWidth',0.5);
cross_area6_2_11 = patch(AT_data.along_track_crevasse.C6_2_11, AT_data.crevasse_pick.P11.crevasse6_2,'r','LineWidth',0.5);
cross_area6_1_14 = patch(AT_data.along_track_crevasse.C6_1_14, AT_data.crevasse_pick.P14.crevasse6_1,'r','LineWidth',0.5);
cross_area6_2_14 = patch(AT_data.along_track_crevasse.C6_2_14, AT_data.crevasse_pick.P14.crevasse6_2,'b','LineWidth',0.5);
cross_area6_18 = patch(AT_data.along_track_crevasse.C6_18, AT_data.crevasse_pick.P18.crevasse6,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 6 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 6.1 P11', 'basal crevasse 6.2 P11', 'basal crevasse 6.1 P14', ...
    'basal crevasse 6.2 P14', 'basal crevasse 6 P18', 'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 7
figure(21)
h1 = plot(AT_data.along_track_crevasse.C7_11, AT_data.crevasse_pick.P11.crevasse7);
hold on
h2 = plot(AT_data.along_track_crevasse.C7_14, AT_data.crevasse_pick.P14.crevasse7);
h3 = plot(AT_data.along_track_crevasse.C7_18, AT_data.crevasse_pick.P18.crevasse7);
cross_area7_11 = patch(AT_data.along_track_crevasse.C7_11, AT_data.crevasse_pick.P11.crevasse7,'r','LineWidth',0.5);
cross_area7_14 = patch(AT_data.along_track_crevasse.C7_14, AT_data.crevasse_pick.P14.crevasse7,'b','LineWidth',0.5);
cross_area7_18 = patch(AT_data.along_track_crevasse.C7_18, AT_data.crevasse_pick.P18.crevasse7,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 7 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 7 P11', 'basal crevasse 7 P14', 'basal crevasse 7 P18', 'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 8
figure(22)
h1 = plot(AT_data.along_track_crevasse.C8_11, AT_data.crevasse_pick.P11.crevasse8);
hold on
h2 = plot(AT_data.along_track_crevasse.C8_14, AT_data.crevasse_pick.P14.crevasse8);
h3 = plot(AT_data.along_track_crevasse.C8_18, AT_data.crevasse_pick.P18.crevasse8);
cross_area8_11 = patch(AT_data.along_track_crevasse.C8_11, AT_data.crevasse_pick.P11.crevasse8,'r','LineWidth',0.5);
cross_area8_14 = patch(AT_data.along_track_crevasse.C8_14, AT_data.crevasse_pick.P14.crevasse8,'b','LineWidth',0.5);
cross_area8_18 = patch(AT_data.along_track_crevasse.C8_18, AT_data.crevasse_pick.P18.crevasse8,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 8 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 8 P11', 'basal crevasse 8 P14', 'basal crevasse 8 P18', 'Location', 'northwest');
alpha(0.3)
hold off

% Crevasse 9
figure(22)
h1 = plot(AT_data.along_track_crevasse.C9_11, AT_data.crevasse_pick.P11.crevasse9);
hold on
h2 = plot(AT_data.along_track_crevasse.C9_14, AT_data.crevasse_pick.P14.crevasse9);
%h3 = plot(AT_data.along_track_crevasse.C9_18, AT_data.crevasse_pick.P18.crevasse9);
cross_area9_11 = patch(AT_data.along_track_crevasse.C9_11, AT_data.crevasse_pick.P11.crevasse9,'r','LineWidth',0.5);
cross_area9_14 = patch(AT_data.along_track_crevasse.C9_14, AT_data.crevasse_pick.P14.crevasse9,'b','LineWidth',0.5);
%cross_area9_18 = patch(AT_data.along_track_crevasse.C9_18, AT_data.crevasse_pick.P18.crevasse9,'g','LineWidth',0.5);
title('Basal Crevasse Base Pick 9 (2011,2014,2018)');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 9 P11', 'basal crevasse 9 P14', 'Location', 'northwest');
alpha(0.3)
hold off
