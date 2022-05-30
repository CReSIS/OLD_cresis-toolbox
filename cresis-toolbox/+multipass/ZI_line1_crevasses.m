%% Basal Crevasse Apex Picking section:
% Crevasse Apex calculation and binary file formation
AT_data.crevasse.apex_bin.P10 = islocalmax(AT_data.interp_data.P10);
AT_data.crevasse.apex_bin.P14A = islocalmax(AT_data.interp_data.P14A);
AT_data.crevasse.apex_bin.P16 = islocalmax(AT_data.interp_data.P16);
AT_data.crevasse.apex_bin.P17 = islocalmax(AT_data.interp_data.P17);
AT_data.crevasse.apex_bin.P18 = islocalmax(AT_data.interp_data.P18);
AT_data.crevasse.apex_bin.P19 = islocalmax(AT_data.interp_data.P19);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
AT_data.crevasse.apex_pic_P10 = AT_data.interp_data.P10.*AT_data.crevasse.apex_bin.P10;
AT_data.crevasse.apex_pic_P10(AT_data.crevasse.apex_pic_P10 == 0) = NaN;
AT_data.crevasse.apex_pic_P14A = AT_data.interp_data.P14A.*AT_data.crevasse.apex_bin.P14A;
AT_data.crevasse.apex_pic_P14A(AT_data.crevasse.apex_pic_P14 == 0) = NaN;
AT_data.crevasse.apex_pic_P16 = AT_data.interp_data.P16.*AT_data.crevasse.apex_bin.P16;
AT_data.crevasse.apex_pic_P16(AT_data.crevasse.apex_pic_P16 == 0) = NaN;
AT_data.crevasse.apex_pic_P17 = AT_data.interp_data.P17.*AT_data.crevasse.apex_bin.P17;
AT_data.crevasse.apex_pic_P17(AT_data.crevasse.apex_pic_P17 == 0) = NaN;
AT_data.crevasse.apex_pic_P18 = AT_data.interp_data.P18.*AT_data.crevasse.apex_bin.P18;
AT_data.crevasse.apex_pic_P18(AT_data.crevasse.apex_pic_P18 == 0) = NaN;
AT_data.crevasse.apex_pic_P19 = AT_data.interp_data.P19.*AT_data.crevasse.apex_bin.P19;
AT_data.crevasse.apex_pic_P19(AT_data.crevasse.apex_pic_P19 == 0) = NaN;

%% Basal Crevasse Wall Base Picking:
% Crevasse wall base calculation and binary file formation
AT_data.crevasse.base_bin.P10 = islocalmin(AT_data.interp_data.P10);
AT_data.crevasse.base_bin.P14A = islocalmin(AT_data.interp_data.P14A);
AT_data.crevasse.base_bin.P16 = islocalmin(AT_data.interp_data.P16);
AT_data.crevasse.base_bin.P17 = islocalmin(AT_data.interp_data.P17);
AT_data.crevasse.base_bin.P18 = islocalmin(AT_data.interp_data.P18);
AT_data.crevasse.base_bin.P19 = islocalmin(AT_data.interp_data.P19);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
AT_data.crevasse.base_pic_P10 = AT_data.interp_data.P10.*AT_data.crevasse.base_bin.P10;
AT_data.crevasse.base_pic_P10(AT_data.crevasse.base_pic_P10 == 0) = NaN;
AT_data.crevasse.base_pic_P14A = AT_data.interp_data.P14A.*AT_data.crevasse.base_bin.P14A;
AT_data.crevasse.base_pic_P14A(AT_data.crevasse.base_pic_P14A == 0) = NaN;
AT_data.crevasse.base_pic_P16 = AT_data.interp_data.P16.*AT_data.crevasse.base_bin.P16;
AT_data.crevasse.base_pic_P16(AT_data.crevasse.base_pic_P16 == 0) = NaN;
AT_data.crevasse.base_pic_P17 = AT_data.interp_data.P17.*AT_data.crevasse.base_bin.P17;
AT_data.crevasse.base_pic_P17(AT_data.crevasse.base_pic_P17 == 0) = NaN;
AT_data.crevasse.base_pic_P18 = AT_data.interp_data.P18.*AT_data.crevasse.base_bin.P18;
AT_data.crevasse.base_pic_P18(AT_data.crevasse.base_pic_P18 == 0) = NaN;
AT_data.crevasse.base_pic_P19 = AT_data.interp_data.P19.*AT_data.crevasse.base_bin.P19;
AT_data.crevasse.base_pic_P19(AT_data.crevasse.base_pic_P19 == 0) = NaN;

%% AFTER MANUAL EDITTING EXPORT AND SAVE AS CSV
% after mannually editting of each profile, save as new variable and export
% variable as a csv 
% DO NO DO IF CSV ALREADY EXISTS
% COMMENT SECTION OUT AFTER FINALS PICKS AS TO NO OVERWRITE EXISTING FILES
crevasse_2010 = AT_data.crevasse.base_pic_P10;
crevasse_2014A = AT_data.crevasse.base_pic_P14A;
crevasse_2016 = AT_data.crevasse.base_pic_P16;
crevasse_2017 = AT_data.crevasse.base_pic_P17;
crevasse_2018 = AT_data.crevasse.base_pic_P18;
crevasse_2019 = AT_data.crevasse.base_pic_P19;
writematrix(crevasse_2010, 'ZI1_Crevasse_2010.csv','Delimiter','comma');
writematrix(crevasse_2014A, 'ZI1_Crevasse_2014A.csv','Delimiter','comma');
writematrix(crevasse_2016, 'ZI1_Crevasse_2016.csv','Delimiter','comma');
writematrix(crevasse_2017, 'ZI1_Crevasse_2017.csv','Delimiter','comma');
writematrix(crevasse_2018, 'ZI1_Crevasse_2018.csv','Delimiter','comma');
writematrix(crevasse_2019, 'ZI1_Crevasse_2019.csv','Delimiter','comma');
%% Reload new csv in after editting
AT_data.crevasse.base_pic_P10 = csvread('ZI1_Crevasse_2010.csv');
AT_data.crevasse.base_pic_P14A = csvread('ZI1_Crevasse_2014A.csv');
AT_data.crevasse.base_pic_P16 = csvread('ZI1_Crevasse_2016.csv');
AT_data.crevasse.base_pic_P17 = csvread('ZI1_Crevasse_2017.csv');
AT_data.crevasse.base_pic_P18 = csvread('ZI1_Crevasse_2018.csv');
AT_data.crevasse.base_pic_P19 = csvread('ZI1_Crevasse_2019.csv');
%% 2010 Basal Crevasse picks

%% 2014 Basal Crevasse picks

%% 2016 Basal Crevasse picks

%% 2017 Basal Crevasse picks

%% 2018 Basal Crevasse Picks

%% 2019 Basal Crevasse picks







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
title('Basal Crevasse Base Picks (2014) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2018', 'Base Pick 2018', 'Location', 'southeast');

crevasse1 = AT_data.interp_data.P11(479:497);
along_track_crevasse1 = AT_data.query_array.P11(479:497)/1e3;
cross_area = cumtrapz(along_track_crevasse1, crevasse1);
figure(15)
h1 = plot(along_track_crevasse1, crevasse1);
cross_area1 = patch(along_track_crevasse1, crevasse1,'b','LineWidth',1.5);
title('Basal Crevasse Base Pick 1 (2011) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 1', 'Location', 'southeast');

