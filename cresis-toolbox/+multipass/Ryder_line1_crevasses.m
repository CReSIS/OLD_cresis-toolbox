%% Basal Crevasse Apex Picking section:
% Crevasse Apex calculation and binary file formation
%AT_data.crevasse.apex_bin.P07 = islocalmax(AT_data.interp_data.P07);
AT_data.crevasse.apex_bin.P11 = islocalmax(AT_data.interp_data.P11);
AT_data.crevasse.apex_bin.P13 = islocalmax(AT_data.interp_data.P13);
AT_data.crevasse.apex_bin.P15 = islocalmax(AT_data.interp_data.P15);
AT_data.crevasse.apex_bin.P19 = islocalmax(AT_data.interp_data.P19);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
%AT_data.crevasse.apex_pic_P07 = AT_data.interp_data.P07.*AT_data.crevasse.apex_bin.P07;
%AT_data.crevasse.apex_pic_P07(AT_data.crevasse.apex_pic_P07 == 0) = NaN;
AT_data.crevasse.apex_pic_P11 = AT_data.interp_data.P11.*AT_data.crevasse.apex_bin.P11;
AT_data.crevasse.apex_pic_P11(AT_data.crevasse.apex_pic_P11 == 0) = NaN;
AT_data.crevasse.apex_pic_P13 = AT_data.interp_data.P13.*AT_data.crevasse.apex_bin.P13;
AT_data.crevasse.apex_pic_P13(AT_data.crevasse.apex_pic_P13 == 0) = NaN;
AT_data.crevasse.apex_pic_P15 = AT_data.interp_data.P15.*AT_data.crevasse.apex_bin.P15;
AT_data.crevasse.apex_pic_P15(AT_data.crevasse.apex_pic_P15 == 0) = NaN;
AT_data.crevasse.apex_pic_P19 = AT_data.interp_data.P19.*AT_data.crevasse.apex_bin.P19;
AT_data.crevasse.apex_pic_P19(AT_data.crevasse.apex_pic_P19 == 0) = NaN;

%% Basal Crevasse Wall Base Picking:
% Crevasse wall base calculation and binary file formation
%AT_data.crevasse.base_bin.P07 = islocalmin(AT_data.interp_data.P07);
AT_data.crevasse.base_bin.P11 = islocalmin(AT_data.interp_data.P11);
AT_data.crevasse.base_bin.P13 = islocalmin(AT_data.interp_data.P13);
AT_data.crevasse.base_bin.P15 = islocalmin(AT_data.interp_data.P15);
AT_data.crevasse.base_bin.P19 = islocalmin(AT_data.interp_data.P19);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
%AT_data.crevasse.base_pic_P07 = AT_data.interp_data.P07.*AT_data.crevasse.base_bin.P07;
AT_data.crevasse.base_pic_P07(AT_data.crevasse.base_pic_P07 == 0) = NaN;
AT_data.crevasse.base_pic_P11 = AT_data.interp_data.P11.*AT_data.crevasse.base_bin.P11;
AT_data.crevasse.base_pic_P11(AT_data.crevasse.base_pic_P11 == 0) = NaN;
AT_data.crevasse.base_pic_P13 = AT_data.interp_data.P13.*AT_data.crevasse.base_bin.P13;
AT_data.crevasse.base_pic_P13(AT_data.crevasse.base_pic_P13 == 0) = NaN;
AT_data.crevasse.base_pic_P15 = AT_data.interp_data.P15.*AT_data.crevasse.base_bin.P15;
AT_data.crevasse.base_pic_P15(AT_data.crevasse.base_pic_P15 == 0) = NaN;
AT_data.crevasse.base_pic_P19 = AT_data.interp_data.P19.*AT_data.crevasse.base_bin.P19;
AT_data.crevasse.base_pic_P19(AT_data.crevasse.base_pic_P19 == 0) = NaN;

%% AFTER MANUAL EDITTING EXPORT AND SAVE AS CSV
% after mannually editting of each profile, save as new variable and export
% variable as a csv 
% DO NO DO IF CSV ALREADY EXISTS
% COMMENT SECTION OUT AFTER FINALS PICKS AS TO NO OVERWRITE EXISTING FILES
%crevasse_2007 = AT_data.crevasse.base_pic_P07;
crevasse_2011 = AT_data.crevasse.base_pic_P11;
crevasse_2013 = AT_data.crevasse.base_pic_P13;
crevasse_2015 = AT_data.crevasse.base_pic_P15;
crevasse_2019 = AT_data.crevasse.base_pic_P19;
%writematrix(crevasse_2007, 'R1_Crevasse_2007.csv','Delimiter','comma');
writematrix(crevasse_2011, 'R1_Crevasse_2011.csv','Delimiter','comma');
writematrix(crevasse_2013, 'R1_Crevasse_2013.csv','Delimiter','comma');
writematrix(crevasse_2015, 'R1_Crevasse_2015.csv','Delimiter','comma');
writematrix(crevasse_2019, 'R1_Crevasse_2019.csv','Delimiter','comma');
%% Reload new csv in after editting
%AT_data.crevasse.base_pic_P07 = csvread('R1_Crevasse_2007.csv');
AT_data.crevasse.base_pic_P11 = csvread('R1_Crevasse_2011.csv');
AT_data.crevasse.base_pic_P13 = csvread('R1_Crevasse_2013.csv');
AT_data.crevasse.base_pic_P15 = csvread('R1_Crevasse_2015.csv');
AT_data.crevasse.base_pic_P19 = csvread('R1_Crevasse_2019.csv');

%% 2007 Basal Crevasse picks

%% 2011 Basal Crevasse picks

%% 2013 Basal Crevasse picks

%% 2015 Basal Crevasse picks

%% 2019 Basal Crevasse Picks







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

