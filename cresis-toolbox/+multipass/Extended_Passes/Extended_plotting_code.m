
%% Annual melt and topography grouping epochs 1-5 line P1
fig1_2 = figure(201);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
ylim([-100 200]);
xlim([0 60]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 1 2011-2014')
legend('2011', '2014', 'Location','southeast');
han = axes(fig1_2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig2_2 = figure(202);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 1 2014-2015')
legend('2014', '2015', 'Location','southeast');
han = axes(fig2_2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig3_2 = figure(203);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 1 2015-2017')
legend('2015', '2017', 'Location','southeast');
han = axes(fig3_2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig4_2 = figure(204);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 1 2017-2018')
legend('2017', '2018', 'Location','southeast');
han = axes(fig2_2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig5_2 = figure(205);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(6).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(6).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 1 2018-2019')
legend('2018', '2019', 'Location','southeast');
han = axes(fig5_2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig6_2 = figure(206);
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(1).Interp_Surf_Corrected(379:end), 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(1).Interp_Bed_Corrected(379:end), 'color','r');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(2).Interp_Surf_Corrected(379:end), 'color','b');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(2).Interp_Bed_Corrected(379:end), 'color','b');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(3).Interp_Surf_Corrected(379:end), 'color','k');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(3).Interp_Bed_Corrected(379:end), 'color','k');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(4).Interp_Surf_Corrected(379:end), 'color','g');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(4).Interp_Bed_Corrected(379:end), 'color','g');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(5).Interp_Surf_Corrected(379:end), 'color','c');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(5).Interp_Bed_Corrected(379:end), 'color','c');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(6).Interp_Surf_Corrected(379:end), 'color','m');
plot(XY_data(1).Interp_AT_Dist(379:end)/1e3, XY_data(6).Interp_Bed_Corrected(379:end), 'color','m');
title('Elevation Profiles 2011-2019');
ylim([-500 100]);
xlim([0 60]);
legend('2011','2014','2015','2017','2018','2019','Location','southeast');
ylabel('elevation (m)');
xlabel('Along Track Distance (km)');

fig7_2 = figure(207);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Melt_SMB_yr);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB_yr);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Melt_SMB_yr);
title('Melt Profiles 2011-2019');
ylim([-200 200]);
xlim([0 60]);
legend('2011-2014','2014-2015','2015-2017','2017-2018','2018-2019','Location','southwest');
ylabel('melt rate (m/yr)');
xlabel('Along Track Distance (km)');
%% P1 Histograms
figure8_2 = figure(208);
histogram(XY_data(1).Melt_SMB_yr(379:3776),50,'BinWidth',2); % was 3960
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P1 Melt Frequency 2011-2014');

figure9_2 = figure(209);
histogram(XY_data(2).Melt_SMB_yr(379:3776),50,'BinWidth',2);
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P1 Melt Frequency 2014-2015');

figure10_2 = figure(210);
histogram(XY_data(3).Melt_SMB_yr(379:3776),50,'BinWidth',2); % was 3840
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P1 Melt Frequency 2015-2017');

figure11_2 = figure(211);
histogram(XY_data(4).Melt_SMB_yr(379:3511),50,'BinWidth',2); %3511
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P1 Melt Frequency 2017-2018');

figure12_2 = figure(212);
histogram(XY_data(5).Melt_SMB_yr(379:3511),50,'BinWidth',2); %3721
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P1 Melt Frequency 2018-2019');
% histogram(XY_data(2).Melt_SMB_yr);
% histogram(XY_data(3).Melt_SMB_yr);
% histogram(XY_data(4).Melt_SMB_yr);
% histogram(XY_data(5).Melt_SMB_yr);
%% Annual melt and topography grouping epochs 1-2 line P2
fig1_3 = figure(301);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
ylim([-50 100]);
xlim([0 60]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 2 2013-2014')
legend('2013', '2014', 'Location','southeast');
han = axes(fig1_3,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig2_3 = figure(302);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-50 100]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
xlim([0 60]);
suptitle('Petermann Line 2 2014-2017')
legend('2014', '2017', 'Location','southeast');
han = axes(fig2_3,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig3_3 = figure(303);
plot(XY_data(1).Interp_AT_Dist(343:end)/1e3, XY_data(1).Interp_Surf_Corrected(343:end), 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist(343:end)/1e3, XY_data(1).Interp_Bed_Corrected(343:end), 'color','r');
plot(XY_data(1).Interp_AT_Dist(343:end)/1e3, XY_data(2).Interp_Surf_Corrected(343:end), 'color','b');
plot(XY_data(1).Interp_AT_Dist(343:end)/1e3, XY_data(2).Interp_Bed_Corrected(343:end), 'color','b');
plot(XY_data(1).Interp_AT_Dist(343:end)/1e3, XY_data(3).Interp_Surf_Corrected(343:end), 'color','k');
plot(XY_data(1).Interp_AT_Dist(343:end)/1e3, XY_data(3).Interp_Bed_Corrected(343:end), 'color','k');
title('Elevation Profiles 2013-2017');
ylim([-500 100]);
xlim([0 60]);
legend('2013','2014','2017','Location','southeast');
ylabel('elevation (m)');
xlabel('Along Track Distance (km)');

fig3_4 = figure(304);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
title('Melt Profiles 2013-2017');
ylim([-50 100]);
xlim([0 60]);
legend('2013-2014','2014-2017','Location','southwest');
ylabel('melt rate (m/yr)');
xlabel('Along Track Distance (km)');

%% P4 figures
fig4_1 = figure(401);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
ylim([-100 200]);
xlim([0 60]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-600 100]);
xlim([0 60]);
suptitle('Petermann Line 4 2010-2011')
legend('2010', '2011', 'Location','southeast');
han = axes(fig4_1,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig4_2 = figure(402);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-600 100]);
xlim([0 60]);
suptitle('Petermann Line 4 2011-2013')
legend('2011', '2013', 'Location','southeast');
han = axes(fig4_2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig4_3 = figure(403);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-600 100]);
xlim([0 60]);
suptitle('Petermann Line 4 2013-2014')
legend('2013', '2014', 'Location','southeast');
han = axes(fig4_3,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

fig4_4 = figure(404);
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB_yr);
ylabel('Basal Melt (m)');
ylim([-100 200]);
xlim([0 60]);
title('Annual Basal Melt');
legend('Time offset thickness','Location','southeast');
subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-600 100]);
xlim([0 60]);
suptitle('Petermann Line 1 2014-2017')
legend('2014', '2017', 'Location','southeast');
han = axes(fig4_4,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');


fig4_5 = figure(405);
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(1).Interp_Surf_Corrected(419:3083), 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(1).Interp_Bed_Corrected(419:3083), 'color','r');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(2).Interp_Surf_Corrected(419:3083), 'color','b');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(2).Interp_Bed_Corrected(419:3083), 'color','b');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(3).Interp_Surf_Corrected(419:3083), 'color','k');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(3).Interp_Bed_Corrected(419:3083), 'color','k');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(4).Interp_Surf_Corrected(419:3083), 'color','g');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(4).Interp_Bed_Corrected(419:3083), 'color','g');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(5).Interp_Surf_Corrected(419:3083), 'color','c');
plot(XY_data(1).Interp_AT_Dist(419:3083)/1e3, XY_data(5).Interp_Bed_Corrected(419:3083), 'color','c');
title('Elevation Profiles 2011-2019');
ylim([-600 100]);
xlim([0 60]);
legend('2010','2011','2013','2014','2017','Location','southeast');
ylabel('elevation (m)');
xlabel('Along Track Distance (km)');

fi4_5 = figure(406);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Melt_SMB_yr);
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB_yr);
title('Melt Profiles 2011-2019');
ylim([-200 200]);
xlim([0 60]);
legend('2010-2011','2011-2013','2013-2014','2014-2017','Location','southwest');
ylabel('melt rate (m/yr)');
xlabel('Along Track Distance (km)');


%% P3 Histograms
figure4_7 = figure(407);
histogram(XY_data(1).Melt_SMB_yr(419:3083),50,'BinWidth',2); % was 3960
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P4 Melt Frequency 2010-2011');

figure4_8 = figure(408);
histogram(XY_data(2).Melt_SMB_yr(419:3083),50,'BinWidth',2);
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P4 Melt Frequency 2011-2013');

figure4_9 = figure(409);
histogram(XY_data(3).Melt_SMB_yr(419:3083),50,'BinWidth',2); % was 3840
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P4 Melt Frequency 2013-2014');

figure4_10 = figure(410);
histogram(XY_data(4).Melt_SMB_yr(419:3083),50,'BinWidth',2); %3511
hold on 
xlabel('Melt Rate');
ylabel('Frequency');
ylim([0 600])
xlim([-75 200])
title('P4 Melt Frequency 2014-2017');