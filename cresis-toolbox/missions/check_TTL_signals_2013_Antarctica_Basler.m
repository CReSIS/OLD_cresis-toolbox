% check_TTL_signals_2013_Antarctica_Basler
%
% Plots the TTL signals.
% Also stores this into a Matlab file.
%
% Author: John Paden

meas_paths = {};

% ===================================================================
% User Settings
% ===================================================================
fn = 'C:\Users\dangermo\Documents\Travel\Antarctica_2013\mcords4_fasttimegain_20130915\CH1.csv';
csv_format = '%f %f %f %f %f';
oscilloscope_bin_rng = []; % Plot 

% ===================================================================
% Automated Section
% ===================================================================
fid = fopen(fn,'r');
C = textscan(fid,csv_format,'HeaderLines',19,'Delimiter',',');
fclose(fid);
tx_signal_time = C{1};
tx_PRI = C{2};
tx_TTL0 = C{3};
tx_DDS = C{4};
clear C;

figure(1); clf;
fig_pos = get(1,'Position');
h_PRI = plot(tx_signal_time*1e6,tx_PRI, 'r');
hold on;
h_TTL0 = plot(tx_signal_time*1e6,tx_TTL0, 'k');
h_DDS = plot(tx_signal_time*1e6,tx_DDS, 'b');
hold off;
grid on;
xlabel('DSO time (us)');
ylabel('volts (V)');
legend([h_PRI h_TTL0 h_DDS], {'PRI', 'TTL0', 'DDS'});

