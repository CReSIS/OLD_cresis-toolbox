%% Plot everything
clearvars -except gRadar tmp
load('PhaseSlope_Data_3Chunk.mat')
cls = physconst('lightspeed');
savedir='/cresis/snfs1/scratch/bmiller/MultiPass_Results_12_1/Official_Paper_Results';

% % Plot each line together
% leg2 = {'Top', 'Mid', 'Bot'}; legger={};
% linsty = {'-','--','.-'}; lincol = {'r','b','g'};
% figure(50)
% for pid = 1:length(fns)
%   colnow = lincol{pid};
%   dispnow = phaseout{pid};
%   for id_ch = 1:numchunks
%     stynow = linsty{id_ch};
%     ynow = dispnow(id_ch,:)./middleind(id_ch);
%     plot(alongx, ynow, [colnow stynow])
%     hold on
%     legger{end+1} = sprintf('%s %s',leg{pid},leg2{id_ch});
%   end
% end
% ylabel('radians by range bine (rad/~)')
% xlabel('Along Track')
% title('Depth Change Across Years')
% legend(legger)
% grid on
% hold off


filt_disp = 301;
%Plot each line together
leg2 = {'Top', 'Mid', 'Bot'}; legger={};
linsty = {'-','--',':'}; lincol = {'r','b','g','c','m','k'};
figure(60)
for pid = fliplr(1:length(fns))
  colnow = lincol{pid};
  dispnow = disp_off{pid};
  %Extra averaging of dispout
  dispnow = fir_dec(dispnow,ones(1,filt_disp)/filt_disp);
  for id_ch = 1:numchunks
    stynow = linsty{id_ch};
%     ynow = dispnow(id_ch,:)./crossy(id_ch);
    ynow = dispnow(id_ch,:)./(cls/(2*BW*sqrt(eps_ice)));
    plot(alongx{pid}, ynow, [colnow stynow])
    hold on
    legger{end+1} = sprintf('%s %s',leg{pid},leg2{id_ch});
  end
end
ylabel('Depth change by Depth (m\cdotm^{-1})')
xlabel('Along Track (km)')
title('Depth Change Across Years')
legend(legger,'Position',[0.8078    0.1196    0.1825    0.7972])
set(findall(gcf,'type','line'),'LineWidth',2)
set(get(gcf,'Children'),'fontsize',14)
set(gcf,'Position',[147         236        1264         563])
grid on
hold off
savefn = '3Chunk_PhaseEstimate_DepthbyDepth';
saveas(gcf,fullfile(savedir,[savefn '.fig']))
saveas(gcf,fullfile(savedir,[savefn '.png']))

%Plot the delta of depth change
%Plot each line together
leg2 = {'Top', 'Mid', 'Bot'}; legger={};
figure(61)
for pid = fliplr(1:length(fns))
  colnow = lincol{pid};
  dispnow = disp_off{pid};
  %Extra averaging of dispout
  dispnow = fir_dec(dispnow,ones(1,filt_disp)/filt_disp);
  for id_ch = 1:numchunks
    stynow = linsty{id_ch};
%     ynow = dispnow(id_ch,:)./crossy(id_ch);
    ynow = dispnow(id_ch,:);
    plot(alongx{pid}, ynow, [colnow stynow])
    hold on
    legger{end+1} = sprintf('%s %s',leg{pid},leg2{id_ch});
  end
end
ylabel('Depth Change (m)')
xlabel('Along Track (km)')
title('Depth Change Across Years')
legend(legger,'Position',[0.8078    0.1196    0.1825    0.7972])
set(findall(gcf,'type','line'),'LineWidth',2)
set(get(gcf,'Children'),'fontsize',14)
set(gcf,'Position',[147         236        1264         563])
grid on
hold off
savefn = '3Chunk_PhaseEstimate_DeltaDepth';
saveas(gcf,fullfile(savedir,[savefn '.fig']))
saveas(gcf,fullfile(savedir,[savefn '.png']))

% %Plot the phase output
% for pid = 1:length(fns)
%   figure(30+pid)
%   imagesc(alongx, crossy, phaseout{pid})
% %   [Con, h] = contourf(phaseout{pid},100);
% %   set(h,'LineColor','none')
%   c = colorbar;
%   c.Label.String = 'Phase (rad)';
%   xlabel('Along Track (m)')
%   ylabel('Depth (m)')
%   title(sprintf('%s to 2014',leg{pid}))
%   caxis([-0.07 0])
% end
% %Plot the displacement output
% for pid = 1:length(fns)
%   figure(40+pid)
%   imagesc(alongx, crossy, disp_off{pid})
%   c = colorbar;
%   c.Label.String = 'Displacement (m)';
%   xlabel('Along Track (km)')
%   ylabel('Depth (m)')
%   title(sprintf('%s to 2014',leg{pid}))
%   caxis([-10 0]*1e-3)
% end

