% Temporal structure of individual cold pools ranked by their change in dD
load('workspace_17_cp_detection_algorithm_11min.mat');
% load('workspace_17_cp_detection_algorithm_1min.mat');
clearvars -except Taf qair wspd rh prec u v t_min_ind ...
                  t_max_ind t_end_ind t1min factor iso_cp_matrix cp_matrix ...
                  dD d18O q_iso DXS ship orange green blue mustard ...
                  Taf_comp2 t_comp2 dD_comp2 d18O_comp2 q_iso_comp2 DXS_comp2...
                  qair_comp2 wspd_comp2 rh_comp2 prec_comp2 u_comp2 v_comp2
i17 = [65,67:78,82:84,86];
Taf_comp2 = Taf_comp2(i17,:);
wspd_comp2 = wspd_comp2(i17,:);
qair_comp2 = qair_comp2(i17,:);
prec_comp2 = prec_comp2(i17,:);
rh_comp2 = rh_comp2(i17,:);
u_comp2 = u_comp2(i17,:);
v_comp2 = v_comp2(i17,:);

j17 = [1,3:14,18:20,22];
d18O_comp2 = d18O_comp2(j17,:);
dD_comp2 = dD_comp2(j17,:);
DXS_comp2 = DXS_comp2(j17,:);

subset = [1:3, 5:8, 11:17]; % cuts to 14 events

%% Identifying the strongest and weakest cold pools based on dD
num_vector = i17;
delta_dD = dD(t_max_ind(num_vector)-factor) - dD(t_min_ind(num_vector)-factor);
cp_vector = [delta_dD, num_vector', (1:length(num_vector))'];
sorted_cp = sortrows(cp_vector,1,'descend');
cp_order = sorted_cp(:,3);
cp_order100 = sorted_cp(:,2);
sz = size(Taf_comp2);
Taf_comp2_sort = NaN(sz);
wspd_comp2_sort = NaN(sz);
qair_comp2_sort = NaN(sz);
prec_comp2_sort = NaN(sz);
d18O_comp2_sort = NaN(sz);
dD_comp2_sort = NaN(sz);
DXS_comp2_sort = NaN(sz);

% sort all 17 cold pools
% Anomalies from values at t_max!!!
for k = 1:length(cp_order)
    Taf_comp2_sort(k,:) = Taf_comp2(cp_order(k),:) -Taf(t_max_ind(cp_order100(k)));
    prec_comp2_sort(k,:)= prec_comp2(cp_order(k),:);
    qair_comp2_sort(k,:)= qair_comp2(cp_order(k),:)-qair(t_max_ind(cp_order100(k)));
    wspd_comp2_sort(k,:)= wspd_comp2(cp_order(k),:)-wspd(t_max_ind(cp_order100(k)));
    DXS_comp2_sort(k,:) = DXS_comp2(cp_order(k),:) -DXS(t_max_ind(cp_order100(k))-factor);
    dD_comp2_sort(k,:)  = dD_comp2(cp_order(k),:)  -dD(t_max_ind(cp_order100(k))-factor);
    d18O_comp2_sort(k,:)= d18O_comp2(cp_order(k),:)-d18O(t_max_ind(cp_order100(k))-factor);
end

% % Anomalies from mean values!!!
% for k = 1:length(cp_order)
%     Taf_comp2_sort(k,:) = Taf_comp2(cp_order(k),:)-Taf(t_max_ind(cp_order100(k)));
%     prec_comp2_sort(k,:)= prec_comp2(cp_order(k),:);
%     qair_comp2_sort(k,:)= qair_comp2(cp_order(k),:)-qairf(t_max_ind(cp_order100(k)));
%     wspd_comp2_sort(k,:)= wspd_comp2(cp_order(k),:)-wspdf(t_max_ind(cp_order100(k)));
%     DXS_comp2_sort(k,:) = DXS_comp2(cp_order(k),:)-DXS(t_max_ind(cp_order100(k))-factor);
%     dD_comp2_sort(k,:)  = dD_comp2(cp_order(k),:) -dD(t_max_ind(cp_order100(k))-factor);
%     d18O_comp2_sort(k,:)= d18O_comp2(cp_order(k),:)-d18O(t_max_ind(cp_order100(k))-factor);
% end

% Eliminating NaN values and duplicating last row and last column for inclusion
% in pcolor function
ncp = length(num_vector);
Taf_comp2_sort(ncp+1,:) = Taf_comp2_sort(ncp,:);
Taf_comp2_sort(:,112) = Taf_comp2_sort(:,111);
prec_comp2_sort(ncp+1,:) = prec_comp2_sort(ncp,:);
prec_comp2_sort(:,112) = prec_comp2_sort(:,111);
qair_comp2_sort(ncp+1,:) = qair_comp2_sort(ncp,:);
qair_comp2_sort(:,112) = qair_comp2_sort(:,111);
wspd_comp2_sort(ncp+1,:) = wspd_comp2_sort(ncp,:);
wspd_comp2_sort(:,112) = wspd_comp2_sort(:,111);
dD_comp2_sort(ncp+1,:) = dD_comp2_sort(ncp,:);
dD_comp2_sort(:,112) = dD_comp2_sort(:,111);
DXS_comp2_sort(ncp+1,:) = DXS_comp2_sort(ncp,:);
DXS_comp2_sort(:,112) = DXS_comp2_sort(:,111);
d18O_comp2_sort(ncp+1,:) = d18O_comp2_sort(ncp,:);
d18O_comp2_sort(:,112) = d18O_comp2_sort(:,111);

%% load colormaps
try
    addpath('/Users/estefania/Documents/Data/cbrewer')
    [tmap]=cbrewer('div', 'BrBG', 12);
    [bw]=cbrewer('seq', 'Greys', 12);
    % Cite as:
    % Charles (2021). cbrewer : colorbrewer schemes for Matlab 
    % (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab), 
    % MATLAB Central File Exchange. Retrieved December 2, 2021.
    % -> Not available anymore, and Simon doesn't have it.
catch
    % Scott Lowe (2025). cbrewer2 (https://github.com/scottclowe/cbrewer2), GitHub. Retrieved March 29, 2025.
    addpath('/Users/deszoeks/Documents/MATLAB/user_tools/graphics/cbrewer2repo/cbrewer2')
    addpath('/Users/deszoeks/Documents/MATLAB/user_tools/graphics/colorspace') % a dependency of cbrewer2
    [tmap]=cbrewer2('div', 'BrBG', 12);
    [bw]=cbrewer2('seq', 'Greys', 12);
end

% Cite as:
% Charles (2021). cbrewer : colorbrewer schemes for Matlab 
% (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab), 
% MATLAB Central File Exchange. Retrieved December 2, 2021.

% load diverging_colormap.mat
%% Plots
% n = length(num_vector) + 1;
n = length(subset) + 1;
st = subset([1:end, end]);
figure;
    scale = [1 1 2 1.7]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
subplot(311)
    pcolor(1:n,t_comp2,dD_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['dD - dD_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'a', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(312)
    pcolor(1:n,t_comp2,Taf_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('T - T_{min} [\circC]')
    text(-3.3, 75, 'b', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-4.5 4.5])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(313)
    pcolor(1:n,t_comp2,qair_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('q - q_{tmin} [g/kg]')
    text(-3.3, 75, 'c', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-3 3])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',32)
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

%%
figure;  
    scale = [1 1 1 1.5]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
subplot(311)
    pcolor(1:n,t_comp2,dD_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['dD - dD_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'd', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(313)
    pcolor(1:n,t_comp2,DXS_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['DXS - DXS_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'e', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-2 2])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(312)
    pcolor(1:n,t_comp2,d18O_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['d^1^8O - d^1^8O_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'f', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-1 1])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off 
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',32)
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

%%
figure;
    scale = [1 1 1 1.5]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
ax1 = subplot(311);
    pcolor(1:n,t_comp2,d18O_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(ax1, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['d^{18}O - d^{18}O_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'd', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-1.2 1.2])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
ax2 = subplot(312);
    pcolor(1:n,t_comp2,wspd_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(ax2, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('U - U_{tmin} [m/s]')
    text(-3.3, 75, 'e', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off 
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
% set(findall(gcf,'-property','Fontsize'),'FontSize',32)
% set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])
% 
% figure;
%     scale = [1 1 1 1.5]; % 1.5 is the scaling of the y-dimension of the figure
%     cx = get(gcf,'Position');
%     set(gcf,'Position',scale.*cx);
%     cx = get(gcf,'PaperPosition');
%     set(gcf,'PaperPosition',scale.*cx);
ax3 = subplot(313)
    pcolor(1:n,t_comp2,prec_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(ax3, bw)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('rain rate [mm/hr]')
    text(-3.3, 75, 'f', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([2e-2 5])
    set(gca,'ColorScale','log')
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',32)
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

% 6-panel plot w/ d18O
figure;
clf;
    scale = [1 1 2 1.7]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
subplot(321)
    pcolor(1:n,t_comp2,dD_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['dD - dD_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'a', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
subplot(323)
    pcolor(1:n,t_comp2,Taf_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('T - T_{min} [\circC]')
    text(-3.3, 75, 'b', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-4.5 4.5])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
subplot(325)
    pcolor(1:n,t_comp2,qair_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('q - q_{tmin} [g/kg]')
    text(-3.3, 75, 'c', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-3 3])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
% right column
ax1 = subplot(322);
    pcolor(1:n,t_comp2,d18O_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(ax1, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title(['d^{18}O - d^{18}O_{tmin} [',char(8240),']'])
    text(-3.3, 75, 'd', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-1.2 1.2])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
ax2 = subplot(324);
    pcolor(1:n,t_comp2,wspd_comp2_sort(st,1:end-1)'); shading flat
    colorbar
    colormap(ax2, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title('U - U_{tmin} [m/s]')
    text(-3.3, 75, 'e', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off 
ax3 = subplot(326);
    pcolor(1:n,t_comp2,prec_comp2_sort(st,1:end-1)'); shading flat
    colorbar('YTick', [0.02 0.1 1 5])
    colormap(ax3, bw)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title('rain rate [mm/hr]')
    text(-3.3, 75, 'f', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([2e-2 5])
    set(gca,'ColorScale','log')
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off

% set figure settings
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',28)
set(findall(gcf,'-property','Fontweight'), 'Fontweight','normal')
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

saveas(gcf, 'rank14cp','epsc')
saveas(gcf, 'rank14cp','svg')
saveas(gcf, 'rank14cp','png')
saveas(gcf, 'rank14cp','pdf')
