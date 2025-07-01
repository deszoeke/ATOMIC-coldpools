% Copyright 2025 Simon P. de Szoeke and Estefanía Quiñones 
% Meléndez.
% 
% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation 
% files (the “Software”), to deal in the Software without 
% restriction, including without limitation the rights to use, 
% copy, modify, merge, publish, distribute, sublicense, and/or 
% sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following 
% conditions:
%
% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.
% Temporal structure of individual cold pools ranked by their change in dD
load('data/workspace_17_cp_detection_algorithm_11min.mat');
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
n = length(subset);

%% Identifying the strongest and weakest cold pools based on dD
num_vector = i17;
% old way:
% delta_dD = dD(t_max_ind(num_vector)-factor) - dD(t_min_ind(num_vector)-factor); 
% new way of ranking by max(dD)
mxdD = NaN(size(num_vector));
for i = 1:length(num_vector)
    nv = num_vector(i);
    mxdD(i) = max(dD((t_max_ind(nv):t_end_ind(nv))-factor));
end
delta_dD = dD(t_max_ind(num_vector)-factor) - mxdD';
% try ranking from dD sort of 14 events in excel
% rankdD_cronOrder = [7, 13, 10, 5, 3, 4, 1, 6, 8, 12, 14, 9, 11, 2];
cp_vector = [delta_dD, num_vector', (1:length(num_vector))'];
sorted_cp = sortrows(cp_vector,1,'descend');
cp_order = sorted_cp(:,3);
cp_order100 = sorted_cp(:,2);
% cp_order = rankdD_cronOrder;
% i14 = i17(subset);
% cp_order100 = i14(rankdD_cronOrder);

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
addpath(genpath("deps"))
[tmap]=cbrewer2('div', 'BrBG', 12);
[bw]=cbrewer2('seq', 'Greys', 12);

%% Plots

% 6-panel plot w/ d18O
figure;
st = subset([1:end, end]);
pnx = -2;
clf;
    scale = [1 1 2 1.7]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
subplot(321)
    pcolor(1:n+1,t_comp2,dD_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['{\delta}D - {\delta}dD_{t0} [',char(8240),']'])
    text(pnx, 75, 'a', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
subplot(323)
    pcolor(1:n+1,t_comp2,Taf_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('T - T_{0} [\circC]')
    text(pnx, 75, 'b', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-4.5 4.5])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
subplot(325)
    pcolor(1:n+1,t_comp2,qair_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('q - q_{t0} [g/kg]')
    text(pnx, 75, 'c', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-3 3])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
% right column
ax1 = subplot(322);
    pcolor(1:n+1,t_comp2,d18O_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(ax1, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title(['{\delta}^{18}O - {\delta}^{18}O_{t0} [',char(8240),']'])
    text(pnx, 75, 'd', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-1.2 1.2])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
ax2 = subplot(324);
    pcolor(1:n+1,t_comp2,wspd_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(ax2, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title('U - U_{t0} [m/s]')
    text(pnx, 75, 'e', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off 
ax3 = subplot(326);
    pcolor(1:n+1,t_comp2,prec_comp2_sort(st,1:end-1)'); shading flat
    cb = colorbar('YTick', [0.02 0.1 1 5], 'TickDirection', 'out');
    colormap(ax3, bw)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title('rain rate [mm/hr]')
    text(pnx, 75, 'f', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    clim([2e-2 5])
    set(gca,'ColorScale','log')
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off

% set figure settings
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
set(findall(gcf,'-property','Fontweight'), 'Fontweight','normal')
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

orient landscape

% saveas(gcf, 'rank14cp2','epsc')
% saveas(gcf, 'rank14cp2','svg')
% saveas(gcf, 'rank14cp2','png')
% saveas(gcf, 'rank14cp2','pdf')


%% plot comp_mean
%{
st = [ 1     2     3     5     6     7     8    11    12    13    14    15    16    17 ];
colors = [[0.8, 0.5, 0]; [0.1, 0.6, 0.2]; [0.1, 0.8, 0.6]; [0.1, 0.6, 1]; [0, 0, 0]; 0.6+[0, 0, 0]];
clf()
ch = 'a';
ax = vert_axes_stack(6);
plot_time_comp(ax(1), t_comp2, Taf_comp2(st,:), [0.8, 0.5, 0])
ylabel(ax(1), 'T (°C)')
text(ax(1), -55, 25, string(char(ch)), 'FontSize',13); ch=ch+1;
plot_time_comp(ax(2), t_comp2, qair_comp2(st,:), [0.1, 0.6, 0.2])
ylabel(ax(2), 'q (g/kg)')
text(ax(2), -55, 15.0, string(char(ch)), 'FontSize',13); ch=ch+1;
plot_time_comp(ax(3), t_comp2, dD_comp2(st,:), [0.1, 0.8, 0.6])
ylabel(ax(3), '\deltaD (‰)')
text(ax(3), -55, -67, string(char(ch)), 'FontSize',13); ch=ch+1;
plot_time_comp(ax(4), t_comp2, d18O_comp2(st,:), [0.1, 0.6, 1])
ylabel(ax(4), '\delta^{18}O (‰)')
text(ax(4), -55, -9.7, string(char(ch)), 'FontSize',13); ch=ch+1;
plot_time_comp(ax(5), t_comp2, DXS_comp2(st,:), 'k')
ylabel(ax(5), 'DXS (‰)')
text(ax(5), -55, 13, string(char(ch)), 'FontSize',13); ch=ch+1;
plot_time_comp(ax(6), t_comp2, prec_comp2(st,:), 0.6.*[1 1 1])
ylabel(ax(6), 'rain (mm/h)')
text(ax(6), -55, 1.1, string(char(ch)), 'FontSize',13);
set(ax(1:5), 'xticklabel', [])

% add strongest 2 events
% 07-Feb-2020 12:29:00
plot(ax(1),t_comp2,  Taf_comp2( 8,:), '--', 'color',[    0.8000    0.5000         0])
plot(ax(2),t_comp2, qair_comp2( 8,:), '--', 'color',[    0.1000    0.6000    0.2000])
plot(ax(3),t_comp2,   dD_comp2( 8,:), '--', 'color',[    0.1000    0.8000    0.6000])
plot(ax(4),t_comp2, d18O_comp2( 8,:), '--', 'color',[    0.1000    0.6000    1.0000])
plot(ax(5),t_comp2,  DXS_comp2( 8,:), '--', 'color',[         0         0         0])
plot(ax(6),t_comp2, max(1e-4,prec_comp2( 8,:)), '--', 'color',[    0.6000    0.6000    0.6000])
% 10-Feb-2020 15:59:00
plot(ax(1),t_comp2,  Taf_comp2(17,:), '-.', 'color',[    0.8000    0.5000         0])
plot(ax(2),t_comp2, qair_comp2(17,:), '-.', 'color',[    0.1000    0.6000    0.2000])
plot(ax(3),t_comp2,   dD_comp2(17,:), '-.', 'color',[    0.1000    0.8000    0.6000])
plot(ax(4),t_comp2, d18O_comp2(17,:), '-.', 'color',[    0.1000    0.6000    1.0000])
plot(ax(5),t_comp2,  DXS_comp2(17,:), '-.', 'color',[         0         0         0])
plot(ax(6),t_comp2, max(1e-4,prec_comp2(17,:)), '-.', 'color',[    0.6000    0.6000    0.6000])

set(ax(6), 'yscale','log', 'ytick',10.^(-3:1), ...
    'yticklabel',["","10^{-2}","","10^0",""], ...
    'ylim',[0.99e-3, 1.01e1])
xlabel(ax(6), 'composite time (minute)')
% set(ax(1), 'XTickLabel', ["","","t_0","t_{min}","t_{end}","",""], 'XAxisLocation','top')
for i=1:6; ax(i).FontSize = 13; end

% composite time gridlines
for i=1:length(ax)
    plot(ax(i), -19.7+[0 0], ax(i).YLim, 'k-')
    plot(ax(i),       [0 0], ax(i).YLim, 'k-')
    plot(ax(i),  31.5+[0 0], ax(i).YLim, 'k-')
end
text(ax(1), -19.7, 27, "t_0", 'horizontalalignment','center', 'fontsize',13, 'clipping','off')
text(ax(1),    0 , 27, "t_{min}", 'horizontalalignment','center', 'fontsize',13, 'clipping','off')
text(ax(1),  31.5, 27, "t_{end}", 'horizontalalignment','center', 'fontsize',13, 'clipping','off')
%}

% fmt = ["epsc"; "svg"; "png"; "pdf"];
% for i=1:length(fmt)
%     saveas(gcf, 'composite_6panel_1col',fmt(i))
% end
