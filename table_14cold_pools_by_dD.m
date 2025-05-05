% Temporal structure of individual cold pools ranked by their change in dD
load('workspace_17_cp_detection_algorithm_11min.mat');
% load('workspace_17_cp_detection_algorithm_1min.mat');
clearvars -except Taf qair wspd rh prec u v t_min_ind ...
                  t_max_ind t_end_ind t1min factor iso_cp_matrix cp_matrix ...
                  dD d18O q_iso DXS ship orange green blue mustard ...
                  Taf_comp2 t_comp2 dD_comp2 d18O_comp2 q_iso_comp2 DXS_comp2...
                  qair_comp2 wspd_comp2 rh_comp2 prec_comp2 u_comp2 v_comp2
% there's no prec variable
prate = ncread('data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc', 'prate'); % mm/h
%[sum(prate)/60] = mm

i17 = [65,67:78,82:84,86]; % chron order
Taf_comp2 = Taf_comp2(i17,:);
wspd_comp2 = wspd_comp2(i17,:);
qair_comp2 = qair_comp2(i17,:);
prec_comp2 = prec_comp2(i17,:);
rh_comp2 = rh_comp2(i17,:);
u_comp2 = u_comp2(i17,:);
v_comp2 = v_comp2(i17,:);

j17 = [1,3:14,18:20,22]; % chron order
d18O_comp2 = d18O_comp2(j17,:);
dD_comp2 = dD_comp2(j17,:);
DXS_comp2 = DXS_comp2(j17,:);

subset = [1:3, 5:8, 11:17]; % cuts to 14 events in dD order

%% sort cold pools by dD
num_vector = i17;
% old way:
% delta_dD = dD(t_max_ind(num_vector)-factor) - dD(t_min_ind(num_vector)-factor); 
% new way of ranking by max(dD)
mxdD = NaN(size(num_vector));
for i = 1:length(num_vector)
    nv = num_vector(i);
    mxdD(i) = max(dD((t_max_ind(nv):t_end_ind(nv))-factor)); % chron order
end
delta_dD = dD(t_max_ind(num_vector)-factor) - mxdD'; % chron order
% try ranking from dD sort of 14 events in excel
cp_vector = [delta_dD, num_vector', (1:length(num_vector))'];
sorted_cp = sortrows(cp_vector,1,'descend');
delta_dD_sort = sorted_cp(:,1); % sorted delta_dD
cp_order = sorted_cp(:,3);    % sorts chron-ordered data by dD
cp_order100 = sorted_cp(:,2); % num_vector=i17 sorted by dD
% cp_order100 = i17(cp_order)
cp_order100_subset = i17(cp_order(subset));


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
    DXS_comp2_sort(k,:) = DXS_comp2(cp_order(k),:) -DXS( t_max_ind(cp_order100(k))-factor);
    dD_comp2_sort(k,:)  = dD_comp2(cp_order(k),:)  -dD(  t_max_ind(cp_order100(k))-factor);
    d18O_comp2_sort(k,:)= d18O_comp2(cp_order(k),:)-d18O(t_max_ind(cp_order100(k))-factor);
end

% plot max differences for each cold pool
% cold pool difference Dcp differences any variable T  T(t_min) - T(t_max)
DcpT = @(T) T(t_min_ind(cp_order100)) - T(t_max_ind(cp_order100));
% Dcpmx = @(q) max(q(t_min_ind(cp_order100)-10:t_min_ind(cp_order100)+10))- max(q(t_max_ind(cp_order100)-10:t_max_ind(cp_order100)+10)); % nope
DiffT = -DcpT(Taf);
Diffq =  DcpT(qair);

q0 = qair(t_max_ind(cp_order100)); % -> 17 cold pools
dD0 = dD(t_max_ind(cp_order100)-factor);
% prec0 = prec(t_max_ind(cp_order100)); % no prec

DiffdD = dD(t_min_ind(cp_order100)-factor) - dD0;

% find max quantities in each cold pool
qmax = NaN(length(cp_order100),1);
pmax = NaN(length(cp_order100),1);
psum = NaN(length(cp_order100),1);
dDmax = NaN(length(cp_order100),1);
for i = 1:length(cp_order100)
    cpi = cp_order100(i);
    qmax(i) = max( qair(t_max_ind(cpi):t_min_ind(cpi)+10) );
    pmax(i) = max( prec_comp2_sort(i,1:80), [], 2);
    % psum(i) = sum( prec_comp2_sort(i, -10<t_comp2<30) ); % composite time index not suitable
    psum(i) = sum(prate(t_max_ind(cpi):t_min_ind(cpi)),'omitnan')/60; % -> mm
    dDmax(i) = max( dD_comp2_sort(i,1:80), [], 2);
end

clf()
subplot(3,4,1)
% plot(DiffT(subset), -delta_dD_sort(subset), '.', 'markersize',10, 'linestyle','none')
% hold on
plot(DiffT(subset), DiffdD(subset), '.', 'markersize',10, 'linestyle','none')
% plot(DiffT(subset), dDmax(subset), '.', 'markersize',10, 'linestyle','none') % composite but same as above
r = corrcoef(DiffT(subset), DiffdD(subset));
text(0.2, 5.5, sprintf('R^2=%0.2f',r(2,1)^2));
xlim([0,4])
ylim([0,7])
ylabel('{\delta}D(t_{min})-{\delta}D(t_0) [10^{-3}]') %, {\delta}D_{max}-{\delta}D(t_0)')
xlabel('-(T_{min} - T_0) [°C]')
text(-2.4, 7, 'a')
axis square

subplot(3,4,5)
plot(DiffT(subset), Diffq(subset), '.', 'markersize',10, 'linestyle','none')
% hold on
% plot(DiffT(subset), qmax(subset)-q0(subset)', '.', 'markersize',10, 'linestyle','none')
r = corrcoef(DiffT(subset), Diffq(subset));
text(0.2, 1.2, sprintf('R^2=%0.2f',r(2,1)^2));
xlim([0,4])
ylim([0,2])
ylabel('q(t_{min})-q(t_0) [g/kg]') %, q_{max}-q(t_0)')
xlabel('-(T_{min} - T_0) [°C]')
text(-2.4, 2, 'b')
axis square

subplot(3,4,9) % no relationship in scatter plot
% plot(DiffT(subset), pmax(subset), '.', 'markersize',10, 'linestyle','none')
semilogy(DiffT(subset), max(0.01,psum(subset)), '.', 'markersize',10, 'linestyle','none')
r = corrcoef(DiffT(subset), psum(subset));
text(1.6, 2e-2, sprintf('R^2=%0.2f',r(2,1)^2));
xlim([0,4])
ylim([0.01,1])
ylabel('rain [mm]')
xlabel('-(T_{min} - T_0) [°C]')
text(-2.4, 1, 'c')
axis square

saveas(gcf, 'table1_fig.eps')
saveas(gcf, 'table1_fig.png')
saveas(gcf, 'table1_fig.pdf')
saveas(gcf, 'table1_fig.svg')

% effectively reads max, cumulatives off Fig 7
maxdD   = max(  dD_comp2_sort(subset, -25<t_comp2<45), [], 2);
cumprec = sum(prec_comp2_sort(subset, -25<t_comp2<45), [], 2);


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
    addpath('/Users/estefania/Documents/Data/cbrewer') % ignore warnings
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

%% Plots
% n = length(num_vector) + 1;
n = length(subset) + 1;
st = subset([1:end, end]); % already used max(dD) above
% st = subset(rankdD_cronOrder([1:end, end]));

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
    colorbar('TickDirection', 'out');
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title(['dD - dD_{t0} [',char(8240),']'])
    text(-3.3, 75, 'a', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
subplot(323)
    pcolor(1:n,t_comp2,Taf_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('T - T_{0} [\circC]')
    text(-3.3, 75, 'b', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-4.5 4.5])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
subplot(325)
    pcolor(1:n,t_comp2,qair_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('time')
    title('q - q_{t0} [g/kg]')
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
    colorbar('TickDirection', 'out');
    colormap(ax1, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title(['d^{18}O - d^{18}O_{t0} [',char(8240),']'])
    text(-3.3, 75, 'd', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-1.2 1.2])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off
ax2 = subplot(324);
    pcolor(1:n,t_comp2,wspd_comp2_sort(st,1:end-1)'); shading flat
    colorbar('TickDirection', 'out');
    colormap(ax2, tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel(' ')
    title('U - U_{t0} [m/s]')
    text(-3.3, 75, 'e', 'fontsize',32, 'clipping','off')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','13','','11','','9','','7','','5','','3','','1'})
    box off 
ax3 = subplot(326);
    pcolor(1:n,t_comp2,prec_comp2_sort(st,1:end-1)'); shading flat
    cb = colorbar('YTick', [0.02 0.1 1 5], 'TickDirection', 'out');
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
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
set(findall(gcf,'-property','Fontweight'), 'Fontweight','normal')
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

orient landscape
% saveas(gcf, 'rank14cp','epsc')
% saveas(gcf, 'rank14cp','svg')
% saveas(gcf, 'rank14cp','png')
% saveas(gcf, 'rank14cp','pdf')

% saveas(gcf, 'rank14cp2','epsc')
% saveas(gcf, 'rank14cp2','svg')
% saveas(gcf, 'rank14cp2','png')
% saveas(gcf, 'rank14cp2','pdf')
