% Temporal structure of individual cold pools ranked by their change in dD
load('workspace_17_cp_detection_algorithm_11min.mat');
% load('workspace_17_cp_detection_algorithm_1min.mat');
clearvars -except Taf qair wspd rh prec u v t_min_ind ...
                  t_max_ind t_end_ind t1min factor iso_cp_matrix cp_matrix ...
                  dD d18O q_iso DXS ship orange green blue mustard ...
                  Taf_comp2 t_comp2 dD_comp2 d18O_comp2 q_iso_comp2 DXS_comp2...
                  qair_comp2 wspd_comp2 rh_comp2 prec_comp2 u_comp2 v_comp2
Taf_comp2 = Taf_comp2([65,67:78,82:84,86],:);
wspd_comp2 = wspd_comp2([65,67:78,82:84,86],:);
qair_comp2 = qair_comp2([65,67:78,82:84,86],:);
prec_comp2 = prec_comp2([65,67:78,82:84,86],:);
rh_comp2 = rh_comp2([65,67:78,82:84,86],:);
u_comp2 = u_comp2([65,67:78,82:84,86],:);
v_comp2 = v_comp2([65,67:78,82:84,86],:);

d18O_comp2 = d18O_comp2([1,3:14,18:20,22],:);
dD_comp2 = dD_comp2([1,3:14,18:20,22],:);
DXS_comp2 = DXS_comp2([1,3:14,18:20,22],:);

%% Identifying the strongest and weakest cold pools based on dD
num_vector = [65,67:78,82:84,86]; % 65:86;

% SPD subset further to 14 valid cold pools based on those included later
% (not sure why 3 were excluded)
subset = [1:3, 5:8, 11:17]; % cuts to 14 events
num_vector = num_vector(subset);

delta_dD = dD(t_max_ind(num_vector)-factor) - dD(t_min_ind(num_vector)-factor);
cp_vector = [delta_dD, num_vector', (1:length(num_vector))'];
sorted_cp = sortrows(cp_vector,1,'descend');
cp_order = sorted_cp(:,3);
cp_order100 = sorted_cp(:,2);
Taf_comp2_sort = 9999*ones(size(Taf_comp2));
dD_comp2_sort = 9999*ones(size(dD_comp2));
wspd_comp2_sort = 9999*ones(size(wspd_comp2));
qair_comp2_sort = 9999*ones(size(qair_comp2));
prec_comp2_sort = 9999*ones(size(prec_comp2));
DXS_comp2_sort = 9999*ones(size(DXS_comp2));
d18O_comp2_sort = 9999*ones(size(d18O_comp2));

% Anomalies from values at t_max!!!
ncp = length(cp_order);
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

% %%
% addpath('/Users/estefania/Documents/Data/cbrewer')
% [tmap]=cbrewer('div', 'BrBG', 12);
% [bw]=cbrewer('seq', 'Greys', 12);
% % Cite as:
% % Charles (2021). cbrewer : colorbrewer schemes for Matlab 
% % (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab), 
% % MATLAB Central File Exchange. Retrieved December 2, 2021.
%
% -> Not available anymore, and Simon doesn't have it. Try
% Scott Lowe (2025). cbrewer2 (https://github.com/scottclowe/cbrewer2), GitHub. Retrieved March 29, 2025.
addpath('/Users/deszoeks/Documents/MATLAB/user_tools/graphics/cbrewer2repo/cbrewer2')

% load diverging_colormap.mat

%% Plots
n = length(num_vector) + 1;
figure;
    scale = [1 1 1 1.5]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
subplot(311)
    pcolor(1:n,t_comp2,dD_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title(['dD - dD_t_m_a_x [',char(8240),']'])
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(312)
    pcolor(1:n,t_comp2,Taf_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title('T - T_m_a_x [\circC]')
    ylim([-30 60])
    caxis([-3 3])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(313)
    pcolor(1:n,t_comp2,qair_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title('q - q_t_m_a_x [g/kg]')
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
    pcolor(1:n,t_comp2,dD_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title(['dD - dD_t_m_a_x [',char(8240),']'])
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(313)
    pcolor(1:n,t_comp2,DXS_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title(['DXS - DXS_t_m_a_x [',char(8240),']'])
    ylim([-30 60])
    caxis([-2 2])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(312)
    pcolor(1:n,t_comp2,d18O_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title(['d^1^8O - d^1^8O_t_m_a_x [',char(8240),']'])
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
subplot(311)
    pcolor(1:n,t_comp2,dD_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title(['dD - dD_t_m_a_x [',char(8240),']'])
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
subplot(312)
    pcolor(1:n,t_comp2,wspd_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(tmap)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title('U - U_t_m_a_x [m/s]')
    ylim([-30 60])
    caxis([-6 6])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off 
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',32)
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])

figure;
    scale = [1 1 1 1.5]; % 1.5 is the scaling of the y-dimension of the figure
    cx = get(gcf,'Position');
    set(gcf,'Position',scale.*cx);
    cx = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',scale.*cx);
subplot(313)
    pcolor(1:n,t_comp2,prec_comp2_sort(:,1:end-1)'); shading flat
    colorbar
    colormap(bw)
    hold on;
    plot([0 22],[0 0],':k','LineWidth',1.5)
    plot([0 22],[-19.7 -19.7],':k','LineWidth',1.5)
    plot([0 22],[31.5 31.5],':k','LineWidth',1.5)    
    ylabel('min rel. to cp')
    title('rain rate [mm/hr]')
    ylim([-30 60])
    caxis([0 1.5])
    set(gca,'TickDir','out'); % The only other option is 'in'
    xticks(1.5:1:18.5)
    xticklabels({'','2','','4','','6','','8','','10','','12','','14','','16','','18'})
    box off
set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
set(findall(gcf,'-property','Fontsize'),'FontSize',32)
set(findall(gcf,'-property','TickLength'),'TickLength',[.02,.1])
