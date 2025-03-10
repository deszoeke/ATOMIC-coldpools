% This code runs with output variables from
% cold_pool_detection_algorithm_recovery_modified.m
% It attempts to replicate Vogel's (2021) Fig. 3
% Last Modified June 16, 2022

% load('workspace_17_cp_detection_algorithm_11min.mat');
% load('workspace_17_cp_detection_algorithm_1min.mat');
% load('workspace_16_cp_detection_algorithm_11min.mat');
% load('workspace_16_cp_detection_algorithm_1min.mat');
load('workspace_cp_detection_algorithm_11min_w_surface_fluxes.mat');
 
clearvars -except Taf qair wspd rh prec u v t_min_ind ...
                  t_max_ind t_end_ind t1min factor iso_cp_matrix cp_matrix ...
                  dD d18O q_iso DXS ship ... % orange green blue mustard ...
                  Taf_comp2 Ta_comp2 t_comp2 dD_comp2 d18O_comp2 q_iso_comp2 DXS_comp2...
                  qair_comp2 wspd_comp2 rh_comp2 prec_comp2 u_comp2 v_comp2 ...
                  sh_comp2 lh_comp2
                  % qs_comp2
                  
% Taf_comp2 = Ta_comp2;

load('mask_extreme_30%_cold_pools.mat')
% load('mask_extreme_15%_cold_pools.mat')
blue = [0 0.4470 0.7410];

%% Plotting the 17 cold pools that overlap with isotopic data
Taf_comp2_17 = Taf_comp2([65,67,69:78,82:84,86],:);
wspd_comp2_17 = wspd_comp2([65,67,69:78,82:84,86],:);
qair_comp2_17 = qair_comp2([65,67,69:78,82:84,86],:);
prec_comp2_17 = prec_comp2([65,67,69:78,82:84,86],:);
rh_comp2_17 = rh_comp2([65,67,69:78,82:84,86],:);
u_comp2_17 = u_comp2([65,67,69:78,82:84,86],:);
v_comp2_17 = v_comp2([65,67,69:78,82:84,86],:);
sh_comp2_17 = sh_comp2([65,67,69:78,82:84,86],:);
lh_comp2_17 = lh_comp2([65,67,69:78,82:84,86],:);

% qs_comp2_17 = qs_comp2([65,67,69:78,82:84,86],:); % CHECK UNITS!!!

Tmax_i = 41; % -18.7 @ 41; t_comp2(41) = 18.7;
             % 30.4 @ 81; t_comp2(81) = 30.4;

%%
subplot(521)
    hold on;
    mask_lines_plotting(t_comp2,Taf_comp2,Taf_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    ylabel('T-T_m_a_x [K]')
    % title('17 cold pools composite')
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-2.5 0.4])
    legend('all','30% strongest','30% weakest','Location','southwest','box','off')
subplot(523)
    % plot(t_comp2,nanmean(q_iso_comp2),'-sk')
    hold on;
    mask_lines_plotting(t_comp2,qair_comp2,qair_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('q-q_t_m_a_x [gkg^-^1]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-.8 .8]) % for 11 min file
    ylim([-.9 1.1])% for 1 min file    
    % legend('Picarro','PSD')
subplot(525)
    hold on;
    mask_lines_plotting(t_comp2,rh_comp2,rh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('RH-RH_t_m_a_x [%]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-5 15])
subplot(522)
    hold on;
    mask_lines_plotting(t_comp2,wspd_comp2,wspd_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('U-U_t_m_a_x [ms^-^1]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-2.5 1.5])
subplot(524) % no anomalies for RR, rain rate
    plot(t_comp2,nanmean(prec_comp2_17),'-k','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(prec_comp2_17)+std(prec_comp2_17,'omitnan')/sqrt(n),'.k','MarkerSize',.2)
    plot(t_comp2,nanmean(prec_comp2_17)-std(prec_comp2_17,'omitnan')/sqrt(n),'.k','MarkerSize',.2)    
    % Strong cold pools
    plot(t_comp2,nanmean(prec_comp2(mask_strong100,:)),'-r','LineWidth',2)
    plot(t_comp2,nanmean(prec_comp2(mask_strong100,:))+std(prec_comp2(mask_strong100,:),'omitnan')/sqrt(n),'.r','MarkerSize',.2)
    plot(t_comp2,nanmean(prec_comp2(mask_strong100,:))-std(prec_comp2(mask_strong100,:),'omitnan')/sqrt(n),'.r','MarkerSize',.2)    
    % Weak cold pools
    plot(t_comp2,nanmean(prec_comp2(mask_weak100,:)),'-','LineWidth',2,'Color',blue)
    plot(t_comp2,nanmean(prec_comp2(mask_weak100,:))+std(prec_comp2(mask_weak100,:),'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)
    plot(t_comp2,nanmean(prec_comp2(mask_weak100,:))-std(prec_comp2(mask_weak100,:),'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)    
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    box off
    ylabel('RR [mmh^-^1]')
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([0 3.5])
% Surface Fluxes %    
subplot(527)
    hold on;
    mask_lines_plotting(t_comp2,sh_comp2,sh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('SH-SH_t_m_a_x [Wm^-^2]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-24 11])
    right_yaxis_plotting(t_comp2,sh_comp2,sh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
subplot(529)
    hold on;
    mask_lines_plotting(t_comp2,lh_comp2,lh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('LH-LH_t_m_a_x [Wm^-^2]')
    box off
    xlabel('time [minutes]')
    xlim([-60 60])
    ylim([-61 54])    
% Isotopes %    
subplot(526)
    hold on;
    mask_lines_plotting(t_comp2,dD_comp2,dD_comp2,mask_strong,mask_weak,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel(['\deltaD-\deltaD_t_m_a_x [',char(8240),']'])
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-2 3])  % for 11 min file
    ylim([-2 3.4])% for 1 min file 
subplot(528)
    hold on;
    mask_lines_plotting(t_comp2,d18O_comp2,d18O_comp2,mask_strong,mask_weak,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel(['\delta^1^8O-\delta^1^8O_t_m_a_x [',char(8240),']'])
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-.4 .5])
    box off
subplot(5,2,10)
    hold on;
    mask_lines_plotting(t_comp2,DXS_comp2,DXS_comp2,mask_strong,mask_weak,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel(['DXS-DXS_t_m_a_x [',char(8240),']'])
    xlabel('time [minutes]')
    xlim([-60 60])
    ylim([-0.8 1.2])
    box off
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])

%% Comparing q_air to qs, both from PSL data set to q_iso
subplot(421)
    hold on;
    mask_lines_plotting(t_comp2,qair_comp2,qair_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    % plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('q-q_t_m_a_x [gkg^-^1]')
    box off
    xlim([-60 60])
    ylim([-.9 1.1])% for 1 min file    
    xlabel('time [minutes]')
subplot(422)
    hold on;
    mask_lines_plotting(t_comp2,qs_comp2/1000,qs_comp2_17/1000,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    % plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('qs-qs_t_m_a_x [gkg^-^1]')
    box off
    xlim([-60 60])
    ylim([-.9 1.1])% for 1 min file    
    xlabel('time [minutes]')
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])

%% Comparing q_air from PSL sensor to q_iso
subplot(421)
    hold on;
    mask_lines_plotting(t_comp2,qair_comp2,qair_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    % plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('q-q_t_m_a_x [gkg^-^1]')
    box off
    xlim([-60 60])
    ylim([-.9 1.1])% for 1 min file    
    title('PSL sensor')
    xlabel('time [minutes]')
subplot(423)
    hold on;
    mask_lines_plotting(t_comp2,q_iso_comp2,q_iso_comp2_17,mask_strong,mask_weak,n,Tmax_i)
    plot([-18.7 -18.7],[-100 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-100 100],':k','LineWidth',.5)
    % plot([30.4 30.4],[-100 100],':k','LineWidth',.5)
    ylabel('qs-qs_t_m_a_x [gkg^-^1]')
    box off
    xlim([-60 60])
    ylim([-.9 1.1])% for 1 min file    
    title('Picarro')
    xlabel('time [minutes]')
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])
