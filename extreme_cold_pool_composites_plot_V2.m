% This code runs with output variables from
% cold_pool_detection_algorithm_recovery_modified.m
% It attempts to replicate Vogel's (2021) Fig. 3
% Last Modified June 16, 2022

% load('workspace_17_cp_detection_algorithm_11min.mat');
% load('workspace_17_cp_detection_algorithm_1min.mat');
% load('workspace_16_cp_detection_algorithm_11min.mat');
% load('workspace_16_cp_detection_algorithm_1min.mat');
load('workspace_cp_detection_algorithm_1min_w_surface_fluxes.mat');
 
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
figure;
subplot(521)
    hold on;
    right_yaxis_plotting(t_comp2,Taf_comp2,Taf_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    ylabel('T [K]')
    % title('17 cold pools composite')
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([23.4 26.5])
    legend('all','30% strongest','30% weakest','Location','southwest','box','off')
subplot(523)
    % plot(t_comp2,nanmean(q_iso_comp2),'-sk')
    hold on;
    right_yaxis_plotting(t_comp2,qair_comp2,qair_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    ylabel('q [gkg^-^1]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([13.6 15.7])
%     ylim([-.8 .8]) % for 11 min file
%     ylim([-.9 1.1])% for 1 min file    
    % legend('Picarro','PSD')
subplot(525)
    hold on;
    right_yaxis_plotting(t_comp2,rh_comp2,rh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    ylabel('RH [%]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([65 85])
subplot(522)
    hold on;
    right_yaxis_plotting(t_comp2,wspd_comp2,wspd_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    ylabel('U [ms^-^1]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([7 11.5])
subplot(524) % no anomalies for RR, rain rate
    % yyaxis right
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
    plot([-18.7 -18.7],[0 4],':k','LineWidth',.5)
%     plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[0 4],':k','LineWidth',.5)
    plot([30.4 30.4],[0 4],':k','LineWidth',.5)
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([0 4])
    % yyaxis left
    % set(gca,'yticklabels',[])
    ylabel('RR [mmh^-^1]')

% Surface Fluxes %    
subplot(527)
    hold on;
    right_yaxis_plotting(t_comp2,sh_comp2,sh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    ylabel('SH [Wm^-^2]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-33 -2])
subplot(529)
    hold on;
    right_yaxis_plotting(t_comp2,lh_comp2,lh_comp2_17,mask_strong100,mask_weak100,n,Tmax_i)
    ylabel('LH [Wm^-^2]')
    box off
    xlabel('time [minutes]')
    xlim([-60 60])
    ylim([-252 -160])
% Isotopes %    
subplot(526)
    hold on;
    right_yaxis_plotting(t_comp2,dD_comp2,dD_comp2,mask_strong,mask_weak,n,Tmax_i)
    ylabel(['\deltaD [',char(8240),']'])
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-73.5 -69])
%     ylim([-2 3])  % for 11 min file
%     ylim([-2 3.4])% for 1 min file 
subplot(528)
    hold on;
    right_yaxis_plotting(t_comp2,d18O_comp2,d18O_comp2,mask_strong,mask_weak,n,Tmax_i)
    ylabel(['\delta^1^8O [',char(8240),']'])
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-10.65 -9.95])
    box off
subplot(5,2,10)
    hold on;
    right_yaxis_plotting(t_comp2,DXS_comp2,DXS_comp2,mask_strong,mask_weak,n,Tmax_i)
    ylabel(['DXS [',char(8240),']'])
    xlabel('time [minutes]')
    xlim([-60 60])
    ylim([7.2 8.5])
    box off
set(findall(gcf,'-property','Fontsize'),'FontSize',10)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])
