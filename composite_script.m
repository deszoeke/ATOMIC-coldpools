%% Using entire RHB air temp timeseries
% Loading 1-min data
load '1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' qair t1min Ta rr sst slp rh
load 'iso_data_1min_intervals_FLAGGED_w_runningmean.mat' d18O dD DXS iso_time
ind = find(t1min>=iso_time(1));
ind = ind(1);
Ta   = Ta(ind:ind+length(iso_time)-1);   % Ta filtered   => 11-min running average
qair = qair(ind:ind+length(iso_time)-1); % qair filtered => 11-min running average
sst  = sst(ind:ind+length(iso_time)-1);  % sst filtered  => 11-min running average
slp  = slp(ind:ind+length(iso_time)-1);  % slp filtered  => 11-min running average
rr   = rr(ind:ind+length(iso_time)-1);
rh   = rh(ind:ind+length(iso_time)-1);

%% Uploading data from COLD POOLS detection algorithm %%
load cold_pool_detection_workspace_17.mat
cp_matrix = cp_matrix(1:15,:); % eliminating last row
load cold_pool_flag_1min_UPDATED_for_iso_times.mat

%% FLAG for cold pool times (indexes)
t_max0_matrix = 9999*ones(length(t_max),length(cp_matrix)); % time
t_min0_matrix = 9999*ones(length(t_max),length(cp_matrix)); % time

for k = 1:length(t_max)
    ii = t_max_ind(k):t_end_ind(k);
    t_max0_matrix(k,1:length(ii)) = (t1min(ii)-t_max(k))./(t_min(k)-t_max(k));
    t_min0_matrix(k,1:length(ii)) = (t1min(ii)-t_min(k))./(t_end(k)-t_min(k));
end
t_max0_matrix(t_max0_matrix==9999) = NaN;
t_min0_matrix(t_min0_matrix==9999) = NaN;
%% Interpolating to obtain a single normalized time vector
cp_total = size(cp_matrix,1);
% For t_min @ 0 and t_end @ 1
t_norm    = 0:0.05:1;
% T_norm  = 0:0.10:3; % Simon says we don't need to normalize time!!!
Taf_norm  = 9999*ones(cp_total,length(t_norm));
Ta_norm   = 9999*ones(cp_total,length(t_norm));
prec_norm = 9999*ones(cp_total,length(t_norm));
qair_norm = 9999*ones(cp_total,length(t_norm));
rh_norm   = 9999*ones(cp_total,length(t_norm));
% wdir_norm = 9999*ones(cp_total,length(t_norm));
% wspd_norm = 9999*ones(cp_total,length(t_norm));
% u_norm    = 9999*ones(cp_total,length(t_norm));
% v_norm    = 9999*ones(cp_total,length(t_norm));
% sh_norm   = 9999*ones(cp_total,length(t_norm));
% lh_norm   = 9999*ones(cp_total,length(t_norm));
% sst_norm   = 9999*ones(cp_total,length(t_norm));
% slp_norm   = 9999*ones(cp_total,length(t_norm));
dD_norm   = 9999*ones(cp_total,length(t_norm));
d18O_norm = 9999*ones(cp_total,length(t_norm));
DXS_norm  = 9999*ones(cp_total,length(t_norm));
% q_iso_norm = 9999*ones(cp_total,length(t_norm));
% mr_norm    = 9999*ones(cp_total,length(t_norm));

for k = 1:cp_total % the last cold pool #17 is empty
    Taf_norm(k,:)  = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),Taf(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    Ta_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),Ta(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    prec_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),rr(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    qair_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),qair(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    rh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),rh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    dD_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),dD(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    d18O_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),d18O(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    DXS_norm(k,:)  = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),DXS(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     wdir_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),wdir(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     wspd_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),wspd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     u_norm(k,:)    = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),u(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     v_norm(k,:)    = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),v(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     sh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),sh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     lh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),lh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
end

%% For t_max @ 0 and t_min @ 1
Taf_norm2  = zeros(cp_total,length(t_norm));
Ta_norm2   = zeros(cp_total,length(t_norm));
prec_norm2 = zeros(cp_total,length(t_norm));
qair_norm2 = zeros(cp_total,length(t_norm));
rh_norm2   = zeros(cp_total,length(t_norm));
% wdir_norm2 = zeros(cp_total,length(t_norm));
% wspd_norm2 = zeros(cp_total,length(t_norm));
% u_norm2    = zeros(cp_total,length(t_norm));
% v_norm2    = zeros(cp_total,length(t_norm));
% sh_norm2   = zeros(cp_total,length(t_norm));
% lh_norm2   = zeros(cp_total,length(t_norm));
dD_norm2   = zeros(cp_total,length(t_norm));
d18O_norm2 = zeros(cp_total,length(t_norm));
DXS_norm2  = zeros(cp_total,length(t_norm));
% q_iso_norm2 = zeros(cp_total,length(t_norm));
% mr_norm2    = zeros(cp_total,length(t_norm));

for k = 1:cp_total % the last cold pool #17 is empty
    Taf_norm2(k,:)  = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),Taf(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    Ta_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),Ta(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    prec_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),rr(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    qair_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),qair(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    rh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),rh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     wdir_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),wdir(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     wspd_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),wspd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     u_norm2(k,:)    = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),u(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     v_norm2(k,:)    = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),v(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     sh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),sh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
%     lh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),lh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    dD_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),dD(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    d18O_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),d18O(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    DXS_norm2(k,:)  = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),DXS(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
end
%% Creating composite time series in normalized time
Taf_comp  = [Taf_norm2 Taf_norm(:,2:end)];
Ta_comp   = [Ta_norm2 Ta_norm(:,2:end)];
prec_comp = [prec_norm2 prec_norm(:,2:end)];
qair_comp = [qair_norm2 qair_norm(:,2:end)];
rh_comp   = [rh_norm2 rh_norm(:,2:end)];
% wdir_comp = [wdir_norm2 wdir_norm(:,2:end)];
% wspd_comp = [wspd_norm2 wspd_norm(:,2:end)];
% u_comp = [u_norm2 u_norm(:,2:end)];
% v_comp = [v_norm2 v_norm(:,2:end)];
% sh_comp   = [sh_norm2 sh_norm(:,2:end)];
% lh_comp   = [lh_norm2 lh_norm(:,2:end)];
dD_comp   = [dD_norm2 dD_norm(:,2:end)];
d18O_comp = [d18O_norm2 d18O_norm(:,2:end)];
DXS_comp  = [DXS_norm2 DXS_norm(:,2:end)];
% q_iso_comp = [q_iso_norm2 q_iso_norm(:,2:end)];
% mr_comp = [mr_norm2 mr_norm(:,2:end)];

% t_norm2 = linspace(-1,0,21);
% t_norm = linspace(0,1,21); 
% t_comp = [t_norm2 t_norm(2:end)];
% t_comp = [-20:1:0,1:1.5:28,30];
t_comp = [linspace(-30.4,0,21) linspace(1,24.9,20)];

% mean(t_min([65,67:68,70:86])-t_max([65,67:68,70:86]))*24*60; % 19.6667 min
mean(t_min-t_max)*24*60; % Previously => 18.6931 min; NOW => 30.4286 min
% mean(t_end([65,67:68,70:86])-t_min([65,67:68,70:86]))*24*60; % 25.7143 min
mean(t_end-t_min)*24*60; % Previously => 30.3762 min; NOW => 24.9286 min
%% Incorporating background data for 40 min before and after t_min
bgwindow = 40; % min
bg_Taf   = zeros(cp_total,bgwindow+1);
bg_Taf2  = zeros(cp_total,bgwindow-10+1);
bg_Ta    = zeros(cp_total,bgwindow+1);
bg_Ta2   = zeros(cp_total,bgwindow-10+1);
bg_rh    = zeros(cp_total,bgwindow+1);
bg_rh2   = zeros(cp_total,bgwindow-10+1);
bg_prec  = zeros(cp_total,bgwindow+1);
bg_prec2 = zeros(cp_total,bgwindow-10+1);
bg_qair  = zeros(cp_total,bgwindow+1);
bg_qair2 = zeros(cp_total,bgwindow-10+1);
% bg_wspd  = zeros(cp_total,bgwindow+1);
% bg_wspd2 = zeros(cp_total,bgwindow-10+1);
% bg_wdir  = zeros(cp_total,bgwindow+1);
% bg_wdir2 = zeros(cp_total,bgwindow-10+1);
% bg_u     = zeros(cp_total,bgwindow+1);
% bg_u2    = zeros(cp_total,bgwindow-10+1);
% bg_v     = zeros(cp_total,bgwindow+1);
% bg_v2    = zeros(cp_total,bgwindow-10+1);
% bg_sh    = zeros(cp_total,bgwindow+1);
% bg_sh2   = zeros(cp_total,bgwindow-10+1);
% bg_lh    = zeros(cp_total,bgwindow+1);
% bg_lh2   = zeros(cp_total,bgwindow-10+1);
bg_dD    = zeros(cp_total,bgwindow+1);
bg_dD2   = zeros(cp_total,bgwindow-10+1);
bg_d18O  = zeros(cp_total,bgwindow+1);
bg_d18O2 = zeros(cp_total,bgwindow-10+1);
bg_DXS   = zeros(cp_total,bgwindow+1);
bg_DXS2  = zeros(cp_total,bgwindow-10+1);

% The first and last cold pools must be removed, no full background data
% available for both
for k = 1:cp_total %
    bg_Taf(k,:)  = Taf(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_Taf2(k,:) = Taf(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_Ta(k,:)   = Ta(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_Ta2(k,:)  = Ta(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_rh(k,:)   = rh(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_rh2(k,:)  = rh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_prec(k,:) = rr(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_prec2(k,:)= rr(t_end_ind(k):t_end_ind(k)+bgwindow-10);
%     bg_wspd(k,:) = wspd(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_wspd2(k,:)= wspd(t_end_ind(k):t_end_ind(k)+bgwindow-10);
%     bg_wdir(k,:) = wdir(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_wdir2(k,:)= wdir(t_end_ind(k):t_end_ind(k)+bgwindow-10);
%     bg_qair(k,:) = qair(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_qair2(k,:)= qair(t_end_ind(k):t_end_ind(k)+bgwindow-10);
%     bg_u(k,:)    = u(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_u2(k,:)   = u(t_end_ind(k):t_end_ind(k)+bgwindow-10);
%     bg_v(k,:)    = v(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_v2(k,:)   = v(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
%     bg_sh(k,:)   = sh(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_sh2(k,:)  = sh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
%     bg_lh(k,:)   = lh(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_lh2(k,:)  = lh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_dD(k,:)   = dD(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_dD2(k,:)  = dD(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_d18O(k,:) = d18O(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_d18O2(k,:)= d18O(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_DXS(k,:)  = DXS(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_DXS2(k,:) = DXS(t_end_ind(k):t_end_ind(k)+bgwindow-10);

end

Taf_comp2  = [bg_Taf(:,1:end-1) Taf_comp(:,1:end-1) bg_Taf2(:,2:end)];
Ta_comp2   = [bg_Ta(:,1:end-1)  Ta_comp(:,1:end-1)  bg_Ta2(:,2:end)];
prec_comp2 = [bg_prec(:,1:end-1) prec_comp(:,1:end-1) bg_prec2(:,2:end)];
qair_comp2 = [bg_qair(:,1:end-1) qair_comp(:,1:end-1) bg_qair2(:,2:end)];
rh_comp2   = [bg_rh(:,1:end-1)   rh_comp(:,1:end-1)   bg_rh2(:,2:end)];
% wspd_comp2 = [bg_wspd(:,1:end-1) wspd_comp(:,1:end-1) bg_wspd2(:,2:end)];
% wdir_comp2 = [bg_wdir(:,1:end-1) wdir_comp(:,1:end-1) bg_wdir2(:,2:end)];
% u_comp2    = [bg_u(:,1:end-1)  u_comp(:,1:end-1)  bg_u2(:,2:end)];
% v_comp2    = [bg_v(:,1:end-1)  v_comp(:,1:end-1)  bg_v2(:,2:end)];
% sh_comp2   = [bg_sh(:,1:end-1) sh_comp(:,1:end-1) bg_sh2(:,2:end)];
% lh_comp2   = [bg_lh(:,1:end-1) lh_comp(:,1:end-1) bg_lh2(:,2:end)];
dD_comp2   = [bg_dD(:,1:end-1)   dD_comp(:,1:end-1) bg_dD2(:,2:end)];
d18O_comp2 = [bg_d18O(:,1:end-1) d18O_comp(:,1:end-1) bg_d18O2(:,2:end)];
DXS_comp2  = [bg_DXS(:,1:end-1)  DXS_comp(:,1:end-1) bg_DXS2(:,2:end)];
% q_iso_comp2 = [bg_q_iso(:,1:end-1) q_iso_comp(:,1:end-1) bg_q_iso2(:,2:end)];
% mr_comp2 = [bg_mr(:,1:end-1) mr_comp(:,1:end-1) bg_mr2(:,2:end)];

t_comp2 = [linspace(-60,-30.4,41) t_comp(2:end-1) linspace(24.9,60,30)]; %[-40:0,1:1.5:28,30:40];

%% Plots
orange = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
mustard = [0.9290 0.6940 0.1250];
green = [0.4196 0.5569 0.1373];

figure;
subplot(321)
    plot(t_comp2,nanmean(Taf_comp2),'-','LineWidth',2,'Color',orange)
    hold on;
    plot(t_comp2,nanmean(Taf_comp2)+std(Taf_comp2,'omitnan')/sqrt(16),'-','Color',orange)
    plot(t_comp2,median(Taf_comp2,'omitnan'),'-o','Color',orange)
    plot(t_comp2,nanmean(Taf_comp2)-std(Taf_comp2,'omitnan')/sqrt(16),'-','Color',orange)
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    box off
    ylabel('T_a [filtered] [\circC]')
    xlim([-60 60])
    ylim([24.6 26.5])
subplot(322)
    plot(t_comp2,nanmean(Ta_comp2),'-','LineWidth',2,'Color',orange)
    hold on;
    plot(t_comp2,nanmean(Ta_comp2)+std(Ta_comp2,'omitnan')/sqrt(16),'-','Color',orange)
    plot(t_comp2,median(Ta_comp2,'omitnan'),'-o','Color',orange)
    plot(t_comp2,nanmean(Ta_comp2)-std(Ta_comp2,'omitnan')/sqrt(16),'-','Color',orange)
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    box off
    ylabel('T_a [\circC]')
    xlim([-60 60])
    ylim([24.6 26.5])
subplot(323)
    plot(t_comp2,nanmean(rh_comp2),'-b','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(rh_comp2)+std(rh_comp2,'omitnan')/sqrt(16),'-b')
    plot(t_comp2,median(rh_comp2,'omitnan'),'-ob')
    plot(t_comp2,nanmean(rh_comp2)-std(rh_comp2,'omitnan')/sqrt(16),'-b')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    set(gca,'xticklabels',[])
    ylabel('RH [%]')
    box off
    xlim([-60 60])
    ylim([67 79])
subplot(324)
    plot(t_comp2,nanmean(dD_comp2),'-b','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(dD_comp2)+std(dD_comp2,'omitnan')/sqrt(16),'-b')
    plot(t_comp2,median(dD_comp2,'omitnan'),'-ob')
    plot(t_comp2,nanmean(dD_comp2)-std(dD_comp2,'omitnan')/sqrt(16),'-b')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel(['\deltaD [',char(8240),']'])
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-74 -70])
    legend('mean','std','median')
subplot(325)
    plot(t_comp2,nanmean(d18O_comp2),'-m','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(d18O_comp2)+std(d18O_comp2,'omitnan')/sqrt(16),'-m')
    plot(t_comp2,median(d18O_comp2,'omitnan'),'-om')
    plot(t_comp2,nanmean(d18O_comp2)-std(d18O_comp2,'omitnan')/sqrt(16),'-m')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel(['\delta^1^8O [',char(8240),']'])
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-10.6 -10.1])
subplot(326)
    plot(t_comp2,nanmean(DXS_comp2),'-','LineWidth',2,'Color',mustard)
    hold on;
    plot(t_comp2,nanmean(DXS_comp2)+std(DXS_comp2,'omitnan')/sqrt(16),'-','Color',mustard)
    plot(t_comp2,median(DXS_comp2,'omitnan'),'-o','Color',mustard)
    plot(t_comp2,nanmean(DXS_comp2)-std(DXS_comp2,'omitnan')/sqrt(16),'-','Color',mustard)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    set(gca,'xticklabels',[])
    ylabel(['d-excess [',char(8240),']'])
    xlim([-60 60])
    ylim([10.4 11.9])
    box off
    set(findall(gcf,'-property','Fontsize'),'FontSize',16)

%%    FLAG for recovery times    %%
%           AND                   %
%  FLAG for peak cold pool times  %
% Defining peak cold pool times as the 5-minute period centered on the minimum temperature time (t_min)
recovery_flag_1min = zeros(size(t1min));
peak_flag_1min = zeros(size(t1min));
for k = 1:length(t_max)
    ii  = t_min_ind(k):t_end_ind(k);
    iii = t_min_ind(k)-3:t_min_ind(k)+2;
    recovery_flag_1min(ii) = 1;
    peak_flag_1min(iii) = 1;
end
% recovery_flag_1min = 1 means value corresponds to a time flagged within a
% cold pool RECOVERY
% recovery_flag_1min = 0 means value corresponds to a time flagged within a
% cold pool FRONT
%% Identifying the strongest and weakest cold pools based on delta_T
num_vector = 1:length(t_max); %
date_vector = t_max;
cp_vector = [delta_T(num_vector)', num_vector', date_vector'];
sorted_cp = sortrows(cp_vector,1);