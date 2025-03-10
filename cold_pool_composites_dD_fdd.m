%% Cold pool composites for dD_dd, fdd, dD_ob %%
load conserved_properties_for_isotopes_1min_1.1to1.3km_linear_interp.mat
% load conserved_properties_for_isotopes_1min_1.1to1.3km_nearest_interp.mat
Ta = TA - 273.15; % air temperature at 17m [in degrees C]
rh = RH;
%% RUNNING COLD POOLS detection algorithm %%
[t_max,t_min,t_max_ind,t_min_ind,t_end,t_end_ind,~,~,delta_T,~,~,Taf] = cold_pool_detection_algorithm(t1min,Ta);
% Ta is filtered within the cold_pool_detection_algorithm function!
% Eliminating the incomplete cold pool because of the emergency in-port
t_max = t_max(:,[1:2,4:end]);
t_min = t_min(:,[1:2,4:end]);
t_end = t_end(:,[1:2,4:end]);
t_max_ind = t_max_ind(:,[1:2,4:end]);
t_end_ind = t_end_ind(:,[1:2,4:end]);
delta_T = delta_T(:,[1:2,4:end]);
%% Identifying the strongest and weakest cold pools based on delta_T
num_vector = 1:length(t_max); %
date_vector = t_max;
cp_vector = [delta_T(num_vector)', num_vector', date_vector'];
sorted_cp = sortrows(cp_vector,1);
%% FLAG for cold pool times (indexes)
cold_pool_flag_1min = zeros(size(t1min));
cp_matrix = zeros(length(t_max),121); % matrix of cold pool indexes
t_max0_matrix = 9999*ones(length(t_max),121); % time
t_min0_matrix = 9999*ones(length(t_max),121); % time

for k = 1:length(t_max)
    ii = t_max_ind(k):t_end_ind(k);
    cp_matrix(k,1:length(ii)) = ii;
    cold_pool_flag_1min(t_max_ind(k):t_end_ind(k)) = 1;
    t_max0_matrix(k,1:length(ii)) = (t1min(ii)-t_max(k))./(t_min(k)-t_max(k));
    t_min0_matrix(k,1:length(ii)) = (t1min(ii)-t_min(k))./(t_end(k)-t_min(k));
end
cp_matrix(cp_matrix==0) = NaN;
t_max0_matrix(t_max0_matrix==9999) = NaN;
t_min0_matrix(t_min0_matrix==9999) = NaN;
%% Interpolating to obtain a single normalized time vector
% For t_min @ 0 and t_end @ 1
cp_total = length(t_max);
t_norm    = 0:0.05:1;
% t_norm  = 0:0.10:3; % Simon says we don't need to normalize time!!!
Taf_norm  = 9999*ones(cp_total,length(t_norm));
Ta_norm   = 9999*ones(cp_total,length(t_norm));
rh_norm   = 9999*ones(cp_total,length(t_norm));
dD_norm   = 9999*ones(cp_total,length(t_norm));
fdd_norm  = 9999*ones(cp_total,length(t_norm));
dDdd_norm = 9999*ones(cp_total,length(t_norm));

for k = 1:length(t_max) % 49:62%[1:68,70:length(t_max)] % [65:68,70:86] %
    Taf_norm(k,:)  = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),Taf(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    Ta_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),Ta(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    rh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),rh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    dD_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),dD_ob(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    fdd_norm(k,:)  = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),fdd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    dDdd_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),dD_dd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
end
%% For t_max @ 0 and t_min @ 1
Taf_norm2  = zeros(cp_total,length(t_norm));
Ta_norm2   = zeros(cp_total,length(t_norm));
rh_norm2   = zeros(cp_total,length(t_norm));
dD_norm2   = zeros(cp_total,length(t_norm));
fdd_norm2  = zeros(cp_total,length(t_norm));
dDdd_norm2 = zeros(cp_total,length(t_norm));

for k = 1:length(t_max) %65:86
    Taf_norm2(k,:)  = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),Taf(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    Ta_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),Ta(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    rh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),rh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    dD_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),dD_ob(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    fdd_norm2(k,:)  = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),fdd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    dDdd_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),dD_dd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
end
%% Creating composite time series in normalized time
Taf_comp  = [Taf_norm2 Taf_norm(:,2:end)];
Ta_comp   = [Ta_norm2 Ta_norm(:,2:end)];
rh_comp   = [rh_norm2 rh_norm(:,2:end)];
dD_comp   = [dD_norm2 dD_norm(:,2:end)];
fdd_comp  = [fdd_norm2 fdd_norm(:,2:end)];
dDdd_comp = [dDdd_norm2 dDdd_norm(:,2:end)];

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
bg_dD    = zeros(cp_total,bgwindow+1);
bg_dD2   = zeros(cp_total,bgwindow-10+1);
bg_fdd   = zeros(cp_total,bgwindow+1);
bg_fdd2  = zeros(cp_total,bgwindow-10+1);
bg_dDdd  = zeros(cp_total,bgwindow+1);
bg_dDdd2 = zeros(cp_total,bgwindow-10+1);
% The first and last cold pools must be removed, no full background data
% available for both
for k = 1:length(t_max) % 49:62 % [65:68,70:86]
    bg_Taf(k,:)  = Taf(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_Taf2(k,:) = Taf(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_Ta(k,:)   = Ta(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_Ta2(k,:)  = Ta(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_rh(k,:)   = rh(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_rh2(k,:)  = rh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_dD(k,:)   = dD_ob(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_dD2(k,:)  = dD_ob(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_fdd(k,:)  = fdd(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_fdd2(k,:) = fdd(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_dDdd(k,:) = dD_dd(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_dDdd2(k,:)= dD_dd(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
end

Taf_comp2  = [bg_Taf(:,1:end-1) Taf_comp(:,1:end-1) bg_Taf2(:,2:end)];
Ta_comp2   = [bg_Ta(:,1:end-1) Ta_comp(:,1:end-1) bg_Ta2(:,2:end)];
rh_comp2   = [bg_rh(:,1:end-1) rh_comp(:,1:end-1) bg_rh2(:,2:end)];
dD_comp2   = [bg_dD(:,1:end-1) dD_comp(:,1:end-1) bg_dD2(:,2:end)];
fdd_comp2  = [bg_fdd(:,1:end-1) fdd_comp(:,1:end-1) bg_fdd2(:,2:end)];
dDdd_comp2 = [bg_dDdd(:,1:end-1) dDdd_comp(:,1:end-1) bg_dDdd2(:,2:end)];

t_comp2 = [linspace(-60,-30.4,41) t_comp(2:end-1) linspace(24.9,60,30)]; %[-40:0,1:1.5:28,30:40];
%% Weigthed means by fdd 
dD_ob_w = mean(dD_comp2.*fdd_comp2,'omitnan')./mean(fdd_comp2,'omitnan');
dD_dd_w = mean(dDdd_comp2.*fdd_comp2,'omitnan')./mean(fdd_comp2,'omitnan');
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
    plot(t_comp2,dD_ob_w,'-k','LineWidth',2)
    hold on;
    plot(t_comp2,dD_ob_w+var(dD_comp2,'omitnan'),'-k')
%     plot(t_comp2,median(dD_comp2,'omitnan'),'-ok')
    plot(t_comp2,dD_ob_w-var(dD_comp2,'omitnan'),'-k')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel(['\deltaD_o_b (weighted) [',char(8240),']'])
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-90 -60])
subplot(325)
    plot(t_comp2,dD_dd_w,'-k','LineWidth',2)
    hold on;
    plot(t_comp2,dD_dd_w+var(dDdd_comp2,'omitnan'),'-k')
%     plot(t_comp2,median(dDdd_comp2,'omitnan'),'-ok')
    plot(t_comp2,dD_dd_w-var(dDdd_comp2,'omitnan'),'-k')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel(['\deltaD_d_d (weighted) [',char(8240),']'])
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-120 -40])
subplot(326)
    plot(t_comp2,nanmean(fdd_comp2),'-k','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(fdd_comp2)+std(fdd_comp2,'omitnan')/sqrt(16),'-k')
    plot(t_comp2,median(fdd_comp2,'omitnan'),'-ok')
    plot(t_comp2,nanmean(fdd_comp2)-std(fdd_comp2,'omitnan')/sqrt(16),'-k')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('fdd')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([.12 .4])
    legend('mean','std','median')