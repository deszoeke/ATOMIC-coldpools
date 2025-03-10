% function [] = q_dD_framework_cold_pools(cp_num)
R_VSMOW = 155.76*10^(-6); % FOR Deuterium
%% Loading datasets
load '1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' qair t1min Ta rr sst slp
load 'iso_data_1min_intervals_FLAGGED_w_runningmean.mat' d18O dD DXS iso_time
ind = find(t1min>=iso_time(1));
ind = ind(1);
Taf   = Ta(ind:ind+length(iso_time)-1);   % Ta filtered   => 11-min running average
qairf = qair(ind:ind+length(iso_time)-1); % qair filtered => 11-min running average
sstf  = sst(ind:ind+length(iso_time)-1);  % sst filtered  => 11-min running average
slpf  = slp(ind:ind+length(iso_time)-1);  % slp filtered  => 11-min running average
rrf   = rr(ind:ind+length(iso_time)-1);
dDf   = dD;   % filtered => 11-min running average
d18Of = d18O; % filtered => 11-min running average
DXSf  = DXS;  % filtered => 11-min running average
time = t1min(ind:ind+length(iso_time)-1);
rrf(rrf<=0) = NaN;

%% Identifying cold pools
del_T = Taf(2:end)-Taf(1:end-1);
cand  = find(del_T<-0.05); % candidates positions
%% for plotting purposes
% figure; 
subplot(511);hold on;
    plot(t1min,Ta);
    plot(time,Taf,'LineWidth',1)
    ylabel('Abs Temp [K]')
    datetick('x','mmm/dd','keeplimits','keepticks')
subplot(512);hold on;
    plot(t1min,qair)
    plot(time,qairf,'LineWidth',1)
    datetick('x','mmm/dd','keeplimits','keepticks')
    ylabel('q [g/kg]')
subplot(513);hold on;
    plot(time,dD)
    plot(time,dDf,'LineWidth',1)
    datetick('x','mmm/dd','keeplimits','keepticks')
    ylabel('\deltaD [permil]')
    yyaxis right
    plot(time,rrf,'-b')
    ylabel('RR [mm/hr]')
subplot(514);hold on;
    plot(time,d18O)
    plot(time,d18Of,'LineWidth',1)
    datetick('x','mmm/dd','keeplimits','keepticks')
    ylabel('\delta^1^8O [permil]')
subplot(515);hold on;
    plot(time,DXS)
    plot(time,DXSf,'LineWidth',1)
    datetick('x','mmm/dd','keeplimits','keepticks')
    ylabel('DXS [permil]')
xlabel('2020')

%% Separating candidates [28 possible cold pools]
c=1;
cand_ind(c) = cand(1);
for i = 2:length(cand)-1
    if cand(i)-cand(i-1)>1 && cand(i+1)-cand(i)==1
        c = c+1;
        cand_ind(c) = cand(i);
    end
end
% Verifying cold pool validity
c = 0;
for i = 1:length(cand_ind)
    if Taf(cand_ind(i)) >= nanmean(Taf)+0.6
        c = c + 1; 
        disp(c)
        cand_ind(i) = 0;
    end
end
cand_ind(cand_ind==0) = [];

% Time passed between possible cold pools
t_pass = time(cand_ind(2:end)) - time(cand_ind(1:end-1)); % in days
t_pass = [1 t_pass*24]; % in hours
cand_ind(t_pass<1) = [];
t_max_offset = 20; % in minutes
% cand_ind(2) = []; % to join candidates 8 and 11 into a single cold pool

%% Identifying t_max: time of the cold pool front onset
c=1;
for k = 1:length(cand_ind)
    for ii = cand_ind(k):-1:1
        if del_T(ii)>0 && ((time(cand_ind(k))-time(ii+1))*24*60)<t_max_offset
            t_max(c) = time(ii+1);
            t_max_ind(c) = ii+1;
            c = c+1;
            break
        elseif del_T(ii)>0 && ((time(cand_ind(k))-time(ii+1))*24*60)>=t_max_offset
            t_max(c) = time(cand_ind(k)-t_max_offset);
            t_max_ind(c) = (cand_ind(k))-t_max_offset;
            c = c+1;
            break
        end
    end
end
% Time passed between t_max times
% t_p = t1min(cand_ind) - t_max; % in days
% t_p = t_p*24*60; % in min
t_p = cand_ind - t_max_ind; % in min
find(t_p>20);

%% Identifying t_min: time of end of cold pool front
c=1;
for k = 1:length(cand_ind)
    for ii = cand_ind(k):1:length(del_T)
        if del_T(ii)>0
            t_min(c) = time(ii);
            t_min_ind(c) = ii;
            c = c+1;
            break
        elseif ii == length(del_T)
            t_min(c) = time(ii);
            t_min_ind(c) = ii;
            break
        end
    end
end
t_min_offset = 20; % within 20min of the previous minimum

%% Identifying t_end: end of cold pool
delta_T = Taf(t_max_ind) - Taf(t_min_ind);
T_min = Taf(t_min_ind);
max_len = 100; % max length of cold pool

[t_end,t_end_ind,end_flag] = t_end_w_max_length(delta_T,T_min,max_len,cand_ind,t_min_ind,Taf,time);
% [t_end,t_end_ind,end_flag] = t_end_w_next_cp(delta_T,T_min,cand_ind,t_min_ind,Taf,time);
%% FLAG for cold pool times (indexes)
t_max(isnan(t_end)) = [];
t_min(isnan(t_end)) = [];
% end_flag(isnan(t_end)) = [];
t_end(isnan(t_end)) = [];

t = time;
cold_pool_flag_1min = zeros(size(t));
cp_matrix = ones(length(t_max),61); % matrix of cold pool indexes
for k = 1:length(t_max)-1
    max_ind(k) = find(t==t_max(k));
    min_ind(k) = find(t==t_min(k));
    end_ind(k) = find(t==t_end(k));
    ii = max_ind(k):end_ind(k);
    cp_matrix(k,1:length(ii)) = ii;
    cold_pool_flag_1min(max_ind(k):end_ind(k)) = 1;
end
cp_matrix(cp_matrix==0) = 700; % Taf(700) = NaN
cp_matrix(cp_matrix==1) = 700; % Taf(700) = NaN

%% for plotting purposes
% figure;
hold on;
for k = 1:length(t_max)-1
    patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1],[.8 .8 .8],'EdgeColor','none')
end
plot(t,cold_pool_flag_1min*16,'-r','Color','r')
ylim([0 1])
datetick('x','mmm/dd','keeplimits','keepticks')
title('cold pool flag'); 
grid on

%% for plotting purposes
% figure;
% subplot(511); 
%     hold on;
%     for k = 1:length(t_max)
%         patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*max(Ta),[.8 .8 .8],'EdgeColor','none')
%     end
%     plot(date_t,Ta)
%     plot(date_t,Taf,'LineWidth',2)
%     ylim([min(Ta) max(Ta)])
%     ylabel('Abs Temp [K]')
%     datetick('x','HH:MM','keeplimits','keepticks')
% subplot(512); 
%     hold on;
%     for k = 1:length(t_max)
%         patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*max(q),[.8 .8 .8],'EdgeColor','none')
%     end
%     plot(date_t,q)
%     plot(date_t,qairf,'LineWidth',2)
%     ylim([min(q) max(q)])
%     ylabel('q [g/kg]')
%     datetick('x','HH:MM','keeplimits','keepticks')
% subplot(513); 
%     hold on;
%     for k = 1:length(t_max)
%         patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*min(dD),[.8 .8 .8],'EdgeColor','none')
%     end
%     plot(date_t,dD)
%     plot(date_t,dDf,'LineWidth',2)
%     ylim([min(dD) max(dD)])
%     ylabel('\deltaD [permil]')
%     datetick('x','HH:MM','keeplimits','keepticks')
% subplot(514); 
%     hold on;
%     for k = 1:length(t_max)
%         patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*min(d18O),[.8 .8 .8],'EdgeColor','none')
%     end
%     plot(date_t,d18O)
%     plot(date_t,d18Of,'LineWidth',2)
%     ylim([min(d18O) max(d18O)])
%     ylabel('\delta^1^8O [permil]')
%     datetick('x','HH:MM','keeplimits','keepticks')
% subplot(515); 
%     hold on;
%     for k = 1:length(t_max)
%         patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*min(DXS),[.8 .8 .8],'EdgeColor','none')
%     end
%     plot(date_t,DXS)
%     plot(date_t,DXSf,'LineWidth',2)
%     ylim([min(DXS) max(DXS)])
%     ylabel('\DXS [permil]')
%     datetick('x','HH:MM','keeplimits','keepticks')
% xlabel('Feb-09-2020')
% set(findall(gcf,'-property','Fontsize'),'FontSize',14)

%% cold pool matrices
for k = 1:length(max_ind)
    dummy = t(max_ind(k):end_ind(k));
    tcold(k,1:length(dummy)) = t(max_ind(k):end_ind(k));
    qcold(k,1:length(dummy)) = qairf(max_ind(k):end_ind(k));
    dDcold(k,1:length(dummy))= dD(max_ind(k):end_ind(k));
   DXScold(k,1:length(dummy))= DXS(max_ind(k):end_ind(k));
  d18Ocold(k,1:length(dummy))= d18O(max_ind(k):end_ind(k));
end
%% Creating new variables where cold pools times are NaN's
qnan = qairf;
qnan(cold_pool_flag_1min==1)= NaN;
dDnan = dD;
dDnan(cold_pool_flag_1min==1)= NaN;
d18Onan = d18O;
d18Onan(cold_pool_flag_1min==1)= NaN;
DXSnan = DXS;
DXSnan(cold_pool_flag_1min==1)= NaN;

for k = 1:length(max_ind)
    dummy = t(max_ind(k)-61:end_ind(k));
    qcoldnan(k,1:length(dummy)) = qnan(max_ind(k)-61:end_ind(k));
    dDcoldnan(k,1:length(dummy))= dDnan(max_ind(k)-61:end_ind(k));
    d18Ocoldnan(k,1:length(dummy))= d18Onan(max_ind(k)-61:end_ind(k));
    DXScoldnan(k,1:length(dummy)) = DXSnan(max_ind(k)-61:end_ind(k));    
end

% 'f' subscript for front
for k = 1:length(max_ind)
    dummy = t(max_ind(k):min_ind(k)-1);
    tcoldf(k,1:length(dummy)) = t(max_ind(k):min_ind(k)-1);
    qcoldf(k,1:length(dummy)) = qairf(max_ind(k):min_ind(k)-1);
    dDcoldf(k,1:length(dummy))= dD(max_ind(k):min_ind(k)-1);
   DXScoldf(k,1:length(dummy))= DXS(max_ind(k):min_ind(k)-1);
  d18Ocoldf(k,1:length(dummy))= d18O(max_ind(k):min_ind(k)-1);
end

% 'w' subscript for wake
for k = 1:length(max_ind)
    dummy = t(min_ind(k):end_ind(k));
    tcoldw(k,1:length(dummy)) = t(min_ind(k):end_ind(k));
    qcoldw(k,1:length(dummy)) = qairf(min_ind(k):end_ind(k));
    dDcoldw(k,1:length(dummy))= dD(min_ind(k):end_ind(k));
   DXScoldw(k,1:length(dummy))= DXS(min_ind(k):end_ind(k));    
  d18Ocoldw(k,1:length(dummy))= d18O(min_ind(k):end_ind(k));
end

%% converting zeros in NaN's
tcoldf(tcoldf==0) = NaN;
qcoldf(qcoldf==0) = NaN;
dDcoldf(dDcoldf==0) = NaN;
DXScoldf(DXScoldf==0) = NaN;
d18Ocoldf(d18Ocoldf==0) = NaN;

tcoldw(tcoldw==0) = NaN;
qcoldw(qcoldw==0) = NaN;
dDcoldw(dDcoldw==0) = NaN;
DXScoldw(DXScoldw==0) = NaN;
d18Ocoldw(d18Ocoldw==0) = NaN;

tcold(tcold==0) = NaN;
qcold(qcold==0) = NaN;
dDcold(dDcold==0) = NaN;
DXScold(DXScold==0) = NaN;
d18Ocold(d18Ocold==0) = NaN;

qcoldnan(qcoldnan==0) = NaN;
dDcoldnan(dDcoldnan==0) = NaN;
DXScoldnan(DXScoldnan==0) = NaN;
d18Ocoldnan(d18Ocoldnan==0) = NaN;

%% Craig-Gordon surface estimates 
load 'RHS_Eq9_MJ79_LIMITED.mat' alpha_e alpha_e_D del_oc del_oc_D
load 'conserved_variables_10minLIMITED.mat' q_en q_surf time_q 
delta_e_D = ((1./alpha_e_D)*(1+del_oc_D) - 1) * 1000;
delta_e_O = ((1./alpha_e)*(1+del_oc) - 1) * 1000;
for k = 1:length(max_ind)
    ind = find(time_q >= t_max(k));
    start_t_ids(k)    = ind(1);
    qsurf_cold(k)     = q_surf(start_t_ids(k));
    dDsurf_cold(k)    = delta_e_D(start_t_ids(k));
    qent_cold(k)      = q_en(start_t_ids(k));
    d18Osurf_cold(k)  = delta_e_O(start_t_ids(k));
    DXSsurf_cold(k)   = dDsurf_cold(k) - (8*d18Osurf_cold(k));
end
%% Identifying centroid for 60-minute data before cold pool onset
% SIMPLE APPROACH: Assuming dD and q at the cold pool onset are the same as dD and q for the centroid
total_cp = length(max_ind); % total number of cold pools
q_cp_mean = qcoldf(1:total_cp,1);
dD_cp_mean = dDcoldf(1:total_cp,1);
d18O_cp_mean = d18Ocoldf(1:total_cp,1);
DXS_cp_mean = DXScoldf(1:total_cp,1);

% CENTROID APPROACH: Using centroid in the q*dD-q space
q_cp_mean = nanmean(qcoldnan(1:total_cp,1:61),2);
dD_cp_mean = nanmean(dDcoldnan(1:total_cp,1:61),2);
d18O_cp_mean = nanmean(d18Ocoldnan(1:total_cp,1:61),2);
DXS_cp_mean = nanmean(DXScoldnan(1:total_cp,1:61),2);

[m1,b1,b]   = mixing_line_slope_yint(qsurf_cold,d18Osurf_cold,q_cp_mean',d18O_cp_mean');
[m2,b2,bb]  = mixing_line_slope_yint(qsurf_cold,dDsurf_cold,q_cp_mean',dD_cp_mean');
[m3,b3,bbb] = mixing_line_slope_yint(qsurf_cold,DXSsurf_cold,q_cp_mean',DXS_cp_mean');

% Y1 = (qsurf_cold.*d18Osurf_cold);
% Y2 = (qsurf_cold.*dDsurf_cold);
% Y3 = (qsurf_cold.*DXSsurf_cold);
% m1 = (q_cp_mean'.*d18O_cp_mean' - Y1)./(q_cp_mean' - qsurf_cold); % slope
% m2 = (q_cp_mean'.*dD_cp_mean' - Y2)./(q_cp_mean' - qsurf_cold); % slope
% m3 = (q_cp_mean'.*DXS_cp_mean' - Y3)./(q_cp_mean' - qsurf_cold); % slope
% b = q_cp_mean'.*d18O_cp_mean' - (m1.*q_cp_mean');
% b1 = Y1 - (m1.*qsurf_cold); % same results for both points: good check
% bb = q_cp_mean'.*dD_cp_mean' - (m2.*q_cp_mean');
% b2 = Y2 - (m2.*qsurf_cold); % same results for both points: good check
% bbb = q_cp_mean'.*DXS_cp_mean' - (m3.*q_cp_mean');
% b3 = Y3 - (m3.*qsurf_cold); % same results for both points: good check
qXd18O_ent = (m1.*qent_cold) + b1;
qXdD_ent   = (m2.*qent_cold) + b2;
qXDXS_ent  = (m3.*qent_cold) + b3;
d18O_ent = qXd18O_ent./qent_cold;
dD_ent = qXdD_ent./qent_cold;
DXS_ent = qXDXS_ent./qent_cold;

%% Incorporating Rayleigh curves in q-dD diagram 
% for peak dD
[dD0p,ind] = max(dDcold,[],2,'omitnan'); % peak: ~-70 permil
    Rvsmow = 155.76e-6; % unitless, deuterium
R0p = Rvsmow * (1 + dD0p*1e-3); % dD0 in permil

for i = 1:length(max_ind)
    q0p(i) = qcold(i,ind(i)); % peak: ~15 g/kg
    T0p(i) = Taf(cp_matrix(i,ind(i))) + 273.15; % peak: ~298 Kelvin
    p0p(i) = slpf(cp_matrix(i,ind(i)))*10^2; % pressure [in Pa]
    Tsurf_cold(i) = sstf(cp_matrix(i,ind(i))) + 273.15; % surface (SST): ~299 Kelvin
    test(i) = dDcold(i,ind(i)); % peak: ~-70 permil
end
% for surface source
dD0 = dDsurf_cold; % surface: ~-64 permil
R0 = Rvsmow * (1 + dD0*1e-3); % dD0 in permil
q0  = qsurf_cold;  % surface: ~ 21 g/kg
T0  = Tsurf_cold;  % surface: ~299 Kelvin (SST)

%% Plots for dD [for plotting purposes]
for i = [2:11,13:16] %[9,14:15] %[8,10,11,13,16] %1:length(max_ind) %[2:7,9,14:15] %1:length(max_ind)-1 %
    % Calculating slope of best fit line for the 60-min data 
    slope = polyfit(qcoldnan(i,1:61),dDcoldnan(i,1:61),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,qcoldnan(i,1:61));
% figure; 
set(gcf, 'Position', get(0, 'Screensize'));
j = plot(qcoldf(i,:),dDcoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(qcoldw(i,:),dDcoldw(i,:),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 60 minutes) in gray dots
jj = scatter(qcoldnan(i,1:61),dDcoldnan(i,1:61),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qcoldnan(i,1:61),y_est,'r-','LineWidth',2)
% Cold pool data (up to 60 minutes) in time-colored dots
plot(qcoldf(i,:),dDcoldf(i,:),'ok','MarkerSize',9);%,'MarkerFaceColor','b')
jjj = scatter(qcold(i,1:61),dDcold(i,1:61),70,(tcold(i,1:61)-tcold(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qcoldw(i,:),dDcoldw(i,:),'sk','MarkerSize',11);%,'MarkerFaceColor','c')

% Adding mixing lines between end member points (entraiment & surface)
Xfit = linspace(qent_cold(i),qsurf_cold(i),10);
Ylinearfit = (m2(i).*Xfit) + b2(i);
Ycurvefit = (Ylinearfit)./Xfit;
j3 = plot(Xfit,Ycurvefit,'--k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
clearvars Xfit Ycurvefit Ylinearfit
% Identifying t_max & t_min
plot(qcoldf(i,1),dDcoldf(i,1),'^k','MarkerFaceColor','k','MarkerSize',12)
plot(qcoldw(i,1),dDcoldw(i,1),'^k','MarkerFaceColor','w','MarkerSize',12)

% % Incorporating Rayleigh curves
% % for peak dD
% [q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0p(i),T0p(i),R0p(i));
% plot(q_rayleigh,delt_rayleigh,'-b','LineWidth',.5)
% % for surface source
% [q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0(i),T0(i),R0(i));
% plot(q_rayleigh,delt_rayleigh,'-m','LineWidth',.5)

% Adding mixing lines between Rayleigh liquid evap curve and surface source
% Adding mixing lines between end member points (hydrometeor & surface)
dD_prec = 10; % in permil
[q_lost,dD_lost,f_lost] = Rayleigh_liquid_evap(p0p(i),q0p(i),T0p(i),dD_prec); % p[Pa],q[g/kg],T[K],dD[permil]
% f = 0.75; % remaining fraction of liquid drop
cmap = b2rcolormap(length(f_lost)+1);
for kk = 1:length(f_lost)
    [m,b,b1] = mixing_line_slope_yint(qsurf_cold(i),dDsurf_cold(i),q_lost(kk),dD_lost(kk));
    Xfit = linspace(q_lost(kk),qsurf_cold(i),500);
    Ylinearfit = (m.*Xfit) + b;
    Ycurvefit = (Ylinearfit)./Xfit;
    j4 = plot(Xfit,Ycurvefit,'--','Color',cmap(kk,:),'LineWidth',2);
    j4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars m b b1 Xfit Ycurvefit Ylinearfit j4
end

% Adding markers for end member points (entraiment & surface)
orange = [0.8500 0.3250 0.0980];
mustard = [0.9290 0.6940 0.1250];
blue = [0 0.4470 0.7410];
plot(qent_cold(i),dD_ent(i),'ok','MarkerFaceColor',orange,'MarkerSize',10)
plot(qsurf_cold(i),dDsurf_cold(i),'ok','MarkerFaceColor',mustard,'MarkerSize',10)

% Colorbar properties
colormap(jet(12))
han = colorbar;
han.Title.String = "minutes since cp onset";
caxis([0 60])

% Figure properties
xlim([8 23]); 
ylim([-84 -64]); % ylim([-90 -64])
xlabel('specific humidity q [g/kg]')
ylabel(['isotopic composition \deltaD [',char(8240),']'])
box on
grid on
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
% legend('background','front','wake','entrainment','surface','t_0','t_m_i_n','Rayleigh peak \deltaD','Rayleigh surface','Location','southwest')

% Including Ta,q,dD timeseries in small plot
inset_plot(i,tcold,Taf,cp_matrix,qairf,dDf,rrf)
% Including liquid drop evaporation in another small plot
axes('Position',[.16 .75 .15 .15]); box on
plot(q_lost,dD_lost,'-r');
for kk = 1:length(f_lost)-1
    hold on;
    plot(q_lost(kk:kk+1),dD_lost(kk:kk+1),'-','Color',cmap(kk,:),'LineWidth',2);
end
xlabel('q_l_o_s_t');ylabel(['\deltaD_l_o_s_t [',char(8240),']'])

% Saving figure in different formats
% saveas(gcf,['model&obs_dD_vs_q_CP#',num2str(i),'_curve_fit_centroid.png'])
% saveas(gcf,['model&obs_dD_vs_q_CP#',num2str(i),'_curve_fit_centroid.fig'])

end

%%
% figure; hold on
for kk = 1:length(f_lost)
    [m,b,b1] = mixing_line_slope_yint(q0p(i),dD0p(i),q_lost(kk),dD_lost(kk));
    Xfit = linspace(q_lost(kk),qsurf_cold(i),2); % q
    Ylinearfit = (m.*Xfit) + b; % q*dD
    Ycurvefit = (Ylinearfit)./Xfit; % dD
    plot(1./Xfit,Ycurvefit,'-o','Color',1-cmap(kk,:))%,'MarkerFaceColor',cmap(kk,:));
    % 1-cmap => (darker color in middle of colorbar)
    % xlim([0 18*10^4]); ylim([-68 8])
    clearvars m b b1 Xfit Ycurvefit Ylinearfit
end
xlabel('1/q')
ylabel(['\deltaD [',char(8240),']'])

clearvars y_est slope q_rayleigh delt_rayleigh q_lost dD_lost m b b1

% end
