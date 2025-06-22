% conserved_properties_for_isotopes.m
%
% Applying the extrapolation and mixing fraction analysis to RHB data in
% order to obtain isotope concentrations for entrained and downdraft end members
% EQM Last modified: Nov 16 2022
% SPdeS tweaked some paths 2025-06-05

addpath('./thermo'); 

%% Air type: observed %%
% q and theta @ 400m
zrf = 400; % [in meters]; reference height
zm = 17;   % [in meters]; measurement level
Rd = 287.04;
Cp = 1005.7;
[th_ob,q_ob,t_adj] = height_adj(zrf,zm,Rd,Cp);
load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' t
pos1 = find(t_adj>=t(1)); % rounding up to closest (in time) PSD surface data point!!!
pos_adj = pos1(1);
clearvars pos1
q_ob = q_ob(pos_adj:pos_adj+length(t)-1);   % keeping 2nd leg data only
th_ob = th_ob(pos_adj:pos_adj+length(t)-1); % keeping 2nd leg data only
t_10min = t;
% Loading iso data %
load('data/iso_data_1min_intervals_FLAGGED_w_runningmean.mat');% RHB Picarro data
load 'data/cold_pool_flag_1min.mat' t1min
ind0 = find(t1min>=iso_time(1));
ind0 = ind0(1);
% Incorporating isotope data %
pos_i = 999999*ones(size(t));
for l = 1:length(t)
    pos2 = find(iso_time>=t(l)); % rounding up to closest isotope surface data point!!!
                                 % try rounding to nearest data point
    pos_i(l) = pos2(1);
    clearvars pos2
end
d18O_ob = d18O(pos_i);
dD_ob   = dD(pos_i);
de_ob   = dD_ob - (8*d18O_ob); % deuterium excess

%% Air type: surface [using properties in equilibrium w/ surface ocean] %%
% q and theta @ 1m
pos = 999999*ones(size(t_10min));
for l = 1:length(t_10min)
    pos1 = find(t1min>=t_10min(l)); % rounding up to closest (in time) PSD surface data point!!!
                                    % try rounding to nearest (in space) data point
    pos(l) = pos1(1);
    clearvars pos1
end
load 'data/1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' sst slp Ta rh
T = Ta(pos(1:end))+273.15; % in degrees Kelvin
SLP = slp(pos(1:end)); % [in hPa]
SST = sst(pos(1:end)); % [in degrees C]
RH  =  rh(pos(1:end));  % in percent
[q_surf,th_surf] = air_types(0,0,0,'surface',SLP,SST);

%% Air type: surface [using Craig-Gordon] %%
% Incorporating isotope data %
h_prime = 0.875; % relative humidity at which we are modeling the iso ratios
filename = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
Tliq = ncread(filename,'tskin'); % liquid temperature at [???]m [in degrees C]
[T_skin,~] = colocate_in_time(filename,Tliq);
% [dD_surf,d18O_surf] = air_types_iso(dD_ob,d18O_ob,'surface',T_skin',RH'/100,h_prime);
[dD_surf,d18O_surf] = air_types_iso(dD_ob,d18O_ob,'surface',SST',RH'/100,h_prime);
% load air_type_surface_CG.mat
% d18O_surf = delta_aO;
% dD_surf = delta_aD;
% de_surf = DXS_a;
de_surf = dD_surf - (8*d18O_surf); % deuterium excess

%% Air type: entrained %%
% q and theta @ 1km
load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' q th h
h_low = [ 800  900 1000 1100 1200];
h_hi  = [1000 1100 1200 1300 1400];
q_ent_1km  = q(h==1e3,:)*1e3; % check units for q_inth, looks like is cg/kg
th_ent_1km = th(h==1e3,:);
for k = 1:length(h_hi)
     q_ent(k,:)  = mean( q(h>=h_low(k) & h<=h_hi(k),:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
    th_ent(k,:) = mean(th(h>=h_low(k) & h<=h_hi(k),:),'omitnan');
    % Incorporating isotope data %
    [dD_ent(k,:),d18O_ent(k,:)] = iso_estimates_from_centroid_approach(q_ob,dD_ob,d18O_ob,q_surf,dD_surf,d18O_surf,q_ent(k,:));
     de_ent(k,:) = dD_ent(k,:) - (8*d18O_ent(k,:)); % deuterium excess
end
2+2;
%% Air type: downdraft %%
% from mean cloud layer
% q and theta @ mean theta_w (wet-bub potential temp) above 1km and below the trade inversion line (6g/kg contour) for each sounding
% Extracting trade inversion height (mixed layer depth)
for k = 12:size(q,2)
    h6 = double(h(q(:,k)*1000<=6));
    if h6(1) <= 6000
        trade_inv(k) = h6(1); % trade inversion
    else
        trade_inv(k) = NaN; % trade inversion
    end
    clearvars h6
end
thw_index = find(h==1e3);
load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' thw
for k = 1:size(thw,2)
    if trade_inv(k) >= 1e3
        trade_index = find(h==trade_inv(k));
        th_dd(k) = nanmean(thw(thw_index:trade_index,k));
    else
        th_dd(k) = NaN; % trade inversion
    end
end
% addpath('C:\Users\quinones\Documents\Data\thermo')
q_dd = qs(1e3*1e2,th_dd-273.15)*1e3; % q_d in g/kg    
    % qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
    % p [Pa], T [degrees C], qs [kg/kg]
    % theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
%% Isotope data for downdraft end member %%
% Isotope concentrations calculated based on the mixing fractions obtained from soundings
for k = 1:length(h_hi)
    [fen(k,:), fss(k,:), fdd(k,:)] = th_q_to_mixfraction(th_ob,q_ob, th_ent(k,:),q_ent(k,:), th_surf,q_surf, th_dd,q_dd);
    dD_dd(k,:) = (q_ob.*dD_ob' - q_ent(k,:).*fen(k,:).*dD_ent(k,:) - q_surf.*fss(k,:).*dD_surf')./(q_dd.*fdd(k,:));
    d18O_dd(k,:) = (q_ob.*d18O_ob' - q_ent(k,:).*fen(k,:).*d18O_ent(k,:) - q_surf.*fss(k,:).*d18O_surf')./(q_dd.*fdd(k,:));
    de_dd(k,:) = dD_dd(k,:) - (8*d18O_dd(k,:)); % deuterium excess
end

fdd(fdd<0)=NaN;
fss(fss<0)=NaN;
fen(fen<0)=NaN;

for k = 1:length(h_hi)
    dD_dd(k,:) = (q_ob.*dD_ob' - q_ent(k,:).*fen(k,:).*dD_ent(k,:) - q_surf.*fss(k,:).*dD_surf')./(q_dd.*fdd(k,:));
    d18O_dd(k,:) = (q_ob.*d18O_ob' - q_ent(k,:).*fen(k,:).*d18O_ent(k,:) - q_surf.*fss(k,:).*d18O_surf')./(q_dd.*fdd(k,:));
    de_dd(k,:) = dD_dd(k,:) - (8*d18O_dd(k,:)); % deuterium excess
end

% fdd_mask = fdd;
% fss_mask = fss;
% fen_mask = fen;
% 
% fdd_mask(iso_cold_pool_flag_10min==0) = NaN;
% fss_mask(iso_cold_pool_flag_10min==0) = NaN;
% fen_mask(iso_cold_pool_flag_10min==0) = NaN;

% X = -1:.1:1;
% Ydd = normpdf(X,mean(fdd,'omitnan'),std(fdd,'omitnan'));
% Yss = normpdf(X,mean(fss,'omitnan'),std(fss,'omitnan'));
% Yen = normpdf(X,mean(fen,'omitnan'),std(fen,'omitnan'));
%% Cold pool times
load 'data/cold_pool_flag_1min.mat' cold_pool_flag_1min
cold_pool_flag_10min = cold_pool_flag_1min(1:10:end);
ind1 = find(t_adj>=t_10min(1));
ind1 = ind1(1);
iso_cold_pool_flag = cold_pool_flag_10min(ind1:ind1+length(t_10min)-1);
% flag = 0; % script outputs plot for outside of cold pool times 
flag = 1; % script outputs cold pool plot
% dummy = 9999*ones(size(iso_cold_pool_flag));
% dummy(1:24:end) = iso_cold_pool_flag(1:24:end);
% dummy(dummy==9999) = NaN;
% sounding_flag = dummy;
2+2;

%% Answering Mampi's question: Do cold pool cases with lower dD exhibit more partial evaporation?
% % Step #1
% q_evap = q_dd - q_ent(4,:);
% % Step #2
% % Given dD_ent, dD_dd, q_ent and q_dd, find dD_evap
% % Linear equation: q_evap.*(1 + dD_evap) = q_dd.*(1 + dD_dd) - q_ent.*(1 + dD_ent);
% % Solve for dD_evap:
% dD_evap = (q_dd.*dD_dd(4,:) - q_ent(4,:).*dD_ent(4,:))./(q_evap);
% % Step #3
% % Estimate fraction that is complete evaporated
% % cf = ???; % completely evaporated fraction;
% filename = 'EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
% dDp = ncread(filename,'dD'); %
% d18Op = ncread(filename,'d18O'); %
% ctime = ncread(filename,'collection_time'); %
% ctime = datenum('01012020','mmddyyyy') + ctime*(1/(3600*24));
% delay_flag = ncread(filename,'delay_flag'); %
% 
% T_oc = T;
% alpha_D = exp( 1158.8*(T_oc.^3./10^12) - 1620.1*(T_oc.^2./10^9)...
%              + 794.84*(T_oc./10^6) - 161.04/10^3 + 2.9992*(10^6./T_oc.^3) ); % for deuterium: D/H ; for 20 C, it gives 1.0844
% alpha_18O = exp( -7.685/10^3 + 6.7123./T_oc - 1.6664*(10^3./T_oc.^2)...
%              + 0.35041*(10^6./T_oc.^3) ); % for oxygen: 18O/16O; for 20 C, it gives 1.0098
% % alpha is the temperature-dependent fractionation factor of the vapour to liquid phase transition.
% 
% dD_l = mean(dDp(delay_flag==0)); % from hydrometeor
% % Calculating Xp,eq (based on Equation 4 Graf et al. 2019)
% dD_p = (1000./alpha_D')  .*((dD_l  ./1000)+ 1 - alpha_D');
% 
% dD_eq_liq = dD_ob; % from hydrometeor??? or observations???
% dD_eq_vap = dD_p; % from the rainwater samples
% % dD_evap = cf.*dD_l + (1-cf).*dD_eq_vap;
% % dD_evap - dD_eq_vap = c.*[dD_l - dD_eq_liq];
% cf = (dD_evap' - dD_eq_vap)./(dD_l - dD_eq_liq);

%% Plots in dD-theta space %%
% figure;
% plot(th_ob,dD_ob,'sk','MarkerSize',5,'MarkerFaceColor','k')
for k = 4 %1:size(dD_ent,1)
    figure;
    plot(th_ob(iso_cold_pool_flag==flag),dD_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',5,'MarkerFaceColor','k')
%     plot(th_ob(iso_cold_pool_flag==flag),d18O_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',5,'MarkerFaceColor','k')
%     plot(th_ob,dD_ob,'sk','MarkerSize',5,'MarkerFaceColor','k')
    hold on;
    scatter(th_ob(iso_cold_pool_flag==flag),dD_ob(iso_cold_pool_flag==flag),13,q_ob(iso_cold_pool_flag==flag),'s','filled')
%     scatter(th_ob(iso_cold_pool_flag==flag),d18O_ob(iso_cold_pool_flag==flag),13,q_ob(iso_cold_pool_flag==flag),'s','filled')
%     scatter(th_ob,dD_ob,13,q_ob,'s','filled')
    colormap(flip(jet(12)))
    caxis([11 17])
    hdl=colorbar;
    ylabel(hdl,'q [g kg^-^1])','FontSize',16,'Rotation',90);
%     scatter(th_ent(k,iso_cold_pool_flag==flag),d18O_ent(k,iso_cold_pool_flag==flag),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(k,1:24:end),d18O_ent(k,1:24:end),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(th_ent(k,:),'omitnan'),mean(d18O_ent(k,:),'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(th_ent(k,iso_cold_pool_flag==flag),'omitnan'),mean(d18O_ent(k,iso_cold_pool_flag==flag),'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
    scatter(th_ent(k,iso_cold_pool_flag==flag),dD_ent(k,iso_cold_pool_flag==flag),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_surf(iso_cold_pool_flag==flag),d18O_surf(iso_cold_pool_flag==flag),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(1:24:end),d18O_surf(1:24:end),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(mean(th_surf,'omitnan'),mean(dD_surf,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
%     scatter(mean(th_surf(iso_cold_pool_flag==flag),'omitnan'),mean(dD_surf(iso_cold_pool_flag==flag),'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
    scatter(th_surf(iso_cold_pool_flag==flag),dD_surf(iso_cold_pool_flag==flag),15,[0.9290 0.6940 0.1250],'filled')
    scatter(th_dd(iso_cold_pool_flag==flag),dD_dd(k,iso_cold_pool_flag==flag),15,[0 0.4470 0.7410],'filled')
%     scatter(th_dd(iso_cold_pool_flag==flag),d18O_dd(k,iso_cold_pool_flag==flag),15,[0 0.4470 0.7410],'filled')
%     scatter(mean(th_dd(iso_cold_pool_flag==flag),'omitnan'),mean(dD_dd(k,iso_cold_pool_flag==flag),'omitnan'),50,[0 0.4470 0.7410],'filled')
    ylim([-92 -62]);
    xlim([290 302.5]);
%     ylabel(['DXS [',char(8240),']'])
    xlabel('\theta [K]')
    ylabel(['\deltaD [',char(8240),']'])
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength', [.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
%     hold on;text(-66,1,'surface','FontSize',15)
%     hold on;text(-75,21.5,'entrained','FontSize',15)
%     hold on;text(-60,11.5,'downdraft','FontSize',15)
%     caxis([297 300])
%     set(gca, 'YDir','reverse')
%     set(gca, 'XDir','reverse')
%     title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
%     text(-58,0,['RMSE = ',num2str(rmse(k)),''],'FontSize',15)
%     saveas(gcf,['RHB-soundings-cps-thvsdD-height-range',num2str(k)],'png')
%     saveas(gcf,['RHB-soundings-cps-thvsdD-height-range',num2str(k)],'fig')
end
% hold on;
% scatter(mean(dD_surf,'omitnan'),mean(th_surf,'omitnan'),50,'r','filled')
% legend('RB SBL obs','0.8-1km','0.9-1.1km','1-1.2km','1.1-1.3km','1.2-1.4km','surface','FontSize',15)
%% Plots in dD-DXS space %%
for k = 4 % 1:size(dD_ent,1)
    figure;
%     plot(dD_ob,de_ob,'sk','MarkerSize',5)
    plot(dD_ob(iso_cold_pool_flag==flag),de_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',5)
    hold on;
%     scatter(dD_ob,de_ob,13,th_ob,'s','filled')
    scatter(dD_ob(iso_cold_pool_flag==flag),de_ob(iso_cold_pool_flag==flag),13,th_ob(iso_cold_pool_flag==flag),'s','filled')
    colormap(jet(15))
    colorbar;
    scatter(mean(dD_ent(k,:),'omitnan'),mean(de_ent(k,:),'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
%     scatter(dD_ent(k,:),de_ent(k,:),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(dD_ent(k,iso_cold_pool_flag==flag),de_ent(k,iso_cold_pool_flag==flag),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(dD_dd(k,iso_cold_pool_flag==flag),'omitnan'),mean(de_dd(k,iso_cold_pool_flag==flag),'omitnan'),50,[0 0.4470 0.7410],'filled')
    scatter(dD_dd(k,:),de_dd(k,:),15,[0 0.4470 0.7410],'filled')
%     scatter(dD_dd(k,iso_cold_pool_flag==flag),de_dd(k,iso_cold_pool_flag==flag),15,[0 0.4470 0.7410],'filled')
    scatter(mean(dD_surf,'omitnan'),mean(de_surf,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
%     scatter(dD_surf,de_surf,15,[0.9290 0.6940 0.1250],'filled')
%     scatter(dD_surf(iso_cold_pool_flag==flag),de_surf(iso_cold_pool_flag==flag),15,[0.9290 0.6940 0.1250],'filled')
    xlim([-84 -56]) 
    ylim([-4 30])
    ylabel(['DXS [',char(8240),']'])
    xlabel(['\deltaD [',char(8240),']'])
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength', [.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
    hold on;text(-66,1,'surface','FontSize',15)
    hold on;text(-75,21.5,'entrained','FontSize',15)
    hold on;text(-60,11.5,'downdraft','FontSize',15)
    caxis([295 300])
    set(gca, 'YDir','reverse')
    set(gca, 'XDir','reverse')
    title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
%     saveas(gcf,['RHB-soundings-DXS-dD-dd-mean-height-range',num2str(k)],'png')
%     saveas(gcf,['RHB-soundings-DXS-dD-dd-mean-height-range',num2str(k)],'fig')
end
%% Plots in theta-q space %%
for k = 4 %1:size(dD_ent,1)
    figure;
    plot(th_ob(iso_cold_pool_flag==flag),q_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',4)
    hold on;
    scatter(th_ob(iso_cold_pool_flag==flag),q_ob(iso_cold_pool_flag==flag),13,dD_ob(iso_cold_pool_flag==flag),'s','filled')
    colormap(flip(jet(16)))
    hdl=colorbar;
    ylabel(hdl,['\deltaD [',char(8240),']'],'FontSize',16,'Rotation',90);
    scatter(th_ent(k,iso_cold_pool_flag==flag),q_ent(k,iso_cold_pool_flag==flag),10,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(k,1:24:end),q_ent(k,1:24:end),10,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(dD_ent(k,:),'omitnan'),mean(de_ent(k,:),'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
    scatter(th_surf(iso_cold_pool_flag==flag),q_surf(iso_cold_pool_flag==flag),10,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(1:24:end),q_surf(1:24:end),10,[0.9290 0.6940 0.1250],'filled')
%     scatter(mean(th_surf,'omitnan'),mean(q_surf,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
    scatter(th_dd(iso_cold_pool_flag==flag),q_dd(iso_cold_pool_flag==flag),10,[0 0.4470 0.7410],'filled')
%     scatter(mean(th_dd,'omitnan'),mean(q_dd,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
    xlim([288 304]); % xlim([288.8 302.2])
    ylim([7 23])
    ylabel('q [g/kg]')
    xlabel('\theta [K]')
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength', [.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
    hold on;text(294.25,21.5,'surface','FontSize',15)
    hold on;text(295,9,'entrained','FontSize',15)
    hold on;text(289.5,11,'downdraft','FontSize',15)
    caxis([-81 -65])
%     title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
%     saveas(gcf,['RHB-soundings-th-q-plain-height-range',num2str(k)],'png')
%     saveas(gcf,['RHB-soundings-th-q-plain-height-range',num2str(k)],'fig')
end
%% Plots to showcase sensitivity analysis %%
% 1st try : q-dD space
time_clips = [301:600;601:900;901:1200;1201:1500;1501:1800;1801:2100];
for k = 1:size(time_clips,1)
    figure;
    title(['From ',datestr(t(time_clips(k,1)),'mmm-dd'),' to ',datestr(t(time_clips(k,end)),'mmm-dd')])
    for c = 1:length(h_hi)
        hold on;
        plot(q_ent(c,time_clips(k,:)),dD_dd(c,time_clips(k,:)))
%         plot(q_dd(time_clips(k,:)),q_ent(c,time_clips(k,:)))
    end
    xlabel('q_e_n_t [g/kg]')
%     xlabel('q_d_d [g/kg]')
    ylabel(['\deltaD_d_d [',char(8240),']'])
%     ylabel(['DXS_e_n_t [',char(8240),']'])
    legend('0.8-1km','0.9-1.1km','1-1.2km','1.1-1.3km','1.2-1.4km')
end

figure;
title('Entire timeseries')
for c = 1:length(h_hi)
    hold on;
    plot(q_ent(c,:),de_ent(c,:))
end
xlabel('q_e_n_t [g/kg]')
% xlabel('q_d_d [g/kg]')
ylabel(['DXS_e_n_t [',char(8240),']'])
legend('0.8-1km','0.9-1.1km','1-1.2km','1.1-1.3km','1.2-1.4km')
%%
figure;
plot(t,fdd(4,:))
datetick('x','mmm-dd','keeplimits','keepticks')
