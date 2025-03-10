load('RHS_Eq9_MJ79_variablesUPDATED.mat') % loads 10-min data!!!
load ('1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat');
clearvars -except alpha_e alpha_e_D del_oc del_oc_D sst slp

%% Cold pool times
load('cold_pool_flag_1min.mat')
flag = 0; % script outputs plot for outside of cold pool times 
% flag = 1; % script outputs cold pool plot

%% Air type: observed
% q and theta @ 400m
zrf = 400; % [in meters]; reference height
zm = 17;   % [in meters]; measurement level
Rd = 287.04;
Cp = 1005.7;
[th_ob,q_ob,t_adj] = height_adj(zrf,zm,Rd,Cp);

cold_pool_flag_10min = cold_pool_flag_1min(1:10:end);
% th_ob(cold_pool_flag_10min==flag) = NaN;
% q_ob(cold_pool_flag_10min==flag)  = NaN;

% addpath('C:\Users\estef\Documents\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files')
load('2nd_leg_sounding_data_10min_linear_interp.mat')

pos_adj = 999999*ones(size(t));
for l = 1:length(t)
%     if l == length(t)
%         break
%     end
    pos1 = find(t_adj>=t(l)); % rounding up to closest (in time) PSD surface data point!!!
                              % try rounding to nearest (in space) data point
    pos_adj(l) = pos1(1);
    clearvars pos1
end

q_ob = q_ob(pos_adj);
th_ob = th_ob(pos_adj);

t_10min = t;

%% Loading iso data
load('iso_data_1min_intervals_FLAGGED_w_runningmean.mat');
ind0 = find(t1min>=iso_time(1));
ind0 = ind0(1);
iso_cold_pool_flag = cold_pool_flag_1min(ind0:ind0+length(dD)-1);
% dD(iso_cold_pool_flag==flag)   = NaN;
% d18O(iso_cold_pool_flag==flag) = NaN;

% Incorporating isotope data %
pos_i = 999999*ones(size(t));
for l = 1:length(t)
%     if l == length(t)
%         break
%     end
    pos2 = find(iso_time>=t(l)); % rounding up to closest isotope surface data point!!!
                                 % try rounding to nearest data point
    pos_i(l) = pos2(1);
    clearvars pos2
end

iso_d18O = d18O(pos_i);
iso_dD = dD(pos_i);
iso_de = iso_dD - (8*iso_d18O); % deuterium excess

%% Loading and replacing previous sounding data
% % Going from 10 min to 4 hours for the next 3 air types
% load('full_214_soundings_Level2_h&p_same_size.mat'); % data at 4hr intervals
% % or load ('sounding_data_Level2_h&p_same_size.mat')
% % or load ('2hPa_radiosonde_data_ATOMIC.m')
% 
% % Using only 2nd leg sounding data
% ind = find(t>=t_10min(1));
% ind = ind(1);
% t = t(ind:end);
% cp_ind = find(t1min>=t(1));
% cp_ind = cp_ind(1);
% cp_end = cp_ind + 240*(length(t)-1);
% cold_pool_flag_4hr = cold_pool_flag_1min(cp_ind:240:cp_end);

%% Air type: surface
% q and theta @ 1m
filename = 'EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]

T = (mean(Ta,'omitnan')+273.15)./100; % in degrees Kelvin
D2H  = 0.98258 - (0.02546./T) + (0.02421./(T^(5/2)));
D18O = 0.96671 + (0.007406/(T^.5)) - (0.004861./(T^3));
h_iso = linspace(0.5,1,5); % range based on ??? data (from above molecular layer to 20m height)
alpha_k_O = (1/D18O)*(1-h_iso)+ h_iso; % alpha_k must be greater than 1; checked!!!
alpha_k_D = (1/D2H)*(1-h_iso)+ h_iso;  % alpha_k must be greater than 1; checked!!!
axa_D = alpha_e_D.*alpha_k_D; % 10-min data!!!
axa_O = alpha_e  .*alpha_k_O; % 10-min data!!!
% Solving for the isotopic concentration del_v0_D
delta_e_D = ((1./axa_D)*(1+del_oc_D) - 1) * 1000; % 10-min data!!!
delta_e_O = ((1./axa_O)*(1+del_oc) - 1) * 1000; % 10-min data!!!

surf_dD = delta_e_D(:,5);   % using values from h=1;
surf_d18O = delta_e_O(:,5); % using values from h=1;
surf_de = surf_dD - (8*surf_d18O); % deuterium excess

% slp(cold_pool_flag_1min==flag) = NaN;
% sst(cold_pool_flag_1min==flag) = NaN;
pos = 999999*ones(size(t_10min));
for l = 1:length(t_10min)
    pos1 = find(t1min>=t_10min(l)); % rounding up to closest (in time) PSD surface data point!!!
                              % try rounding to nearest (in space) data point
    pos(l) = pos1(1);
    clearvars pos1
end
SLP = slp(pos(1:end)); % [in hPa]
T = sst(pos(1:end)); % [in degrees C]
addpath('C:\Users\quinones\Documents\Data\thermo')
q_surf = qs(SLP'*100,T')*1e3; % in g/kg; slp must be in Pa and T in degrees C
q_surf = q_surf';
th_surf = (T + 273.15).*(1e5./(SLP*1e2)).^(Rd/Cp); % (Potential Temp in degrees K)

%% cold_pool_flag_10min for 2nd leg only %%
pos_2nd_leg = find(t_adj==t_10min(1));
iso_cold_pool_flag_10min = cold_pool_flag_10min(pos_2nd_leg:pos_2nd_leg+length(t_10min)-1);
% surf_dD(iso_cold_pool_flag_10min==flag) = NaN;
% surf_d18O(iso_cold_pool_flag_10min==flag) = NaN;
% surf_de(iso_cold_pool_flag_10min==flag) = NaN;

%% Air type: entrained
% q and theta @ 1km
q_ent  = q(h==1e3,:)*1e3; % check units for q_inth, looks like is cg/kg
th_ent = th(h==1e3,:);

% q_ent = mean(q(h>=800 & h<=1200,:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_ent = mean(th(h>=800 & h<=1200,:),'omitnan');

% q_ent = mean(q(h>=1000 & h<=1400,:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_ent = mean(th(h>=1000 & h<=1400,:),'omitnan');

% q_ent = mean(q(h>=1100 & h<=1400,:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_ent = mean(th(h>=1100 & h<=1400,:),'omitnan');

% q_ent = mean(q(h>=1100 & h<=1300,:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_ent = mean(th(h>=1100 & h<=1300,:),'omitnan');

[ent_dD,ent_d18O] = iso_estimates_from_centroid_approach(q_ob,iso_dD,iso_d18O,q_surf,surf_dD,surf_d18O,q_ent);
 ent_de = ent_dD - (8*ent_d18O); % deuterium excess

%% Air type: downdraft from mean cloud layer 
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
for k = 1:size(thw,2)
    if trade_inv(k) >= 1e3
        trade_index = find(h==trade_inv(k));
        th_d(k) = nanmean(thw(thw_index:trade_index,k));
    else
        th_d(k) = NaN; % trade inversion
    end
end
addpath('C:\Users\quinones\Documents\Data\thermo')
q_d = qs(1e3*1e2,th_d-273.15)*1e3; % q_d in g/kg    
%   qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
%   p [Pa], T [degrees C], qs [kg/kg]
%   theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008

%% dd_iso's calculated based on the mixing fractions obtained from soundings
[fen, fss, fdd] = th_q_to_mixfraction(th_ob,q_ob, th_ent,q_ent, th_surf,q_surf, th_d,q_d);

fdd(fdd<0)=NaN;
fss(fss<0)=NaN;
fen(fen<0)=NaN;

fdd_mask = fdd;
fss_mask = fss;
fen_mask = fen;

fdd_mask(iso_cold_pool_flag_10min==0) = NaN;
fss_mask(iso_cold_pool_flag_10min==0) = NaN;
fen_mask(iso_cold_pool_flag_10min==0) = NaN;

X = -1:.1:1;
Ydd = normpdf(X,mean(fdd,'omitnan'),std(fdd,'omitnan'));
Yss = normpdf(X,mean(fss,'omitnan'),std(fss,'omitnan'));
Yen = normpdf(X,mean(fen,'omitnan'),std(fen,'omitnan'));

%% Mixing fractions plot
time = t_10min;
% figure;
color = 'g'; %[0 0.4470 0.7410]; %[0.8500 0.3250 0.0980]; %'k'; %[0 0.4470 0.7410]; %
subplot(311)
    plot(time,fdd,'Color',color)
    hold on;
    plot(time,zeros(size(time)),'k')
%     plot(time,fdd_mask,'-r')
    ylabel('dd fraction')
    ylim([-0.2 0.8])
    datetick('x','mm/dd','keeplimits','keepticks')
%     legend('ambient','cold pools')
    legend('@1km','mean [0.8 1.2] km','mean [1.0 1.4] km','mean [1.1 1.4] km','mean [1.1 1.3] km')

subplot(312)
    plot(time,fen,'Color',color)
    hold on;
    plot(time,ones(size(time)),'--k')
%     plot(time,fen_mask,'-r')
    ylabel('ent fraction')
    ylim([0 1.5])
    datetick('x','mm/dd','keeplimits','keepticks')

subplot(313)
    plot(time,fss,'Color',color)
    hold on;
    plot(time,zeros(size(time)),'k')
%     plot(time,fss_mask,'-r')
    ylabel('srf fraction')
    ylim([-0.5 0.5])
    datetick('x','mm/dd','keeplimits','keepticks')

set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','Fontsize'),'FontSize',13)
set(findall(gcf,'-property','XLim'),'XLim',[datenum('01/28/2020','mm/dd/yyyy') datenum('02/11/2020','mm/dd/yyyy')])
%%
dd_dD = (q_ob.*iso_dD' - q_ent.*fen.*ent_dD - q_surf.*fss.*surf_dD')./(q_d.*fdd);
dd_d18O = (q_ob.*iso_d18O' - q_ent.*fen.*ent_d18O - q_surf.*fss.*surf_d18O')./(q_d.*fdd);
dd_de = dd_dD - (8*dd_d18O); % deuterium excess
% dd_d18O(dd_d18O>0) = NaN;
% dd_dD(dd_dD>0) = NaN;
% dd_de(dd_dD>0) = NaN;

%% Plots in dD-DXS space %%
figure;
plot(iso_dD(iso_cold_pool_flag_10min==flag),iso_de(iso_cold_pool_flag_10min==flag),'sk','MarkerSize',5)
hold on;
scatter(iso_dD(iso_cold_pool_flag_10min==flag),iso_de(iso_cold_pool_flag_10min==flag),13,th_ob(iso_cold_pool_flag_10min==flag),'s','filled')
colormap(jet(15))
colorbar;
hold on;
scatter(mean(ent_dD,'omitnan'),mean(ent_de,'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
hold on;
scatter(mean(dd_dD(iso_cold_pool_flag_10min==flag),'omitnan'),mean(dd_de(iso_cold_pool_flag_10min==flag),'omitnan'),50,[0 0.4470 0.7410],'filled')
hold on;
scatter(mean(surf_dD,'omitnan'),mean(surf_de,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
xlim([-85 -55]) 
ylim([-5 22])
ylabel(['DXS [',char(8240),']'])
xlabel(['\deltaD [',char(8240),']'])
set(findall(gcf,'-property','Fontsize'),'FontSize',30)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05 .05])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
box on
hold on;text(-69,1,'surface','FontSize',15)
hold on;text(-68,8,'entrained','FontSize',15)
hold on;text(-80,8,'downdraft','FontSize',15)
% legend('observed')
axis square
caxis([295 300])
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
grid on

%% Plotting isotope timeseries for entrained model only %%
% figure;
color = 'g'; %[0 0.4470 0.7410]; %[0.8500 0.3250 0.0980]; %'k'; %[0 0.4470 0.7410]; %
subplot(311)
plot(time,ent_d18O,'Color',color)
hold on;
datetick('x','mm/dd','keeplimits','keepticks')
ylabel(['\delta^1^8O_e_n_t [',char(8240),']'])
legend('@1km','mean [0.8 1.2] km','mean [1.0 1.4] km','mean [1.1 1.4] km','mean [1.1 1.3] km')

subplot(312)
plot(time,ent_dD,'Color',color)
hold on;
ylabel(['\deltaD_e_n_t [',char(8240),']'])
datetick('x','mm/dd','keeplimits','keepticks')

subplot(313)
plot(time,ent_de,'Color',color)
hold on;
ylabel(['DXS_e_n_t [',char(8240),']'])
datetick('x','mm/dd','keeplimits','keepticks')

set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','Fontsize'),'FontSize',13)

%% Plotting isotope timeseries for entrained model %%
figure;
color = [0.8500 0.3250 0.0980];
subplot(311)
plot(time,ent_d18O,'Color',color)
hold on;
plot(time,surf_d18O,'Color',[0.9290 0.6940 0.1250])
plot(time,iso_d18O,'k')
plot(time,dd_d18O,'Color',[0 0.4470 0.7410])
datetick('x','mm/dd','keeplimits','keepticks')
ylabel(['\delta^1^8O [',char(8240),']'])
ylim([-20 0])
legend('entrained','surface','observed','downdraft')

subplot(312)
plot(time,ent_dD,'Color',color)
hold on;
plot(time,surf_dD,'Color',[0.9290 0.6940 0.1250])
plot(time,iso_dD,'k')
plot(time,dd_dD,'Color',[0 0.4470 0.7410])
ylabel(['\deltaD [',char(8240),']'])
datetick('x','mm/dd','keeplimits','keepticks')

subplot(313)
plot(time,ent_de,'Color',color)
hold on;
plot(time,surf_de,'Color',[0.9290 0.6940 0.1250])
plot(time,iso_de,'k')
plot(time,dd_de,'Color',[0 0.4470 0.7410])
ylabel(['DXS [',char(8240),']'])
datetick('x','mm/dd','keeplimits','keepticks')

set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','Fontsize'),'FontSize',13)
