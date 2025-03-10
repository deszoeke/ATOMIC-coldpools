filename = 'EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
ncinfo(filename);
time_p = ncread(filename,'collection_time'); % time 
time_p = time_p/3600/24 + datenum('20200101','yyyymmdd');

delay_flag = ncread(filename,'delay_flag');
ind = 7:11; % indices of precip samples coinciding with iso data
d18O_p = ncread(filename,'d18O');
dD_p   = ncread(filename,'dD');

filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
time_T = ncread(filename,'time');
time_T = time_T/3600/24 + datenum('20200101','yyyymmdd');
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]

% Loading water vapor data [Xv]
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
time_v = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_v = datenum('01012020','mmddyyyy') + time_v*(1/(3600*24));

prec_mask = zeros(size(time_p)); % precip samples mask
for k = 1:length(time_p)-1
    if time_v(1)<=time_p(k)
        prec_mask(k) = find(time_v==time_p(k));
    end
end
prec_mask(prec_mask==0) = 1;
prec_mask(10) = prec_mask(10)+1;

%% Calculating equilibrium alphas based on Ta %%
factor = 24480;
T_oc = 273.15 + Ta(prec_mask+factor); % surface air temperature in Kelvin => from PSD instrument at (17m)
% T_oc = 273.15 + 20; % for sanity check; At 20 C, α is 1.0850 for 2H1H^16O/1H_2^16O and 
                                                      % 1.0098 for 1H_2^18O/1H_2^16O (Majoube, 1971b).
T_oc(isnan(T_oc)) = 273.15 + Ta(41342);
alpha_D = exp( 1158.8*(T_oc.^3./10^12) - 1620.1*(T_oc.^2./10^9)...
             + 794.84*(T_oc./10^6) - 161.04/10^3 + 2.9992*(10^6./T_oc.^3) ); % for deuterium: D/H ; for 20 C, it gives 1.0844
alpha_18O = exp( -7.685/10^3 + 6.7123./T_oc - 1.6664*(10^3./T_oc.^2)...
             + 0.35041*(10^6./T_oc.^3) ); % for oxygen: 18O/16O; for 20 C, it gives 1.0098
% alpha is the temperature-dependent fractionation factor of the vapour to liquid phase transition.

% Calculating Xp,eq (based on Equation 4 Graf et al. 2019)
d18O_eqv = (1000./alpha_18O').*((d18O_p./1000)+ 1 - alpha_18O');
  dD_eqv = (1000./alpha_D')  .*((dD_p  ./1000)+ 1 - alpha_D');
 DXS_eqv = dD_eqv - 8*d18O_eqv; % deuterium excess

% d18O_p = alpha_18O.*(d18O_p+1)-1; % input delta in unitless (mult. by 1000), not in permil "units"
% dD_p   = alpha_D.*(dD_p+1)-1;
DXS_p  = dD_p - 8.*d18O_p; % deuterium excess

%% Figure %%
hold on;
% yyaxis right
yyaxis left
scatter(time_p(ind),d18O_p(ind),'r','filled')
scatter(time_p(ind),d18O_eqv(ind),'b','filled')
% ylabel(['\deltaD_p [',char(8240),']'])
ylim([-12 -7])
% ylabel(['\deltaD_v [',char(8240),']'])
% xlim([datenum('Jan/26/2020 14:20') datenum('Feb/10/2020 23:00')])

hold on;
% yyaxis right
yyaxis left
scatter(time_p(ind),dD_p(ind),'r','filled')
scatter(time_p(ind),dD_eqv(ind),'b','filled')
% ylabel(['\delta^1^8O_p [',char(8240),']'])
ylim([-82 -52])
% ylabel(['\delta^1^8O_v [',char(8240),']'])
% xlim([datenum('Jan/26/2020 14:20') datenum('Feb/10/2020 23:00')])

hold on;
% yyaxis right
yyaxis left
scatter(time_p(ind),DXS_p(ind),'r','filled')
scatter(time_p(ind),DXS_eqv(ind),'b','filled')
% ylabel(['DXS_p [',char(8240),']'])
ylim([-2 17])
% yyaxis left
% ylabel(['DXS_v [',char(8240),']'])
% xlim([datenum('Jan/26/2020 14:20') datenum('Feb/10/2020 23:00')])

% zoom files
xlim([datenum('Feb/04/2020 14:20') datenum('Feb/10/2020 23:00')])

% zoomX2 files
xlim([datenum('Feb/10/2020 14:40') datenum('Feb/10/2020 18:00')])
