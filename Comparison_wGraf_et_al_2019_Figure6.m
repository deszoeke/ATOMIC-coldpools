%% Comparison with Graf et al. 2019 Figure 6 (a) %%
% Loading precip samples [Xl = Xp]
filename = 'EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
ncinfo(filename);
time_p = ncread(filename,'collection_time'); % time 
time_p = time_p/3600/24 + datenum('20200101','yyyymmdd');
delay_flag = ncread(filename,'delay_flag');
% ind = 9:11; % indices of precip samples coinciding with iso data
ind = 7:11;

d18O_p = ncread(filename,'d18O');
dD_p   = ncread(filename,'dD');
DXS_p  = dD_p - 8*d18O_p; % deuterium excess

% Loading water vapor data [Xv]
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
dD_v = ncread(filename,'dD'); %
d18O_v = ncread(filename,'d18O'); %
time_v = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_v = datenum('01012020','mmddyyyy') + time_v*(1/(3600*24));
inlet_flag = ncread(filename,'inlet_flag'); %

filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
time_rdir = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_rdir = datenum('01012020','mmddyyyy') + time_rdir/3600/24;
rdir_o = ncread(filename,'rdir'); % original/"raw" variable
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]

% rdir = zeros(size(time_v));
% for k = 1:length(time_v)
%     rdir(k) = rdir_o(time_rdir==time_v(k));
% end
% ship_flag = zeros(size(rdir));
% ship_flag(rdir>-135 & rdir>45) = 1; % 1 = bad wind dir
% 
%   dD_v(ship_flag==1 | inlet_flag==1) = NaN;
% d18O_v(ship_flag==1 | inlet_flag==1) = NaN;
 
DXS_v  = dD_v - 8*d18O_v; % deuterium excess
 
prec_mask = zeros(size(time_p)); % precip samples mask
for k = 1:length(time_p)-1
    if time_v(1)<=time_p(k)
        prec_mask(k) = find(time_v==time_p(k));
    end
end
prec_mask(prec_mask==0) = 1;
prec_mask(10) = prec_mask(10)+1;

%% Calculating equilibrium alphas based on Ta
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
d18O_p_eq = (1000./alpha_18O').*((d18O_p./1000)+ 1 - alpha_18O');
  dD_p_eq = (1000./alpha_D')  .*((dD_p  ./1000)+ 1 - alpha_D');
 DXS_p_eq = dD_p_eq - 8*d18O_p_eq; % deuterium excess

% Calculating deltas
Delta_delta = dD_p_eq - dD_v(prec_mask);   % Deltaδ = δ^2H_p,eq − δ^2H_v,sfc (Equation 5 Graf et al. 2019)
    Delta_d = DXS_p_eq - DXS_v(prec_mask); % d_p,eq − d_v,sfc; (Equation 6 Graf et al. 2019)

%% Plot
figure;
plot([-22 19],[0 0],'-k','LineWidth',1)
hold on;
plot([0 0],[-20 6],'-k','LineWidth',1)
scatter(Delta_delta(prec_mask>1),Delta_d(prec_mask>1),45,time_p(prec_mask>1),'filled')
labels = {'02/06','02/07','02/08','02/09','02/10','02/11'}; % explicitly indicate the labels 
h = colorbar('Ticks',datenum('Feb/06/2020'):1:datenum('Feb/11/2020'),'TickLabels',labels);
colormap(jet(10))
caxis([datenum('Feb/06/2020') datenum('Feb/11/2020')])
xlim([-22 19])
ylim([-20 6])
% set(gca,'xticks',-20:5:15)
% set(gca,'yticks',-20:5:5)
text(3,-17,'Evaporation dominates')
text(-21,3,'Cloud signal dominates')
text(3,3,'Equilibration dominates')
ylabel(['\Deltad [',char(8240),']'])
xlabel(['\Delta\delta [',char(8240),']'])
box on
grid on