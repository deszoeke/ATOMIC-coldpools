filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
time_T = ncread(filename,'time');
time_T = time_T/3600/24 + datenum('20200101','yyyymmdd');
T_skin = ncread(filename,'tskin'); % liquid temperature at [???]m [in degrees C]
RH  = ncread(filename,'rhair')/100;
% T_L = ncread(filename,'tsea'); % liquid temperature at [???]m [in degrees C]
% T_L = mean(ncread(filename,'tskin'),'omitnan'); % liquid temperature at [???]m [in degrees C]

% Co-locating RH & T_L variable in time %
load '2nd_leg_sounding_data_10min_linear_interp.mat' t
pos_i = 999999*ones(size(t));
for l = 1:length(t)
    pos2 = find(time_T>=t(l)); % rounding up to closest isotope surface data point!!!
                                    % try rounding to nearest data point
    pos_i(l) = pos2(1);
    clearvars pos2
end
time_RH = time_T(pos_i);
h_ship = RH(pos_i);
T_L    = T_skin(pos_i);
% load h_ship.mat
% h_ob = 0.75; % input RH in decimal point; unitless
h_ob = h_ship; % input RH in decimal point; unitless

% Loading water vapor data [Xv] %
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
time_v = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_v = datenum('01012020','mmddyyyy') + time_v'*(1/(3600*24));
d18O_v = ncread(filename,'d18O');
dD_v   = ncread(filename,'dD');
ship_flag = ncread(filename,'ship_flag');
inlet_flag = ncread(filename,'inlet_flag');
% d18O_a = d18O_ob; % mean(d18O_v(ship_flag==0 & inlet_flag==0),'omitnan');
% dHDO_a   = dD_ob; % mean(dD_v(ship_flag==0 & inlet_flag==0),'omitnan');
load deltas_observed_endmembers.mat
d18O_a = d18O_ob'; % for mixing diagram purposes (DXS vs dD plot w/ endmembers) 
dHDO_a =   dD_ob'; % for mixing diagram purposes (DXS vs dD plot w/ endmembers) 

% Calculating d_E's %
[d18O_e, d_18OE_approx] = delta_E_eqm(h_ob,'18O/16O',d18O_a,T_L);
[dHDO_e, d_HDOE_approx] = delta_E_eqm(h_ob,  'D/H'  ,dHDO_a,T_L);

% figure;
% plot(d_epsHDO,1-h,'.r')
% hold on;
% plot(d_eps18,1-h,'.b')
% ylabel('1-h')
% xlabel('\Delta\epsilon')
% legend('HDO','^1^8O')

%% Solving for delta_a
h = 0.875; % 0:0.05:1; % 0.42; % 0.75; % 
n = 0.5; 
% n = 1/2 for an open water body under natural conditions
theta = 0.88; % scaling factor/weighting term
% 0.50 for evaporation in the eastern Mediterranean Sea
% 0.88 for the North American Great Lakes

% Diffusivity ratios from Hellmann and Harvey (2020)
T = (T_L+274.15)/100; % degrees K
DrHDO  = 0.98258 - (0.02546./T) + (0.02421./(T.^(5/2)));
Dr18O  = 0.96671 + (0.007406./(T.^(1/2))) - (0.004861./(T.^(3)));

nCdHDO = ((DrHDO).^(-n)) - 1; % C_d = molecular diffusivity ratio
nCd18O = ((Dr18O).^(-n)) - 1; % C_d = molecular diffusivity ratio  

d_epsHDO = (1-h)*theta*nCdHDO; % delta_epsilon; unitless
d_eps18O = (1-h)*theta*nCd18O; % delta_epsilon; unitless

% Calculating alphas %
% alpha is the temperature-dependent fractionation factor of the vapour to liquid phase transition
T_oc = (T_L+274.15); % degrees K
alpha_D = exp( 1158.8*(T_oc.^3./10^12) - 1620.1*(T_oc.^2./10^9)...
             + 794.84*(T_oc./10^6) - 161.04/10^3 + 2.9992*(10^6./T_oc.^3));
             % for deuterium: D/H ; for 20 C, it gives 1.0844
alpha_O = exp( -7.685/10^3 + 6.7123./T_oc - 1.6664*(10^3./T_oc.^2)...
             + 0.35041*(10^6./T_oc.^3) ); 
             % for oxygen: 18O/16O; for 20 C, it gives 1.0098
dHDO_L = 7; % in permil from Sebastian's Meteor Values
d18O_L = 1; % in permil from Sebastian's Meteor Values
% dHDO_L = 5; % from Gat 1996
% d18O_L = 1; % from Gat 1996
eps_star_O = 1 - 1./alpha_O; % unitless
eps_star_D = 1 - 1./alpha_D; % unitless

delta_aO = ((alpha_O) .* d18O_L - d18O_e.*((1-h)+d_eps18O) - eps_star_O.*10^3 - d_eps18O.*10^3)./(h);
delta_aD = ((alpha_D) .* dHDO_L - dHDO_e.*((1-h)+d_epsHDO) - eps_star_D.*10^3 - d_epsHDO.*10^3)./(h);
% delta_a0 = ((1./alpha_D) .* dHDO_L - dHDO_e.*((1-h)+0) - eps_star_D - 0)./(h);
DXS_a = delta_aD - 8*delta_aO;

%% Plot h versus h*deltas %%
offset = delta_aO - 57;

figure;
% plot(h_ship.*delta_aO,h_ship,'.b')
plot(h_ship.*offset,h_ship,'.b')
hold on;
plot(h_ship.*delta_aD,h_ship,'.r')
% legend('18O','HDO')
legend('18O w/offset','HDO')
xlabel('h*\delta')
ylabel('h')
title('Diffusivity ratios from Hellmann & Harvey (2020)')
text(-50,0.7,'offset = - 57 permil')
% title('Diffusivity ratios from  Merlivat (1978)')

%% Plot 1-h versus deltas %%
offset = delta_aO - 57;

figure;
plot(delta_aO,(1-h_ship),'.b')
% plot(offset,(1-h_ship),'.b')
hold on;
plot(delta_aD,(1-h_ship),'.r')
legend('18O','HDO')
% legend('18O w/offset','HDO')
xlabel('\delta')
ylabel('1 - h')
title('Diffusivity ratios from Hellmann & Harvey (2020)')
text(-68,0.3,'offset = - 57 permil')
