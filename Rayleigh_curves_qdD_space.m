%% {q,dD} pair distributions

% Loading iso data: RHB Picarro data
load 'iso_data_1min_intervals_FLAGGED_w_runningmean.mat' d18O dD DXS iso_time
% variables: d18O dD DXS factor iso_time

% Using q from the Picarro
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
q_iso = ncread(filename,'q'); % specific humidity; ratio of the mass of water vapor to the mass of moist air

% Using q from the PSD instruments
load '1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' qair t1min Ta
% variables: qair rh rr ship slp sst t1min Ta u v wdir wspd
ind0 = find(t1min==iso_time(1));
% Incorporating isotope data %
pos = (ind0:1:ind0+length(iso_time)-1)';
Ta = Ta(pos)'+273.15; % in degrees Kelvin
q_psd = qair(pos)'; % in g/kg
q_psd_ppm = q_psd*1000;
% filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';

%% Conversion from ppmv to g/kg
T = 25 + 273.15; % in Kelvin; assuming constant temp
rh = 80; % relative humidity in percent
P = 1013.25*100; % pressure in Pa
    Rd    = 287.04; % [J K-1 kg-1] dry air
    Rv    = 461.50; % [J K-1 kg-1] moist air/vapor
addpath C:\Users\quinones\Documents\Data\thermo
Pv = rh/100*es(T-273.15,P); %T [degrees C], P [Pa], es [Pa]
Pd = P - Pv;
rho_d = Pd/(Rd*T); % density of dry air in kg/m^3
rho_v = Pv/(Rv*T); % density of vapor in kg/m^3
% cf = rho_v/(10^3*rho_d); % conversion factor from ppmv to g/kg

%% Conversion from g/kg to ppmv
cf = (10*rho_d)/rho_v; % conversion factor from g/kg to ppmv

%% Rayleigh curve
% alpha equilibrium functions
% equilibrium fractionation ratio for deuterium alpha(T[Kelvin]) = [HDO]/[H2O]"
    Rvsmow = 155.76e-6; % unitless, deuterium
T0 = 25 + 273.15; % in Kelvin; assuming constant temp
rh0 = 80; % relative humidity in percent
q0 = (rh0/100 * 2.541e6 * exp(-5415.0 / T0) * 18/29)*1000; % in g/kg
% q0 = (rh0/100 * 2.541e6 * exp(-5415.0 / T0) * 18/29)*10^6; % in ppm
dD0 = -80; % in per mil
R0 = Rvsmow * (1 + dD0*1e-3); % dD0 in permil
q_rayleigh = 0:0.1:23; % in g/kg
% q_rayleigh = 0:100:20000; % in ppm
% fraction = q_iso./q0; % m/m0;
fraction = q_rayleigh./q0; % m/m0;
    Cp    = 1005.7; % [J K-1 kg-1] specific heat dry air
    L     = 2.5e6; % [J/kg]=[m2/s2]  Latent Heat of Vaporization of Water
addpath C:\Users\quinones\Documents\Data\thermo
p0 = 1013.25*100; % pressure in Pa
qsat = qs(p0,T0-273.15); % p [Pa], T [degrees C], qsat [kg/kg] (unitless)
c = (Rd*T0 + L*qsat)./(L*(Rd/Rv) - Cp*T0);
% temp changes
T = T0.*fraction.^c; % in Kelvin
alphae_d = exp( 1158.8e-12 .*T.^3 - 1620.1e-9 .*T.^2 + 794.84e-6 .*T - 161.04e-3 + 2.9992e6./T.^3 );
rayleigh_step = R0.*fraction.^(alphae_d-1);
delt_rayleigh = 1e3*(rayleigh_step/Rvsmow - 1); % permil

%% Plots
figure;
plot(q_iso*cf,dD,'.')
ylabel(['\deltaD [',char(8240),']'])
xlabel('q_P_i_c_a_r_r_o [ppmv]')
hold on; plot(q_rayleigh*cf,delt_rayleigh)
xlim([0 14500])   
ylim([-500 -30])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','Fontsize'),'FontSize',14)

figure;
plot(q_iso,dD,'.')
ylabel(['\deltaD [',char(8240),']'])
xlabel('q_P_i_c_a_r_r_o [gkg^-^1]')
hold on; plot(q_rayleigh,delt_rayleigh)
xlim([0 14.5])   
ylim([-500 -30])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','Fontsize'),'FontSize',14)
    
figure;
subplot(142)
    plot(q_iso*cf,dD,'.')
    ylabel(['\deltaD [',char(8240),']'])
    xlabel('q_P_i_c_a_r_r_o [ppmv]')
    hold on; plot(q_rayleigh*cf,delt_rayleigh)
subplot(144)
    plot(q_psd*cf,dD,'.')
    ylabel(['\deltaD [',char(8240),']'])
    xlabel('q_P_S_D [ppmv]')
    hold on; plot(q_rayleigh*cf,delt_rayleigh)
    legend('all','Rayleigh')
subplot(141)
    plot(q_iso,dD,'.')
    ylabel(['\deltaD [',char(8240),']'])
    xlabel('q_P_i_c_a_r_r_o [gkg^-^1]')
    hold on; plot(q_rayleigh,delt_rayleigh)
subplot(143)
    plot(q_psd,dD,'.')
    ylabel(['\deltaD [',char(8240),']'])
    xlabel('q_P_S_D [gkg^-^1]')
    hold on; plot(q_rayleigh,delt_rayleigh)
    ylim([-1000 0])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','Fontsize'),'FontSize',14)

%% Including Rayleigh curves in q-dD plots for indivisual cold pools
% Rayleigh curve for condensation from cold pool's peak (highest dD point)
%% CP#13 onset on 08-Feb-2020 10:24:00
q0 = 15.1966; % peak: 15.1966; surface: 21.6714;
T0 = 25 + 273.15; % peak: 25 + 273.15; surface: 27 + 273.15; % in Kelvin;
dD0 = -70.6094; % peak: -70.6094; surface: -64.9727; % in per mil
    Rvsmow = 155.76e-6; % unitless, deuterium
R0 = Rvsmow * (1 + dD0*1e-3); % dD0 in permil
[q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0,T0,R0);

hold on; plot(q_rayleigh,delt_rayleigh)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)

%% CP#20 onset on 09-Feb-2020 15:59:00
q0 = 21.615; % peak: 15.3526; surface: 21.615;
T0 = 26.9 + 273.15; % peak: 25.5 + 273.15; surface: 26.9 + 273.15; % in Kelvin;
dD0 = -64.9643; % peak: -68.5895; surface: -64.9643 % in per mil
    Rvsmow = 155.76e-6; % unitless, deuterium
R0 = Rvsmow * (1 + dD0*1e-3); % dD0 in permil
[q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0,T0,R0);

hold on; plot(q_rayleigh,delt_rayleigh)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)

% Rayleigh curve for evaporation from surface
