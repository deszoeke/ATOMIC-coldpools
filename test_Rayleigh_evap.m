% Testing viability of including Rayleigh evap. of liquid in q-dD diagram

% addpath(".\thermo")
p = 980 * 10^2; % near-surface pressure in Pa
qv = 15; % units: [g/kg]; mass of water vapor per mass of dry air
T = 294; % temp of dry air
    % temp that Simon uses: 294 K (still too warm for the SBL)
    % temp used initially:  298 K (too warm for the SBL)
% T0 = Twet_eqm(p); % use wetbulb temp (not theta) for liquid water [in K]
T0 = Twet(T-273.15,qv/10^3,p)+273.15; % (Tdry[C], qv[kg/kg], p[Pa])
dD0 = 5; % dD of liquid water in permil
    % dD of liquid water Simon uses:     5 permil
    % dD used based on Adriana's paper: 10 permil
    Rvsmow = 155.76e-6; % unitless, deuterium
    Rd = 287.04; % [J K-1 kg-1] dry air
R0 = Rvsmow * (1 + dD0*1e-3); % use dD0 in permil
[delt_rayleigh,f,R] = Rayleigh_curve_evap(T0,R0);
% q_rain = 2;% similar but not the same as a mixing ratio; like a specific humidity for rainwater;
             % seems confusing since q is only defined for water vapor!!!
             % is a very small value!!!
             % q_rain = rho_rain/(rho_air + rho_rain)
             % rho_air = rho_d + rho_v
R_lost = (R0 - f.*R)./(1 - f);
dD_lost = 1e3*(R_lost/Rvsmow - 1); % permil
Tv = T.*(1 + 0.61.*qv/10^3); % in Kelvin?; virtual temperature
rho_rain = 2.64*10^-5; % units: [2.64*10^-5 kg/kg or kg/m^3] value provided by Simon's code from integrating approach
rho_air = p./(Rd.*Tv); % has to be ~1 [kg/m^3]
q_rain = rho_rain./(rho_air+rho_rain);
q_lost = q_rain.*(1-f);

% figure;
% subplot(221);plot(1-f,dD_lost);xlabel('1-f');ylabel(['\deltaD_l_o_s_t [',char(8240),']'])
% subplot(222);plot(dD_lost,1-f);ylabel('1-f');xlabel(['\deltaD_l_o_s_t [',char(8240),']'])
% subplot(223);plot(f,dD_lost);xlabel('f');ylabel(['\deltaD_l_o_s_t [',char(8240),']'])
% subplot(224);plot(dD_lost,f);ylabel('f');xlabel(['\deltaD_l_o_s_t [',char(8240),']'])

figure;
plot(q_lost,dD_lost);xlabel('q_l_o_s_t');ylabel(['\deltaD_l_o_s_t [',char(8240),']'])
