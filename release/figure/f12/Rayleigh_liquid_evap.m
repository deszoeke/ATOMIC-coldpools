function [q_lost,dD_lost,f_lost] = Rayleigh_liquid_evap(p,qv,T,dD0)
% p[Pa],qv[g/kg],T[K],dD0[permil]
% addpath C:\Users\quinones\Documents\Data\thermo
% p = 980 * 10^2; % near-surface pressure in Pa
% qv = 15; % units: [g/kg]; mass of water vapor per mass of dry air
% T = 294; % temp of dry air in Kelvin
    % temp that Simon uses: 294 K (still too warm for the SBL)
    % temp used initially:  298 K (too warm for the SBL)
% T0 = Twet_eqm(p); % use wetbulb temp (not theta) for liquid water [in K]
T0 = Twet(T-273.15,qv/10^3,p)+273.15; % (Tdry[C], qv[kg/kg], p[Pa])
% dD0 = 5; % dD of liquid water in permil
    % dD of liquid water Simon uses:     5 permil
    % dD used based on Adriana's paper: 10 permil
    Rvsmow = 155.76e-6; % unitless, deuterium
    Rd = 287.04; % [J K-1 kg-1] dry air
R0 = Rvsmow * (1 + dD0*1e-3); % use dD0 in permil
[~,f,R] = Rayleigh_curve_evap(T0,R0);
% q_rain = 2;% similar but not the same as a mixing ratio; like a specific humidity for rainwater;
             % seems confusing since q is only defined for water vapor!!!
             % is a very small value!!!
             % q_rain = rho_rain/(rho_air + rho_rain)
             % rho_air = rho_d + rho_v
R_lost = (R0 - f.*R)./(1 - f);
dD_lost = 1e3*(R_lost/Rvsmow - 1); % permil
Tv = T.*(1 + 0.61.*qv/10^3); % in Kelvin?; virtual temperature
    rho_rain = 2.64*10^-5; % units: [2.64*10^-5 kg/m^3] value provided by Simon's code from integrating approach
rho_air = p./(Rd.*Tv); % has to be ~1 [kg/m^3]
q_rain = rho_rain./(rho_air+rho_rain); % in kg/kg
%% !!! Changed units for q_lost !!! from kg/kg to g/kg !!!
q_lost = q_rain.*(1-f).*1e3; % in g/kg
f_lost = 1-f;