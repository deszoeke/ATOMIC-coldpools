function [q_rayleigh,delt_rayleigh]= Rayleigh_curve_cond(q0,T0,R0)
% Reservoirs: water vapor
    Rvsmow = 155.76e-6; % unitless, deuterium
    frac = 1 - 0.5.^(0:10);    % increasing between 0 and 0.999
%     frac2 = 1 + 0.5.^(0:10); % decreasing between 2 and 1.001
%     frac1 = [frac2 flip(frac)]; % decreasing between 2 and 0   
        Rd    = 287.04; % [J K-1 kg-1] dry air
        Rv    = 461.50; % [J K-1 kg-1] moist air/vapor
        Cp    = 1005.7; % [J K-1 kg-1] specific heat dry air
        L     = 2.5e6; % [J/kg]=[m2/s2]  Latent Heat of Vaporization of Water
    addpath C:\Users\quinones\Documents\Data\thermo
    p0 = 1013.25*100; % pressure in Pa
    qsat = qs(p0,T0-273.15); % p [Pa], T [degrees C], qsat [kg/kg] (unitless)
    c = (Rd*T0 + L*qsat)./(L*(Rd/Rv) - Cp*T0);
    T = T0.*frac.^c; % in Kelvin
    alphae_d = exp( 1158.8e-12 .*T.^3 - 1620.1e-9 .*T.^2 + 794.84e-6 .*T - 161.04e-3 + 2.9992e6./T.^3 );
    q_rayleigh = q0.*frac;  
    rayleigh_step = R0.*frac.^(alphae_d-1);
    delt_rayleigh = 1e3*(rayleigh_step/Rvsmow - 1); % permil


