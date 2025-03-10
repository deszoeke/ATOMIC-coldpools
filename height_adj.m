function [th_ob,q_ob,time_height_adj] = height_adj(zrf,zm,Rd,Cp)
    k =  0.4; % von Karman constant
    filename = 'EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
    qa = ncread(filename,'qair'); % air specific humidity from PSL RH at height ztq [g/kg]
    Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]
    qstar = ncread(filename,'qstar'); % specific humidity scaling parameter ['g/kg']
    Tstar = ncread(filename,'tstar'); % temperature scaling parameter [K]
    slp = ncread(filename,'psealevel'); % atmospheric pressure at sea level [in mbar = hPa]
    slp = slp.*100; % in Pa
    L = ncread(filename,'MO_length'); % Monin Obukhov length scale [m] all negative values; valid range: [-500,0]
    time = ncread(filename,'time');
    time_height_adj = time/3600/24 + datenum('20200101','yyyymmdd');
    rdir = ncread(filename,'rdir'); % original/"raw" variable
    
    % Using our ship_flag, not the one from the PSD surface data, and only for the variables it has an impact on: wind variables[???]
    ship = zeros(size(rdir));
    ship(rdir>-135 & rdir>45) = 1; % 1 = bad wind direction
    Ta(ship==1)= NaN;   % air temperature at 17m [in degrees C]
    qa(ship==1)= NaN;   % specific humidity; units: g/kg
    qstar(ship==1)= NaN;    
    Tstar(ship==1)= NaN;    
    L(ship==1)= NaN;    
    slp(ship==1)= NaN;  % atmospheric pressure at sea level; units: mbar = hPa

    g = 9.8; % [m/s^2]
    H = Rd*(Ta+273.15)./g; % scale height [in m]
    % Computing pressure at zm level using hypsometric equation
    pm = slp.*exp(-zm./H); % [in Pa]
    
%   addpath('C:\Users\estef\Documents\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files')
    % Option #1 => Using 3 term equation on surface air temp
    th_m = (Ta + 273.15).*(1e5./pm).^(Rd/Cp); % (Potential Temp in degrees K)
    thstar = (Tstar).*(1e5./slp).^(Rd/Cp); % (Potential Temp in degrees K) % find reference!!!
    th_rf = th_m + (thstar./k).*(log(zrf/zm)-psi_T(zrf./L)+psi_T(zm./L));
    qrf = qa + (qstar./k).*(log(zrf/zm)-psi_T(zrf./L)+psi_T(zm./L));
    th_ob = th_rf; % (Potential Temp in degrees K)
    q_ob = qrf;

    % Option #2 => NOT using empirical functions (psi/phi) on surface air temp
    % th_rf = th_m' + (thstar./k).*(log(zrf/zm));
    % qrf = q_m' + (qstar./k).*(log(zrf/zm));
    % th_ob = th_rf; % (Potential Temp in degrees K)
    % q_ob = qrf;
    
end
