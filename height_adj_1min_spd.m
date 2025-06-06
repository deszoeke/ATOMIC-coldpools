function [th_ob,q_ob,time_height_adj] = height_adj_1min_spd(zrf)
% [th_ob,q_ob,time_height_adj] = height_adj_1min_spd(zrf)
% load observations and adjust with flux-gradient similarity to height zrf

k =  0.4; % von Karman constant
Rd = 287.0400;
Cp = 1.0057e+03;

% 10-min flux variables
filename10 = 'data/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
qstar10 = ncread(filename10,'qstar'); % specific humidity scaling parameter ['g/kg']
Tstar10 = ncread(filename10,'tstar'); % temperature scaling parameter [K]
slp10 = 100*ncread(filename10,'psealevel'); % atmospheric pressure at sea level [hPa --> Pa]
L10 = ncread(filename10,'MO_length'); % Monin Obukhov length scale [m] all negative values; valid range: [-500,0]
time10 = ncread(filename10,'time')/3600/24 + datenum('20200101','yyyymmdd'); % in 10-min resolution

% 1-min met variables
filename = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
%           data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc
qa = ncread(filename,'qair'); % air specific humidity from PSL RH at height ztq [g/kg]
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]
slp = 100*ncread(filename,'psealevel'); % atmospheric pressure at sea level [hPa --> Pa]
time_height_adj = ncread(filename,'time')/3600/24 + datenum('20200101','yyyymmdd'); % in 1-min resolution
rdir = ncread(filename,'rdir'); % original/"raw" variable
zm = ncread(filename,'ztq'); zm=zm(1);

% Using our ship_flag, not the one from the PSD surface data, and only for the variables it has an impact on: wind variables[???]
% bad wind direction starboard-aft half
maskout = rdir<-135 | rdir>45; % 1 = bad wind direction % spd fixed 20250606
Ta(maskout)= NaN;   % air temperature at 17m [in degrees C]
qa(maskout)= NaN;   % specific humidity; units: g/kg
% ^ sufficient
% qstar(ship==1)= NaN;    
% Tstar(ship==1)= NaN;    
% L(ship==1)= NaN;    
% slp(ship==1)= NaN;  % atmospheric pressure at sea level; units: mbar = hPa

g = 9.8; % [m/s^2]
H = Rd*(Ta+273.15)./g; % scale height [in m]
% Computing pressure at zm level using hypsometric equation
pm = slp.*exp(-zm./H); % [in Pa]

% interp thetastar, qstar, and L to 1 min
% thstar = interp1(time10, Tstar10.*(1e5./slp10).^(Rd/Cp) ,time_height_adj,'previous');
% qstar  = interp1(time10, qstar10                        ,time_height_adj,'previous');
% L      = interp1(time10, L10                            ,time_height_adj,'previous');
% smoother
dnmid = time10 + 5/60/24;
thstar = interp1(dnmid, Tstar10.*(1e5./slp10).^(Rd/Cp) ,time_height_adj,'linear');
qstar  = interp1(dnmid, qstar10                        ,time_height_adj,'linear');
L      = interp1(dnmid, L10                            ,time_height_adj,'linear');

% flux-gradient similarity adjustment
th_m = (Ta + 273.15).*(1e5./pm).^(Rd/Cp); % (Potential Temp in degrees K)
% th_ob = th_m + (thstar./k).*(log(zrf/zm)-psi_T(zrf./L)+psi_T(zm./L));
% q_ob  = qa + (qstar./k).*(log(zrf/zm)-psi_T(zrf./L)+psi_T(zm./L));
th_ob = th_m + adj_t(tstar,L, zm,zrf);
q_ob  = qa   + adj_t(qstar,L, zm,zrf);

end

function dt = adj_t(star,L, zm,zrf)
% dt = adj_t(star,L, zm,zrf)
% adjust thermal variable from measurement height zm to reference height zm
% Trf = Tm + adj_t(Tstar,L, zm,zrf)
    k=0.4;
    dt = (star./k).*(log(zrf/zm)-psi_T(zrf./L)+psi_T(zm./L));
end
