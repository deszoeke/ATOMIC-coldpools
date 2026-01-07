% Copyright 2025 Simon P. de Szoeke and Estefanía Quiñones 
% Meléndez.
% 
% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation 
% files (the “Software”), to deal in the Software without 
% restriction, including without limitation the rights to use, 
% copy, modify, merge, publish, distribute, sublicense, and/or 
% sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following 
% conditions:
%
% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.
%% cold pool detection algorith following Vogel et al. 2020 and Vogel (2017) %%
% Last updated on 02/15/2023
% By E. Quiñones Meléndez (quinones@oregonstate.edu)

% Using entire RHB air temp timeseries
% Loading 10-min data instead of 11-min running average
% filename = 'EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
% timef = ncread(filename,'time'); % time of the flux data
% % timef = timef(2575:4715)'; ----%%% 2nd leg data only %%%----
% timef = timef/3600/24 + datenum('20200101','yyyymmdd');

% Using entire RHB air temp timeseries
% Loading 1-min data
filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
t1min = ncread(filename,'time'); % time of the met data in 'seconds since 2020-01-01 00:00:00UTC'
t1min = t1min/3600/24 + datenum('20200101','yyyymmdd');
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]
rdir = ncread(filename,'rdir'); % original/"raw" variable
% Using our ship_flag, not the one from the PSD surface data, and only for the variables it has an impact on: wind variables[???]
ship = zeros(size(rdir));
ship(rdir<-135 & rdir>45) = 1; % 1 = bad wind direction
    % ship = ncread(filename,'ship_contamination_flag'); % ship contamination index for quality control
    %        % Value of 0 implies no significant maneuver nor bad wind direction during the average and is good data.
    %        % Interpolated from 10min data to 1min. This flag is a combination of several criteria including relative
    %        % wind direction (to eliminate unsuitable wind out of - 90/+90 deg sector), ship maneuvers (standard 
    %        % deviation of heading and ship's speed), and reasonable limits on certain other variables, such as ship motion correction.
rh = ncread(filename,'rhair'); % relative humidity; units: percentage
qair = ncread(filename,'qair'); % specific humidity; units: g/kg
rr = ncread(filename,'prate'); % rain rate; units: mm/hr
slp = ncread(filename,'psealevel'); % atmospheric pressure at sea level; units: mbar = hPa
wdir = ncread(filename,'wdir'); % true wind direction, deg
wspd = ncread(filename,'wspd'); % true wind speed; units: m/s
u = wspd.*cosd(270-wdir); % eastward wind speed, m/s
v = wspd.*sind(270-wdir); % northward wind speed, m/s
sst = ncread(filename,'tskin'); % skin-level sea surface temperature [in degrees C]

filename = 'EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
sh = ncread(filename,'hs_bulk'); % 10min resolution; sensible heat flux COARE 3.6 bulk model; surface_downward_sensible_heat_flux; units: W/m^2; range:[-52.3062,  11.8809]
lh = ncread(filename,'hl_bulk'); % 10min resolution; latent heat flux COARE 3.6 bulk model;   surface_downward_latent_heat_flux;   units: W/m^2; range:[-397.0974,-58.9983]
t10min = ncread(filename,'time'); %
t10min = t10min/3600/24 + datenum('20200101','yyyymmdd');
sh = interp1(t10min,sh,t1min); % interpolated to 1min
lh = interp1(t10min,lh,t1min); % interpolated to 1min

% tsea = ncread(filename,'tsea'); % sea water temperature [in degrees C] 
% lat = ncread(filename,'lat');
% lon = ncread(filename,'lon');

%% Using the ship flag
rh(ship==1)= NaN;   % relative humidity; units: percentage
qair(ship==1)= NaN; % specific humidity; units: g/kg
rr(ship==1)= NaN;   % rain rate; units: mm/hr
wdir(ship==1)= NaN; % true wind direction, deg
wspd(ship==1)= NaN; % wind speed; units: m/s
u(ship==1)= NaN;    % eastward wind speed, m/s
v(ship==1)= NaN;    % northward wind speed, m/s
sst(ship==1)= NaN;  % skin-level sea surface temperature [in degrees C]
slp(ship==1)= NaN;  % atmospheric pressure at sea level; units: mbar = hPa
sh(ship==1)= NaN;   % sensible heat flux COARE 3.6 bulk model; surface_downward_sensible_heat_flux; units: W/m^2; range:[-52.3062,  11.8809]
lh(ship==1)= NaN;   % latent heat flux COARE 3.6 bulk model;   surface_downward_latent_heat_flux;   units: W/m^2; range:[-397.0974,-58.9983]

%% Generating 11-min running averages

% qair = movmean(qair,11,'omitnan'); % qair filtered => 11-min running average
% wspd = movmean(wspd,11,'omitnan'); % wspd filtered => 11-min running average
% wdir = movmean(wdir,11,'omitnan'); % wdir filtered => 11-min running average
% rh = movmean(rh,11,'omitnan');     % rh filtered => 11-min running average
% rr = movmean(rr,11,'omitnan');     % prec filtered => 11-min running average
% u = movmean(u,11,'omitnan');       % u filtered => 11-min running average
% v = movmean(v,11,'omitnan');       % v filtered => 11-min running average
% sst = movmean(sst,11,'omitnan');   % skin-level sea surface temperature [in degrees C]
% slp = movmean(slp,11,'omitnan');   % atmospheric pressure at sea level; units: mbar = hPa
% sh = movmean(sh,11,'omitnan');     % sh filtered => 11-min running average
% lh = movmean(lh,11,'omitnan');     % lh filtered => 11-min running average

%% RUNNING COLD POOLS detection algorithm %%
[t_max,t_min,t_max_ind,t_min_ind,t_end,t_end_ind,~,~,delta_T,T_max,T_min,Taf] = cold_pool_detection_algorithm(t1min,Ta);
% Ta is filtered within the cold_pool_detection_algorithm function!

%%    FLAG for recovery times    %%
%           AND                   %
%  FLAG for peak cold pool times  %
% Defining peak cold pool times as the 5-minute period centered on the
% minimum temperature time (t_min)
recovery_flag_1min = zeros(size(t1min));
peak_flag_1min = zeros(size(t1min));
for k = 1:length(t_max)
    ii  = t_min_ind(k):t_end_ind(k);
    iii = t_min_ind(k)-3:t_min_ind(k)+2;
    recovery_flag_1min(ii) = 1;
    peak_flag_1min(iii) = 1;
end
% recovery_flag_1min = 1 means value corresponds to a time flagged within a
% cold pool RECOVERY
% recovery_flag_1min = 0 means value corresponds to a time flagged within a
% cold pool FRONT

% for plotting purposes
figure;
plot(t1min(recovery_flag_1min==1),Ta(recovery_flag_1min==1),'.b')
hold on;
plot(t1min(recovery_flag_1min==0),Ta(recovery_flag_1min==0),'.k')
plot(t1min(peak_flag_1min==1),Ta(peak_flag_1min==1),'.r')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool recovery flag'); ylim([22 28])

% for plotting purposes
figure;
plot(t1min,Ta,'k')
hold on;
ylim([22 28])
yyaxis right
plot(t1min,recovery_flag_1min,'-r','Color','r')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool recovery flag'); ylim([-1 2])

figure;
plot(t1min,Ta,'k')
hold on;
ylim([22 28])
yyaxis right
plot(t1min,peak_flag_1min,'-r','Color','r')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool peak flag'); ylim([0 1])

% plotting all the temp that correspond to cold pool times
figure;
plot(t1min(cp_matrix(~isnan(cp_matrix))),Taf(cp_matrix(~isnan(cp_matrix))),'ok')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool temp')

2+2;
%% Identifying the strongest and weakest cold pools based on delta_T
% num_vector = 1:length(t_max); %
% date_vector = t_max;
% cp_vector = [delta_T(num_vector)', num_vector', date_vector'];
% sorted_cp = sortrows(cp_vector,1);

%% FLAG for cold pool times (indexes)
cold_pool_flag_1min = zeros(size(t1min));
cp_matrix = zeros(length(t_max),121); % matrix of cold pool indexes
t_max0_matrix = 9999*ones(length(t_max),121); % time
t_min0_matrix = 9999*ones(length(t_max),121); % time

for k = 1:length(t_max)
    ii = t_max_ind(k):t_end_ind(k);
    cp_matrix(k,1:length(ii)) = ii;
    cold_pool_flag_1min(t_max_ind(k):t_end_ind(k)) = 1;
    t_max0_matrix(k,1:length(ii)) = (t1min(ii)-t_max(k))./(t_min(k)-t_max(k));
    t_min0_matrix(k,1:length(ii)) = (t1min(ii)-t_min(k))./(t_end(k)-t_min(k));
end
cp_matrix(cp_matrix==0) = NaN;
t_max0_matrix(t_max0_matrix==9999) = NaN;
t_min0_matrix(t_min0_matrix==9999) = NaN;

% for plotting purposes
figure;
hold on;
plot(t1min,cold_pool_flag_1min,'-r','Color','r')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool flag'); ylim([-1 2])

% Plotting all the temp that correspond to cold pool times
figure;
plot(t1min(cp_matrix(~isnan(cp_matrix))),Taf(cp_matrix(~isnan(cp_matrix))),'ok')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool temp')
    
%% Interpolating to obtain a single normalized time vector
% For t_min @ 0 and t_end @ 1
t_norm    = 0:0.05:1;
% T_norm  = 0:0.10:3; % Simon says we don't need to normalize time!!!
Taf_norm  = 9999*ones(101,length(t_norm));
Ta_norm   = 9999*ones(101,length(t_norm));
prec_norm = 9999*ones(101,length(t_norm));
wdir_norm = 9999*ones(101,length(t_norm));
qair_norm = 9999*ones(101,length(t_norm));
rh_norm   = 9999*ones(101,length(t_norm));
wspd_norm = 9999*ones(101,length(t_norm));
u_norm    = 9999*ones(101,length(t_norm));
v_norm    = 9999*ones(101,length(t_norm));
sh_norm   = 9999*ones(101,length(t_norm));
lh_norm   = 9999*ones(101,length(t_norm));

Taf = Ta;

for k = 1:length(t_max) % 49:62%[1:68,70:length(t_max)] % [65:68,70:86] %
    Taf_norm(k,:)  = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),Taf(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    Ta_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),Ta(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    prec_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),rr(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    wdir_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),wdir(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    qair_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),qair(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    rh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),rh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    wspd_norm(k,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),wspd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    u_norm(k,:)    = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),u(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    v_norm(k,:)    = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),v(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    sh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),sh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    lh_norm(k,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),lh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
end
%%
% Taf_norm(69,:)=NaN;
% prec_norm(69,:)=NaN;
% qair_norm(69,:)=NaN;
% wdir_norm(69,:)=NaN;
% rh_norm(69,:)=NaN;
% wspd_norm(69,:)=NaN;
% u_norm(69,:)=NaN;
% v_norm(69,:)=NaN;

% For plotting purposes
figure;
for k = 1:length(t_max) % 49:62 % 65:86 %
    hold on;
    plot(t_norm,Taf_norm(k,:))
end
    ylabel('Temp [\circC]')
    title('101 RHB cold pools; t_m_i_n @ 0; recovery @ 1')
    xlabel('normalized time [by (t - t_m_i_n)/(t_r_e_c - t_m_i_n)]')
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on

%% For t_max @ 0 and t_min @ 1
Taf_norm2  = zeros(101,length(t_norm));
Ta_norm2   = zeros(101,length(t_norm));
prec_norm2 = zeros(101,length(t_norm));
qair_norm2 = zeros(101,length(t_norm));
rh_norm2   = zeros(101,length(t_norm));
wdir_norm2 = zeros(101,length(t_norm));
wspd_norm2 = zeros(101,length(t_norm));
u_norm2    = zeros(101,length(t_norm));
v_norm2    = zeros(101,length(t_norm));
sh_norm2   = zeros(101,length(t_norm));
lh_norm2   = zeros(101,length(t_norm));

for k = 1:length(t_max) %65:86
    Taf_norm2(k,:)  = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),Taf(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    Ta_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),Ta(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    prec_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),rr(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    qair_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),qair(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    rh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),rh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    wdir_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),wdir(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    wspd_norm2(k,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),wspd(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    u_norm2(k,:)    = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),u(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    v_norm2(k,:)    = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),v(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    sh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),sh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
    lh_norm2(k,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),lh(cp_matrix(k,~isnan(cp_matrix(k,:)))),t_norm);
end
% Taf_norm2(69,:)=NaN;
% prec_norm2(69,:)=NaN;
% qair_norm2(69,:)=NaN;
% wdir_norm2(69,:)=NaN;
% rh_norm2(69,:)=NaN;
% wspd_norm2(69,:)=NaN;
% u_norm2(69,:)=NaN;
% v_norm2(69,:)=NaN;

% For plotting purposes
figure;
for k = 1:length(t_max) % 49:62 % 65:86 % 1:length(t_max)
    hold on;
    plot(t_norm,Taf_norm2(k,:))
end
    ylabel('Temp [\circC]')
    title('101 RHB cold pools; t_f_r_o_n_t @ 0; t_m_i_n @ 1')
    xlabel('normalized time [by (t - t_f_r)/(t_m_i_n - t_f_r)]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on
    
%% Creating composite time series in normalized time
Taf_comp = [Taf_norm2 Taf_norm(:,2:end)];
Ta_comp  = [Ta_norm2 Ta_norm(:,2:end)];
prec_comp = [prec_norm2 prec_norm(:,2:end)];
qair_comp = [qair_norm2 qair_norm(:,2:end)];
rh_comp   = [rh_norm2 rh_norm(:,2:end)];
wdir_comp = [wdir_norm2 wdir_norm(:,2:end)];
wspd_comp = [wspd_norm2 wspd_norm(:,2:end)];
u_comp = [u_norm2 u_norm(:,2:end)];
v_comp = [v_norm2 v_norm(:,2:end)];
sh_comp   = [sh_norm2 sh_norm(:,2:end)];
lh_comp   = [lh_norm2 lh_norm(:,2:end)];

% t_norm2 = linspace(-1,0,21);
% t_norm = linspace(0,1,21); 
% t_comp = [t_norm2 t_norm(2:end)];
% t_comp = [-20:1:0,1:1.5:28,30];
t_comp = [linspace(-18.7,0,21) linspace(1,30.4,20)];

mean(t_min([65,67:68,70:86])-t_max([65,67:68,70:86]))*24*60; % 19.6667 min
mean(t_min-t_max)*24*60; % 18.6931 min
mean(t_end([65,67:68,70:86])-t_min([65,67:68,70:86]))*24*60; % 25.7143 min
mean(t_end-t_min)*24*60; % 30.3762 min

figure;
for k = 1:length(t_max) % 49:62 % 65:86 % 1:length(t_max)
    hold on;
    plot(t_comp,Taf_comp(k,:))
end
    ylabel('T_a [\circC]')
    title('101 RHB cold pools; t_f_r_o_n_t @ -20; t_m_i_n @ 0; t_r_e_c @ 30')
    xlabel('time [minutes]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on
    xlim([t_comp(1) t_comp(end)])
    
%% Incorporating background data for 40 min before and after t_min
bgwindow = 40; % min
bg_Taf   = zeros(21,bgwindow+1);
bg_Taf2  = zeros(21,bgwindow-10+1);
bg_Ta    = zeros(21,bgwindow+1);
bg_Ta2   = zeros(21,bgwindow-10+1);
bg_prec  = zeros(21,bgwindow+1);
bg_prec2 = zeros(21,bgwindow-10+1);
bg_wspd  = zeros(21,bgwindow+1);
bg_wspd2 = zeros(21,bgwindow-10+1);
bg_wdir  = zeros(21,bgwindow+1);
bg_wdir2 = zeros(21,bgwindow-10+1);
bg_qair  = zeros(21,bgwindow+1);
bg_qair2 = zeros(21,bgwindow-10+1);
bg_rh    = zeros(21,bgwindow+1);
bg_rh2   = zeros(21,bgwindow-10+1);
bg_u     = zeros(21,bgwindow+1);
bg_u2    = zeros(21,bgwindow-10+1);
bg_v     = zeros(21,bgwindow+1);
bg_v2    = zeros(21,bgwindow-10+1);
bg_sh    = zeros(21,bgwindow+1);
bg_sh2   = zeros(21,bgwindow-10+1);
bg_lh    = zeros(21,bgwindow+1);
bg_lh2   = zeros(21,bgwindow-10+1);
% The first and last cold pools must be removed, no full background data
% available for both
for k = 2:length(t_max)-1 % 49:62 % [65:68,70:86]
    bg_Taf(k,:)  = Taf(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_Taf2(k,:) = Taf(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_Ta(k,:)   = Ta(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_Ta2(k,:)  = Ta(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_prec(k,:) = rr(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_prec2(k,:)= rr(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_wspd(k,:) = wspd(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_wspd2(k,:)= wspd(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_wdir(k,:) = wdir(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_wdir2(k,:)= wdir(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_qair(k,:) = qair(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_qair2(k,:)= qair(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_rh(k,:)   = rh(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_rh2(k,:)  = rh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_u(k,:)    = u(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_u2(k,:)   = u(t_end_ind(k):t_end_ind(k)+bgwindow-10);
    bg_v(k,:)    = v(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_v2(k,:)   = v(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_sh(k,:)   = sh(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_sh2(k,:)  = sh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
    bg_lh(k,:)   = lh(t_max_ind(k)-bgwindow:t_max_ind(k));
    bg_lh2(k,:)  = lh(t_end_ind(k):t_end_ind(k)+bgwindow-10);    
end
% bg_Taf(5,:) = NaN;
% bg_Taf2(5,:)= NaN;
% bg_prec(5,:) = NaN;
% bg_prec2(5,:)= NaN;
% bg_qair(5,:) = NaN;
% bg_qair2(5,:)= NaN;
% bg_wspd(5,:) = NaN;
% bg_wspd2(5,:)= NaN;
% bg_wdir(5,:) = NaN;
% bg_wdir2(5,:)= NaN;
% bg_rh(5,:) = NaN;
% bg_rh2(5,:)= NaN;
% bg_u(5,:) = NaN;
% bg_u2(5,:)= NaN;
% bg_v(5,:) = NaN;
% bg_v2(5,:)= NaN;

Taf_comp2 = [bg_Taf(:,1:end-1) Taf_comp(1:end-1,:) bg_Taf2(:,2:end)];
Ta_comp2  = [bg_Ta(:,1:end-1) Ta_comp(1:end-1,:) bg_Ta2(:,2:end)];
prec_comp2 = [bg_prec(:,1:end-1) prec_comp(1:end-1,:) bg_prec2(:,2:end)];
qair_comp2 = [bg_qair(:,1:end-1) qair_comp(1:end-1,:) bg_qair2(:,2:end)];
wspd_comp2 = [bg_wspd(:,1:end-1) wspd_comp(1:end-1,:) bg_wspd2(:,2:end)];
rh_comp2   = [bg_rh(:,1:end-1) rh_comp(1:end-1,:) bg_rh2(:,2:end)];
wdir_comp2 = [bg_wdir(:,1:end-1) wdir_comp(1:end-1,:) bg_wdir2(:,2:end)];
u_comp2 = [bg_u(:,1:end-1) u_comp(1:end-1,:) bg_u2(:,2:end)];
v_comp2 = [bg_v(:,1:end-1) v_comp(1:end-1,:) bg_v2(:,2:end)];
sh_comp2   = [bg_sh(:,1:end-1) sh_comp(1:end-1,:) bg_sh2(:,2:end)];
lh_comp2   = [bg_lh(:,1:end-1) lh_comp(1:end-1,:) bg_lh2(:,2:end)];

t_comp2 = [linspace(-60,-18.7,41) t_comp(2:end-1) linspace(30.4,60,31)]; %[-40:0,1:1.5:28,30:40];

figure;
for k = 2:length(t_max)-1 % 49:62 % 65:86 % 1:length(t_max)
    hold on;
    plot(t_comp2,Taf_comp2(k,:),'-');% ':','Color',[0.8500 0.3250 0.0980])
end
    ylabel('T_a [\circC]')
    title('99 RHB cold pools; t_f_r_o_n_t @ -18.7; t_m_i_n @ 0; t_r_e_c @ 30.4')
    xlabel('time [minutes]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on
   
clearvars -except bgwindow t_comp t_max0_matrix t_min0_matrix cp_matrix t1min t_comp2 Taf_comp2 Ta_comp2 prec_comp2 qair_comp2 wspd_comp2 rh_comp2 sh_comp2 lh_comp2 wdir_comp2 qs_comp2 u_comp2 v_comp2 t_max t_min t_max_ind t_min_ind t_end t_end_ind Taf qair wspd rh prec u v
    
%% Incorporating background data for 60 min before and after t_min
% bgwindow = 40; % min
% bg_Taf  = zeros(21,bgwindow+1);
% bg_Taf2 = zeros(21,bgwindow*(3/4)+1);
% for k = [65:68,70:86]
%     bg_Taf(k-64,:) = Taf(t_max_ind(k)-bgwindow:t_max_ind(k));
%     bg_Taf2(k-64,:) = Taf(t_end_ind(k):t_end_ind(k)+bgwindow*(3/4));
% end
% bg_Taf(5,:) = NaN;
% bg_Taf2(5,:)= NaN;
% 
% Taf_comp3 = [bg_Taf(:,1:end-1) Taf_comp(65:86,:) bg_Taf2(:,2:end)];
% t_comp3 = [-60:0,1:1.5:28,30:60];
% 
% figure;
% for k = 65:86%1:length(t_max)
%     hold on;
%     plot(t_comp3,Taf_comp3(k-64,:),'-');% ':','Color',[0.8500 0.3250 0.0980])
% end
%     ylabel('Temp [\circC]')
%     title('21 RHB cold pools; t_f_r_o_n_t @ -20; t_m_i_n @ 0; t_r_e_c @ 30')
%     xlabel('time [minutes]');
%     set(findall(gcf,'-property','Fontsize'),'FontSize',20)
%     set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
%     box on
%     grid on

%% Normalizing time for the iso data
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
dD = ncread(filename,'dD'); %
d18O = ncread(filename,'d18O'); %
time = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time = datenum('01012020','mmddyyyy') + time*(1/(3600*24));
inlet_flag = ncread(filename,'inlet_flag'); %

filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
time_rr = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_rr = datenum('01012020','mmddyyyy') + time_rr/3600/24;
rdir_o = ncread(filename,'rdir'); % original/"raw" variable

rdir = zeros(size(time));
for k = 1:length(time)
    rdir(k) = rdir_o(time_rr==time(k));
end
ship_flag = zeros(size(rdir));
ship_flag(rdir<-135 & rdir>45) = 1; % 1 = bad wind dir
% ship_flag = ncread(filename,'ship_flag'); %

dD(ship_flag==1 | inlet_flag==1) = NaN;
d18O(ship_flag==1 | inlet_flag==1) = NaN;

clearvars filename inlet_flag ship_flag rdir rdir_o time_rr

% dD   = movmean(dD,11,'omitnan');   % dD filtered => 11-min running average % in per mil
% d18O = movmean(d18O,11,'omitnan'); % d18O filtered => 11-min running average % in per mil

% mr = movmean(mr_avg_1min,11,'omitnan');     % mr filtered => 11-min running average % in g/kg
% q_iso = (mr/1000)./(1+(mr/1000))*1e3;       % air specific humidity [g/kg] from iso data
DXS = dD - 8*d18O;

% timei(1) == t1min(24481);
factor = 24480; % conversion factor between met-sea flux data index and iso data index
iso_cp_matrix_trim = cp_matrix([65,67,69:78,82:84,86],:) - factor;
iso_cp_matrix = cp_matrix(65:86,:) - factor;

% for plotting purposes
figure;
plot(time(iso_cp_matrix_trim(~isnan(iso_cp_matrix_trim))),dD(iso_cp_matrix_trim(~isnan(iso_cp_matrix_trim))),'ok')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool dD')

figure;
for k = [65,67,69:78,82:84,86] % 2:length(t_max)-1 % 65:86 % 1:length(t_max)
    hold on;
    plot(t_max0_matrix(k,~isnan(cp_matrix(k,:))),dD(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))))
end
    ylabel('dD [permil]')
    title('16 RHB cold pools; t_f_r_o_n_t @ 0; t_m_i_n @ 1')
    xlabel('normalized time [by (t - t_f_r)/(t_m_i_n - t_f_r)]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on
    
%% Interpolating to the same normalized time vector
% For t_min @ 0 and t_end @ 1
t_norm = 0:0.05:1;
dD_norm = zeros(19,length(t_norm));
d18O_norm = zeros(19,length(t_norm));
DXS_norm = zeros(19,length(t_norm));
q_iso_norm = zeros(19,length(t_norm));
mr_norm   = zeros(19,length(t_norm));

for k = [65,67,69:78,82:84,86]
    dD_norm(k-64,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),dD(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
    d18O_norm(k-64,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),d18O(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
    DXS_norm(k-64,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),DXS(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
%     q_iso_norm(k-64,:) = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),q_iso(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
%     mr_norm(k-64,:)   = interp1(t_min0_matrix(k,~isnan(cp_matrix(k,:))),mr(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
end
% dD_norm(5,:)=NaN;
% d18O_norm(5,:)=NaN;
% DXS_norm(5,:)=NaN;
% q_iso_norm(5,:)=NaN;
% mr_norm(5,:)=NaN;

% For plotting purposes
figure;
for k = [65,67,69:78,82:84,86] % 1:length(t_max)
    hold on;
    plot(t_norm,dD_norm(k-64,:))
end
    ylabel('dD [permil]')
    title('16 RHB cold pools; t_m_i_n @ 0; recovery @ 1')
    xlabel('normalized time [by (t - t_m_i_n)/(t_r_e_c - t_m_i_n)]')
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on

% For t_max @ 0 and t_min @ 1
dD_norm2 = zeros(19,length(t_norm));
d18O_norm2 = zeros(19,length(t_norm));
DXS_norm2 = zeros(19,length(t_norm));
q_iso_norm2 = zeros(19,length(t_norm));
mr_norm2   = zeros(19,length(t_norm));

for k = [65,67,69:78,82:84,86]
    dD_norm2(k-64,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),dD(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
    d18O_norm2(k-64,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),d18O(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
    DXS_norm2(k-64,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),DXS(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
%     q_iso_norm2(k-64,:) = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),q_iso(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
%     mr_norm2(k-64,:)   = interp1(t_max0_matrix(k,~isnan(cp_matrix(k,:))),mr(iso_cp_matrix(k-64,~isnan(iso_cp_matrix(k-64,:)))),t_norm);
end
% dD_norm2(5,:)=NaN;
% d18O_norm2(5,:)=NaN;
% DXS_norm2(5,:)=NaN;
% q_iso_norm2(5,:)=NaN;
% mr_norm2(5,:)=NaN;

% For plotting purposes
figure;
for k = [65,67,69:78,82:84,86] % 1:length(t_max)
    hold on;
    plot(t_norm,dD_norm2(k-64,:))
end
    ylabel('dD [permil]')
    title('16 RHB cold pools; t_f_r_o_n_t @ 0; t_m_i_n @ 1')
    xlabel('normalized time [by (t - t_f_r)/(t_m_i_n - t_f_r)]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on

%% Composite dD time series in normalized time
dD_comp = [dD_norm2 dD_norm(:,2:end)];
d18O_comp = [d18O_norm2 d18O_norm(:,2:end)];
DXS_comp = [DXS_norm2 DXS_norm(:,2:end)];
% q_iso_comp = [q_iso_norm2 q_iso_norm(:,2:end)];
% mr_comp = [mr_norm2 mr_norm(:,2:end)];

figure;
for k = [65,67,69:78,82:84,86] % 1:length(t_max)
    hold on;
    plot(t_comp,dD_comp(k-64,:))
end
    ylabel('\deltaD [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    title('16 RHB cold pools; t_f_r_o_n_t @ -18.7; t_m_i_n @ 0; t_r_e_c @ 30.4')
    xlabel('time [minutes]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on
    xlim([t_comp(1) t_comp(end)])
    
%% Incorporating background iso data
bg_dD  = zeros(21,bgwindow+1);
% bg_dD2 = zeros(21,bgwindow/2+1); % for 40 min before and after t_min
bg_dD2 = zeros(21,bgwindow*(3/4)+1); % for 60 min before and after t_min
bg_d18O  = zeros(21,bgwindow+1);
bg_d18O2 = zeros(21,bgwindow*(3/4)+1); % for 60 min before and after t_min
bg_DXS  = zeros(21,bgwindow+1);
bg_DXS2 = zeros(21,bgwindow*(3/4)+1); % for 60 min before and after t_min
bg_q_iso  = zeros(21,bgwindow+1);
bg_q_iso2 = zeros(21,bgwindow*(3/4)+1); % for 60 min before and after t_min
bg_mr  = zeros(21,bgwindow+1);
bg_mr2 = zeros(21,bgwindow*(3/4)+1); % for 60 min before and after t_min

for k = [65,67,69:78,82:84,86]
    bg_dD(k-64,:) = dD(t_max_ind(k)-factor-bgwindow:t_max_ind(k)-factor);
%     bg_dD2(k-64,:) = dD(t_end_ind(k)-factor:t_end_ind(k)-factor+bgwindow/2); % for 40 min before and after t_min
    bg_dD2(k-64,:) = dD(t_end_ind(k)-factor:t_end_ind(k)-factor+bgwindow*(3/4)); % for 60 min before and after t_min
    bg_d18O(k-64,:) = d18O(t_max_ind(k)-factor-bgwindow:t_max_ind(k)-factor);
    bg_d18O2(k-64,:) = d18O(t_end_ind(k)-factor:t_end_ind(k)-factor+bgwindow*(3/4)); % for 60 min before and after t_min
    bg_DXS(k-64,:) = DXS(t_max_ind(k)-factor-bgwindow:t_max_ind(k)-factor);
    bg_DXS2(k-64,:) = DXS(t_end_ind(k)-factor:t_end_ind(k)-factor+bgwindow*(3/4)); % for 60 min before and after t_min
%     bg_mr(k-64,:) = mr(t_max_ind(k)-factor-bgwindow:t_max_ind(k)-factor);
%     bg_mr2(k-64,:) = mr(t_end_ind(k)-factor:t_end_ind(k)-factor+bgwindow*(3/4)); % for 60 min before and after t_min
%     bg_q_iso(k-64,:) = q_iso(t_max_ind(k)-factor-bgwindow:t_max_ind(k)-factor);
%     bg_q_iso2(k-64,:) = q_iso(t_end_ind(k)-factor:t_end_ind(k)-factor+bgwindow*(3/4)); % for 60 min before and after t_min
end

bg_dD2(bg_dD2==0) = NaN;
bg_d18O2(bg_d18O2==0)= NaN;
bg_DXS2(bg_DXS2==0)= NaN;
% bg_mr2(bg_mr2==0)= NaN;
% bg_q_iso2(bg_q_iso2==0)= NaN;

dD_comp2 = [bg_dD(:,1:end-1) dD_comp bg_dD2(:,2:end)];
d18O_comp2 = [bg_d18O(:,1:end-1) d18O_comp bg_d18O2(:,2:end)];
DXS_comp2 = [bg_DXS(:,1:end-1) DXS_comp bg_DXS2(:,2:end)];
% q_iso_comp2 = [bg_q_iso(:,1:end-1) q_iso_comp bg_q_iso2(:,2:end)];
% mr_comp2 = [bg_mr(:,1:end-1) mr_comp bg_mr2(:,2:end)];

dD_comp2(dD_comp2==0) = NaN;
d18O_comp2(d18O_comp2==0)= NaN;
DXS_comp2(DXS_comp2==0)= NaN;
% bg_mr2(bg_mr2==0)= NaN;
% bg_q_iso2(bg_q_iso2==0)= NaN;

% dD_comp3 = [bg_dD(:,1:end-1) dD_comp bg_dD2(:,2:end)];
% t_comp3 = t_comp2;
% Taf_comp3 = Taf_comp2;

figure;
for k = [65,67,69:78,82:84,86] % 1:length(t_max)
    hold on;
    plot(t_comp2,dD_comp2(k-64,:),'-');%':','Color',[.5 .5 .5])
end
    ylabel('\deltaD [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    title('16 RHB cold pools; t_f_r_o_n_t @ -18.7; t_m_i_n @ 0; t_r_e_c @ 30.4')
    xlabel('time [minutes]');
    set(findall(gcf,'-property','Fontsize'),'FontSize',20)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
    box on
    grid on

%% Subtracting the mean from each individual cold pool
for k = 1:19
    dev_dD(k,:) = dD_comp2(k,:)-nanmean(dD_comp2(k,:)); % deviation/anomaly from the mean
    dev_Taf(k,:) = Taf_comp2(k,:)-nanmean(Taf_comp2(k,:)); % deviation/anomaly from the mean
    dev_d18O(k,:) = d18O_comp2(k,:)-nanmean(d18O_comp2(k,:)); % deviation/anomaly from the mean
%     dev_q_iso(k,:) = q_iso_comp2(k,:)-nanmean(q_iso_comp2(k,:)); % deviation/anomaly from the mean
    dev_DXS(k,:) = DXS_comp2(k,:)-nanmean(DXS_comp2(k,:)); % deviation/anomaly from the mean
end
%%
orange = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
mustard = [0.9290 0.6940 0.1250];
green = [0.4196 0.5569 0.1373];

Taf_comp2 = Taf_comp2([65,67,69:78,82:84,86],:);
wspd_comp2 = wspd_comp2([65,67,69:78,82:84,86],:);
qair_comp2 = qair_comp2([65,67,69:78,82:84,86],:);
prec_comp2 = prec_comp2([65,67,69:78,82:84,86],:);
rh_comp2 = rh_comp2([65,67,69:78,82:84,86],:);
u_comp2 = u_comp2([65,67,69:78,82:84,86],:);
v_comp2 = v_comp2([65,67,69:78,82:84,86],:);
dD_comp2 = dD_comp2([1,3,5:14,18:20,22],:);
d18O_comp2 = d18O_comp2([1,3,5:14,18:20,22],:);
DXS_comp2 = DXS_comp2([1,3,5:14,18:20,22],:);

%%
figure;
subplot(521)
    plot(t_comp2,nanmean(Taf_comp2),'-','LineWidth',2,'Color',orange)
    hold on;
    plot(t_comp2,nanmean(Taf_comp2)+std(Taf_comp2,'omitnan')/sqrt(16),'-','Color',orange)
    plot(t_comp2,median(Taf_comp2,'omitnan'),'-o','Color',orange)
    plot(t_comp2,nanmean(Taf_comp2)-std(Taf_comp2,'omitnan')/sqrt(16),'-','Color',orange)
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    ylabel('T_a [\circC]')
    xlim([-60 60])
    ylim([24.5 26.3])
subplot(522)
    plot(t_comp2,nanmean(wspd_comp2),'Color',green,'LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(wspd_comp2)+std(wspd_comp2,'omitnan')/sqrt(16),'Color',green)
    plot(t_comp2,median(wspd_comp2,'omitnan'),'-o','Color',green)
    plot(t_comp2,nanmean(wspd_comp2)-std(wspd_comp2,'omitnan')/sqrt(16),'Color',green)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('U [m/s]')
    box off
    xlim([-60 60])
    ylim([8 11.1])
subplot(523)
%     plot(t_comp2,nanmean(q_iso_comp2),'-sk')
    hold on;
    plot(t_comp2,nanmean(qair_comp2),'-','LineWidth',2,'Color',blue)
    plot(t_comp2,nanmean(qair_comp2)+std(qair_comp2,'omitnan')/sqrt(16),'-','Color',blue)
    plot(t_comp2,median(qair_comp2,'omitnan'),'-o','Color',blue)
    plot(t_comp2,nanmean(qair_comp2)-std(qair_comp2,'omitnan')/sqrt(16),'-','Color',blue)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('q [g/kg]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([14.2 15.4])
%     legend('Picarro','PSD')
subplot(524)
    plot(t_comp2,nanmean(dD_comp2),'-k','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(dD_comp2)+std(dD_comp2,'omitnan')/sqrt(16),'-k')
    plot(t_comp2,median(dD_comp2,'omitnan'),'-ok')
    plot(t_comp2,nanmean(dD_comp2)-std(dD_comp2,'omitnan')/sqrt(16),'-k')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('\deltaD [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-73.3 -69.7])
subplot(525)
    plot(t_comp2,nanmean(prec_comp2),'-r','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(prec_comp2)+std(prec_comp2,'omitnan')/sqrt(16),'-r')
    plot(t_comp2,median(prec_comp2,'omitnan'),'-or')
    plot(t_comp2,nanmean(prec_comp2)-std(prec_comp2,'omitnan')/sqrt(16),'-r')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    box off
    ylabel('rain [mm/hr]')
    xlim([-60 60])
    ylim([0 2.5])
    set(gca,'xticklabels',[])
subplot(526)
    plot(t_comp2,nanmean(d18O_comp2),'-m','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(d18O_comp2)+std(d18O_comp2,'omitnan')/sqrt(16),'-m')
    plot(t_comp2,median(d18O_comp2,'omitnan'),'-om')
    plot(t_comp2,nanmean(d18O_comp2)-std(d18O_comp2,'omitnan')/sqrt(16),'-m')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('\delta^1^8O [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    set(gca,'xticklabels',[])
    xlim([-60 60])
    ylim([-10.6 -10.1])
subplot(527)
    plot(t_comp2,nanmean(rh_comp2),'-b','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(rh_comp2)+std(rh_comp2,'omitnan')/sqrt(16),'-b')
    plot(t_comp2,median(rh_comp2,'omitnan'),'-ob')
    plot(t_comp2,nanmean(rh_comp2)-std(rh_comp2,'omitnan')/sqrt(16),'-b')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    set(gca,'xticklabels',[])
    ylabel('RH [%]')
    box off
    xlim([-60 60])
    ylim([68 80])
subplot(528)
    plot(t_comp2,nanmean(DXS_comp2),'-','LineWidth',2,'Color',mustard)
    hold on;
    plot(t_comp2,nanmean(DXS_comp2)+std(DXS_comp2,'omitnan')/sqrt(16),'-','Color',mustard)
    plot(t_comp2,median(DXS_comp2,'omitnan'),'-o','Color',mustard)
    plot(t_comp2,nanmean(DXS_comp2)-std(DXS_comp2,'omitnan')/sqrt(16),'-','Color',mustard)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    set(gca,'xticklabels',[])
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    xlim([-60 60])
    ylim([10.4 11.9])
    box off
    set(findall(gcf,'-property','Fontsize'),'FontSize',16)
% Including wind components %
subplot(529)
    plot(t_comp2,nanmean(u_comp2),'-c','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(u_comp2)+std(u_comp2,'omitnan')/sqrt(16),'-c')
    plot(t_comp2,median(u_comp2,'omitnan'),'-oc')
    plot(t_comp2,nanmean(u_comp2)-std(u_comp2,'omitnan')/sqrt(16),'-c')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('E-W wind [m/s]')
    box off
    xlim([-60 60])
    ylim([-10.5 -7])
    xlabel('time [minutes]')
subplot(5,2,10)
    plot(t_comp2,nanmean(v_comp2),'-g','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(v_comp2)+std(v_comp2,'omitnan')/sqrt(16),'-g')
    plot(t_comp2,median(v_comp2,'omitnan'),'-og')
    plot(t_comp2,nanmean(v_comp2)-std(v_comp2,'omitnan')/sqrt(16),'-g')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('N-S wind [m/s]')
    xlabel('time [minutes]')
    xlim([-60 60])
    ylim([-4.2 -1.6])
    box off
    set(findall(gcf,'-property','Fontsize'),'FontSize',16)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])
    
%% Plots
addpath('C:\Users\estef\Documents\MATLAB\addaxis6')
orange = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
green = [0.4196 0.5569 0.1373];
mustard = [0.9290 0.6940 0.1250];

figure;
% subplot(211)
    plot(t_comp2,nanmean(DXS_comp2),'-k','LineWidth',2)
    hold on;
    wdir(ship==1) = NaN;
    plot(t1min,wdir,'-m','LineWidth',2);%,'Color',[0.5 0.5 0.5]) 
    ylim([10 120])
%     plot(t_comp2,nanmean(dD_comp2)+std(dD_comp2,'omitnan'),'-k')
%     plot(t_comp2,median(dD_comp2,'omitnan'),'-ok')
%     plot(t_comp2,nanmean(dD_comp2)-std(dD_comp2,'omitnan'),'-k')    
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box on
    grid on
% subplot(212)
    addaxis(t_comp2,nanmean(q_iso_comp2),'-','LineWidth',2,'Color',blue)
%     hold on;
%     plot(t_comp2,nanmean(Taf_comp2)+std(Taf_comp2,'omitnan'),'-','Color',orange)
%     plot(t_comp2,nanmean(Taf_comp2)-std(Taf_comp2,'omitnan'),'-','Color',orange)
%     plot(t_comp2,median(Taf_comp2,'omitnan'),'-o','Color',orange)
    addaxislabel(2,'specific humiity [g/kg]')
    xlabel('time [minutes]')
    box on
    grid on
%     ylim([24.9 26.05])
    addaxis(t_comp2,nanmean(Taf_comp2)*-1,'-','LineWidth',2,'Color',orange)
    addaxislabel(3,'(-1)*T_a [\circC]')
% legend('mean')%,'std','median')
set(findall(gcf,'-property','Fontsize'),'FontSize',26)

%% Timeseries plot for the WIW poster %%
figure;
subplot(212)
    yyaxis right
    plot(t_avg_1min,dD,'-k','LineWidth',2)
    hold on;
    plot(t1min,cold_pool_flag_1min*-80,'-r')
    xticklabels([29:1:31,1:11])
    xlim([datenum('Jan/29/2020') datenum('Feb/11/2020')])
    ylim([-80 -65])
    box on
    grid on
    xlabel('time [minutes]')
    ylabel('\deltaD [^{\fontsize{10}o}/{\fontsize{10}oo}]')
subplot(211)
    plot(t1min,Taf,'-','LineWidth',2,'Color',orange)
    hold on;
    plot(t1min,cold_pool_flag_1min*28,'-k')
    ylim([23 28])
    set(gca,'xticklabels',[])
    xlim([datenum('Jan/29/2020') datenum('Feb/11/2020')])
    ylabel('T_a [\circC]')
    box on
    grid on
set(findall(gcf,'-property','Fontsize'),'FontSize',20)

%% Identifying percentage of cold pools with dD increase
dD_up_flag = 9999*ones(1,size(dD_comp,1));
dD_up_t = 9999*ones(1,size(dD_comp,1));
dD_up_ind = 9999*ones(1,size(dD_comp,1));

for k = 1:size(dD_comp,1)
    dD_base = dD(t_max_ind(k+64)-factor);
    for ii = t_min_ind(k+64)-factor%:1:t_end_ind(k+64)-factor
        if dD(ii) > dD_base
            dD_up_t(k) = t1min(ii+factor);
            dD_up_ind(k) = ii;
            dD_up_flag(k) = 1; % 1 = during cold pool times dD increases "significantly"
            break
        end
        dD_up_t(k) = t1min(ii+factor);
        dD_up_ind(k) = ii;
        dD_up_flag(k) = 0; % 0 = during cold pool times dD does NOT increases "significantly"
    end
end

%% Criterion for flagging dD INCREASE
del_dD = dD(2:end)-dD(1:end-1);
figure;
histogram(del_dD)
xlim([-0.4 0.4])
xlabel('\DeltadD [change in dD]')
cand = find(del_D<-0.05); % candidates positions

for k = 1:size(dD_comp,1)
    test(k) = dD_comp(k,21) - dD_comp(k,1);
end
figure;
histogram(test,'BinEdges',[-3:1:8])
xlabel('\DeltadD [dD(t_m_i_n) - dD(t_m_a_x)]')
ylabel('# of cold pools')
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
grid on

%% Vogel paper plots
figure;
subplot(211)
plot(t1min,Ta,'-k','LineWidth',.5)
hold on;
plot(t1min,Taf,'-k','LineWidth',1.5)
datetick('x','HHAM mm/dd','keeplimits')
ylabel('T [\circC]')
legend('1 min data','11-point centered running mean')

subplot(212)
plot(t1min(2:end),del_T,'-k','LineWidth',1.5)
hold on;
plot(t1min,-.05*ones(size(t1min)),'--k')
hold on;
plot(t1min,nanmean(del_T/10)*ones(size(t1min)),'-k','LineWidth',.5)
datetick('x','HHAM mm/dd','keeplimits')
ylabel('\deltaT [\circC]')
xlabel('time [UTC]')
