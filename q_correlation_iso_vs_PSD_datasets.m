% Correlation between the two measurements of q during ATOMIC %
% q_iso vs q_PSD: correlation = 0.9532!!!
% Loading datasets
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
t_iso = ncread(filename,'time'); % time of the met data in 'seconds since 2020-01-01 00:00:00UTC'
t_iso = t_iso/3600/24 + datenum('20200101','yyyymmdd');
q_isoo  = ncread(filename,'q');
q_iso = movmean(q_isoo,11,'omitnan'); % 11-min running average
figure; plot(t_iso,q_isoo)%,t_iso,q_iso)

load '1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' qair t1min % 11-min running averages
filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
t1mino = ncread(filename,'time'); % time of the met data in 'seconds since 2020-01-01 00:00:00UTC'
t1mino = t1mino/3600/24 + datenum('20200101','yyyymmdd');
qairo = ncread(filename,'qair'); % specific humidity; units: g/kg
ind = find(t1min>=t_iso(1));
ind = ind(1);
qairo = qairo(ind:ind+length(t_iso)-1); % qair raw for iso times
qair  = qair(ind:ind+length(t_iso)-1);  % 11-min running average for iso times
hold on; plot(t1mino(ind:ind+length(t_iso)-1),qairo)%,t1min(ind:ind+length(t_iso)-1),qair)

q_isoo(isnan(qairo)) = NaN;
qairo(isnan(q_isoo)) = NaN;

c = corr(qairo',q_isoo,'rows','complete');
% c = corr(qairo',q_isoo,'rows','pairwise');
% c = 0.9532 !!!