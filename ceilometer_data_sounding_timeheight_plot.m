%% Adding ceilometer data to sounding time-height plot %%
filename = 'EUREC4A_ATOMIC_RonBrown_Ceilometer_15s_20200109-20200212_v1.1.nc';
t_cb   = ncread(filename,'time')/3600/24 + datenum('20200101','yyyymmdd');
cb = ncread(filename,'first_cloud_base');
open('q_thw_4hr_soundings.fig')
hold on;
plot(t_cb,cb,'.k')

% 11-min running mean

% Raw mean
nanmean(cb); % = 1.2281 km % 1.228091452609265e+03 % in m
