%% Plots following Gat 1996 %%
% Figure 4 => marine evaporation process in terms of the d-diagram
% Loading precip data [Xp]
filename = 'EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
time_p = ncread(filename,'collection_time'); % time 
time_p = time_p/3600/24 + datenum('20200101','yyyymmdd');
delay_flag = ncread(filename,'delay_flag');
% ind = 7:11; % indices of precip samples coinciding with iso data
d18O_p = ncread(filename,'d18O');
dD_p   = ncread(filename,'dD');
d18Ovsmow = 0;
dDvsmow = 0;
% Loading water vapor data [Xv]
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
time_v = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_v = datenum('01012020','mmddyyyy') + time_v*(1/(3600*24));
d18O_v = ncread(filename,'d18O');
dD_v   = ncread(filename,'dD');
ship_flag = ncread(filename,'ship_flag');
inlet_flag = ncread(filename,'inlet_flag');
d18O_a = mean(d18O_v(ship_flag==0 & inlet_flag==0),'omitnan');
dD_a   = mean(dD_v(ship_flag==0 & inlet_flag==0),'omitnan');
d18O_L = 1; % from Gat 1996
dD_L   = 5; % from Gat 1996
d18O_seqv = mean(d18O_surf,'omitnan'); % surface in equilibrium vapor
  dD_seqv = mean(dD_surf,'omitnan');   % surface in equilibrium vapor
d18O_peqv = d18O_eqv; % precip. in equilibrium vapor
  dD_peqv = dD_eqv; % precip. in equilibrium vapor

% Slope of the evaporation line: S_E
y = mx + b;
% epsilon_D = ;
% epsilon_18O = ;
h = .75; % mean relative humidity
m = (h*(dD_a-dD_L)+epsilon_D)/(h*(d18O_a-d18O_L)+epsilon_18O); % S_E

%%
figure; hold on;
plot(d18O_a,dD_a,'ok','MarkerfaceColor','black')
plot(d18O_p(delay_flag==0 & d18O_p<0),dD_p(delay_flag==0 & d18O_p<0),'ok')
plot(d18Ovsmow,dDvsmow,'ok','MarkerfaceColor','black')
text(0,0,'surface water')
plot(d18O_peqv(delay_flag==0 & d18O_p<0),dD_peqv(delay_flag==0 & d18O_p<0),'ok','MarkerSize',4); % evaporated raindrop in vapor equilibrium
plot(d18O_seqv,dD_seqv,'ok','MarkerfaceColor','black'); % evaporated raindrop in vapor equilibrium
plot(d18O_L,dD_L,'ok','MarkerfaceColor','black')

plot(da)
plot(delta_epsilon)
plot(dp); % precipitation
plot(d18Osurf,dDsurf)
xlabel(['\delta^1^8O [',char(8240),']'])
ylabel(['\deltaD [',char(8240),']'])
ylim([-90 30])
xlim([-16 3])
box on

%% Simon's suggested plot
% Plot with (horizontal) delta_D and delta_18O vs. 1-h (vertical), with an offset, 
% to be determined, to make the deltas of the two isotopes lie atop each other. 
% The slope with respect to 1-h is expected to be the same in the turbulent layer
% (around h=0.75), and diverge to the two equilibrium deltas (likewise offset on 
% these axes) at (1-h)=0. The difference of the divergence is due to the difference
% between the kinetic Δϵ's, and in fact this difference has to be in the offset to 
% make the lines of dD and d18O lie on each other at h~=0.75.
h = 0:0.09:1;
figure;
plot(dD_p,1-h,'.r')
hold on;
plot(d18O_p,1-h,'.b')

