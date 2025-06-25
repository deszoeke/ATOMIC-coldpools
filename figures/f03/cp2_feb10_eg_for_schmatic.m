%% Ron H. Brown isotope analizer timeseries
filename = 'data/EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
dD = ncread(filename,'dD'); %
d18O = ncread(filename,'d18O'); %
dD_o = dD; % original/"raw" variable
d18O_o = dD; % original/"raw" variable
time = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time = datenum('01012020','mmddyyyy') + time*(1/(3600*24));
dxs = dD - 8*d18O;
dxs_o = dD_o - 8*d18O_o;
inlet_flag = ncread(filename,'inlet_flag'); %
dD(inlet_flag==1)   = NaN; %-------| poor
d18O(inlet_flag==1) = NaN; %       | data
dxs(inlet_flag==1)  = NaN; %-------| removed

dDs = movmean(dD,11,'omitnan');
d18Os = movmean(d18O,11,'omitnan');
dxss = dDs - 8*d18Os;

%% Met data
filename = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
tair = ncread(filename,'tair'); % C
qair = ncread(filename,'qair'); % kg/kg
time_rr = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_rr = datenum('01012020','mmddyyyy') + time_rr*(1/(3600*24));
rr_o = ncread(filename,'prate'); % original/"raw" variable; rain rate in mm/hr
% We need rr in the same size as the iso data
% Removing rain data outside of iso data times
rr = zeros(size(time));
for k = 1:length(time)
    rr(k) = rr_o(time_rr==time(k));
end
% rr(ship_flag==1) = NaN;
rr(inlet_flag==1) = NaN;

%% Precip data
filename = 'data/EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
dDp = ncread(filename,'dD'); %
d18Op = ncread(filename,'d18O'); %
ctime = ncread(filename,'collection_time'); %
ctime = datenum('01012020','mmddyyyy') + ctime*(1/(3600*24));
delay_flag = ncread(filename,'delay_flag'); %
dxsp = dDp - 8*d18Op;

%%
dD_rr = dDs; %-------|
d18O_rr = d18Os; %   | filtered
dxs_rr = dxss; %-----| data

dD_rr(isnan(rr)) = NaN;
d18O_rr(isnan(rr)) = NaN;
dxs_rr(isnan(rr)) = NaN;

dD_rr(rr<=0.1) = NaN;
d18O_rr(rr<=0.1) = NaN;
dxs_rr(rr<=0.1) = NaN;

tmin = datenum(2020,2,10,16,20,0); % 12, 17, or 20, depending on the metric
% t0   = tmin - 19.7/(24*60);
% tend = tmin + 31.5/(24*60);
t0   = tmin - 21/(24*60); % from 11-min filtered T
tend = tmin + 21/(24*60);

% plot vertical lines
% slicex = @(x) plot([x,x], gca.YLim, 'k-')
slicet = @(T) arrayfun(@(x) plot([x,x], gca().YLim, 'k-'), T);

fsz = 18;

%% Iso timeseries plot
figure;

clf;
ax = vert_axes_stack(5);
axes(ax(1)); hold on;
plot(ax(1), time_rr, tair, 'k.', 'markersize',12); ylim([22, 26]);
slicet([t0, tmin, tend])
ylabel(ax(1), 'T_{air} (\circ{C})')
axes(ax(2)); hold on;
plot(ax(2), time_rr, qair, 'k.', 'markersize',12); ylim([13, 16]);
slicet([t0, tmin, tend])
ylabel(ax(2), 'q_{air} (g/kg)')
axes(ax(3)); hold on;
plot(ax(3), time, dD,'k.','MarkerSize',12); ylim([-71, -64]);
slicet([t0, tmin, tend])
ylabel(ax(3), ['\delta{D} (', char(8240) ')'] )
yticks(ax(3), -71:-64);
ytl = get(ax(3),'yticklabel');
set(ax(3), 'yticklabel',{'','-70','','','','','-65',''})
axes(ax(4)); hold on;
bar(ax(4), time_rr,rr_o); ylim([0 15]);
slicet([t0, tmin, tend])
ylabel(ax(4), 'rain (mm/hr)')
xlabel('2020 Feb 10')
text(ax(1), t0  , 27, "\it{t}_0", 'fontsize',fsz, 'clipping','off')
text(ax(1), tmin, 27, "\it{t}_{min}", 'fontsize',fsz, 'clipping','off')
text(ax(1), tend, 27, "\it{t}_{end}", 'fontsize',fsz, 'clipping','off')

for i = 1:4
    a = ax(i);
    a.FontSize = fsz;
    a.Box = true;
    xlim(a, [datenum(2020,2,10,15,40,0), datenum(2020,2,10,17,10,0)])
    datetick(a, 'x', 'HH:MM', 'keeplimits')
end
for i = 1:3
    a = ax(i);
    a.XTickLabel={''};
end
text(ax(1), gca().XLim(1)+2/(24*60),  25  , "a", 'fontsize',fsz )
text(ax(2), gca().XLim(1)+2/(24*60),  15.5, "b", 'fontsize',fsz )
text(ax(3), gca().XLim(1)+2/(24*60), -65  , "c", 'fontsize',fsz )
text(ax(4), gca().XLim(1)+2/(24*60),  11  , "d", 'fontsize',fsz )
text(ax(1), tmin-15/(24*60),  22.7  , "front", 'fontsize',fsz )
text(ax(1), tmin+6/(24*60) ,  22.7  , "wake" , 'fontsize',fsz )

ax(2).YAxisLocation = 'right';
ax(4).YAxisLocation = 'right';
ax(5).Visible=false;

fmt = ["png", "epsc", "pdf", "svg"];
% for f = fmt
%     saveas(gcf, "cp2_fig3", f)
% end
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
