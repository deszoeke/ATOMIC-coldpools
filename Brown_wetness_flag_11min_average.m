%% Ron H. Brown isotope analizer timeseries
% Loading Meteor data for comparison [OPTIONAL]
filename = 'data/M161_Picarro_Level1.nc';
d18O_Meteor = ncread(filename,'delta_18_16'); %
dD_Meteor   = ncread(filename,'delta_D'); %
time_Meteor = datenum('01012020','mmddyyyy') + ncread(filename,'time')*(1/(3600*24)); % from 'seconds since 2020-01-01 00:00:00' to MatLab time
dxs_Meteor = dD_Meteor - 8*d18O_Meteor;

filename = 'data/EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
dD = ncread(filename,'dD'); %
d18O = ncread(filename,'d18O'); %
dD_o = ncread(filename,'dD'); % original/"raw" variable
d18O_o = ncread(filename,'d18O'); % original/"raw" variable
time = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time = datenum('01012020','mmddyyyy') + time*(1/(3600*24));
dxs = dD - 8*d18O;
dxs_o = dD_o - 8*d18O_o;

% ship_flag = ncread(filename,'ship_flag'); %
% dD(ship_flag==1)   = NaN; %-------| poor
% d18O(ship_flag==1) = NaN; %       | data
% dxs(ship_flag==1)  = NaN; %-------| removed

inlet_flag = ncread(filename,'inlet_flag'); %
dD(inlet_flag==1)   = NaN; %-------| poor
d18O(inlet_flag==1) = NaN; %       | data
dxs(inlet_flag==1)  = NaN; %-------| removed

dDs = movmean(dD,11,'omitnan');
d18Os = movmean(d18O,11,'omitnan');
dxss = dDs - 8*d18Os;

%% Rain data
filename = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
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

%% Iso timeseries plot
figure;
subplot(311)
    plot(time_Meteor,d18O_Meteor,'o','MarkerEdgeColor',[202/260 204/260 206/260],'MarkerFaceColor',[202/260 204/260 206/260],'MarkerSize',2)
    hold on;
%     plot(time,d18O_o,'or','MarkerFaceColor','r','MarkerSize',2)
    hold on;
    plot(time,d18O,'ok','MarkerFaceColor','k','MarkerSize',1)
    hold on;
    plot(time,d18Os,'ob','MarkerFaceColor','b','MarkerSize',1)
    datetick('x','mm/dd')
    hold on;
    plot(time,d18O_rr,'og','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',3)
%     hold on;
%     plot(ctime,d18Op,'or','MarkerFaceColor','r')
%     legend('R/V Meteor','Observations; quality_p_o_o_r','Rain-Free Observations; quality_g_o_o_d','11-min running average for quality_g_o_o_d observations','Rain Measured >= 0.1 mm/hr for quality_g_o_o_d observations')
    legend('R/V Meteor','Rain-Free Observations; inlet_g_o_o_d','11-min running average for inlet_g_o_o_d observations','Rain Measured >= 0.1 mm/hr for inlet_g_o_o_d observations')
%     title('\delta^1^8O')
    ylabel('\delta^1^8O_v [‰]')
    xlabel('Time (mm/dd)')
    xlim([time(1) time(end)])
%     ylim([-12 max(d18Op)]) % when including precip data
    ylim([-12 -7])
    grid on
subplot(312)
    plot(time_Meteor,dD_Meteor,'o','MarkerEdgeColor',[202/260 204/260 206/260],'MarkerFaceColor',[202/260 204/260 206/260],'MarkerSize',2)
    hold on;
%     plot(time,dD_o,'or','MarkerFaceColor','r','MarkerSize',2)
    datetick('x','mm/dd')
    hold on;
    plot(time,dD,'ok','MarkerFaceColor','k','MarkerSize',1)
    hold on;
    plot(time,dDs,'ob','MarkerFaceColor','b','MarkerSize',1)
    hold on;
    plot(time,dD_rr,'og','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',3)
%     hold on;
%     plot(ctime,dDp,'or','MarkerFaceColor','r')
%     title('\deltaD')
    ylabel('\deltaD_v [‰]')
    xlabel('Time (mm/dd)')
    xlim([time(1) time(end)])
%     ylim([-80 max(dDp)]) % when including precip data
    ylim([-80 -50])
    grid on
subplot(313)
    plot(time_Meteor,dxs_Meteor,'o','MarkerEdgeColor',[202/260 204/260 206/260],'MarkerFaceColor',[202/260 204/260 206/260],'MarkerSize',2)
    hold on;
%     plot(time,dxs_o,'or','MarkerFaceColor','r','MarkerSize',2)
    datetick('x','mm/dd')
    hold on;
    plot(time,dxs,'ok','MarkerFaceColor','k','MarkerSize',1)
    hold on;
    plot(time,dxss,'ob','MarkerFaceColor','b','MarkerSize',1)
    hold on;
    plot(time,dxs_rr,'og','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',3)
%     hold on;
%     plot(ctime,dxsp,'or','MarkerFaceColor','r')
%     title('D-Excess')
    ylabel('DXS_v [‰]')
    xlabel('Time (mm/dd)')
    xlim([time(1) time(end)])
%     ylim([min(dxsp) 22]) % when including precip data
    ylim([5 22])
    grid on
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% Rain rate timeseries plot
figure;
    bar(time_rr,rr_o,'r')
    datetick('x','mm/dd')
    hold on;
    bar(time,rr,'b')
%     hold on;
%     plot(time,rrs,'oc','MarkerFaceColor','c','MarkerSize',3)
%     title('rain rate')
    ylabel('rain rate (mm/hr)')
    xlabel('Time')
%     xlim([time(1) time(end)])
    grid on
% Zoom-in to Feb 10 cold pool time frame [OPTIONAL]
    xlim([datenum('02/10/2020 15:40','mm/dd/yyyy HH:MM') datenum('02/10/2020 17:00','mm/dd/yyyy HH:MM')])
    datetick('x','HH:MM')
    ylim([0 15])
% set(findall(gcf,'-property','FontSize'),'FontSize',34)
