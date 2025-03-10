% Script to plot timeseries figure with cold pools times shaded and min T
% times dotted
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
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]
rh = ncread(filename,'rhair'); % relative humidity; units: percentage
qair = ncread(filename,'qair'); % specific humidity; units: g/kg
rr = ncread(filename,'prate'); % rain rate; units: mm/hr
wspd = ncread(filename,'wspd'); % true wind speed; units: m/s
% u = wspd.*cosd(270-wdir); % eastward wind speed, m/s
% v = wspd.*sind(270-wdir); % northward wind speed, m/s

% Ta = movmean(Ta,11,'omitnan'); % Tair filtered => 11-min running average
% qair = movmean(qair,11,'omitnan'); % qair filtered => 11-min running average

rdir = zeros(size(time));
for k = 1:length(time)
    rdir(k) = rdir_o(time_rr==time(k));
end
ship_flag = zeros(size(rdir));
ship_flag(rdir>-135 & rdir>45) = 1; % 1 = bad wind dir

dD(ship_flag==1 | inlet_flag==1) = NaN;
d18O(ship_flag==1 | inlet_flag==1) = NaN;

% clearvars filename inlet_flag ship_flag rdir rdir_o time_rr
% dD   = movmean(dD,11,'omitnan');   % dD filtered => 11-min running average % in per mil
% d18O = movmean(d18O,11,'omitnan'); % d18O filtered => 11-min running average % in per mil
DXS = dD - 8*d18O;

%%
load('cold_pool_flag_1min.mat')
load('workspace_100_cp_detection_algorithm_1min.mat')

%% Timeseries plot for the WIW poster %%
figure;
subplot(411)
    % variables for shading
    yvalues = [22 22 28 28];
    hold on;
    for k = [65,67,69:78,82:84,86]
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time_rr,Ta,'-k','LineWidth',2)
    hold on;
    scatter(time_rr(t_min_ind([65,67,69:78,82:84,86])),Ta(t_min_ind([65,67,69:78,82:84,86])),'r','filled')
    ylim([22.5 28])
    xticks(datenum('Jan/26/2020'):1:datenum('Feb/11/2020'))
    set(gca,'xticklabels',[])
    xlim([datenum('Jan/26/2020 4PM') datenum('Feb/11/2020')])
    ylabel('T_a [\circC]')
    box on
subplot(412)
    % variables for shading
    yvalues = [11 11 18 18];
    hold on;
    for k = [65,67,69:78,82:84,86]
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time_rr,qair,'-k','LineWidth',2)
    hold on;
    scatter(time_rr(t_min_ind([65,67,69:78,82:84,86])),qair(t_min_ind([65,67,69:78,82:84,86])),'r','filled')
    ylim([11 17.5])
    xticks(datenum('Jan/26/2020'):1:datenum('Feb/11/2020'))
    set(gca,'xticklabels',[])
    xlim([datenum('Jan/26/2020 4PM') datenum('Feb/11/2020')])
    ylabel('T_a [\circC]')
    box on
subplot(413)
    % variables for shading
    yvalues = [-80 -80 -64 -64];
    hold on;
    for k = [65,67,69:78,82:84,86]
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time,dD,'-k','LineWidth',2)
    hold on;
    scatter(time_rr(t_min_ind([65,67,69:78,82:84,86])),dD(t_min_ind([65,67,69:78,82:84,86])-factor),'r','filled')
    xticks(datenum('Jan/26/2020'):1:datenum('Feb/11/2020'))
    set(gca,'xticklabels',[])
    xlim([datenum('Jan/26/2020 4PM') datenum('Feb/11/2020')])
    ylim([-80 -65])
    box on
    ylabel(['\deltaD [',char(8240),']'])
subplot(414)
    % variables for shading
    yvalues = [0 0 15 15];
    hold on;
    for k = [65,67,69:78,82:84,86]
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time_rr,rr,'-k','LineWidth',2)
    hold on;
    scatter(time_rr(t_min_ind([65,67,69:78,82:84,86])),rr(t_min_ind([65,67,69:78,82:84,86])),'r','filled')
    ylim([0 15])
    xticks(datenum('Jan/26/2020'):1:datenum('Feb/11/2020'))
    xticklabels(26:42)
    xlim([datenum('Jan/26/2020 4PM') datenum('Feb/11/2020')])
    ylabel('RR [mm/hr]')
    xlabel('DOY 2020')
    box on
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
set(findall(gcf,'-property','TickDir'),'TickDir','both')
