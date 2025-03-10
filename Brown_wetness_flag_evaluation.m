%% Iso data
filename = 'M161_Picarro_Level1.nc';
d18O_Meteor = ncread(filename,'delta_18_16'); %
dD_Meteor   = ncread(filename,'delta_D'); %
time_Meteor = datenum('01012020','mmddyyyy') + ncread(filename,'time')*(1/(3600*24)); % from 'seconds since 2020-01-01 00:00:00' to MatLab time
dxs_Meteor = dD_Meteor - 8*d18O_Meteor;

filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
time_rr = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time_rr = datenum('01012020','mmddyyyy') + time_rr/3600/24;
rr = ncread(filename,'prate'); % rain rate in mm/hr
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]
rdir_o = ncread(filename,'rdir'); % original/"raw" variable

filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
dD = ncread(filename,'dD'); %
d18O = ncread(filename,'d18O'); %
time = ncread(filename,'time'); % in 'seconds since 2020-01-01 00:00:00'
time = datenum('01012020','mmddyyyy') + time*(1/(3600*24));

rdir = zeros(size(time));
for k = 1:length(time)
    rdir(k) = rdir_o(time_rr==time(k));
end
ship_flag = zeros(size(rdir));
ship_flag(rdir>-135 & rdir>45) = 1; % 1 = bad wind dir
% ship_flag = ncread(filename,'ship_flag'); %
inlet_flag = ncread(filename,'inlet_flag'); %

dDs = dD(ship_flag==0 & inlet_flag==0);
d18Os = d18O(ship_flag==0 & inlet_flag==0);
dxss = dDs - 8*d18Os;
times = time(ship_flag==0 & inlet_flag==0);

% dDs = dD(ship_flag==0 | ship_flag==2);
% d18Os = d18O(ship_flag==0 | ship_flag==2);
% dxss = dDs - 8*d18Os;
% times = time(ship_flag==0 | ship_flag==2);

% dDs = dD(inlet_flag==0);
% d18Os = d18O(inlet_flag==0);
% dxss = dDs - 8*d18Os;
% times = time(inlet_flag==0);

dD = dD(ship_flag==1 | inlet_flag==1); %-------|
d18O = d18O(ship_flag==1 | inlet_flag==1); %   | poor
dxs = dD - 8*d18O; %                           | data
time = time(ship_flag==1 | inlet_flag==1); %---|

% dD = dD(inlet_flag==1); %-------|
% d18O = d18O(inlet_flag==1); %   | poor
% dxs = dD - 8*d18O; %            | data
% time = time(inlet_flag==1); %---|

filename = 'EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
dDp = ncread(filename,'dD'); %
d18Op = ncread(filename,'d18O'); %
ctime = ncread(filename,'collection_time'); %
ctime = datenum('01012020','mmddyyyy') + ctime*(1/(3600*24));
delay_flag = ncread(filename,'delay_flag'); %
dxsp = dDp - 8*d18Op;

time_ind = zeros(size(time));
for k = 1:length(time)
    time_ind(k) = find(time_rr==time(k));
end

time_inds = zeros(size(times));
for k = 1:length(times)
    time_inds(k) = find(time_rr==times(k));
end

dD_rr = dD; %-------|
d18O_rr = d18O; %   | poor
dxs_rr = dxs; %-----| data

rr_2 = rr(time_ind);
dD_rr(isnan(rr_2)) = NaN;
d18O_rr(isnan(rr_2)) = NaN;
dxs_rr(isnan(rr_2)) = NaN;

dD_rr(rr_2<=0.1) = NaN;
d18O_rr(rr_2<=0.1) = NaN;
dxs_rr(rr_2<=0.1) = NaN;

dD_rrs = dDs;
d18O_rrs = d18Os;
dxs_rrs = dxss;

rr_s = rr(time_inds);
dD_rrs(isnan(rr_s)) = NaN;
d18O_rrs(isnan(rr_s)) = NaN;
dxs_rrs(isnan(rr_s)) = NaN;

dD_rrs(rr_s<=0.1) = NaN;
d18O_rrs(rr_s<=0.1) = NaN;
dxs_rrs(rr_s<=0.1) = NaN;

%% INCLUDING COLD POOLS %%
[t_max,t_min,t_max_ind,t_min_ind,t_end,t_end_ind,end_flag,cp_ind,delta_T,T_max,T_min,Taf] = cold_pool_detection_algorithm(time_rr,Ta);
%%
figure;
subplot(311)
    % variables for shading
    yvalues = [-12 -12 -7 -7];
    hold on;
    for k = 1:length(cp_ind)
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time_Meteor,d18O_Meteor,'o','MarkerEdgeColor',[202/260 204/260 206/260],'MarkerFaceColor',[202/260 204/260 206/260],'MarkerSize',2)
    hold on;
    plot(times,d18Os,'ob','MarkerFaceColor','b','MarkerSize',2)
    hold on;
    plot(time,d18O,'ok','MarkerFaceColor','b','MarkerSize',3)
    hold on;
    plot(time_rr(time_inds),d18O_rrs,'ok','MarkerFaceColor','g')
    hold on;
    plot(time_rr(time_ind),d18O_rr,'ok','MarkerFaceColor','r')
%     hold on;
%     plot(ctime,d18Op,'or','MarkerFaceColor','r')
    legend('R/V Meteor','Rain-Free Observation; quality_g_o_o_d','Rain-Free Observation; quality_p_o_o_r','Rain Measured >= 0.1 mm/hr for quality_g_o_o_d observations','Rain Measured >= 0.1 mm/hr for quality_p_o_o_r observations','Precip sample')
%     legend('R/V Meteor','Rain-Free Observation; inlet_g_o_o_d','Rain-Free Observation; inlet_r_e_v','Rain Measured >= 0.1 mm/hr for inlet_g_o_o_d observations','Rain Measured >= 0.1 mm/hr for inlet_r_e_v observations','Precip sample')
    title('\delta^1^8O')
    ylabel('\delta^1^8O (per mil)')
    xlabel('Time (mm/dd)')
    xlim([times(1) times(end)])
    datetick('x','mm/dd','keeplimits')    
%     ylim([-12 max(d18Op)])
    ylim([-12 -7])
    grid on
    box on
subplot(312)
    % variables for shading
    yvalues = [-80 -80 -50 -50];
    hold on;
    for k = 1:length(cp_ind)
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time_Meteor,dD_Meteor,'o','MarkerEdgeColor',[202/260 204/260 206/260],'MarkerFaceColor',[202/260 204/260 206/260],'MarkerSize',2)
    hold on;
    plot(times,dDs,'om','MarkerFaceColor','m','MarkerSize',2)
    hold on;
    plot(time,dD,'ok','MarkerFaceColor','m','MarkerSize',3)
    hold on;
    plot(time_rr(time_inds),dD_rrs,'ok','MarkerFaceColor','g')
    hold on;
    plot(time_rr(time_ind),dD_rr,'ok','MarkerFaceColor','r')
%     hold on;
%     plot(ctime,dDp,'or','MarkerFaceColor','r')
    title('\deltaD')
    ylabel('\deltaD (per mil)')
    xlabel('Time (mm/dd)')
    xlim([times(1) times(end)])
    datetick('x','mm/dd','keeplimits')    
%     ylim([-80 max(dDp)])
    ylim([-80 -50])
    grid on
    box on
subplot(313)
    % variables for shading
    yvalues = [5 5 22 22];
    hold on;
    for k = 1:length(cp_ind)
    xvalues = [time_rr(t_max_ind(k)) time_rr(t_end_ind(k)) time_rr(t_end_ind(k)) time_rr(t_max_ind(k))];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars xvalues j
    end
    plot(time_Meteor,dxs_Meteor,'o','MarkerEdgeColor',[202/260 204/260 206/260],'MarkerFaceColor',[202/260 204/260 206/260],'MarkerSize',2)
    hold on;
    plot(times,dxss,'oc','MarkerFaceColor','c','MarkerSize',2)
    hold on;
    plot(time,dxs,'ok','MarkerFaceColor','c','MarkerSize',3)
    hold on;
    plot(time_rr(time_inds),dxs_rrs,'ok','MarkerFaceColor','g')
    hold on;
    plot(time_rr(time_ind),dxs_rr,'ok','MarkerFaceColor','r')
%     hold on;
%     plot(ctime,dxsp,'or','MarkerFaceColor','r')
    title('D-Excess')
    ylabel('\DeltaD (per mil)')
    xlabel('Time (mm/dd)')
    xlim([times(1) times(end)])
    datetick('x','mm/dd','keeplimits')    
%     ylim([min(dxsp) 22])
    ylim([5 22])
    grid on
    box on
set(findall(gcf,'-property','FontSize'),'FontSize',34)

%% Scatter plots
figure;
subplot(221)
    plot(rr_s,d18O_rrs,'ob','MarkerFaceColor','b')
    hold on;
    plot(rr_2,d18O_rr,'ok','MarkerFaceColor','r')    
    title('\delta^1^8O')
    ylabel('\delta^1^8O (per mil)')
    xlabel('rain rate (mm/hr)')
    ylim([-12 -7])
    set(gca,'XScale','log')
    legend('Iso Observation; quality_g_o_o_d','Iso Observation; quality_p_o_o_r');
%     legend('Iso Observation; inlet_g_o_o_d','Iso Observation; inlet_r_e_v');
    grid on
subplot(222)
    plot(rr_s,dD_rrs,'om','MarkerFaceColor','m')
    hold on;
    plot(rr_2,dD_rr,'ok','MarkerFaceColor','r')    
    title('\deltaD')
    ylabel('\deltaD (per mil)')
    xlabel('rain rate (mm/hr)')
    ylim([-80 -50])
    set(gca,'XScale','log')
    grid on
subplot(223)
    plot(rr_s,dxs_rrs,'oc','MarkerFaceColor','c')
    hold on;
    plot(rr_2,dxs_rr,'ok','MarkerFaceColor','r')    
    title('D-Excess')
    ylabel('\DeltaD (per mil)')
    xlabel('rain rate (mm/hr)')
    ylim([5 22])
    set(gca,'XScale','log')
    grid on
set(findall(gcf,'-property','FontSize'),'FontSize',34)
