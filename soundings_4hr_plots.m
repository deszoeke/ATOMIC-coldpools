% load full_214_soundings_Level2_h&p_same_size.mat
% load 2nd_leg_sounding_data_10min_linear_interp.mat
%% Data directly from the sounding netcdf file %%
% Ron Brown data %
filename = ('EUREC4A_RonBrown_Vaisala-RS_L2_v3.0.0.nc');
ascent_flag = ncread(filename,'ascent_flag'); % ascending = 1
h = ncread(filename,'alt'); % 'geopotential height' in meters
h = double(h);
q1 = ncread(filename,'q'); % 'specific humidity' in 'kg/kg'
q1 = q1(:,ascent_flag==1).*10^3; % in g/kg
rh1 = ncread(filename,'rh'); % 'relative_humidity' in decimal point
rh1 = rh1(:,ascent_flag==1).*100; % in percent
t1 = ncread(filename,'launch_time'); % 'time at which the sounding started' in 'seconds since 2020-01-01'
t1 = t1(ascent_flag==1)/3600/24 + datenum('20200101','yyyymmdd'); % in MatLab time
u1 = ncread(filename,'u');% eastward_wind 'u-component of the wind' in m/s
u1 = u1(:,ascent_flag==1);
v1 = ncread(filename,'v');% northward_wind 'v-component of the wind' in m/s
v1 = v1(:,ascent_flag==1);
th1 = ncread(filename,'theta'); % 'potential temperature' air pot temp in K
th1 = th1(:,ascent_flag==1);
ta1 = ncread(filename,'ta'); % 'dry bulb temperature' air temp in K
ta1 = ta1(:,ascent_flag==1);
% BCO data %
filename = ('EUREC4A_BCO_Vaisala-RS_L2_v3.0.0.nc');
ascent_flag2 = ncread(filename,'ascent_flag'); % ascending = 1
q2 = ncread(filename,'q'); % 'specific humidity' in 'kg/kg'
q2 = q2(:,ascent_flag2==1).*10^3; % in g/kg
rh2 = ncread(filename,'rh'); % 'relative_humidity' in decimal point
rh2 = rh2(:,ascent_flag2==1).*100; % in percent
t2 = ncread(filename,'launch_time'); % 'time at which the sounding started' in 'seconds since 2020-01-01'
t2 = t2(ascent_flag2==1)/3600/24 + datenum('20200101','yyyymmdd'); % in MatLab time
u2 = ncread(filename,'u');% eastward_wind 'u-component of the wind' in m/s
u2 = u2(:,ascent_flag2==1);
v2 = ncread(filename,'v');% northward_wind 'v-component of the wind' in m/s
v2 = v2(:,ascent_flag2==1);
th2 = ncread(filename,'theta'); % 'potential temperature' air pot temp in K
th2 = th2(:,ascent_flag2==1);
ta2 = ncread(filename,'ta'); % 'dry bulb temperature' air temp in K
ta2 = ta2(:,ascent_flag2==1);
%% Stitching BCO and RHB soundings together %%
d = t1(2:end)-t1(1:end-1); % difference in days
h8 = 8/24; % 8 hours in day units
ind = find(d>h8);
ind_inport1 = find(t2>=t1(ind(1)) & t2<=t1(ind(1)+1));
ind_inport2 = find(t2>=t1(ind(2)) & t2<=t1(ind(2)+1));

t = [t1(1:ind(1));t2(ind_inport1(2:end));t1(ind(1)+1:ind(2));t2(ind_inport2);t1(ind(2)+1:end)];
ta = [ta1(:,1:ind(1)),ta2(:,ind_inport1(2:end)),ta1(:,ind(1)+1:ind(2)),ta2(:,ind_inport2),ta1(:,ind(2)+1:end)];
th = [th1(:,1:ind(1)),th2(:,ind_inport1(2:end)),th1(:,ind(1)+1:ind(2)),th2(:,ind_inport2),th1(:,ind(2)+1:end)];
rh = [rh1(:,1:ind(1)),rh2(:,ind_inport1(2:end)),rh1(:,ind(1)+1:ind(2)),rh2(:,ind_inport2),rh1(:,ind(2)+1:end)];
q = [q1(:,1:ind(1)),q2(:,ind_inport1(2:end)),q1(:,ind(1)+1:ind(2)),q2(:,ind_inport2),q1(:,ind(2)+1:end)];
u = [u1(:,1:ind(1)),u2(:,ind_inport1(2:end)),u1(:,ind(1)+1:ind(2)),u2(:,ind_inport2),u1(:,ind(2)+1:end)];
v = [v1(:,1:ind(1)),v2(:,ind_inport1(2:end)),v1(:,ind(1)+1:ind(2)),v2(:,ind_inport2),v1(:,ind(2)+1:end)];

%% Loading data to plot inversion height %%
for k = 1:size(q,2)
    h6 = h(q(:,k)<=6);
    if h6(1) <= 6000
        trade_inv(k) = h6(1); % trade inversion height
    else
        trade_inv(k) = NaN; % trade inversion height
    end
    clearvars h6
end
%% Loading cold pool times %%
load 'cold_pool_detection_workspace_17.mat' t_min
cp_flag([1,12,17]) = 1; % cold pools without iso data
flag = 0; % cold pools with iso data
% flag = 1; % cold pools without iso data
h_ent = 1.2*ones(size(cp_flag(cp_flag==flag)));
% h_dd  = 1.1*ones(size(cp_flag(cp_flag==flag)));

%% Height plots for soundings preceding FEB 10 cold pool  %%
target = 'Feb 10 2020 12:00'; % target date
ind = find(t1>=datenum(target));
ind = ind(1:2);
% h = h./10^3;
figure; 
subplot(121); hold on;
    plot(q1(:,ind(1)),h);
    plot(rh1(:,ind(1)),h);
    plot(ta1(:,ind(1)),h);
    plot(th1(:,ind(1)),h);
    plot(u1(:,ind(1)),h);
    plot(v1(:,ind(1)),h);
    ylabel('Height [km]')
    ylim([0 5])
subplot(122); hold on;
    plot(q1(:,ind(2)),h);
    plot(rh1(:,ind(2)),h);
    plot(ta1(:,ind(2)),h);
    plot(th1(:,ind(2)),h);
    plot(u1(:,ind(2)),h);
    plot(v1(:,ind(2)),h);
    ylabel('Height [km]')
    ylim([0 5])
%% "Theta" and "q" time-height plots using contourf %%
figure;
    cmap = b2rcolormap(25);
subplot(2,1,1)
    [~,p] = contourf(t,h/1000,q);
    set(p,'LineColor', 'none');
    hold on; plot(t,trade_inv/1000,'-r','LineWidth',1)
    % plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/29/2020','mm/dd/yyyy') datenum('02/10/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/29/2020','mm/dd/yyyy'):2:datenum('02/10/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('specific humidity (g kg^-^1)')
    colorbar;
    caxis([0 18])
    ylabel('height [km]')
    ylim([0 4])
subplot(2,1,2)
    [~,p] = contourf(t,h/1000,th);
    set(p,'LineColor', 'none');
    hold on; plot(t,trade_inv/1000,'-r','LineWidth',1)
    % plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/29/2020','mm/dd/yyyy') datenum('02/10/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/29/2020','mm/dd/yyyy'):2:datenum('02/10/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('potential temperature (K)')
    colorbar;
%     caxis([280 300])
    ylabel('height [km]')
    set(gca,'TickDir','out');
%     xlabel('2020 date')
    ylim([0 4])
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'Type','axes'),'LineWidth',2)

%% "Theta" and "q" time-height plots using pcolor & contours %%
qmin = round(min(min(q))); % in g/kg
qmax = max(max(q)); % in g/kg
thmin = round(min(min(th))); % in K
thmax = max(max(th)); % in K

step = 5; % contour step
figure;
    cmap = b2rcolormap(25);
subplot(2,1,1)
    p = pcolor(t,h/1000,q);
    set(p,'EdgeColor', 'none');
    hold on;
    contour(t,h/1000,q,qmin:step:qmax,'Color','k');
    plot(t,trade_inv/1000,'Color',[0.7 0 0],'LineWidth',1)
    % plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/08/2020','mm/dd/yyyy') datenum('02/12/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/08/2020','mm/dd/yyyy'):5:datenum('02/13/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('specific humidity (g kg^-^1)')
    colorbar;
    caxis([0 20])
    ylabel('height [km]')
    ylim([0 4])
    % Adding labels for atmo layers
    hold on;
    text(datenum('Jan 08 2020 8PM'),3.8,'free troposphere','FontSize',18,'FontWeight','bold','Color','w')
    text(datenum('Jan 08 2020 8PM'),1.8,'trade cumulus','FontSize',16,'FontWeight','bold','Color','k')
    text(datenum('Jan 08 2020 8PM'),0.3,'subcloud mixed layer','FontSize',18,'FontWeight','bold','Color','k')
    text(datenum('Feb 05 2020'),1.7,'6g/kg','FontSize',18,'FontWeight','bold','Color',[0.7 0 0])% ,'r')

subplot(2,1,2)
    p = pcolor(t,h/1000,th);
    set(p,'EdgeColor', 'none');
    hold on;
    contour(t,h/1000,th,thmin:step:thmax,'Color','k');
    plot(t,trade_inv/1000,'Color',[0.7 0 0],'LineWidth',1)
    % plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/08/2020','mm/dd/yyyy') datenum('02/12/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/08/2020','mm/dd/yyyy'):5:datenum('02/13/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('potential temperature (K)')
    colorbar;
    caxis([295 315])
    ylabel('height [km]')
    set(gca,'TickDir','out');
%     xlabel('2020 date')
    ylim([0 4])
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'Type','axes'),'LineWidth',2)

%% "ta" and "rh" time-height plots using pcolor & contours %%
rhmin = round(min(min(rh))); % in percent
rhmax = max(max(rh)); % in percent
Tmin = round(min(min(ta))); % in K
Tmax = max(max(ta)); % in K
step = 5; % contour step
figure;
    cmap = b2rcolormap(21);
% ta
subplot(2,1,1)
    p = pcolor(t,h/1000,ta);
    set(p,'EdgeColor', 'none');
    hold on;
    contour(t,h/1000,ta,Tmin:step:Tmax,'Color','k');
    title('temperature (K)')
%     plot(t,trade_inv/1000,'Color',[0.7 0 0],'LineWidth',1)
    xlim([datenum('01/08/2020','mm/dd/yyyy') datenum('02/12/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/08/2020','mm/dd/yyyy'):5:datenum('02/13/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    colorbar;
    caxis([280 300])
    ylabel('height [km]')
    ylim([0 4])

% rh
subplot(2,1,2)
    p = pcolor(t,h/1000,rh);
    set(p,'EdgeColor', 'none');
    hold on;
    contour(t,h/1000,rh,rhmin:step*4:rhmax,'Color','k');
    title('relative humidity (%)')
%     plot(t,trade_inv/1000,'Color',[0.7 0 0],'LineWidth',1)
    xlim([datenum('01/08/2020','mm/dd/yyyy') datenum('02/12/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/08/2020','mm/dd/yyyy'):5:datenum('02/13/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    colorbar;
    caxis([0 100])
    ylabel('height [km]')
    ylim([0 4])
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'Type','axes'),'LineWidth',2)

%% "u" and "v" time-height plots using pcolor & contours %%
umin = round(min(min(u))); % in m/s
umax = max(max(u)); % in m/s
vmin = round(min(min(v))); % in m/s
vmax = max(max(v)); % in m/s
step = 5; % contour step
figure;
    cmap = b2rcolormap(25);
% u
subplot(2,1,1)
    p = pcolor(t,h/1000,u);
    set(p,'EdgeColor', 'none');
    hold on;
    contour(t,h/1000,u,umin:step:15,'Color','k');
    contour(t,h/1000,u,15:step:umax,'Color','w');
    title('zonal wind (m s^-^1)')
    plot(t,trade_inv/1000,'Color',[0.7 0 0],'LineWidth',1)
    xlim([datenum('01/08/2020','mm/dd/yyyy') datenum('02/12/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/08/2020','mm/dd/yyyy'):5:datenum('02/13/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    colorbar;
    caxis([-15 15])
    ylabel('height [km]')
    ylim([0 16])

% v
subplot(2,1,2)
    p = pcolor(t,h/1000,v);
    set(p,'EdgeColor', 'none');
    hold on;
    contour(t,h/1000,v,vmin:step:15,'Color','k');
    contour(t,h/1000,v,15:step:vmax,'Color','w');
    title('meridional wind (m s^-^1)')
    plot(t,trade_inv/1000,'Color',[0.7 0 0],'LineWidth',1)
    xlim([datenum('01/08/2020','mm/dd/yyyy') datenum('02/12/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/08/2020','mm/dd/yyyy'):5:datenum('02/13/2020','mm/dd/yyyy'))
    datetick('x','mmm/dd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    colorbar;
    caxis([-15 15])
    ylabel('height [km]')
    ylim([0 16])
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'Type','axes'),'LineWidth',2)
