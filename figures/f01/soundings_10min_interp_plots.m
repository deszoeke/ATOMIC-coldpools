load data/full_214_soundings_Level2_h_p_same_size.mat
% load 2nd_leg_sounding_data_10min_linear_interp.mat
%% Loading data to plot inversion height %%
for k = 12:size(q,2)
    h6 = double(h(q(:,k)*1000<=6));
    if h6(1) <= 6000
        trade_inv(k) = h6(1); % trade inversion height
    else
        trade_inv(k) = NaN; % trade inversion height
    end
    clearvars h6
end
%% Loading cold pool times %%
load data/cold_pool_detection_workspace_17.mat t_min
cp_flag([1,12,17]) = 1; % cold pools without iso data
flag = 0; % cold pools with iso data
% flag = 1; % cold pools without iso data
h_ent = 1.2*ones(size(cp_flag(cp_flag==flag)));
% h_dd  = 1.1*ones(size(cp_flag(cp_flag==flag)));

%% ceilometer data %%
filename = 'data/EUREC4A_ATOMIC_RonBrown_Ceilometer_15s_20200109-20200212_v1.1.nc';
t_cb   = ncread(filename,'time')/3600/24 + datenum('20200101','yyyymmdd');
cb = ncread(filename,'first_cloud_base');

%% "Theta" and "q" time-height plots using pcolor & contours %%
c = -54:3:54;
qmin = round(min(min(q))*1000); % in g/kg
qmax = max(max(q))*1000; % in g/kg
Tmin = round(min(min(thw))); % in K
Tmax = max(max(thw)); % in K
step = 5; % contour step
toff = 2/24; % to center pcolor
figure;
clf;
    cmap = b2rcolormap(25);
subplot(4,1,1)
    p = pcolor(t-toff,h/1000,q*1000);
    set(p,'EdgeColor', 'none');
    % set(p, 'Rasterization', 'on');
    hold on;
    contour(t,h/1000,q*1000,qmin:step:qmax,'Color','k');
    plot(t(12:end),trade_inv(12:end)/1000,'-r','LineWidth',1)
    plot(t_cb,cb/1e3,'.', 'color',0.7+[0 0 0], 'markersize',2)
    plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/29/2020','mm/dd/yyyy') datenum('02/10/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/29/2020','mm/dd/yyyy'):2:datenum('02/10/2020','mm/dd/yyyy'))
    datetick('x','mmmdd','keeplimits','keepticks')
    xticklabels('')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('specific humidity (g kg^-^1)', 'fontweight','normal')
    colorbar;
    caxis([0 20])
    % ylabel('height [km]')
    ylim([0 4])
    text(datenum('Feb 05 2020'),1.7,'6g/kg','FontSize',18,'Color',[0.7 0 0])% ,'r') % humicline
    text(datenum(2020,1,28,0,0,0), 4.7, 'a', 'clipping','off')

subplot(4,1,2)
    p = pcolor(t-toff,h/1000,thw);
    set(p,'EdgeColor', 'none');
    % set(p, 'Rasterization', 'on');
    hold on;
    contour(t,h/1000,thw,Tmin:step:Tmax,'Color','k');
    plot(t(12:end),trade_inv(12:end)/1000,'-r','LineWidth',1)
    plot(t_cb,cb/1e3,'.', 'color',0.7+[0 0 0], 'markersize',2) % ceilo cloud base
    plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/29/2020','mm/dd/yyyy') datenum('02/10/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/29/2020','mm/dd/yyyy'):2:datenum('02/10/2020','mm/dd/yyyy'))
    datetick('x','mmmdd','keeplimits','keepticks')
    xticklabels('')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('wet-bulb potential temperature (K)', 'fontweight','normal')
    colorbar;
    caxis([280 300])
    ylabel('height (km)')
    set(gca,'TickDir','out');
    ylim([0 4])
    text(datenum('Jan-29-2020 03:00'),3.8,'free troposphere','FontSize',22, 'Color','k')
    text(datenum('Jan-29-2020 03:00'),1.6,'trade cumulus','FontSize',22, 'Color','k')
    text(datenum('Jan-29-2020 03:00'),0.3,'subcloud mixed layer','FontSize',22, 'Color','k')
    text(datenum(2020,1,28,0,0,0), 4.7, 'b', 'clipping','off')

subplot(4,1,3)
    p = pcolor(t-toff,h/1000,th);
    set(p,'EdgeColor', 'none');
    % set(p, 'Rasterization', 'on');
    hold on;
    contour(t,h/1000,th,300:5:320,'Color','k');
    plot(t(12:end),trade_inv(12:end)/1000,'-r','LineWidth',1)
    plot(t_cb,cb/1e3,'.', 'color',0.7+[0 0 0], 'markersize',2) % ceilo cloud base
    plot(t_min(cp_flag==flag),h_ent,'k*')
    % plot(t_min(cp_flag==flag),h_dd,'k*')
    xlim([datenum('01/29/2020','mm/dd/yyyy') datenum('02/10/2020 06PM','mm/dd/yyyy HHPM')])
    xticks(datenum('01/29/2020','mm/dd/yyyy'):2:datenum('02/10/2020','mm/dd/yyyy'))
    datetick('x','mmmdd','keeplimits','keepticks')
    set(gca,'TickDir','out','XMinorTick','on');
    colormap(cmap)
    title('potential temperature (K)', 'fontweight','normal')
    colorbar;
    caxis([297 320])
    % ylabel('height [km]')
    set(gca,'TickDir','out');
    xlabel('2020')
    ylim([0 4])
    text(datenum(2020,1,28,0,0,0), 4.7, 'c', 'clipping','off')

% finish plot
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
set(findall(gcf,'Type','axes'),'LineWidth',2)
set(gcf,'renderer','painters')
orient tall;

fmt = ["png", "epsc", "pdf", "svg"];
for i = 1:length(fmt)
    f = strcat('-d', fmt(i));
    print(gcf, "timeheight_q_thw_th", f)
end
