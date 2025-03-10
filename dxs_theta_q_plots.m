% Scatter plots %%
figure;
subplot(131)
scatter(d18O,dD,10,iso_time-iso_time(1),'filled')
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
colorbar
colormap(parula(16))
axis square

figure;
subplot(132)
scatter(d18O,dD,10,q_iso,'filled')
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
colorbar
colormap(jet)
axis square

% figure;
subplot(133)
scatter(d18O,dD,10,DXS,'filled')
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
colorbar
colormap(jet)
axis square

figure;
plot(q_ob,iso_de,'.b')
ylabel(['d-excess [',char(8240),']'])
xlabel('specific humidity [g/kg]')

%% Code needed to create line plots %%
flag_dxs = 9999*ones(size(DXS));
% Separating the data based on its DXS
for k = 1:length(iso_time)
    if DXS(k)<=10
        flag_dxs(k) = 1;
    elseif DXS(k)>10 && DXS(k)<=12
        flag_dxs(k) = 2;
    elseif DXS(k)>12 && DXS(k)<=14
        flag_dxs(k) = 3;
    elseif DXS(k)>14 && DXS(k)<=16
        flag_dxs(k) = 4;
    elseif DXS(k)>16
        flag_dxs(k) = 5;    
    else
        flag_dxs(k) = NaN;
    end
end
figure;
plot(iso_time,flag_dxs)

%% Code needed to create line plots %%
flag_t =  9999*ones(size(iso_time));
load 'cursor_info.mat'
c = 0;
for k = [1,3:5,7:10] % 1:10 % 
    c = c+1;
    index(c) = cursor_info(k).DataIndex;
end
time_cursor = iso_time(index);
% Separating the data based on its DXS and continous days
flag_t(1:index(end)) = 1;
flag_t(index(end)  :index(end-1)) = 2;
flag_t(index(end-1):index(end-2)) = 3;
flag_t(index(end-2):index(end-3)) = 4;
flag_t(index(end-3):index(end-4)) = 5;    
flag_t(index(end-4):index(end-5)) = 6;
flag_t(index(end-5):index(end-6)) = 7;
flag_t(index(end-6):index(end-7)) = 8;
flag_t(index(end-7):end) = 9;
% flag_t(index(end-7):index(end-8)) = 9;
% flag_t(index(end-8):index(end-9)) = 10;
% flag_t(index(end-9):end) = 11;

figure;
plot(iso_time,flag_t)

%% Generating line plots %%
figure;
co = colormap(parula(length(time_cursor)+1));
% subplot(131);
hold on;
for k = 1:length(time_cursor)+1
%     figure;
    plot(d18O(flag_t==k),dD(flag_t==k),'Color',co(k,:))
end
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
axis square

%% Only plotting cold pools %%
load 'cold_pool_flag_1min.mat'
load 'recovery&peak_flags_1min_full_timeseries.mat'
ind0 = find(t1min==iso_time(1));
pos  = ind0:1:ind0+length(iso_time)-1;
flag_cp = cold_pool_flag_1min(pos)';

figure; hold on;
scatter(d18O(flag_cp==1),dD(flag_cp==1),10,DXS(flag_cp==1),'filled')
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
axis square
h = colorbar;
clim([datenum('01/26/2020','mm/dd/yyyy') datenum('02/11/2020','mm/dd/yyyy')])
colormap(jet)
datetick(h,'y','mmm/dd','keeplimits','keepticks')

% Excluding cold pool recovery %
flag_rec = recovery_flag_1min(pos)';

figure; hold on;
scatter(d18O(flag_rec==0 & flag_cp==1),dD(flag_rec==0 & flag_cp==1),10,iso_time(flag_rec==0 & flag_cp==1),'filled')
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
axis square
h = colorbar;
clim([datenum('01/26/2020','mm/dd/yyyy') datenum('02/11/2020','mm/dd/yyyy')])
colormap(jet)
datetick(h,'y','mmm/dd','keeplimits','keepticks')

% Only plotting peak cold pool times %
flag_peak = peak_flag_1min(pos)';

figure;
scatter(d18O(flag_peak==1),dD(flag_peak==1),10,DXS(flag_peak==1),'filled')
ylabel(['\deltaD [',char(8240),']'])
xlabel(['\delta^1^8O [',char(8240),']'])
axis square
h = colorbar;
clim([datenum('01/26/2020','mm/dd/yyyy') datenum('02/11/2020','mm/dd/yyyy')])
colormap(jet)
datetick(h,'y','mmm/dd','keeplimits','keepticks')

%% DXS dD Ta timeseries %%
figure;
plot(time,dD-8*d18O,'-b')
hold on;
plot(time,Ta_full((2575:4715))-10,'-k')
plot(time,dD+80,'-m')
xlim([datenum('Jan/26/2020') datenum('Feb/11/2020 1PM')])
datetick('x','mmm/dd','keepticks','keeplimits')
ylabel('\delta_e_x_c_e_s_s [permil]')
ylim([0 20])
set(findall(gcf,'Type','axes'),'LineWidth',2)
set(findall(gcf,'-property','Fontsize'),'FontSize',20)

% Times of (12) cold pool soundings in both ATOMIC RHB legs:
% Identified by Simon de Szoeke

cpt = ['20200111_1051'
   '20200112_1055'
   '20200118_1848'
   '20200119_1049'
   '20200121_0249'
   '20200121_0641'
   '20200129_1437'
   '20200203_2240'
   '20200207_1509'
   '20200208_1439'
   '20200211_0243'
   '20200212_1041'];
cold_pools_times = datenum(cpt,'yyyymmdd_HHMM');
yy = [0 20];
xx = meshgrid(cold_pools_times,yy);
hold on;
plot(xx(:,6:11),yy,'-r')
legend('DXS','Ta (Ta - 10\circC)','dD (dD + 80 permil)')

%% DXS dD Ta timeseries in separate panels %%
figure;
subplot(414);
plot(time,dD-8*d18O,'-b')
xlim([datenum('Jan/29/2020') datenum('Feb/09/2020 7AM')])
datetick('x','mmm/dd','keepticks','keeplimits')
ylabel('\delta_e_x_c_e_s_s')
ylim([5 15])
% hold on;
% plot(xx(:,6:11),yy,'-r')
% set(gca,'XMinorTick','on')
grid on

subplot(411);
plot(time,Ta_full((2575:4715)),'-k')
xlim([datenum('Jan/29/2020') datenum('Feb/09/2020 7AM')])
datetick('x','mmm/dd','keepticks','keeplimits')
ylabel('Ta [\circC]')
% hold on;
% plot(xx(:,6:11),yy+17,'-r')
ylim([23 29])
% set(gca,'XMinorTick','on')
grid on

subplot(412);
plot(time,dD,'-m')
xlim([datenum('Jan/29/2020') datenum('Feb/09/2020 7AM')])
datetick('x','mmm/dd','keepticks','keeplimits')
ylabel('\deltaD')
% hold on;
% plot(xx(:,6:11),yy-80,'-r')
ylim([-80 -65])
% set(gca,'XMinorTick','on')
grid on

subplot(413);
plot(time,d18O,'-r')
xlim([datenum('Jan/29/2020') datenum('Feb/09/2020 7AM')])
datetick('x','mmm/dd','keepticks','keeplimits')
ylabel('\delta^1^8O')
ylim([-11.5 -9])
% hold on;
% plot(xx(:,6:11),yy,'-r')
% set(gca,'XMinorTick','on')
grid on

set(findall(gcf,'Type','axes'),'LineWidth',2)
set(findall(gcf,'-property','Fontsize'),'FontSize',20)

%% Dd timeseries colored with DXS values
figure;
hold on;
scatter(t,dD,10,iso_de,'filled')
xlim([datenum('Jan/29/2020') datenum('Feb/09/2020 7AM')])
datetick('x','mmm/dd','keepticks','keeplimits')
ylabel('dD [permil]')
% ylim([5 15])
hold on;
plot(xx(:,6:11),yy-77,'-r')
set(findall(gcf,'Type','axes'),'LineWidth',2)
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
box on
c = colorbar;
set(gca,'XMinorTick','on','YMinorTick','off')
grid on

%% DXS vs dD
figure;
plot(iso_de,dD,'ok','MarkerFaceColor','k','MarkerSize',2);
set(findall(gcf,'Type','axes'),'LineWidth',2)
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
xlabel('\delta_e_x_c_e_s_s [permil]')
ylabel('dD [permil]')
box on
grid on
axis square

%% Does DXS has a certain trend in the cold pool events?

% Choosing 3 cold pool events:
%    '20200129_1437'
%    '20200207_1509'
%    '20200208_1439'
% Using 5-hour windows for all three events

figure;
subplot(331);
plot(t,iso_de,'-b')
xlim([datenum('Jan/29/2020 12PM') datenum('Jan/29/2020 5PM')])
datetick('x','mmm/dd','keeplimits')
ylabel('\delta_e_x_c_e_s_s [permil]')
ylim([5 15])
hold on;
plot(xx(:,6:11),yy,'-r')
grid on

subplot(334);
plot(time,Ta_full((2575:4715)),'-k')
xlim([datenum('Jan/29/2020 12PM') datenum('Jan/29/2020 5PM')])
datetick('x','mmm/dd','keeplimits')
ylabel('Ta [\circC]')
hold on;
plot(xx(:,6:11),yy+17,'-r')
ylim([23 29])
grid on

subplot(337);
plot(time,dD,'-m')
xlim([datenum('Jan/29/2020 12PM') datenum('Jan/29/2020 5PM')])
datetick('x','mmm/dd','keeplimits')
ylabel('dD [permil]')
xlabel('Jan/29')
hold on;
plot(xx(:,6:11),yy-80,'-r')
ylim([-80 -65])
grid on

subplot(332);
plot(t,iso_de,'-b')
xlim([datenum('Feb/07/2020 11AM') datenum('Feb/07/2020 4PM')])
datetick('x','mmm/dd','keeplimits')
ylim([5 15])
hold on;
plot(xx(:,6:11),yy,'-r')
grid on

subplot(335);
plot(time,Ta_full((2575:4715)),'-k')
xlim([datenum('Feb/07/2020 11AM') datenum('Feb/07/2020 4PM')])
datetick('x','mmm/dd','keeplimits')
hold on;
plot(xx(:,6:11),yy+17,'-r')
ylim([23 29])
grid on

subplot(338);
plot(time,dD,'-m')
xlim([datenum('Feb/07/2020 11AM') datenum('Feb/07/2020 4PM')])
datetick('x','mmm/dd','keeplimits')
xlabel('Feb/07')
hold on;
plot(xx(:,6:11),yy-80,'-r')
ylim([-80 -65])
grid on

subplot(333);
plot(t,iso_de,'-b')
xlim([datenum('Feb/08/2020 11AM') datenum('Feb/08/2020 4PM')])
datetick('x','mmm/dd','keeplimits')
ylim([5 15])
hold on;
plot(xx(:,6:11),yy,'-r')
grid on

subplot(336);
plot(time,Ta_full((2575:4715)),'-k')
xlim([datenum('Feb/08/2020 11AM') datenum('Feb/08/2020 4PM')])
datetick('x','mmm/dd','keeplimits')
hold on;
plot(xx(:,6:11),yy+17,'-r')
ylim([23 29])
grid on

subplot(339);
plot(time,dD,'-m')
xlim([datenum('Feb/08/2020 11AM') datenum('Feb/08/2020 4PM')])
datetick('x','mmm/dd','keeplimits')
xlabel('Feb/08')
hold on;
plot(xx(:,6:11),yy-80,'-r')
ylim([-80 -65])
grid on

set(findall(gcf,'Type','axes'),'LineWidth',2)
set(findall(gcf,'-property','Fontsize'),'FontSize',20)

%% 4 panels - cold pool plots
open iso&Ta_1min_data_timeseries.fig
subplot(221);
yyaxis left
plot(t_avg_1min,d18O_avg_1min,'.b')
xlim([datenum('Jan/29/2020 11AM') datenum('Jan/29/2020 6PM')])
datetick('x','HHPM','keeplimits')
ylabel('\delta^1^8O [permil]')
ylim([-11.5 -9])
yyaxis right
ylim([22 27])
title('Jan/29/2020')

subplot(222);
yyaxis left
plot(t_avg_1min,d18O_avg_1min,'.b')
xlim([datenum('Feb/07/2020 8AM') datenum('Feb/07/2020 3PM')])
datetick('x','HHPM','keeplimits')
ylabel('\delta^1^8O [permil]')
ylim([-11.5 -9])
yyaxis right
ylim([22 27])
title('Feb/07/2020')

subplot(223);
yyaxis left
plot(t_avg_1min,d18O_avg_1min,'.b')
xlim([datenum('Feb/08/2020 12PM') datenum('Feb/08/2020 7PM')])
datetick('x','HHPM','keeplimits')
ylabel('\delta^1^8O [permil]')
ylim([-11.5 -9])
yyaxis right
ylim([22 27])
title('Feb/08/2020')

subplot(224);
yyaxis left
plot(t_avg_1min,d18O_avg_1min,'.b')
xlim([datenum('Feb/10/2020 11AM') datenum('Feb/10/2020 6PM')])
datetick('x','HHPM','keeplimits')
ylabel('\delta^1^8O [permil]')
ylim([-11.5 -9])
yyaxis right
ylim([22 27])
title('Feb/10/2020')
