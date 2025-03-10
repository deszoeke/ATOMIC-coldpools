%% pdf from P3 data
addpath 'D:\BackUp_Personal_Laptop_March2022\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files\Dean''s P3 files with flags'
addpath 'C:\Users\quinones\Documents\Data\thermo'
% adri_data = ncinfo('P3_ATOMIC_RF01_Isotopes1Hz_V0.1.nc');
% Variables: Time dD d18O
dean_flags = ncinfo('P3_with_vertlocflags_20200210.nc');
% Variables: Time dD d18O qs thetaEs thetav RH Ta ws wd alt
%%
% date = ['0131'; '0203'; '0204'; '0205'; '0209'; '0210']; % dates that coincide with cold pools
date = ['0117';'0119';'0123';'0124';'0131';'0203';'0204';'0205';'0209';'0210';'0211']; % all flight dates
n = length(date);

qair_P3  = zeros(10000,n);
qs_P3    = zeros(10000,n);
dD_P3    = zeros(10000,n);
Ta_P3    = zeros(10000,n);
d18O_P3  = zeros(10000,n);
alt_P3   = zeros(10000,n);
tsecs_P3 = zeros(10000,n);
theta_P3 = zeros(10000,n);
thw_P3   = zeros(10000,n); % theta_w
p_P3 = zeros(10000,n); 

for ii = 1:n
disp(ii)
% temp and qs variables
filename = ['P3_with_vertlocflags_2020',date(ii,:),'.nc'];

% mmr variable
% filename = ['EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-1Hz_2020',date(ii,:),'_v1.1.nc'];

% press and isotope ratios
% filename = ['EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-Isotope-Ratios-1Hz_2020',date(ii,:),'_v1.1.nc'];

tsecs = ncread(filename,'secs_since_2020Jan01'); % time of the P3 data
% tsecs = ncread(filename,'time'); % time of the P3 data
tsecs = tsecs/3600/24 + datenum('20200101','yyyymmdd');
% t1min =
% t10min =
Ta = ncread(filename,'Ta'); % air temperature [in degrees K]
% rh = ncread(filename,'RH'); % relative humidity; units: percentage
qs = ncread(filename,'qs'); % saturated specific humidity; units: kg/kg
qair = ncread(filename,'mmr'); % humidity mixing ratio; ~specific humidity; units: g/kg
% ratio of the mass of water vapor to the mass of dry air
% water vapor mass (dry) mixing ratio
% wspd = ncread(filename,'ws'); % wind speed; units: m/s
dD = ncread(filename,'dD'); % 
d18O = ncread(filename,'d18O'); %
alt =  ncread(filename,'alt'); %
p = ncread(filename,'press'); % corrected static pressure: hPa

% th  = theta(p,Ta); % theta(p[Pa],Temp[K])
Rd = 287.04; % units???
Cp = 1005.7; % units???
th = (Ta).*(1e5./(p*100)).^(Rd/Cp); % (Potential Temp in degrees K)

thw = theta_w(p*100,Ta,qs); % theta_w(p[Pa],Temp[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
% thw = theta_w_eqm(p*100,th,qair/1000); % theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008

Ta   = Ta(~isnan(dD));
qair = qair(~isnan(dD));
qs   = qs(~isnan(dD));
d18O = d18O(~isnan(dD));
tsecs = tsecs(~isnan(dD));
alt  = alt(~isnan(dD));
th   = th(~isnan(dD));
thw  = thw(~isnan(dD));
p = p(~isnan(dD)); 
dD   = dD(~isnan(dD));

hlevel = 3000; % in meters
% hlevel = 400; % in meters

Ta = Ta(alt<=hlevel);
qair = qair(alt<=hlevel);
qs = qs(alt<=hlevel);
d18O = d18O(alt<=hlevel);
tsecs = tsecs(alt<=hlevel);
dD  = dD(alt<=hlevel);
th = th(alt<=hlevel);
thw = thw(alt<=hlevel);
p = p(alt<=hlevel);
alt = alt(alt<=hlevel);

% [N,XEDGES,YEDGES] = histcounts2(qair,dD,'Normalization','pdf');
% Probability density function estimate. The height 
%  of each bar is, (number of observations in bin)
%  / (total number of observations * area of bin).
%  The volume of each bar is the relative number of
%  observations, and the sum of the bar volumes 
%  is less than or equal to 1.

% figure;
% contour(XEDGES(2:end),YEDGES(2:end),N')
% colormap(jet(12))
% colorbar
% xlabel('specific humidity [g/kg]')
% ylabel('dD [per mil]')
% set(findall(gcf,'-property','Fontsize'),'FontSize',14)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
% title(['P3 flight on ',datestr(tsecs(1))])
% caxis([0.0 0.3])
% xlim([7 19]); ylim([-80 -60])

xx = length(qair);
qair_P3(1:xx,ii) = qair;
qs_P3(1:xx,ii) = qs;
dD_P3(1:xx,ii) = dD;
Ta_P3(1:xx,ii) = Ta;
d18O_P3(1:xx,ii) = d18O;
alt_P3(1:xx,ii) = alt;
tsecs_P3(1:xx,ii) = tsecs;
theta_P3(1:xx,ii) = th;
thw_P3(1:xx,ii) = thw;
p_P3(1:xx,ii) = p;

clearvars qair dD d18O alt th tsecs Ta thw qs p
end

qair_P3(qair_P3==0) = NaN;
qs_P3(qs_P3==0)     = NaN;
dD_P3(dD_P3==0)     = NaN;
d18O_P3(d18O_P3==0) = NaN;
Ta_P3(Ta_P3==0)     = NaN;
tsecs_P3(tsecs_P3==0) = NaN;
alt_P3(alt_P3<=0)   = NaN;
theta_P3(theta_P3==0) = NaN;
thw_P3(thw_P3==0)   = NaN;
p_P3(p_P3==0)       = NaN;

%% For including in weakest/strongest cold pools [dot] plots
% average quantities below 1000m from 6 relevant P3 flights
m_qair_P3 = nanmean(nanmean(qair_P3));
m_dD_P3 = nanmean(nanmean(dD_P3));
m_Ta_P3 = nanmean(nanmean(Ta_P3));
m_d18O_P3 = nanmean(nanmean(d18O_P3));

%% PDF plot of full flights below the hlevel

[N,XEDGES,YEDGES] = histcounts2(qair_P3,dD_P3,'Normalization','pdf');
figure;
contour(XEDGES(2:end),YEDGES(2:end),N')
colormap(jet(20))
colorbar
xlabel('specific humidity [g/kg]')
ylabel('dD [per mil]')
set(findall(gcf,'-property','Fontsize'),'FontSize',14)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
title(['Jan 31 to Feb 10 - P3 flights below ',num2str(hlevel),'m'])
caxis([0.0 0.1])
xlim([12 18]); ylim([-76 -65])

[N,XEDGES,YEDGES] = histcounts2(qair_P3,Ta_P3,'Normalization','pdf');
figure;
contour(XEDGES(2:end),YEDGES(2:end),N')
colormap(jet(20))
colorbar
xlabel('specific humidity [g/kg]')
ylabel('air temperature [\circC]')
set(findall(gcf,'-property','Fontsize'),'FontSize',14)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
title(['Jan 31 to Feb 10 - P3 flights below ',num2str(hlevel),'m'])
caxis([0.0 0.1])
xlim([12 18]); ylim([18 31])

%% Plots scripted on 10/10/2022 %%
% dD_P3 = dD_P3 - 8*d18O_P3;
figure;
sample_id = 1:length(Ta_P3);
h = hlevel*ones(size(Ta_P3));
subplot(211)
m = mean(dD_P3(:,5),'omitnan');
plot(sample_id,dD_P3(:,5))
hold on; plot(sample_id,mean(dD_P3(:,5),'omitnan')*ones(size(Ta_P3)))
ylabel(['\deltaD [',char(8240),']'])
xlabel('sample ID')
title(['1km height level flight 0209 mean = ' num2str(m) ''])

subplot(212)
m = mean(dD_P3(:,6),'omitnan');
plot(sample_id,dD_P3(:,6))
hold on; plot(sample_id,mean(dD_P3(:,6),'omitnan')*ones(size(Ta_P3)))
ylabel(['\deltaD [',char(8240),']'])
text();
xlabel('sample ID')
title(['1km height level flight 0210 mean = ' num2str(m) ''])

%% Time-height Plots scripted on 10/13/2022 %%
% generated to determine at which level depleted isotopes are present!!!

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],dD_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = '\deltaD [‰]';
caxis([-100 -60])
colormap(parula(8))
ylabel('height [m]')
xlabel('date in 2020')

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],d18O_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = '\delta^1^8O [‰]';
caxis([-16 -8])
colormap(parula(8))
ylabel('height [m]')
xlabel('date in 2020')

figure;
DXS_P3 = dD_P3 - 8*d18O_P3;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],DXS_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = 'DXS [‰]';
caxis([8 16])
colormap(parula(8))
ylabel('height [m]')
xlabel('date in 2020')

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],theta_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = '\theta [K]';
caxis([295 320])
% colormap(cmap)
ylabel('height [m]')
xlabel('date in 2020')

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],thw_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = '\theta_w [K]';
% caxis([-16 -8])
% colormap(parula(8))
ylabel('height [m]')
xlabel('date in 2020')

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],Ta_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = 'T_a [K]';
caxis([280 304])
h.YTick = 280:4:312;
% colormap(parula(6))
cmap = b2rcolormap(13);
ylabel('height [m]')
xlabel('date in 2020')

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],qair_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = 'mmr [g/kg]';
caxis([0 20])
% h.YTick = 280:4:312;
% colormap(parula(6))
cmap = b2rcolormap(11);
ylabel('height [m]')
xlabel('date in 2020')

figure;
for ii = 1:n
    hold on;
    scatter(tsecs_P3(:,ii),alt_P3(:,ii),[],qs_P3(:,ii),'filled','s')
end
datetick('x','mm/dd','keeplimits','keepticks')
h = colorbar;
h.Title.String = 'qs [g/kg]';
caxis([5 25])
% h.YTick = 280:4:312;
% colormap(parula(6))
cmap = b2rcolormap(11);
ylabel('height [m]')
xlabel('date in 2020')

%% Plotting P3 profiles - script edited on 10/17/2022 %%
for ii = 1:n
    figure;
    subplot(151);plot(theta_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\theta [K]'); ylabel('alt [m]'); xlim([295 315])
    subplot(152);plot(p_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('pressure [hPa]'); xlim([700 1020])
    subplot(153);plot(qair_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('q [g/kg]'); xlim([0 21])
    title(['P3 flight on ',datestr(tsecs_P3(1,ii))])
    subplot(154);plot(dD_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\deltaD [‰]'); xlim([-102 -40])
    subplot(155);plot(d18O_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\delta^1^8O [‰]'); xlim([-21 0])   
    saveas(gcf,['P3-profiles-p',num2str(ii)],'fig')
    saveas(gcf,['P3-profiles-p',num2str(ii)],'png')
end
%%
for ii = 1:n
    figure;
    subplot(141);plot(thw_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\theta_w [K]'); ylabel('alt [m]'); xlim([290 310]);
    subplot(142);plot(qair_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('q [g/kg]'); xlim([0 21])
    title(['P3 flight on ',datestr(tsecs_P3(1,ii))])
    subplot(143);plot(dD_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\deltaD [‰]'); xlim([-102 -40])
    subplot(144);plot(d18O_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\delta^1^8O [‰]'); xlim([-21 0])   
    saveas(gcf,['P3-profiles-thw',num2str(ii)],'fig')
    saveas(gcf,['P3-profiles-thw',num2str(ii)],'png')
end

%% Plotting all P3 profiles in single plot - script edited on 10/17/2022 %%
figure;
for ii = 1:n
    subplot(141); hold on;
    plot(thw_P3(:,ii),alt_P3(:,ii),'Color',[211/255 211/255 211/255])
    xlabel('\theta_w [K]'); ylabel('alt [m]');
    subplot(142); hold on;
    plot(qair_P3(:,ii),alt_P3(:,ii),'Color',[211/255 211/255 211/255])
    xlabel('q [g/kg]');
    title('P3 flights')
    subplot(143); hold on;
    plot(dD_P3(:,ii),alt_P3(:,ii),'Color',[211/255 211/255 211/255])
    xlabel('\deltaD [‰]');
    subplot(144); hold on;
    plot(d18O_P3(:,ii),alt_P3(:,ii),'Color',[211/255 211/255 211/255])
    xlabel('\delta^1^8O [‰]');
end
saveas(gcf,'P3-all-profiles-thw','fig')
saveas(gcf,'P3-all-profiles-thw','png')

%%
figure;
for ii = 1:n
    subplot(151); hold on;
    plot(theta_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\theta [K]'); ylabel('alt [m]');
    subplot(152); hold on;
    plot(thw_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\theta_w [K]');
    subplot(153); hold on;
    plot(qair_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('q [g/kg]');
    title('P3 flights')
    subplot(154); hold on;
    plot(dD_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\deltaD [‰]');
    subplot(155); hold on;
    plot(d18O_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\delta^1^8O [‰]');   
end
saveas(gcf,'P3-all-profiles-thetas','fig')
saveas(gcf,'P3-all-profiles-thetas','png')

%%
figure;
for ii = 1:n
    subplot(151); hold on;
    plot(theta_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\theta [K]'); ylabel('alt [m]');
    subplot(152); hold on;
    plot(p_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('pressure [hPa]');   
    subplot(153); hold on;
    plot(qair_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('q [g/kg]');
    title('P3 flights')
    subplot(154); hold on;
    plot(dD_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\deltaD [‰]');
    subplot(155); hold on;
    plot(d18O_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\delta^1^8O [‰]');   
end
saveas(gcf,'P3-all-profiles-p','fig')
saveas(gcf,'P3-all-profiles-p','png')

% saveas(gcf,'P3-all-flights-mmr','fig')
% saveas(gcf,'P3-all-flights-mmr','png')

%% Computing means at different altitudes for a quick look plot %%
count = 0;
m_alt_P3 = 5:10:2995;
for k = 5:10:2995
    count = count + 1;
    m_thw_P3(count)  = mean(thw_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
    m_qair_P3(count) = mean(qair_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
    m_dD_P3(count)   = mean(dD_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
    m_d18O_P3(count) = mean(d18O_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
end
subplot(141); hold on; plot(m_thw_P3,m_alt_P3,'k','LineWidth',2)
subplot(142); hold on; plot(m_qair_P3,m_alt_P3,'k','LineWidth',2)
subplot(143); hold on; plot(m_dD_P3,m_alt_P3,'k','LineWidth',2)
subplot(144); hold on; plot(m_d18O_P3,m_alt_P3,'k','LineWidth',2)