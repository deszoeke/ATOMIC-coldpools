%% pdf from P3 data
% addpath 'D:\BackUp_Personal_Laptop_March2022\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files\Dean''s P3 files with flags'
% addpath 'C:\Users\quinones\Documents\Data\thermo'
% adri_data = ncinfo('P3_ATOMIC_RF01_Isotopes1Hz_V0.1.nc');
adri_data = ncinfo('data/EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-1Hz_20200117_v1.1.nc');
p3_data = ncinfo('data/EUREC4A_ATOMIC_P3_Flight-Level_20200117_v1.1.nc');
% Variables: Time dD d18O
dean_flags = ncinfo('data/P3_with_vertlocflags_20200210.nc');
% Variables: Time dD d18O qs thetaEs thetav RH Ta ws wd alt
%% Separating profiles for each flight %%
% date = ['0131'; '0203'; '0204'; '0205'; '0209'; '0210']; % dates that coincide with cold pools
date = ['0117';'0119';'0123';'0124';'0131';'0203';'0204';'0205';'0209';'0210';'0211']; % all flight dates
n = length(date);
hlevel = 2000; % in meters
count = 0;

qair_P3  = zeros(1000,n);
qs_P3    = zeros(1000,n);
dD_P3    = zeros(1000,n);
Ta_P3    = zeros(1000,n);
d18O_P3  = zeros(1000,n);
alt_P3   = zeros(1000,n);
tsecs_P3 = zeros(1000,n);
theta_P3 = zeros(1000,n);
thw_P3   = zeros(1000,n); % theta_w
p_P3     = zeros(1000,n); 

for ii = 1:n
    disp(ii)
    % temp and qs variables
    filename = ['P3_with_vertlocflags_2020',date(ii,:),'.nc'];
    flag_prf = ncread(filename,'flag_prf'); % Flag for vertical profiling. When in a profile, flag value is the profile number (1,2,...) for this flight. Flag is NAN otherwise.
    pn = max(flag_prf); % total number of profiles in this flight
    tsecs = ncread(filename,'secs_since_2020Jan01'); % time of the P3 data
    tsecs = tsecs/3600/24 + datenum('20200101','yyyymmdd');
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

    for kk = 1:pn
        % Skipping profiles for which trade inversion is below 2km
        % h6 = alt(qair<=6);
        dum = alt(flag_prf==kk);
        dummy = dum(2:end);
        h6 = dummy(abs(diff(Ta(flag_prf==kk)))>=.5);
        if isempty(h6)
            count = count + 1;
            qair_P3(1:length(qair(flag_prf==kk)),count) = qair(flag_prf==kk);
            qs_P3(1:length(qair(flag_prf==kk)),count)   = qs(flag_prf==kk);
            dD_P3(1:length(qair(flag_prf==kk)),count)   = dD(flag_prf==kk);
            Ta_P3(1:length(qair(flag_prf==kk)),count)   = Ta(flag_prf==kk);
            d18O_P3(1:length(qair(flag_prf==kk)),count) = d18O(flag_prf==kk);
            alt_P3(1:length(qair(flag_prf==kk)),count)  = alt(flag_prf==kk);
            tsecs_P3(1:length(qair(flag_prf==kk)),count) = tsecs(flag_prf==kk);
            theta_P3(1:length(qair(flag_prf==kk)),count) = th(flag_prf==kk);
            thw_P3(1:length(qair(flag_prf==kk)),count)  = thw(flag_prf==kk);
            p_P3(1:length(qair(flag_prf==kk)),count)    = p(flag_prf==kk);
        elseif h6(1) <= 1500
            continue
        else
            count = count + 1;
            qair_P3(1:length(qair(flag_prf==kk)),count) = qair(flag_prf==kk);
            qs_P3(1:length(qair(flag_prf==kk)),count)   = qs(flag_prf==kk);
            dD_P3(1:length(qair(flag_prf==kk)),count)   = dD(flag_prf==kk);
            Ta_P3(1:length(qair(flag_prf==kk)),count)   = Ta(flag_prf==kk);
            d18O_P3(1:length(qair(flag_prf==kk)),count) = d18O(flag_prf==kk);
            alt_P3(1:length(qair(flag_prf==kk)),count)  = alt(flag_prf==kk);
            tsecs_P3(1:length(qair(flag_prf==kk)),count) = tsecs(flag_prf==kk);
            theta_P3(1:length(qair(flag_prf==kk)),count) = th(flag_prf==kk);
            thw_P3(1:length(qair(flag_prf==kk)),count)  = thw(flag_prf==kk);
            p_P3(1:length(qair(flag_prf==kk)),count)    = p(flag_prf==kk);
        end
    end

    clearvars qair dD d18O alt th tsecs Ta thw qs p pn kk flag_prf h6 dummy dum
end

qair_P3(qair_P3==0) = NaN;
qs_P3(qs_P3==0)     = NaN;
dD_P3(dD_P3==0)     = NaN;
d18O_P3(d18O_P3==0) = NaN;
Ta_P3(Ta_P3==0)     = NaN;
tsecs_P3(tsecs_P3==0) = NaN;
theta_P3(theta_P3==0) = NaN;
thw_P3(thw_P3==0)     = NaN;
p_P3(p_P3==0)       = NaN;
alt_P3(alt_P3<=0)   = NaN;

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

%% Plots scripted on 10/13/2022 %%
% time-height plots to determine at which level depleted isotopes are
% present!!!

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
cmap = b2rcolormap(25);
colormap(cmap)
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

%% Plotting P3 profiles -  script edited on 10/17/2022 %%
for ii = 1:n
    figure;
    subplot(151);plot(theta_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\theta [K]'); ylabel('alt [m]')
    subplot(152);plot(p_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('pressure [hPa]');
    subplot(153);plot(qair_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('q [g/kg]');
    title(['P3 flight on ',datestr(tsecs_P3(1,ii))])
    subplot(154);plot(dD_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\deltaD [‰]');
    subplot(155);plot(d18O_P3(:,ii),alt_P3(:,ii),'k')
    xlabel('\delta^1^8O [‰]');
    saveas(gcf,['P3-profiles-p',num2str(ii)],'fig')
    saveas(gcf,['P3-profiles-p',num2str(ii)],'png')
end

%% Plotting flag P3 profiles in single plot - script edited on 10/18/2022 %%
figure;
for ii = 1:count
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
set(findall(gcf,'-property','YLim'),'YLim',[0 hlevel])
% saveas(gcf,'P3-flag-profiles-thw_zoom','fig')
% saveas(gcf,'P3-flag-profiles-thw_zoom','png')
% xlim([-125 -20])
% xlim([-17 -5])

%% Calculating rates between different height levels - script edited on 06/09/2023
% Variables: dD and d18O
% Level: between 0.1 and 1.1 km
    % For dD
p1 = [mean(dD_P3(alt_P3>=90 & alt_P3<=110),'omitnan') mean(alt_P3(alt_P3>=90 & alt_P3<=110),'omitnan')];
p2 = [mean(dD_P3(alt_P3>=1090 & alt_P3<=1110),'omitnan') mean(alt_P3(alt_P3>=1090 & alt_P3<=1110),'omitnan')];
rise = p2(2) - p1(2);   %  999.5158m = ~1km
run  = p2(1) - p1(1);   % -8.0215 permil
rate_dD = run/(rise/10^3); % -8.0215 permil/km
    % For d18O
p1 = [mean(d18O_P3(alt_P3>=90 & alt_P3<=110),'omitnan') mean(alt_P3(alt_P3>=90 & alt_P3<=110),'omitnan')];
p2 = [mean(d18O_P3(alt_P3>=1090 & alt_P3<=1110),'omitnan') mean(alt_P3(alt_P3>=1090 & alt_P3<=1110),'omitnan')];
rise = p2(2) - p1(2);   %  999.5158m = ~1km
run  = p2(1) - p1(1);   % -0.7463 permil
rate_d18O = run/(rise/10^3); % -0.7463 permil/km
% Level: between 1.1 and 1.3 km (values should remain constant)
    % Corroborating constant value by calculating means:
m_bar = mean(d18O_P3(alt_P3>=1100 & alt_P3<=1300),'omitnan'); % -11.5388
m_bar = mean(  dD_P3(alt_P3>=1100 & alt_P3<=1300),'omitnan'); % -78.9368
% Level: above trade inversion (1.3km) and up to 3 km
    % For dD
p1 = [mean(dD_P3(alt_P3>=1290 & alt_P3<=1310),'omitnan') mean(alt_P3(alt_P3>=1290 & alt_P3<=1310),'omitnan')];
p2 = [mean(dD_P3(alt_P3>=2990 & alt_P3<=3010),'omitnan') mean(alt_P3(alt_P3>=2990 & alt_P3<=3010),'omitnan')];
rise = p2(2) - p1(2);   %  1701.74m = ~1.7km
run  = p2(1) - p1(1);   % -111.1105 permil
rate_dD = run/(rise/10^3); % -65.2920 permil/km
    % For d18O
p1 = [mean(d18O_P3(alt_P3>=1290 & alt_P3<=1310),'omitnan') mean(alt_P3(alt_P3>=1290 & alt_P3<=1310),'omitnan')];
p2 = [mean(d18O_P3(alt_P3>=2990 & alt_P3<=3010),'omitnan') mean(alt_P3(alt_P3>=2990 & alt_P3<=3010),'omitnan')];
rise = p2(2) - p1(2);   %  1701.74m = ~1.7km
run  = p2(1) - p1(1);   % -10.2029 permil
rate_d18O = run/(rise/10^3); % -5.9955 permil/km

%%
figure;
for ii = 1:count
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
set(findall(gcf,'-property','YLim'),'YLim',[0 hlevel])
saveas(gcf,'P3-flag-profiles-thetas','fig')
saveas(gcf,'P3-flag-profiles-thetas','png')

%%
figure;
for ii = 1:count
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
set(findall(gcf,'-property','YLim'),'YLim',[0 hlevel])
saveas(gcf,'P3-flag-profiles-p','fig')
saveas(gcf,'P3-flag-profiles-p','png')

%% Computing means at different altitudes for a quick look plot %%
counter = 0;
m_alt_P3 = 5:10:2995;
for k = 5:10:2995
    counter = counter + 1;
    m_thw_P3(counter)  = mean(thw_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
    m_qair_P3(counter) = mean(qair_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
    m_dD_P3(counter)   = mean(dD_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
    m_d18O_P3(counter) = mean(d18O_P3(alt_P3>=k-24 & alt_P3<=k+25),'omitnan');
end
subplot(141); hold on; plot(m_thw_P3,m_alt_P3,'k','LineWidth',2)
subplot(142); hold on; plot(m_qair_P3,m_alt_P3,'k','LineWidth',2)
subplot(143); hold on; plot(m_dD_P3,m_alt_P3,'k','LineWidth',2)
subplot(144); hold on; plot(m_d18O_P3,m_alt_P3,'k','LineWidth',2)