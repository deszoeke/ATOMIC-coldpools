%% Comparing P3 profiles isotopic ratios (dD_ent,d18O_ent) with modeled (extrapolation approach) isotopic ratios 
% Ratios directly from P3 profiles
addpath 'D:\BackUp_Personal_Laptop_March2022\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files\Dean''s P3 files with flags'
date = ['0117';'0119';'0123';'0124';'0131';'0203';'0204';'0205';'0209';'0210';'0211']; % all flight dates
n = length(date);
count = 0;
q_P3     = 99999*ones(1000,n);
dD_P3    = 99999*ones(1000,n);
Ta_P3    = 99999*ones(1000,n);
d18O_P3  = 99999*ones(1000,n);
alt_P3   = 99999*ones(1000,n);
tsecs_P3 = 99999*ones(1000,n);
theta_P3 = 99999*ones(1000,n);
thw_P3   = 99999*ones(1000,n); % theta_w
p_P3     = 99999*ones(1000,n); 
counter = 0;

for ii = 1:n
    disp(ii)
    filename = ['P3_with_vertlocflags_2020',date(ii,:),'.nc'];
    flag_prf = ncread(filename,'flag_prf'); % Flag for vertical profiling. When in a profile, flag value is the profile number (1,2,...) for this flight. Flag is NAN otherwise.
    pn = max(flag_prf); % total number of profiles in this flight
    t_fl = ncread(filename,'secs_since_2020Jan01'); % time of the P3 data in seconds
    t_fl = t_fl/3600/24 + datenum('20200101','yyyymmdd');
    Ta_fl = ncread(filename,'Ta'); % air temperature [in degrees K]
    P_fl = ncread(filename,'press'); % corrected static pressure: hPa
    RH_fl = ncread(filename,'RH'); % relative humidity; units: percentage
    q_fl= ncread(filename,'mmr'); % humidity mixing ratio; units: g/kg
    % ratio of the mass of water vapor to the mass of dry air; water vapor mass (dry) mixing ratio
    e = 2.71828; % Euler number
    T = Ta_fl - 273.15;  % air temperature [in degrees C]
    % q_fl = (6.112 .* e.^((17.67.*T)./(T+243.5)) .* RH_fl/100 .* 622.18)./(P_fl-(6.112 .* e.^((17.67.*T)./(T+243.5)) .* RH_fl/100 .* 0.37782));
    % specific humidity; units: g/kg; P [hPa]; T [C]; RH [%];
    qs_fl = ncread(filename,'qs'); % saturated specific humidity; units: kg/kg
    Rd = 287.04; % units???
    Cp = 1005.7; % units???
    th_fl = (Ta_fl).*(1e5./(P_fl*100)).^(Rd/Cp); % (Potential Temp in degrees K); theta(p[Pa],Temp[K])
    addpath('C:\Users\quinones\Documents\Data\thermo')
%   qsat_fl = qs(1e3*1e2,th_fl-273.15); % qsat in kg/kg; ERROR: gives really moist qsat for high alt, when it should be really dry up there    
    thw_fl = theta_w(P_fl*100,Ta_fl,qs_fl); % theta_w(p[Pa],Temp[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008; % theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
    dD_fl = ncread(filename,'dD'); % 
    d18O_fl = ncread(filename,'d18O'); %
    alt_fl =  ncread(filename,'alt'); %
    disp(pn)
    for kk = 1:pn
        dum = alt_fl(flag_prf==kk);
        dummy = dum(2:end);
        T5 = dummy(abs(diff(Ta_fl(flag_prf==kk)))>=.5);
        if isempty(T5)
            count = count + 1;
            q    = q_fl(flag_prf==kk);
            dD   = dD_fl(flag_prf==kk);
            Tair   = Ta_fl(flag_prf==kk);
            d18O = d18O_fl(flag_prf==kk);
            alt  = alt_fl(flag_prf==kk);
            tsecs = t_fl(flag_prf==kk);
            theta = th_fl(flag_prf==kk);
            thw  = thw_fl(flag_prf==kk);
            p    = P_fl(flag_prf==kk);
        elseif T5(1) <= 1500
            continue
        else
            count = count + 1;
            q    = q_fl(flag_prf==kk);
            dD   = dD_fl(flag_prf==kk);
            Tair   = Ta_fl(flag_prf==kk);
            d18O = d18O_fl(flag_prf==kk);
            alt  = alt_fl(flag_prf==kk);
            tsecs = t_fl(flag_prf==kk);
            theta = th_fl(flag_prf==kk);
            thw  = thw_fl(flag_prf==kk);
            p    = P_fl(flag_prf==kk);
        end

        %% Air type: observed %%
        % q and theta @ 400m
        %     qs = ncread(filename,'qs'); % saturated specific humidity; units: kg/kg
        %     wspd = ncread(filename,'ws'); % wind speed; units: m/s    
        zrf = 400; % [in meters]; reference height
        dist = 10; % window
        q_ob(count) = mean(q(alt>=zrf-dist & alt<=zrf+dist),'omitnan');
        th_ob(count) = mean(theta(alt>=zrf-dist & alt<=zrf+dist),'omitnan');
        tsecs_ob(count) = mean(tsecs(alt>=zrf-dist & alt<=zrf+dist),'omitnan');
        d18O_ob(count) = mean(d18O(alt>=zrf-dist & alt<=zrf+dist),'omitnan');
        dD_ob(count)   = mean(dD(alt>=zrf-dist & alt<=zrf+dist),'omitnan');
        de_ob(count)   = dD_ob(count) - (8*d18O_ob(count)); % deuterium excess
        %% Air type: surface %%
        % q and theta @ 1m [mean from the PSD surface data set]
        load '1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' sst slp Ta
        SLP = mean(slp,'omitnan'); % [in hPa]
        SST = mean(sst,'omitnan'); % [in degrees C]
%         [q_surf(count),th_surf(count)] = air_types(0,0,0,'surface',SLP,SST);
        Ts = mean(Ta,'omitnan')+273.15; % in degrees Kelvin
%         [dD_surf(count),d18O_surf(count)] = air_types_iso(0,0,'surface',Ts);
%         de_surf(count) = dD_surf(count) - (8*d18O_surf(count)); % deuterium excess
        %% Air type: entrained %%
        % q and theta @ 1km
        h_low = [ 800  900 1000 1100 1200];
        h_hi  = [1000 1100 1200 1300 1400];
        hrf = 1e3;
         q_ent_1km(count) = mean(    q(alt>=hrf-20 & alt<=hrf+20),'omitnan'); % check units for q_inth, looks like is cg/kg
        th_ent_1km(count) = mean(theta(alt>=hrf-20 & alt<=hrf+20),'omitnan');
        for k = 1:length(h_hi)
             q_ent(k,count) = mean( q(alt>=h_low(k) & alt<=h_hi(k),:),'omitnan'); % check units for q_inth, looks like is cg/kg
            th_ent(k,count) = mean(theta(alt>=h_low(k) & alt<=h_hi(k),:),'omitnan');
            % Incorporating isotope data %
            dD_ent(k,count) = mean(  dD(alt>=h_low(k) & alt<=h_hi(k),:),'omitnan');
          d18O_ent(k,count) = mean(d18O(alt>=h_low(k) & alt<=h_hi(k),:),'omitnan');
            de_ent(k,count) = dD_ent(k,count) - (8*d18O_ent(k,count)); % deuterium excess
        end
        %% Air type: downdraft %%
        % from mean cloud layer
        % q and theta @ mean theta_w (wet-bub potential temp) above 1km and below the trade inversion line (6g/kg contour) for each sounding
        % Extracting trade inversion height (mixed layer depth)
        h6 = double(alt(q<=6));
        if isempty(h6)
            disp(['[',num2str(count),']'])
            count = count - 1;
            continue
        end
        if h6(1) <= 6000
            trade_inv(count) = h6(1); % trade inversion
        else
            trade_inv(count) = NaN; % trade inversion
        end
        thw_index = find(alt>=hrf-10 & alt<=trade_inv(count));
        if isempty(thw_index)
            disp(['[',num2str(count),']'])
            count = count - 1;
            continue
        end
        thw_index = thw_index(1);
        hrf_alt(count) = alt(thw_index);
        if hrf_alt(count) == trade_inv(count)
            disp(['[',num2str(count),']'])
            count = count - 1;
            continue
        end
        if trade_inv(count) >= 1e3
            counter = counter + 1;
            trade_index = find(alt==trade_inv(count));
            th_dd(count) = nanmean(thw(thw_index:trade_index));
        else
            th_dd(count)   = NaN;
        end
        clearvars h6 thw_index trade_index
        addpath('C:\Users\quinones\Documents\Data\thermo')
        q_dd(count) = qs(1e3*1e2,th_dd(count)-273.15)*1e3; % q_d in g/kg    
        % qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
        % p [Pa], T [degrees C], qs [kg/kg]
    end
    q_P3(1:length(q),count)       = q;
    dD_P3(1:length(q),count)      = dD;
    Ta_P3(1:length(q),count)      = Tair;
    d18O_P3(1:length(q),count)    = d18O;
    alt_P3(1:length(q),count)     = alt;
    tsecs_P3(1:length(q),count)   = tsecs;
    theta_P3(1:length(q),count)   = theta;
    thw_P3(1:length(q),count)     = thw;
    p_P3(1:length(q),count)       = p;

    clearvars q dD Tair d18O alt tsecs theta thw p    
end

q_P3(q_P3==99999)       = NaN;
dD_P3(dD_P3==99999)     = NaN;
d18O_P3(d18O_P3==99999) = NaN;
Ta_P3(Ta_P3==99999)     = NaN;
tsecs_P3(tsecs_P3==99999) = NaN;
theta_P3(theta_P3==99999) = NaN;
thw_P3(thw_P3==99999)   = NaN;
p_P3(p_P3==99999)       = NaN;
alt_P3(alt_P3==99999)   = NaN;
q_dd(length(q_dd)+1)    = NaN;
th_dd(length(th_dd)+1)  = NaN;

2+2;
%% Plots in dD-theta space %%
figure;
plot(dD_ob,th_ob,'sk','MarkerSize',5,'MarkerFaceColor','k')
for k = 1:size(dD_ent,1)
%     figure;
%     plot(dD_ob,th_ob,'sk','MarkerSize',5,'MarkerFaceColor','k')
%     hold on;
%     scatter(dD_ob,de_ob,13,th_ob,'s','filled')
%     colormap(jet(15))
%     colorbar;
    hold on;
    scatter(mean(dD_ent_mf(k,:),'omitnan'),mean(th_ent(k,:),'omitnan'),50,'filled')
%     scatter(dD_ent_mf(k,:),de_ent_mf(k,:),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(dD_dd_mf(k,:),'omitnan'),mean(th_dd_mf(k,:),'omitnan'),50,[0 0.4470 0.7410],'filled')
%     scatter(dD_dd_mf(k,1:24:end),de_dd_mf(k,1:24:end),15,[0 0.4470 0.7410],'filled')
    scatter(mean(dD_ent(k,:),'omitnan'),mean(th_ent(k,:),'omitnan'),25,'black','filled')
%     scatter(dD_ent(k,:),de_ent(k,:),15,'black','filled')
%     scatter(mean(dD_surf,'omitnan'),mean(th_surf,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
%     scatter(dD_surf,de_surf,15,[0.9290 0.6940 0.1250],'filled')
%     scatter(dD_surf,th_surf,50,'r','filled')
%     xlim([-84 -56]) 
%     ylim([-4 42])
%     ylabel(['DXS [',char(8240),']'])
    ylabel('\theta [K]')
    xlabel(['\deltaD [',char(8240),']'])
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
%     hold on;text(-66,1,'surface','FontSize',15)
    hold on;text(-75,21.5,'entrained','FontSize',15)
%     hold on;text(-60,11.5,'downdraft','FontSize',15)
%     caxis([297 300])
    legend('P3 SBL obs','0.8-1km','','0.9-1.1km','','1-1.2km','','1.1-1.3km','','1.2-1.4km','P3 layer obs','surface','FontSize',15)
%     set(gca, 'YDir','reverse')
%     set(gca, 'XDir','reverse')
%     title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
%     text(-58,0,['RMSE = ',num2str(rmse(k)),''],'FontSize',15)
%     saveas(gcf,['P3-28-profiles-DXS-dD-dd-height-range',num2str(k)],'png')
%     saveas(gcf,['P3-28-profiles-DXS-dD-dd-height-range',num2str(k)],'fig')
end
hold on;
scatter(mean(dD_surf,'omitnan'),mean(th_surf,'omitnan'),50,'r','filled')
legend('P3 SBL obs','0.8-1km','','0.9-1.1km','','1-1.2km','','1.1-1.3km','','1.2-1.4km','P3 layer obs','surface','FontSize',15)

%% Plots in th-q space %%
for k = 1:size(dD_ent,1)
    figure;
    plot(th_ob,q_ob,'sk','MarkerSize',5)
    hold on;
    scatter(th_ob,q_ob,13,dD_ob,'s','filled')
%     colormap(flip(jet(18))) % for DXS % for comparison with RHB data
    colormap(flip(jet(16))) % FOR dD  % for comparison with RHB data
    colorbar;
    hold on;
    scatter(th_ent(k,:),q_ent(k,:),10,[0.8500 0.3250 0.0980],'filled')
    hold on;
    scatter(th_surf,q_surf,10,[0.9290 0.6940 0.1250],'filled')
    hold on;
    scatter(th_dd,q_dd,10,[0 0.4470 0.7410],'filled')
    xlim([289 302]) % for comparison with RHB data
%     xlim([296 302])
    ylim([8 22])
    ylabel('q [g/kg]')
    xlabel('\theta [K]') 
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
    hold on;text(297,21,'surface','FontSize',15)
    hold on;text(295,10,'entrained','FontSize',15)
    hold on;text(291,14,'downdraft','FontSize',15)
    caxis([-81 -65]) % FOR dD % for comparison with RHB data
%     caxis([-79 -67]) % FOR dD 
%     caxis([5 14]) % FOR DXS % for comparison with RHB data
    title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
    saveas(gcf,['P3-60-profiles-dD-height-range',num2str(k)],'png')
    saveas(gcf,['P3-60-profiles-dD-height-range',num2str(k)],'fig')
end

%% Plotting timeseries %%
fl_num = 1:length(q_ob);
figure;
subplot(511)
    plot(fl_num,q_ob,'k',fl_num,q_ent(4,:),fl_num,q_surf)
    legend('observed','entrained','surface')
    ylabel('q [g/kg]')
subplot(512)
    plot(fl_num,th_ob,'k',fl_num,th_ent(4,:),fl_num,th_surf)
    ylabel('\theta [K]')
subplot(513)
    plot(fl_num,dD_ob,'k',fl_num,dD_ent(4,:),fl_num,dD_surf)
    ylabel(['\deltaD [',char(8240),']'])
subplot(514)
    plot(fl_num,d18O_ob,'k',fl_num,d18O_ent(4,:),fl_num,d18O_surf)
    ylabel(['\delta^1^8O [',char(8240),']'])
subplot(515)
    plot(fl_num,de_ob,'k',fl_num,de_ent(4,:),fl_num,de_surf)
    ylabel(['DXS [',char(8240),']'])
    xlabel('In-house P3 profile #')
set(findall(gcf,'-property','Fontsize'),'FontSize',15)
set(findall(gcf,'-property','TickLength'),'TickLength',[.01 .01])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)

%% Ratios (dD_ent,d18O_ent) based on extrapolation approach using exclusively P3 data!
for k = 1:length(h_hi)
    [dD_ent_mf(k,:),d18O_ent_mf(k,:)] = iso_estimates_extrapolation_approach(q_ob,dD_ob,d18O_ob,q_surf,dD_surf,d18O_surf,q_ent(k,:));
    % RMSE %
    rmse_dD(k)   = sqrt(sum((dD_ent(k,:)-dD_ent_mf(k,:)).^2,'omitnan')/sum(~isnan(dD_ent_mf(k,:))));
    rmse_d18O(k) = sqrt(sum((d18O_ent(k,:)-d18O_ent_mf(k,:)).^2,'omitnan')/sum(~isnan(dD_ent_mf(k,:))));
end
de_ent_mf = dD_ent_mf - (8*d18O_ent_mf); % deuterium excess
for k = 1:length(h_hi)
    rmse_de(k) = sqrt(sum((de_ent(k,:)-de_ent_mf(k,:)).^2,'omitnan')/sum(~isnan(de_ent_mf(k,:))));
end

%% Plotting comparison between ent timeseries %%
% figure('units','normalized','outerposition',[0 0 1 1])
for k = 1:size(dD_ent,1)
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(311)
        hold on;
        plot(tsecs_ob,dD_ent(k,:),'.k',tsecs_ob,dD_ent_mf(k,:),'o',MarkerFaceColor='auto')
        ylabel(['\deltaD_e_n_t[',char(8240),']'])
        datetick('x','mmm-dd','keeplimits','keepticks')
        ylim([-94 -69])
        title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
        legend('from P3 profile','from mixing fractions','orientation','horizontal','location','southeast')
        text([datenum('01/15/2020','mm/dd/yyyy')],-67,['RMSE = ',num2str(rmse_dD(k)),''],'FontSize',15)
    subplot(312)
        hold on;
        plot(tsecs_ob,d18O_ent(k,:),'.k',tsecs_ob,d18O_ent_mf(k,:),'o',MarkerFaceColor='auto')
        ylabel(['\delta^1^8O_e_n_t [',char(8240),']'])
        datetick('x','mmm-dd','keeplimits','keepticks')
        ylim([-17 -9.5])
        text([datenum('01/15/2020','mm/dd/yyyy')],-9,['RMSE = ',num2str(rmse_d18O(k)),''],'FontSize',15)
    subplot(313)
        hold on;
        plot(tsecs_ob,de_ent(k,:),'.k',tsecs_ob,de_ent_mf(k,:),'o',MarkerFaceColor='auto')
        ylabel(['DXS_e_n_t [',char(8240),']'])
        datetick('x','mmm-dd','keeplimits','keepticks')
        xlabel('date in 2020')
        ylim([6 50])
        text([datenum('01/15/2020','mm/dd/yyyy')],53,['RMSE = ',num2str(rmse_de(k)),''],'FontSize',15)
    set(findall(gcf,'-property','Fontsize'),'FontSize',15)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.01 .01])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    saveas(gcf,['P3-28-profiles-ent-timeseries-comparison',num2str(k)],'png')
    saveas(gcf,['P3-28-profiles-ent-timeseries-comparison',num2str(k)],'fig')
end
% legend('[0.8 1.0] km','mixing fractions',...
%        '[0.9 1.1] km','mixing fractions',...
%        '[1.0 1.2] km','mixing fractions',...
%        '[1.1 1.3] km','mixing fractions',...
%        '[1.2 1.4] km','mixing fractions')
% saveas(gcf,'P3-28-profiles-ent-timeseries-comparison2','png')
% saveas(gcf,'P3-28-profiles-ent-timeseries-comparison2','fig')

%% Plotting mean diff and RMSD %%
figure;
subplot(311)
    plot(h_low,mean(dD_ent,2),'.k',h_low,mean(dD_ent_mf,2,'omitnan'),'o',MarkerFaceColor='auto')
    hold on; plot(h_low,mean(dD_ent,2),'.k');
    ylabel(['\deltaD_e_n_t[',char(8240),']'])
subplot(312)
    plot(h_low,mean(d18O_ent,2),'.k',h_low,mean(d18O_ent_mf,2,'omitnan'),'o',MarkerFaceColor='auto')
    ylabel(['\delta^1^8O_e_n_t [',char(8240),']'])
subplot(313)
    plot(h_low,mean(de_ent,2),'.k',h_low,mean(de_ent_mf,2,'omitnan'),'o',MarkerFaceColor='auto')
    ylabel(['DXS_e_n_t [',char(8240),']'])
    xlabel('height [m]')
figure;
subplot(311)
    plot(h_low,mean(dD_ent-dD_ent_mf,2,'omitnan'),'.k')
    hold on; plot(h_low,[0,0,0,0,0],'-k')
    ylabel(['\deltaD_e_n_t[',char(8240),']'])
subplot(312)
    plot(h_low,mean(d18O_ent-d18O_ent_mf,2,'omitnan'),'.k')
    ylabel(['\delta^1^8O_e_n_t [',char(8240),']'])
    ylim([0 2])
subplot(313)
    plot(h_low,mean(de_ent-de_ent_mf,2,'omitnan'),'.k')
    ylabel(['DXS_e_n_t [',char(8240),']'])
    xlabel('height [m]')
% RMSD %
rmsd = rms(dD_ent-dD_ent_mf,2,'omitnan');
rmsd = rms(d18O_ent-d18O_ent_mf,2,'omitnan');
rmsd = rms(de_ent-de_ent_mf,2,'omitnan');
figure;
subplot(311)
    plot(h_low,rms(dD_ent-dD_ent_mf,2,'omitnan'),'.k')
    ylabel(['\deltaD_e_n_t[',char(8240),']'])
subplot(312)
    plot(h_low,rms(d18O_ent-d18O_ent_mf,2,'omitnan'),'.k')
    ylabel(['\delta^1^8O_e_n_t [',char(8240),']'])
%     ylim([0 2])
subplot(313)
    plot(h_low,rms(de_ent-de_ent_mf,2,'omitnan'),'.k')
    ylabel(['DXS_e_n_t [',char(8240),']'])
    xlabel('height [m]')

%% Ratios (dD_dd,d18O_dd) based on mixing fractions using exclusively P3 data!
% Isotope concentrations calculated based on mixing fractions obtained from P3 profiles
for k = 1:length(h_hi)
    [fen(k,:), fss(k,:), fdd(k,:)] = th_q_to_mixfraction(th_ob,q_ob, th_ent(k,:),q_ent(k,:), th_surf,q_surf, th_dd,q_dd);
end
fdd(fdd<0)=NaN;
fss(fss<0)=NaN;
fen(fen<0)=NaN;
for k = 1:length(h_hi)
    dD_dd_mf(k,:) = (q_ob.*dD_ob - q_ent(k,:).*fen(k,:).*dD_ent_mf(k,:) - q_surf.*fss(k,:).*dD_surf)./(q_dd.*fdd(k,:));
    d18O_dd_mf(k,:) = (q_ob.*d18O_ob - q_ent(k,:).*fen(k,:).*d18O_ent_mf(k,:) - q_surf.*fss(k,:).*d18O_surf)./(q_dd.*fdd(k,:));
    de_dd_mf(k,:) = dD_dd_mf(k,:) - (8*d18O_dd_mf(k,:)); % deuterium excess
end
%% Plotting mean diff and RMSD %%
figure;
subplot(311)
    plot(h_low,mean(dD_ent,2),'.k',h_low,mean(dD_ent_mf,2,'omitnan'),'o',MarkerFaceColor='auto')
    hold on; plot(h_low,mean(dD_ent,2),'.k');
    ylabel(['mean [',char(8240),']'])
    title(['\deltaD_e_n_t[',char(8240),']'])
subplot(312)
    plot(h_low,mean(dD_ent-dD_ent_mf,2,'omitnan'),'.k')
    hold on; plot(h_low,[0,0,0,0,0],'-k')
    ylabel(['obs - est [',char(8240),']'])
subplot(313)
    plot(h_low,rms(dD_ent-dD_ent_mf,2,'omitnan'),'.k')
    ylabel(['RMSD [',char(8240),']'])
    xlabel('height [m]')
figure;
subplot(311)
    plot(h_low,mean(d18O_ent,2),'.k',h_low,mean(d18O_ent_mf,2,'omitnan'),'o',MarkerFaceColor='auto')
    ylabel(['mean [',char(8240),']'])
    title(['\delta^1^8O_e_n_t [',char(8240),']'])
subplot(312)
    plot(h_low,mean(d18O_ent-d18O_ent_mf,2,'omitnan'),'.k')
    hold on; plot(h_low,[0,0,0,0,0],'-k')
    ylabel(['obs - est [',char(8240),']'])
    ylim([0 2])
subplot(313)
    plot(h_low,rms(d18O_ent-d18O_ent_mf,2,'omitnan'),'.k')
    ylabel(['RMSD [',char(8240),']'])
    xlabel('height [m]')
figure;
subplot(311)
    plot(h_low,mean(de_ent,2),'.k',h_low,mean(de_ent_mf,2,'omitnan'),'o',MarkerFaceColor='auto')
    ylabel(['mean [',char(8240),']'])
    title(['DXS_e_n_t [',char(8240),']'])
subplot(312)
    plot(h_low,mean(de_ent-de_ent_mf,2,'omitnan'),'.k')
    hold on; plot(h_low,[0,0,0,0,0],'-k')
    ylabel(['obs - est [',char(8240),']'])
    xlabel('height [m]')
subplot(313)
    plot(h_low,rms(de_ent-de_ent_mf,2,'omitnan'),'.k')
    ylabel(['RMSD [',char(8240),']'])
    xlabel('height [m]')

%% Plotting dd timeseries %%
for k = 1:size(dD_ent,1)
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(311)
        plot(tsecs_ob,dD_dd_mf(k,:),'o',MarkerFaceColor='auto')
        ylabel(['\deltaD_d_d [',char(8240),']'])
        datetick('x','mmm-dd','keeplimits','keepticks')
        title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
    subplot(312)
        plot(tsecs_ob,d18O_dd_mf(k,:),'o',MarkerFaceColor='auto')
        ylabel(['\delta^1^8O_d_d [',char(8240),']'])
        datetick('x','mmm-dd','keeplimits','keepticks')
    subplot(313)
        plot(tsecs_ob,de_dd_mf(k,:),'o',MarkerFaceColor='auto')
        ylabel(['DXS_d_d [',char(8240),']'])
        datetick('x','mmm-dd','keeplimits','keepticks')
        xlabel('date in 2020')
    set(findall(gcf,'-property','Fontsize'),'FontSize',15)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.01 .01])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    saveas(gcf,['P3-28-profiles-dd-timeseries',num2str(k)],'png')
    saveas(gcf,['P3-28-profiles-dd-timeseries',num2str(k)],'fig')
end
%% Plotting comparison between surf timeseries %%
figure;
subplot(311)
    plot(tsecs_ob,dD_surf,'ok',tsecs_ob,dD_surf_mf(4,:),'o',MarkerFaceColor='auto')
    ylabel(['\deltaD_s_u_r_f [',char(8240),']'])
    datetick('x','mmm-dd','keeplimits','keepticks')
    legend('from P3 profile','from mixing fractions')
subplot(312)
    plot(tsecs_ob,d18O_surf,'ok',tsecs_ob,d18O_surf_mf(4,:),'o',MarkerFaceColor='auto')
    ylabel(['\delta^1^8O_s_u_r_f [',char(8240),']'])
    datetick('x','mmm-dd','keeplimits','keepticks')
subplot(313)
    plot(tsecs_ob,de_surf,'ok',tsecs_ob,de_surf_mf(4,:),'o',MarkerFaceColor='auto')
    ylabel(['DXS_s_u_r_f [',char(8240),']'])
    datetick('x','mmm-dd','keeplimits','keepticks')
    xlabel('date in 2020')
set(findall(gcf,'-property','Fontsize'),'FontSize',15)
set(findall(gcf,'-property','TickLength'),'TickLength',[.01 .01])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
saveas(gcf,'P3-28-profiles-surf-timeseries-comparison','png')
saveas(gcf,'P3-28-profiles-surf-timeseries-comparison','fig')

%% Including data points in time-height plots from Oct 13 2022 %%
% open('time_height_P3_iso_data.fig') % dD plot
% open('time_height_P3_iso_data_18O.fig')
% open('time_height_P3_iso_data_DXS.fig')
hold on;
plot(tsecs_ob(1:end-1),trade_inv,'ok',tsecs_ob(1:end-1),hrf_alt,'o',MarkerFaceColor='auto')
ylim([0 3000])
