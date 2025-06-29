% conserved_properties_for_isotopes.m
%
% Applying the extrapolation and mixing fraction analysis to RHB data in
% order to obtain isotope concentrations for entrained and downdraft end members
% EQM Last modified: Nov 16 2022
% SPdeS tweaked some paths 2025-06-05

% addpath('./thermo'); 

%% Air type: observed %%
% q and theta @ 400m
zrf = 400; % [in meters]; reference height
zm = 17;   % [in meters]; measurement level
Rd = 287.04;
Cp = 1005.7;
[th_ob,q_ob,t_adj] = height_adj(zrf,zm,Rd,Cp);
load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' t
pos1 = find(t_adj>=t(1)); % rounding up to closest (in time) PSD surface data point!!!
pos_adj = pos1(1);
clearvars pos1
q_ob = q_ob(pos_adj:pos_adj+length(t)-1);   % keeping 2nd leg data only
th_ob = th_ob(pos_adj:pos_adj+length(t)-1); % keeping 2nd leg data only
t_10min = t;
% Loading iso data %
load('data/iso_data_1min_intervals_FLAGGED_w_runningmean.mat');% RHB Picarro data
load 'data/cold_pool_flag_1min.mat' t1min
ind0 = find(t1min>=iso_time(1));
ind0 = ind0(1);
% Incorporating isotope data %
pos_i = 999999*ones(size(t));
for l = 1:length(t)
    pos2 = find(iso_time>=t(l)); % rounding up to closest isotope surface data point!!!
                                 % try rounding to nearest data point
    pos_i(l) = pos2(1);
    clearvars pos2
end
d18O_ob = d18O(pos_i);
dD_ob   = dD(pos_i);
de_ob   = dD_ob - (8*d18O_ob); % deuterium excess

%% Air type: surface [using properties in equilibrium w/ surface ocean] %%
% q and theta @ 1m
pos = 999999*ones(size(t_10min));
for l = 1:length(t_10min)
    pos1 = find(t1min>=t_10min(l)); % rounding up to closest (in time) PSD surface data point!!!
                                    % try rounding to nearest (in space) data point
    pos(l) = pos1(1);
    clearvars pos1
end
load 'data/1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' sst slp Ta rh
T = Ta(pos(1:end))+273.15; % in degrees Kelvin
SLP = slp(pos(1:end)); % [in hPa]
SST = sst(pos(1:end)); % [in degrees C]
RH  =  rh(pos(1:end));  % in percent
[q_surf,th_surf] = air_types(0,0,0,'surface',SLP,SST);

%% Air type: surface [using Craig-Gordon] %%
% Incorporating isotope data %
h_prime = 0.875; % relative humidity at which we are modeling the iso ratios
filename = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
Tliq = ncread(filename,'tskin'); % liquid temperature at [???]m [in degrees C]
[T_skin,~] = colocate_in_time(filename,Tliq);
[dD_surf,d18O_surf] = air_types_iso(dD_ob,d18O_ob,'surface',SST',RH'/100,h_prime);
de_surf = dD_surf - (8*d18O_surf); % deuterium excess

%% Air type: entrained %%
% q and theta @ 1km
load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' q th h
h_low = [ 800  900 1000 1100 1200];
h_hi  = [1000 1100 1200 1300 1400];
q_ent_1km  = q(h==1e3,:)*1e3; % check units for q_inth, looks like is cg/kg
th_ent_1km = th(h==1e3,:);
for k = 1:length(h_hi)
     q_ent(k,:)  = mean( q(h>=h_low(k) & h<=h_hi(k),:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
    th_ent(k,:) = mean(th(h>=h_low(k) & h<=h_hi(k),:),'omitnan');
    % Incorporating isotope data %
    [dD_ent(k,:),d18O_ent(k,:)] = iso_estimates_from_centroid_approach(q_ob,dD_ob,d18O_ob,q_surf,dD_surf,d18O_surf,q_ent(k,:));
     de_ent(k,:) = dD_ent(k,:) - (8*d18O_ent(k,:)); % deuterium excess
end
2+2;
%% Air type: downdraft %%
% from mean cloud layer
% q and theta @ mean theta_w (wet-bub potential temp) above 1km and below the trade inversion line (6g/kg contour) for each sounding
% Extracting trade inversion height (mixed layer depth)
for k = 12:size(q,2)
    h6 = double(h(q(:,k)*1000<=6));
    if h6(1) <= 6000
        trade_inv(k) = h6(1); % trade inversion
    else
        trade_inv(k) = NaN; % trade inversion
    end
    clearvars h6
end
thw_index = find(h==1e3);
load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' thw
for k = 1:size(thw,2)
    if trade_inv(k) >= 1e3
        trade_index = find(h==trade_inv(k));
        th_dd(k) = nanmean(thw(thw_index:trade_index,k));
    else
        th_dd(k) = NaN; % trade inversion
    end
end
q_dd = qs(1e3*1e2,th_dd-273.15)*1e3; % q_d in g/kg    
%% Isotope data for downdraft end member %%
% Isotope concentrations calculated based on the mixing fractions obtained from soundings
for k = 1:length(h_hi)
    [fen(k,:), fss(k,:), fdd(k,:)] = th_q_to_mixfraction(th_ob,q_ob, th_ent(k,:),q_ent(k,:), th_surf,q_surf, th_dd,q_dd);
    dD_dd(k,:) = (q_ob.*dD_ob' - q_ent(k,:).*fen(k,:).*dD_ent(k,:) - q_surf.*fss(k,:).*dD_surf')./(q_dd.*fdd(k,:));
    d18O_dd(k,:) = (q_ob.*d18O_ob' - q_ent(k,:).*fen(k,:).*d18O_ent(k,:) - q_surf.*fss(k,:).*d18O_surf')./(q_dd.*fdd(k,:));
    de_dd(k,:) = dD_dd(k,:) - (8*d18O_dd(k,:)); % deuterium excess
end

fdd(fdd<0)=NaN;
fss(fss<0)=NaN;
fen(fen<0)=NaN;

for k = 1:length(h_hi)
    dD_dd(k,:) = (q_ob.*dD_ob' - q_ent(k,:).*fen(k,:).*dD_ent(k,:) - q_surf.*fss(k,:).*dD_surf')./(q_dd.*fdd(k,:));
    d18O_dd(k,:) = (q_ob.*d18O_ob' - q_ent(k,:).*fen(k,:).*d18O_ent(k,:) - q_surf.*fss(k,:).*d18O_surf')./(q_dd.*fdd(k,:));
    de_dd(k,:) = dD_dd(k,:) - (8*d18O_dd(k,:)); % deuterium excess
end


%% Cold pool times
load 'data/cold_pool_flag_1min.mat' cold_pool_flag_1min
cold_pool_flag_10min = cold_pool_flag_1min(1:10:end);
ind1 = find(t_adj>=t_10min(1));
ind1 = ind1(1);
iso_cold_pool_flag = cold_pool_flag_10min(ind1:ind1+length(t_10min)-1);
% flag = 0; % script outputs plot for outside of cold pool times 

k = 4; % level of entrainment end member

%% 4-panel plot
figure;

clf(); hold on

%% Plots in theta-q space %%
for ipanel = 1:2
    ax(ipanel) = subplot(2,2,ipanel, 'align', 'fontsize',18); hold on;
    flag = mod(ipanel,2); % cold pool plot
    if flag
        titlestr = "cold pool";
        ii = (iso_cold_pool_flag==flag);
    else
        titlestr = "environment";
        ff = find(iso_cold_pool_flag==flag);
        ii = ff(1:10:end);
    end

    plot(th_ob(iso_cold_pool_flag==flag),q_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',4)
    hold on;
    scatter(th_ob(iso_cold_pool_flag==flag),q_ob(iso_cold_pool_flag==flag),13,dD_ob(iso_cold_pool_flag==flag),'s','filled')
    colormap(flip(jet(16)))
    hdl=colorbar;
    if ~flag
        ylabel(hdl,'q (g kg^-^1)','FontSize',18,'Rotation',90);
    else
        set(hdl, 'visible','off')
    end
    ylabel(hdl,['\deltaD (',char(8240),')'],'FontSize',18,'Rotation',90);
    scatter(th_ent(k,ii),q_ent(k,ii),10,[0.8500 0.3250 0.0980],'filled')
    scatter(th_surf(ii),q_surf(ii),10,[0.9290 0.6940 0.1250],'filled')
    scatter(th_dd(ii),q_dd(ii),10,[0 0.4470 0.7410],'filled')
    xlim([289 303]); % xlim([288.8 302.2])
    ylim([7 23])
    axis square
    box on
    ylabel('q (g kg^{-1})')
    % xlabel('\theta (K)')
    clim([-81 -65])
    title(titlestr, 'fontweight','normal')
    if ~flag
        text(292,21.5,'surface','FontSize',15)
        text(292,9,'entrained','FontSize',15)
        text(289.3,15.8,'downdraft','FontSize',15)
    end
end
% set(findall(gcf,'-property','Fontsize'),'FontSize',18)
% set(findall(gcf,'-property','TickLength'),'TickLength', [.05 .05])
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1)

%% Plots in dD-theta space
for ipanel = 3:4
    ax(ipanel) = subplot(2,2,ipanel, 'align', 'fontsize',18); hold on;
    flag = mod(ipanel,2); % cold pool plot
    if flag
        titlestr = "cold pool";
        ii = (iso_cold_pool_flag==flag);
    else
        titlestr = "environment";
        ff = find(iso_cold_pool_flag==flag);
        ii = ff(1:10:end);
    end

    plot(th_ob(ii),dD_ob(ii),'sk','MarkerSize',5,'MarkerFaceColor','k')
    hold on;
    scatter(th_ob(ii),dD_ob(ii),13,q_ob(ii),'s','filled')
    colormap(flip(jet(12)))
    clim([12 16])
    hdl=colorbar;
    if ~flag
        ylabel(hdl,'q (g kg^-^1)','FontSize',16,'Rotation',90);
    else
        set(hdl, 'visible','off')
    end
    scatter(th_ent(k,ii),dD_ent(k,ii),15,[0.8500 0.3250 0.0980],'filled')
    scatter(th_surf(ii),dD_surf(ii),15,[0.9290 0.6940 0.1250],'filled')
    scatter(th_dd(ii),dD_dd(k,ii),15,[0 0.4470 0.7410],'filled')
    ylim([-92 -62]);
    xlim([289 303]);
    xlabel('\theta (K)')
    ylabel(['\deltaD (',char(8240),')'])
    % set(findall(gcf,'-property','Fontsize'),'FontSize',14)
    set(findall(gcf,'-property','TickLength'),'TickLength', [.05 .05])
    % set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
    % title(titlestr, 'fontweight','normal')
end

text(ax(1), 290,21, "a", 'fontsize',18)
text(ax(2), 290,21, "c", 'fontsize',18)
text(ax(3), 290,-66, "b", 'fontsize',18)
text(ax(4), 290,-66, "d", 'fontsize',18)

% save in several formats
fmt = ["epsc"; "svg"; "png"; "pdf"];
% for i = 1:length(fmt)
%     saveas(gcf, "con_prop_dD_4panel", fmt(i))
% end

%% Plots in dD-DXS space %%
for k = 4 % 1:size(dD_ent,1)
    figure;
%     plot(dD_ob,de_ob,'sk','MarkerSize',5)
    plot(dD_ob(iso_cold_pool_flag==flag),de_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',5)
    hold on;
%     scatter(dD_ob,de_ob,13,th_ob,'s','filled')
    scatter(dD_ob(iso_cold_pool_flag==flag),de_ob(iso_cold_pool_flag==flag),13,th_ob(iso_cold_pool_flag==flag),'s','filled')
    colormap(jet(15))
    colorbar;
    scatter(mean(dD_ent(k,:),'omitnan'),mean(de_ent(k,:),'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
%     scatter(dD_ent(k,:),de_ent(k,:),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(dD_ent(k,iso_cold_pool_flag==flag),de_ent(k,iso_cold_pool_flag==flag),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(dD_dd(k,iso_cold_pool_flag==flag),'omitnan'),mean(de_dd(k,iso_cold_pool_flag==flag),'omitnan'),50,[0 0.4470 0.7410],'filled')
    scatter(dD_dd(k,:),de_dd(k,:),15,[0 0.4470 0.7410],'filled')
%     scatter(dD_dd(k,iso_cold_pool_flag==flag),de_dd(k,iso_cold_pool_flag==flag),15,[0 0.4470 0.7410],'filled')
    scatter(mean(dD_surf,'omitnan'),mean(de_surf,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
%     scatter(dD_surf,de_surf,15,[0.9290 0.6940 0.1250],'filled')
%     scatter(dD_surf(iso_cold_pool_flag==flag),de_surf(iso_cold_pool_flag==flag),15,[0.9290 0.6940 0.1250],'filled')
    xlim([-84 -56]) 
    ylim([-4 30])
    ylabel(['DXS [',char(8240),']'])
    xlabel(['\deltaD [',char(8240),']'])
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength', [.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
    hold on;text(-66,1,'surface','FontSize',15)
    hold on;text(-75,21.5,'entrained','FontSize',15)
    hold on;text(-60,11.5,'downdraft','FontSize',15)
    caxis([295 300])
    set(gca, 'YDir','reverse')
    set(gca, 'XDir','reverse')
    title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
%     saveas(gcf,['RHB-soundings-DXS-dD-dd-mean-height-range',num2str(k)],'png')
%     saveas(gcf,['RHB-soundings-DXS-dD-dd-mean-height-range',num2str(k)],'fig')
end

%% Plots to showcase sensitivity analysis %%
% 1st try : q-dD space
time_clips = [301:600;601:900;901:1200;1201:1500;1501:1800;1801:2100];
for k = 1:size(time_clips,1)
    figure;
    title(['From ',datestr(t(time_clips(k,1)),'mmm-dd'),' to ',datestr(t(time_clips(k,end)),'mmm-dd')])
    for c = 1:length(h_hi)
        hold on;
        plot(q_ent(c,time_clips(k,:)),dD_dd(c,time_clips(k,:)))
%         plot(q_dd(time_clips(k,:)),q_ent(c,time_clips(k,:)))
    end
    xlabel('q_e_n_t [g/kg]')
%     xlabel('q_d_d [g/kg]')
    ylabel(['\deltaD_d_d [',char(8240),']'])
%     ylabel(['DXS_e_n_t [',char(8240),']'])
    legend('0.8-1km','0.9-1.1km','1-1.2km','1.1-1.3km','1.2-1.4km')
end

figure;
title('Entire timeseries')
for c = 1:length(h_hi)
    hold on;
    plot(q_ent(c,:),de_ent(c,:))
end
xlabel('q_e_n_t [g/kg]')
% xlabel('q_d_d [g/kg]')
ylabel(['DXS_e_n_t [',char(8240),']'])
legend('0.8-1km','0.9-1.1km','1-1.2km','1.1-1.3km','1.2-1.4km')
%%
figure;
plot(t,fdd(4,:))
datetick('x','mmm-dd','keeplimits','keepticks')
