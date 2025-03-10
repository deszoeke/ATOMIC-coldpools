%% Air type: observed
% q and theta @ 400m
zrf = 400; % [in meters]; reference height
zm = 17;   % [in meters]; measurement level
Rd = 287.04;
Cp = 1005.7;
[th_ob,q_ob,t_adj] = height_adj(zrf,zm,Rd,Cp);

% addpath('C:\Users\estef\Documents\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files')
load('2nd_leg_sounding_data_10min_linear_interp.mat')

pos_adj = 999999*ones(size(t));
for l = 1:length(t)
    if l == length(t)
        break
    end
    pos1 = find(t_adj>=t(l)); % rounding up to closest (in time) PSD surface data point!!!
                              % try rounding to nearest (in space) data point
    pos_adj(l) = pos1(1);
    clearvars pos1
end

q_ob = q_ob(pos_adj(1:end-1));
th_ob = th_ob(pos_adj(1:end-1));

t_10min = t;

%% Loading and replacing previous sounding data
% Going from 10 min to 4 hours for the next 3 air types
load('full_214_soundings_Level2_h&p_same_size.mat'); % data at 4hr intervals
% or load ('sounding_data_Level2_h&p_same_size.mat')
% or load ('2hPa_radiosonde_data_ATOMIC.m')

% Using only 2nd leg sounding data
ind = find(t>=t_10min(1));
t = t(ind(1):end);

%% Air type: entrained
% q and theta @ 1km
q_en1km  = q(h==1e3,ind(1):end-1)*1e3; % check units for q_inth, looks like is cg/kg
th_en1km = th(h==1e3,ind(1):end-1);
% % q and theta for layer between 0.8 and 1.2 km
% q_en  = mean(q(h>=800 & h<=1200 ,ind(1):end),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_en = mean(th(h>=800 & h<=1200,ind(1):end),'omitnan');
% % q and theta for layer between 1 and 1.4 km
% q_en  = mean(q(h>=1000 & h<=1400 ,ind(1):end),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_en = mean(th(h>=1000 & h<=1400,ind(1):end),'omitnan');
% % q and theta for layer between 1.1 and 1.4 km
% q_en  = mean(q(h>=1100 & h<=1400 ,ind(1):end),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
% th_en = mean(th(h>=1100 & h<=1400,ind(1):end),'omitnan');
% q and theta for layer between 1.1 and 1.3 km
q_en  = mean(q(h>=1100 & h<=1300 ,ind(1):end-1),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
th_en = mean(th(h>=1100 & h<=1300,ind(1):end-1),'omitnan');

%% Air type: downdraft from within cloud layer 
% q and theta @ minimum theta_w (wet-bub potential temp) within the cloud layer, defined as area/region below the trade inversion line (6g/kg contour) for each sounding
%     % Extracting trade inversion height (mixed layer depth)
%     for k = 1:size(q,2)
%         h6 = double(h(q(:,k)*1000<=6));
%         if h6(1) <= 6000
%             trade_inv(k) = h6(1); % trade inversion
%         else
%             trade_inv(k) = NaN; % trade inversion
%         end
%     end
    
%     for k = 1:size(thw,2)
%         trade_index = find(h==trade_inv(k));
%         th_d(k) = nanmin(thw(1:trade_index,k));
%     end
%     addpath('/Users/estefania/Documents/Research Year 2020-2021/Recovery_PersonalLaptop_09022020/OneDrive/Documents/ATOMIC_Files/RHB raw files/thermo')
%     addpath('C:\Users\estef\Documents\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files\thermo')

%     q_d = qs(1e3*1e2,th_d-273.15)*1e3; % q_d in g/kg    
%     qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
%     p [Pa], T [degrees C], qs [kg/kg]
%     theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
%     th_d = th_d(7:end-1);
%     q_d = q_d(7:end-1);
    
%% Air type: downdraft from mean cloud layer 
% % q and theta @ mean theta_w (wet-bub potential temp) above 1km and below the trade inversion line (6g/kg contour) for each sounding
    % Extracting trade inversion height (mixed layer depth)
    for k = 1:size(q,2)
        h6 = double(h(q(:,k)*1000<=6));
        if h6(1) <= 6000
            trade_inv(k) = h6(1); % trade inversion
        else
            trade_inv(k) = NaN; % trade inversion
        end
        clearvars h6
    end
    thw_index = find(h==1e3);
    for k = 1:size(thw,2)
        if trade_inv(k) >= 1e3
            trade_index = find(h==trade_inv(k));
            th_d(k) = nanmean(thw(thw_index:trade_index,k));
        else
            th_d(k) = NaN; % trade inversion
        end
    end
%   addpath('/Users/estefania/Documents/Research Year 2020-2021/Recovery_PersonalLaptop_09022020/OneDrive/Documents/ATOMIC_Files/RHB raw files/thermo')
    addpath('C:\Users\quinones\Documents\Data\thermo')
    q_d = qs(1e3*1e2,th_d-273.15)*1e3; % q_d in g/kg    
%   qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
%   p [Pa], T [degrees C], qs [kg/kg]
%   theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
    th_d = th_d(ind(1):end-1);
    q_d = q_d(ind(1):end-1);
    
%% Air type: downdraft from 1km
% % q and theta @ 1km
% addpath('/Users/estefania/Documents/Research Year 2020-2021/Recovery_PersonalLaptop_09022020/OneDrive/Documents/ATOMIC_Files/RHB raw files/thermo')
%     theta = th(h==1e3,7:end-1); % potential temp in K
%     % p is already in Pa
%     % q is already in kg/kg
%     th_d = theta_w_eqm(p(h==1e3,7:end-1),theta,q(h==1e3,7:end-1));
%     q_d = qs(1e3*1e2,th_d-273.15)*1e3; % q_d in g/kg
% qs(p,T) is saturation specific humidity based on Wexler's formula for es
% with enhancement factor (see es.m).
% p [Pa], T [degrees C], qs [kg/kg]
% theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008

%% Air type: free troposphere
% q and theta @ minimum theta_w (wet-bub potential temp) pressure level for each sounding
%     addpath('/Users/estefania/Documents/Research Year 2020-2021/Recovery_PersonalLaptop_09022020/OneDrive/Documents/ATOMIC_Files/RHB raw files/thermo')
%     [min_thw,min_idx] = min(theta_w_eqm(p(:,7:end-1),th(:,7:end-1),q(:,7:end-1)),[],1); % min_idx gives the indices
%     th_d = min_thw; % in K
%     q_d = qs(1e3*1e2,th_d-273.15)*1e3; % q_d in g/kg
% qs(p,T) is saturation specific humidity based on Wexler's formula for es
% with enhancement factor (see es.m).
% p [Pa], T [degrees C], qs [kg/kg]
% theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008

%% Air type: surface
% q and theta @ 1m(?)
load ('1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat');
% p = 1000; % mb or hPa
% qv = qs;
% Pi = (p/100000).^((Rd/Cp)*(1-0.28*qv));
% th = th_int(:,p_int==1000); % in K
% Temp = th.*Pi;
% qsat = qs(p*100,temp); % (saturation specific humidity in g/kg)

pos = 999999*ones(size(t));
for l = 1:length(t)
    if l == length(t)
        break
    end
    pos1 = find(t1min>=t(l)); % rounding up to closest (in time) PSD surface data point!!!
                              % try rounding to nearest (in space) data point
    pos(l) = pos1(1);
    clearvars pos1
end
SLP = slp(pos(1:end-1)); % [in hPa]
T = sst(pos(1:end-1)); % [in degrees C]
q_surf = qs(SLP'*100,T')*1e3; % in g/kg; slp must be in Pa and T in degrees C
% q_test = qs_full(pos(7:end-1)); % [in g/kg]
th_surf = (T' + 273.15).*(1e5./(SLP'*1e2)).^(Rd/Cp); % (Potential Temp in degrees K)
q_surf = q_surf';
th_surf = th_surf';

%% Loading iso data
load('iso_data_1min_intervals_FLAGGED_w_runningmean.mat');

% Incorporating isotope data %
pos_i = 9999*ones(size(t_10min));
for l = 1:length(t_10min)
%     if l == 203
%         break
%     end
    pos2 = find(iso_time>=t_10min(l)); % rounding up to closest isotope surface data point!!!
                                                  % try rounding to nearest data point
    pos_i(l) = pos2(1);
end
pos_i(pos_i==9999) = 1;
iso_d18O = d18O(pos_i);
iso_dD = dD(pos_i);
% iso_mr = mr(pos_i);

iso_d18O = iso_d18O(1:end-1);
iso_dD = iso_dD(1:end-1);
% iso_mr = iso_mr(7:end-1);

% iso_d18O(iso_d18O==iso_d18O(1)) = NaN;
% iso_dD(iso_dD==iso_dD(1)) = NaN;
% iso_mr(iso_mr==iso_mr(1)) = NaN;
iso_de = iso_dD - 8*iso_d18O; % deuterium excess

%% Plots %%
figure;
scatter(th_ob,q_ob,15,'k')
hold on;
scatter(th_en,q_en,10,[0.8500 0.3250 0.0980],'filled')
hold on;
scatter(th_d,q_d,10,[0 0.4470 0.7410],'filled')
hold on;
scatter(th_surf,q_surf,10,[0.9290 0.6940 0.1250],'filled')
xlim([288.8 302.2])
ylim([7 23])
ylabel('q [g/kg]')
xlabel('\theta [K]')
set(findall(gcf,'-property','Fontsize'),'FontSize',30)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05 .05])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
box on
hold on;text(294.25,21.5,'surface','FontSize',15)
hold on;text(295,9,'entrained','FontSize',15)
hold on;text(289.5,11,'downdraft','FontSize',15)
% legend('observed')
axis square

%% With iso data
figure;
scatter(th_en,q_en,25,[0.8500 0.3250 0.0980],'filled')
hold on;
scatter(th_ob,q_ob,30,'k')
hold on;
scatter(th_d,q_d,25,[0 0.4470 0.7410],'filled')
hold on;
scatter(th_surf,q_surf,25,[0.9290 0.6940 0.1250],'filled')
hold on;
scatter(th_ob,q_ob,25,iso_de,'filled')
colormap(flip(jet(15)))
han = colorbar;
han.Title.String = "d_e_x_c_e_s_s [per mil]";
xlim([288.8 302.2])
ylim([7 23])
ylabel('q [g/kg]')
xlabel('\theta [K]')
set(findall(gcf,'-property','Fontsize'),'FontSize',30)
box on
legend('entrained','observed','downdraft','surface')

%% Saving data in .mat file %%
time_soundings = t(7:end-1)';
time_PSD_surface_data = t1min(pos(7:end-1)); % only applies to q_surf and th_surf
clearvars -except q_d q_en q_ob q_surf time_soundings time_PSD_surface_data th_d th_en th_ob th_surf iso_mr iso_dD iso_d18O

%% Mix fraction plot %%
% load('conserved_variables_with_iso_data_full_214_soundings.mat')

[fen, fss, fdd] = th_q_to_mixfraction(th_ob,q_ob, th_en,q_en, th_surf,q_surf, th_d,q_d);
% [fen, fss, fdd] = th_q_to_mixfraction(thml,qml, then,qen, thss,qss, thdd,qdd)
% This function diagnoses BL air as a 3-part mixture: fen + fss + fdd = 1
% for each timestamp (each row), 207 in this case.
% fen = entrained air fraction
% fss = fraction of air in equilibrium with the sea surface
% fdd = saturated downdraft air fraction
% The mixing fractions of these 3 end members are inverted algebraically
% from observed thml (th_ob), qml (q_ob) as in de Szoeke (2018).

%% Interpolate sounding end members to the 4 hour PSL met/flux times
time_PSD = time_PSD_surface_data(~isnan(fdd));
fdd = fdd(~isnan(fdd));
fss = fss(~isnan(fss));
fen = fen(~isnan(fen));

fdd_4hr = interp1(time_PSD,fdd,time_PSD_surface_data);
fss_4hr = interp1(time_PSD,fss,time_PSD_surface_data);
fen_4hr = interp1(time_PSD,fen,time_PSD_surface_data);

%% Interpolate sounding end members to the 10 min PSL met/flux times
fdd_10min = interp1(time_PSD,fdd,t1min);
fss_10min = interp1(time_PSD,fss,t1min);
fen_10min = interp1(time_PSD,fen,t1min);

%% Plotting q and theta for entrained air
color = 'g'; %[0.8500 0.3250 0.0980]; %'k'; %[0 0.4470 0.7410]; %

% figure;
subplot(311)
hold on;
plot(t,q_en,'Color',color)
datetick('x','mm/dd','keeplimits','keepticks')
ylabel('q_e_n_t [g/kg]')
legend('@1km','mean [0.8 1.2] km','mean [1.0 1.4] km','mean [1.1 1.4] km','mean [1.1 1.3] km')

subplot(312)
hold on;
plot(t,th_en,'-.','Color',color)
datetick('x','mm/dd','keeplimits','keepticks')
ylabel('\theta_e_n_t [K]')

subplot(313)
hold on;
plot(t,q_en-q_en1km,'Color',color)
plot(t,th_en-th_en1km,'-.','Color',color)
datetick('x','mm/dd','keeplimits','keepticks')
ylabel('diff [mean - 1km]')

%% scatter(th_en,q_en,15,'r','filled')
% hold on;
% scatter(th_ob,q_ob,15,'k','filled')
% hold on;
scatter(fdd,fss,30,'k');
hold on;
scatter(fdd,fss,25,iso_de,'filled');
colormap(flip(jet(15)))
han = colorbar;
han.YLabel.String = "d_e_x_c_e_s_s [per mil]";
hold on;
plot([0 1],[0 0],'-k');
hold on;
plot([0 0],[0 1],'-k');
hold on;
plot([1 0],[0 1],'-k');
% hold on;
% contourf(fdd,fss,fen)
set(gca,'Xdir','reverse')
% scatter(th_surf,q_surf,15,'m','filled')
xlim([-0.2 1])
ylim([-0.2 1])
ylabel('surface fraction ~ q [g/kg]')
xlabel('downdraft fraction ~ \theta [K]')
set(findall(gcf,'-property','Fontsize'),'FontSize',30)
box on
% legend('entrained','observed','downdraft','surface')
title('BL air mixture fraction')

%% Plotting for iso data timeframe only (no 1st leg sounding data)
% BL conserved properties plot
figure;
scatter(th_en(~isnan(iso_de)),q_en(~isnan(iso_de)),25,[0 0.4470 0.7410],'filled')
hold on;
scatter(th_ob(~isnan(iso_de)),q_ob(~isnan(iso_de)),30,'k')
hold on;
scatter( th_d(~isnan(iso_de)), q_d(~isnan(iso_de)),25,[0.9290 0.6940 0.1250],'filled')
hold on;
scatter(th_surf(~isnan(iso_de)),q_surf(~isnan(iso_de)),25,[0.8500 0.3250 0.0980],'filled')
hold on;
scatter(th_ob(~isnan(iso_de)),q_ob(~isnan(iso_de)),25,iso_de(~isnan(iso_de)),'filled')
colormap(flip(jet(15)))
han = colorbar;
han.Title.String = "d_e_x_c_e_s_s [per mil]";
xlim([285 305])
ylim([5 25])
ylabel('specific humidity q [g/kg]')
xlabel('potential temperature \theta [K]')
set(findall(gcf,'-property','Fontsize'),'FontSize',30)
box on
legend('entrained','observed','downdraft','surface')
title('BL conserved properties')

%% BL air mixture fraction plot
figure;
scatter(fdd(~isnan(iso_de)),fss(~isnan(iso_de)),30,'k');
hold on;
scatter(fdd(~isnan(iso_de)),fss(~isnan(iso_de)),25,iso_d18O(~isnan(iso_de)),'filled');
colormap(flip(jet(15)))
han = colorbar;
han.YLabel.String = "d18O [per mil]";
hold on;
plot([0 1],[0 0],'-k');
hold on;
plot([0 0],[0 1],'-k');
hold on;
plot([1 0],[0 1],'-k');
set(gca,'Xdir','reverse')
xlim([-0.2 1])
ylim([-0.2 1])
ylabel('surface fraction ~ q [g/kg]')
xlabel('downdraft fraction ~ \theta [K]')
set(findall(gcf,'-property','Fontsize'),'FontSize',30)
box on
axis square
title('BL air mixture fraction')

%% Area plots for BL air mix fraction %%
% converting time into seconds since 2020 01 01
% time_in_sec = (time_PSD_surface_data - datenum('2020-01-01 00:00:00','yyyy-mm-dd HH:MM:SS'))*3600*24;
% time_PSD_surface_data(isnan(fss)==1) = [];
% fdd(isnan(fdd)==1) = [];
% fen(isnan(fen)==1) = [];
% fss(isnan(fss)==1) = [];
num_days = 9+4/24:4/24:43+12/24;

figure;
% subplot(2,1,1);
% % area(time_soundings,[fen;fss;fdd]'); hold on;
% % t3 = time_soundings;
% % patch([t3(91) t3(91) t3(92) t3(92)],[0 1 1 0],'w')
% % patch([t3(129) t3(129) t3(130) t3(130)],[0 1 1 0],'w')
bar(t1min,[fen_10min;fss_10min;fdd_10min]',1,'stacked')
datetick('x','mmm/dd','keepticks','keeplimits')
xlabel('2020 date')
ylabel('BL air mixture fraction')
legend({'entrained fraction','surface fraction','downdraft fraction'})
set(findall(gcf,'-property','Fontsize'),'FontSize',15)
set(findall(gcf,'Type','axes'),'LineWidth',2)
% xlim([datenum('Jan/09/2020') datenum('Feb/12/2020 1PM')])

% color for area plot ---- RGB Triplet ---- %
% blue                  [0 0.4470 0.7410]       entrained
% mustard               [0.9290 0.6940 0.1250]  downdraft
% orange                [0.8500 0.3250 0.0980]  surface