% Copyright 2025 Simon P. de Szoeke and Estefanía Quiñones 
% Meléndez.
% 
% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation 
% files (the “Software”), to deal in the Software without 
% restriction, including without limitation the rights to use, 
% copy, modify, merge, publish, distribute, sublicense, and/or 
% sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following 
% conditions:
%
% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.
% time subsetting and interpolation not correct

%% Air type: observed
% q and theta @ 400m
zrf = 400; % [in meters]; reference height
zm = 17;   % [in meters]; measurement level
Rd = 287.04;
Cp = 1005.7;

% addpath('C:\Users\estef\Documents\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB raw files')
load('data/2nd_leg_sounding_data_10min_linear_interp.mat') % sounding data interpolated to 10 min
% 10 min time base
t_10min = t; % (1:end-1);

[th_ob_,q_ob_,t_adj] = height_adj(zrf,zm,Rd,Cp); % all PSD 10 min times, length 5040
% interp t_adj -> t_10min
intit = @(x) interp1(t_adj, x, t_10min, 'nearest');
th_ob = intit(th_ob_);
q_ob  = intit( q_ob_);

%% Load 4-h sounding data
% Going from 10 min to 4 hours for the next 3 air types
load('data/full_214_soundings_Level2_h&p_same_size.mat'); % data at 4hr intervals
% or load ('sounding_data_Level2_h&p_same_size.mat')
% or load ('2hPa_radiosonde_data_ATOMIC.m')

% Using only 2nd leg sounding data
ind = find(t >= t_10min(1), 1, 'first');
t_4h = t(ind:end); % length 100

%% Air type: entrained
% q and theta @ 1km
q_en1km  = q( h==1e3, ind:end)*1e3; % check units for q_inth, looks like is cg/kg
th_en1km = th(h==1e3, ind:end);
% vertically average
q_en  = mean(q(h>=1100 & h<=1300 ,ind:end),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
th_en = mean(th(h>=1100 & h<=1300,ind:end),'omitnan');
    
%% Air type: downdraft from mean cloud layer 
% % q and theta @ mean theta_w (wet-bulb potential temp) above 1km and below the trade inversion line (6g/kg contour) for each sounding
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
th_d_ = NaN + zeros(1,size(thw,2));
for k = 1:size(thw,2)
    if trade_inv(k) >= 1e3
        trade_index = find(h==trade_inv(k));
        th_d_(k) = nanmean(thw(thw_index:trade_index,k));
    else
        th_d_(k) = NaN; % trade inversion
    end
end
%   addpath('/Users/estefania/Documents/Research Year 2020-2021/Recovery_PersonalLaptop_09022020/OneDrive/Documents/ATOMIC_Files/RHB raw files/thermo')
%   addpath('C:\Users\quinones\Documents\Data\thermo')

q_d_ = qs(1e3*1e2,th_d_-273.15)*1e3; % q_d in g/kg    
%   qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
%   p [Pa], T [degrees C], qs [kg/kg]
%   theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
th_d = th_d_(ind:end);
q_d = q_d_(ind:end);
% upsample sounding sounding-derived data to 10-min
intit = @(x) interp1(t_4h, x, t_10min, 'linear');
th_d_10 = intit(th_d);
q_d_10  = intit( q_d);
th_en_10 = intit( th_en );
q_en_10  = intit(  q_en );

%% Air type: surface
% q and theta @ 1m(?)
load ('data/1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat');
% interpolate to t_10min
intit = @(x) interp1(t1min, x, t_10min(1:end), 'nearest');
SLP = intit(slp)';
T = intit(sst)';
q_surf = qs(SLP'*100,T')*1e3;
th_surf = (T + 273.15).*(1e5./(SLP*1e2)).^(Rd/Cp);

%% Loading iso data
load('data/iso_data_1min_intervals_FLAGGED_w_runningmean.mat');
intit = @(x) interp1(iso_time, x, t_10min(1:end), 'nearest');
iso_d18O = intit( d18O );
iso_dD   = intit( dD   );
iso_de = iso_dD - 8*iso_d18O; % deuterium excess

%% 
F = load("data/cold_pool_flag_1min_UPDATED_for_iso_times.mat");
cold_pool_flag_10min = logical(interp1(F.t1min,F.cold_pool_flag_1min, t_10min, 'nearest'));

%% Plots %%
figure;

%% color q-theta plot with iso data
clf; subplot(1,1,1); hold on;
scatter(th_en,q_en,25,[0.8500 0.3250 0.0980],'filled')
scatter(th_d,q_d,25,[0 0.4470 0.7410],'filled')
scatter(th_surf,q_surf,25,[0.9290 0.6940 0.1250],'filled')
% scatter(th_ob,q_ob,18,'k')
scatter(th_ob,q_ob,15,iso_de,'filled')
plot(th_ob( cold_pool_flag_10min),q_ob( cold_pool_flag_10min),'ko')
clim([9 16]); b2rcolormap(14);
han = colorbar;
han.Label.String = ['deuterium excess (' char(8240) ')'];
xlim([288.8 302.2])
ylim([7 23])
ylabel('q [g/kg]')
xlabel('\theta [K]')
set(findall(gcf,'-property','Fontsize'),'FontSize',18)
box on; axis square
% legend('entrained','observed','downdraft','surface')

%% Saving data in .mat file %%
time_soundings = t_4h(7:end-1)';
% time_PSD_surface_data = t1min(pos(7:end-1)); % only applies to q_surf and th_surf
clearvars -except q_d q_en q_ob q_surf time_soundings time_PSD_surface_data th_d th_en th_ob th_surf iso_mr iso_dD iso_d18O

%% Mix fraction plot %%
% load('conserved_variables_with_iso_data_full_214_soundings.mat')

[fen, fss, fdd] = th_q_to_mixfraction(th_ob(:),q_ob(:), th_en_10(:),q_en_10(:), th_surf(:),q_surf(:), th_d_10(:),q_d_10(:));
% This function diagnoses BL air as a 3-part mixture: fen + fss + fdd = 1
% for each timestamp (each row), 207 in this case.
% fen = entrained air fraction
% fss = fraction of air in equilibrium with the sea surface
% fdd = saturated downdraft air fraction
% The mixing fractions of these 3 end members are inverted algebraically
% from observed thml (th_ob), qml (q_ob) as in de Szoeke (2018).

%% sounding derived end members already interpolated to 10-min
fdd_10min = fdd;
fss_10min = fss;
fen_10min = fen;

%% Plot q and theta for entrained air
color = 'b'; %[0.8500 0.3250 0.0980]; %'k'; %[0 0.4470 0.7410]; %

% figure;
subplot(311)
hold on;
plot(t_4h,q_en,'Color',color)
datetick('x','mm/dd','keeplimits','keepticks')
ylabel('q_e_n_t [g/kg]')
legend('@1km','mean [0.8 1.2] km','mean [1.0 1.4] km','mean [1.1 1.4] km','mean [1.1 1.3] km')

subplot(312)
hold on;
plot(t_4h,th_en,'-.','Color',color)
datetick('x','mm/dd','keeplimits','keepticks')
ylabel('\theta_e_n_t [K]')

subplot(313)
hold on;
plot(t_4h,q_en-q_en1km,'Color',color)
plot(t_4h,th_en-th_en1km,'-.','Color',color)
datetick('x','mm/dd','keeplimits','keepticks')
ylabel('diff [mean - 1km]')

%% Plot Fig 9
clrs = colororder();
figure

% scatter(th_en,q_en,15,'r','filled')
clf;
subplot(1,2,2, 'align','fontsize',18)
scatter(fdd,fss,20,'k');
hold on;
scatter(fdd,fss,15,iso_dD,'filled');
colormap(flip(jet(15)))
han = colorbar;
han.YLabel.String = ['\deltaD (', char(8240), ')'];
hold on;
plot([0 1],[0 0],'-k');
hold on;
plot([0 0],[0 1],'-k');
hold on;
plot([1 0],[0 1],'-k');
plot(1,0, 'ko', 'markerfacecolor',clrs(1,:), 'markersize',10)
plot(0,0, 'ko', 'markerfacecolor',clrs(2,:), 'markersize',10)
plot(0,1, 'ko', 'markerfacecolor',clrs(3,:), 'markersize',10)
set(gca,'Xdir','reverse')
xlim([-0.2 1.1])
ylim([-0.2 1.1])
ylabel('surface fraction')
xlabel('downdraft fraction')
set(findall(gcf,'-property','Fontsize'),'FontSize',18)
box on; axis square;
% legend('entrained','observed','downdraft','surface')
text(0.95, 0.95, 'b', 'FontSize',18)
title('\theta-q mixture fraction', 'fontweight','normal')

subplot(1,2,1, 'align','fontsize',18); hold on;
plot(th_ob(~cold_pool_flag_10min),q_ob(~cold_pool_flag_10min),'.', 'color',0.7+[0 0 0])
plot(th_ob( cold_pool_flag_10min),q_ob( cold_pool_flag_10min),'ko')
scatter(th_en,q_en,10,[0.8500 0.3250 0.0980],'filled')
scatter(th_d,q_d,10,[0 0.4470 0.7410],'filled')
scatter(th_surf,q_surf,10,[0.9290 0.6940 0.1250],'filled')
xlim([288.8 302.2])
ylim([7 23])
ylabel('specific humididty q (g/kg)')
xlabel('potential temperature \theta (K)')
hdl = colorbar(); hdl.Visible = "off";
set(findall(gcf,'-property','Fontsize'),'FontSize',18)
set(findall(gcf,'-property','TickLength'),'TickLength',[.05 .05])
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
text(293.2,21.5,'surface','FontSize',18)
text(294,8.75,'entrained','FontSize',18)
text(289.5,11,'downdraft','FontSize',18)
text(290, 21, 'a', 'FontSize',18)
% legend('observed')
box on; axis square

% print Fig 9
orient landscape
fmt = ["epsc"; "svg"; "png"; "pdf"];
% for i = 1:length(fmt)
%     saveas(gcf, "con_prop_thq_ab", fmt(i))
% end

% END working part of script
%{
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
%}