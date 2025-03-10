% Adjusting the code extracted from the conserved_properties_for_isotopes.m
% file to obtain results in 1-min resolution
% Created on: Jan 25 2023
% Last modified: Jan 25 2023
%% Air type: observed %%
% q and theta @ 400m
zrf = 400; % [in meters]; reference height
zm = 17;   % [in meters]; measurement level
Rd = 287.04;
Cp = 1005.7;
[th_ob,q_ob,t_adj] = height_adj_1min(zrf,zm,Rd,Cp);
load '2nd_leg_sounding_data_10min_linear_interp.mat' t
pos1 = find(t_adj>=t(1)); % rounding up to closest (in time) PSD surface data point!!!
pos2 = find(t_adj>=t(end));
pos1 = pos1(1); % start of 2nd leg
pos2 = pos2(1); %   end of 2nd leg
q_ob = q_ob(pos1:pos2);   % keeping 2nd leg data only
th_ob = th_ob(pos1:pos2); % keeping 2nd leg data only
time  = t_adj(pos1:pos2); % 1-min resolution
clearvars t t_adj pos1 pos2
% Loading iso data %
load('iso_data_1min_intervals_FLAGGED_w_runningmean.mat'); % RHB Picarro data
ind0 = find(iso_time>=time(1));
ind0 = ind0(1);
% Incorporating isotope data %
pos_i = ind0:1:ind0+length(time)-1;
d18O_ob = d18O(pos_i);
dD_ob   = dD(pos_i);
de_ob   = dD_ob - (8*d18O_ob); % deuterium excess
clearvars dD d18O DXS iso_time pos_i ind0 factor
%% Air type: surface [using properties in equilibrium w/ surface ocean] %%
% q and theta @ 1m
load '1min_res_PSD_surface_variables_FLAGGED_w_runningmean.mat' sst slp Ta rh t1min
ind0 = find(t1min==time(1));
% Incorporating isotope data %
pos = ind0:1:ind0+length(time)-1;
TA = Ta(pos(1:end))+273.15; % in degrees Kelvin
SLP = slp(pos(1:end)); % [in hPa]
SST = sst(pos(1:end)); % [in degrees C]
RH  =  rh(pos(1:end));  % in percent
[q_surf,th_surf] = air_types(0,0,0,'surface',SLP,SST);
clearvars rh slp sst t1min Ta
%% Air type: surface [using Craig-Gordon] %%
% Incorporating isotope data %
h_prime = 0.875; % relative humidity at which we are modeling the iso ratios
filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
Tskin = ncread(filename,'tskin'); % liquid temperature at [???]m [in degrees C]
Tskin = Tskin(pos);
% [dD_surf,d18O_surf] = air_types_iso(dD_ob,d18O_ob,'surface',T_skin',RH'/100,h_prime);
[dD_surf,d18O_surf] = air_types_iso(dD_ob,d18O_ob,'surface',SST',RH'/100,h_prime);
% load air_type_surface_CG.mat
% d18O_surf = delta_aO;
% dD_surf = delta_aD;
% de_surf = DXS_a;
de_surf = dD_surf - (8*d18O_surf); % deuterium excess
%% Air type: entrained %%
% q and theta @ 1km
load '2nd_leg_sounding_data_10min_linear_interp.mat' q th h t
h_low = 1100;
h_hi  = 1300;
 q_rf = mean( q(h>=h_low & h<=h_hi,:),'omitnan')*1e3; % check units for q_inth, looks like is cg/kg
th_rf = mean(th(h>=h_low & h<=h_hi,:),'omitnan');
% Adjusting time resolution from 10-min to 1-min
th_ent = interp1(t,th_rf,time,'nearest');
q_ent  = interp1(t, q_rf,time,'nearest');
% Incorporating isotope data %
[dD_ent,d18O_ent] = iso_estimates_from_centroid_approach(q_ob,dD_ob,d18O_ob,q_surf,dD_surf,d18O_surf,q_ent);
de_ent = dD_ent - (8*d18O_ent); % deuterium excess
clearvars th_rf q_rf
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
load '2nd_leg_sounding_data_10min_linear_interp.mat' thw
for k = 1:size(thw,2)
    if trade_inv(k) >= 1e3
        trade_index = find(h==trade_inv(k));
        th_rf(k) = nanmean(thw(thw_index:trade_index,k));
    else
        th_rf(k) = NaN; % trade inversion
    end
end
addpath('C:\Users\quinones\Documents\Data\thermo')
q_rf = qs(1e3*1e2,th_rf-273.15)*1e3; % q_d in g/kg    
    % qs(p,T) is saturation specific humidity based on Wexler's formula for es with enhancement factor (see es.m).
    % p [Pa], T [degrees C], qs [kg/kg]
    % theta_w_eqm(p[Pa],theta[K],qv[kg/kg]) from eqn 3.8 of Davies-Jones 2008
% Adjusting time resolution from 10-min to 1-min
th_dd = interp1(t,th_rf,time,'nearest');
 q_dd = interp1(t, q_rf,time,'nearest');
 clearvars th_rf q_rf q th thw t h trade_inv
%% Isotope data for downdraft end member %%
% Isotope concentrations calculated based on the mixing fractions obtained from soundings
[fen, fss, fdd] = th_q_to_mixfraction(th_ob,q_ob, th_ent,q_ent, th_surf,q_surf, th_dd,q_dd);
% dD_dd = (q_ob.*dD_ob' - q_ent.*fen.*dD_ent - q_surf.*fss.*dD_surf')./(q_dd.*fdd);
% d18O_dd = (q_ob.*d18O_ob' - q_ent.*fen.*d18O_ent - q_surf.*fss.*d18O_surf')./(q_dd.*fdd);
% de_dd = dD_dd - (8*d18O_dd); % deuterium excess

fdd(fdd<0)=NaN;
fss(fss<0)=NaN;
fen(fen<0)=NaN;

% This model for downdraft isotopic composition has been replaced by hydrometeor evaporation model
% to avoid the use of mixing fractions to compute dD of downdraft source
% dD_dd = (q_ob.*dD_ob' - q_ent.*fen.*dD_ent - q_surf.*fss.*dD_surf')./(q_dd.*fdd);
% d18O_dd = (q_ob.*d18O_ob' - q_ent.*fen.*d18O_ent - q_surf.*fss.*d18O_surf')./(q_dd.*fdd);
% de_dd = dD_dd - (8*d18O_dd); % deuterium excess
2+2;
%% Cold pool times
load 'cold_pool_flag_1min.mat'
iso_cold_pool_flag = cold_pool_flag_1min(pos);
load 'recovery&peak_flags_1min_full_timeseries.mat'
recovery_flag = recovery_flag_1min(pos);
peak_flag = peak_flag_1min(pos);
flag = 0; % script outputs plot for outside of cold pool times 
flag = 1; % script outputs cold pool plot

% for plotting purposes
t = t1min(pos);
figure;
plot(t(recovery_flag==1),TA(recovery_flag==1)-273.15,'.b')
hold on;
plot(t(recovery_flag==0),TA(recovery_flag==0)-273.15,'.k')
plot(t(peak_flag==1),TA(peak_flag==1)-273.15,'.r')
datetick('x','mm/dd','keeplimits','keepticks')
title('cold pool recovery flag'); ylim([22 28])
%% Slope => partial derivatives with respect to air temperature %%
sfdd   = diff(fdd)./diff(TA);
dsfdd  = diff(dfdd)./diff(dTA);
sdD_ob = diff(dD_ob')./diff(TA);
dsdD_ob = diff(ddD_ob')./diff(dTA);

hold on;plot(time(2:end),sdD_ob,'*k')

sfdd   = diff(fdd)./diff(th_ob);
dsfdd  = diff(dfdd)./diff(dth_ob);
sdD_ob = diff(dD_ob')./diff(th_ob);
dsdD_ob = diff(ddD_ob')./diff(dth_ob);

%% Correlation plots %%
% dTA = detrend(TA,'omitnan');
% dfdd = detrend(fdd,'omitnan');
% Detrend by differencing
% dTA = diff(TA);
% dfdd = diff(fdd);
% Detrend by linear regression
dTA  = detrend_hours_linear(TA,time,6);
dfdd = detrend_hours_linear(fdd,time,6);

[RHO,PVAL] = corr(dTA',-1*dfdd','rows','complete'); % RHO = 0.7036; diff=> 0.2080; lin=> 0.8894(6hrs); 0.7654(24hrs)
% The relationship between two variables is generally considered strong 
% when their r value is larger than 0.7. The correlation r measures the 
% strength of the linear relationship between two quantitative variables.
[RHO,PVAL] = corr(TA',-1*fdd','rows','complete'); % RHO = 0.5559;

figure; hold on;
plot(time,TA,'.')
yyaxis right
plot(time,-1*fdd,'.')
datetick('x','mm/dd','keeplimits','keepticks')
legend('T_a_i_r','-f_d_d')

figure; hold on;
% plot(time,dTA)
% plot(time(2:end),dTA)
plot(time(1:length(dTA)),dTA)
yyaxis right
% plot(time,-1*dfdd)
% plot(time(2:end),-1*dfdd)
plot(time(1:length(dfdd)),-1*dfdd)
datetick('x','mm/dd','keeplimits','keepticks')
legend('detrend(T_a_i_r)','detrend(-f_d_d)')

% dth_ob = detrend(th_ob,'omitnan');
% dth_ob = diff(th_ob);
dth_ob  = detrend_hours_linear(th_ob,time,6);
[RHO,PVAL] = corr(th_ob',-1*fdd','rows','complete');   % RHO = 0.4843;
[RHO,PVAL] = corr(dth_ob',-1*dfdd','rows','complete'); % RHO = 0.6776; diff=> 0.9607; lin=> 0.9498(6hrs); 0.5914(24hrs)

figure; hold on;
plot(time,th_ob)
yyaxis right
plot(time,-1*fdd)
datetick('x','mm/dd','keeplimits','keepticks')
legend('\theta_o_b','-f_d_d')

figure; hold on;
% plot(time,dth_ob)
% plot(time(2:end),dth_ob)
plot(time(1:length(dth_ob)),dth_ob)
yyaxis right
% plot(time,-1*dfdd)
% plot(time(2:end),-1*dfdd)
plot(time(1:length(dfdd)),-1*dfdd)
datetick('x','mm/dd','keeplimits','keepticks')
legend('detrend(\theta_o_b)','detrend(-f_d_d)')

ddD_ob = detrend_hours_linear(dD_ob',time,6);
% [RHO,PVAL] = corr(TA',-1*dD_ob,'rows','complete'); % RHO = 0.0062
% [RHO,PVAL] = corr(dTA',ddD_ob,'rows','complete');  % RHO = 0.0702
% 
% figure; hold on;
% plot(time,TA)
% yyaxis right
% plot(time,-1*dD_ob)
% datetick('x','mm/dd','keeplimits','keepticks')
% legend('T_a_i_r','-dD_o_b')
% 
% figure; hold on;
% plot(time,dTA)
% yyaxis right
% plot(time,ddD_ob)
% datetick('x','mm/dd','keeplimits','keepticks')
% legend('detrend(T_a_i_r)','detrend(dD_o_b)')
% 
% [RHO,PVAL] = corr(th_ob',-1*dD_ob,'rows','complete'); % RHO = 0.0599
% [RHO,PVAL] = corr(dth_ob',ddD_ob,'rows','complete');  % RHO = 0.0477
% 
% figure; hold on;
% plot(time,th_ob)
% yyaxis right
% plot(time,-1*dD_ob)
% datetick('x','mm/dd','keeplimits','keepticks')
% legend('\theta_o_b','-dD_o_b')
% 
% figure; hold on;
% plot(time,dth_ob)
% yyaxis right
% plot(time,ddD_ob)
% datetick('x','mm/dd','keeplimits','keepticks')
% legend('detrend(\theta_o_b)','detrend(dD_o_b)')
% 
% dfss = detrend(fss,'omitnan');
% [RHO,PVAL] = corr(th_ob',fss','rows','complete');   % RHO = 0.2477
% [RHO,PVAL] = corr(dth_ob',dfss','rows','complete'); % RHO = 0.1789
% [RHO,PVAL] = corr(TA',fss','rows','complete');      % RHO = 0.2067
% [RHO,PVAL] = corr(dTA',dfss','rows','complete');    % RHO = 0.1566
% 
% dfen = detrend(fen,'omitnan');
% [RHO,PVAL] = corr(th_ob',fen','rows','complete');   % RHO = 0.1029
% [RHO,PVAL] = corr(dth_ob',dfen','rows','complete'); % RHO = 0.2632
% [RHO,PVAL] = corr(TA',fen','rows','complete');      % RHO = 0.1681
% [RHO,PVAL] = corr(dTA',dfen','rows','complete');    % RHO = 0.2901
% 
% ddD_dd = detrend(dD_dd,'omitnan');
% [RHO,PVAL] = corr(th_ob',dD_dd','rows','complete');   % RHO = 0.0065
% [RHO,PVAL] = corr(dth_ob',ddD_dd','rows','complete'); % RHO = 0.0413
% [RHO,PVAL] = corr(TA',dD_dd','rows','complete');      % RHO = 0.0257
% [RHO,PVAL] = corr(dTA',ddD_dd','rows','complete');    % RHO = 0.0507
% 
% plot(time(dD_dd<100),dD_dd(dD_dd<100))
% 
[RHO,PVAL] = corr(dD_ob,fdd','rows','complete');       % RHO = 0.0179
[RHO,PVAL] = corr(ddD_ob',dfdd','rows','complete');     % RHO = 0.1619
[RHO,PVAL] = corr(dD_dd',-1*fdd','rows','complete');   % RHO = 0.0567
[RHO,PVAL] = corr(ddD_dd',-1*dfdd','rows','complete'); % RHO = 0.0224

%% Cumulative density function (distribution)
Y = dD_dd(dD_dd<-64 & dD_dd>-85);%(iso_cold_pool_flag==flag);
X = Y(~isnan(Y)); 
figure;
plot(sort(X),(1:length(X))./length(X),'-')
hold on;
% plot([0 0.7],[.5 .5])
% text(0.67,0.515,'mean')
% plot([0 0.7],[.95 .95])
% text(0.64,0.965,'95th percentile')
% plot([0 0.7],[.87 .87])
% text(0.65,0.87,'87th percentile')
% plot([.285 .285],[0 1])
% plot([.26 .26],[0 1])
% plot([.198 .198],[0 1])
% plot([.406 .406],[0 1])
Y7 = dD_dd(fdd>=0.285);
X7 = Y7(~isnan(Y7)); 
plot(sort(X7),(1:length(X7))./length(X7),'-')
YY = dD_dd(iso_cold_pool_flag==flag);
XX = YY(~isnan(YY));
plot(sort(XX),(1:length(XX))./length(XX),'-')
Y3 = dD_dd(iso_cold_pool_flag==flag & fdd>=0.285);
X3 = Y3(~isnan(Y3)); 
plot(sort(X3),(1:length(X3))./length(X3),'-')
xlabel('dD_d_d')
ylabel('cummulative density function')
legend('dd estimate','fdd>=0.285','cold pools','cold pools & fdd>=0.285')

Y4 = dD_ob;
X4 = Y4(~isnan(Y4));
figure;
plot(sort(X4),(1:length(X4))./length(X4),'--')
hold on;
Y8 = dD_ob(fdd>=0.285);
X8 = Y8(~isnan(Y8));
plot(sort(X8),(1:length(X8))./length(X8),'--')
Y5 = dD_ob(iso_cold_pool_flag==flag);
X5 = Y5(~isnan(Y5)); 
plot(sort(X5),(1:length(X5))./length(X5),'--')
Y6 = dD_ob(iso_cold_pool_flag==flag & fdd>=0.285);
X6 = Y6(~isnan(Y6)); 
plot(sort(X6),(1:length(X6))./length(X6),'--')
xlabel('dD_o_b')
ylabel('cummulative density function')
legend('analyzer','fdd>=0.285','cold pools','cold pools & fdd>=0.285')
2+2;

figure; hold on;
plot(sort(X4),cumsum(1:length(X4))./length(X4))
plot(sort(X8),cumsum(1:length(X8))./length(X8))
plot(sort(X5),cumsum(1:length(X5))./length(X5))
plot(sort(X6),cumsum(1:length(X6))./length(X6))

figure; hold on;
plot(sort(X),cumsum(1:length(X))./length(X))
plot(sort(X7),cumsum(1:length(X7))./length(X7))
plot(sort(XX),cumsum(1:length(XX))./length(XX))
plot(sort(X3),cumsum(1:length(X3))./length(X3))

%% Plots in dD-theta space %%
% figure;
% plot(th_ob,dD_ob,'sk','MarkerSize',5,'MarkerFaceColor','k')
    figure;
%     plot(th_ob(iso_cold_pool_flag==flag & fdd>=0.285),dD_ob(iso_cold_pool_flag==flag & fdd>=0.285),'sk','MarkerSize',5,'MarkerFaceColor','k')
%     plot(th_ob(iso_cold_pool_flag==flag),dD_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',5,'MarkerFaceColor','k')
%     plot(th_ob(iso_cold_pool_flag==flag),d18O_ob(iso_cold_pool_flag==flag),'sk','MarkerSize',5,'MarkerFaceColor','k')
%     plot(th_ob,dD_ob,'sk','MarkerSize',5,'MarkerFaceColor','k')
%     plot(th_ob(iso_cold_pool_flag==flag & recovery_flag==0),dD_ob(iso_cold_pool_flag==flag & recovery_flag==0),'sk','MarkerSize',5,'MarkerFaceColor','k')
%       plot(th_ob(recovery_flag==1),dD_ob(recovery_flag==1),'sk','MarkerSize',5,'MarkerFaceColor','k')
      plot(th_ob(peak_flag==1),dD_ob(peak_flag==1),'sk','MarkerSize',5,'MarkerFaceColor','k')
    hold on;
%     scatter(th_ob(iso_cold_pool_flag==flag & fdd>=0.285),dD_ob(iso_cold_pool_flag==flag & fdd>=0.285),13,q_ob(iso_cold_pool_flag==flag & fdd>=0.285),'s','filled')
%     scatter(th_ob(iso_cold_pool_flag==flag),dD_ob(iso_cold_pool_flag==flag),13,q_ob(iso_cold_pool_flag==flag),'s','filled')
%     scatter(th_ob(iso_cold_pool_flag==flag & recovery_flag==0),dD_ob(iso_cold_pool_flag==flag & recovery_flag==0),13,q_ob(iso_cold_pool_flag==flag & recovery_flag==0),'s','filled')
%     scatter(th_ob(recovery_flag==1),dD_ob(recovery_flag==1),13,q_ob(recovery_flag==1),'s','filled')
    scatter(th_ob(peak_flag==1),dD_ob(peak_flag==1),13,q_ob(peak_flag==1),'s','filled')
%     scatter(th_ob(iso_cold_pool_flag==flag),d18O_ob(iso_cold_pool_flag==flag),13,q_ob(iso_cold_pool_flag==flag),'s','filled')
%     scatter(th_ob,dD_ob,13,q_ob,'s','filled')
    colormap(flip(jet(12)))
    caxis([11 17])
    hdl=colorbar;
    ylabel(hdl,'q [g kg^-^1])','FontSize',16,'Rotation',90);
%     scatter(th_ent ,dD_ent ,15,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(th_ent,'omitnan'),mean(d18O_ent,'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
%     scatter(mean(th_ent(iso_cold_pool_flag==flag),'omitnan'),mean(d18O_ent(iso_cold_pool_flag==flag),'omitnan'),50,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(iso_cold_pool_flag==flag & fdd>=0.285),dD_ent(iso_cold_pool_flag==flag & fdd>=0.285),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(iso_cold_pool_flag==flag),dD_ent(iso_cold_pool_flag==flag),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(iso_cold_pool_flag==flag & recovery_flag==0),dD_ent(iso_cold_pool_flag==flag & recovery_flag==0),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(recovery_flag==1),dD_ent(recovery_flag==1),15,[0.8500 0.3250 0.0980],'filled')
    scatter(th_ent(peak_flag==1),dD_ent(peak_flag==1),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_ent(iso_cold_pool_flag==flag),d18O_ent(iso_cold_pool_flag==flag),15,[0.8500 0.3250 0.0980],'filled')
%     scatter(th_surf ,dD_surf ,15,[0.9290 0.6940 0.1250],'filled')
%     scatter(mean(th_surf,'omitnan'),mean(dD_surf,'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
%     scatter(mean(th_surf(iso_cold_pool_flag==flag),'omitnan'),mean(dD_surf(iso_cold_pool_flag==flag),'omitnan'),50,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(iso_cold_pool_flag==flag),dD_surf(iso_cold_pool_flag==flag),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(iso_cold_pool_flag==flag & recovery_flag==0),dD_surf(iso_cold_pool_flag==flag & recovery_flag==0),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(recovery_flag==1),dD_surf(recovery_flag==1),15,[0.9290 0.6940 0.1250],'filled')
    scatter(th_surf(peak_flag==1),dD_surf(peak_flag==1),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(iso_cold_pool_flag==flag),d18O_surf(iso_cold_pool_flag==flag),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_surf(iso_cold_pool_flag==flag & fdd>=0.285),dD_surf(iso_cold_pool_flag==flag & fdd>=0.285),15,[0.9290 0.6940 0.1250],'filled')
%     scatter(th_dd ,dD_dd ,15,[0 0.4470 0.7410],'filled')
%     scatter(th_dd(iso_cold_pool_flag==flag & fdd>=0.285),dD_dd(iso_cold_pool_flag==flag & fdd>=0.285),15,[0 0.4470 0.7410],'filled')
%     scatter(th_dd(iso_cold_pool_flag==flag),dD_dd(iso_cold_pool_flag==flag),15,[0 0.4470 0.7410],'filled')
%     scatter(th_dd(iso_cold_pool_flag==flag & recovery_flag==0),dD_dd(iso_cold_pool_flag==flag  & recovery_flag==0),15,[0 0.4470 0.7410],'filled')
%     scatter(th_dd(recovery_flag==1),dD_dd(recovery_flag==1),15,[0 0.4470 0.7410],'filled')
    scatter(th_dd(peak_flag==1),dD_dd(peak_flag==1),15,[0 0.4470 0.7410],'filled')
%     scatter(th_dd(iso_cold_pool_flag==flag),d18O_dd(iso_cold_pool_flag==flag),15,[0 0.4470 0.7410],'filled')
%     scatter(mean(th_dd(iso_cold_pool_flag==flag),'omitnan'),mean(dD_dd(iso_cold_pool_flag==flag),'omitnan'),50,[0 0.4470 0.7410],'filled')
%     ylim([-96 -61]); % xlim([-84 -56])
    xlim([290 302.5]);
    ylim([-90 -60]);
%     ylabel(['DXS [',char(8240),']'])
    xlabel('\theta [K]')
    ylabel(['\deltaD [',char(8240),']'])
    set(findall(gcf,'-property','Fontsize'),'FontSize',30)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.05 .05])
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
    box on
    axis square
%     hold on;text(-66,1,'surface','FontSize',15)
%     hold on;text(-75,21.5,'entrained','FontSize',15)
%     hold on;text(-60,11.5,'downdraft','FontSize',15)
%     caxis([297 300])
%     set(gca, 'YDir','reverse')
%     set(gca, 'XDir','reverse')
%     title(['Height range [',num2str(h_low(k)/1000),':',num2str(h_hi(k)/1000),'] km'],'FontSize',17)
%     text(-58,0,['RMSE = ',num2str(rmse(k)),''],'FontSize',15)
%     saveas(gcf,['RHB-soundings-cps-thvsdD-height-range',num2str(k)],'png')
%     saveas(gcf,['RHB-soundings-cps-thvsdD-height-range',num2str(k)],'fig')
% hold on;
% scatter(mean(dD_surf,'omitnan'),mean(th_surf,'omitnan'),50,'r','filled')
% legend('RB SBL obs','0.8-1km','0.9-1.1km','1-1.2km','1.1-1.3km','1.2-1.4km','surface','FontSize',15)