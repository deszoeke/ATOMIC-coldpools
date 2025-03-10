%% Constellation plots for cold pool events with iso data %%
% Last updated on 01/25/2022
% By E. Quinones Melendez (quinones@oregonstate.edu)

load('2nd_leg_cold_pool&conserved_variables_10min.mat')
clearvars -except tcold time qsurf_cold

% addpath('C:\Users\estef\Documents\Research Year 2020-2021\Recovery_PersonalLaptop_09022020\OneDrive\Documents\ATOMIC_Files\RHB processed files')
filename = 'EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
timef = ncread(filename,'time'); % time of the flux data
timef = timef(2575:4715)'; % 2nd leg data only
timef = timef/3600/24 + datenum('20200101','yyyymmdd');

qa = ncread(filename,'qair'); % air specific humidity from PSL RH at height ztq [g/kg]
q_m = qa(2575:4715)'; % 2nd leg data only
Ta = ncread(filename,'tair'); % air temperature at 17m [in degrees C]

%% q and Dd for cold pool events %%
% Loading iso data
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
timep = ncread(filename,'time'); % times for each 1 min average;'seconds since 2020-01-01 00:00:00'
dD_avg_1min = ncread(filename,'dD'); % in permil ('1e-3');'hydrogen isotope ratio in water vapor (units permil)';'2H/1H ratio in water vapor expressed relative to Vienna Standard Mean Ocean Water'
d18O_avg_1min = ncread(filename,'d18O'); % in permil ('1e-3');'oxygen isotope ratio in water vapor (units permil)';'18O/16O ratio in water vapor expressed relative to Vienna Standard Mean Ocean Water'
mr = ncread(filename,'mmr'); % in g/kg; 'water vapor mass (DRY) mixing ratio'; humidity_mixing_ratio;'ratio of the mass of water vapor to the mass of dry air'
vmr = ncread(filename,'vmr'); % in ppmv (parts per mass volume);'water vapor volume (WET) mixing ratio'; mole_fraction_of_water_vapor_in_air; 'ratio of moles water vapor to moles moist air'
ship_flag = ncread(filename,'ship_flag'); %'measurement quality concern flag based on the ship contamination flag';[0=good,1=poor,2=in port]
inlet_flag = ncread(filename,'inlet_flag'); %'measurement quality concern flag based on the inlet blower direction';[0=normal,1=reversed]
q_iso = ncread(filename,'q'); % air specific humidity [g/kg];'ratio of the mass of water vapor to the mass of moist air'
% q_iso_test = (mr/1000)./(1+(mr/1000))*1e3; % air specific humidity [g/kg] from iso data

% Eliminating NaN's from isodata %
timep(isnan(dD_avg_1min)==1)=[];
mr(isnan(dD_avg_1min)==1)=[];
vmr(isnan(dD_avg_1min)==1)=[];
q_iso(isnan(dD_avg_1min)==1)=[];
ship_flag(isnan(dD_avg_1min)==1)=[];
inlet_flag(isnan(dD_avg_1min)==1)=[];
dD_avg_1min(isnan(dD_avg_1min)==1)=[];
d18O(isnan(d18O_avg_1min)==1)=[];

% Incorporating isotope data %
timep = timep/3600/24 + datenum('20200101','yyyymmdd');
d_ex_1min = dD_avg_1min - 8*d18O_avg_1min; % deuterium excess
% dD_avg_1min(d_ex_1min<0) = NaN;
% q_iso(d_ex_1min<0) = NaN;
% d_ex_1min(d_ex_1min<0) = NaN;

%% Identifying cold pool times (indexes) in isotope timeseries
load('2nd_leg_cold_pool_times.mat')
% Restricting cold pool timeseries to the first 62 min
for k = [2,5,21]
     t_end2(k) = t_max2(k)+(1/24/60)*62;
end

for k = 1:22
    max_ind(k) = find(round(timep,9)==round(t_max2(k),9));
    min_ind(k) = find(round(timep,9)==round(t_min2(k),9));
    end_ind(k) = find(round(timep,9)==round(t_end2(k),9));
end

%% Including q-dD points for surface and entrainment end members
load('10min_res_PSD_wind_surface_variables.mat')
c = 0;
for k = 1:22
    for c = 1:10
    if ~isempty(find((time_full)==tcold(k,c)))
        start_t_ind(k,1) = find((time_full)==tcold(k,c));
    end
    end
end
start_t_ind(19,1) = 4455;
start_t_ind(21,1) = 4689;

for k = 1:22
    qsurf_cold_test(k,:) = qs_full(start_t_ind(k,1)-1:start_t_ind(k,1)+5);
    tcold_10min(k,:)= time_full(start_t_ind(k,1)-1:start_t_ind(k,1)+5);
end

% From 10-min data for 2nd leg only!!!
for k = 1:22
    start_t_ids(k) = find(time==tcold_10min(k,1));
end

% load('k_fractionation_coeff_at_PSD_surface_data_times.mat')
% load('k_H2Ofractionation_coeff_at_iso_data_times.mat')
load('RHS_Eq9_MJ79_variablesUPDATED.mat')
load('conserved_variables_10min_intervals.mat')
% Kinetic effect correction => new fractionation factor: 1/(alpha_e*alpha_k)
% Valid for an undersaturated open system
  % Diffusivities (D) adopted from Merlivat 1978:
    % D_H18O/D_H16O = 0.9723 ± 0.0007
    % D_HD16O/D_H216O = 0.9755 ± 0.0009
% D18O = 0.9723; % D'/D % Hellmann's and Harvey's # 0.96671
% D2H  = 0.9755; % D'/D % Hellmann's and Harvey's # 0.98258
  % Diffusivities (D) adopted from Hellmann and Harvey 2020:
    % D_H18O/D_H16O = 
    % D_HD16O/D_H216O = 
T = (nanmean(Ta)+273.15)./100; % in degrees Kelvin
D2H  = 0.98258 - (0.02546./T) + (0.02421./(T^(5/2)));
D18O = 0.96671 + (0.007406/(T^.5)) - (0.004861./(T^3));

% Not using the turbulent exponent (n) => same as assuming n = 1?
% h = 1; % yields the same result as using alpha_e alone
h = linspace(0.5,1,5); % range based on ??? data (from above molecular layer to 20m height)
alpha_k_O = (1/D18O)*(1-h)+ h; % alpha_k must be greater than 1; checked!!!
alpha_k_D = (1/D2H)*(1-h)+ h;  % alpha_k must be greater than 1; checked!!!
axa_D = alpha_e_D.*alpha_k_D; 
axa_O = alpha_e  .*alpha_k_O; 
% Solving for the isotopic concentration del_v0_D
delta_e_D = ((1./axa_D)*(1+del_oc_D) - 1) * 1000;
delta_e_O = ((1./axa_O)*(1+del_oc) - 1) * 1000;
% delta_e_D = ((1./alpha_e_D)*(1+del_oc_D) - 1) * 1000;
% delta_e_O = ((1./alpha_e)*(1+del_oc) - 1) * 1000;

for k = 1:22
    qsurf_cold_test(k,:) = q_surf(start_t_ids(k):start_t_ids(k)+6);
    dDsurf_cold(k,:)     = delta_e_D(start_t_ids(k):start_t_ids(k)+6,5);
    qent_cold(k,:)       = q_en(start_t_ids(k):start_t_ids(k)+6);
    d18Osurf_cold(k,:)   = delta_e_O(start_t_ids(k):start_t_ids(k)+6,5);
    DXSsurf_cold(k,:)    = dDsurf_cold(k,:) - (8*d18Osurf_cold(k,:));
end

%%
for k = 1:length(max_ind)
    dummy = timep(max_ind(k)-90:end_ind(k));
    tcold(k,1:length(dummy)) = timep(max_ind(k)-90:end_ind(k));
    qcold(k,1:length(dummy)) = q_iso(max_ind(k)-90:end_ind(k));
    dDcold(k,1:length(dummy))= dD_avg_1min(max_ind(k)-90:end_ind(k));
    DXScold(k,1:length(dummy))= d_ex_1min(max_ind(k)-90:end_ind(k));
    d18Ocold(k,1:length(dummy))= d18O_avg_1min(max_ind(k)-90:end_ind(k));
    cold_pool_mask_1min(max_ind(k):end_ind(k)) = 1;
end

%% Creating new variables where cold pools times are NaN's
qnan = q_iso;
qnan(cold_pool_mask_1min==1)= NaN;
dDnan = dD_avg_1min;
dDnan(cold_pool_mask_1min==1)= NaN;
d18Onan = d18O_avg_1min;
d18Onan(cold_pool_mask_1min==1)= NaN;
DXSnan = d_ex_1min;
DXSnan(cold_pool_mask_1min==1)= NaN;

for k = 1:length(max_ind)
    dummy = timep(max_ind(k)-90:end_ind(k));
    qcoldnan(k,1:length(dummy)) = qnan(max_ind(k)-90:end_ind(k));
    dDcoldnan(k,1:length(dummy))= dDnan(max_ind(k)-90:end_ind(k));
    d18Ocoldnan(k,1:length(dummy))= d18Onan(max_ind(k)-90:end_ind(k));
    DXScoldnan(k,1:length(dummy))= DXSnan(max_ind(k)-90:end_ind(k));    
end
% Restricting cold pool timeseries to the first 60 min after the cold pool front ends
qcoldnan(:,147:end) = [];
dDcoldnan(:,147:end) = [];
d18Ocoldnan(:,147:end) = [];
DXScoldnan(:,147:end) = [];

% 'f' subscript for front
for k = 1:length(max_ind)
    dummy = timep(max_ind(k):min_ind(k));
    tcoldf(k,1:length(dummy)) = timep(max_ind(k):min_ind(k));
    qcoldf(k,1:length(dummy)) = q_iso(max_ind(k):min_ind(k));
    dDcoldf(k,1:length(dummy))= dD_avg_1min(max_ind(k):min_ind(k));
    DXScoldf(k,1:length(dummy))= d_ex_1min(max_ind(k):min_ind(k));
    d18Ocoldf(k,1:length(dummy))= d18O_avg_1min(max_ind(k):min_ind(k));
end

% 'w' subscript for wake
for k = 1:length(max_ind)
    dummy = timep(min_ind(k):end_ind(k));
    tcoldw(k,1:length(dummy)) = timep(min_ind(k):end_ind(k));
    qcoldw(k,1:length(dummy)) = q_iso(min_ind(k):end_ind(k));
    dDcoldw(k,1:length(dummy))= dD_avg_1min(min_ind(k):end_ind(k));
    DXScoldw(k,1:length(dummy))= d_ex_1min(min_ind(k):end_ind(k)); 
    d18Ocoldw(k,1:length(dummy))= d18O_avg_1min(min_ind(k):end_ind(k));    
end

%% converting zeros in NaN's
tcoldf(tcoldf==0) = NaN;
qcoldf(qcoldf==0) = NaN;
dDcoldf(dDcoldf==0) = NaN;
DXScoldf(DXScoldf==0) = NaN;
d18Ocoldf(d18Ocoldf==0) = NaN;

tcoldw(tcoldw==0) = NaN;
qcoldw(qcoldw==0) = NaN;
dDcoldw(dDcoldw==0) = NaN;
DXScoldw(DXScoldw==0) = NaN;
d18Ocoldw(d18Ocoldw==0) = NaN;

tcold(tcold==0) = NaN;
qcold(qcold==0) = NaN;
dDcold(dDcold==0) = NaN;
DXScold(DXScold==0) = NaN;
d18Ocold(d18Ocold==0) = NaN;

qcoldnan(qcoldnan==0) = NaN;
% qscoldnan(qscoldnan==0) = NaN;
dDcoldnan(dDcoldnan==0) = NaN;
DXScoldnan(DXScoldnan==0) = NaN;
d18Ocoldnan(d18Ocoldnan==0) = NaN;

%%
% tcold(1,44:51) = NaN;
% qcold(1,44:51) = NaN;
% dDcold(1,44:51) = NaN;
% DXScold(1,44:51) = NaN;

%% restricting wake timeseries to the first 60 min after the cold pool front ends
tcoldw(:,61:end) = [];
qcoldw(:,61:end) = [];
dDcoldw(:,61:end) = [];
DXScoldw(:,61:end) = [];
d18Ocoldw(:,61:end) = [];

for k = 1:22
    dum = find(isnan(tcoldf(k,:)));
    if isempty(dum) == 1
        tcoldw(k,14:end) = NaN;
        qcoldw(k,14:end) = NaN;
        dDcoldw(k,14:end) = NaN;
        DXScoldw(k,14:end) = NaN;
        d18Ocoldw(k,14:end) = NaN;
    else
        ind = 60-dum(1);
        tcoldw(k,ind+1:end) = NaN;
        qcoldw(k,ind+1:end) = NaN;
        dDcoldw(k,ind+1:end) = NaN;
        DXScoldw(k,ind+1:end) = NaN;
        d18Ocoldw(k,ind+1:end) = NaN;
    end
end

%% Identifying centroid for 90-minute data before cold pool onset
% SIMPLE APPROACH: Assuming dD and q at the cold pool onset are the same as dD and q for the centroid
q_cp_mean = qcoldf(:,1);
dD_cp_mean = (qcoldf(:,1).*dDcoldf(:,1))./10^3;
d18O_cp_mean = (qcoldf(:,1).*d18Ocoldf(:,1))./10^3;
DXS_cp_mean = (qcoldf(:,1).*DXScoldf(:,1))./10^3;

% CENTROID APPROACH: Using centroid in the q*dD-q space
q_cp_mean = nanmean(qcoldnan(:,1:90),2);
% qs_cp_mean = nanmean(qscoldnan(:,1:90),2);
dD_cp_mean = nanmean((qcoldnan(:,1:90).*dDcoldnan(:,1:90)),2)./10^3;
d18O_cp_mean = nanmean((qcoldnan(:,1:90).*d18Ocoldnan(:,1:90)),2)./10^3;
DXS_cp_mean = nanmean((qcoldnan(:,1:90).*DXScoldnan(:,1:90)),2)./10^3;

% Re-defining qsurf_cold in terms of varying RH
         % (20m = roughly where the Picarro sampled)

Y1 = (qsurf_cold(:,1).*d18Osurf_cold(:,1))./10^3;
Y2 = (qsurf_cold(:,1).*dDsurf_cold(:,1))./10^3;
Y3 = (qsurf_cold(:,1).*DXSsurf_cold(:,1))./10^3;
m1 = (d18O_cp_mean - Y1)./(q_cp_mean - qsurf_cold(:,1)); % slope
m2 = (dD_cp_mean - Y2)./(q_cp_mean - qsurf_cold(:,1)); % slope
m3 = (DXS_cp_mean - Y3)./(q_cp_mean - qsurf_cold(:,1)); % slope
b = d18O_cp_mean - (m1.*q_cp_mean);
b1 = Y1 - (m1.*qsurf_cold(:,1)); % same results for both points: good check
bb = dD_cp_mean - (m2.*q_cp_mean);
b2 = Y2 - (m2.*qsurf_cold(:,1)); % same results for both points: good check
bbb = DXS_cp_mean - (m3.*q_cp_mean);
b3 = Y3 - (m3.*qsurf_cold(:,1)); % same results for both points: good check
qXd18O_ent = (m1.*qent_cold(:,1)) + b1;
qXdD_ent   = (m2.*qent_cold(:,1)) + b2;
qXDXS_ent  = (m3.*qent_cold(:,1)) + b3;
d18O_ent = (qXd18O_ent.*10^3)./qent_cold(:,1);
dD_ent = (qXdD_ent.*10^3)./qent_cold(:,1);
DXS_ent = (qXDXS_ent.*10^3)./qent_cold(:,1);

%% Plots for dD
for i = 13 %[8,13,20] %1:22
    % Calculating slope of best fit line for the 90-min data 
    slope = polyfit(qcoldnan(i,1:90),dDcoldnan(i,1:90),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,qcoldnan(i,1:90));
figure; set(gcf, 'Position', get(0, 'Screensize'));
j = plot(qcoldf(i,:),dDcoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(qcoldw(i,:),dDcoldw(i,:),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 90 minutes) in gray dots
jj = scatter(qcoldnan(i,1:90),dDcoldnan(i,1:90),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
jjj = scatter(qcoldf(i,:),dDcoldf(i,:),70,(tcoldf(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
j4 = scatter(qcoldw(i,:),dDcoldw(i,:),70,(tcoldw(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
j4.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([8 23]); ylim([-78 -64])
% ylim([-95 -45]) % ylim([6 12.5])

Xfit = linspace(qent_cold(i,1),qsurf_cold(i,1),10);
Ylinearfit = (m2(i,1).*Xfit) + b2(i,1);
Ycurvefit = (Ylinearfit.*10^3)./Xfit;

% Adding mixing lines and end member points (entraiment & surface)
plot(qent_cold(i,1),dD_ent(i),'ok','MarkerFaceColor',[0.2, 0, 0],'MarkerSize',10)
j3 = plot(Xfit,Ycurvefit,'--k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qsurf_cold(i,1),dDsurf_cold(i,1),'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)
% plot(nanmean(q_ob),nanmean(dD),'ok','MarkerFaceColor','w','MarkerSize',10)

colormap(jet(12))
han = colorbar;
han.Title.String = "minutes since cp onset";
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
caxis([0 60])

hold on;
plot(qcoldf(i,:),dDcoldf(i,:),'ok','MarkerSize',8);%,'MarkerFaceColor','b')
plot(qcoldw(i,:),dDcoldw(i,:),'sk','MarkerSize',10);%,'MarkerFaceColor','c')
plot(qcoldf(i,1),dDcoldf(i,1),'^k','MarkerFaceColor','k','MarkerSize',12)
plot(qcoldw(i,1),dDcoldw(i,1),'^k','MarkerFaceColor','w','MarkerSize',12)
plot(qcoldnan(i,1:90),y_est,'r-','LineWidth',2)

% LIMITS FOR ZOOM-IN plots
% ylim([min(dDcoldnan(i,1:90))-.25 max(dDcoldnan(i,1:90))+3.75])
xlabel('specific humidity q [g/kg]')
ylabel('dD [per mil]')

clearvars y_est slope
box on
grid on
legend('entrainment','surface','cp front','cp wake','front onset','wake onset','Location','northwest')

title(['Cold Pool #',num2str(i),'; onset on ',datestr(tcoldf(i,1))])
% saveas(gcf,['dD_vs_q_CP#',num2str(i),'_range_RH.png'])
end

%% Plots for d18O
for i = 13 %[8,10,11,13:17] %1:22
    % Calculating slope of best fit line for the 90-min data 
    slope = polyfit(qcoldnan(i,1:90),d18Ocoldnan(i,1:90),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,qcoldnan(i,1:90));
figure; set(gcf, 'Position', get(0, 'Screensize'));
j = plot(qcoldf(i,:),d18Ocoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(qcoldw(i,:),d18Ocoldw(i,:),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 90 minutes) in gray dots
jj = scatter(qcoldnan(i,1:90),d18Ocoldnan(i,1:90),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
jjj = scatter(qcoldf(i,:),d18Ocoldf(i,:),70,(tcoldf(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
j4 = scatter(qcoldw(i,:),d18Ocoldw(i,:),70,(tcoldw(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
j4.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([8 23]) ; %ylim([-78 -64])
% ylim([-95 -45]) % ylim([6 12.5])

Xfit = linspace(qent_cold(i,1),qsurf_cold(i,1),10);
Ylinearfit = (m1(i,1).*Xfit) + b1(i,1);
Ycurvefit = (Ylinearfit.*10^3)./Xfit;

% Adding mixing lines and end member points (entraiment & surface)
plot(qent_cold(i,1),d18O_ent(i),'ok','MarkerFaceColor',[0.2, 0, 0],'MarkerSize',10)
j3 = plot(Xfit,Ycurvefit,'--k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qsurf_cold(i,1),d18Osurf_cold(i,1),'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)

colormap(jet(12))
han = colorbar;
han.Title.String = "minutes since cp onset";
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
caxis([0 60])

plot(qcoldf(i,:),d18Ocoldf(i,:),'ok','MarkerSize',8);%,'MarkerFaceColor','b')
plot(qcoldw(i,:),d18Ocoldw(i,:),'sk','MarkerSize',10);%,'MarkerFaceColor','c')
plot(qcoldf(i,1),d18Ocoldf(i,1),'^k','MarkerFaceColor','k','MarkerSize',12)
plot(qcoldw(i,1),d18Ocoldw(i,1),'^k','MarkerFaceColor','w','MarkerSize',12)
plot(qcoldnan(i,1:90),y_est,'r-','LineWidth',2)

% xlim([285 305])
% ylim([min(d18Ocoldnan(i,1:90))-.25 max(d18Ocoldnan(i,1:90))+3.75])
xlabel('specific humidity q [g/kg]')
ylabel('d18O [per mil]')

clearvars y_est slope
box on
grid on
legend('entrainment','surface','cp front','cp wake','front onset','wake onset','Location','northwest')

title(['Cold Pool #',num2str(i),'; onset on ',datestr(tcoldf(i,1))])
% saveas(gcf,['d18O_vs_q_CP#',num2str(i),'_new_ff.png'])
end

%% Plots for DXS
for i = 13 %[8,10,11,13:17] %1:22
    % Calculating slope of best fit line for the 90-min data 
    slope = polyfit(qcoldnan(i,1:90),DXScoldnan(i,1:90),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,qcoldnan(i,1:90));
figure; set(gcf, 'Position', get(0, 'Screensize'));
j = plot(qcoldf(i,:),DXScoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(qcoldw(i,:),DXScoldw(i,:),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 90 minutes) in gray dots
jj = scatter(qcoldnan(i,1:90),DXScoldnan(i,1:90),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
jjj = scatter(qcoldf(i,:),DXScoldf(i,:),70,(tcoldf(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
j4 = scatter(qcoldw(i,:),DXScoldw(i,:),70,(tcoldw(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
j4.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([8 23]) % xlim([-78 -64])
% ylim([-95 -45]) % ylim([6 12.5])

Xfit = linspace(qent_cold(i,1),qsurf_cold(i,1),10);
Ylinearfit = (m3(i,1).*Xfit) + b3(i,1);
Ycurvefit = (Ylinearfit.*10^3)./Xfit;

% Adding mixing lines and end member points (entraiment & surface)
plot(qent_cold(i,1),DXS_ent(i),'ok','MarkerFaceColor',[0.2, 0, 0],'MarkerSize',10)
j3 = plot(Xfit,Ycurvefit,'--k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qsurf_cold(i,1),DXSsurf_cold(i,1),'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)

colormap(jet(12))
han = colorbar;
han.Title.String = "minutes since cp onset";
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
caxis([0 60])

plot(qcoldf(i,:),DXScoldf(i,:),'ok','MarkerSize',8);%,'MarkerFaceColor','b')
plot(qcoldw(i,:),DXScoldw(i,:),'sk','MarkerSize',10);%,'MarkerFaceColor','c')
plot(qcoldf(i,1),DXScoldf(i,1),'^k','MarkerFaceColor','k','MarkerSize',12)
plot(qcoldw(i,1),DXScoldw(i,1),'^k','MarkerFaceColor','w','MarkerSize',12)
plot(qcoldnan(i,1:90),y_est,'r-','LineWidth',2)

% xlim([285 305])
% ylim([min(DXScoldnan(i,1:90))-.25 max(DXScoldnan(i,1:90))+3.75])
xlabel('specific humidity q [g/kg]')
ylabel('DXS [per mil]')

clearvars y_est slope
box on
grid on
legend('entrainment','surface','cp front','cp wake','front onset','wake onset','Location','northeast')

title(['Cold Pool #',num2str(i),'; onset on ',datestr(tcoldw(i,1))])
% saveas(gcf,['DXS_vs_q_CP#',num2str(i),'_new_ff.png'])
end