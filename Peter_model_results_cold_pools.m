load ColdPoolTimeseries_9Feb_hg.mat
time = out.time(:,1); % 1-min data in days starting at Feb 09 2020
date_t = time + datenum('dec-31-2019','mmm-dd-yyyy'); % in MatLab time 
R_VSMOW = 155.76*10^(-6); % FOR Deuterium

rD = out.HDOV_LOW(:,1); % HDO water vapor mass mixing ratio
r = out.QV_LOW(:,1); % WATER VAPOR MASS MIXING RATIO
dD = ((rD./r)-1)*10^3;
rO = out.O18V_LOW(:,1); % O18 water vapor mass mixing ratio
d18O = ((rO./r)-1)*10^3;
DXS = dD - 8*d18O;

q = r*10^3;
x = out.x(:,1);
y = out.y(:,1);
Ta = out.TABS_LOW(:,1);

%% Identifying cold pools
Taf = movmean(Ta,11,'omitnan'); % Ta filtered => 11-min running average
qairf = movmean(q(:,1),11,'omitnan'); % qair filtered => 11-min running average
dDf = movmean(dD(:,1),11,'omitnan'); % wspd filtered => 11-min running average
d18Of = movmean(d18O(:,1),11,'omitnan'); % qs filtered => 11-min running average
DXSf  = movmean(DXS(:,1),11,'omitnan'); % qs filtered => 11-min running average

figure; 
subplot(511);hold on;
    plot(out.time(:,1),Ta(:,1));
    plot(out.time(:,1),Taf(:,1),'LineWidth',2)
    ylabel('Abs Temp [K]')
subplot(512);hold on;
    plot(out.time(:,1),q(:,1))
    plot(out.time(:,1),qairf(:,1),'LineWidth',2)
    ylabel('q [g/kg]')
subplot(513);hold on;
    plot(out.time(:,1),dD(:,1))
    plot(out.time(:,1),dDf(:,1),'LineWidth',2)
    ylabel('\deltaD [permil]')
subplot(514);hold on;
    plot(out.time(:,1),d18O(:,1))
    plot(out.time(:,1),d18Of(:,1),'LineWidth',2)
    ylabel('\delta^1^8O [permil]')
subplot(515);hold on;
    plot(out.time(:,1),DXS(:,1))
    plot(out.time(:,1),DXSf(:,1),'LineWidth',2)
    ylabel('DXS [permil]')
xlabel('time in DOY')

del_T = Taf(2:end)-Taf(1:end-1);
cand = find(del_T<-0.05); % candidates positions

%% Separating candidates [5 possible cold pools]
c=1;
cand_ind(c) = cand(1);
for i = 2:length(cand)-1
    if cand(i)-cand(i-1)>1 && cand(i+1)-cand(i)==1
        c = c+1;
        cand_ind(c) = cand(i);
    end
end

% Time passed between possible cold pools
t_pass = time(cand_ind(2:end)) - time(cand_ind(1:end-1)); % in days
t_pass = t_pass*24; % in hours
t_max_offset = 20; % in minutes
cand_ind(2) = []; % to join candidates 8 and 11 into a single cold pool

%% Identifying t_max: time of the cold pool front onset
c=1;
for k = 1:length(cand_ind)
    for ii = cand_ind(k):-1:1
        if del_T(ii)>0 && ((time(cand_ind(k))-time(ii+1))*24*60)<t_max_offset
            t_max(c) = time(ii+1);
            t_max_ind(c) = ii+1;
            c = c+1;
            break
        elseif del_T(ii)>0 && ((time(cand_ind(k))-time(ii+1))*24*60)>=t_max_offset
            t_max(c) = time(cand_ind(k)-t_max_offset);
            t_max_ind(c) = (cand_ind(k))-t_max_offset;
            c = c+1;
            break
        end
    end
end
% Time passed between t_max times
% t_p = t1min(cand_ind) - t_max; % in days
% t_p = t_p*24*60; % in min
t_p = cand_ind - t_max_ind; % in min
find(t_p>20);

%% Identifying t_min: time of end of cold pool front
c=1;
for k = 1:length(cand_ind)
    for ii = cand_ind(k):1:length(del_T)
        if del_T(ii)>0
            t_min(c) = time(ii);
            t_min_ind(c) = ii;
            c = c+1;
            break
        elseif ii == length(del_T)
            t_min(c) = time(ii);
            t_min_ind(c) = ii;
            break
        end
    end
end
t_min_offset = 20; % within 20min of the previous minimum

%% Identifying t_end: end of cold pool
delta_T = Taf(t_max_ind) - Taf(t_min_ind);
T_min = Taf(t_min_ind);
T_max = Taf(t_max_ind);

for k = 1:length(cand_ind)-1
    for ii = t_min_ind(k):1:cand_ind(k+1)
%         if Taf(ii) > T_min(k) + delta_T(k)/exp(1) % for recovery defined as when ~1/3 of delta_T is reached; 1/3 = 0.3333; 1/exp(1) = 0.3679 = e^(-1); 
        if Taf(ii) > T_min(k) + delta_T(k)*(1-(1/exp(1))) % for recovery defined as when ~2/3 of delta_T is reached; 2/3 = 0.6667; 1-(1/exp(1)) = 0.6321 = 1 - e^(-1);
            t_end(k) = time(ii);
            t_end_ind(k) = ii;
            end_flag(k) = 1; % 1 = cold pool ends on its own; recovered cold pool
            break
        end
        t_end(k) = time(ii-1);
        t_end_ind(k) = ii-1;
        end_flag(k) = 0; % cold pool ends at the onset of the next cold pool
    end
end
k = length(cand_ind);
for ii = t_min_ind(k):1:length(time)
%     if Taf(ii) > T_min(k) + delta_T(k)/exp(1) % for recovery defined as when ~1/3 of delta_T is reached; 1/3 = 0.3333; 1/exp(1) = 0.3679 = e^(-1); 
    if Taf(ii) > T_min(k) + delta_T(k)*(1-(1/exp(1))) % for recovery defined as when ~2/3 of delta_T is reached; 2/3 = 0.6667; 1-(1/exp(1)) = 0.6321 = 1 - e^(-1);
        t_end(k) = time(ii);
        t_end_ind(k) = ii;
        end_flag(k) = 1; % 1 = cold pool ends on its own; recovered cold pool
        break
    end
    t_end(k) = time(ii-1);
    t_end_ind(k) = ii-1;
    end_flag(k) = 0; % cold pool ends at the onset of the next cold pool
end

%% FLAG for cold pool times (indexes)
t = time(:,1);
cold_pool_flag_1min = zeros(size(t));
cp_matrix = zeros(length(t_max),61); % matrix of cold pool indexes
for k = 1:length(t_max)
    max_ind(k) = find(t==t_max(k));
    min_ind(k) = find(t==t_min(k));
    end_ind(k) = find(t==t_end(k));
    ii = max_ind(k):end_ind(k);
    cp_matrix(k,1:length(ii)) = ii;
    cold_pool_flag_1min(max_ind(k):end_ind(k)) = 1;
end
cp_matrix(cp_matrix==0) = NaN;
t = t + datenum('20191231','yyyymmdd');

% for plotting purposes
figure;
hold on;
for k = 1:length(t_max)
    patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 2 2 -1],[.8 .8 .8],'EdgeColor','none')
end
plot(t,cold_pool_flag_1min,'-r','Color','r')
ylim([0 1])
datetick('x','HH:MM','keeplimits','keepticks')
title('cold pool flag'); 
xlabel('Feb-09-2020')
grid on

%%
figure;
subplot(511); 
    hold on;
    for k = 1:length(t_max)
        patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*max(Ta(:,1)),[.8 .8 .8],'EdgeColor','none')
    end
    plot(date_t(:,1),Ta(:,1))
    plot(date_t(:,1),Taf(:,1),'LineWidth',2)
    ylim([min(Ta(:,1)) max(Ta(:,1))])
    ylabel('Abs Temp [K]')
    datetick('x','HH:MM','keeplimits','keepticks')
subplot(512); 
    hold on;
    for k = 1:length(t_max)
        patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*max(q(:,1)),[.8 .8 .8],'EdgeColor','none')
    end
    plot(date_t(:,1),q(:,1))
    plot(date_t(:,1),qairf(:,1),'LineWidth',2)
    ylim([min(q(:,1)) max(q(:,1))])
    ylabel('q [g/kg]')
    datetick('x','HH:MM','keeplimits','keepticks')
subplot(513); 
    hold on;
    for k = 1:length(t_max)
        patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*min(dD(:,1)),[.8 .8 .8],'EdgeColor','none')
    end
    plot(date_t(:,1),dD(:,1))
    plot(date_t(:,1),dDf(:,1),'LineWidth',2)
    ylim([min(dD(:,1)) max(dD(:,1))])
    ylabel('\deltaD [permil]')
    datetick('x','HH:MM','keeplimits','keepticks')
subplot(514); 
    hold on;
    for k = 1:length(t_max)
        patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*min(d18O(:,1)),[.8 .8 .8],'EdgeColor','none')
    end
    plot(date_t(:,1),d18O(:,1))
    plot(date_t(:,1),d18Of(:,1),'LineWidth',2)
    ylim([min(d18O(:,1)) max(d18O(:,1))])
    ylabel('\delta^1^8O [permil]')
    datetick('x','HH:MM','keeplimits','keepticks')
subplot(515); 
    hold on;
    for k = 1:length(t_max)
        patch([t(max_ind(k)) t(max_ind(k)) t(end_ind(k)) t(end_ind(k))],[-1 1 1 -1].*min(DXS(:,1)),[.8 .8 .8],'EdgeColor','none')
    end
    plot(date_t(:,1),DXS(:,1))
    plot(date_t(:,1),DXSf(:,1),'LineWidth',2)
    ylim([min(DXS(:,1)) max(DXS(:,1))])
    ylabel('\DXS [permil]')
    datetick('x','HH:MM','keeplimits','keepticks')
xlabel('Feb-09-2020')
set(findall(gcf,'-property','Fontsize'),'FontSize',14)

%%
k = 1;
    dummy = t(max_ind(k)-6:end_ind(k));
    tcold(k,1:length(dummy)) = t(max_ind(k)-6:end_ind(k));
    qcold(k,1:length(dummy)) = q(max_ind(k)-6:end_ind(k));
    dDcold(k,1:length(dummy))= dD(max_ind(k)-6:end_ind(k));
    DXScold(k,1:length(dummy))= DXS(max_ind(k)-6:end_ind(k));
    d18Ocold(k,1:length(dummy))= d18O(max_ind(k)-6:end_ind(k));

for k = 2:length(max_ind)
    dummy = t(max_ind(k)-61:end_ind(k));
    tcold(k,1:length(dummy)) = t(max_ind(k)-61:end_ind(k));
    qcold(k,1:length(dummy)) = q(max_ind(k)-61:end_ind(k));
    dDcold(k,1:length(dummy))= dD(max_ind(k)-61:end_ind(k));
   DXScold(k,1:length(dummy))= DXS(max_ind(k)-61:end_ind(k));
  d18Ocold(k,1:length(dummy))= d18O(max_ind(k)-61:end_ind(k));
end
%% Creating new variables where cold pools times are NaN's
qnan = q(:,1);
qnan(cold_pool_flag_1min==1)= NaN;
dDnan = dD(:,1);
dDnan(cold_pool_flag_1min==1)= NaN;
d18Onan = d18O(:,1);
d18Onan(cold_pool_flag_1min==1)= NaN;
DXSnan = DXS(:,1);
DXSnan(cold_pool_flag_1min==1)= NaN;

for k = 2:length(max_ind)
    dummy = t(max_ind(k)-61:end_ind(k));
    qcoldnan(k,1:length(dummy)) = qnan(max_ind(k)-61:end_ind(k));
    dDcoldnan(k,1:length(dummy))= dDnan(max_ind(k)-61:end_ind(k));
    d18Ocoldnan(k,1:length(dummy))= d18Onan(max_ind(k)-61:end_ind(k));
    DXScoldnan(k,1:length(dummy)) = DXSnan(max_ind(k)-61:end_ind(k));    
end

% 'f' subscript for front
for k = 1:length(max_ind)
    dummy = t(max_ind(k):min_ind(k));
    tcoldf(k,1:length(dummy)) = t(max_ind(k):min_ind(k));
    qcoldf(k,1:length(dummy)) = q(max_ind(k):min_ind(k));
    dDcoldf(k,1:length(dummy))= dD(max_ind(k):min_ind(k));
   DXScoldf(k,1:length(dummy))= DXS(max_ind(k):min_ind(k));
  d18Ocoldf(k,1:length(dummy))= d18O(max_ind(k):min_ind(k));
end

% 'w' subscript for wake
for k = 1:length(max_ind)
    dummy = t(min_ind(k):end_ind(k));
    tcoldw(k,1:length(dummy)) = t(min_ind(k):end_ind(k));
    qcoldw(k,1:length(dummy)) = q(min_ind(k):end_ind(k));
    dDcoldw(k,1:length(dummy))= dD(min_ind(k):end_ind(k));
   DXScoldw(k,1:length(dummy))= DXS(min_ind(k):end_ind(k));    
  d18Ocoldw(k,1:length(dummy))= d18O(min_ind(k):end_ind(k));
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
dDcoldnan(dDcoldnan==0) = NaN;
DXScoldnan(DXScoldnan==0) = NaN;
d18Ocoldnan(d18Ocoldnan==0) = NaN;

%% 
load('RHS_Eq9_MJ79_LIMITED.mat')
load('conserved_variables_10minLIMITED.mat')

k = 19; % all model cols pools using the same surface properties
delta_e_D = ((1./alpha_e_D)*(1+del_oc_D) - 1) * 1000;
delta_e_O = ((1./alpha_e)*(1+del_oc) - 1) * 1000;

qsurf_cold      = q_surf(start_t_ids(k));
dDsurf_cold     = delta_e_D(start_t_ids(k));
qent_cold       = q_en(start_t_ids(k));
d18Osurf_cold   = delta_e_O(start_t_ids(k));
DXSsurf_cold    = dDsurf_cold - (8*d18Osurf_cold);

%% Identifying centroid for 90-minute data before cold pool onset
% SIMPLE APPROACH: Assuming dD and q at the cold pool onset are the same as dD and q for the centroid
q_cp_mean = qcoldf(:,1);
dD_cp_mean = (qcoldf(:,1).*dDcoldf(:,1))./10^3;
d18O_cp_mean = (qcoldf(:,1).*d18Ocoldf(:,1))./10^3;
DXS_cp_mean = (qcoldf(:,1).*DXScoldf(:,1))./10^3;

% CENTROID APPROACH: Using centroid in the q*dD-q space
q_cp_mean = nanmean(qcoldnan(:,1:61),2);
dD_cp_mean = nanmean((qcoldnan(:,1:61).*dDcoldnan(:,1:61)),2)./10^3;
d18O_cp_mean = nanmean((qcoldnan(:,1:61).*d18Ocoldnan(:,1:61)),2)./10^3;
DXS_cp_mean = nanmean((qcoldnan(:,1:61).*DXScoldnan(:,1:61)),2)./10^3;

Y1 = (qsurf_cold.*d18Osurf_cold)./10^3;
Y2 = (qsurf_cold.*dDsurf_cold)./10^3;
Y3 = (qsurf_cold(:,1).*DXSsurf_cold)./10^3;
m1 = (d18O_cp_mean - Y1)./(q_cp_mean - qsurf_cold); % slope
m2 = (dD_cp_mean - Y2)./(q_cp_mean - qsurf_cold); % slope
m3 = (DXS_cp_mean - Y3)./(q_cp_mean - qsurf_cold); % slope
b = d18O_cp_mean - (m1.*q_cp_mean);
b1 = Y1 - (m1.*qsurf_cold); % same results for both points: good check
bb = dD_cp_mean - (m2.*q_cp_mean);
b2 = Y2 - (m2.*qsurf_cold); % same results for both points: good check
bbb = DXS_cp_mean - (m3.*q_cp_mean);
b3 = Y3 - (m3.*qsurf_cold); % same results for both points: good check
qXd18O_ent = (m1.*qent_cold) + b1;
qXdD_ent   = (m2.*qent_cold) + b2;
qXDXS_ent  = (m3.*qent_cold) + b3;
d18O_ent = (qXd18O_ent.*10^3)./qent_cold;
dD_ent = (qXdD_ent.*10^3)./qent_cold;
DXS_ent = (qXDXS_ent.*10^3)./qent_cold;

%% Plots for dD
for i = 3
    % Calculating slope of best fit line for the 60-min data 
    slope = polyfit(qcoldnan(i,1:61),dDcoldnan(i,1:61),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,qcoldnan(i,1:61));
figure; set(gcf, 'Position', get(0, 'Screensize'));
j = plot(qcoldf(i,:),dDcoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(qcoldw(i,:),dDcoldw(i,:),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 60 minutes) in gray dots
jj = scatter(qcoldnan(i,1:61),dDcoldnan(i,1:61),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
jjj = scatter(qcoldf(i,:),dDcoldf(i,:),70,(tcoldf(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
j4 = scatter(qcoldw(i,:),dDcoldw(i,:),70,(tcoldw(i,:)-tcoldf(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
j4.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([8 23]); %ylim([-90 -64])

Xfit = linspace(qent_cold,qsurf_cold,10);
Ylinearfit = (m2(i).*Xfit) + b2(i);
Ycurvefit = (Ylinearfit.*10^3)./Xfit;

% Adding mixing lines and end member points (entraiment & surface)
plot(qent_cold,dD_ent(i),'ok','MarkerFaceColor',[0.2, 0, 0],'MarkerSize',10)
j3 = plot(Xfit,Ycurvefit,'-k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qsurf_cold,dDsurf_cold,'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)
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
plot(qcoldnan(i,1:61),y_est,'r-','LineWidth',2)

% xlim([285 305])
% ylim([min(dDcold(i,1:90))-.25 max(dDcold(i,1:90))+3.75])
xlabel('specific humidity q [g/kg]')
ylabel('dD [per mil]')
% ylabel('DXS [per mil]')
% xlabel('dD [per mil]')

clearvars y_est slope
box on
grid on
legend('entrainment','surface','cp front','cp wake','front onset','wake onset','Location','northwest')

title(['Model cp #',num2str(i),'; onset on ',datestr(tcoldf(i,1))])
saveas(gcf,['model&obs_dD_vs_q_CP#',num2str(i),'_curve_fit_centroid.png'])
end
