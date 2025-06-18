%% mixing fractions in isotopic space q-q*dD %%
% [fen, fss, fdd] = th_q_to_mixfraction(th_ob,q_ob, th_en,q_en, th_surf,q_surf, th_d,q_d);
% figure; hold on; % for all cps in single plot
% set(gcf, 'Position', [1 1 750 650]);

q_dD_framework_cold_pools

% Allocating "slope" variables
b=1:17;
c=1:61;
slope = 9999*ones(size(b));
int = 9999*ones(size(b));
Error_slope = 9999*ones(size(b));
% allocating variables
fen = 9999*ones(size(c));
fss = 9999*ones(size(c));
fee = 9999*ones(size(c));
t_cp = 9999*ones(size(c));

%% for iso mixing fraction timeseries plot
    % iso_fen = 9999*ones(size(b));
    % iso_fdd = 9999*ones(size(b));
    % iso_fss = 9999*ones(size(b));

for i = [2:11,13:16] %[7,8] %[9,14:15] %[8,10,11,13,16] %1:length(max_ind) %[2:7,9,14:15] %1:length(max_ind)-1 %
    q_bg  = qcoldnan(i,1:61);
    dD_bg = q_bg.*dDcoldnan(i,1:61);
    q_cp  = qcold(i,1:61);
    dD_cp = q_cp.*dDcold(i,1:61);
    q_en  = qent_cold(i)*ones(1,61);
    dD_en = q_en.*(dD_ent(i)*ones(1,61));
    q_surf  = qsurf_cold(i)*ones(1,61);
    dD_surf = q_surf.*(dDsurf_cold(i)*ones(1,61));
    t_cp(i,:) = tcold(i,1:61);
    
    dD_prec = 10; % in permil
    [q_lost,dD_lost,f_lost] = Rayleigh_liquid_evap(p0p(i),q0p(i),T0p(i),dD_prec); % p[Pa],q[g/kg],T[K],dD[permil]
    f = 0.25; % remaining fraction of liquid drop
    q_ev  = q_lost(f_lost==f)*ones(1,61);
    dD_ev = q_ev.*(dD_lost(f_lost==f)*ones(1,61));
    [fen, fss, fee] = x_y_to_mixfraction(dD_cp,q_cp, dD_en,q_en, dD_surf,q_surf, dD_ev,q_ev);
    [fen_bg, fss_bg, fee_bg] = x_y_to_mixfraction(dD_bg,q_bg, dD_en,q_en, dD_surf,q_surf, dD_ev,q_ev);
    fdd=1-fen-fss;
    fdd_bg = 1-fen_bg-fss_bg; %?

% for iso mixing fraction timeseries plot
% iso_fen(i) = fen(1);
% iso_fss(i) = fss(1);
% iso_fdd(i) = fdd(1);
   
    % Calculating the mass of water going into the mixing fractions
    % This methodology only applies to when fen == 0, correct?
    % Otherwise the linear relationship does not hold anymore, right?
    ind = find(fen>0 & fen<0.1); % for when fen == 0
    dummy = 1:length(fen); % for the entire dataset, for plotting limits
    % [fw1,fw2] = mass_mixing_test(ind,fee,fss,q_ev,q_surf);
    [fwd1,fwd2] = mass_mixing_test(dummy,fee,fss,q_ev,q_surf);
    [fw1,fw2,fw3,q] = mass_mixing(dummy,fee,fss,fen,q_ev,q_surf,q_en);
    %     fa1 = fee(ind);
    %     fa2 = fss(ind);
    %     q1 = q_ev(ind);
    %     q2 = q_surf(ind);
    %     q = fa1.*q1 + fa2.*q2;
    %     fw1 = (fa1.*q1)./q;
    %     fw2 = (fa2.*q2)./q;

    % Calculating linear fit in mixing fraction space
    X = fss(~isnan(fss));
    Y = fee(~isnan(fee));
    % N = 1; % linear polynomial
    % [P,S,MU] = polyfit(X,Y,N);
    % Xeval = 0:0.1:1;
    % cloud_fit = polyval(P,Xeval,S,MU);
    % slope(i) = (cloud_fit(end)-cloud_fit(1))./(Xeval(end)-Xeval(1));
    % cloud_fit2 = slope(i).*Xeval + cloud_fit(1);
    %     n = length(X);
    %     slope(i) = (n*sum(X.*Y)-(sum(X).*sum(Y)))./((n*sum(X.^2))-(sum(X)).^2);
    %     int(i) = ((sum(X.^2)*sum(Y))-(sum(X)*sum(X.*Y)))./((n*sum(X.^2))-(sum(X)).^2);
    %     S = sqrt(sum((Y-(slope(i).*X)-int(i)).^2)./(n-2));
    %     Error_slope(i) = S.*sqrt(n./((n*sum(X.^2))-(sum(X).^2)));
    %     disp((Error_slope(i)/slope(i))*100)
    mdl = fitlm(X,Y);
    Error_slope(i) = mdl.Rsquared.Adjusted; % R-squared
    slope(i) = mdl.Coefficients.Estimate(2);
%     max_fen(i) = max(fen,[],'omitnan');
%     max_fss(i) = max(fss,[],'omitnan');    
%     max_fee(i) = max(fee,[],'omitnan');
%     min_fen(i) = min(fen,[],'omitnan');
%     min_fss(i) = min(fss,[],'omitnan');    
%     min_fee(i) = min(fee,[],'omitnan');
    Xeval = 0:0.1:1;
    cloud_fit2 = slope(i).*Xeval + mdl.Coefficients.Estimate(1);

% Plots for isotopes vapor mixing fractions %%
    % Using specific humidity
    % set(gcf, 'Position', [1 1 900 500]);
    figure; hold on;
    j = scatter(fss_bg,fee_bg,15,[.5 .5 .5],'filled');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    j = scatter(fss,fee,30,'k');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % scatter(fss,fee,25,dD_cp,'filled');
    j = scatter(fss,fee,25,dD_cp./q_cp,'filled'); % color coded by dD
    % j = scatter(fss,fee,25,(t_cp-t_cp(1))*24*60,'filled'); % COLOR CODED BY TIME
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    colormap(flip(jet(15)))
    han = colorbar;
    % han.YLabel.String = 'q*\deltaD';% [',char(8240),']'];
    han.YLabel.String = ['\deltaD [',char(8240),']'];
%     han.YLabel.String = ('minutes since cp onset');
    j = plot([0 1],[0 0],'-k');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    j = plot([0 0],[0 1],'-k');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    j = plot([1 0],[0 1],'-k');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % set(gca,'Xdir','reverse')
    plot(Xeval,cloud_fit2,'-','LineWidth',1)
    text(0.45,0.8,['slope = ',num2str(round(slope(i),2)),''])
    xlim([-0.2 1])
    ylim([-0.2 1])
    xlabel('surface fraction')
    ylabel('evaporatee fraction')
    text(0.6,0.9,['CP#',num2str(i),''])
    set(findall(gcf,'-property','Fontsize'),'FontSize',18)
    box on
    title('BL air mixture fraction')
    axis square
    caxis([min(dD_cp./q_cp) max(dD_cp./q_cp)])
    scatter(0,1,55,'k')% fee=1; fss = 0 and fen = 0
    scatter(1,0,55,[.5 .5 .5],'filled')% fss=1; fee = 0 and fen = 0
    scatter(0,0,55,'k','filled')% fen=1; fee = 0 and fss = 0
    %legend('linear fit','Location','southeast')
% Saving figure in different formats
%     saveas(gcf,['mixing_fractions_qdD_vs_q_CP#',num2str(i),'_wslope.png'])
%     saveas(gcf,['mixing_fractions_qdD_vs_q_CP#',num2str(i),'_wslope.fig'])
end

%% Allocating variables for calculating cold pool age estimates in a separate script
% Script name: cold_pool_age_estimates.m
% Variables saved:
mf_fee(i,1:61) = fee; % mf stands for mixing fractions
mf_fen(i,1:61) = fen;
mf_fss(i,1:61) = fss;
mf_fee_bg(i,1:61) = fee_bg; % mf stands for mixing fractions
mf_fen_bg(i,1:61) = fen_bg;
mf_fss_bg(i,1:61) = fss_bg;
mf_dD_cp(i,1:61) = dD_cp./q_cp;
slope(slope==9999) = NaN;
slope = slope(1:16);
mf_t_cp = t_cp;
% Xeval;
Yeval(i,1:11) = cloud_fit2;
% end

% Saving variables for calculating cold pool age estimates in a separate script
% Script name: cold_pool_age_estimates.m
save('mixing_fractions_vars.mat','mf_t_cp','mf_fee','mf_fen','mf_fss','mf_fee_bg','mf_fen_bg','mf_fss_bg','mf_dD_cp','slope','Xeval','Yeval');

%% Calculating linear fit in mixing fraction space
    X = fw2(~isnan(fw2));
    Y = fw1(~isnan(fw1));
    N = 1; % linear polynomial
    [P,S,MU] = polyfit(X,Y,N);
    Xeval = 0:0.1:1;
    cloud_fit = polyval(P,Xeval,S,MU);
    slope(i) = (cloud_fit(end)-cloud_fit(1))./(Xeval(end)-Xeval(1));
    cloud_fit2 = slope(i).*Xeval + cloud_fit(1);

%% Plots for conserved properties air mixing fractions %%
yellow = [0.9290 0.6940 0.1250];
orange = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
% Using fen,fss,fdd from conserved properties
% PLOT #1
%   Same axes as water mixing fractions
    figure; hold on;
    scatter(fss,fdd,5,[.5 .5 .5],'filled')
%   Highlighting cold pools
    scatter(fss(iso_cold_pool_flag==1),fdd(iso_cold_pool_flag==1),5,'k','filled')
    xlabel('surface fraction ~ q')
    ylabel('downdraft fraction ~\theta')
    plot([0 1],[0 0],'-k');
    plot([0 0],[0 1],'-k');
    plot([1 0],[0 1],'-k');
    scatter(0,1,55,blue,'filled')% fdd=1; fss = 0 and fen = 0
    scatter(1,0,55,yellow,'filled')% fss=1; fdd = 0 and fen = 0
    scatter(0,0,55,orange,'filled')% fen=1; fdd = 0 and fss = 0
    ylim([-.2 1])
    xlim([-.2 1])
    axis square
    box on
    set(findall(gcf,'-property','Fontsize'),'FontSize',16)
    legend('non-cold pools','cold pools','FontSize',12);legend boxoff 
% PLOT #2
%   Same axes as conserved properties plot
    figure; hold on;
    scatter(fdd,fss,5,[.5 .5 .5],'filled')
%   Highlighting cold pools
    scatter(fdd(iso_cold_pool_flag==1),fss(iso_cold_pool_flag==1),5,'k','filled')
    ylabel('surface fraction ~ q')
    xlabel('downdraft fraction ~\theta')
    plot([0 1],[0 0],'-k');
    plot([0 0],[0 1],'-k');
    plot([1 0],[0 1],'-k');
    scatter(1,0,55,blue,'filled')% fdd=1; fss = 0 and fen = 0
    scatter(0,1,55,yellow,'filled')% fss=1; fdd = 0 and fen = 0
    scatter(0,0,55,orange,'filled')% fen=1; fdd = 0 and fss = 0
    ylim([-.2 1])
    xlim([-.2 1])
    axis square
    box on
    set(findall(gcf,'-property','Fontsize'),'FontSize',16)
    legend('non-cold pools','cold pools','Location','northwest','FontSize',12); legend boxoff 
    set(gca, 'XDir','reverse')
%% Plots for water mixing fractions %%
%  Using mass
%     figure; hold on;
%     set(gcf, 'Position', [1 1 750 650]);
    % scatter(fa2,fa1,15,[.5 .5 .5],'filled')
    scatter(fw2,fw1,30,'k');
    % scatter(fw2,fw1,25,dD_cp(ind),'filled');
    % scatter(fw2,fw1,25,dD_cp(ind)./q_cp(ind),'filled'); % for 2 member plots 
    scatter(fw2,fw1,25,dD_cp./q_cp,'filled'); % for 3 member plots 
    colormap(flip(jet(15)))
    han = colorbar;
    % han.YLabel.String = 'q*\deltaD';% [',char(8240),']'];
    han.YLabel.String = ['\deltaD [',char(8240),']'];
    hold on;
    plot([0 1],[0 0],'-k');
    hold on;
    plot([0 0],[0 1],'-k');
    hold on;
    plot([1 0],[0 1],'-k');
    % set(gca,'Xdir','reverse')
    plot(Xeval,cloud_fit2,'-','LineWidth',1)
    text(0.45,2*10^(-4),['slope = ',num2str(round(slope(i),2)),''])
%     xlim([min(fwd2) max(fwd2)])
%     ylim([min(fwd1) max(fwd1)])
    ylim([-.5e-4 3.5e-4])
    xlabel('surface fraction')
    ylabel('evaporee fraction')
    text(0.5,3*10^(-4),['CP#',num2str(i),'']) % for 3 member plots
    % text(0.9999,2*10^(-4),['CP#',num2str(i),'']) % for 2 member plots
    set(findall(gcf,'-property','Fontsize'),'FontSize',18)
    box on
    title('BL water mixture fraction')
    axis square
    caxis([min(dD_cp./q_cp) max(dD_cp./q_cp)])

    % Saving figure in different formats
%     saveas(gcf,['mixing_fractions_qdD_vs_q_CP#',num2str(i),'_water.png'])
%     saveas(gcf,['mixing_fractions_qdD_vs_q_CP#',num2str(i),'_water.fig'])

%%
% end
%% Plots: slope vs mixing fraction
slope(slope==9999) = NaN;
Error_slope(Error_slope==9999) = NaN;
% int(int==9999) = NaN;
max_fen(max_fen==0) = NaN;
max_fss(max_fss==0) = NaN;
max_fee(max_fee==0) = NaN;
min_fen(min_fen==0) = NaN;
min_fss(min_fss==0) = NaN;
min_fee(min_fee==0) = NaN;

%% slope vs max(mixing fraction)
% mfr = mixing fraction
% max(mfr) = min(1-mfr) ??? ;
% For fen
    figure; hold on;
    scatter(max_fen,slope,25,Error_slope,'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('max(f_e_n)')
    ylabel('slope in q-q*\deltaD space')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([0.5 1],[.2 .2],'--k')

% For fss
    figure; hold on;
    scatter(max_fss,slope,25,Error_slope,'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('max(f_s_s)')
    ylabel('slope in q-q*\deltaD space')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([0.2 1],[.2 .2],'--k')

% For fee
    figure; hold on;
    scatter(max_fee,slope,25,Error_slope,'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('max(f_e_e)')
    ylabel('slope in q-q*\deltaD space')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([0 .4],[.2 .2],'--k')
    
%% slope vs min(mixing fraction)
% mfr = mixing fraction
% min(mfr) = max(1-mfr) ??? ;
% For fen
    figure; hold on;
    scatter(min_fen,slope,25,Error_slope,'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('min(f_e_n)')
    ylabel('slope in q-q*\deltaD space')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([-.2 .8],[.2 .2],'--k')

% For fss
    figure; hold on;
    scatter(min_fss,slope,25,Error_slope,'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('min(f_s_s)')
    ylabel('slope in q-q*\deltaD space')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([0 .5],[.2 .2],'--k')

% For fee
    figure; hold on;
    scatter(min_fee,slope,25,Error_slope,'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('min(f_e_e)')
    ylabel('slope in q-q*\deltaD space')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([-.1 0],[.2 .2],'--k')

%% Plots for max temp vs max dD
% At time of min temp
Y_max = dD(t_min_ind(1:16));  % dD at time of min temp
X_max = T_min(1:16);          % COLDEST TEMP
Z_max = rrf(t_min_ind(1:16)); % rain rate at time of min temp

% During entire cold pool
% Y_max = max(dDcold,[],2); % peak dD
% X_max = T_min(1:16);      % COLDEST TEMP
% for jj = 1:length(b)
%     Z_max(jj) = max(rrf(max_ind(jj):end_ind(jj)));  % peak rain rate
% end
% Z_max = % rain rate at time of peak dD

figure; hold on;
plot(X_max,Y_max,'ok','MarkerFaceColor','w')   % at time of min temp
% plot(X_max,Y_max,'dk','MarkerFaceColor','k')   % peak dD
scatter(X_max,Y_max,25,Z_max,'filled')    % for RainRate
h = colorbar;
cmap = colormap(jet(9));
colormap(flip(cmap(1:8,:)))
xlabel('T_m_i_n')
ylabel('\deltaD')
ylabel(h,'RR [mm/hr]','Rotation',270)
axis square
caxis([0 8])
% ylim([0 0.8])

%% Plots for max change in temp vs max change in dD
% At time of min temp
Y_max = dD(t_min_ind(1:16))-dD(t_max_ind(1:16));  % dD anomaly at time of min temp
X_max = T_min(1:16)-Taf(t_max_ind(1:16));         % COLDEST TEMP anomaly
Z_max = t_max(1:16);                              % cold pool date (chronological order)

% During entire cold pool
Y_max2 = max(dDcold,[],2) - dD(t_max_ind(1:16));   % peak dD anomaly
X_max2 = T_min(1:16)-Taf(t_max_ind(1:16));         % COLDEST TEMP anomaly

figure; hold on;
plot(X_max,Y_max,'ok','MarkerFaceColor','w')     % at time of min temp
plot(X_max2,Y_max2,'dk','MarkerFaceColor','k')   % peak dD
scatter(X_max,Y_max,25,Z_max,'filled')           % for date
scatter(X_max2,Y_max2,25,Z_max,'filled')         % for date
h = colorbar;
cmap = colormap(jet(14));
xlabel('\DeltaT_m_i_n')
ylabel('\Delta\deltaD')
ylabel(h,'cold pool start time','Rotation',270)
axis square
h.TickLabels = datestr(h.Ticks,'dd/mmm');
caxis([datenum('01/28/2020','mm/dd/yyyy') datenum('02/11/2020','mm/dd/yyyy')])