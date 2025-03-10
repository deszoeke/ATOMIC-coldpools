%% mixing fractions in isotopic space q-q*dD %%
% [fen, fss, fdd] = th_q_to_mixfraction(th_ob,q_ob, th_en,q_en, th_surf,q_surf, th_d,q_d);
% figure; hold on; % for all cps in single plot
% set(gcf, 'Position', [1 1 750 650]);

% Allocating "slope" variables
b=1:17;
slope = 9999*ones(size(b));
int = 9999*ones(size(b));
Error_slope = 9999*ones(size(b));

for i = [2:11,13:16] %[9,14:15] %[8,10,11,13,16] %1:length(max_ind) %[2:7,9,14:15] %1:length(max_ind)-1 %
    q_bg  = qcoldnan(i,1:61);
    dD_bg = q_bg.*dDcoldnan(i,1:61);
    q_cp  = qcold(i,1:61);
    dD_cp = q_cp.*dDcold(i,1:61);
    q_en  = qent_cold(i)*ones(1,61);
    dD_en = q_en.*(dD_ent(i)*ones(1,61));
    q_surf  = qsurf_cold(i)*ones(1,61);
    dD_surf = q_surf.*(dDsurf_cold(i)*ones(1,61));
    
    dD_prec = 10; % in permil
    [q_lost,dD_lost,f_lost] = Rayleigh_liquid_evap(p0p(i),q0p(i),T0p(i),dD_prec); % p[Pa],q[g/kg],T[K],dD[permil]
    f = 0.25; % remaining fraction of liquid drop
    q_ev  = q_lost(f_lost==f)*ones(1,61);
    dD_ev = q_ev.*(dD_lost(f_lost==f)*ones(1,61));
    [fen, fss, fee] = x_y_to_mixfraction(dD_cp,q_cp, dD_en,q_en, dD_surf,q_surf, dD_ev,q_ev);
    [fen_bg, fss_bg, fee_bg] = x_y_to_mixfraction(dD_bg,q_bg, dD_en,q_en, dD_surf,q_surf, dD_ev,q_ev);

    % Calculating the mass of water going into the mixing fractions
    % This methodology only applies to when fen == 0, correct?
    % Otherwise the linear relationship does not hold anymore, right?
    ind = find(fen>0 & fen<0.1); % for when fen == 0
    dummy = 1:length(fen); % for the entire dataset, for plotting limits
    % [fw1,fw2] = mass_mixing_test(ind,fev,fss,q_ev,q_surf);
    [fwd1,fwd2] = mass_mixing_test(dummy,fev,fss,q_ev,q_surf);
    [fw1,fw2,fw3,q] = mass_mixing(dummy,fev,fss,fen,q_ev,q_surf,q_en);
    %     fa1 = fev(ind);
    %     fa2 = fss(ind);
    %     q1 = q_ev(ind);
    %     q2 = q_surf(ind);
    %     q = fa1.*q1 + fa2.*q2;
    %     fw1 = (fa1.*q1)./q;
    %     fw2 = (fa2.*q2)./q;

    % Calculating linear fit in mixing fraction space
    X = fss(~isnan(fss));
    Y = fee(~isnan(fee));
%     N = 1; % linear polynomial
    mdl = fitlm(X,Y);
    Error_slope(i) = mdl.Rsquared.Adjusted; % R-squared
    slope(i) = mdl.Coefficients.Estimate(2);
    max_fen(i) = max(fen,[],'omitnan');
    max_fss(i) = max(fss,[],'omitnan');    
    max_fee(i) = max(fee,[],'omitnan');
    min_fen(i) = min(fen,[],'omitnan');
    min_fss(i) = min(fss,[],'omitnan');    
    min_fee(i) = min(fee,[],'omitnan');
end

%% Plots for air mixing fractions %%
%     % Using specific humidity
%     % set(gcf, 'Position', [1 1 900 500]);
%     figure; hold on;
%     j = scatter(fss_bg,fev_bg,15,[.5 .5 .5],'filled');
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     j = scatter(fss,fev,30,'k');
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     % scatter(fss,fev,25,dD_cp,'filled');
%     j = scatter(fss,fev,25,dD_cp./q_cp,'filled');
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     colormap(flip(jet(15)))
%     han = colorbar;
%     % han.YLabel.String = 'q*\deltaD';% [',char(8240),']'];
%     han.YLabel.String = ['\deltaD [',char(8240),']'];
%     j = plot([0 1],[0 0],'-k');
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     j = plot([0 0],[0 1],'-k');
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     j = plot([1 0],[0 1],'-k');
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     % set(gca,'Xdir','reverse')
%     plot(Xeval,cloud_fit2,'-','LineWidth',1)
%     text(0.45,0.8,['slope = ',num2str(round(slope,2)),''])
%     legend('linear fit','Location','southeast')
%     xlim([-0.2 1])
%     ylim([-0.2 1])
%     xlabel('surface fraction')
%     ylabel('evaporatee fraction')
%     text(0.6,0.9,['CP#',num2str(i),''])
%     set(findall(gcf,'-property','Fontsize'),'FontSize',18)
%     box on
%     title('BL air mixture fraction')
%     axis square
%     caxis([min(dD_cp./q_cp) max(dD_cp./q_cp)])

    % Saving figure in different formats
%     saveas(gcf,['mixing_fractions_qdD_vs_q_CP#',num2str(i),'_wslope.png'])
%     saveas(gcf,['mixing_fractions_qdD_vs_q_CP#',num2str(i),'_wslope.fig'])
    
%% Calculating linear fit in mixing fraction space
    X = fw2(~isnan(fw2));
    Y = fw1(~isnan(fw1));
    N = 1; % linear polynomial
    [P,S,MU] = polyfit(X,Y,N);
    Xeval = 0:0.1:1;
    cloud_fit = polyval(P,Xeval,S,MU);
    slope = (cloud_fit(end)-cloud_fit(1))./(Xeval(end)-Xeval(1));
    cloud_fit2 = slope.*Xeval + cloud_fit(1);

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
    text(0.45,2*10^(-4),['slope = ',num2str(round(slope,2)),''])
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

% end

%% For fee (color-coded by CP number)
for k = [4,11,13]
    figure; hold on;
    scatter(max_fee(k),slope(k),25,k,'filled')
    h = colorbar;
    cmap = colormap(jet(14));
    xlabel('max(f_h_y'')')
    ylabel('slope (in f_h_y''-f_s_S'')')
    title(['CP ',num2str(k)])
%     ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([2 16])
    plot([0 .4],[.2 .2],'--k')
    ylim([-.2 .8])
end

%% For fss (color-coded by R_squared)
    figure; hold on;
    scatter(max_fss([2:11,13:16]),slope([2:11,13:16]),25,Error_slope([2:11,13:16]),'filled')
    h = colorbar;
    cmap = colormap(jet(20));
    xlabel('max(f_s_s)')
    ylabel('slope (in f_h_y''-f_s_S'')')
    ylabel(h,'R^2','Rotation',270)
    axis square
    caxis([0 1])
    plot([0 .8],[.2 .2],'--k')