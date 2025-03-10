% Code runs with output from "q_dD_FRAMEWORK_COLD_POOLS.m" script
%% Plots in 1/q-dD space (straight mixing lines)
for i = [2:11,13:16] %[9,14:15] %[8,10,11,13,16] %1:length(max_ind) %[2:7,9,14:15] %1:length(max_ind)-1 %
    % Calculating slope of best fit line for the 60-min data 
    slope = polyfit(1./qcoldnan(i,1:61),dDcoldnan(i,1:61),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,1./qcoldnan(i,1:61));
figure; 
set(gcf, 'Position', get(0, 'Screensize'));
j = plot(1./qcoldf(i,:),dDcoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(1./qcoldw(i,:),dDcoldw(i,:),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 60 minutes) in gray dots
jj = scatter(1./qcoldnan(i,1:61),dDcoldnan(i,1:61),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(1./qcoldnan(i,1:61),y_est,'r-','LineWidth',2)
% Cold pool data (up to 60 minutes) in time-colored dots
plot(1./qcoldf(i,:),dDcoldf(i,:),'ok','MarkerSize',9);%,'MarkerFaceColor','b')
jjj = scatter(1./qcold(i,1:61),dDcold(i,1:61),70,(tcold(i,1:61)-tcold(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(1./qcoldw(i,:),dDcoldw(i,:),'sk','MarkerSize',11);%,'MarkerFaceColor','c')

% Adding mixing lines and end member points (entraiment & surface)
Xfit = linspace(qent_cold(i),qsurf_cold(i),10);
Ylinearfit = (m2(i).*Xfit) + b2(i); % in units of q*dD
Ycurvefit = (Ylinearfit)./Xfit; % in units of dD
plot(1./qent_cold(i),dD_ent(i),'ok','MarkerFaceColor',[0.2, 0, 0],'MarkerSize',10)
j3 = plot(1./Xfit,Ycurvefit,'--k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(1./qsurf_cold(i),dDsurf_cold(i),'ok','MarkerFaceColor',[.5 .5 .5],'MarkerSize',10)
clearvars Xfit Ycurvefit Ylinearfit
% Identifying t_max & t_min
plot(1./qcoldf(i,1),dDcoldf(i,1),'^k','MarkerFaceColor','k','MarkerSize',12)
plot(1./qcoldw(i,1),dDcoldw(i,1),'^k','MarkerFaceColor','w','MarkerSize',12)

% Incorporating Rayleigh curves
% % for peak dD
% [q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0p(i),T0p(i),R0p(i));
% plot(q_rayleigh,q_rayleigh.*delt_rayleigh,'-b','LineWidth',.5)
% % for surface source
% [q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0(i),T0(i),R0(i));
% plot(q_rayleigh,q_rayleigh.*delt_rayleigh,'-m','LineWidth',.5)

% Adding mixing lines between Rayleigh liquid evap curve and peak dD
dD_prec = 15.6; % in permil
[q_lost,dD_lost,f_lost] = Rayleigh_liquid_evap(p0p(i),q0p(i),T0p(i),dD_prec); % p[Pa],q[g/kg],T[K],dD[permil]
% f = 0.75; % remaining fraction of liquid drop
cmap = b2rcolormap(length(f_lost)+1);
for kk = 1:length(f_lost)
    [m,b,b1] = mixing_line_slope_yint(q0p(i),dD0p(i),q_lost(kk),dD_lost(kk));
    Xfit = linspace(q_lost(kk),qsurf_cold(i),500);
    Ylinearfit = (m.*Xfit) + b; % in units of q*dD
    Ycurvefit = (Ylinearfit)./Xfit; % in units of dD
    j4 = plot(1./Xfit,Ycurvefit,'--','Color',cmap(kk,:),'LineWidth',2);
    j4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    clearvars m b b1 Xfit Ycurvefit Ylinearfit j4
end

% Colorbar properties
colormap(jet(12))
han = colorbar;
han.Title.String = "minutes since cp onset";
caxis([0 60])

% Figure properties
xlim([1/23 1/8]); % xlim([8 23]); 
ylim([-84 -64])
xlabel('1/q [kg/g]')
ylabel(['\deltaD [',char(8240),']'])
box on
grid on
set(findall(gcf,'-property','Fontsize'),'FontSize',20)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
legend('background','front','wake','entrainment','surface','t_0','t_m_i_n','Location','southeast')
    %...,'Rayleigh peak \deltaD','Rayleigh surface','Location','southwest')

% Including Ta,q,dD timeseries in small plot
% inset_plot(i,tcold,Taf,cp_matrix,qairf,dDf,rrf)
% Including liquid drop evaporation in another small plot
axes('Position',[.16 .16 .15 .15]); box on
plot(1./q_lost,dD_lost,'-r');
for kk = 1:length(f_lost)-1
    hold on;
    plot(1./q_lost(kk:kk+1),dD_lost(kk:kk+1),'-','Color',cmap(kk,:),'LineWidth',2);
end
xlabel('1/q_l_o_s_t');ylabel(['\deltaD_l_o_s_t [',char(8240),']'])
title(['CP #',num2str(i),'; onset on ',datestr(tcold(i,1))])

% Saving figure in different formats
saveas(gcf,['obs_dD_vs_q1_CP#',num2str(i),'_curve_fit_centroid.png'])
saveas(gcf,['obs_dD_vs_q1_CP#',num2str(i),'_curve_fit_centroid.fig'])

% figure; hold on
% for kk = 1:length(f_lost)
%     [m,b,b1] = mixing_line_slope_yint(q0p(i),dD0p(i),q_lost(kk),dD_lost(kk));
%     Xfit = linspace(q_lost(kk),qsurf_cold(i),2); % q
%     Ylinearfit = (m.*Xfit) + b; % q*dD
%     Ycurvefit = (Ylinearfit)./Xfit; % dD
%     plot(1./Xfit,Ycurvefit,'-o','Color',1-cmap(kk,:))%,'MarkerFaceColor',cmap(kk,:));
%     % 1-cmap => (darker color in middle of colorbar)
%     xlim([0 18*10^4]); % ylim([-68 8])
%     clearvars m b b1 Xfit Ycurvefit Ylinearfit
% end
% xlabel('1/q')
% ylabel(['\deltaD [',char(8240),']'])

clearvars y_est slope q_rayleigh delt_rayleigh q_lost dD_lost m b b1

end
