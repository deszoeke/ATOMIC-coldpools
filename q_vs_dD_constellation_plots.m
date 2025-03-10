% Code runs with output from "q_dD_FRAMEWORK_COLD_POOLS.m" script
%% Plots for dD [for plotting purposes]
for i = 8 %[2:11,13:16] %[9,14:15] %[8,10,11,13,16] %1:length(max_ind) %[2:7,9,14:15] %1:length(max_ind)-1 %
    % Calculating slope of best fit line for the 60-min data 
    slope = polyfit(qcoldnan(i,1:61),dDcoldnan(i,1:61),1);
    disp(['Equation is y = ' num2str(slope(1)) '*x + ' num2str(slope(2))])
    y_est = polyval(slope,qcoldnan(i,1:61));
figure; 
set(gcf, 'Position', get(0, 'Screensize'));
j = plot(qcoldf(i,:),dDcoldf(i,:),'--k');
j.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold on;
j2 = plot(qcoldw(i,1:61),dDcoldw(i,1:61),'--k');
j2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Background data (previous 60 minutes) in gray dots
scatter(qcoldnan(i,1:61),dDcoldnan(i,1:61),15,[.5 .5 .5],'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jj = plot(qcoldnan(i,1:61),y_est,'r-','LineWidth',2);
jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Cold pool data (up to 60 minutes) in time-colored dots
plot(qcoldf(i,:),dDcoldf(i,:),'ok','MarkerSize',9);%,'MarkerFaceColor','b')
jjj = scatter(qcold(i,1:61),dDcold(i,1:61),70,(tcold(i,1:61)-tcold(i,1))*1440,'filled'); %-datenum('01312020','mmddyyyy'),'filled')
jjj.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(qcoldw(i,1:61),dDcoldw(i,1:61),'sk','MarkerSize',11);%,'MarkerFaceColor','c')

% Adding mixing lines between end member points (entraiment & surface)
Xfit = linspace(qent_cold(i),qsurf_cold(i),10);
Ylinearfit = (m2(i).*Xfit) + b2(i);
Ycurvefit = (Ylinearfit)./Xfit;
j3 = plot(Xfit,Ycurvefit,'--k','LineWidth',2);
j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
clearvars Xfit Ycurvefit Ylinearfit
% Identifying t_max & t_min
plot(qcoldf(i,1),dDcoldf(i,1),'^k','MarkerFaceColor','k','MarkerSize',12)
plot(qcoldw(i,1),dDcoldw(i,1),'^k','MarkerFaceColor','w','MarkerSize',12)

% % Incorporating Rayleigh curves
% % for peak dD
% [q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0p(i),T0p(i),R0p(i));
% plot(q_rayleigh,delt_rayleigh,'-b','LineWidth',.5)
% % for surface source
% [q_rayleigh,delt_rayleigh] = Rayleigh_curve_cond(q0(i),T0(i),R0(i));
% plot(q_rayleigh,delt_rayleigh,'-m','LineWidth',.5)

% Adding markers for end member points (entraiment & surface)
plot(qent_cold(i),dD_ent(i),'ok','MarkerFaceColor',orange,'MarkerSize',10)
plot(qsurf_cold(i),dDsurf_cold(i),'ok','MarkerFaceColor',mustard,'MarkerSize',10)

% Colorbar properties
colormap(jet(12))
han = colorbar;
han.Label.String = "minutes since t_0";
caxis([0 60])

% Figure properties
xlim([8 23]); 
% ylim([-84 -64]); % ylim([-90 -64])
xlabel('specific humidity, q [g/kg]')
ylabel(['\deltaD [',char(8240),']'])
box on
grid on
set(findall(gcf,'-property','Fontsize'),'FontSize',28)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
legend('background','front','wake','t_0','t_m_i_n','entrainment','surface','Location','northwest')
title(['onset on ',datestr(tcold(i,1))])

% Saving figure in different formats
% saveas(gcf,['obs_dD_vs_q_CP#',num2str(i),'_ent_surf_sources.png'])
% saveas(gcf,['obs_dD_vs_q_CP#',num2str(i),'_ent_surf_sources.fig'])

end