% This code runs with output variables from
% cold_pool_detection_algorithm_recovery_modified.m
% Last Modified May 24, 2022

orange = [0.8500 0.3250 0.0980];
blue = [0 0.4470 0.7410];
mustard = [0.9290 0.6940 0.1250];
green = [0.4196 0.5569 0.1373];

% Plotting the entire cold pool matrix (99 in total)
Taf_comp2 = Taf_comp2(2:end,:);
wspd_comp2 = wspd_comp2(2:end,:);
qair_comp2 = qair_comp2(2:end,:);
prec_comp2 = prec_comp2(2:end,:);
rh_comp2 = rh_comp2(2:end,:);
u_comp2 = u_comp2(2:end,:);
v_comp2 = v_comp2(2:end,:);
sh_comp2 = sh_comp2(2:end,:);
lh_comp2 = lh_comp2(2:end,:);

%% Plotting the entire cold pool matrix (99 in total)
figure;
subplot(721)
    plot(t_comp2,nanmean(Taf_comp2),'-','LineWidth',2,'Color',orange)
    hold on;
    plot(t_comp2,nanmean(Taf_comp2)+std(Taf_comp2,'omitnan')/sqrt(99),'-','Color',orange)
    plot(t_comp2,median(Taf_comp2,'omitnan'),'-o','Color',orange)
    plot(t_comp2,nanmean(Taf_comp2)-std(Taf_comp2,'omitnan')/sqrt(99),'-','Color',orange)
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    ylabel('T_a [\circC]')
    title('99 cold pools composite')
    xlim([-60 60])
    ylim([24.3 26.3])
subplot(723)
    plot(t_comp2,nanmean(wspd_comp2),'Color',green,'LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(wspd_comp2)+std(wspd_comp2,'omitnan')/sqrt(99),'Color',green)
    plot(t_comp2,median(wspd_comp2,'omitnan'),'-o','Color',green)
    plot(t_comp2,nanmean(wspd_comp2)-std(wspd_comp2,'omitnan')/sqrt(99),'Color',green)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('U [m/s]')
    box off
    xlim([-60 60])
    ylim([8 11.3])
subplot(725)
%     plot(t_comp2,nanmean(q_iso_comp2),'-sk')
    hold on;
    plot(t_comp2,nanmean(qair_comp2),'-','LineWidth',2,'Color',blue)
    plot(t_comp2,nanmean(qair_comp2)+std(qair_comp2,'omitnan')/sqrt(99),'-','Color',blue)
    plot(t_comp2,median(qair_comp2,'omitnan'),'-o','Color',blue)
    plot(t_comp2,nanmean(qair_comp2)-std(qair_comp2,'omitnan')/sqrt(99),'-','Color',blue)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('q [g/kg]')
    box off
    xlim([-60 60])
    ylim([14.2 15.4])
%     legend('Picarro','PSD')
subplot(727)
    plot(t_comp2,nanmean(prec_comp2),'-r','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(prec_comp2)+std(prec_comp2,'omitnan')/sqrt(99),'-r')
    plot(t_comp2,median(prec_comp2,'omitnan'),'-or')
    plot(t_comp2,nanmean(prec_comp2)-std(prec_comp2,'omitnan')/sqrt(99),'-r')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    box off
    ylabel('rain [mm/hr]')
    xlim([-60 60])
    ylim([0 3.5])
subplot(729)
    plot(t_comp2,nanmean(rh_comp2),'-b','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(rh_comp2)+std(rh_comp2,'omitnan')/sqrt(99),'-b')
    plot(t_comp2,median(rh_comp2,'omitnan'),'-ob')
    plot(t_comp2,nanmean(rh_comp2)-std(rh_comp2,'omitnan')/sqrt(99),'-b')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('RH [%]')
    box off
    xlim([-60 60])
    ylim([68 80])
% Including wind components %
    subplot(7,2,11)
        plot(t_comp2,nanmean(u_comp2),'-c','LineWidth',2)
        hold on;
        plot(t_comp2,nanmean(u_comp2)+std(u_comp2,'omitnan')/sqrt(99),'-c')
        plot(t_comp2,median(u_comp2,'omitnan'),'-oc')
        plot(t_comp2,nanmean(u_comp2)-std(u_comp2,'omitnan')/sqrt(99),'-c')    
        plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
        plot([0 0],[-100 100],':k','LineWidth',1.5)
        plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
        ylabel('E-W wind [m/s]')
        box off
        xlim([-60 60])
        ylim([-10.5 -7])
    subplot(7,2,13)
        plot(t_comp2,nanmean(v_comp2),'-g','LineWidth',2)
        hold on;
        plot(t_comp2,nanmean(v_comp2)+std(v_comp2,'omitnan')/sqrt(99),'-g')
        plot(t_comp2,median(v_comp2,'omitnan'),'-og')
        plot(t_comp2,nanmean(v_comp2)-std(v_comp2,'omitnan')/sqrt(99),'-g')    
        plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
        plot([0 0],[-100 100],':k','LineWidth',1.5)
        plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
        ylabel('N-S wind [m/s]')
        xlabel('time [minutes]')
        xlim([-60 60])
        ylim([-5.2 -1.6])
        box off
% Including surface fluxes %
% subplot(7,2,11)
%     plot(t_comp2,nanmean(sh_comp2),'-c','LineWidth',2)
%     hold on;
%     plot(t_comp2,nanmean(sh_comp2)+std(sh_comp2,'omitnan')/sqrt(99),'-c')
%     plot(t_comp2,median(sh_comp2,'omitnan'),'-oc')
%     plot(t_comp2,nanmean(sh_comp2)-std(sh_comp2,'omitnan')/sqrt(99),'-c')    
%     plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
%     plot([0 0],[-100 100],':k','LineWidth',1.5)
%     plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
%     ylabel('SH [Wm^-^2]')
%     box off
%     xlim([-60 60])
%     ylim([-23 0])
% subplot(7,2,13)
%     plot(t_comp2,nanmean(lh_comp2),'-g','LineWidth',2)
%     hold on;
%     plot(t_comp2,nanmean(lh_comp2)+std(lh_comp2,'omitnan')/sqrt(99),'-g')
%     plot(t_comp2,median(lh_comp2,'omitnan'),'-og')
%     plot(t_comp2,nanmean(lh_comp2)-std(lh_comp2,'omitnan')/sqrt(99),'-g')    
%     plot([-19.7 -19.7],[-230 0],':k','LineWidth',1.5)
%     plot([0 0],[-230 0],':k','LineWidth',1.5)
%     plot([31.5 31.5],[-230 0],':k','LineWidth',1.5)
%     ylabel('LH [Wm^-^2]')
%     xlabel('time [minutes]')
%     xlim([-60 60])
%     ylim([-225 -180])
%     box off
    set(findall(gcf,'-property','Fontsize'),'FontSize',14)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])  
    
%% Plotting the 16 cold pools that overlap with isotopic data
% already subset to 15!
% Taf_comp2 = Taf_comp2([65,67,69:78,82:84,86]-1,:);
% wspd_comp2 = wspd_comp2([65,67,69:78,82:84,86]-1,:);
% qair_comp2 = qair_comp2([65,67,69:78,82:84,86]-1,:);
% prec_comp2 = prec_comp2([65,67,69:78,82:84,86]-1,:);
% rh_comp2 = rh_comp2([65,67,69:78,82:84,86]-1,:);
% u_comp2 = u_comp2([65,67,69:78,82:84,86]-1,:);
% v_comp2 = v_comp2([65,67,69:78,82:84,86]-1,:);
% sh_comp2 = sh_comp2([65,67,69:78,82:84,86]-1,:);
% lh_comp2 = lh_comp2([65,67,69:78,82:84,86]-1,:);

%%
subplot(722)
    plot(t_comp2,nanmean(Taf_comp2),'-','LineWidth',2,'Color',orange)
    hold on;
    plot(t_comp2,nanmean(Taf_comp2)+std(Taf_comp2,'omitnan')/sqrt(17),'-','Color',orange)
    plot(t_comp2,median(Taf_comp2,'omitnan'),'-o','Color',orange)
    plot(t_comp2,nanmean(Taf_comp2)-std(Taf_comp2,'omitnan')/sqrt(17),'-','Color',orange)
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('DXS [^{\fontsize{10}o}/{\fontsize{10}oo}]')
    box off
    ylabel('T_a [\circC]')
    title('16 cold pools composite')
    xlim([-60 60])
    ylim([24.3 26.3])
subplot(724)
    plot(t_comp2,nanmean(wspd_comp2),'Color',green,'LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(wspd_comp2)+std(wspd_comp2,'omitnan')/sqrt(17),'Color',green)
    plot(t_comp2,median(wspd_comp2,'omitnan'),'-o','Color',green)
    plot(t_comp2,nanmean(wspd_comp2)-std(wspd_comp2,'omitnan')/sqrt(17),'Color',green)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('U [m/s]')
    box off
    xlim([-60 60])
    ylim([8 11.3])
subplot(726)
%     plot(t_comp2,nanmean(q_iso_comp2),'-sk')
    hold on;
    plot(t_comp2,nanmean(qair_comp2),'-','LineWidth',2,'Color',blue)
    plot(t_comp2,nanmean(qair_comp2)+std(qair_comp2,'omitnan')/sqrt(17),'-','Color',blue)
    plot(t_comp2,median(qair_comp2,'omitnan'),'-o','Color',blue)
    plot(t_comp2,nanmean(qair_comp2)-std(qair_comp2,'omitnan')/sqrt(17),'-','Color',blue)    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('q [g/kg]')
    box off
    xlim([-60 60])
    ylim([14.2 15.4])
%     legend('Picarro','PSD')
subplot(728)
    plot(t_comp2,nanmean(prec_comp2),'-r','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(prec_comp2)+std(prec_comp2,'omitnan')/sqrt(17),'-r')
    plot(t_comp2,median(prec_comp2,'omitnan'),'-or')
    plot(t_comp2,nanmean(prec_comp2)-std(prec_comp2,'omitnan')/sqrt(17),'-r')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    box off
    ylabel('rain [mm/hr]')
    xlim([-60 60])
    ylim([0 3.5])
subplot(7,2,10)
    plot(t_comp2,nanmean(rh_comp2),'-b','LineWidth',2)
    hold on;
    plot(t_comp2,nanmean(rh_comp2)+std(rh_comp2,'omitnan')/sqrt(17),'-b')
    plot(t_comp2,median(rh_comp2,'omitnan'),'-ob')
    plot(t_comp2,nanmean(rh_comp2)-std(rh_comp2,'omitnan')/sqrt(17),'-b')    
    plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
    plot([0 0],[-100 100],':k','LineWidth',1.5)
    plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
    ylabel('RH [%]')
    box off
    xlim([-60 60])
    ylim([68 80])
% Including wind components %
    subplot(7,2,12)
        plot(t_comp2,nanmean(u_comp2),'-c','LineWidth',2)
        hold on;
        plot(t_comp2,nanmean(u_comp2)+std(u_comp2,'omitnan')/sqrt(17),'-c')
        plot(t_comp2,median(u_comp2,'omitnan'),'-oc')
        plot(t_comp2,nanmean(u_comp2)-std(u_comp2,'omitnan')/sqrt(17),'-c')    
        plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
        plot([0 0],[-100 100],':k','LineWidth',1.5)
        plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
        ylabel('E-W wind [m/s]')
        box off
        xlim([-60 60])
        ylim([-10.5 -7])
    subplot(7,2,14)
        plot(t_comp2,nanmean(v_comp2),'-g','LineWidth',2)
        hold on;
        plot(t_comp2,nanmean(v_comp2)+std(v_comp2,'omitnan')/sqrt(17),'-g')
        plot(t_comp2,median(v_comp2,'omitnan'),'-og')
        plot(t_comp2,nanmean(v_comp2)-std(v_comp2,'omitnan')/sqrt(17),'-g')    
        plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
        plot([0 0],[-100 100],':k','LineWidth',1.5)
        plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
        ylabel('N-S wind [m/s]')
        xlabel('time [minutes]')
        xlim([-60 60])
        ylim([-5.2 -1.6])
        box off
% Including surface fluxes %
% subplot(7,2,12)
%     plot(t_comp2,nanmean(sh_comp2),'-c','LineWidth',2)
%     hold on;
%     plot(t_comp2,nanmean(sh_comp2)+std(sh_comp2,'omitnan')/sqrt(99),'-c')
%     plot(t_comp2,median(sh_comp2,'omitnan'),'-oc')
%     plot(t_comp2,nanmean(sh_comp2)-std(sh_comp2,'omitnan')/sqrt(99),'-c')    
%     plot([-19.7 -19.7],[-100 100],':k','LineWidth',1.5)
%     plot([0 0],[-100 100],':k','LineWidth',1.5)
%     plot([31.5 31.5],[-100 100],':k','LineWidth',1.5)
%     ylabel('SH [Wm^-^2]')
%     box off
%     xlim([-60 60])
%     ylim([-23 0])
% subplot(7,2,14)
%     plot(t_comp2,nanmean(lh_comp2),'-g','LineWidth',2)
%     hold on;
%     plot(t_comp2,nanmean(lh_comp2)+std(lh_comp2,'omitnan')/sqrt(99),'-g')
%     plot(t_comp2,median(lh_comp2,'omitnan'),'-og')
%     plot(t_comp2,nanmean(lh_comp2)-std(lh_comp2,'omitnan')/sqrt(99),'-g')    
%     plot([-19.7 -19.7],[-230 0],':k','LineWidth',1.5)
%     plot([0 0],[-230 0],':k','LineWidth',1.5)
%     plot([31.5 31.5],[-230 0],':k','LineWidth',1.5)
%     ylabel('LH [Wm^-^2]')
%     xlabel('time [minutes]')
%     xlim([-60 60])
%     ylim([-225 -180])
%     box off
    set(findall(gcf,'-property','Fontsize'),'FontSize',14)
    set(findall(gcf,'-property','TickLength'),'TickLength',[.05,.1])