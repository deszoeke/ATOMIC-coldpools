function conserved_plot_iso(iso,th_ob,q_ob,th_en,q_en,th_d,q_d,th_surf,q_surf)
    figure;
    scatter(th_ob,q_ob,15,'k','s')
    hold on;
    scatter(th_en,q_en,10,[0.8500 0.3250 0.0980],'filled')
    hold on;
    scatter(th_d,q_d,10,[0 0.4470 0.7410],'filled')
    hold on;
    scatter(th_surf,q_surf,10,[0.9290 0.6940 0.1250],'filled')
    hold on;
    scatter(th_ob,q_ob,13,iso,'s','filled')
    colormap(flip(jet(18)))
    colorbar;
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
end
