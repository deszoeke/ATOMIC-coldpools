function mask_lines_plotting(t_comp2,var,var_17,mask_strong100,mask_weak100,n,Tmax_i)
    blue = [0 0.4470 0.7410];
    anomaly = var_17-var_17(:,Tmax_i);
    anomaly_strong = var(mask_strong100,:)-var(mask_strong100,Tmax_i);
    anomaly_weak = var(mask_weak100,:)-var(mask_weak100,Tmax_i);

    plot(t_comp2,nanmean(anomaly),'-k','LineWidth',2)
    hold on;
    j = plot(t_comp2,nanmean(anomaly)+std(anomaly,'omitnan')/sqrt(17),'.k','MarkerSize',.2);
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    jj = plot(t_comp2,nanmean(anomaly)-std(anomaly,'omitnan')/sqrt(17),'.k','MarkerSize',.2);  
    jj.Annotation.LegendInformation.IconDisplayStyle = 'off';

    plot(t_comp2,nanmean(anomaly_strong),'-r','LineWidth',2)
    j3 = plot(t_comp2,nanmean(anomaly_strong)+std(anomaly_strong,'omitnan')/sqrt(n),'.r','MarkerSize',.2);
    j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    j4 = plot(t_comp2,nanmean(anomaly_strong)-std(anomaly_strong,'omitnan')/sqrt(n),'.r','MarkerSize',.2);
    j4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(t_comp2,nanmean(anomaly_weak),'-','LineWidth',2,'Color',blue)
    plot(t_comp2,nanmean(anomaly_weak)+std(anomaly_weak,'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)
    plot(t_comp2,nanmean(anomaly_weak)-std(anomaly_weak,'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)       
end