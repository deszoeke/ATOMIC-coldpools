function right_yaxis_plotting(t_comp2,var,var_17,mask_strong100,mask_weak100,n,Tmax_i)
    blue = [0 0.4470 0.7410];
    anomaly = var_17-var_17(:,Tmax_i); 
    anomaly = nanmean(anomaly)+nanmean(var_17(:,Tmax_i)); % anomaly centered in actual value range (not centered in ZERO)
    anomaly_strong = var(mask_strong100,:)-var(mask_strong100,Tmax_i);
    anomaly_strong = nanmean(anomaly_strong)+nanmean(var_17(:,Tmax_i)); % anomaly centered in actual value range (not centered in ZERO)
    anomaly_weak = var(mask_weak100,:)-var(mask_weak100,Tmax_i);
    anomaly_weak = nanmean(anomaly_weak)+nanmean(var_17(:,Tmax_i)); % anomaly centered in actual value range (not centered in ZERO)

    % yyaxis right
    plot(t_comp2,anomaly,'-k','LineWidth',2)
    hold on;
    j = plot(t_comp2,anomaly+std(anomaly,'omitnan')/sqrt(17),'.k','MarkerSize',.2);
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
    jj = plot(t_comp2,anomaly-std(anomaly,'omitnan')/sqrt(17),'.k','MarkerSize',.2);  
    jj.Annotation.LegendInformation.IconDisplayStyle = 'off';

    plot(t_comp2,anomaly_strong,'-r','LineWidth',2)
    j3 = plot(t_comp2,anomaly_strong+std(anomaly_strong,'omitnan')/sqrt(n),'.r','MarkerSize',.2);
    j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
    j4 = plot(t_comp2,anomaly_strong-std(anomaly_strong,'omitnan')/sqrt(n),'.r','MarkerSize',.2);
    j4.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(t_comp2,anomaly_weak,'-','LineWidth',2,'Color',blue)
    plot(t_comp2,anomaly_weak+std(anomaly_weak,'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)
    plot(t_comp2,anomaly_weak-std(anomaly_weak,'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)       
    
    plot([-18.7 -18.7],[-300 100],':k','LineWidth',.5)
    plot([-60 60],[0 0],':k','LineWidth',.5)
    plot([0 0],[-300 100],':k','LineWidth',.5)
    plot([30.4 30.4],[-300 100],':k','LineWidth',.5)

    % yyaxis left
    % set(gca,'yticklabels',[])
end

% function mask_lines_plotting(t_comp2,var,var_17,mask_strong100,mask_weak100,n,Tmax_i)
%     blue = [0 0.4470 0.7410];
%     anomaly = var_17-var_17(:,Tmax_i);
%     anomaly_strong = var(mask_strong100,:)-var(mask_strong100,Tmax_i);
%     anomaly_weak = var(mask_weak100,:)-var(mask_weak100,Tmax_i);
% 
%     plot(t_comp2,nanmean(anomaly),'-k','LineWidth',2)
%     hold on;
%     j = plot(t_comp2,nanmean(anomaly)+std(anomaly,'omitnan')/sqrt(17),'.k','MarkerSize',.2);
%     j.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     jj = plot(t_comp2,nanmean(anomaly)-std(anomaly,'omitnan')/sqrt(17),'.k','MarkerSize',.2);  
%     jj.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
%     plot(t_comp2,nanmean(anomaly_strong),'-r','LineWidth',2)
%     j3 = plot(t_comp2,nanmean(anomaly_strong)+std(anomaly_strong,'omitnan')/sqrt(n),'.r','MarkerSize',.2);
%     j3.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     j4 = plot(t_comp2,nanmean(anomaly_strong)-std(anomaly_strong,'omitnan')/sqrt(n),'.r','MarkerSize',.2);
%     j4.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     plot(t_comp2,nanmean(anomaly_weak),'-','LineWidth',2,'Color',blue)
%     plot(t_comp2,nanmean(anomaly_weak)+std(anomaly_weak,'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)
%     plot(t_comp2,nanmean(anomaly_weak)-std(anomaly_weak,'omitnan')/sqrt(n),'.','Color',blue,'MarkerSize',.2)       
% end