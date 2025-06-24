function plot_time_comp(ax, t_comp2, y, color)
    me = mean(y, 'omitnan');
    se = std(y, 'omitnan') / sqrt(size(y,1));

    plot(ax, t_comp2, me, '-', 'linewidth',1.7, 'color',color);
    hold(ax, 'on');
    plot(ax, t_comp2, me + se, '-', 'linewidth',0.7, 'color',color);
    plot(ax, t_comp2, me - se, '-', 'linewidth',0.7, 'color',color);
    plot(ax, t_comp2, median(y, 'omitnan'), '.', 'color',color);
end