function ax = vert_axes_stack(n)
% n = 6; % number of axes

margins = 0.15; % top and bottom margins as fraction of figure
gap = 0.02; % gap between axes
height = (1 - 2 * margins - (n - 1) * gap) / n;

for k = 1:n
    bottom = 1 - margins - k * height - (k - 1) * gap;
    ax(k) = axes('Position', [0.16, bottom, 0.64, height], 'fontsize',18);
    
    % Example plot for each axis
    % plot(rand(10,1));
    
    % Optional: remove x-axis labels except bottom
    if k < n
        ax(k).XTickLabel = [];
    end
end

end
