function ax = vert_axes_stack(n)
% n = 6; % number of axes

margins = 0.08; % top and bottom margins as fraction of figure
gap = 0.02; % gap between axes
height = (1 - 2 * margins - (n - 1) * gap) / n;

for k = 1:n
    bottom = 1 - margins - k * height - (k - 1) * gap;
    ax(k) = axes('Position', [0.1, bottom, 0.8, height], 'fontsize',14);
    
    % Example plot for each axis
    % plot(rand(10,1));
    
    % Optional: remove x-axis labels except bottom
    if k < n
        ax(k).XTickLabel = [];
    end
end

end