% plot_W_iso_ages.m

lambda_mass_flux; % -> lambda_W
regress_transform_cp;

rankdD_cronOrder = [7, 13, 10, 5, 3, 4, 1, 6, 8, 12, 14, 9, 11, 2];

clf()
subplot(3,3,1, 'fontsize',12)
hold on
j = 0; % chronological index
n = zeros(14,1);
for i = [2:11, 13:16]
    j = j + 1;
    tdd = t_age_vec{i}(1);
    p = proj{i};
    [~, imx] = max(p);
    ii = imx:find(isfinite(p) & p>=0, 1, 'last');
    % p(p < 0) = NaN;
    x = -log( p(ii) );
    y = lambda_W(i) * 60 * (tdd + (0:length(x)-1)');
    n(j) = sum(isfinite(p(ii)));
    if n(j)>1
        plot(x,y, '.-')
    else
        plot(x,y,'.')
    end
    text(x(end)+0.1, y(end), sprintf("%i", rankdD_cronOrder(j)))
end
axis tight
ylim([-0.04, 0.4])
xlim([-0.8, 8])
xlabel({'nondimensional','isotope age -log(g''/g''_{dd})'})
ylabel({'nondimensional','water vapor age \lambda_W(t-t_{dd})'})
axis square

% saveas(gcf, 'W_iso_ages.eps', 'epsc')
% saveas(gcf, 'W_iso_ages.svg')
% saveas(gcf, 'W_iso_ages.png')

% saveas(gcf, 'W_iso_ages2.eps', 'epsc')
% saveas(gcf, 'W_iso_ages2.svg')
% saveas(gcf, 'W_iso_ages2.png')

% CP 3 and 9 have no data