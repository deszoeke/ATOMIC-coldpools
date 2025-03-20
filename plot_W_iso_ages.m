% plot_W_iso_ages.m

lambda_mass_flux; % -> lambda_W
regress_transform_cp;

clf()
subplot(2,2,1, 'fontsize',12)
hold on
for i = [2:11, 13:16]
    tdd = t_age_vec{i}(1);
    p = proj{i};
    [~, imx] = max(p);
    ii = imx:find(isfinite(p) & p>=0, 1, 'last');
    p(p < 0) = NaN;
    x = -log( p(ii) );
    y = lambda_W(i) * 60 * (tdd + (0:length(x)-1)');

    plot(x,y, '.-')
    text(x(end)+0.1,y(end), sprintf("%i", i))
end
axis tight
ylim([0, 0.4])
xlim([0, 8])
xlabel('nondimensional isotope age -log(g''/g''_{dd})')
ylabel('nondimensional water vapor age \lambda_W(t-t_{dd})')
axis square

saveas(gcf, 'W_iso_ages.eps', 'epsc')
saveas(gcf, 'W_iso_ages.svg')
saveas(gcf, 'W_iso_ages.png')