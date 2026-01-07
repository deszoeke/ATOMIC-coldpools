% Copyright 2025 Simon P. de Szoeke and Estefanía Quiñones 
% Meléndez.
% 
% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation 
% files (the “Software”), to deal in the Software without 
% restriction, including without limitation the rights to use, 
% copy, modify, merge, publish, distribute, sublicense, and/or 
% sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following 
% conditions:
%
% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.
% plot_W_iso_ages.m

lambda_mass_flux; % -> lambda_W
regress_transform_cp;

rankdD_cronOrder = [7, 13, 10, 5, 3, 4, 1, 6, 8, 12, 14, 9, 11, 2];

clf()
subplot(1,1,1, 'fontsize',14)
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
        plot(x,y, '.-', 'LineWidth',1.4, 'MarkerSize',7)
    else
        plot(x,y,'.')
    end
    text(x(end)+0.1, y(end)+0.01, sprintf("%i", rankdD_cronOrder(j)), 'FontSize',12)
end
axis tight
ylim([-0.02, 0.4])
xlim([-0.4, 8])
plot([-0.4, 8], [0, 0], 'k-', 'LineWidth',0.2)
plot([0, 0], [-0.02, 0.4], 'k-', 'LineWidth',0.2)
plot([0, 0.4],[0, 0.4], 'k--', 'LineWidth',0.2)
xlabel({'nondimensional','isotope age -log(g''/g''_{dd})'})
ylabel({'nondimensional','water vapor age \lambda_{W,SBL}(t-t_{dd})'})
axis square

% saveas(gcf, 'W_iso_ages2.eps', 'epsc')
% saveas(gcf, 'W_iso_ages2.svg')
% saveas(gcf, 'W_iso_ages2.png')

% CP 9 and 10 have no data