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
% eigcov.m
%
% Given vectors x and y, find the principal direction of variability,
% assuming the error is the same in x and y. Find the intersection xd,yd
% of the line in the principal direction through the centroid (xm, ym),
% with the diagonal x+y=1, and xd at its y=0 intercept.
% SPdeS

load data/mixing_fractions_vars
% plot(mf_fss', mf_fee', '.', linestyle="none")
% hold on; plot([0, 1, 0, 0], [0, 0, 1, 0], "k")
% axis equal
ncp  = size(mf_fee, 1);
proj = cell(ncp,1);
t_age_vec = cell(ncp,1);
slope = NaN(ncp,1);
xd    = NaN(ncp,1);
yd    = NaN(ncp,1);
xb    = NaN(ncp,1);
lambda_iso = NaN(ncp,1);
dt_dd_src  = NaN(ncp,1);

% start test plot
plot([0 1 0 0], [0 0 1 0], "k")
hold on

for cpi = 2:size(mf_fee,1) % each cold pool
    % start at the isotope peak and regress just the recovery
    [~, imx] = max(mf_fss(cpi,:) + mf_fee(cpi,:)); % typical slope is ~1
    [p,b,c,d,e] = cp_mf_regress_transform(mf_fss(cpi,imx:end), mf_fee(cpi,imx:end));
    proj{cpi} = p;
    slope(cpi) = b; 
    xd(cpi) = c;
    yd(cpi) = d;
    xb(cpi) = e; % y-intercept (no downdraft)

    p(p<0) = NaN;
    % start at the isotope peak
    [~, imx] = max(p);
    pp = p(imx:end);
    fii = find(isfinite(pp));
    Pf = polyfit(fii-1, log(pp(fii)), 1);
    % isotope decay rate, 1/minute
    lambda_iso(cpi) = -Pf(1); 
    % extrapolate to the evaporative downdraft origin time
    % time since the downdraft source when the peak is observed at the ship
    % dt_dd_src(cpi) = -log(pp(1)) / lambda_iso(cpi); % minutes
    dt_dd_src(cpi) = -Pf(2) / lambda_iso(cpi); % minutes
    
    t_age_vec{cpi} = dt_dd_src(cpi) + (0:length(pp));

    % add test plot
    plot(mf_fss(cpi,imx:end),mf_fee(cpi,imx:end),'.-')
    plot(mf_fss(cpi,imx),mf_fee(cpi,imx),'o')
    plot([xb(cpi),xd(cpi)],[0,yd(cpi)])
    
    % subplot(2,1,2)
    % plot(t_age_vec{cpi}(fii), log(pp(fii)))

end


anom = @(x) x-mean(x, "omitmissing");

% functions (must come at end of script)

function [xm,ym, v1,v2, d1,d2, slope1] = eigcov(x,y)
% eigcov.m
%
% Given vectors x and y, find the principal direction of variability,
% assuming the error is the same in x and y. Find the intersection xd,yd
% of the line in the principal direction through the centroid (xm, ym),
% with the diagonal x+y=1, and xd at its y=0 intercept.
% SPdeS
    xm = mean(x, "omitmissing");
    ym = mean(y, "omitmissing");
    [V,D] = eig(cov(x-xm, y-ym));
    % 2 eigenvalues
    [S, I] = sort(diag(D), 'descend');
    d1 = S(1);
    d2 = S(2);
    % 2 eigenvectors
    v1 = V(:,I(1));
    v2 = V(:,I(2));

    slope1 = v1(2)/v1(1);
end

function [proj, sl1, xd,yd, xb] = cp_mf_regress_transform(x,y)
    % y-intercept, ybase=0
    xbase = @(xm,ym,sl) xm-ym/sl;
    
    % diagonal intercept
    xdiag = @(xm,ym,sl) (sl*xm+1-ym)/(1+sl);
    ydiag = @(xm,ym,sl) ym+sl*(1-ym-xm)/(1+sl);

    ii = isfinite(x) & isfinite(y);
    [xm,ym, v1,v2, d1,d2, sl1] = eigcov(x(ii), y(ii));

    xd = xdiag(xm,ym,sl1);
    yd = ydiag(xm,ym,sl1);
    xb = xbase(xm,ym,sl1);
    yb = 0; 

    [~, imx] = max(abs(v1));
    vn = v1 / sqrt(v1'*v1) * sign(v1(imx)); % unit normal

    % projection of x,y onto principal eigenvector
    proj = [x-xb; y-0]'*vn / sqrt((xd-xb)^2 + yd^2);
    % test: proj=0 at xb,0 and proj=1 at xd,yd
end
