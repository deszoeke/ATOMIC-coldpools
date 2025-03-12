% eigcov.m
%
% Given vectors x and y, find the principal direction of variability,
% assuming the error is the same in x and y. Find the intersection xd,yd
% of the line in the principal direction through the centroid (mx, my),
% with the diagonal x+y=1, and xd at its y=0 intercept.
% SPdeS

[mx,my, v1,v2, d1,d2, sl1] = eigcov(x,y);

xd = xdiag(mx,my,sl1);
yd = ydiag(mx,my,sl1);
xb = xbase(mx,my,sl1);

anom = @(x) x-mean(x, "omitmissing");

% y-intercept, ybase=0
xbase = @(xm,ym,sl) xm-ym/sl;

% diagonal intercept
xdiag = @(xm,ym,sl) (sl*xm+1-ym)/(1+sl);
ydiag = @(xm,ym,sl) ym+sl*(1-ym-xm)/(1+sl);

function [mx,my, v1,v2, d1,d2, slope1] = eigcov(x,y)
    mx = mean(x, "omitmissing");
    my = mean(y, "omitmissing");
    [V,D] = eig(cov(x-xm, y-ym));
    % 2 eigenvectors
    v1 = V(:,1);
    v2 = V(:,2);
    % 2 eigenvalues
    d1 = D(1,1);
    d2 = D(2,2);
    slope1 = v1(2)/v1(1);
end