function plotlims(x,y, tx,ty, t0,t1, varargin)
% plot (x, y) with tx,ty in the interval [t0, t1]
    yy = ( (t0 <= ty) & (ty <= t1) ); % indexes y
    xx = ( (t0 <= tx) & (tx <= t1) ); % indexes x
    plot(x(xx), y(yy), varargin{:});
end