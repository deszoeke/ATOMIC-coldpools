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
function plot_time_comp(ax, t_comp2, y, color)

    % background mean for each event
    bkm = @(T) mean(T(:,t_comp2<-19.7),2);
    % standard deviation of background mean among events
    sbk = @(T) std(bkm(T));

    me = mean(y, 'omitnan');
    % se = sqrt( std(y, 'omitnan').^2 - sbk(y)^2 ); % can be negative
    se = std(y-bkm(y), 'omitnan');

    plot(ax, t_comp2, me, '-', 'linewidth',1.4, 'color',color);
    hold(ax, 'on');
    plot(ax, t_comp2, me + se, '-', 'linewidth',0.5, 'color',color);
    plot(ax, t_comp2, me - se, '-', 'linewidth',0.5, 'color',color);
    plot(ax, t_comp2, median(y, 'omitnan'), '.', 'color',color);
end