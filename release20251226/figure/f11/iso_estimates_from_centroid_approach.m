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
function [dD_mf1,d18O_mf1] = iso_estimates_from_centroid_approach(q,dD,d18O,qsurf,dDsurf,d18Osurf,q_mf1)
    % CENTROID APPROACH: Using centroid in the q*dD-q space
    % Allocating variables
    aver_win = 10; % 50; % 90; % averaging window for sensitivity purposes
       q_cp_mean = 999999*ones(size(qsurf));
    d18O_cp_mean = 999999*ones(size(qsurf));
      dD_cp_mean = 999999*ones(size(qsurf));
    for k = 1:length(qsurf)
        if k>=aver_win && k<=length(qsurf)-aver_win
            bg = (k-(aver_win-1):k);
            q_cp_mean(k) = mean(q(bg),'omitnan'); % bg = background
            dD_cp_mean(k) = mean((q(bg)'.*dD(bg)),'omitnan')./10^3;
            d18O_cp_mean(k) = mean((q(bg)'.*d18O(bg)),'omitnan')./10^3;
        end
        clearvars bg
    end
       q_cp_mean(q_cp_mean   ==999999) = NaN;
    d18O_cp_mean(d18O_cp_mean==999999) = NaN;
      dD_cp_mean(dD_cp_mean  ==999999) = NaN;

    Y1 = (qsurf.*d18Osurf')./10^3;
    Y2 = (qsurf.*dDsurf')./10^3;
    m1 = (d18O_cp_mean - Y1)./(q_cp_mean - qsurf); % slope
    m2 = (dD_cp_mean - Y2)./(q_cp_mean - qsurf); % slope
    b = d18O_cp_mean - (m1.*q_cp_mean);
    b1 = Y1 - (m1.*qsurf); % same results for both points: good check
    bb = dD_cp_mean - (m2.*q_cp_mean);
    b2 = Y2 - (m2.*qsurf); % same results for both points: good check
    qXd18O_mf1 = (m1.*q_mf1) + b1;
    qXdD_mf1   = (m2.*q_mf1) + b2;
    d18O_mf1 = (qXd18O_mf1.*10^3)./q_mf1;
    dD_mf1 = (qXdD_mf1.*10^3)./q_mf1;        
end