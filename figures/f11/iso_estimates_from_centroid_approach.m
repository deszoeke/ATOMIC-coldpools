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