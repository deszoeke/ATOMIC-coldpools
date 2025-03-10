function [dD_mf1,d18O_mf1] = iso_estimates_from_centroid_approach(qsurf,dDsurf,d18Osurf,q_mf1)
    % CENTROID APPROACH: Using centroid in the q*dD-q space
    for k = 1:length(qsurf)
        if k>=90 && k<=length(qsurf)-90
            bg = (k-89:k);
            q_cp_mean(k) = nanmean(q(bg)); % bg = background
            dD_cp_mean(k) = nanmean((q(bg).*dD(bg)))./10^3;
            d18O_cp_mean(k) = nanmean((q(bg).*d18O(bg)))./10^3;
        end
        clearvars bg
    end
    Y1 = (qsurf.*d18Osurf)./10^3;
    Y2 = (qsurf.*dDsurf)./10^3;
    m1 = (d18O_cp_mean - Y1)./(q_cp_mean - qsurf); % slope
    m2 = (dD_cp_mean - Y2)./(q_cp_mean - qsurf); % slope
    b = d18O_cp_mean - (m1.*q_cp_mean);
    b1 = Y1 - (m1.*qsurf); % same results for both points: good check
    bb = dD_cp_mean - (m2.*q_cp_mean);
    b2 = Y2 - (m2.*qsurf); % same results for both points: good check
    qXd18O_mf1 = (m1.*q_mf1) + b1;
    qXdD_mf1   = (m2.*q_mf1) + b2;
    d18O_mf1(k) = (qXd18O_mf1.*10^3)./q_mf1;
    dD_mf1(k) = (qXdD_mf1.*10^3)./q_mf1;        
end