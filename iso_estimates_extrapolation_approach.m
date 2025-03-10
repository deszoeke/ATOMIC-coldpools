function [dD_ent_mf,d18O_ent_mf] = iso_estimates_extrapolation_approach(q_ob,dD_ob,d18O_ob,q_surf,dD_surf,d18O_surf,q_ent)
    % EXTRAPOLATION APPROACH: Using single point in q*dD-q space
    % Allocating variables
    X1 = (q_ob.*d18O_ob)./10^3; % sp = single point
    X2 = (q_ob.*dD_ob)./10^3;

    Y1 = (q_surf.*d18O_surf)./10^3;
    Y2 = (q_surf.*dD_surf)./10^3;
    m1 = (X1 - Y1)./(q_ob - q_surf); % slope
    m2 = (X2 - Y2)./(q_ob - q_surf); % slope
    b  = X1 - (m1.*q_ob);
    b1 = Y1 - (m1.*q_surf); % same results for both points: good check
    bb = X2 - (m2.*q_ob);
    b2 = Y2 - (m2.*q_surf); % same results for both points: good check
    qXd18O_ent = (m1.*q_ent) + b1;
    qXdD_ent   = (m2.*q_ent) + b2;
    d18O_ent_mf = (qXd18O_ent.*10^3)./q_ent;
    dD_ent_mf = (qXdD_ent.*10^3)./q_ent;        
end
