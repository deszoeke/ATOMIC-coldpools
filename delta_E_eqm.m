function [d_E, d_E_approx] = delta_E_eqm(h_ob,isotopologue,d_a,T_L)
% d_a = -18.12; % input d_a in permil
% h_ob = 0.433; % input RH in decimal point; unitless
% isotopologue = '18O/16O'; % input isotopologue in string format
% T_L = 13; % input liquid temperature in degrees Celsius

% Diffusivity ratios from Hellmann and Harvey (2020)
T = (T_L+274.15)/100; % degrees K
DrHDO  = 0.98258 - (0.02546./T) + (0.02421./(T.^(5/2)));
Dr18O  = 0.96671 + (0.007406./(T.^(1/2))) - (0.004861./(T.^(3)));
T_oc = (T_L+274.15); % degrees K

% Diffusivity ratios from Merlivat (1978)
% DrHDO = 0.9757; % 𝐷HD^16O/D𝐻_2^16O
% Dr18O = 0.9727; % 𝐷H_2^18O/D𝐻_2^16O

n   = 0.5; 
% n = 1/2 for an open water body under natural conditions
theta = 0.88; % scaling factor/weighting term
% 0.50 for evaporation in the eastern Mediterranean Sea
% 0.88 for the North American Great Lakes
% alpha is the temperature-dependent fractionation factor of the vapour to liquid phase transition.

switch isotopologue 
    case '18O/16O'
    alpha_star = exp( -7.685/10^3 + 6.7123./T_oc - 1.6664*(10^3./T_oc.^2)...
                 + 0.35041*(10^6./T_oc.^3) ); 
                 % for oxygen: 18O/16O; for 20 C, it gives 1.0098
    % eps_k = 18.9; % in permil
    nCd = ((Dr18O).^(-n)) - 1; % C_d = molecular diffusivity ratio  
    d_L = 1; % in permil from Sebastian's Meteor Values
    % d18O_L = 1; % from Gat 1996

    case 'D/H'
    alpha_star = exp( 1158.8*(T_oc.^3./10^12) - 1620.1*(T_oc.^2./10^9)...
                 + 794.84*(T_oc./10^6) - 161.04/10^3 + 2.9992*(10^6./T_oc.^3));
                 % for deuterium: D/H ; for 20 C, it gives 1.0844
    % eps_k = [???];
    nCd = ((DrHDO).^(-n)) - 1; 
    d_L = 7; % in permil from Sebastian's Meteor Values
    % dD_L = 5; % from Gat 1996

    otherwise
    disp('ERROR - invalid entry for isotopologue variable - valid options are ''18O/16O'' or ''H/D')
end

% Calculating delta_epsilon using equation 6b from Gat 1996
d_eps = (1-h_ob).*theta.*nCd; % delta_epsilon; unitless
eps_star = 1 - 1./alpha_star; % unitless
% d_E_Xiao = ((d_L./alpha_star) - h_ob.*d_a - eps_star.*10^3 - eps_k.*(1-h_ob))./((1-h_ob) + (1-h_ob).*(eps_k/10^3));
d_E = (alpha_star.*d_L - h_ob.*d_a - eps_star.*10^3 - d_eps.*10^3)./((1-h_ob) + d_eps); % from Gat's 1996 Craig-Gordon model
d_E_approx = (d_L - h_ob.*d_a - eps_star.*10^3 - d_eps.*10^3)./(1-h_ob); % from Gat's 1996 Craig-Gordon model
end