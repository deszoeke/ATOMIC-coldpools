function [dD_type,d18O_type] = air_types_iso(dD,d18O,type,var1,var2,var3)

switch type 
    case 'surface'
        % disp('surface')
        T_L = var1;
        rh  = var2;
        % Calculating d_E's %
        [d18O_e, ~] = delta_E_eqm(rh,'18O/16O',d18O,T_L);
        [dHDO_e, ~] = delta_E_eqm(rh,  'D/H'  ,dD,T_L);
        % Solving for delta_a %
        h = var3; 
        n = 0.5; 
        % n = 1/2 for an open water body under natural conditions
        theta = 0.88; % scaling factor/weighting term
        % 0.50 for evaporation in the eastern Mediterranean Sea
        % 0.88 for the North American Great Lakes
        
        % Diffusivity ratios from Hellmann and Harvey (2020)
        T = (T_L+274.15)/100; % degrees K
        DrHDO  = 0.98258 - (0.02546./T) + (0.02421./(T.^(5/2)));
        Dr18O  = 0.96671 + (0.007406./(T.^(1/2))) - (0.004861./(T.^(3)));
        
        nCdHDO = ((DrHDO).^(-n)) - 1; % C_d = molecular diffusivity ratio
        nCd18O = ((Dr18O).^(-n)) - 1; % C_d = molecular diffusivity ratio  
        
        d_epsHDO = (1-h)*theta*nCdHDO; % delta_epsilon; unitless
        d_eps18O = (1-h)*theta*nCd18O; % delta_epsilon; unitless
        
        % Calculating alphas %
        % alpha is the temperature-dependent fractionation factor of the vapour to liquid phase transition
        T_oc = (T_L+274.15); % degrees K
        alpha_D = exp( 1158.8*(T_oc.^3./10^12) - 1620.1*(T_oc.^2./10^9)...
                     + 794.84*(T_oc./10^6) - 161.04/10^3 + 2.9992*(10^6./T_oc.^3));
                     % for deuterium: D/H ; for 20 C, it gives 1.0844
        alpha_O = exp( -7.685/10^3 + 6.7123./T_oc - 1.6664*(10^3./T_oc.^2)...
                     + 0.35041*(10^6./T_oc.^3) ); 
                     % for oxygen: 18O/16O; for 20 C, it gives 1.0098
        dHDO_L = 7; % in permil from Sebastian's Meteor Values
        d18O_L = 1; % in permil from Sebastian's Meteor Values
        % dHDO_L = 5; % from Gat 1996
        % d18O_L = 1; % from Gat 1996
        eps_star_O = 1 - 1./alpha_O; % unitless
        eps_star_D = 1 - 1./alpha_D; % unitless
        
        delta_aO = ((alpha_O) .* d18O_L - d18O_e.*((1-h)+d_eps18O) - eps_star_O.*10^3 - d_eps18O.*10^3)./(h);
        delta_aD = ((alpha_D) .* dHDO_L - dHDO_e.*((1-h)+d_epsHDO) - eps_star_D.*10^3 - d_epsHDO.*10^3)./(h);
        % delta_a0 = ((1./alpha_D) .* dHDO_L - dHDO_e.*((1-h)+0) - eps_star_D - 0)./(h);
        
        dD_type   = delta_aD;
        d18O_type = delta_aO;
        % DXS_a = delta_aD - 8*delta_aO;

%         load 'RHS_Eq9_MJ79_variablesUPDATED.mat' alpha_e_D alpha_e del_oc_D del_oc % loads 10-min data!!!
%         T = var1;
%         D2H  = 0.98258 - (0.02546./T) + (0.02421./(T.^(5/2)));
%         D18O = 0.96671 + (0.007406./(T.^.5)) - (0.004861./(T.^3));
%         h_iso = linspace(0.5,1,5); % range based on ??? data (from above molecular layer to 20m height)
%         % If we use h = 1, alpha_k_O and alpha_k_D are both 1 and the
%         % dependency on Temp is lost!!! That can't be right!!!
%         alpha_k_O = (1./D18O).*(1-h_iso(end))+ h_iso(end); % alpha_k must be greater than 1; check!!!
%         alpha_k_D = (1./D2H).*(1-h_iso(end)) + h_iso(end); % alpha_k must be greater than 1; check!!!
%         axa_D = alpha_e_D.*alpha_k_D'; % 10-min data!!!
%         axa_O = alpha_e  .*alpha_k_O'; % 10-min data!!!
%         % Solving for the isotopic concentration del_v0_D
%         delta_e_D = ((1./axa_D).*(1+del_oc_D) - 1) * 1000; % 10-min data!!!
%         delta_e_O = ((1./axa_O).*(1+del_oc) - 1) * 1000; % 10-min data!!!
%         dD_type = delta_e_D;   % using values from h=1;
%         d18O_type = delta_e_O; % using values from h=1;
    case 'downdraft'
        disp('dd')
    case 'entrained'
        disp('entrained')
    otherwise 
        disp('ERROR: invalid type entry')
end
