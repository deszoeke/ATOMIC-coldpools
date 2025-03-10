function alpha = alpha_e(isotopologue, T)
% Returns equilibrium fractionation coefficient for either D/H or 18O/16O. 
% The coefficient is computed from the passed temperature(s)
% 
% Inputs:
% isotopologue: str
%     Either 'D' or '18O'.
% T: scalar or array-like
%     Temperature(s) in kelvin.


switch upper(isotopologue)
  case {'2H', 'D'}
    alpha = exp( 1158.8e-12 .*T.^3 - 1620.1e-9 .*T.^2 ...
                  + 794.84e-6 .*T - 161.04e-3 + 2.9992e6./T.^3 );
  case '18O'
    alpha = exp( -7.685e-3 + 6.7123 ./T - 1.6664e3 ./T.^2 ...
                  + 0.35041e6 ./T.^3 );
end
end