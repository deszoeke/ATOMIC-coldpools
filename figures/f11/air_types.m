function [q_type,th_type] = air_types(q,theta,h,type,var1,var2,var3)
Rd = 287.04;
Cp = 1005.7;

switch type 
    case 'surface'
        % var2 must be in degrees C
%         disp('surface')
        SLP = var1; % [in hPa]
        T = var2; % [in degrees C]
        % addpath('C:\Users\quinones\Documents\Data\thermo')
        q_type = qs(SLP*100,T)*1e3; % in g/kg; slp must be in Pa and T in degrees C
        th_type = (T + 273.15).*(1e5./(SLP*1e2)).^(Rd/Cp); % (Potential Temp in degrees K)
    case 'downdraft'
        disp('dd')
    case 'entrained'
        disp('entrained')
    otherwise 
        disp('ERROR: invalid type entry')
end