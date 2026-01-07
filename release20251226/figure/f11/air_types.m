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