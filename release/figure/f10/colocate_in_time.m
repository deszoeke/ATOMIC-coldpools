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
function [new_var,new_time] = colocate_in_time(filename,old_var)
    % filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
    time_T = ncread(filename,'time');
    time_T = time_T/3600/24 + datenum('20200101','yyyymmdd');
    % T_L = ncread(filename,'tskin'); % liquid temperature at [???]m [in degrees C]
    % RH  = ncread(filename,'rhair')/100;
    
    % Co-locating RH & T_L variable in time %
    load 'data/2nd_leg_sounding_data_10min_linear_interp.mat' t
    pos_i = 999999*ones(size(t));
    for l = 1:length(t)
        pos2 = find(time_T>=t(l)); % rounding up to closest isotope surface data point!!!
                                        % try rounding to nearest data point
        pos_i(l) = pos2(1);
        clearvars pos2
    end
    time_RH = time_T(pos_i);
    % new_RH = RH(pos_i);
    new_var  = old_var(pos_i);
    new_time = time_RH;
end
