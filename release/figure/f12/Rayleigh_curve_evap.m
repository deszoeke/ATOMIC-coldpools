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
function [delt_rayleigh,frac,rayleigh_step]= Rayleigh_curve_evap(T0,R0)
% Reservoirs: liquid water drops
    Rvsmow = 155.76e-6; % unitless, deuterium
    frac = 0:0.05:1; %  20 values in linear dist of fraction
%   frac = 0:0.01:1; % 100 values in linear dist of fraction

% The following distribution of "f" values did not work as well as expected    
%     frac1 = 0.5 - 0.5.^(1:10); % increasing between 0 and 0.499
%     frac2 = 0.5 + 0.5.^(1:10); % decreasing between 1 and 0.501
%     frac = [frac1 0.5 flip(frac2)]; % increasing between 0 and 1

% The following distribution of "f" values work well    
%     frac1 = 1   - 0.5.^(2:10); % increasing between 0.75 and 0.999
%     frac2 =       0.5.^(1:10); % decreasing between 0.5 and 0.001
%     frac = [flip(frac2) frac1];% increasing between 0 and 0.999
% liquid water temp doesn't vary the same way as air temp when vapor is the reservoir, 
% so we use constant liquid water temp in this function
    T = T0; % in Kelvin
% alpha must be less than one for this reservoir (liquid water);
% so we calculate rayleigh_step using 1/alphae_d!!!
    alphae_d = 1./exp( 1158.8e-12 .*T.^3 - 1620.1e-9 .*T.^2 + 794.84e-6 .*T - 161.04e-3 + 2.9992e6./T.^3 );
    rayleigh_step = R0.*frac.^(alphae_d-1);
    delt_rayleigh = 1e3*(rayleigh_step/Rvsmow - 1); % permil