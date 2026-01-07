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
function [m,b,bb] = mixing_line_slope_yint(q0,dD0,q1,dD1)
Y0 = q0.*dD0;
Y1 = q1.*dD1;
m = (Y1 - Y0)./(q1 - q0); % slope
bb = Y1 - (m.*q1); % y-intercept
b  = Y0 - (m.*q0); % same results for both points: good check

%% Test
% q0  =   1.1327e-05; %0.1; % g/kg
% dD0 = -44.6691; %-60; % permil
% q1  =  14.4619; %15;
% dD1 = -64.7553; %-66;
% figure; plot([q0 q1],[dD0 dD1],'ok','MarkerFaceColor','k')
% figure; plot([q0 q1],[q0 q1].*[dD0 dD1],'ok','MarkerFaceColor','k')
% Y0 = q0.*dD0;
% Y1 = q1.*dD1;
% m = (Y1 - Y0)./(q1 - q0); % slope
% bb = Y1 - (m.*q1); % y-intercept
% b  = Y0 - (m.*q0); % same results for both points: good check
% Xfit = linspace(q0,q1,10);
% Ylinearfit = (m.*Xfit) + b;
% hold on; plot(Xfit,Ylinearfit,'or')
% Ycurvefit = (Ylinearfit)./Xfit;
% hold on; plot(Xfit,Ycurvefit,'or')
