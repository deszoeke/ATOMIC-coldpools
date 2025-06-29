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
