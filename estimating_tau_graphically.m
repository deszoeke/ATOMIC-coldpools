function [tau,tau_ss,tau_ee] = estimating_tau_graphically(i)
% using mixing fractions
load 'mixing_fractions_vars.mat' mf_fss mf_fee mf_fen mf_t_cp slope
pk_ind = find(mf_fen(i,:)==min(mf_fen(i,:)));
in_y = (mf_t_cp(i,pk_ind:end)-mf_t_cp(i,pk_ind))*24*60;
% y = log(in_y)./std(in_y,'omitnan'); % normalized by std. dev. before taking the log
y = log(in_y); % log of time in minutes

% With projection
fss_prime = mf_fss(i,:) - mf_fss(i,1);
% % with extrapolation
% extra_fss = max(fss_prime):0.01:fss_peak;
% fss_prime = [fss_prime extra_fss];
fee_prime = slope(i).*fss_prime;
% Redifining negative values in projection to avoid complex numbers when 
% evaluating the natural log
fee_prime(fee_prime<0) = NaN;
fss_prime(fss_prime<0) = NaN;
% Dividing by the mixing fraction at peak dD normalizes the projection
in_f = fee_prime(pk_ind:end)./(fee_prime(pk_ind)); % gives same values as fss_prime after normalizing
% in_f = fss_prime(pk_ind:end)./(fss_prime(pk_ind));
% flipping the axes!!!
% NOW THIS SECTION GIVES LAMBDA INSTEAD OF TAU!!!
x = log(in_f);

% Eliminating NaNs
y(isnan(x)) = [];
x(isnan(x)) = [];
x = x(2:end);
y = y(2:end);

% Using polyfit to calculate slope
p = polyfit(y(2:end),x(2:end),1);
f = polyval(p,y(2:end));
tau = -p(1);
disp(tau)
% lambda = -p(1);
% tau = 1/lambda;

% Using covariance matrix to calculate slope
% tau = sum(y-mean(y,'omitnan').*x-mean(x,'omitnan'))/sum(y-mean(y,'omitnan').^2);
% disp(tau)

figure; hold on;
plot(y,x,'o',y(2:end),f,'--r')
title(['tau = ',num2str(tau)])

%% Without projection
in_x = (mf_t_cp(i,pk_ind:end)-mf_t_cp(i,pk_ind))*24*60;
% x = log(in_x)./std(in_x,'omitnan'); % normalized by std. dev. before taking the log
x = log(in_x);

in_fee = mf_fee(i,pk_ind:end)./mf_fee(i,pk_ind);
in_fss = mf_fss(i,pk_ind:end)./mf_fss(i,pk_ind);

% Redifining negative values to avoid complex numbers when evaluating the natural log
in_fee(in_fee<0) = NaN;
in_fss(in_fss<0) = NaN;

yee = log(in_fee);
yss = log(in_fss);

% Eliminating NaNs
x(isnan(yee)) = [];
yss(isnan(yee)) = [];
yee(isnan(yee)) = [];

% Linear fit in log space
win = 15; % averaging window start point
pee = polyfit(x(win:end),yee(win:end),1);
fee = polyval(pee,x(win:end));
pss = polyfit(x(win:end),yss(win:end),1);
fss = polyval(pss,x(win:end));

lambda_ee = -pee(1);
tau_ee = 1/lambda_ee;
lambda_ss = -pss(1);
tau_ss = 1/lambda_ss;

% figure; hold on;
% plot(x,yee,'o',x(3:end),fee,'--r')
% legend('data','linear fit')
% title(['tau_e_e = ',num2str(tau_ee)])

% figure; hold on;
% plot(x,yss,'o',x(3:end),fss,'--r')
% legend('data','linear fit')
% title(['tau_s_s = ',num2str(tau_ss)])
