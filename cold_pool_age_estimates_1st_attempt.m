% cold_pool_age_estimates_1st_attempt.m
% Script to calculate cold pool age estimates
% Created by Estefanía Quiñones Meléndez
% Last modified on Nov 13, 2023
load 'mixing_fractions_vars.mat'
% Projection unto the slope line
i = 8;
O = [Xeval(1),Yeval(i,1)]; % origin
v = [Xeval(2),Yeval(i,2)]; % vector along the slope line
len_v = sqrt(v(1).^2 + v(2).^2);
u = v./len_v; % unit vector along the slope line
len_u = sqrt(u(1).^2 + u(2).^2);

x = [mf_fss(i,40),mf_fee(i,40)]; % dot to be projected unto line
proj = (sum(x.*v)/(norm(v)^2))*v;
% time = 
% equation = 
A = v;
B = x;
% C = proj;
% z represents time progression:
for i = 16
    xprime = mf_fss(i,:);
    xdd_ind = find(mf_dD_cp==max(mf_dD_cp(i,:))); % mf_fss @ fen = 0 (actually, @ max(dD))
    xdd = mf_fss(xdd_ind);
    x = xprime - xdd;
    %     y = yprime - ydd; % goes unused
    z = sqrt(1+(slope(i)).^2).*x;
end

%% Skipping the projection
% Normalizing cold pool times
% fron_ind = find(mf_t_cp(i,:)<=tcoldw(i,1))
% Plotting all decay behaviour in single plot
figure; hold on;
t_n = 1:61; % normalized time
plot(t_n,exp(-t_n),'--b','LineWidth',1.5)
plot(t_n,exp(-0.2.*t_n),'--g','LineWidth',1.5)
plot(t_n,exp(-0.1.*t_n),'--k','LineWidth',1.5)
plot(t_n,exp(-0.05.*t_n),'--r','LineWidth',1.5)
for i = [2:11,13:16]
    wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
    l = length(wake_ind);
    % plot(1:l,mf_fee(i,wake_ind)) %./max(mf_fee(i,wake_ind))
    plot(1:l,mf_fss(i,wake_ind)) %./max(mf_fee(i,wake_ind))
    % plot(1:l,mf_fen(i,wake_ind)) %./max(mf_fee(i,wake_ind))
end
% legend('\tau = 1 min','\tau = 5 min','\tau = 10 min','\tau = 20 min')
legend('\lambda = 1','\lambda = 0.2','\lambda = 0.1','\lambda = 0.05')

%% Model
U = 5; % m/s
C = 1870; % [J K-1 kg-1]; turbulent exchange coefficient/eddy difusivity
E_ent = slope(i); % flux by entrainment
E_sfc = -C.*U.*(mf_fen(i,:)-mf_fss(i,:)); % flux by surface
lambda_ent = E_ent./(mf_fen(i,wake_ind(1)) - mf_fen(i,1));
lambda_sfc = E_sfc./(mf_fss(i,wake_ind(1)) - mf_fss(i,1));
lambda = lambda_ent + lambda_sfc; % decay coefficient
tau = 1./lambda;
tau = mean(tau);
n = mf_fen(i,wake_ind(1)) - mf_fen(i,1); % numerator;
d = mf_fen(i,:) - mf_fen(i,1); % denominator
l = log(n./d); % logarithm
cp_age_fen = tau.*l; % time/age from the model
cp_age_fee = tau.*log((mf_fee(i,wake_ind(1)) - mf_fee(i,1))./(mf_fee(i,:) - mf_fee(i,1))); % time/age from the model
cp_age_fss = tau.*log((mf_fss(i,wake_ind(1)) - mf_fss(i,1))./(mf_fss(i,:) - mf_fss(i,1))); % time/age from the model
test = log((mf_fen(i,wake_ind(1)) - mf_fen(i,1))./(mf_fen(i,:) - mf_fen(i,1))); % time/age from the model
% Use the projection instead of the quotient to do the normalization along
% the diagonal line
% For plotting purposes
figure; hold on
plot(mf_t_cp(i,:),cp_age_fen,'-ob','LineWidth',1.5)
plot(mf_t_cp(i,:),cp_age_fee,'-ok','LineWidth',1.5)
plot(mf_t_cp(i,:),cp_age_fss,'-or','LineWidth',1.5)
legend('from f_e_n','from f_e_e','from f_s_s')
ylabel('cold pool age')
% xlabel('time in minutes')
plot(t_n,exp(-0.15.*t_n),'--g','LineWidth',1.5)
legend('/lambda = 1','/lambda = 0.05','/lambda = 0.1')
t_ob % time of observation

%% CODE BELOW NOT WORKING AS EXPECTED => FIX!!!
% figure; hold on;
% for i = [2:11,13:16]
%     fron_ind = find(mf_t_cp(i,:)<=tcoldw(i,1));
%     l = length(fron_ind);
%     plot(1:l,mf_fee(i,fron_ind))
% end

