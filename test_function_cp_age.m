% test_function_cp_age.m
% Script to test & use cold_pool_age_estimates.m
% Created by Estefanía Quiñones Meléndez
% Last modified on Feb 28, 2024
i = 8;
% i = 16;
[tau,cp_age_f,~,pk_ind,f_peak] = cold_pool_age_estimates(i);
% INTRODUCE NEW ARGUMENT: f_var => should be f_peak or f_intersection
% for 8 and 16 is when f_peak will be used
%% For plotting purposes => 1:1 plots
mf_t = (mf_t_cp(i,:) - mf_t_cp(i,pk_ind))*24*60; % time in minutes
figure; hold on
plot(mf_t(pk_ind:end)./tau,cp_age_f(pk_ind:end)./tau,'-ok','LineWidth',1.5)
ylabel('vapor non-dimensional age')
xlabel('clock non-dimensional age')
title(['CP #',num2str(i),''])
grid on
%%
i = 7;
ideal = 8;
[tau,cp_age_f,~,pk_ind] = cold_pool_ages(i,ideal,f_peak);
%% For plotting purposes => 1:1 plots
% for negative tau
mf_t = (mf_t_cp(i,:) - mf_t_cp(i,pk_ind))*24*60; % time in minutes
figure; hold on
plot(mf_t(pk_ind:end)./-tau,cp_age_f(pk_ind:end)./-tau,'-ok','LineWidth',1.5)
ylabel('vapor non-dimensional age')
xlabel('clock non-dimensional age')
title(['CP #',num2str(i),''])
grid on

%% For plotting purposes => in mixing fraction space
% figure; hold on;
ideal = 16;
[~,~,~,~,f_peak] = cold_pool_age_estimates(ideal);
load 'mixing_fractions_vars.mat' mf_fss mf_fee Xeval Yeval slope
for i = [2:10,15:16]
    figure; hold on;
    [~,~,~,~,~,n,d,fss_prime,fee_prime] = cold_pool_ages(i,ideal,f_peak);
    scatter(mf_fss(i,:),mf_fee(i,:),5,n./d,'filled')
    plot([0 1],[1 0],'-k')
    plot([0 0],[1 0],'-k')
    plot([0 1],[0 0],'-k')
    plot([0 1],slope(i).*[0 1],'-r')
    plot(Xeval,Yeval(i,:),'--r')
    xlim([-.2 1])
    scatter(fss_prime,fee_prime,5,n./d,'filled')
    title(['CP #',num2str(i),''])
    colorbar
end

%% Testing tau estimates
mf_tau = zeros(1,16);
mf_tau(mf_tau==0) = NaN;
for i = 1:16 % [2:10,16] %
    [tau,~,~] = estimating_tau_graphically(i);
    % Saving tau for each CP
    mf_tau(i) = tau;
end
%% Without projection
mf_tau_ss = zeros(1,16);
mf_tau_ss(mf_tau_ss==0) = NaN;
for i = 1:16 % [2:10,16] %
    [~,tau_ss,~] = estimating_tau_graphically(i);
    % Saving tau for each CP
    mf_tau_ss(i) = tau_ss;
end
mf_tau_ee = zeros(1,16);
mf_tau_ee(mf_tau_ee==0) = NaN;
for i = 1:16 % [2:10,16] %
    [~,~,tau_ee] = estimating_tau_graphically(i);
    % Saving tau for each CP
    mf_tau_ee(i) = tau_ee;
end
%% Comparing different tau's
id = 1:16;
figure; hold on;
plot(id,mf_tau,'ob','MarkerFaceColor','b')
plot(id,mf_tau_proj,'*r')
plot(id,mf_tau_ss,'sk','MarkerFaceColor','k')
 xlabel('cold pool ID')
ylabel('time scale (\tau) [min]')
legend('based on fluxes','based on isotopes (projection)','based on isotopes')
