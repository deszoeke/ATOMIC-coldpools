function [tau,cp_age_fee,cp_age_fss,pk_ind,fss_peak,fee_peak]=cold_pool_age_estimates(i)
% Script to calculate cold pool age estimates
% Created by Estefanía Quiñones Meléndez
% Last modified on Dec 11, 2023
    %% Plotting all decay behaviour in single plot
    load 'mixing_fractions_vars.mat' slope mf_t_cp mf_fen mf_fss mf_fee tcoldw
%     figure; hold on;
%     t_n = 1:61; % normalized time
%     plot(t_n,exp(-t_n),'--b','LineWidth',1.5)
%     plot(t_n,exp(-0.2.*t_n),'--g','LineWidth',1.5)
%     plot(t_n,exp(-0.1.*t_n),'--k','LineWidth',1.5)
%     plot(t_n,exp(-0.05.*t_n),'--r','LineWidth',1.5)
    for k = [2:11,13:16]
        wake_ind = find(mf_t_cp(k,:)>=tcoldw(k,1));
        l = length(wake_ind);
    %     plot(1:l,mf_fee(i,wake_ind)) %./max(mf_fee(i,wake_ind))
    %     plot(1:l,mf_fss(i,wake_ind)) %./max(mf_fee(i,wake_ind))
    %     plot(1:l,mf_fen(k,wake_ind)) %./max(mf_fee(i,wake_ind))
    end
    % legend('\tau = 1 min','\tau = 5 min','\tau = 10 min','\tau = 20 min')
%     legend('\lambda = 1','\lambda = 0.2','\lambda = 0.1','\lambda = 0.05')
    % title('f_e_e')
    % title('f_s_s')
%     title('f_e_n')
%% Input surface flux by hand to test code
% for reference:
wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
E_sfc = +180.5783; % units W/m^2
Lv    = 2.5e6; % units [J/kg]=[m^2/s^2]; Latent Heat of Vaporization of Water
E_sfc = E_sfc./Lv; % units [kg m^-2 s^-1]
% Changing units for the flux
% h = 700; % units in meters; SBL height
% rho = 1; % units in kg/m^3; air density
%% Using surface flux from PSD dataset

% E_sfc = ; % density units [kg m^-2 s^-1]

E_ent = slope(i)*10^-4; % flux by entrainment

%% Craig-Gordon surface estimates
% Adjust this loop to get qsurf_cold and qent_cold for each cold pool timestep,
% not just for the start of the cold pool
% load 'RHS_Eq9_MJ79_LIMITED.mat' alpha_e alpha_e_D del_oc del_oc_D
load 'conserved_variables_10minLIMITED.mat' q_en q_surf time_q 
% delta_e_D = ((1./alpha_e_D)*(1+del_oc_D) - 1) * 1000;
l = length(mf_t_cp);
if i == 16
    l = 21;
end
for k = 1:l
    % for reference: wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
    ind = find(time_q >= tcoldw(i,k));
    start_t_ids    = ind(1);
    qsurf_cold(k)     = q_surf(start_t_ids);
    % dDsurf_cold(i,k)= delta_e_D(start_t_ids);
    qent_cold(k)      = q_en(start_t_ids);
end
%% Model
% Using the same tau for every cp %
% q's are in g/kg units!
% wake index (wake_ind) corresponds to  T_min  time
% peak index  (pk_ind)  corresponds to peak dD time
pk_ind = find(mf_fen(i,:)==min(mf_fen(i,:)));
% 1 index corresponds to the background/enviroment time
% lambda_ent = E_ent./(qent_cold(i) - qent_cold(i-1));
% lambda_sfc = E_sfc./(qsurf_cold(i) - qsurf_cold(i-1));
if i == 16
    wake_ind(1) = 21;
end
lambda_ent = E_ent./(qent_cold(wake_ind(1)) - qent_cold(1));
lambda_sfc = E_sfc./(qsurf_cold(wake_ind(1)) - qsurf_cold(1));
lambda = lambda_ent + lambda_sfc; % decay coefficient
tau = 1./lambda;
tau = tau./60;
if tau < 0
    tau = -tau;
end
%% With projection
fss_prime = mf_fss(i,:) - mf_fss(i,1);
fee_prime = slope(i).*fss_prime;
fee_peak = max(fee_prime);
fss_peak = max(fss_prime);
% Redifining negative values in projection to avoid complex numbers when 
% evaluating the natural log
fee_prime(fee_prime<0) = NaN;
fss_prime(fss_prime<0) = NaN;
% Dividing by the mixing fraction at peak dD normalizes the projection
in_fee = fee_prime./(fee_prime(wake_ind(1))); % gives same values as fss_prime after normalizing
% in_fss = fss_prime./(fss_prime(wake_ind(1))); %
n = in_fee - in_fee(1); % numerator
d = in_fee(pk_ind) - in_fee(1); % denominator
cp_age_fee = tau.*log(n./d); % time/age from the model
% cp_age_fss = tau.*log(n./d); % time/age from the model
cp_age_fss = tau.*log(in_fee);
% Adjusted
if i == 16
cp_age_fee = -cp_age_fee; % time/age from the model
cp_age_fss = -cp_age_fss; % time/age from the model
end
% cp_age_fee = -10^(4).*tau.*log(in_fee); % time/age from the model
% cp_age_fss = -10^(4).*tau.*log(in_fss); % time/age from the model
% use the same f values from to define the new triangle

%% For plotting purposes
% figure; hold on;
% plot(mf_t_cp(i,:),mf_fss(i,:),'-')
% % plot(mf_t_cp(i,:),in_fss,'-.')
% plot(mf_t_cp(i,:),in_fee,'-.')
% plot(mf_t_cp(i,:),mf_fee(i,:),'-')
% legend('f_s_s','f''','f_e_e')
% % legend('f_s_s','f''_s_s','f''_e_e','f_e_e')
% title(['CP #',num2str(i),''])
% ylabel('f')
% datetick('x','HH:MM PM','keeplimits','keepticks')

%% For plotting purposes
% figure; hold on
% plot(mf_t_cp(i,:),cp_age_fee,'-ok','LineWidth',1.5)
% % plot(mf_t_cp(i,wake_ind(1):end),cp_age_fee(wake_ind(1):end),'-ok','LineWidth',1.5)
% % plot(mf_t_cp(i,:),cp_age_fss,'-or','LineWidth',1.5)
% % legend('from f_e_e','from f_s_s')
% ylabel('cold pool age (in minutes)')
% datetick('x','HH:MM PM','keeplimits','keepticks')
% title(['CP #',num2str(i),''])
% legend('from f')
% grid on

%% Trim cold pool to start at the peak dD time
% Use NaN's instead of skipping the values as input to the log function!!!
% To avoid issues when plotting
% Calculate the error (absolute error): is it less than 10% 

%% Calculating the surface flux
% % CODE BELOW NOT WORKING AS EXPECTED => FIX!!!
% U = 5; % m/s
% C = 1870; % [J K-1 kg-1]; turbulent exchange coefficient/eddy difusivity
% E_sfc = -C.*U.*(q_en(i,:)-q_ss(i,:)); % flux by surface
%% Model without projection
% n = mf_fen(i,wake_ind(1)) - mf_fen(i,1); % numerator;
% d = mf_fen(i,:) - mf_fen(i,1); % denominator
% l = log(n./d); % logarithm
% cp_age_fen = tau.*l; % time/age from the model
% cp_age_fee = tau.*log((mf_fee(i,wake_ind(1)) - mf_fee(i,1))./(mf_fee(i,:) - mf_fee(i,1))); % time/age from the model
% cp_age_fss = tau.*log((mf_fss(i,wake_ind(1)) - mf_fss(i,1))./(mf_fss(i,:) - mf_fss(i,1))); % time/age from the model
% test = log((mf_fen(i,wake_ind(1)) - mf_fen(i,1))./(mf_fen(i,:) - mf_fen(i,1))); % time/age from the model
% % Use the projection instead of the quotient to do the normalization along the diagonal line
%% For plotting purposes: ZOOM-IN
% figure; hold on
% plot(mf_t_cp(i,pk_ind:end),cp_age_fee(pk_ind:end),'-ok','LineWidth',1.5)
% % plot(mf_t_cp(i,wake_ind:end),cp_age_fee(wake_ind(1):end),'-ok','LineWidth',1.5)
% % plot(mf_t_cp(i,:),cp_age_fen,'-ob','LineWidth',1.5)
% % plot(mf_t_cp(i,:),cp_age_fee,'-ok','LineWidth',1.5)
% % plot(mf_t_cp(i,:),cp_age_fss,'-or','LineWidth',1.5)
% % legend('from f_e_n','from f_e_e','from f_s_s')
% ylabel('cold pool age (in minutes)')
% datetick('x','HH:MM PM','keeplimits','keepticks')
% title(['CP #',num2str(i),''])
% legend('from f','Location','southeast')
% grid on

% Adding decay model line for comparison
% plot(t_n,exp(-0.15.*t_n),'--g','LineWidth',1.5)
% legend('/lambda = 1','/lambda = 0.05','/lambda = 0.1')
% t_ob % time of observation
end