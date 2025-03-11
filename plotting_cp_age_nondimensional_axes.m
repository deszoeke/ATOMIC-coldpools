% plotting_cp_age_nondimensional_axes.m
% Script to plot cold pool ages in non-dimensional space
% Created by Estefanía Quiñones Meléndez
% Last modified on Mar 11, 2024
%% Selecting intercept point in mixing fraction space
ideal = 8;
% ideal = 16;
[~,~,~,~,f_peak] = cold_pool_age_estimates(ideal);
disp(f_peak)
load 'mixing_fractions_vars.mat' mf_t_cp slope
% Add tau as a variable that's saved in the loop
%% Plot
% Pre-allocating tau variable
mf_tau = zeros(1,16);
mf_tau(mf_tau==0) = NaN;
figure; hold on;
for i = [2:11,13:16] % [8,16] % 
    [tau,cp_age_f,~,pk_ind,wake_ind] = cold_pool_ages(i,ideal,f_peak);
    disp(['CP #',num2str(i),'tau =',num2str(tau),''])
% Saving tau for each CP
    mf_tau(i) = tau; % mixing fraction tau = lambda_i*(t-t0)
    % SPdeS 2025-2-25
    % Dt = t-t0 is estimated by dimensionalizing the nondimensional age from
    % the (isotope) mixing fraction. Then the nondimensional age from the 
    % mass fluxes is computed as lambda_W*Dt.
    % For each cold pool we need the lambda_W, and either Dt or lambda_iso.
    % lambda_W from the fluxes and lambda_iso from the isotope decay are
    % computed in cold_pool_ages(). Simon's not sure what tau is. It could
    % be 1/lambda, or it could be nondimensional age Dt/lambda. Rework
    % inside of cold_pool_ages assuming Dt = the same const for flux (W) and isotope
    % (iso) methods. The time of the peak is irrelevant to the age Dt
    % except for the ideal cold pools where the time of the peak
    % is assumed to correspond to the time of the evaporative downdraft.
% Defining wall-clock non-dimesional age reference point
    % Referenced to peak dD
    mf_t = (mf_t_cp(i,:) - mf_t_cp(i,pk_ind))*24*60; % time in minutes
    % Adding offset so that it follows 1:1 line
    % [offset hardcoded]
    mf_t = mf_t + cp_age_f(pk_ind);
    
    % Referenced to baseline (t_max: before cold pool starts)
    % NOT an option, ALMOST follows 1:1 line, CP#8 & #16 DON'T START AT (0,0)
    % mf_t = (mf_t_cp(i,:) - mf_t_cp(i,1))*24*60; % time in minutes
    
    % Referenced to T_min (t_min)
    % NOT an option, ROUGHLY follows 1:1 line, CP#8 & #16 START close to (0,0)
    % mf_t = (mf_t_cp(i,:) - mf_t_cp(i,wake_ind(1)))*24*60; % time in minutes
    
    % Referenced to full recovery (t_end)
    % NOT an option, doesn't follow 1:1 line
    %     mf_t = mf_t_cp(i,:);
    %     mf_t(isnan(mf_t_cp(i,:))) = [];
    %     cp_age_f(isnan(mf_t_cp(i,:))) = [];
    %     mf_t = (mf_t - mf_t(end))*24*60; % time in minutes

% Plotting dimesional axes
    % plot(mf_t(pk_ind:end),cp_age_f(pk_ind:end),'-o','LineWidth',1.5)
%     plot(mf_t(pk_ind:end),cp_age_f(pk_ind:end),'-o','LineWidth',1.5)
% Plotting non-dimesional axes
    plot(mf_t(pk_ind:end)./-tau,cp_age_f(pk_ind:end)./-tau,'-o','LineWidth',1.5)
%     text(mf_t(pk_ind),cp_age_f(pk_ind),['CP#',num2str(i),'']);
end
    title(['All cold pools (using CP#',num2str(ideal),' as ideal)'])
    ylabel('vapor cold pool age')
    xlabel('observation time (normalized)')
    grid on
    axis equal
%     plot(0:80,0:80,'-k') % adding 1:1 line for reference
    plot(0:10,0:10,'-k') % adding 1:1 line for reference
    %     legend('CP#2','CP#3','CP#4','CP#5','CP#6','CP#7','CP#8','CP#9','CP#10',...
    %         'CP#11','CP#13','CP#14','CP#15','CP#16')
%% Plot
% Replacing zeros with NaNs
mf_tau_16(mf_tau_16==0) = NaN;
mf_tau_8(mf_tau_8==0) = NaN;
% Comparing tau for all CPs
figure; hold on;
plot(slope,-mf_tau,'.','MarkerSize',10)
xlabel('slope in mixing fraction space')
ylabel('time scale (\tau) [min]')
% Comparing tau for all CPs
load 'dD_&_Tmin_anomaly_vars.mat' T_min
figure; hold on;
plot(T_min,-mf_tau,'.','MarkerSize',10)
ylabel('time scale (\tau) [min]')
xlabel('T_m_i_n [^oC]')
