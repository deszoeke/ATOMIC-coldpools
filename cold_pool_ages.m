function [tau,cp_age_fee,cp_age_fss,pk_ind,wake_ind,n,d,fss_prime,fee_prime]=cold_pool_ages(i,ideal,fss_peak)
% Script to calculate cold pool age estimates for non-ideal cases
% Based on estimates from CP#8 & CP#16 analysis
% Created by Estefanía Quiñones Meléndez
% Last modified on Mar 06, 2024
load 'mixing_fractions_vars.mat' slope mf_t_cp mf_fen mf_fss tcoldw tcold
%% Input surface flux by hand to test code
% for reference:
wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
% E_sfc = +180.5783; % units W/m^2
% Lv    = 2.5e6; % units [J/kg]=[m^2/s^2]; Latent Heat of Vaporization of Water
% E_sfc = E_sfc./Lv; % units [kg m^-2 s^-1]
% % Changing units for the flux
% % h = 700; % units in meters; SBL height
% % rho = 1; % units in kg/m^3; air density
% E_ent = slope(i)*10^-4; % units [kg m^-2 s^-1]; flux by entrainment

%% Using surface flux from PSD dataset
% loading surface fluxes
load 'mixing_fractions_vars' mf_flux
% creating surface flux matrix
%     filename = ('EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc');
%     hl_bulk = ncread(filename,'hl_bulk'); % units W/m^2, surface_downward_latent_heat_flux
%     % hs_bulk = ncread(filename,'hs_bulk'); % units W/m^2, surface_downward_sensible_heat_flux
%     time = datenum('01012020','mmddyyyy') + ncread(filename,'time')/3600/24;
%     mf_flux = zeros(16,61);
%     mf_flux(mf_flux==0) = NaN;
%     for k = [2:11,13:16]
%         for j = 1:61
%             ind = find(time>=mf_t_cp(k,j));
%                 if isempty(ind)
%                     continue
%                 end
%             mf_flux(k,j) = hl_bulk(ind(1));
%         end
%     end
E_sfc = -mean(mf_flux(i,:),'omitnan'); % units W/m^2
Lv    = 2.5e6; % units [J/kg]=[m^2/s^2]; Latent Heat of Vaporization of Water
% h     = 700; % units in meters; SBL height
% rho   = 1; % units in kg/m^3; air density

% E_sfc = E_sfc./Lv./h./rho;
E_sfc = E_sfc./Lv; % units [kg m^-2 s^-1]; flux by surface

E_ent = slope(i)*10^-4; % units [kg m^-2 s^-1]; flux by entrainment

%% Craig-Gordon surface estimates
% Adjust this loop to get qsurf_cold and qent_cold for each cold pool timestep,
% not just for the start of the cold pool
% load 'RHS_Eq9_MJ79_LIMITED.mat' alpha_e alpha_e_D del_oc del_oc_D
load 'conserved_variables_10minLIMITED.mat' q_en q_surf time_q 
% delta_e_D = ((1./alpha_e_D)*(1+del_oc_D) - 1) * 1000;
mf_t = mf_t_cp(i,:);
mf_t(isnan(mf_t)) = [];
l = length(mf_t);
if i == 16
    l = 21;
end
for k = 1:l
    % for reference: wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
    ind = find(time_q >= tcold(i,k));
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
lambda = lambda_ent + lambda_sfc; % decay coefficient = lambda_W from fluxes (SPdeS)
tau = 1./lambda;
tau = tau./60;
if tau > 0
    tau = -tau;
end
disp(tau)
%% With projection
fss_prime = mf_fss(i,:) - mf_fss(i,1);
% % With extrapolation
% extra_fss = max(fss_prime):0.01:fss_peak;
% fss_prime = [fss_prime extra_fss];
fee_prime = slope(ideal).*fss_prime;
% Redifining negative values in projection to avoid complex numbers when 
% evaluating the natural log
fee_prime(fee_prime<0) = NaN;
fss_prime(fss_prime<0) = NaN;
% Dividing by the mixing fraction at peak dD normalizes the projection
in_fee = fee_prime./(fee_prime(wake_ind(1))); % gives same values as fss_prime after normalizing
in_fss = fss_prime./(fss_prime(wake_ind(1)));
n = in_fee - in_fee(1); % numerator
d = fss_peak./(fss_prime(wake_ind(1))) - in_fee(1); % denominator
cp_age_fee = tau.*log(n./d); % time/age from the model
% cp_age_fss = tau.*log(n./d); % time/age from the model
cp_age_fss = tau.*log(in_fee);

% SPdeS: I'm not sure about using this tau=-1/lambda_W from the fluxes here. 
% I get confused by tau. And it's the wrong tau (tau_W). So let's avoid tau.
% There should be another lambda, lambda_iso from the slope of log(in_fee) vs time.
% I don't think it's calculated here, but you have computed it and we need it here.
% From that get Dt = cp_age_fss (corrected) = -log(in_fee) / lambda_iso
% From this, set Dt at the peak, and then increment the time from the peak
% with the regular clock time axis. Dt is a vector. It's the same
% time vector whether we're thinking of a flux adjustment or an isotope
% adjustment, because the cold pool downdraft happened at some definite
% time, and the clock shows some definite later time. It's just the
% difference of those times.
% Then the two different nondimensional times on the plot are 
% Dt/lambda_W vs. Dt/lambda_iso.
% Using Dt/lambda_iso will make a pretty boring plot, with
% straight lines with slopes lambda_iso/lambda_W. But we actually have the
% isotope nondimensional time = -log( in_fee ) independently!
% We used it to generate the slope that we solved for \lambda_iso. 
% Plot it here again and plot two estimates of nondimensional time
% -log( in_fee ) vs. Dt/lambda_W.

% Adjusted
% cp_age_fee = -10^(4).*tau.*log(in_fee); % time/age from the model
% cp_age_fss = -10^(4).*tau.*log(in_fss); % time/age from the model
% use the same f values from to define the new triangle

%% For plotting purposes
% Constructing a new time vector for the extrapolated values
mf_t = mf_t_cp(i,:);
% 1:length(in_fee);
% mf_t(1:length(mf_t_cp(i,:))) = mf_t_cp(i,:);

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
end