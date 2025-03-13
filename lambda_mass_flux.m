% Estimate Lambda using surface flux from PSD dataset
% 
% Adjust this to average qsurf_cold and qent_cold over the whole wake,
% not just for the start of the cold pool.
%
% Assume that the cold pool wake recovery starts at the peak, 
% and ends when Stef NaNs it out.
%
% Might want to go 1-2 e-foldings, i.e. longer
% for 40 min cold pools, shorter for 1-min cold pools
% balance intensity dynamic range with representativeness
% --if you sample too long, other processes not representative
% of the wake start to dominate.
%
% mass flux velocity, m/s
%  W_src = E_src / (rho * Lv * (q-q_src))
% see eqns 9-13 in the paper

Lv    = 2.4e6; % units [J/kg]=[m^2/s^2]; Latent Heat of Vaporization of Water
h     = 700; % units in meters; SBL height
rho   = 1.1; % units in kg/m^3; air density

% load mixing fractions and surface fluxes
load mixing_fractions_vars.mat slope mf_t_cp mf_fen mf_fss tcoldw tcold mf_flux
load conserved_variables_10minLIMITED.mat q_en q_surf time_q

% need to reconstruct
ncp = size(mf_fss,1);
nt  = size(mf_fss,2);

% allocate variables for each cold pool
% q in mf coordinates
qsurf_cold = NaN(size(mf_fss));
qent_cold  = NaN(size(mf_fss));
for i = 1:ncp
    for k = 1:nt % k indexes time in one cold pool
        % interpolates 10 min nearest neighbor to 1 min
        % for reference: 
        % wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
        start_t_ids = find(time_q <= mf_t_cp(i,k), 'last');
        qsurf_cold(i,k)     = 1e-3 * q_surf(start_t_ids);
        qent_cold(i,k)      = 1e-3 * q_en(start_t_ids); % interpolated from soundings
    end
end
% qsurf_cold, qent_cold in kg/kg

% compute flux decay constants for each cold pool
lambda_ent = NaN(ncp,1);
lambda_sfc = NaN(ncp,1);
lambda_W   = NaN(ncp,1);
for i = 1:size(mf_fss,1)
    % wake index (wake_ind) corresponds to  T_min  time
    % peak index  (pk_ind)  corresponds to peak dD time
    [~, pk_ind] = min( mf_fen(i,:) );
    ind = pk_ind:find(isfinite(mf_fen(i,:), 'last')); % range, includes NaNs
    if i == 16 % 16 has 2 cold pool peaks and recoveries
        % wake_ind = 21;
        % wake_ind is unused. should this hack be for pk_ind? 
        % it doesn't look right. 
        % the first CP16 is at ind=14:23, second at 28:42
        % the indices will change if NaNs are spliced out so we don't do that
        ind = 14:23;
    end

    % kinematic q flux
    wq_sfc = -mean(mf_flux(i,ind),'omitnan') ./ (rho * Lv); % units W/m^2 -> kg/kg m s^-1 flux by surface
    % wq_ent should be multiplied by wq_sfc, E_sfc right???
    wq_ent = slope(i)*10^-4; % units [kg/kg m s^-1]; flux by entrainment
    % Is slope computed right? Dimensions of W were not right.
    % Check dimensions. rho and h were not used???

    % q-q_src in kg/kg units now
    % mass flux, m/s
    W_ent = wq_ent / mean(qent_cold(i,ind) - qent_cold(i,ind), 'omitnan');
    W_sfc = wq_sfc / mean(qsurf_cold(i,ind) - qsurf_cold(i,ind), 'omitnan');
    lambda_ent(i) = W_end / h;
    lambda_sfc(i) = W_sfc / h;
    lambda_W(i) = lambda_ent(i) + lambda_sfc(i); % decay coefficient
    % decay coefficient = lambda_W from fluxes (SPdeS)
end
