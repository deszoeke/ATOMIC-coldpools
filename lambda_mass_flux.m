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
load mixing_fractions_vars.mat slope mf_t_cp mf_fen mf_fss mf_fee tcoldw tcold mf_flux
load conserved_variables_10minLIMITED.mat q_en q_surf time_q
H = load('conserved_variables_10min&mixing_fractions.mat'); % q_ob, time=time_q

% need to reconstruct
ncp = size(mf_fss,1);
nt  = size(mf_fss,2);

% allocate variables for each cold pool
% q in mf coordinates
qsurf_cold = NaN(size(mf_fss));
qent_cold  = NaN(size(mf_fss));
qbl_cold   = NaN(size(mf_fss));
for i = [2:11, 13:ncp]
    for k = find(isfinite(mf_t_cp(i,:))) % k indexes time in one cold pool
        % interpolates 10 min nearest neighbor to 1 min
        % for reference: 
        % wake_ind = find(mf_t_cp(i,:)>=tcoldw(i,1));
        start_t_ids = find(time_q <= mf_t_cp(i,k), 1, 'last');
        qsurf_cold(i,k)     = 1e-3 * q_surf(start_t_ids);
        qent_cold(i,k)      = 1e-3 * q_en(start_t_ids); % interpolated from soundings
        qbl_cold(i,k)       = 1e-3 * H.q_ob(start_t_ids);
    end
end
% qsurf_cold, qent_cold in kg/kg

% compute flux decay constants for each cold pool
lambda_ent = NaN(ncp,1);
lambda_sfc = NaN(ncp,1);
lambda_W   = NaN(ncp,1);
omg        = NaN(ncp,1);
for i = [2:11, 13:ncp]
    % wake index (wake_ind) corresponds to  T_min  time
    % peak index  (pk_ind)  corresponds to peak dD time
    [~, pk_ind] = min( mf_fen(i,:) );
    ind = pk_ind:find(isfinite(mf_fen(i,:)), 1, 'last'); % range, includes NaNs
    if i == 16 % 16 has 2 cold pool peaks and recoveries
        % wake_ind = 21;
        % wake_ind is unused. should this hack be for pk_ind? 
        % it doesn't look right. 
        % the first CP16 is at ind=14:23, second at 28:42
        % the indices will change if NaNs are spliced out so we don't do that
        ind = 14:23;
    end

    % omg = W_sfc/(W_sfc+W_ent) is y-intercept of the mixing fraction (no downdraft)
    [~,~,~,~,omg(i)] = cp_mf_regress_transform(mf_fss(i,ind), mf_fee(i,ind));

    % average delta qs over wake
    % omg(i) = mean(('omitnan')
    % kinematic fluxes wq
    wq_sfc = -mf_flux(i,ind) ./ (rho * Lv); % units W/m^2 -> kg/kg m s^-1 flux by surface
    % wq_ent should be multiplied by wq_sfc, E_sfc right???
    % wq_ent = slope(i)*10^-4; % units [kg/kg m s^-1]; flux by entrainment
    % Is slope computed right? Dimensions of W were not right.
    % Check dimensions. rho and h were not used???

    % q-q_src in kg/kg units now
    % mass flux, m/s
    W_sfc = mean(wq_sfc ./ (qsurf_cold(i,ind) - qbl_cold(i,ind)), 'omitnan');
    W_ent = ((1-omg(i))/omg(i)) * W_sfc;
    lambda_sfc(i) = W_sfc / h;
    lambda_ent(i) = W_ent / h;
    % decay coefficient from fluxes in wake, in absence of downdraft, 1/s
    lambda_W(i) = lambda_ent(i) + lambda_sfc(i);
end

%% function library

function [xm,ym, v1,v2, d1,d2, slope1] = eigcov(x,y)
% eigcov.m
%
% Given vectors x and y, find the principal direction of variability,
% assuming the error is the same in x and y. Find the intersection xd,yd
% of the line in the principal direction through the centroid (xm, ym),
% with the diagonal x+y=1, and xd at its y=0 intercept.
% SPdeS
    xm = mean(x, "omitmissing");
    ym = mean(y, "omitmissing");
    [V,D] = eig(cov(x-xm, y-ym));
    % 2 eigenvalues
    [S, I] = sort(diag(D), 'descend');
    d1 = S(1);
    d2 = S(2);
    % 2 eigenvectors
    v1 = V(:,I(1));
    v2 = V(:,I(2));

    slope1 = v1(2)/v1(1);
end

function [proj, sl1, xd,yd, xb] = cp_mf_regress_transform(x,y)
    % y-intercept, ybase=0
    xbase = @(xm,ym,sl) xm-ym/sl;
    
    % diagonal intercept
    xdiag = @(xm,ym,sl) (sl*xm+1-ym)/(1+sl);
    ydiag = @(xm,ym,sl) ym+sl*(1-ym-xm)/(1+sl);

    ii = isfinite(x) & isfinite(y);
    [xm,ym, v1,v2, d1,d2, sl1] = eigcov(x(ii), y(ii));

    xd = xdiag(xm,ym,sl1);
    yd = ydiag(xm,ym,sl1);
    xb = xbase(xm,ym,sl1);
    yb = 0; 

    [~, imx] = max(abs(v1));
    vn = v1 / sqrt(v1'*v1) * sign(v1(imx)); % unit normal

    % projection of x,y onto principal eigenvector
    proj = [x-xb; y-0]'*vn / sqrt((xd-xb)^2 + yd^2);
    % test: proj=0 at xb,0 and proj=1 at xd,yd
end
