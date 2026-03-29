% calc_W_ent_res.m

lambda_mass_flux; % -> lambda_W, mf_flux, W_sfc
regress_transform_cp; % -> lambda_iso [1/min]

% Re-solve for mass flux W_ent as a residual for h = 350 m
% using eqn. 6.
ncp = size(mf_fss,1);
nt  = size(mf_fss,2);
W_ent_r = NaN(ncp,1);
W_sfc = NaN(ncp,1);

hconst = 350; % m

for i = [2:11, 13:ncp]
    [~, pk_ind] = min( mf_fen(i,:) );
    ind = pk_ind:find(isfinite(mf_fen(i,:)), 1, 'last'); % range, inc
    
    wq_sfc = -mf_flux(i,ind) ./ (rho * Lv); % units W/m^2 -> kg/kg m
    W_sfc(i) = mean(wq_sfc ./ (qsurf_cold(i,ind) - qbl_cold(i,ind)), 'omitnan');
    W_ent_r(i) = (lambda_iso(i)/60) * hconst - W_sfc(i);
end

rankdD_cronOrder = [7, 13, 10, 5, 3, 4, 1, 6, 8, 12, 14, 9, 11, 2];
W_ent_r( rankdD_cronOrder )
