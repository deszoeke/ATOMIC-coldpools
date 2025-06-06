% plot_q_dD.m
% 2025-06-05 SPdeS

% addpath('~/Documents/MATLAB/user_tools/thermo/')
% addpath('./thermo') % symbolic link <- ~/Documents/MATLAB/user_tools/thermo/

% get CG surface, entrained, and downdraft end members
kp = allchild(groot); % keep existing figures
conserved_properties_for_isotopes;
for x = allchild(groot)
    if ~ismember(x,kp)
        close(x) % close other opened figures
    end
end
dt_10min = datetime(t_10min, 'ConvertFrom','datenum');
% dD_ent(4,:)
% q_ent(4,:)
% dD_surf
% q_surf


% mixing fractions and mass flux data
lambda_mass_flux; % mf_t_cp, mf_fen
cp_dt = datetime(mf_t_cp, 'ConvertFrom','datenum'); % -> datetime
H.dt = datetime(H.time, 'ConvertFrom','datenum');

% regress_transform_cp; % don't need

% isotope data
isofile = 'data/EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
isoinfo = ncinfo(isofile);
isotime = datetime(2020,1,1,0,0,0) + seconds(ncread(isofile, 'time'));
dD = ncread(isofile, 'dD');
qP = ncread(isofile, 'q'); % Picarro

% met data
metfile = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
metinfo = ncinfo(metfile);
mettime = datetime(2020,1,1,0,0,0) + seconds(ncread(metfile, 'time'));
qair = ncread(metfile, 'qair');
[th_adj,q_adj,time_adj] = height_adj_1min_spd(zrf);


rankdD_cronOrder = [7, 13, 10, 5, 3, 4, 1, 6, 8, 12, 14, 9, 11, 2];
tonset = [datetime(2020,2,7,12,29,0); datetime(2020,2,10,15,59,0)];

% isotope helper functions
delta = @(qi,q) qi./q - 1;
mixt = @(q1,q2,n) q1:(q2-q1)/n:q2;
% qi = @(q,di) (1+di) .* q;
% test
% 1e3 * delta( mixt(qien,qisf,20), mixt(qen,qsf,20) )

clf(); subplot(1,1,1, 'fontsize',14); hold on
cmap = colororder();
j = 0; % chronological index
for i = [8] %, 16] % ranks 1 and 2
    j = j + 1;
    [~, pk_ind] = min( mf_fen(i,:) );
    ind = pk_ind:find(isfinite(mf_fen(i,:)), 1, 'last'); % range, includes NaNs
    
    % background: median of 60 min prior
    bb = ( (tonset(j)-minutes(60) <= dt_10min) & (dt_10min < tonset(j)) ); % indexes q_surf, q_ent, dD_surf, dD_ent
    q_ent_b  = mean( q_ent(4,bb) );
    q_surf_b = mean( q_surf( bb) );
    qi_ent  = mean( q_ent( 4,bb).*dD_ent( 4,bb) , 'omitnan');
    qi_surf = mean( q_surf(  bb).*dD_surf(  bb) , 'omitnan');
    qx = mixt(q_ent_b,q_surf_b,20);
    dy = delta( mixt(qi_ent,qi_surf,20), mixt(q_ent_b,q_surf_b,20) );
    plot( qx, dy, 'k-')
    plot( qx(1), dy(1), 'ko', 'MarkerEdgeColor','k', 'MarkerFaceColor','k') % entrainment
    plot( qx(end), dy(end), 'ko', 'MarkerEdgeColor','k', 'MarkerFaceColor',0.7+[0, 0, 0]) % surface

    bb1 = ( (tonset(j)-minutes(60) <= isotime) & (isotime < tonset(j)) ); % indexes isotime, q, dD
    % plot(qP(bb1)-0.4, dD(bb1), '.', 'color','k');
    mm1 = ( (tonset(j)-minutes(60) <= mettime) & (mettime < tonset(j)) );
    plot(qair(mm1), dD(bb1), '.', 'color','k');
    % onset where T drops
    ii0 = ( (tonset(j) <= isotime) & (isotime <= cp_dt(i,ind(1))) ); % indexes isotime, q, dD
    mm0 = ( (tonset(j) <= mettime) & (mettime <= cp_dt(i,ind(1))) ); % indexes qair from met tower
    % plot(qP(ii0)-0.4, dD(ii0), '.-', 'color',cmap(1+2*(j-1),:));
    plot(qair(mm0), dD(ii0), '.-', 'color',cmap(1+2*(j-1),:));
    % wake where T recovers
    iir = ( (cp_dt(i,ind(1)) <= isotime) & (isotime <= cp_dt(i,ind(1))+hours(1)) ); % indexes isotime, q, dD
    mmr = ( (cp_dt(i,ind(1)) <= mettime) & (mettime <= cp_dt(i,ind(1))+hours(1)) ); % indexes qair from met tower
    % plot(qP(iir)-0.4, dD(iir), '.-', 'color',cmap(2+2*(j-1),:));
    plot(qair(mmr), dD(iir), '.', 'color',cmap(2+2*(j-1),:));
    title(string(tonset(j)), 'fontweight','normal')
end
axis tight
xlabel('q (10^{-3} kg kg^{-1})')
ylabel(['\deltaD (', char(8240), ')'])
ylim([-76, -64]); xlim([11, 22])

% saveas(gcf, 'W_iso_ages.eps', 'epsc')
% saveas(gcf, 'W_iso_ages.svg')
% saveas(gcf, 'W_iso_ages.png')

saveas(gcf, 'W_iso_ages2.eps', 'epsc')
saveas(gcf, 'W_iso_ages2.svg')
saveas(gcf, 'W_iso_ages2.png')

% CP 9 and 10 have no data