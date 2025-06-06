% plot_q_dD.m
% 2025-06-05 SPdeS

% addpath('~/Documents/MATLAB/user_tools/thermo/')
% addpath('./thermo') % symbolic link <- ~/Documents/MATLAB/user_tools/thermo/

zm = 17;
zrf = 400;

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

% 10-min flux variables for height adjustments
filename10 = 'data/EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc';
qstar10 = ncread(filename10,'qstar'); % specific humidity scaling parameter ['g/kg']
Tstar10 = ncread(filename10,'tstar'); % temperature scaling parameter [K]
slp10 = 100*ncread(filename10,'psealevel'); % atmospheric pressure at sea level [hPa --> Pa]
L10 = ncread(filename10,'MO_length'); % Monin Obukhov length scale [m] all negative values; valid range: [-500,0]
datenum10 = ncread(filename10,'time')/3600/24 + datenum('20200101','yyyymmdd'); % in 10-min resolution
datetime10 = datetime(datenum10, 'ConvertFrom','datenum');
dq10 = adj_t(qstar10,L10, 17,400);
dt10mid = datetime10 + minutes(5);

% mixing fractions and mass flux data
lambda_mass_flux; % mf_t_cp, mf_fen
cp_dt = datetime(mf_t_cp, 'ConvertFrom','datenum'); % -> datetime % array has no background before cold pools
H.dt = datetime(H.time, 'ConvertFrom','datenum'); % 10-min
% for extrapolating dD_en from background and surface end member
H.x = H.fen./(H.fen+H.fss); % for extrapolating entrainment source

mf_dt_cp = datetime(mf_t_cp, 'ConvertFrom', 'datenum');
% regress_transform_cp; % don't need

% isotope data
isofile = 'data/EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
isoinfo = ncinfo(isofile);
isotime = datetime(2020,1,1,0,0,0) + seconds(ncread(isofile, 'time'));
isodatenum = datenum(isotime);
dD = ncread(isofile, 'dD');
qP = ncread(isofile, 'q'); % Picarro
qPadj = qP + interp1(dt10mid,dq10, isotime);
qiPadj 

% met data
metfile = 'data/EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
metinfo = ncinfo(metfile);
mettime = datetime(2020,1,1,0,0,0) + seconds(ncread(metfile, 'time')); % datetime
qair = ncread(metfile, 'qair');
qadj = qair + interp1(dt10mid,dq10, mettime, 'linear');

% [th_adj,q_adj,time_adj] = height_adj_1min_spd(zrf);
% dt_adj = datetime(time_adj, 'ConvertFrom','datenum'); % -> datetime


rankdD_cronOrder = [7, 13, 10, 5, 3, 4, 1, 6, 8, 12, 14, 9, 11, 2];
tonset = [datetime(2020,2,7,12,29,0); datetime(2020,2,10,15,59,0)];

% isotope helper functions
delta = @(qi,q) qi./q - 1;
mixt = @(q1,q2,n) q1:(q2-q1)/n:q2;
% qi = @(q,di) (1+di) .* q;
% test
% 1e3 * delta( mixt(qien,qisf,20), mixt(qen,qsf,20) )
pic = false;

clf(); subplot(1,1,1, 'fontsize',14); hold on
cmap = colororder();
% j = 0; % chronological index
% for i = [8] %, 16] % ranks 1 and 2
%     j = j + 1;
j=1; i=8;
    [~, pk_ind] = min( mf_fen(i,:) );
    ind = pk_ind:find(isfinite(mf_fen(i,:)), 1, 'last'); % range, includes NaNs
    
    % background: median of 60 min prior
    bb = ( (tonset(j)-minutes(60) <= dt_10min) & (dt_10min < tonset(j)) ); % indexes q_surf, q_ent, dD_surf, dD_ent
    q_ent_b  = mean( q_ent(4,bb) );
    q_surf_b = mean( q_surf( bb) );
    qi_ent  = mean( q_ent( 4,bb).*dD_ent( 4,bb) , 'omitnan');
    qi_surf = mean( q_surf(  bb).*dD_surf(  bb) , 'omitnan');
    % re-extrapolate qi_ent
    inbg = @(dt) ( (tonset(j)-minutes(60))<=dt & dt<tonset(j) );
    qi_bg = mean( qPadj(inbg(isotime)) , 'omitnan');
    x_bg  = mean( H.x(inbg(H.dt)) , 'omitnan');
    qi_en = (qi_bg - (1-x_bg)*qi_surf) / x_bg; % extrapolate entrainment source from background and surface

    qx = mixt(q_ent_b,q_surf_b,20);
    % dy = delta( mixt(qi_ent,qi_surf,20), mixt(q_ent_b,q_surf_b,20) );
    dy = delta( mixt(qi_en,qi_surf,20), mixt(q_ent_b,q_surf_b,20) );
    plot( qx, dy, 'k-')
    plot( qx(1), dy(1), 'ko', 'MarkerEdgeColor','k', 'MarkerFaceColor','k') % entrainment
    plot( qx(end), dy(end), 'ko', 'MarkerEdgeColor','k', 'MarkerFaceColor',0.7+[0, 0, 0]) % surface

    % plot background prior to cold pool
    plotlims(qair,dD, mettime,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k')
    % onset where T drops
    plotlims(qair, dD, mettime,isotime, tonset(j),cp_dt(i,ind(1)), '.-', 'color',cmap(1+2*(j-1),:));
    % wake where T recovers
    plotlims(qair, dD, mettime,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2+2*(j-1),:));
    
    % % plot background prior to cold pool
    % plotlims(q_adj,dD, dt_adj,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k')
    % % onset where T drops
    % plotlims(q_adj, dD, dt_adj,isotime, tonset(j),cp_dt(i,ind(1)), '.-', 'color',cmap(1+2*(j-1),:));
    % % wake where T recovers
    % plotlims(q_adj, dD, dt_adj,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2+2*(j-1),:));

    if pic % plot Picarro q, with offset
        plotlims(qPadj, dD, isotime,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k');
        plotlims(qPadj, dD, isotime,isotime, tonset(j),cp_dt(i,ind(1)), '.-', 'color',cmap(1+2*(j-1),:));
        plotlims(qPadj, dD, isotime,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2+2*(j-1),:));
    end

    title(string(tonset(j)), 'fontweight','normal')
% end

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

function plotlims(x,y, tx,ty, t0,t1, varargin)
% plot (x, y) with tx,ty in the interval [t0, t1]
    yy = ( (t0 <= isotime) & (isotime <= t1) ); % indexes y
    xx = ( (t0 <= mettime) & (mettime <= t1) ); % indexes x
    plot(x(xx), y(yy), varargin{:});
end