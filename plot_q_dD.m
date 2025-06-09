% plot_q_dD.m
% 2025-06-05 SPdeS

% addpath('~/Documents/MATLAB/user_tools/thermo/')
% addpath('./thermo') % symbolic link <- ~/Documents/MATLAB/user_tools/thermo/

Rvsmow = 155.76e-6; % unitless, deuterium
zm = 17.;
zrf = 400.;

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
dq10 = adj_t(qstar10,L10, 17,400);qi_e
dt10mid = datetime10 + minutes(5);

% mixing fractions and mass flux data
lambda_mass_flux; % mf_t_cp, mf_fen
cp_dt = datetime(mf_t_cp, 'ConvertFrom','datenum'); % -> datetime % array has no background before cold pools
H.dt = datetime(H.time, 'ConvertFrom','datenum'); % 10-min

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
% qiPadj 

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
delta = @(qi,q) qi./(Rvsmow*q) - 1;
qi = @(q,di) Rvsmow*(1+di) .* q; % di unitless (not permil), qi has units of q, 
mean_delta = @(q,d) 1e3*mean( delta( qi(q,1e-3*d), q ), 'omitnan'); % d permil input and output
mixt = @(q1,q2,n) q1:(q2-q1)/n:q2; % represents mixture between q1 and q2
xweight = @(s,a,e) (a-s)/(e-s); % e fraction of a in e-a----s

% test
% 1e3 * delta( mixt(qien,qisf,20), mixt(qen,qsf,20) )

cmap = colororder();
cmjt = hsv(7);
cmjt = cmjt([1,2,3,5,6,7],:);
pic = false;
picadj = true;
psl = false;
psladj = false;
k_ent = 4; % height to get entrainment end member
ms = 8; % markersize

%{
% color wake obs by time
for hi = 1:6
    m0 = 10*(hi-1);
    % inc = @(t) tonset(j)+minutes(m0) <= t & t <= tonset(j)+minutes(m0+10);
    inc = @(t) cp_dt(i,ind(1))+minutes(m0) <= t & t <= cp_dt(i,ind(1))+minutes(m0+10);
    plot(qPadj(inc(isotime)), dD(inc(isotime)), 'o', ...
        'markersize',ms, 'markerfacecolor',cmjt(i,:), 'markeredgecolor','k')
end
%}

clf();
is = [8, 16]; % ranks 1 and 2
for j = 1:2  % chronological index
    i = is(j);
    subplot(2,2,j, 'fontsize',14); hold on
    % j=1; i=8;
    % j=2; i=16;
    [~, pk_ind] = min( mf_fen(i,:) );
    ind = pk_ind:find(isfinite(mf_fen(i,:)), 1, 'last'); % range, includes NaNs
    
    % background: median of 60 min prior
    inbg = @(dt) ( (tonset(j)-minutes(60))<=dt & (dt<tonset(j)) );
    bb = inbg(dt_10min)'; % indexes q_surf, q_ent, dD_surf, dD_ent
    bb = bb & isfinite(dD_surf);
    q_surf_b = mean( q_surf( bb) );
    qi_surf_b = mean( qi(q_surf(  bb        ),1e-3.*dD_surf(      bb)) , 'omitnan');
    qi_adj_b  = mean( qi(qPadj(inbg(isotime)),1e-3.*dD(inbg(isotime))) , 'omitnan');
    q_ent_b  = mean( q_ent(4,bb) );
    q_adj_b = mean(qPadj(inbg(isotime)),'omitnan');
    % reextrapolate dD_ent from background centroid
    x = xweight(q_surf_b, q_adj_b, q_ent_b);
    qi_ent_b = (qi_adj_b - (1-x)*qi_surf_b) / x; % extrapolate using moisture
    dD_ent = 1e3*delta(qi_ent_b, q_ent_b); % permil
    plotqd(q_surf_b, qi_surf_b, q_ent_b, qi_ent_b, 20)

    if pic % plot Picarro q
        plotlims(qP, dD, isotime,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k');
        plotlims(qP, dD, isotime,isotime, tonset(j),cp_dt(i,ind(1)), '.-', 'color',cmap(1+2*(j-1),:), 'markersize',ms);
        plotlims(qP, dD, isotime,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2+2*(j-1),:));
    end
    if picadj % plot adjusted Picarro q
        plotlims(qPadj, dD, isotime,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k');
        plotlims(qPadj, dD, isotime,isotime, tonset(j),cp_dt(i,ind(2)), '.-', 'color',cmap(1,:), 'markersize',ms);
        plotlims(qPadj, dD, isotime,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2,:));
    end
    if psl % plot unadjusted PSL q
        plotlims(qair, dD, mettime,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k');                          % background prior to cold pool
        plotlims(qair, dD, mettime,isotime, tonset(j),cp_dt(i,ind(1)), '.-', 'color',cmap(1+2*(j-1),:), 'markersize',ms);                % onset where T drops
        plotlims(qair, dD, mettime,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2+2*(j-1),:)); % wake where T recovers
    end
    if psladj % plot adjusted PSL q
        plotlims(q_adj, dD, dt_adj,isotime, tonset(j)-minutes(60),tonset(j), '.', 'color','k');
        plotlims(q_adj, dD, dt_adj,isotime, tonset(j),cp_dt(i,ind(1)), '.-', 'color',cmap(1+2*(j-1),:), 'markersize',ms);
        plotlims(q_adj, dD, dt_adj,isotime, cp_dt(i,ind(1)),cp_dt(i,ind(1))+hours(1), '.-', 'color',cmap(2+2*(j-1),:));
    end

    text(11, -65, char(uint8('a')-1+j), 'fontsize',14)
    title(string(tonset(j)), 'fontweight','normal')
    axis square
    ax = gca();
    ax.XAxis.MinorTick       = 'on';
    ax.XAxis.MinorTickValues = 10:1:22;
    ax.TickLength=[0.04, 0.04];
    xlabel('q (10^{-3} kg kg^{-1})')
    ylabel(['\deltaD (', char(8240), ')'])
    ylim([-76, -64]); xlim([10, 22])
end

saveas(gcf, 'q-dD_cps.eps', 'epsc')
saveas(gcf, 'q-dD_cps.svg')
saveas(gcf, 'q-dD_cps.png')
saveas(gcf, 'q-dD_cps.pdf')

% CP 9 and 10 have no data
