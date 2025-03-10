%% data availability plot: isotopes %%
% RHB data
filename = 'EUREC4A_ATOMIC_RonBrown_Isotope-Analyzer_1min_20200126-20200210_v1.0.nc';
t_vapor  = ncread(filename,'time')/3600/24 + datenum('20200101','yyyymmdd');
vapor    = 5*ones(size(t_vapor));
% P-3 data
date = ['0117';'0119';'0123';'0124';'0131';'0203';'0204';'0205';'0209';'0210';'0211']; % all flight dates
n = length(date);
t_p3 = 9999*ones(1,2*n);
for k = 1:n
    filename = ['EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-1Hz_2020',num2str(date(k,1:4)),'_v1.1.nc'];
    tt = ncread(filename,'time')/3600/24 + datenum('20200101','yyyymmdd');
    t_p3(2*k-1) = tt(1);
    t_p3(2*k)   = tt(end);
end
p3       = 4*ones(size(t_p3));
% Precip data
filename = 'EUREC4A_ATOMIC_RonBrown_Precipitation-Isotope-Ratios_20200105-20200212_v1.0.nc';
t_prec   = ncread(filename,'collection_time')/3600/24 + datenum('20200101','yyyymmdd');
prec     = 3*ones(size(t_prec));
% CTD casts data
data  = readtable('RonBrown_CTD-Isotope-Ratios.csv');
t_ctd = datenum(data.Var1);
ctd   = 2*ones(size(t_ctd));
% Sea surface data
data = readtable('RonBrown_Surface-Isotope-Ratios.csv');
t_surf = datenum(data.Var1);
surf   = 1*ones(size(t_surf));
%% Plotting
figure; hold on;
t = t_vapor;
y = vapor;
yvalues = [y(1)-0.5 y(1)-0.5 y(1)+0.5 y(1)+0.5];
xvalues = [t(1) t(end) t(end) t(1)];
j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
t = t_p3;
y = p3;
for k = 1:n
    yvalues = [y(1)-0.5 y(1)-0.5 y(1)+0.5 y(1)+0.5];
    xvalues = [t(2*k-1) t(2*k) t(2*k) t(2*k-1)];
    j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
    j.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
plot(t_prec,prec,'.k',t_ctd,ctd,'.k',t_surf,surf,'.k','MarkerSize',16)
datetick('x','mmm/dd','keeplimits','keepticks')
yticks([1 2 3 4 5])
yticklabels({'Sea surface','CTD casts','Rainwater','P-3 vapor','RHB vapor'})
box on
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
plot([t_prec(1)-1 t_surf(end)+2],[1.5 1.5],'-k')
plot([t_prec(1)-1 t_surf(end)+2],[2.5 2.5],'-k')
plot([t_prec(1)-1 t_surf(end)+2],[3.5 3.5],'-k')
plot([t_prec(1)-1 t_surf(end)+2],[4.5 4.5],'-k')
legend('continous sampling','discrete sampling','Location','northoutside')
ylim([0.5 5.5])
xlim([t_prec(1)-.9 t_surf(end)+2])
ax = gca;
ax.XGrid = 'on';
saveas(gcf,'data_availability_isotopes.fig')
saveas(gcf,'data_availability_isotopes.png')

%% data availability plot: traditional meteorology %%
% Radiosonde data
filename = 'EUREC4A_RonBrown_Vaisala-RS_L2_v3.0.0.nc';
t_rad    = ncread(filename,'launch_time')/3600/24 + datenum('20200101','yyyymmdd');
rad      = 5*ones(size(t_rad));
    % % Ceilometer data
    % % filename = 'EUREC4A_ATOMIC_RonBrown_Ceilometer_1hour_20200109-20200212_v1.1.nc';
    % % t_ceil   = ncread(filename,'time')/3600/24 + datenum('20200101','yyyymmdd');
    % % ceil     = 3*ones(size(t_ceil));
% Surface data [PSD data & Disdrometer rain data % Ceilometer data]
filename = 'EUREC4A_ATOMIC_RonBrown_1min_nav_met_sea_20200109-20200212_v1.3.nc';
t_psd    = double(ncread(filename,'time'))/3600/24 + datenum('20200101','yyyymmdd');
psd      = 4*ones(size(t_psd));

%% Plotting
figure; hold on;
t = t_psd;
y = psd;
yvalues = [y(1)-0.5 y(1)-0.5 y(1)+0.5 y(1)+0.5];
xvalues = [t(1) t(end) t(end) t(1)];
patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
% t = t_;
% y = ;
% yvalues = [y(1)-0.5 y(1)-0.5 y(1)+0.5 y(1)+0.5];
% xvalues = [t(1) t(end) t(end) t(1)];
% j = patch(xvalues, yvalues, [0.9 0.9 0.9],'LineStyle','none');
% j.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_rad,rad,'.k','MarkerSize',16)
plot([t_prec(1)-1 t_surf(end)+2],[3.5 3.5],'-k')
plot([t_prec(1)-1 t_surf(end)+2],[4.5 4.5],'-k')
plot([t_prec(1)-1 t_surf(end)+2],[5.5 5.5],'-k')
legend('continous sampling','discrete sampling','Location','northoutside')
datetick('x','mmm/dd','keeplimits','keepticks')
yticks([4 5])
yticklabels({'Surface','Radiosondes'})
box on
set(findall(gcf,'-property','Fontsize'),'FontSize',16)
ylim([0.5 5.5])
xlim([t_prec(1)-.9 t_surf(end)+2])
ax = gca;
ax.XGrid = 'on';
saveas(gcf,'data_availability_trad-met.fig')
saveas(gcf,'data_availability_trad-met.png')

%% NEXT STEP: good quality data availability plot %%
% In house ship contamination flag
% color bad data in red