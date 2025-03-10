addpath('C:\Users\quinones\Documents\Data\thermo')
thermo_constants

Rvsmow = 1; % will work relatively for defining R-delta reln., but not absolute

delta = @(R) 1e3*(R/Rvsmow - 1); % permil
Rati = @(d) Rvsmow * (1 + d*1e-3); % d in permil

% 0) givens
% initial conditions from subcloud layer air
R0 = Rati(-70);
q0 = 15e-3; % kg/kg
Ta = 298; % K
pa = 1e5; % Pa
T0 = Tlcl(pa, Ta, q0);
% LCL is initial temperature to start condensing

% final target
qfin = 1e-5; % ~1% of initial or nearly zero

% 1) Discretize q into nstep log-spaced steps
nstep = 100;
qs = logspace(log10(q0), log10(qfin), nstep+1)'; % decreasing (cumulative)
% qs = linspace(q0, qfin, nstep+1)'; % decreasing (cumulative)
fh = qs(2:end)./qs(1:end-1); % fractional step (const. if log-spaced)

% initialize T arrays
Ts = T0 + zeros(nstep+1, 1);
Th = T0 + zeros(nstep,   1);
Rvrs = zeros(nstep, 1); % Rayleigh ratio=R(i+1)/R(i) of steps

% 2) Iterate through qs steps, compute T and half-step Th
% Ts(1) = T0 % already
for istep = 1:nstep
    [Th(istep), Ts(istep+1)] = condense_T_Dlog_step( Ts(istep), qs(istep), qs(istep+1) );
    Rvrs(istep) = rayleigh_step(fh(istep), Th(istep));
end

% 3) Vector cumulative product to get R throughout Rayleigh process
fractions = qs./q0;
Rvs = cumprod( [R0; Rvrs] ); % remaining vapor isotope ratio
% ^Only works for steps evenly weighted with constant in Dlog(q)?
% Rvs = cumprod( [R0; fh.*Rvrs] ) ./ fractions; % doesn't seem to matter.

% isotope ratio of resulting liquid condensed in each step
Rls = (1-fh.*Rvrs) ./ (1-fh);
% cumulative liquid
ql = q0 - qs;
% diff(ql) = -diff(qs)

% cumulative liquid isotope ratio
% Rl = cumsum(Rls.*diff(ql))./cumsum(diff(ql));  % wrong

Rl_remove = (R0 - fractions.*Rvs) ./ (1-fractions); % total R of liquid removed
Rl_remove_step = NaN(nstep, 1);
for i = 1:nstep
    Rl_remove_step(i) = (Rvs(i) - fh(i).*Rvs(i+1)) ./ (1-fh(i));
end

% diagnostic plots
halfstep = sqrt(qs(2)/qs(1));

%{
plot(1-fractions, delta(Rvs));
ylabel('vapor \deltaD (permil)')
xlabel('fraction condensed')
clf
plot((1-fractions(1:end-1))*halfstep, Th)
ylabel('temperature (K)')
xlabel('midpoint fraction condensed')
clf
plot((1-fractions(1:end-1))*halfstep, alpha_e('D', Th))
ylabel('alpha_{eD}')
xlabel('midpoint fraction condensed')
clf
plot(Th, alpha_e('D', Th))
%}

clf
plot(delta(Rvs), Ts, '.-');
hold on
axis ij
plot(delta(Rl_remove_step), Th, '.-')
coi = get(gca, 'ColorOrderIndex');
plot(delta(Rl_remove(2:end)), Ts(2:end), '.-')
set(gca, 'ColorOrderIndex', coi)
plot(delta(alpha_e('D', T0)*R0), T0, 'o')
plot(ql*1e3, Ts ,'-')
plot(delta([R0 R0]),[292 195], 'k--')
plot([0 0],[292 195], 'k-')
set(gca, 'fontsize',16)
xlabel('\deltaD (permil)')
ylabel('temperature (K), q (g/kg)')
l = flipud(get(gca, 'children'));
legend(l([1:3 5]), {'vapor remains', 'local liquid removed', 'total liquid removed', 'q_{liquid}'}, 'location','southwest')
legend boxoff
axis tight
xlim([-800, 10]); ylim([200 293])
xlim([-70, 15]); ylim([260 293])
xlim([-20, 12]); ylim([260 293])

% saveas(gcf, 'liquid_dD.png')
% saveas(gcf, 'liquid_dD.svg')

function [Th, T1] = condense_T_Dlog_step( T0, q0, q1 )
% Take a full moist adiabatic condensation step from q0,T0 -> q1,T1
% Compute Th of the half step [to get alpha_e(Th)] and T1 at the end of the
% step.

% constants not inherited from lexical or calling scope
Cp = 1.0057e3;
Rv = 461.5;
Rd = 287.04;
Rd_Rv = Rd/Rv;

cee = @(T, q) (Rd.*T + Lv(T).*q) ./ (Lv(T).*Rd_Rv - Cp.*T);

% before step
q1_q0 = q1/q0;
Dlogq = log(q1_q0);

% take half step
Dlogqh = Dlogq/2;
expDlogqh = exp(Dlogqh);
DlogTh = cee(T0,q0) * Dlogqh;

% evaluate Th, qh, alphah at half step
qh = q0 * expDlogqh;
Th = T0 * exp(DlogTh);
% alphah = alpha_e(isotopologue, Th);

% full step, integrated, needed for the next step
T1 = T0 * q1_q0^cee(Th, qh);
end

function fnate = rayleigh_step(fraction, Th)
% Rayleigh fractionate vapor isotope ratio by the returned factor, 
% evaluated if fraction evaporates at temperature Th. Equilibrium
% fractionation is evaluated at Th.

fnate = fraction.^(alpha_eD(Th)-1); % =R1/R0
end

function a = alpha_eD(T)
% Deuterium equilibrium fractionation for liquid/vapor
a = exp( 1158.8e-12 .*T.^3 - 1620.1e-9 .*T.^2 + ...
         794.84e-6 .*T - 161.04e-3 + 2.9992e6./T.^3 );
end

function a = alpha_e18O(T)
% oxygen-18 equilibrium fractionation for liquid/vapor
a = exp( -7.685e-3 + 6.7123 ./T - 1.6664e3 ./T.^2 ...
         + 0.35041e6 ./T.^3 );
end