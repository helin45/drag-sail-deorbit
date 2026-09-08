P_ROOT = fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(P_ROOT,'lib'))); cd(P_ROOT);  % repo root + lib on path
%% =========================================================================
%  MONTE CARLO SENSITIVITY ANALYSIS
%  Parameters: theta0, dphi0, rho_scale
%  Fixed: phi=85, alpha=0.75, h=450km
%  Duration: 1 orbital period per sample
% =========================================================================
clear; clc; close all; clear functions

load('sam_results.mat')
load('inertia_tensors.mat','inertiaDB');

mu = 3.986004418e14; Re = 6.3781e6;

%% Fixed simulation parameters (matching ROADM/SAM sweep)
year_ep  = 2025; month_ep = 1; day_ep = 15; hour_utc = 12;
f107a = 150; f107 = 150; ap = 15;
ode_reltol = 1e-6; ode_abstol = 1e-8;
d_theta_deg = 2.0;

adbsat_base = fullfile(project_root(), 'adbsat_processed');

%% Reference case: phi=85, alpha=0.75, h=450km
phi_sel   = 85;
accom_str = '0p75';
alt       = 450;

a_orb  = Re + alt*1e3;
T_orb  = 2*pi*sqrt(a_orb^3/mu);
v_circ = sqrt(mu/a_orb);

I_idx   = find(inertiaDB.angles_deg == phi_sel);
Iyy_val = inertiaDB.I(2,2,I_idx);
g_idx   = find(ref_geom.phi_deg == phi_sel);

fname   = fullfile(adbsat_base, sprintf('%dkm/%ddeg_CLL_accom_%s.mat', alt, phi_sel, accom_str));
raw     = load(fname);
Cm_sort = raw.aedb.aero.Cm_BY(:); Cm_sort(isnan(Cm_sort)) = 0;
Cd_sort = -raw.aedb.aero.Cf_wX(:); Cd_sort(isnan(Cd_sort)) = 0;

clear sails_c
sails_c(1).angle    = phi_sel;
sails_c(1).alpha    = (-179.5:1:179.5)';
sails_c(1).Cm_pitch = Cm_sort;
sails_c(1).Cd       = Cd_sort;
sails_c(1).Iyy      = Iyy_val;
sails_c(1).A_ref    = ref_geom.A_ref(g_idx);
sails_c(1).L_ref    = ref_geom.L_ref(g_idx);

p.a_orb  = a_orb;   p.inc    = deg2rad(0);
p.RAAN   = deg2rad(90);
p.T_orb  = T_orb;   p.v_circ = v_circ;
p.year   = year_ep; p.month  = month_ep;
p.day    = day_ep;  p.hour   = hour_utc;
p.f107a  = f107a;   p.f107   = f107; p.ap = ap;
p.Iyy    = Iyy_val; p.sails  = sails_c; p.n_sails = 1;

%% =========================================================================
%  MONTE CARLO SAMPLE DISTRIBUTIONS
%  theta0:     uniform [-180, 180] deg  (full range of initial pitch angles)
%  dphi0:      uniform [-2, 2] deg/s    (initial pitch rate uncertainty)
%  rho_scale:  uniform [0.85, 1.15]     (±15% density uncertainty)
% =========================================================================

N_mc = 500; rng(42);

theta0_samples  = deg2rad(-180 + 360*rand(N_mc,1));   % full angular range
dphi0_samples   = deg2rad(-2   +   4*rand(N_mc,1));   % ±2 deg/s
rho_samples     = 0.85 + 0.30*rand(N_mc,1);           % ±15% density

t_vec = linspace(0, T_orb, 3000)';
Cd_mc = zeros(N_mc,1);
opts  = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);

fprintf('Running %d MC samples (1 orbit each)...\n', N_mc);
for mi = 1:N_mc
    p.rho_scale = rho_samples(mi);
    set_sails(sails_c);

    [t_out, X_out] = ode113(@(t,X) attitude_ode(t,X,p), t_vec, ...
                             [theta0_samples(mi); dphi0_samples(mi)], opts);

    Cd_mc(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sails_c, d_theta_deg);

    if mod(mi,50) == 0
        fprintf('  %d/%d | Cd_bar = %.4f\n', mi, N_mc, Cd_mc(mi));
    end
end

save('mc_results_new.mat', 'Cd_mc', 'theta0_samples', 'dphi0_samples', ...
     'rho_samples', 'N_mc', 'phi_sel', 'alt', 'accom_str');
fprintf('\nSaved mc_results_new.mat\n');

%% =========================================================================
%  STATISTICS
% =========================================================================

fprintf('\n=== MC Results (phi=%d, alpha=%s, h=%d km) ===\n', phi_sel, accom_str, alt)
fprintf('N            = %d\n',  N_mc)
fprintf('Mean Cd_bar  = %.6f\n', mean(Cd_mc))
fprintf('Std          = %.6f\n', std(Cd_mc))
fprintf('CV           = %.4f%%\n', std(Cd_mc)/mean(Cd_mc)*100)
fprintf('95%% CI       = [%.6f, %.6f]\n', ...
    mean(Cd_mc) - 1.96*std(Cd_mc)/sqrt(N_mc), ...
    mean(Cd_mc) + 1.96*std(Cd_mc)/sqrt(N_mc))
fprintf('Min / Max    = %.6f / %.6f\n', min(Cd_mc), max(Cd_mc))
fprintf('Range        = %.6f\n', max(Cd_mc) - min(Cd_mc))
fprintf('IQR          = %.6f\n', iqr(Cd_mc))

%% Spearman rank correlations
rank_theta0 = tiedrank(theta0_samples);
rank_dphi0  = tiedrank(dphi0_samples);
rank_rho    = tiedrank(rho_samples);
rank_CD     = tiedrank(Cd_mc);

tmp = corrcoef(rank_theta0, rank_CD); r_theta0_sp = tmp(1,2);
tmp = corrcoef(rank_dphi0,  rank_CD); r_dphi0_sp  = tmp(1,2);
tmp = corrcoef(rank_rho,    rank_CD); r_rho_sp    = tmp(1,2);

fprintf('\nSpearman rank correlations with Cd_bar:\n')
fprintf('  theta0    : rho = %.4f\n', r_theta0_sp)
fprintf('  dphi0     : rho = %.4f\n', r_dphi0_sp)
fprintf('  rho_scale : rho = %.4f\n', r_rho_sp)

%% Pearson correlations
tmp = corrcoef(theta0_samples, Cd_mc); r_theta0 = tmp(1,2);
tmp = corrcoef(dphi0_samples,  Cd_mc); r_dphi0  = tmp(1,2);
tmp = corrcoef(rho_samples,    Cd_mc); r_rho    = tmp(1,2);

fprintf('\nPearson correlations with Cd_bar:\n')
fprintf('  theta0    : r = %.4f\n', r_theta0)
fprintf('  dphi0     : r = %.4f\n', r_dphi0)
fprintf('  rho_scale : r = %.4f\n', r_rho)

%% Quadratic fit vs dphi0
p2 = polyfit(rad2deg(dphi0_samples), Cd_mc, 2);
fprintf('\nQuadratic fit Cd_bar = a*dphi0^2 + b*dphi0 + c:\n')
fprintf('  a = %.6f\n  b = %.6f\n  c = %.6f\n', p2(1), p2(2), p2(3))
fprintf('  Vertex at dphi0 = %.4f deg/s\n', -p2(2)/(2*p2(1)))
fprintf('  Cd_bar at vertex = %.6f\n', polyval(p2, -p2(2)/(2*p2(1))))

%% Percentiles
fprintf('\nPercentiles of Cd_bar:\n')
pcts = [5 10 25 50 75 90 95];
vals = prctile(Cd_mc, pcts);
for i = 1:length(pcts)
    fprintf('  P%02d = %.6f\n', pcts(i), vals(i))
end

%% =========================================================================
%  PLOTS
% =========================================================================

alt_color   = [0.10 0.72 0.40];   % 450 km green
cmap_rho    = parula(N_mc);

figure('Position',[50 50 1400 900]);

%% 1. Histogram
subplot(2,3,1)
histogram(Cd_mc, 30, 'FaceColor', [0.20 0.45 0.85], 'EdgeColor','none', 'Normalization','pdf')
hold on
x_fit = linspace(min(Cd_mc), max(Cd_mc), 200);
plot(x_fit, normpdf(x_fit, mean(Cd_mc), std(Cd_mc)), 'r-', 'LineWidth', 2)
xline(mean(Cd_mc),        'r-',  'LineWidth', 1.5, 'DisplayName','Mean')
xline(prctile(Cd_mc, 5),  'k--', 'LineWidth', 1.2, 'DisplayName','5th/95th pct')
xline(prctile(Cd_mc, 95), 'k--', 'LineWidth', 1.2, 'HandleVisibility','off')
xlabel('$\bar{C}_D$','Interpreter','latex')
ylabel('PDF','Interpreter','latex')
title('MC Distribution','Interpreter','latex')
legend({'Samples','Normal fit','Mean','5th/95th pct'},'Interpreter','latex','Location','best')
grid on; box on

%% 2. Cd vs theta0
subplot(2,3,2)
scatter(rad2deg(theta0_samples), Cd_mc, 15, rho_samples, 'filled')
cb = colorbar; cb.Label.String = '$s_\rho$'; cb.Label.Interpreter = 'latex';
xlabel('$\theta_0$ (deg)','Interpreter','latex')
ylabel('$\bar{C}_D$','Interpreter','latex')
title('$\bar{C}_D$ vs $\theta_0$','Interpreter','latex')
grid on; box on

%% 3. Cd vs dphi0 with quadratic fit
subplot(2,3,3)
scatter(rad2deg(dphi0_samples), Cd_mc, 15, rho_samples, 'filled')
hold on
x_q = linspace(-2, 2, 200);
plot(x_q, polyval(p2, x_q), 'r-', 'LineWidth', 2)
cb = colorbar; cb.Label.String = '$s_\rho$'; cb.Label.Interpreter = 'latex';
xlabel('$\dot{\theta}_0$ (deg s$^{-1}$)','Interpreter','latex')
ylabel('$\bar{C}_D$','Interpreter','latex')
title('$\bar{C}_D$ vs $\dot{\theta}_0$','Interpreter','latex')
legend({'Samples','Quadratic fit'},'Interpreter','latex','Location','best')
grid on; box on

%% 4. Cd vs rho_scale
subplot(2,3,4)
scatter(rho_samples, Cd_mc, 15, rad2deg(theta0_samples), 'filled')
cb = colorbar; cb.Label.String = '$\theta_0$ (deg)'; cb.Label.Interpreter = 'latex';
xlabel('$s_\rho$','Interpreter','latex')
ylabel('$\bar{C}_D$','Interpreter','latex')
title('$\bar{C}_D$ vs $s_\rho$','Interpreter','latex')
grid on; box on

%% 5. Spearman bar chart
subplot(2,3,5)
spearman = [r_theta0_sp, r_dphi0_sp, r_rho_sp];
pearson  = [r_theta0,    r_dphi0,    r_rho];
b = bar([pearson; spearman]', 'grouped');
b(1).FaceColor = [0.20 0.45 0.85];
b(2).FaceColor = [0.85 0.30 0.10];
set(gca, 'XTickLabel', {'$\theta_0$','$\dot{\theta}_0$','$s_\rho$'}, ...
    'TickLabelInterpreter','latex')
ylabel('Correlation coefficient','Interpreter','latex')
title('Pearson vs Spearman','Interpreter','latex')
legend({'Pearson','Spearman'},'Interpreter','latex','Location','best')
yline(0,'k-','LineWidth',1)
grid on; box on

%% 6. Empirical CDF
subplot(2,3,6)
CD_sort = sort(Cd_mc);
cdf_vals = (1:N_mc)/N_mc;
plot(CD_sort, cdf_vals, 'b-', 'LineWidth', 2)
hold on
xline(mean(Cd_mc),        'r-',  'LineWidth', 1.5, 'DisplayName','Mean')
xline(prctile(Cd_mc, 5),  'k--', 'LineWidth', 1.2, 'DisplayName','5th/95th pct')
xline(prctile(Cd_mc, 95), 'k--', 'LineWidth', 1.2, 'HandleVisibility','off')
xlabel('$\bar{C}_D$','Interpreter','latex')
ylabel('CDF','Interpreter','latex')
title('Empirical CDF','Interpreter','latex')
legend('Interpreter','latex','Location','best')
grid on; box on

sgtitle(sprintf('Monte Carlo Sensitivity Analysis ($\\phi=%d^\\circ$, $\\alpha=0.75$, $h=%d$~km, $N=%d$)', ...
        phi_sel, alt, N_mc), 'Interpreter','latex','FontSize',12)

exportgraphics(gcf,'mc_sensitivity_new.png','Resolution',300)
fprintf('Saved mc_sensitivity_new.png\n')

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function CD_bar = compute_CDbar(t_out, phi_deg_t, sails_c, d_theta_deg)
    edges   = (-180 : d_theta_deg : 180);
    centres = edges(1:end-1) + d_theta_deg/2;
    n_bins  = length(centres);
    dt_k    = diff(t_out); dt_k(end+1) = dt_k(end);
    Delta_t = zeros(1, n_bins);
    for k = 1:length(t_out)
        wrapped = wrapTo180(phi_deg_t(k));
        bi = find(edges <= wrapped, 1, 'last');
        if ~isempty(bi) && bi <= n_bins
            Delta_t(bi) = Delta_t(bi) + dt_k(k);
        end
    end
    pdf_SAM = Delta_t / (sum(dt_k) * d_theta_deg);
    Cd_bins = zeros(1, n_bins);
    for bi = 1:n_bins
        Cd_bins(bi) = interp1(sails_c(1).alpha, sails_c(1).Cd, ...
                              centres(bi), 'pchip', 'extrap');
    end
    CD_bar = sum(Cd_bins .* pdf_SAM) * d_theta_deg;
end

function dX = attitude_ode(t, X, p)
    phi  = X(1); dphi = X(2);
    [~, tau_aero] = compute_rho_torque(t, phi, p);
    n_orb  = 2*pi/p.T_orb;
    tau_gg = -1.5*n_orb^2*p.Iyy*sin(2*phi);
    dX     = [dphi; (tau_aero + tau_gg)/p.Iyy];
end

function [rho, tau_total] = compute_rho_torque(t, phi, p)
    n_orb = 2*pi/p.T_orb;
    nu    = mod(n_orb*t, 2*pi);
    [r_eci, ~] = kep2eci_circ(p.a_orb, p.inc, p.RAAN, nu);
    utc_dt  = datetime(p.year,p.month,p.day,p.hour,0,0) + seconds(t);
    utc_vec = [utc_dt.Year utc_dt.Month utc_dt.Day ...
               utc_dt.Hour utc_dt.Minute utc_dt.Second];
    lla   = eci2lla(r_eci', utc_vec);
    doy   = day(utc_dt,'dayofyear');
    utc_s = hour(utc_dt)*3600 + minute(utc_dt)*60 + second(utc_dt);
    [~, rho_arr] = atmosnrlmsise00(lla(3), lla(1), lla(2), ...
                   p.year, doy, utc_s, p.f107a, p.f107, p.ap);
    rho = rho_arr(6) * p.rho_scale;
    q   = 0.5*rho*p.v_circ^2;
    sails     = set_sails();
    phi_deg   = rad2deg(phi);
    tau_total = 0;
    for s = 1:length(sails)
        Cm_s      = interp1(sails(s).alpha, sails(s).Cm_pitch, ...
                            phi_deg, 'pchip', 'extrap');
        tau_total = tau_total + q*sails(s).A_ref*sails(s).L_ref*Cm_s;
    end
end

function sails = set_sails(new_sails)
    persistent stored_sails;
    if nargin > 0; stored_sails = new_sails; end
    sails = stored_sails;
end

function [r, v] = kep2eci_circ(a, inc, RAAN, nu)
    mu_ = 3.986004418e14;
    r_pf = a*[cos(nu);sin(nu);0];
    v_pf = sqrt(mu_/a)*[-sin(nu);cos(nu);0];
    cO=cos(RAAN); sO=sin(RAAN); ci=cos(inc); si=sin(inc);
    R = [cO,-sO*ci,sO*si; sO,cO*ci,-cO*si; 0,si,ci];
    r = R*r_pf; v = R*v_pf;
end
