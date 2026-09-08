%% =========================================================================
%  montecarlo.m  --  Monte Carlo sensitivity of the orientation-averaged
%                    drag coefficient Cd_bar
%
%  Samples the initial pitch angle theta0, the initial pitch rate dtheta0
%  and an atmospheric-density scale factor s_rho simultaneously, propagates
%  the 1-DOF pitch model for one orbit per sample, and reduces each run to
%  Cd_bar via the SAM residence-time average. Reports Pearson / Spearman
%  correlations, percentiles and a quadratic fit vs dtheta0, and a
%  six-panel summary figure.
%
%  Everything is controlled by the CFG block below (geometry, altitude,
%  accommodation, sample count, sampling ranges, solver). Set
%  CFG.per_phi_barchart = true to also sweep phi = 45:5:85 with MC error
%  bars (writes mc_barchart.mat / .png).
%
%  Outputs: mc_results.mat, mc_workspace.mat, mc_extended_analysis.png
%  Requires adbsat_processed/ + ref_geometry.mat + inertia_tensors.mat.
%
%  (Consolidates the former montecarlo_multiparam / montecarlo_sensitivity /
%  montecarlo_barchart scripts.)
% =========================================================================
clear; clc; close all; clear functions

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'lib')));
addpath(genpath(fullfile(here, 'adbsat_processed')));
cd(here);

%% ============================ CFG ====================================
CFG.phi_deg          = 75;             % apex half-angle for the main run
CFG.alt_km           = 450;
CFG.accom_str        = '0p75';
CFG.N                = 200;            % number of MC samples
CFG.seed             = 42;
CFG.theta0_range_deg   = [-10  10];    % uniform
CFG.dtheta0_range_dps  = [ -2   2];    % uniform  (deg/s)
CFG.rho_scale_range    = [0.85 1.15];  % uniform
CFG.solver             = 'ode45';      % 'ode45' or 'ode113'
CFG.n_steps            = 3000;
CFG.ode_reltol         = 1e-8;
CFG.ode_abstol         = 1e-10;
CFG.d_theta_deg        = 2.0;
CFG.per_phi_barchart   = false;        % also do the phi = 45:5:85 bar chart
CFG.barchart_phi_vec   = 45:5:85;
CFG.barchart_N         = 100;

% epoch / space weather (matches ROADM)
CFG.year = 2025; CFG.month = 6; CFG.day = 21; CFG.hour = 10;
CFG.f107a = 150; CFG.f107 = 150; CFG.ap = 15;
CFG.inc_deg = 98; CFG.RAAN_deg = 90;
%% ====================================================================

mu = 3.986004418e14; Re = 6.3781e6;

load('ref_geometry.mat',    'ref_geom');
load('inertia_tensors.mat', 'inertiaDB');
adbsat_base = fullfile(here, 'adbsat_processed');

a_orb  = Re + CFG.alt_km*1e3;
T_orb  = 2*pi*sqrt(a_orb^3/mu);
v_circ = sqrt(mu/a_orb);
doy0   = day(datetime(CFG.year, CFG.month, CFG.day), 'dayofyear');

[sail, p] = build_case(CFG.phi_deg, CFG.alt_km, CFG.accom_str, ...
                       adbsat_base, ref_geom, inertiaDB, CFG, ...
                       a_orb, T_orb, v_circ, doy0);

odefun = str2func(CFG.solver);
opts   = odeset('RelTol', CFG.ode_reltol, 'AbsTol', CFG.ode_abstol);
t_vec  = linspace(0, T_orb, CFG.n_steps)';

%% ----- sample --------------------------------------------------------
rng(CFG.seed);
N = CFG.N;
theta0_samples = deg2rad(CFG.theta0_range_deg(1) + diff(CFG.theta0_range_deg)*rand(N,1));
dphi0_samples  = deg2rad(CFG.dtheta0_range_dps(1) + diff(CFG.dtheta0_range_dps)*rand(N,1));
rho_samples    = CFG.rho_scale_range(1) + diff(CFG.rho_scale_range)*rand(N,1);

CD_mc = zeros(N,1);
fprintf('Running %d MC samples (phi=%d, alpha=%s, h=%d km)...\n', ...
        N, CFG.phi_deg, CFG.accom_str, CFG.alt_km);
for mi = 1:N
    p.rho_scale = rho_samples(mi);
    set_sails(sail);
    [t_out, X_out] = odefun(@(t,X) attitude_ode(t,X,p), t_vec, ...
                            [theta0_samples(mi); dphi0_samples(mi)], opts);
    CD_mc(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail, CFG.d_theta_deg);
    if mod(mi,20) == 0
        fprintf('  %d/%d | Cd_bar=%.4f\n', mi, N, CD_mc(mi));
    end
end

save('mc_workspace.mat');   % validation/ scripts load this
fprintf('Saved mc_workspace.mat\n');

%% ----- statistics ---------------------------------------------------
fprintf('\n=== MC Results ===\n');
fprintf('N            = %d\n', N);
fprintf('Mean Cd_bar  = %.6f\n', mean(CD_mc));
fprintf('Std          = %.6f\n', std(CD_mc));
fprintf('CV           = %.4f%%\n', std(CD_mc)/mean(CD_mc)*100);
fprintf('95%% CI       = [%.6f, %.6f]\n', ...
    mean(CD_mc)-1.96*std(CD_mc)/sqrt(N), mean(CD_mc)+1.96*std(CD_mc)/sqrt(N));
fprintf('Min / Max    = %.6f / %.6f\n', min(CD_mc), max(CD_mc));

tmp = corrcoef(theta0_samples, CD_mc); r_phi0  = tmp(1,2);
tmp = corrcoef(dphi0_samples,  CD_mc); r_dphi0 = tmp(1,2);
tmp = corrcoef(rho_samples,    CD_mc); r_rho   = tmp(1,2);
tmp = corrcoef(tiedrank(theta0_samples), tiedrank(CD_mc)); r_phi0_sp  = tmp(1,2);
tmp = corrcoef(tiedrank(dphi0_samples),  tiedrank(CD_mc)); r_dphi0_sp = tmp(1,2);
tmp = corrcoef(tiedrank(rho_samples),    tiedrank(CD_mc)); r_rho_sp   = tmp(1,2);
fprintf('\nPearson  : theta0 %.4f | dtheta0 %.4f | s_rho %.4f\n', r_phi0, r_dphi0, r_rho);
fprintf('Spearman : theta0 %.4f | dtheta0 %.4f | s_rho %.4f\n', r_phi0_sp, r_dphi0_sp, r_rho_sp);

pcts = [5 10 25 50 75 90 95];
vals = prctile(CD_mc, pcts);
fprintf('\nPercentiles of Cd_bar:\n');
for i = 1:numel(pcts); fprintf('  P%02d = %.6f\n', pcts(i), vals(i)); end

p2 = polyfit(rad2deg(dphi0_samples), CD_mc, 2);
fprintf('\nQuadratic fit vs dtheta0: a=%.6f b=%.6f c=%.6f | vertex %.4f deg/s\n', ...
        p2(1), p2(2), p2(3), -p2(2)/(2*p2(1)));

%% ----- figure -----------------------------------------------------
figure('Position',[50 50 1400 900]);

subplot(2,3,1);
histogram(CD_mc, 30, 'FaceColor',[0.2 0.4 0.8], 'EdgeColor','none', 'Normalization','pdf');
hold on;
x_fit = linspace(min(CD_mc), max(CD_mc), 200);
plot(x_fit, normpdf(x_fit, mean(CD_mc), std(CD_mc)), 'r-', 'LineWidth', 2);
xline(mean(CD_mc), 'r-', 'LineWidth', 2);
xline(prctile(CD_mc,5),  'k--', 'LineWidth', 1.2);
xline(prctile(CD_mc,95), 'k--', 'LineWidth', 1.2);
xlabel('$\bar{C}_D$','Interpreter','latex'); ylabel('PDF','Interpreter','latex');
title('MC distribution','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',11,'Box','on');

subplot(2,3,2);
scatter(rad2deg(theta0_samples), CD_mc, 20, rho_samples, 'filled');
cb = colorbar; cb.Label.String = '$s_\rho$'; cb.Label.Interpreter = 'latex';
xlabel('$\theta_0$ (deg)','Interpreter','latex'); ylabel('$\bar{C}_D$','Interpreter','latex');
title('$\bar{C}_D$ vs $\theta_0$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',11,'Box','on');

subplot(2,3,3);
scatter(rad2deg(dphi0_samples), CD_mc, 20, rho_samples, 'filled'); hold on;
x_q = linspace(CFG.dtheta0_range_dps(1), CFG.dtheta0_range_dps(2), 200);
plot(x_q, polyval(p2, x_q), 'r-', 'LineWidth', 2);
cb = colorbar; cb.Label.String = '$s_\rho$'; cb.Label.Interpreter = 'latex';
xlabel('$\dot{\theta}_0$ (deg/s)','Interpreter','latex'); ylabel('$\bar{C}_D$','Interpreter','latex');
title('$\bar{C}_D$ vs $\dot{\theta}_0$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',11,'Box','on');

subplot(2,3,4);
scatter(rho_samples, CD_mc, 20, rad2deg(theta0_samples), 'filled');
cb = colorbar; cb.Label.String = '$\theta_0$ (deg)'; cb.Label.Interpreter = 'latex';
xlabel('$s_\rho$','Interpreter','latex'); ylabel('$\bar{C}_D$','Interpreter','latex');
title('$\bar{C}_D$ vs $s_\rho$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',11,'Box','on');

subplot(2,3,5);
b = bar([[r_phi0 r_dphi0 r_rho]; [r_phi0_sp r_dphi0_sp r_rho_sp]]', 'grouped');
b(1).FaceColor = [0.2 0.4 0.8]; b(2).FaceColor = [0.8 0.3 0.2];
set(gca,'XTickLabel',{'$\theta_0$','$\dot{\theta}_0$','$s_\rho$'}, ...
    'TickLabelInterpreter','latex','FontSize',11,'Box','on');
ylabel('Correlation','Interpreter','latex'); yline(0,'k-');
legend({'Pearson','Spearman'},'Interpreter','latex','Location','best');
title('Pearson vs Spearman','Interpreter','latex');

subplot(2,3,6);
plot(sort(CD_mc), (1:N)/N, 'b-', 'LineWidth', 2); hold on;
xline(mean(CD_mc), 'r-', 'LineWidth', 1.5);
xline(prctile(CD_mc,5),  'k--', 'LineWidth', 1.2);
xline(prctile(CD_mc,95), 'k--', 'LineWidth', 1.2);
xlabel('$\bar{C}_D$','Interpreter','latex'); ylabel('CDF','Interpreter','latex');
title('Empirical CDF','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',11,'Box','on');

sgtitle(sprintf(['Monte Carlo sensitivity ($\\phi=%d^\\circ$, $\\alpha=%s$, ' ...
        '$h=%d$~km, $N=%d$)'], CFG.phi_deg, strrep(CFG.accom_str,'p','.'), ...
        CFG.alt_km, N), 'Interpreter','latex','FontSize',13);
exportgraphics(gcf, 'mc_extended_analysis.png', 'Resolution', 300);
fprintf('Saved mc_extended_analysis.png\n');

save('mc_results.mat', 'CD_mc', 'theta0_samples', 'dphi0_samples', 'rho_samples', ...
     'N', 'r_phi0', 'r_dphi0', 'r_rho', 'r_phi0_sp', 'r_dphi0_sp', 'r_rho_sp', ...
     'p2', 'pcts', 'vals', 'CFG');
fprintf('Saved mc_results.mat\n');

%% ----- optional: per-phi bar chart -------------------------------
if CFG.per_phi_barchart
    pv  = CFG.barchart_phi_vec;
    Nb  = CFG.barchart_N;
    rng(CFG.seed);
    dphi0_b = deg2rad(CFG.dtheta0_range_dps(1) + diff(CFG.dtheta0_range_dps)*rand(Nb,1));
    rho_b   = CFG.rho_scale_range(1) + diff(CFG.rho_scale_range)*rand(Nb,1);
    theta0_fixed = deg2rad(5);

    CD_nom = zeros(numel(pv),1); CD_lo = CD_nom; CD_hi = CD_nom;
    CD_mn  = CD_nom;             CD_sd = CD_nom;

    for k = 1:numel(pv)
        [sail_k, p_k] = build_case(pv(k), CFG.alt_km, CFG.accom_str, ...
                                   adbsat_base, ref_geom, inertiaDB, CFG, ...
                                   a_orb, T_orb, v_circ, doy0);
        set_sails(sail_k);
        p_k.rho_scale = 1.0;
        [t_out,X_out] = odefun(@(t,X) attitude_ode(t,X,p_k), t_vec, [theta0_fixed;0], opts);
        CD_nom(k) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_k, CFG.d_theta_deg);

        cd_k = zeros(Nb,1);
        for mi = 1:Nb
            p_k.rho_scale = rho_b(mi);
            set_sails(sail_k);
            [t_out,X_out] = odefun(@(t,X) attitude_ode(t,X,p_k), t_vec, ...
                                   [theta0_fixed; dphi0_b(mi)], opts);
            cd_k(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_k, CFG.d_theta_deg);
        end
        CD_lo(k) = min(cd_k); CD_hi(k) = max(cd_k);
        CD_mn(k) = mean(cd_k); CD_sd(k) = std(cd_k);
        fprintf('  phi=%d | nom=%.4f | [%.4f, %.4f]\n', pv(k), CD_nom(k), CD_lo(k), CD_hi(k));
    end

    save('mc_barchart.mat', 'pv', 'CD_nom', 'CD_lo', 'CD_hi', 'CD_mn', 'CD_sd');
    fig = figure('Position',[100 100 800 500]);
    bar(pv, CD_nom, 0.6, 'FaceColor',[0.2 0.4 0.8], 'EdgeColor','none'); hold on;
    errorbar(pv, CD_nom, CD_nom-CD_lo, CD_hi-CD_nom, 'k.', 'LineWidth', 1.5, 'CapSize', 8);
    xlabel('$\phi$ (deg)','Interpreter','latex'); ylabel('$\bar{C}_D$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex','FontSize',12,'Box','on','XTick',pv,'YGrid','on');
    exportgraphics(fig, 'mc_barchart.png', 'Resolution', 300);
    fprintf('Saved mc_barchart.mat / mc_barchart.png\n');
end

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================
function [sail, p] = build_case(phi_deg, alt_km, accom_str, adbsat_base, ...
                                ref_geom, inertiaDB, CFG, a_orb, T_orb, v_circ, doy0)
    fname = fullfile(adbsat_base, sprintf('%dkm', alt_km), ...
                     sprintf('%ddeg_CLL_accom_%s.mat', phi_deg, accom_str));
    if ~isfile(fname); error('ADM file not found: %s', fname); end
    raw = load(fname);
    Cm  = raw.aedb.aero.Cm_BY(:);  Cm(isnan(Cm)) = 0;
    Cd  = -raw.aedb.aero.Cf_wX(:); Cd(isnan(Cd)) = 0;

    I_idx = find(inertiaDB.angles_deg == phi_deg);
    g_idx = find(ref_geom.phi_deg == phi_deg);

    sail(1).angle    = phi_deg;
    sail(1).alpha    = (-179.5:1:179.5)';
    sail(1).Cm_pitch = Cm;
    sail(1).Cd       = Cd;
    sail(1).Iyy      = inertiaDB.I(2,2,I_idx);
    sail(1).A_ref    = ref_geom.A_ref(g_idx);
    sail(1).L_ref    = ref_geom.L_ref(g_idx);

    p.a_orb  = a_orb;                p.inc   = deg2rad(CFG.inc_deg);
    p.RAAN   = deg2rad(CFG.RAAN_deg); p.T_orb = T_orb;
    p.v_circ = v_circ;               p.doy0  = doy0;
    p.year   = CFG.year;             p.month = CFG.month;
    p.day    = CFG.day;              p.hour  = CFG.hour;
    p.f107a  = CFG.f107a;            p.f107  = CFG.f107;   p.ap = CFG.ap;
    p.Iyy    = sail(1).Iyy;          p.sails = sail;       p.n_sails = 1;
    p.rho_scale = 1.0;
end

function CD_bar = compute_CDbar(t_out, phi_deg_t, sail, d_theta_deg)
    edges   = (-180:d_theta_deg:180);
    centres = edges(1:end-1) + d_theta_deg/2;
    n_bins  = numel(centres);
    dt_k = diff(t_out); dt_k(end+1) = dt_k(end);
    Delta_t = zeros(1, n_bins);
    for k = 1:numel(t_out)
        bi = find(edges <= phi_deg_t(k), 1, 'last');
        if ~isempty(bi) && bi <= n_bins
            Delta_t(bi) = Delta_t(bi) + dt_k(k);
        end
    end
    pdf_SAM = Delta_t / (sum(dt_k) * d_theta_deg);
    Cd_bins = zeros(1, n_bins);
    for bi = 1:n_bins
        Cd_bins(bi) = interp1(sail(1).alpha, sail(1).Cd, centres(bi), 'pchip', 'extrap');
    end
    CD_bar = sum(Cd_bins .* pdf_SAM) * d_theta_deg;
end

function dX = attitude_ode(t, X, p)
    phi = X(1); dphi = X(2);
    [~, tau_aero] = compute_rho_torque(t, phi, p);
    n_orb  = 2*pi/p.T_orb;
    tau_gg = -1.5*n_orb^2*p.Iyy*sin(2*phi);
    dX = [dphi; (tau_aero + tau_gg)/p.Iyy];
end

function [rho, tau_total] = compute_rho_torque(t, phi, p)
    n_orb = 2*pi/p.T_orb;
    nu    = mod(n_orb*t, 2*pi);
    [r_eci, ~] = kep2eci_circ(p.a_orb, p.inc, p.RAAN, nu);
    utc_dt  = datetime(p.year, p.month, p.day, p.hour, 0, 0) + seconds(t);
    utc_vec = [utc_dt.Year, utc_dt.Month, utc_dt.Day, ...
               utc_dt.Hour, utc_dt.Minute, utc_dt.Second];
    lla = eci2lla(r_eci', utc_vec);
    doy   = p.doy0 + t/86400;
    utc_s = mod(p.hour*3600 + t, 86400);
    [~, rho_arr] = atmosnrlmsise00(lla(3), lla(1), lla(2), ...
                   p.year, doy, utc_s, p.f107a, p.f107, p.ap);
    rho = rho_arr(6) * p.rho_scale;
    q   = 0.5 * rho * p.v_circ^2;
    sails = set_sails();
    phi_deg = rad2deg(phi);
    tau_total = 0;
    for s = 1:numel(sails)
        Cm_s = interp1(sails(s).alpha, sails(s).Cm_pitch, phi_deg, 'pchip', 'extrap');
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
    r_pf = a*[cos(nu); sin(nu); 0];
    v_pf = sqrt(mu_/a)*[-sin(nu); cos(nu); 0];
    cO = cos(RAAN); sO = sin(RAAN);
    ci = cos(inc);  si = sin(inc);
    R = [cO,-sO*ci,sO*si; sO,cO*ci,-cO*si; 0,si,ci];
    r = R*r_pf; v = R*v_pf;
end
