P_ROOT = fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(P_ROOT,'lib'))); cd(P_ROOT);  % repo root + lib on path
clear; clc; close all; clear functions
 
load('sam_results.mat')
load('inertia_tensors.mat','inertiaDB');
 
%% Constants
mu = 3.986004418e14;
Re = 6.3781e6;
 
% FIX 3: January epoch, 4th year of typical 11-year solar cycle
% Year 4 corresponds to rising phase toward solar maximum
% F10.7 ~150 sfu is reasonable; Roberts uses a "typical" cycle so keep 150
% but shift to January for correct density distribution
year_ep  = 2025;
month_ep = 1;       % <-- FIX 3: January (was June)
day_ep   = 15;      % mid-January
 
f107a = 150; f107 = 150; ap = 15;
 
ode_reltol = 1e-6;
ode_abstol = 1e-8;
 
%% Simulation duration (MULTI-DAY)
t_final = 2 * 24 * 3600;   % 2 days
tspan   = [0 t_final];
 
%% Configuration
phi_sel = 75;
alt     = 450;
accom   = '0p75';
 
adbsat_base = fullfile(project_root(), 'adbsat_processed');
 
I_idx   = find(inertiaDB.angles_deg == phi_sel);
Iyy_val = inertiaDB.I(2,2,I_idx);
g_idx   = find(ref_geom.phi_deg == phi_sel);
 
fname = fullfile(adbsat_base, sprintf('%dkm',alt), ...
                 sprintf('%ddeg_CLL_accom_%s.mat',phi_sel,accom));
 
raw = load(fname);
 
Cm_sort = raw.aedb.aero.Cm_BY(:);
Cm_sort(isnan(Cm_sort)) = 0;
 
Cd_sort = -raw.aedb.aero.Cf_wX(:);
Cd_sort(isnan(Cd_sort)) = 0;
 
%% Sail struct
clear sails_c
sails_c(1).angle   = phi_sel;
sails_c(1).alpha   = (-179.5:1:179.5)';
sails_c(1).Cm_pitch = Cm_sort;
sails_c(1).Cd      = Cd_sort;
sails_c(1).Iyy     = Iyy_val;
sails_c(1).A_ref   = ref_geom.A_ref(g_idx);
sails_c(1).L_ref   = ref_geom.L_ref(g_idx);
 
%% Orbit
a_orb  = Re + alt*1e3;
T_orb  = 2*pi*sqrt(a_orb^3/mu);
v_circ = sqrt(mu/a_orb);
 
%% Initial conditions matching Roberts & Harkness 2007
% phi0 = 2.8 rad (~160 deg), dphi0 = 1 rev/orbit in inertial frame
phi0_test  = deg2rad(160);
dphi0_test = 2*pi / T_orb;   % 1 rev/orbit (rad/s)
 
%% Run BOTH cases
cases = [0, 12];  % midnight and midday
 
results = struct;
 
for k = 1:length(cases)
 
    hour_utc = cases(k);
 
    p.a_orb  = a_orb;
    p.inc    = deg2rad(0);   % FIX 1: equatorial orbit (was 98 deg SSO)
    p.RAAN   = deg2rad(90);
    p.T_orb  = T_orb;
    p.v_circ = v_circ;
 
    p.year  = year_ep;
    p.month = month_ep;
    p.day   = day_ep;
    p.hour  = hour_utc;
 
    p.f107a = f107a;
    p.f107  = f107;
    p.ap    = ap;
 
    p.Iyy   = Iyy_val;
    p.sails = sails_c;
    p.n_sails = 1;
 
    set_sails(sails_c);
 
    opts = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
 
    [t_out, X_out] = ode113(@(t,X) attitude_ode(t,X,p), ...
                            tspan, [phi0_test; dphi0_test], opts);
 
    results(k).hour  = hour_utc;
    results(k).t     = t_out;
    results(k).theta = rad2deg(X_out(:,1));
    results(k).omega = rad2deg(X_out(:,2));
 
    % FIX 2: wrapped angle of attack for plotting (mod into [-180, 180])
    results(k).theta_wrapped = wrapTo180(results(k).theta);
 
end
 
%% ===== PLOTS =====
 
for k = 1:length(results)
 
    t             = results(k).t;
    theta_raw     = results(k).theta;
    theta_wrapped = results(k).theta_wrapped;
    omega         = results(k).omega;
 
    hour_label = results(k).hour;
 
    %% Full duration: raw (accumulated) angle and rate
    figure('Name', sprintf('Full duration, hour=%d', hour_label));
    subplot(1,2,1)
    plot(t/3600, theta_raw, 'LineWidth', 1)
    xlabel('Time (hours)', 'Interpreter', 'latex')
    ylabel('$\theta$ (deg, accumulated)', 'Interpreter', 'latex')
    title(sprintf('Pitch -- accumulated (hour = %d)', hour_label), 'Interpreter', 'latex')
    grid on
 
    subplot(1,2,2)
    plot(t/3600, omega, 'LineWidth', 1)
    xlabel('Time (hours)', 'Interpreter', 'latex')
    ylabel('$\dot{\theta}$ (deg s$^{-1}$)', 'Interpreter', 'latex')
    title(sprintf('Pitch rate (hour = %d)', hour_label), 'Interpreter', 'latex')
    grid on
 
    %% Full duration: WRAPPED angle of attack (comparable to Roberts Fig. 8)
    figure('Name', sprintf('Wrapped AoA, hour=%d', hour_label));
    plot(t/3600, theta_wrapped, 'LineWidth', 0.8)
    xlabel('Time (hours)', 'Interpreter', 'latex')
    ylabel('$\theta_\mathrm{wrapped}$ (deg)', 'Interpreter', 'latex')
    title(sprintf('Angle of attack -- wrapped to $[-180, 180]$ (hour = %d)', hour_label), ...
          'Interpreter', 'latex')
    ylim([-180 180])
    grid on
 
    %% 1-hour zoom: wrapped angle, first hour (comparable to Roberts detail panels)
    idx_zoom = t <= 3600;
 
    figure('Name', sprintf('1-hr zoom wrapped, hour=%d', hour_label));
    plot(t(idx_zoom)/60, theta_wrapped(idx_zoom), 'LineWidth', 1)
    xlabel('Time (min)', 'Interpreter', 'latex')
    ylabel('$\theta_\mathrm{wrapped}$ (deg)', 'Interpreter', 'latex')
    title(sprintf('Zoom -- first hour, wrapped (hour = %d)', hour_label), ...
          'Interpreter', 'latex')
    grid on
 
end
 
%% ===== FUNCTIONS =====
 
function dX = attitude_ode(t,X,p)
 
    phi  = X(1);
    dphi = X(2);
 
    [~, tau_aero] = compute_rho_torque(t, phi, p);
 
    n_orb  = 2*pi / p.T_orb;
    tau_gg = -1.5 * n_orb^2 * p.Iyy * sin(2*phi);
 
    dX = [dphi; (tau_aero + tau_gg) / p.Iyy];
 
end
 
function [rho, tau_total] = compute_rho_torque(t, phi, p)
 
    n_orb = 2*pi / p.T_orb;
    nu    = mod(n_orb*t, 2*pi);
 
    [r_eci, ~] = kep2eci_circ(p.a_orb, p.inc, p.RAAN, nu);
 
    %% Proper time handling
    utc_dt = datetime(p.year, p.month, p.day, p.hour, 0, 0) + seconds(t);
 
    utc_vec = [utc_dt.Year, utc_dt.Month, utc_dt.Day, ...
               utc_dt.Hour, utc_dt.Minute, utc_dt.Second];
 
    lla = eci2lla(r_eci', utc_vec);
 
    doy   = day(utc_dt, 'dayofyear');
    utc_s = hour(utc_dt)*3600 + minute(utc_dt)*60 + second(utc_dt);
 
    [~, rho_arr] = atmosnrlmsise00( ...
        lla(3), lla(1), lla(2), ...
        p.year, doy, utc_s, ...
        p.f107a, p.f107, p.ap);
 
    rho = rho_arr(6);
 
    q = 0.5 * rho * p.v_circ^2;
 
    sails = set_sails();
 
    phi_deg = rad2deg(phi);
 
    tau_total = 0;
 
    for s = 1:length(sails)
 
        Cm_s = interp1(sails(s).alpha, sails(s).Cm_pitch, ...
                       phi_deg, 'pchip', 'extrap');
 
        tau_total = tau_total + q * sails(s).A_ref * ...
                                   sails(s).L_ref * Cm_s;
    end
 
end
 
function sails = set_sails(new_sails)
    persistent stored_sails;
    if nargin > 0
        stored_sails = new_sails;
    end
    sails = stored_sails;
end
 
function [r, v] = kep2eci_circ(a, inc, RAAN, nu)
 
    mu_ = 3.986004418e14;
 
    r_pf = a * [cos(nu); sin(nu); 0];
    v_pf = sqrt(mu_/a) * [-sin(nu); cos(nu); 0];
 
    cO = cos(RAAN); sO = sin(RAAN);
    ci = cos(inc);  si = sin(inc);
 
    R = [cO, -sO*ci,  sO*si;
         sO,  cO*ci, -cO*si;
          0,  si,     ci];
 
    r = R * r_pf;
    v = R * v_pf;
 
end