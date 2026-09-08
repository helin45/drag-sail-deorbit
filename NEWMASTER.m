%% =========================================================================
%  DRAG SAIL ATTITUDE DYNAMICS + SAM — FULL PARAMETER SWEEP
%  Equatorial orbit, January epoch, inertiaDB Iyy
%  phi: 45:5:85 deg | accom: 0.50:0.05:1.00 | alt: 350/450/650 km
% =========================================================================
clear; clc; close all;
clear functions
 
%% =========================================================================
%  SECTION 1 — CONSTANTS & SETTINGS
% =========================================================================
 
mu = 3.986004418e14;
Re = 6.3781e6;
 
% Equatorial orbit, January epoch (matched to Roberts & Harkness 2007)
year_ep  = 2025;
month_ep = 1;
day_ep   = 15;
hour_utc = 0;          % midnight deployment as primary case
 
f107a = 150; f107 = 150; ap = 15;
 
ode_reltol = 1e-6;
ode_abstol = 1e-8;
 
% SAM settings
d_theta_deg = 2.0;     % bin width (deg)
 
% Attitude initial conditions (Roberts & Harkness 2007)
phi0  = deg2rad(160);  % 2.8 rad
dphi0 = 0.0;           % will be set per orbit below
 
% Simulation: one full orbit for SAM (attitude converges within ~1 orbit)
% Extend if you want multi-orbit averaging — change multiplier below
n_orb_sim = 1;
 
%% =========================================================================
%  SECTION 2 — LOAD GEOMETRY & INERTIA
% =========================================================================
 
load('ref_geometry.mat',   'ref_geom');
load('inertia_tensors.mat','inertiaDB');
 
adbsat_base = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';
 
%% =========================================================================
%  SECTION 3 — PARAMETER SPACE
% =========================================================================
 
alt_vec    = [350, 450, 650];
phi_vec    = 45:5:85;                % apex half-angles (deg)
accom_vec  = {'0p50','0p55','0p60','0p65','0p70','0p75', ...
              '0p80','0p85','0p90','0p95','1p00'};
accom_vals = 0.50:0.05:1.00;
 
n_alts  = length(alt_vec);
n_phi   = length(phi_vec);
n_accom = length(accom_vec);
 
% Pre-allocate outputs
Cd_bar_3D = zeros(n_alts, n_accom, n_phi);   % [alt x accom x phi]
ref_results = struct();
 
%% =========================================================================
%  SECTION 4 — MAIN SWEEP
% =========================================================================
 
for alt_i = 1:n_alts
    alt     = alt_vec(alt_i);
    alt_dir = fullfile(adbsat_base, sprintf('%dkm', alt));
 
    a_orb  = Re + alt*1e3;
    T_orb  = 2*pi * sqrt(a_orb^3/mu);
    v_circ = sqrt(mu/a_orb);
    t_span = [0, n_orb_sim * T_orb];
 
    fprintf('\n====== Altitude: %d km (T_orb = %.1f min) ======\n', ...
            alt, T_orb/60);
 
    for ac = 1:n_accom
        accom_str = accom_vec{ac};
        fprintf('  Accommodation: %s\n', accom_str);
 
        for si = 1:n_phi
            phi_sel = phi_vec(si);
 
            %% --- Load aero database ---
            fname = fullfile(alt_dir, ...
                    sprintf('%ddeg_CLL_accom_%s.mat', phi_sel, accom_str));
            if ~isfile(fname)
                error('File not found: %s', fname);
            end
            raw  = load(fname);
            aero = raw.aedb.aero;
 
            Cm_sort = aero.Cm_BY(:); Cm_sort(isnan(Cm_sort)) = 0;
            Cd_sort = -aero.Cf_wX(:); Cd_sort(isnan(Cd_sort)) = 0;
 
            %% --- Inertia from database (your model, not Roberts sphere) ---
            I_idx   = find(inertiaDB.angles_deg == phi_sel);
            Iyy_val = inertiaDB.I(2,2,I_idx);
 
            g_idx = find(ref_geom.phi_deg == phi_sel);
 
            %% --- Build sail struct ---
            clear sails_c
            sails_c(1).angle     = phi_sel;
            sails_c(1).alpha     = (-179.5:1:179.5)';
            sails_c(1).Cm_pitch  = Cm_sort;
            sails_c(1).Cd        = Cd_sort;
            sails_c(1).Iyy       = Iyy_val;
            sails_c(1).A_ref     = ref_geom.A_ref(g_idx);
            sails_c(1).L_ref     = ref_geom.L_ref(g_idx);
 
            %% --- ODE parameters ---
            p.a_orb  = a_orb;
            p.inc    = deg2rad(0);       % equatorial
            p.RAAN   = deg2rad(90);
            p.T_orb  = T_orb;
            p.v_circ = v_circ;
 
            p.year   = year_ep;
            p.month  = month_ep;
            p.day    = day_ep;
            p.hour   = hour_utc;
 
            p.f107a  = f107a;
            p.f107   = f107;
            p.ap     = ap;
 
            p.Iyy    = Iyy_val;
            p.sails  = sails_c;
            p.n_sails = 1;
 
            % Initial pitch rate: 1 rev/orbit (Roberts & Harkness 2007)
            dphi0_val = 2*pi / T_orb;
 
            set_sails(sails_c);
 
            opts = odeset('RelTol', ode_reltol, 'AbsTol', ode_abstol);
 
            [t_out, X_out] = ode113(@(t,X) attitude_ode(t,X,p), ...
                                    t_span, [phi0; dphi0_val], opts);
 
            phi_rad   = X_out(:,1);
            phi_deg_t = rad2deg(phi_rad);
 
            %% --- SAM: Statistical Averaging Model ---
            edges   = (-180 : d_theta_deg : 180);
            centres = edges(1:end-1) + d_theta_deg/2;
            n_bins  = length(centres);
 
            dt_k        = diff(t_out);
            dt_k(end+1) = dt_k(end);   % pad to match length
 
            Delta_t = zeros(1, n_bins);
            for k = 1:length(t_out)
                bi = find(edges <= wrapTo180(phi_deg_t(k)), 1, 'last');
                if ~isempty(bi) && bi <= n_bins
                    Delta_t(bi) = Delta_t(bi) + dt_k(k);
                end
            end
 
            T_total = sum(dt_k);
            pdf_SAM = Delta_t / (T_total * d_theta_deg);
 
            %% --- Cd_bar from SAM ---
            Cd_bins = zeros(1, n_bins);
            for bi = 1:n_bins
                Cd_bins(bi) = interp1(sails_c(1).alpha, sails_c(1).Cd, ...
                                      centres(bi), 'pchip', 'extrap');
            end
            Cd_bar = sum(Cd_bins .* pdf_SAM) * d_theta_deg;
            Cd_bar_3D(alt_i, ac, si) = Cd_bar;
 
            fprintf('    phi=%2ddeg | accom=%s | Cd_bar=%.4f | Iyy=%.5f kg.m2\n', ...
                    phi_sel, accom_str, Cd_bar, Iyy_val);
 
            %% --- Store full reference results ---
            ref_results(alt_i, ac, si).label        = sprintf('%dkm_%ddeg_%s', ...
                                                       alt, phi_sel, accom_str);
            ref_results(alt_i, ac, si).alt_km       = alt;
            ref_results(alt_i, ac, si).phi_deg_sail = phi_sel;
            ref_results(alt_i, ac, si).accom        = accom_vals(ac);
            ref_results(alt_i, ac, si).Iyy          = Iyy_val;
            ref_results(alt_i, ac, si).t            = t_out;
            ref_results(alt_i, ac, si).phi_deg      = phi_deg_t;
            ref_results(alt_i, ac, si).phi_wrapped  = wrapTo180(phi_deg_t);
            ref_results(alt_i, ac, si).omega_deg_s  = rad2deg(X_out(:,2));
            ref_results(alt_i, ac, si).centres      = centres;
            ref_results(alt_i, ac, si).pdf_SAM      = pdf_SAM;
            ref_results(alt_i, ac, si).Cd_bar       = Cd_bar;
            ref_results(alt_i, ac, si).T_orb_min    = T_orb/60;
 
            %% --- Periodic save: every sail angle within alt+accom loop ---
            save('attitude_SAM_sweep_partial.mat', ...
                 'Cd_bar_3D', 'ref_results', ...
                 'alt_vec', 'accom_vec', 'accom_vals', 'phi_vec', ...
                 'inertiaDB', 'ref_geom');
 
        end % phi loop
    end % accom loop
end % alt loop
 
%% =========================================================================
%  SECTION 5 — FINAL SAVE
% =========================================================================
 
save('attitude_SAM_sweep_FINAL.mat', ...
     'Cd_bar_3D', 'ref_results', ...
     'alt_vec', 'accom_vec', 'accom_vals', 'phi_vec', ...
     'inertiaDB', 'ref_geom');
 
fprintf('\n=== Sweep complete ===\n');
fprintf('Saved: attitude_SAM_sweep_FINAL.mat\n');
fprintf('Cd_bar_3D dimensions: [alt x accom x phi] = %s\n', ...
        mat2str(size(Cd_bar_3D)));
fprintf('alt_vec    = %s\n', mat2str(alt_vec));
fprintf('accom_vals = %s\n', mat2str(accom_vals));
fprintf('phi_vec    = %s deg\n', mat2str(phi_vec));
 
%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================
 
function dX = attitude_ode(t, X, p)
    phi  = X(1);
    dphi = X(2);
    [~, tau_aero] = compute_rho_torque(t, phi, p);
    n_orb  = 2*pi / p.T_orb;
    tau_gg = -1.5 * n_orb^2 * p.Iyy * sin(2*phi);
    dX     = [dphi; (tau_aero + tau_gg) / p.Iyy];
end
 
function [rho, tau_total] = compute_rho_torque(t, phi, p)
 
    n_orb = 2*pi / p.T_orb;
    nu    = mod(n_orb*t, 2*pi);
 
    [r_eci, ~] = kep2eci_circ(p.a_orb, p.inc, p.RAAN, nu);
 
    utc_dt  = datetime(p.year, p.month, p.day, p.hour, 0, 0) + seconds(t);
    utc_vec = [utc_dt.Year, utc_dt.Month, utc_dt.Day, ...
               utc_dt.Hour, utc_dt.Minute, utc_dt.Second];
 
    lla   = eci2lla(r_eci', utc_vec);
    doy   = day(utc_dt, 'dayofyear');
    utc_s = hour(utc_dt)*3600 + minute(utc_dt)*60 + second(utc_dt);
 
    [~, rho_arr] = atmosnrlmsise00( ...
        lla(3), lla(1), lla(2), ...
        p.year, doy, utc_s, ...
        p.f107a, p.f107, p.ap);
 
    rho = rho_arr(6);
    q   = 0.5 * rho * p.v_circ^2;
 
    sails     = set_sails();
    phi_deg   = rad2deg(phi);
    tau_total = 0;
 
    for s = 1:length(sails)
        Cm_s      = interp1(sails(s).alpha, sails(s).Cm_pitch, ...
                            phi_deg, 'pchip', 'extrap');
        tau_total = tau_total + q * sails(s).A_ref * sails(s).L_ref * Cm_s;
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