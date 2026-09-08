%% =========================================================================
%  ROADM.m  --  Reduced-Order Attitude Dynamics Model
%
%  Stage 2 of the pipeline (ADM -> ROADM -> SAM).
%
%  Propagates the single-degree-of-freedom pitch equation of motion for
%  every (apex half-angle phi, accommodation coefficient alpha, altitude h)
%  combination. Aerodynamic torque at each timestep interpolates the
%  pitching-moment coefficient C_M(theta) from the ADM at the instantaneous
%  pitch angle and scales it by the orbit-varying dynamic pressure from the
%  NRLMSISE-00 density model; gravity-gradient torque is added.
%
%  Output:  roadm_results.mat
%      roadm      : struct [n_alt x n_accom x n_phi] with the pitch angle
%                   time history theta(t) and rate theta_dot(t) for each case
%      sail_db    : struct [n_accom x n_phi] with the ADM C_D(theta) /
%                   C_M(theta) curves and reference geometry / inertia
%      plus alt_vec, accom_vec, accom_vals, phi_vec, inertiaDB, ref_geom
%      and a meta struct recording the epoch / solver settings.
%
%  SAM.m consumes roadm_results.mat. Requires adbsat_processed/ (ADM output,
%  not in this repo) and inertia_tensors.mat + ref_geometry.mat (run
%  setup_databases).
%
%  Canonical settings from the former MASTERCODE.m: June epoch, 98 deg
%  inclination, ode45, zero initial pitch rate. NOTE: N_ORBITS below is 1
%  (as in MASTERCODE.m). The dissertation headline runs used 15 orbits --
%  set N_ORBITS = 15 to reproduce those.
% =========================================================================
clear; clc; close all;
clear functions

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'lib')));
addpath(genpath(fullfile(here, 'adbsat_processed')));
cd(here);

%% ----- SECTION 1 : constants & settings ---------------------------------
mu = 3.986004418e14;
Re = 6.3781e6;

year_ep  = 2025;   month_ep = 6;   day_ep = 21;   hour_utc = 10;
f107a    = 150;    f107     = 150; ap     = 15;

N_ORBITS   = 1;                 % <-- set to 15 for the dissertation headline runs
n_steps    = 3000 * N_ORBITS;   % keep ~3000 samples per orbit
ode_reltol = 1e-8;
ode_abstol = 1e-10;

phi0  = deg2rad(5);             % initial pitch angle
dphi0 = 0.0;                    % initial pitch rate

inc_deg  = 98.0;               % Sun-synchronous
RAAN_deg = 90.0;

%% ----- SECTION 2 : databases ------------------------------------------
if ~isfile('inertia_tensors.mat')
    error('inertia_tensors.mat not found -- run setup_databases inertia');
end
if ~isfile('ref_geometry.mat')
    error('ref_geometry.mat not found -- run setup_databases geometry');
end
load('inertia_tensors.mat', 'inertiaDB');
load('ref_geometry.mat',    'ref_geom');

adbsat_base = fullfile(here, 'adbsat_processed');

%% ----- SECTION 3 : parameter space ----------------------------------
alt_vec    = [350, 450, 650];
accom_vec  = {'0p50','0p55','0p60','0p65','0p70','0p75','0p80','0p85','0p90','0p95','1p00'};
accom_vals = 0.50:0.05:1.00;
phi_vec    = [45 50 55 60 65 70 75 80 85];

n_alts  = numel(alt_vec);
n_accom = numel(accom_vec);
n_phi   = numel(phi_vec);

roadm   = struct('alt_km',{},'accom',{},'phi_deg',{}, ...
                 't',{},'theta_deg',{},'theta_dot_deg',{},'T_orb_s',{},'Iyy',{});
roadm(n_alts,n_accom,n_phi).alt_km = [];   % preallocate
sail_db = struct('phi_deg',{},'accom',{},'alpha_deg',{},'Cd',{},'Cm_pitch',{}, ...
                 'A_ref',{},'L_ref',{},'Iyy',{});
sail_db(n_accom,n_phi).phi_deg = [];

%% ----- SECTION 4 : propagation sweep -------------------------------
opts = odeset('RelTol', ode_reltol, 'AbsTol', ode_abstol);

for alt_i = 1:n_alts
    alt     = alt_vec(alt_i);
    alt_dir = fullfile(adbsat_base, sprintf('%dkm', alt));
    fprintf('\n====== Altitude: %d km ======\n', alt);

    a_orb  = Re + alt*1e3;
    T_orb  = 2*pi * sqrt(a_orb^3 / mu);
    v_circ = sqrt(mu / a_orb);
    t_vec  = linspace(0, N_ORBITS*T_orb, n_steps)';
    doy0   = day(datetime(year_ep, month_ep, day_ep), 'dayofyear');
    jd0    = juliandate(datetime(year_ep, month_ep, day_ep, hour_utc, 0, 0));

    for ac = 1:n_accom
        accom_str = accom_vec{ac};
        fprintf('  Accommodation: %s\n', accom_str);

        for si = 1:n_phi
            phi = phi_vec(si);

            fname = fullfile(alt_dir, sprintf('%ddeg_CLL_accom_%s.mat', phi, accom_str));
            if ~isfile(fname)
                error('ADM file not found: %s', fname);
            end
            raw  = load(fname);
            aero = raw.aedb.aero;

            Cm_sort = aero.Cm_BY(:);  Cm_sort(isnan(Cm_sort)) = 0;
            Cd_sort = -aero.Cf_wX(:); Cd_sort(isnan(Cd_sort)) = 0;

            I_idx   = find(inertiaDB.angles_deg == phi);
            Iyy_val = inertiaDB.I(2, 2, I_idx);
            g_idx   = find(ref_geom.phi_deg == phi);

            clear sails_c
            sails_c(1).angle    = phi;
            sails_c(1).alpha    = (-179.5:1:179.5)';
            sails_c(1).Cm_pitch = Cm_sort;
            sails_c(1).Cd       = Cd_sort;
            sails_c(1).Iyy      = Iyy_val;
            sails_c(1).A_ref    = ref_geom.A_ref(g_idx);
            sails_c(1).L_ref    = ref_geom.L_ref(g_idx);

            p.a_orb  = a_orb;            p.inc    = deg2rad(inc_deg);
            p.RAAN   = deg2rad(RAAN_deg); p.T_orb = T_orb;
            p.v_circ = v_circ;           p.jd0    = jd0;
            p.doy0   = doy0;             p.year   = year_ep;
            p.hour   = hour_utc;         p.f107a  = f107a;
            p.f107   = f107;             p.ap     = ap;
            p.Iyy    = Iyy_val;          p.sails  = sails_c;
            p.n_sails = 1;              p.month  = month_ep;
            p.day    = day_ep;

            set_sails(sails_c);
            [t_out, X_out] = ode45(@(t,X) attitude_ode(t, X, p), ...
                                   t_vec, [phi0; dphi0], opts);

            roadm(alt_i,ac,si).alt_km        = alt;
            roadm(alt_i,ac,si).accom         = accom_vals(ac);
            roadm(alt_i,ac,si).phi_deg       = phi;
            roadm(alt_i,ac,si).t             = t_out;
            roadm(alt_i,ac,si).theta_deg     = rad2deg(X_out(:,1));
            roadm(alt_i,ac,si).theta_dot_deg = rad2deg(X_out(:,2));
            roadm(alt_i,ac,si).T_orb_s       = T_orb;
            roadm(alt_i,ac,si).Iyy           = Iyy_val;

            if alt_i == 1
                sail_db(ac,si).phi_deg   = phi;
                sail_db(ac,si).accom     = accom_vals(ac);
                sail_db(ac,si).alpha_deg = sails_c(1).alpha;
                sail_db(ac,si).Cd        = Cd_sort;
                sail_db(ac,si).Cm_pitch  = Cm_sort;
                sail_db(ac,si).A_ref     = sails_c(1).A_ref;
                sail_db(ac,si).L_ref     = sails_c(1).L_ref;
                sail_db(ac,si).Iyy       = Iyy_val;
            end

            fprintf('    phi=%2ddeg | accom=%s | done (%d steps)\n', ...
                    phi, accom_str, numel(t_out));

            if mod(si, 3) == 0
                save('roadm_results_partial.mat', 'roadm', 'sail_db', ...
                     'alt_vec', 'accom_vec', 'accom_vals', 'phi_vec', ...
                     'inertiaDB', 'ref_geom', '-v7.3');
            end
        end
    end
end

%% ----- SECTION 5 : save --------------------------------------------
meta = struct('N_ORBITS', N_ORBITS, 'n_steps', n_steps, ...
              'year', year_ep, 'month', month_ep, 'day', day_ep, 'hour', hour_utc, ...
              'f107a', f107a, 'f107', f107, 'ap', ap, ...
              'inc_deg', inc_deg, 'RAAN_deg', RAAN_deg, ...
              'phi0_deg', rad2deg(phi0), 'dphi0_dps', rad2deg(dphi0), ...
              'ode_reltol', ode_reltol, 'ode_abstol', ode_abstol, 'solver', 'ode45');

save('roadm_results.mat', 'roadm', 'sail_db', ...
     'alt_vec', 'accom_vec', 'accom_vals', 'phi_vec', ...
     'inertiaDB', 'ref_geom', 'meta', '-v7.3');

fprintf('\nSaved roadm_results.mat  (%d x %d x %d cases)\n', n_alts, n_accom, n_phi);
fprintf('Next: run SAM\n');

%% =========================================================================
%  LOCAL FUNCTIONS  (canonical versions, from MASTERCODE.m)
% =========================================================================
function dX = attitude_ode(t, X, p)
    phi  = X(1);
    dphi = X(2);
    [~, tau_aero] = compute_rho_torque(t, phi, p);
    n_orb  = 2*pi / p.T_orb;
    tau_gg = -1.5 * n_orb^2 * p.Iyy * sin(2*phi);
    ddphi  = (tau_aero + tau_gg) / p.Iyy;
    dX     = [dphi; ddphi];
end

function [rho, tau_total] = compute_rho_torque(t, phi, p)
    n_orb = 2*pi / p.T_orb;
    nu    = mod(n_orb * t, 2*pi);

    [r_eci, ~] = kep2eci_circ(p.a_orb, p.inc, p.RAAN, nu);

    utc_dt  = datetime(p.year, p.month, p.day, p.hour, 0, 0) + seconds(t);
    utc_vec = [utc_dt.Year, utc_dt.Month, utc_dt.Day, ...
               utc_dt.Hour, utc_dt.Minute, utc_dt.Second];
    lla   = eci2lla(r_eci', utc_vec);
    lat_d = lla(1);
    lon_d = lla(2);
    alt_m = lla(3);

    doy   = p.doy0 + t / 86400;
    utc_s = mod(p.hour * 3600 + t, 86400);

    [~, rho_arr] = atmosnrlmsise00(alt_m, lat_d, lon_d, ...
                                    p.year, doy, utc_s, ...
                                    p.f107a, p.f107, p.ap);
    rho = rho_arr(6);

    q = 0.5 * rho * p.v_circ^2;

    sails   = set_sails();
    n_sails = length(sails);

    phi_deg   = rad2deg(phi);
    tau_total = 0;
    for s = 1:n_sails
        Cm_s = interp1(sails(s).alpha, sails(s).Cm_pitch, ...
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
    R = [ cO, -sO*ci,  sO*si;
          sO,  cO*ci, -cO*si;
           0,  si,     ci   ];
    r = R * r_pf;
    v = R * v_pf;
end
