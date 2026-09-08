clear; clc; close all; clear functions

load('drag_sail_attitude_results.mat')
load('inertia_tensors.mat','inertiaDB');

%% Constants
mu = 3.986004418e14;
Re = 6.3781e6;
omega_E = 7.2921159e-5;   % Earth rotation [rad/s]

year_ep  = 2025;
month_ep = 1;
day_ep   = 15;

f107a = 150; f107 = 150; ap = 15;

ode_reltol = 1e-6;
ode_abstol = 1e-8;

%% Simulation duration
t_final = 2 * 24 * 3600;
tspan   = [0 t_final];

%% Configuration
phi_sel = 75;
alt     = 450;
accom   = '0p75';

adbsat_base = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';

%% Roberts sphere inertia
m_host    = 10;
rho_body  = 760;
r_sphere  = (3*m_host / (4*pi*rho_body))^(1/3);
Iyy_val   = 0.4 * m_host * r_sphere^2;

g_idx = find(ref_geom.phi_deg == phi_sel);

fname = fullfile(adbsat_base, sprintf('%dkm',alt), ...
                 sprintf('%ddeg_CLL_accom_%s.mat',phi_sel,accom));

raw = load(fname);

Cm_sort = raw.aedb.aero.Cm_BY(:); Cm_sort(isnan(Cm_sort))=0;
Cd_sort = -raw.aedb.aero.Cf_wX(:); Cd_sort(isnan(Cd_sort))=0;

%% Sail
clear sails_c
sails_c(1).angle    = phi_sel;
sails_c(1).alpha    = (-179.5:1:179.5)';
sails_c(1).Cm_pitch = Cm_sort;
sails_c(1).Cd       = Cd_sort;
sails_c(1).Iyy      = Iyy_val;
sails_c(1).A_ref    = ref_geom.A_ref(g_idx);
sails_c(1).L_ref    = ref_geom.L_ref(g_idx);

%% Orbit
a_orb  = Re + alt*1e3;
T_orb  = 2*pi*sqrt(a_orb^3/mu);

%% Initial conditions
phi0_test  = deg2rad(160);
dphi0_test = 2*pi / T_orb;

%% Cases
cases = [0, 12];
results = struct;

for k = 1:length(cases)

    hour_utc = cases(k);

    p.a_orb  = a_orb;
    p.inc    = deg2rad(0);
    p.RAAN   = deg2rad(90);
    p.T_orb  = T_orb;

    p.year  = year_ep;
    p.month = month_ep;
    p.day   = day_ep;
    p.hour  = hour_utc;

    p.f107a = f107a;
    p.f107  = f107;
    p.ap    = ap;

    p.Iyy   = Iyy_val;
    p.sails = sails_c;

    set_sails(sails_c);

    opts = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);

    [t_out, X_out] = ode113(@(t,X) attitude_ode(t,X,p,omega_E), ...
                            tspan, [phi0_test; dphi0_test], opts);

    results(k).hour          = hour_utc;
    results(k).t             = t_out;
    results(k).theta         = rad2deg(X_out(:,1));
    results(k).omega         = rad2deg(X_out(:,2));
    results(k).theta_wrapped = wrapTo180(rad2deg(X_out(:,1)));

end

%% ===== FUNCTIONS =====

function dX = attitude_ode(t,X,p,omega_E)

    phi  = X(1);
    dphi = X(2);

    [~,tau_aero] = compute_rho_torque(t,phi,p,omega_E);

    n_orb  = 2*pi/p.T_orb;
    tau_gg = -1.5*n_orb^2*p.Iyy*sin(2*phi);

    dX = [dphi; (tau_aero + tau_gg)/p.Iyy];

end

function [rho,tau_total] = compute_rho_torque(t,phi,p,omega_E)

    n_orb = 2*pi/p.T_orb;
    nu    = mod(n_orb*t,2*pi);

    [r_eci, v_eci] = kep2eci_circ(p.a_orb,p.inc,p.RAAN,nu);

    %% Time
    utc_dt = datetime(p.year,p.month,p.day,p.hour,0,0) + seconds(t);
    utc_vec = [utc_dt.Year, utc_dt.Month, utc_dt.Day,...
               utc_dt.Hour, utc_dt.Minute, utc_dt.Second];

    lla = eci2lla(r_eci',utc_vec);

    doy   = day(utc_dt,'dayofyear');
    utc_s = hour(utc_dt)*3600 + minute(utc_dt)*60 + second(utc_dt);

    [~,rho_arr] = atmosnrlmsise00(...
        lla(3), lla(1), lla(2),...
        p.year, doy, utc_s,...
        p.f107a, p.f107, p.ap);

    rho = rho_arr(6);

    %% FIX: relative velocity
    omega_vec = [0;0;omega_E];
    v_atm = cross(omega_vec, r_eci);
    v_rel = v_eci - v_atm;

    q = 0.5 * rho * norm(v_rel)^2;

    sails = set_sails();
    phi_deg = rad2deg(phi);

    tau_total = 0;

    for s = 1:length(sails)
        Cm_s = interp1(sails(s).alpha, sails(s).Cm_pitch,...
                       phi_deg,'pchip','extrap');

        tau_total = tau_total + q * sails(s).A_ref *...
                                   sails(s).L_ref * Cm_s;
    end
end

function sails = set_sails(new_sails)
    persistent stored_sails;
    if nargin>0
        stored_sails = new_sails;
    end
    sails = stored_sails;
end

function [r,v] = kep2eci_circ(a,inc,RAAN,nu)

    mu_ = 3.986004418e14;

    r_pf = a*[cos(nu);sin(nu);0];
    v_pf = sqrt(mu_/a)*[-sin(nu);cos(nu);0];

    cO=cos(RAAN); sO=sin(RAAN);
    ci=cos(inc);  si=sin(inc);

    R=[cO,-sO*ci,sO*si;
       sO, cO*ci,-cO*si;
       0,  si,    ci];

    r=R*r_pf;
    v=R*v_pf;
end