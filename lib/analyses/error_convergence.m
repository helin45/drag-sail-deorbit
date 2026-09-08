%% =========================================================================
%  ERROR ANALYSIS
%  1. n_steps convergence
%  2. Bin width (d_theta) sensitivity
%  3. Monte Carlo over dphi0
%  4. Atmospheric density scaling
% =========================================================================
clear; clc; close all; clear functions

addpath(genpath(fullfile(project_root(), 'adbsat_processed')));
addpath(genpath(fullfile(project_root(), 'adbsat_processed')));
addpath(genpath(fullfile(project_root(), 'adbsat_processed')));

%% Shared constants
mu=3.986004418e14; Re=6.3781e6;
year_ep=2025; month_ep=6; day_ep=21; hour_utc=10;
f107a=150; f107=150; ap=15;
ode_reltol=1e-8; ode_abstol=1e-10;
phi0=deg2rad(5); dphi0=0.0;
adbsat_base = fullfile(project_root(), 'adbsat_processed');

load('ref_geometry.mat','ref_geom');
load('inertia_tensors.mat','inertiaDB');

% Representative case: phi=75, alpha=0.75, 450 km
phi_sel=75; accom_str='0p75'; alt=450;
a_orb=Re+alt*1e3; T_orb=2*pi*sqrt(a_orb^3/mu); v_circ=sqrt(mu/a_orb);
doy0=day(datetime(year_ep,month_ep,day_ep),'dayofyear');
I_idx=find(inertiaDB.angles_deg==phi_sel); Iyy_val=inertiaDB.I(2,2,I_idx);
g_idx=find(ref_geom.phi_deg==phi_sel);
fname=fullfile(adbsat_base,sprintf('%dkm/%ddeg_CLL_accom_%s.mat',alt,phi_sel,accom_str));
raw=load(fname);
Cm_sort=raw.aedb.aero.Cm_BY(:); Cm_sort(isnan(Cm_sort))=0;
Cd_sort=-raw.aedb.aero.Cf_wX(:); Cd_sort(isnan(Cd_sort))=0;
sails_c(1).angle=phi_sel; sails_c(1).alpha=(-179.5:1:179.5)';
sails_c(1).Cm_pitch=Cm_sort; sails_c(1).Cd=Cd_sort;
sails_c(1).Iyy=Iyy_val; sails_c(1).A_ref=ref_geom.A_ref(g_idx); sails_c(1).L_ref=ref_geom.L_ref(g_idx);
p.a_orb=a_orb; p.inc=deg2rad(98); p.RAAN=deg2rad(90); p.T_orb=T_orb;
p.v_circ=v_circ; p.doy0=doy0; p.year=year_ep; p.hour=hour_utc;
p.f107a=f107a; p.f107=f107; p.ap=ap; p.Iyy=Iyy_val;
p.sails=sails_c; p.n_sails=1; p.month=month_ep; p.day=day_ep;
set_sails(sails_c);
opts=odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);

%% =========================================================================
%  1. n_steps CONVERGENCE
% =========================================================================
fprintf('\n--- 1. n_steps convergence ---\n');
n_steps_vec=[500 1000 2000 3000 5000 10000];
CD_nsteps=zeros(size(n_steps_vec));
for ni=1:length(n_steps_vec)
    t_vec=linspace(0,T_orb,n_steps_vec(ni))';
    [t_out,X_out]=ode45(@(t,X) attitude_ode(t,X,p),t_vec,[phi0;dphi0],opts);
    CD_nsteps(ni)=compute_CDbar(t_out,rad2deg(X_out(:,1)),sails_c,2.0);
    fprintf('  n_steps=%5d | Cd_bar=%.6f\n',n_steps_vec(ni),CD_nsteps(ni));
end
rel_err_nsteps=abs(CD_nsteps-CD_nsteps(end))./CD_nsteps(end)*100;
fprintf('Relative error vs n_steps=10000:\n');
for ni=1:length(n_steps_vec)
    fprintf('  n_steps=%5d | rel_err=%.4f%%\n',n_steps_vec(ni),rel_err_nsteps(ni));
end

%% =========================================================================
%  2. BIN WIDTH SENSITIVITY
% =========================================================================
fprintf('\n--- 2. Bin width sensitivity ---\n');
d_theta_vec=[0.5 1.0 2.0 4.0 8.0];
t_vec=linspace(0,T_orb,3000)';
[t_out,X_out]=ode45(@(t,X) attitude_ode(t,X,p),t_vec,[phi0;dphi0],opts);
phi_deg_t=rad2deg(X_out(:,1));
CD_bins=zeros(size(d_theta_vec));
for di=1:length(d_theta_vec)
    CD_bins(di)=compute_CDbar(t_out,phi_deg_t,sails_c,d_theta_vec(di));
    fprintf('  d_theta=%.1f deg | Cd_bar=%.6f\n',d_theta_vec(di),CD_bins(di));
end
rel_err_bins=abs(CD_bins-CD_bins(1))./CD_bins(1)*100;
fprintf('Relative error vs d_theta=0.5 deg:\n');
for di=1:length(d_theta_vec)
    fprintf('  d_theta=%.1f | rel_err=%.4f%%\n',d_theta_vec(di),rel_err_bins(di));
end

%% =========================================================================
%  3. MONTE CARLO over dphi0
% =========================================================================
fprintf('\n--- 3. Monte Carlo over initial pitch rate ---\n');
N_mc=100;
dphi0_max=deg2rad(2);  % adjust to physically justified max rad/s
rng(42);
dphi0_samples=dphi0_max*(2*rand(N_mc,1)-1);  % uniform [-max, +max]
t_vec=linspace(0,T_orb,3000)';
CD_mc=zeros(N_mc,1);
for mi=1:N_mc
    [t_out,X_out]=ode45(@(t,X) attitude_ode(t,X,p),t_vec,[phi0;dphi0_samples(mi)],opts);
    CD_mc(mi)=compute_CDbar(t_out,rad2deg(X_out(:,1)),sails_c,2.0);
    if mod(mi,10)==0; fprintf('  MC sample %d/%d\n',mi,N_mc); end
end
fprintf('MC results (dphi0 uniform +/-%.4f rad/s):\n',dphi0_max);
fprintf('  Mean Cd_bar = %.6f\n',mean(CD_mc));
fprintf('  Std         = %.6f\n',std(CD_mc));
fprintf('  CV          = %.4f%%\n',std(CD_mc)/mean(CD_mc)*100);
fprintf('  Min/Max     = %.6f / %.6f\n',min(CD_mc),max(CD_mc));

%% =========================================================================
%  4. ATMOSPHERIC DENSITY SCALING
% =========================================================================
fprintf('\n--- 4. Atmospheric density scaling ---\n');
rho_scales=[0.85 1.00 1.15];
t_vec=linspace(0,T_orb,3000)';
CD_rho=zeros(size(rho_scales));
for ri=1:length(rho_scales)
    p.rho_scale=rho_scales(ri);
    set_sails(sails_c);
    [t_out,X_out]=ode45(@(t,X) attitude_ode_scaled(t,X,p),t_vec,[phi0;dphi0],opts);
    CD_rho(ri)=compute_CDbar(t_out,rad2deg(X_out(:,1)),sails_c,2.0);
    fprintf('  rho_scale=%.2f | Cd_bar=%.6f\n',rho_scales(ri),CD_rho(ri));
end
fprintf('Sensitivity to +/-15%% density: delta_Cd_bar = %.6f (%.4f%%)\n',...
    abs(CD_rho(3)-CD_rho(1)), abs(CD_rho(3)-CD_rho(1))/CD_rho(2)*100);

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================
function CD_bar=compute_CDbar(t_out,phi_deg_t,sails_c,d_theta_deg)
    edges=(-180:d_theta_deg:180);
    centres=edges(1:end-1)+d_theta_deg/2;
    n_bins=length(centres);
    dt_k=diff(t_out); dt_k(end+1)=dt_k(end);
    Delta_t=zeros(1,n_bins);
    for k=1:length(t_out)
        bi=find(edges<=phi_deg_t(k),1,'last');
        if ~isempty(bi)&&bi<=n_bins
            Delta_t(bi)=Delta_t(bi)+dt_k(k);
        end
    end
    T_total=sum(dt_k);
    pdf_SAM=Delta_t/(T_total*d_theta_deg);
    Cd_bins=zeros(1,n_bins);
    for bi=1:n_bins
        Cd_bins(bi)=interp1(sails_c(1).alpha,sails_c(1).Cd,centres(bi),'pchip','extrap');
    end
    CD_bar=sum(Cd_bins.*pdf_SAM)*d_theta_deg;
end

function dX=attitude_ode(t,X,p)
    phi=X(1); dphi=X(2);
    [~,tau_aero]=compute_rho_torque(t,phi,p,1.0);
    n_orb=2*pi/p.T_orb;
    tau_gg=-1.5*n_orb^2*p.Iyy*sin(2*phi);
    dX=[dphi;(tau_aero+tau_gg)/p.Iyy];
end

function dX=attitude_ode_scaled(t,X,p)
    phi=X(1); dphi=X(2);
    [~,tau_aero]=compute_rho_torque(t,phi,p,p.rho_scale);
    n_orb=2*pi/p.T_orb;
    tau_gg=-1.5*n_orb^2*p.Iyy*sin(2*phi);
    dX=[dphi;(tau_aero+tau_gg)/p.Iyy];
end

function [rho,tau_total]=compute_rho_torque(t,phi,p,rho_scale)
    n_orb=2*pi/p.T_orb; nu=mod(n_orb*t,2*pi);
    [r_eci,~]=kep2eci_circ(p.a_orb,p.inc,p.RAAN,nu);
    utc_dt=datetime(p.year,p.month,p.day,p.hour,0,0)+seconds(t);
    utc_vec=[utc_dt.Year,utc_dt.Month,utc_dt.Day,utc_dt.Hour,utc_dt.Minute,utc_dt.Second];
    lla=eci2lla(r_eci',utc_vec);
    doy=p.doy0+t/86400; utc_s=mod(p.hour*3600+t,86400);
    [~,rho_arr]=atmosnrlmsise00(lla(3),lla(1),lla(2),p.year,doy,utc_s,p.f107a,p.f107,p.ap);
    rho=rho_arr(6)*rho_scale;
    q=0.5*rho*p.v_circ^2;
    sails=set_sails(); tau_total=0;
    phi_deg=rad2deg(phi);
    for s=1:length(sails)
        Cm_s=interp1(sails(s).alpha,sails(s).Cm_pitch,phi_deg,'pchip','extrap');
        tau_total=tau_total+q*sails(s).A_ref*sails(s).L_ref*Cm_s;
    end
end

function sails=set_sails(new_sails)
    persistent stored_sails;
    if nargin>0; stored_sails=new_sails; end
    sails=stored_sails;
end

function [r,v]=kep2eci_circ(a,inc,RAAN,nu)
    mu_=3.986004418e14;
    r_pf=a*[cos(nu);sin(nu);0]; v_pf=sqrt(mu_/a)*[-sin(nu);cos(nu);0];
    cO=cos(RAAN);sO=sin(RAAN);ci=cos(inc);si=sin(inc);
    R=[cO,-sO*ci,sO*si;sO,cO*ci,-cO*si;0,si,ci];
    r=R*r_pf; v=R*v_pf;
end