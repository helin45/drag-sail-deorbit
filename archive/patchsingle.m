% Load
load('drag_sail_attitude_results.mat')
clear functions

% Constants
mu=3.986004418e14; Re=6.3781e6;
year_ep=2025; month_ep=6; day_ep=21; hour_utc=10;
f107a=150; f107=150; ap=15;
n_steps=3000; ode_reltol=1e-8; ode_abstol=1e-10;
phi0=deg2rad(5); dphi0=0.0; d_theta_deg=2.0;
load('inertia_tensors.mat','inertiaDB');

% Target combination
alt_i=2; ac=1; si=4;   % 450km, 0p50, 60deg
alt=450; phi=60; accom_str='0p50';

fname = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed/450km/60deg_CLL_accom_0p50.mat';
raw=load(fname); aero=raw.aedb.aero;
Cm_sort=aero.Cm_BY(:); Cm_sort(isnan(Cm_sort))=0;
Cd_sort=-aero.Cf_wX(:); Cd_sort(isnan(Cd_sort))=0;

I_idx=find(inertiaDB.angles_deg==phi); Iyy_val=inertiaDB.I(2,2,I_idx);
g_idx=find(ref_geom.phi_deg==phi);

clear sails_c
sails_c(1).angle=phi; sails_c(1).alpha=(-179.5:1:179.5)';
sails_c(1).Cm_pitch=Cm_sort; sails_c(1).Cd=Cd_sort;
sails_c(1).Iyy=Iyy_val; sails_c(1).A_ref=ref_geom.A_ref(g_idx);
sails_c(1).L_ref=ref_geom.L_ref(g_idx);

a_orb=Re+alt*1e3; T_orb=2*pi*sqrt(a_orb^3/mu); v_circ=sqrt(mu/a_orb);
t_vec=linspace(0,T_orb,n_steps)';
jd0=juliandate(datetime(year_ep,month_ep,day_ep,hour_utc,0,0));
doy0=day(datetime(year_ep,month_ep,day_ep),'dayofyear');

p.a_orb=a_orb; p.inc=deg2rad(98); p.RAAN=deg2rad(90); p.T_orb=T_orb;
p.v_circ=v_circ; p.jd0=jd0; p.doy0=doy0; p.year=year_ep;
p.hour=hour_utc; p.f107a=f107a; p.f107=f107; p.ap=ap;
p.Iyy=Iyy_val; p.sails=sails_c; p.n_sails=1;
p.month=month_ep; p.day=day_ep;

set_sails(sails_c);
opts=odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
[t_out,X_out]=ode45(@(t,X)attitude_ode(t,X,p),t_vec,[phi0;dphi0],opts);
phi_deg_t=rad2deg(X_out(:,1));

edges=(-180:d_theta_deg:180); centres=edges(1:end-1)+d_theta_deg/2;
n_bins=length(centres); dt_k=diff(t_out); dt_k(end+1)=dt_k(end);
Delta_t=zeros(1,n_bins);
for k=1:length(t_out)
    bi=find(edges<=phi_deg_t(k),1,'last');
    if ~isempty(bi)&&bi<=n_bins; Delta_t(bi)=Delta_t(bi)+dt_k(k); end
end
pdf_SAM=Delta_t/(sum(dt_k)*d_theta_deg);
Cd_bins=zeros(1,n_bins);
for bi=1:n_bins
    Cd_bins(bi)=interp1(sails_c(1).alpha,sails_c(1).Cd,centres(bi),'pchip','extrap');
end
Cd_bar=sum(Cd_bins.*pdf_SAM)*d_theta_deg;
Cd_bar_3D(alt_i,ac,si)=Cd_bar;
fprintf('Patched: 450km | 0p50 | 60deg | Cd_bar=%.4f\n',Cd_bar);

save('drag_sail_attitude_results.mat','Cd_bar_3D','ref_results',...
     'alt_vec','accom_vec','accom_vals','phi_vec','inertiaDB','ref_geom');
fprintf('Saved\n');

function dX=attitude_ode(t,X,p)
    phi=X(1); dphi=X(2);
    [~,tau_aero]=compute_rho_torque(t,phi,p);
    n_orb=2*pi/p.T_orb;
    tau_gg=-1.5*n_orb^2*p.Iyy*sin(2*phi);
    dX=[dphi;(tau_aero+tau_gg)/p.Iyy];
end
function [rho,tau_total]=compute_rho_torque(t,phi,p)
    n_orb=2*pi/p.T_orb; nu=mod(n_orb*t,2*pi);
    [r_eci,~]=kep2eci_circ(p.a_orb,p.inc,p.RAAN,nu);
    utc_dt=datetime(p.year,p.month,p.day,p.hour,0,0)+seconds(t);
    utc_vec=[utc_dt.Year,utc_dt.Month,utc_dt.Day,utc_dt.Hour,utc_dt.Minute,utc_dt.Second];
    lla=eci2lla(r_eci',utc_vec);
    doy=p.doy0+t/86400; utc_s=mod(p.hour*3600+t,86400);
    [~,rho_arr]=atmosnrlmsise00(lla(3),lla(1),lla(2),p.year,doy,utc_s,p.f107a,p.f107,p.ap);
    rho=rho_arr(6); q=0.5*rho*p.v_circ^2;
    sails=set_sails(); n_sails=length(sails);
    phi_deg=rad2deg(phi); tau_total=0;
    for s=1:n_sails
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