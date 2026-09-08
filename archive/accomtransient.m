clear functions
load('drag_sail_attitude_results.mat')

mu=3.986004418e14; Re=6.3781e6;
year_ep=2025; month_ep=6; day_ep=21; hour_utc=10;
f107a=150; f107=150; ap=15;
n_steps=3000; ode_reltol=1e-8; ode_abstol=1e-10;
phi0=deg2rad(5); dphi0=0.0;
load('inertia_tensors.mat','inertiaDB');
adbsat_base = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';

% Fixed phi=75deg, all altitudes, selected accommodation values
phi_sel    = 75;
si_sel     = find(phi_vec == phi_sel);
accom_sel  = {'0p50','0p75','1p00'};
accom_vals_sel = [0.50 0.75 1.00];
clrs = [0.2 0.2 0.8; 0.2 0.7 0.4; 0.8 0.2 0.2];
alt_labels = {'350 km','450 km','650 km'};

I_idx   = find(inertiaDB.angles_deg == phi_sel);
Iyy_val = inertiaDB.I(2,2,I_idx);
g_idx   = find(ref_geom.phi_deg == phi_sel);

fig = figure('Position',[60 60 1200 800]);
for alt_i = 1:3
    alt    = alt_vec(alt_i);
    alt_dir = fullfile(adbsat_base, sprintf('%dkm', alt));
    a_orb  = Re + alt*1e3;
    T_orb  = 2*pi*sqrt(a_orb^3/mu);
    v_circ = sqrt(mu/a_orb);
    t_vec  = linspace(0,T_orb,n_steps)';
    doy0   = day(datetime(year_ep,month_ep,day_ep),'dayofyear');

    subplot(3,1,alt_i);
    hold on;

    for ai = 1:length(accom_sel)
        fname = fullfile(alt_dir, sprintf('%ddeg_CLL_accom_%s.mat', phi_sel, accom_sel{ai}));
        raw   = load(fname);
        Cm_sort = raw.aedb.aero.Cm_BY(:); Cm_sort(isnan(Cm_sort))=0;
        Cd_sort = -raw.aedb.aero.Cf_wX(:); Cd_sort(isnan(Cd_sort))=0;

        clear sails_c
        sails_c(1).angle    = phi_sel;
        sails_c(1).alpha    = (-179.5:1:179.5)';
        sails_c(1).Cm_pitch = Cm_sort;
        sails_c(1).Cd       = Cd_sort;
        sails_c(1).Iyy      = Iyy_val;
        sails_c(1).A_ref    = ref_geom.A_ref(g_idx);
        sails_c(1).L_ref    = ref_geom.L_ref(g_idx);

        p.a_orb=a_orb; p.inc=deg2rad(98); p.RAAN=deg2rad(90); p.T_orb=T_orb;
        p.v_circ=v_circ; p.doy0=doy0; p.year=year_ep;
        p.hour=hour_utc; p.f107a=f107a; p.f107=f107; p.ap=ap;
        p.Iyy=Iyy_val; p.sails=sails_c; p.n_sails=1;
        p.month=month_ep; p.day=day_ep;

        set_sails(sails_c);
        opts = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
        [t_out,X_out] = ode45(@(t,X) attitude_ode(t,X,p), t_vec, [phi0;dphi0], opts);

        plot(t_out/60, rad2deg(X_out(:,1)), 'Color',clrs(ai,:), 'LineWidth',1.2, ...
             'DisplayName', sprintf('$\\alpha = %.2f$', accom_vals_sel(ai)));
    end

    yline(0,'k--','HandleVisibility','off');
    ylabel('$\theta$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i == 3
        xlabel('Time (min)','Interpreter','latex');
        lg = legend('Location','best'); set(lg,'Interpreter','latex');
    end
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex',...
        'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
end
sgtitle(sprintf('Effect of $\\alpha$ on pitch transient ($\\phi = %d^\\circ$)', phi_sel),...
        'Interpreter','latex','FontSize',12);

savefig(fig,'plot_accom_effect_transient.fig');
exportgraphics(fig,'plot_accom_effect_transient.png','Resolution',300);
fprintf('Saved\n');

%% Local functions
function dX = attitude_ode(t,X,p)
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
    sails=set_sails(); phi_deg=rad2deg(phi); tau_total=0;
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