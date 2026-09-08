P_ROOT = fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(P_ROOT,'lib'))); cd(P_ROOT);  % repo root + lib on path
load('sam_results.mat')
clear functions

mu=3.986004418e14; Re=6.3781e6;
year_ep=2025; month_ep=6; day_ep=21; hour_utc=10;
f107a=150; f107=150; ap=15;
n_steps=3000; ode_reltol=1e-8; ode_abstol=1e-10;
phi0=deg2rad(5); dphi0=0.0; d_theta_deg=2.0;
load('inertia_tensors.mat','inertiaDB');

adbsat_base = fullfile(project_root(), 'adbsat_processed');

% Selected accommodation coefficients to compare
accom_sel      = {'0p50','0p75','1p00'};
accom_vals_sel = [0.50, 0.75, 1.00];
clrs_accom_sel = [0.2 0.2 0.8; 0.2 0.7 0.4; 0.8 0.2 0.2];

% Fixed phi=75deg, all altitudes
phi     = 75;
si_sel  = find(phi_vec == phi);
I_idx   = find(inertiaDB.angles_deg == phi);
Iyy_val = inertiaDB.I(2,2,I_idx);
g_idx   = find(ref_geom.phi_deg == phi);

pdf_store   = cell(3, length(accom_sel));
cent_store  = cell(3, length(accom_sel));
alt_labels  = {'350 km', '450 km', '650 km'};

for alt_i = 1:3
    alt     = alt_vec(alt_i);
    alt_dir = fullfile(adbsat_base, sprintf('%dkm', alt));
    a_orb   = Re + alt*1e3;
    T_orb   = 2*pi*sqrt(a_orb^3/mu);
    v_circ  = sqrt(mu/a_orb);
    t_vec   = linspace(0,T_orb,n_steps)';
    jd0     = juliandate(datetime(year_ep,month_ep,day_ep,hour_utc,0,0));
    doy0    = day(datetime(year_ep,month_ep,day_ep),'dayofyear');

    for ai = 1:length(accom_sel)
        fname = fullfile(alt_dir, sprintf('%ddeg_CLL_accom_%s.mat', phi, accom_sel{ai}));
        raw   = load(fname);
        Cm_sort = raw.aedb.aero.Cm_BY(:); Cm_sort(isnan(Cm_sort))=0;
        Cd_sort = -raw.aedb.aero.Cf_wX(:); Cd_sort(isnan(Cd_sort))=0;

        clear sails_c
        sails_c(1).angle    = phi;
        sails_c(1).alpha    = (-179.5:1:179.5)';
        sails_c(1).Cm_pitch = Cm_sort;
        sails_c(1).Cd       = Cd_sort;
        sails_c(1).Iyy      = Iyy_val;
        sails_c(1).A_ref    = ref_geom.A_ref(g_idx);
        sails_c(1).L_ref    = ref_geom.L_ref(g_idx);

        p.a_orb=a_orb; p.inc=deg2rad(98); p.RAAN=deg2rad(90); p.T_orb=T_orb;
        p.v_circ=v_circ; p.jd0=jd0; p.doy0=doy0; p.year=year_ep;
        p.hour=hour_utc; p.f107a=f107a; p.f107=f107; p.ap=ap;
        p.Iyy=Iyy_val; p.sails=sails_c; p.n_sails=1;
        p.month=month_ep; p.day=day_ep;

        set_sails(sails_c);
        opts = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
        [t_out,X_out] = ode45(@(t,X) attitude_ode(t,X,p), t_vec, [phi0;dphi0], opts);

        phi_deg_t = rad2deg(X_out(:,1));
        edges   = (-180:d_theta_deg:180);
        centres = edges(1:end-1)+d_theta_deg/2;
        n_bins  = length(centres);
        dt_k    = diff(t_out); dt_k(end+1)=dt_k(end);

        Delta_t = zeros(1,n_bins);
        for k = 1:length(t_out)
            bi = find(edges<=phi_deg_t(k),1,'last');
            if ~isempty(bi)&&bi<=n_bins
                Delta_t(bi)=Delta_t(bi)+dt_k(k);
            end
        end
        pdf_SAM = Delta_t/(sum(dt_k)*d_theta_deg);
        pdf_store{alt_i,ai}  = pdf_SAM;
        cent_store{alt_i,ai} = centres;
        fprintf('Done: %dkm | phi=%ddeg | accom=%s\n', alt, phi, accom_sel{ai});
    end
end

%% Plot 1: PDF bar subplots — all altitudes x selected accom (phi=75deg)
fig1 = figure('Position',[60 60 1200 800]);
for alt_i = 1:3
    for ai = 1:length(accom_sel)
        subplot(3,3,(alt_i-1)*3+ai);
        bar(cent_store{alt_i,ai}, pdf_store{alt_i,ai}, 1, ...
            'FaceColor',clrs_accom_sel(ai,:),'EdgeColor','none','FaceAlpha',0.8);
        xlim([-20 20]);
        ttl = sprintf('%s, $\\alpha=%.2f$', alt_labels{alt_i}, accom_vals_sel(ai));
        title(ttl,'Interpreter','latex','FontSize',9);
        if ai==1; ylabel('$p(\theta)$ (deg$^{-1}$)','Interpreter','latex'); end
        if alt_i==3; xlabel('$\theta$ (deg)','Interpreter','latex'); end
        set(gca,'Box','on','LineWidth',1.0,'FontSize',9,'TickLabelInterpreter','latex');
    end
end
savefig(fig1,'plot_PDF_accom_effect_bar.fig');
exportgraphics(fig1,'plot_PDF_accom_effect_bar.png','Resolution',300);
fprintf('Saved PDF accom bar plot\n');

%% Plot 2: Scatter PDF zoomed -10 to 10, log scale — all accom overlaid per altitude
fig2 = figure('Position',[60 500 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for ai = 1:length(accom_sel)
        c    = cent_store{alt_i,ai};
        p_   = pdf_store{alt_i,ai};
        mask = c>=-10 & c<=10 & p_>0;
        scatter(c(mask), p_(mask), 30, 'filled', ...
                'MarkerFaceColor',clrs_accom_sel(ai,:), ...
                'MarkerFaceAlpha',0.8, ...
                'DisplayName',sprintf('$\\alpha=%.2f$',accom_vals_sel(ai)));
        hold on;
    end
    set(gca,'YScale','log');
    xlim([-10 10]);
    xlabel('$\theta$ (deg)','Interpreter','latex');
    ylabel('$p(\theta)$ (deg$^{-1}$)','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i==3
        lg = legend('Location','best','FontSize',9);
        set(lg,'Interpreter','latex');
    end
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex', ...
        'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
end
savefig(fig2,'plot_PDF_accom_effect_log.fig');
exportgraphics(fig2,'plot_PDF_accom_effect_log.png','Resolution',300);
fprintf('Saved PDF accom log plot\n');

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