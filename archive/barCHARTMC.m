%% mc_barchart.m
cd('/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE')
load('mc_workspace.mat');

t_vec = linspace(0,T_orb,3000)';
opts  = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
N_mc  = 100;
rng(42);

phi_vec    = [45 50 55 60 65 70 75 80 85];
accom_str  = '0p75';
alt        = 450;
adbsat_base= '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';

% Sample dphi0 and rho together
dphi0_samples = deg2rad(-2 + 4*rand(N_mc,1));
rho_samples   = 0.85 + 0.30*rand(N_mc,1);
phi0_fixed    = deg2rad(5);

CD_nominal = zeros(length(phi_vec),1);
CD_min     = zeros(length(phi_vec),1);
CD_max     = zeros(length(phi_vec),1);
CD_mean    = zeros(length(phi_vec),1);
CD_std     = zeros(length(phi_vec),1);

for si = 1:length(phi_vec)
    phi = phi_vec(si);
    fprintf('Processing phi=%d deg...\n', phi);

    fname = fullfile(adbsat_base, sprintf('%dkm/%ddeg_CLL_accom_%s.mat', alt, phi, accom_str));
    raw   = load(fname);
    Cm_s  = raw.aedb.aero.Cm_BY(:); Cm_s(isnan(Cm_s))=0;
    Cd_s  = -raw.aedb.aero.Cf_wX(:); Cd_s(isnan(Cd_s))=0;
    sail_alpha = (-179.5:1:179.5)';

    I_idx    = find(inertiaDB.angles_deg==phi);
    Iyy_val  = inertiaDB.I(2,2,I_idx);
    g_idx    = find(ref_geom.phi_deg==phi);
    sail_Aref= ref_geom.A_ref(g_idx);
    sail_Lref= ref_geom.L_ref(g_idx);

    % Nominal: dphi0=0, rho_scale=1.0
    p_nom = build_params(p, 1.0, sail_alpha, Cm_s, Cd_s, Iyy_val, sail_Aref, sail_Lref);
    [t_out,X_out] = ode45(@(t,X) attitude_ode_direct(t,X,p_nom), t_vec, [phi0_fixed; 0], opts);
    CD_nominal(si) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_alpha, Cd_s, 2.0);

    % MC: vary dphi0 and rho simultaneously
    CD_mc_phi = zeros(N_mc,1);
    for mi = 1:N_mc
        p_loc = build_params(p, rho_samples(mi), sail_alpha, Cm_s, Cd_s, Iyy_val, sail_Aref, sail_Lref);
        [t_out,X_out] = ode45(@(t,X) attitude_ode_direct(t,X,p_loc), t_vec, ...
                              [phi0_fixed; dphi0_samples(mi)], opts);
        CD_mc_phi(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_alpha, Cd_s, 2.0);
    end

    CD_min(si)  = min(CD_mc_phi);
    CD_max(si)  = max(CD_mc_phi);
    CD_mean(si) = mean(CD_mc_phi);
    CD_std(si)  = std(CD_mc_phi);

    fprintf('  phi=%d | Nominal=%.4f | Min=%.4f | Max=%.4f | Mean=%.4f\n',...
        phi, CD_nominal(si), CD_min(si), CD_max(si), CD_mean(si));

    save('mc_barchart_partial.mat','CD_nominal','CD_min','CD_max','CD_mean','CD_std','phi_vec');
end

save('mc_barchart.mat','CD_nominal','CD_min','CD_max','CD_mean','CD_std','phi_vec');
fprintf('Saved mc_barchart.mat\n');

%% Plot
err_low = CD_nominal - CD_min;
err_hi  = CD_max     - CD_nominal;

fig = figure('Position',[100 100 800 500]);
b = bar(phi_vec, CD_nominal, 0.6, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none');
hold on;
errorbar(phi_vec, CD_nominal, err_low, err_hi, 'k.', 'LineWidth', 1.5, 'CapSize', 8);
xlabel('$\phi$ (deg)','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',12,'Box','on','LineWidth',1.2,...
    'XTick',phi_vec,'YGrid','on','GridAlpha',0.3);
exportgraphics(fig,'mc_barchart.png','Resolution',300);
fprintf('Saved mc_barchart.png\n');

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================
function p_out = build_params(p_in, rho_scale, alpha, Cm, Cd, Iyy, Aref, Lref)
    p_out            = p_in;
    p_out.rho_scale  = rho_scale;
    p_out.sail_alpha = alpha;
    p_out.sail_Cm    = Cm;
    p_out.sail_Cd    = Cd;
    p_out.sail_Iyy   = Iyy;
    p_out.sail_Aref  = Aref;
    p_out.sail_Lref  = Lref;
end

function CD_bar = compute_CDbar(t_out, phi_deg_t, sail_alpha, sail_Cd, d_theta_deg)
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
        Cd_bins(bi)=interp1(sail_alpha, sail_Cd, centres(bi),'pchip','extrap');
    end
    CD_bar=sum(Cd_bins.*pdf_SAM)*d_theta_deg;
end

function dX = attitude_ode_direct(t, X, p)
    phi=X(1); dphi=X(2);
    [~,tau_aero] = compute_rho_torque_direct(t, phi, p);
    n_orb  = 2*pi/p.T_orb;
    tau_gg = -1.5*n_orb^2*p.sail_Iyy*sin(2*phi);
    dX     = [dphi; (tau_aero+tau_gg)/p.sail_Iyy];
end

function [rho, tau_total] = compute_rho_torque_direct(t, phi, p)
    n_orb=2*pi/p.T_orb; nu=mod(n_orb*t,2*pi);
    [r_eci,~]=kep2eci_circ(p.a_orb,p.inc,p.RAAN,nu);
    utc_dt=datetime(p.year,p.month,p.day,p.hour,0,0)+seconds(t);
    utc_vec=[utc_dt.Year,utc_dt.Month,utc_dt.Day,utc_dt.Hour,utc_dt.Minute,utc_dt.Second];
    lla=eci2lla(r_eci',utc_vec);
    doy=p.doy0+t/86400; utc_s=mod(p.hour*3600+t,86400);
    [~,rho_arr]=atmosnrlmsise00(lla(3),lla(1),lla(2),p.year,doy,utc_s,p.f107a,p.f107,p.ap);
    rho=rho_arr(6)*p.rho_scale;
    q=0.5*rho*p.v_circ^2;
    phi_deg=rad2deg(phi);
    Cm_s=interp1(p.sail_alpha, p.sail_Cm, phi_deg,'pchip','extrap');
    tau_total=q*p.sail_Aref*p.sail_Lref*Cm_s;
end

function [r,v]=kep2eci_circ(a,inc,RAAN,nu)
    mu_=3.986004418e14;
    r_pf=a*[cos(nu);sin(nu);0]; v_pf=sqrt(mu_/a)*[-sin(nu);cos(nu);0];
    cO=cos(RAAN);sO=sin(RAAN);ci=cos(inc);si=sin(inc);
    R=[cO,-sO*ci,sO*si;sO,cO*ci,-cO*si;0,si,ci];
    r=R*r_pf; v=R*v_pf;
end