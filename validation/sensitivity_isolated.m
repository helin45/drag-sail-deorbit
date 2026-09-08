
cd(fileparts(fileparts(mfilename('fullpath'))));  % run from repo root (loads mc_workspace.mat)
load('mc_workspace.mat');
t_vec = linspace(0,T_orb,3000)';
opts  = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
N_test = 50;
rng(42);
rho_scales = [0.85, 1.00, 1.15];
rho_labels = {'$s_\rho = 0.85$','$s_\rho = 1.00$','$s_\rho = 1.15$'};
clrs       = [0.2 0.4 0.8; 0.2 0.7 0.2; 0.8 0.2 0.2];

sail_alpha = sails_c(1).alpha;
sail_Cm    = sails_c(1).Cm_pitch;
sail_Cd    = sails_c(1).Cd;
sail_Iyy   = sails_c(1).Iyy;
sail_Aref  = sails_c(1).A_ref;
sail_Lref  = sails_c(1).L_ref;

%% =========================================================================
%  TEST 1: Isolated theta0, dphi0=0, 3 densities
% =========================================================================
phi0_test = deg2rad(-10 + 20*rand(N_test,1));
CD_t1     = zeros(N_test,3);

for ri = 1:3
    rs    = rho_scales(ri);
    p_loc = build_params(p, rs, sail_alpha, sail_Cm, sail_Cd, sail_Iyy, sail_Aref, sail_Lref);
    CD_tmp = zeros(N_test,1);
    parfor mi = 1:N_test
        [t_out,X_out] = ode45(@(t,X) attitude_ode_direct(t,X,p_loc), t_vec, ...
                              [phi0_test(mi); 0], opts);
        CD_tmp(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_alpha, sail_Cd, 2.0);
    end
    CD_t1(:,ri) = CD_tmp;
    fprintf('Test1 | s_rho=%.2f | Mean=%.4f | CV=%.4f%%\n', rs, mean(CD_tmp), std(CD_tmp)/mean(CD_tmp)*100);
    save('isolated_test1_partial.mat','CD_t1','phi0_test','rho_scales');
    fprintf('  Saved test1 partial (ri=%d)\n', ri);
end
save('isolated_test1.mat','CD_t1','phi0_test','rho_scales');
fprintf('Saved isolated_test1.mat\n');

%% =========================================================================
%  TEST 2: Isolated density, theta0=phi0 fixed, dphi0=0
% =========================================================================

phi0_fixed = deg2rad(5);  % explicitly define for parfor workers
rho_test = 0.85 + 0.30*rand(N_test,1);
CD_t2    = zeros(N_test,1);
parfor mi = 1:N_test
    p_loc = build_params(p, rho_test(mi), sail_alpha, sail_Cm, sail_Cd, sail_Iyy, sail_Aref, sail_Lref);
    [t_out,X_out] = ode45(@(t,X) attitude_ode_direct(t,X,p_loc), t_vec, [phi0_fixed; 0], opts);
    CD_t2(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_alpha, sail_Cd, 2.0);
end

%% =========================================================================
%  TEST 3: Isolated dphi0, theta0=phi0 fixed, 3 densities
% =========================================================================

dphi0_test = deg2rad(-2 + 4*rand(N_test,1));
CD_t3      = zeros(N_test,3);
for ri = 1:3
    rs    = rho_scales(ri);
    p_loc = build_params(p, rs, sail_alpha, sail_Cm, sail_Cd, sail_Iyy, sail_Aref, sail_Lref);
    CD_tmp = zeros(N_test,1);
    parfor mi = 1:N_test
        [t_out,X_out] = ode45(@(t,X) attitude_ode_direct(t,X,p_loc), t_vec, ...
                              [phi0_fixed; dphi0_test(mi)], opts);
        CD_tmp(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sail_alpha, sail_Cd, 2.0);
    end
    CD_t3(:,ri) = CD_tmp;
    fprintf('Test3 | s_rho=%.2f | Mean=%.4f | CV=%.4f%%\n', rs, mean(CD_tmp), std(CD_tmp)/mean(CD_tmp)*100);
    save('isolated_test3_partial.mat','CD_t3','dphi0_test','rho_scales');
    fprintf('  Saved test3 partial (ri=%d)\n', ri);
end
%% plotting
f1 = figure('Position',[100 100 600 450]);
hold on;
for ri=1:3
    scatter(rad2deg(phi0_test), CD_t1(:,ri), 40, 'filled',...
        'MarkerFaceColor',clrs(ri,:),'DisplayName',rho_labels{ri});
end
xlabel('$\theta_0$ (deg)','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
legend('Interpreter','latex','Location','best');
set(gca,'TickLabelInterpreter','latex','FontSize',12,'Box','on','LineWidth',1.2,...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.15);
exportgraphics(f1,'isolated_theta0_3rho.png','Resolution',300);
fprintf('Saved isolated_theta0_3rho.png\n');

f2 = figure('Position',[100 100 600 450]);
scatter(rho_test, CD_t2, 40, [0.2 0.4 0.8], 'filled');
xlabel('$s_\rho$','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',12,'Box','on','LineWidth',1.2,...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.15);
exportgraphics(f2,'isolated_rho.png','Resolution',300);
fprintf('Saved isolated_rho.png\n');

f3 = figure('Position',[100 100 600 450]);
hold on;
for ri=1:3
    scatter(rad2deg(dphi0_test), CD_t3(:,ri), 40, 'filled',...
        'MarkerFaceColor',clrs(ri,:),'DisplayName',rho_labels{ri});
end
xlabel('$\dot{\theta}_0$ (deg/s)','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
legend('Interpreter','latex','Location','south');
set(gca,'TickLabelInterpreter','latex','FontSize',12,'Box','on','LineWidth',1.2,...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.15);
exportgraphics(f3,'isolated_dphi0_3rho.png','Resolution',300);
fprintf('Saved isolated_dphi0_3rho.png\n');

%% Final save
save('isolated_tests_all.mat','CD_t1','CD_t2','CD_t3',...
     'phi0_test','rho_test','dphi0_test','rho_scales');
fprintf('Saved isolated_tests_all.mat\n');

%% functions
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