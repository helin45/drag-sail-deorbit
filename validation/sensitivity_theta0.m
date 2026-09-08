%% test_theta0.m
cd(fileparts(fileparts(mfilename('fullpath'))));  % run from repo root (loads mc_workspace.mat)
load('mc_workspace.mat');

N_test = 50;
rng(42);
phi0_test = deg2rad(-10 + 20*rand(N_test,1));
CD_test = zeros(N_test,1);
p.rho_scale = 1.0;
set_sails(sails_c);
opts = odeset('RelTol',ode_reltol,'AbsTol',ode_abstol);
t_vec = linspace(0,T_orb,3000)';

for mi = 1:N_test
    [t_out,X_out] = ode45(@(t,X) attitude_ode(t,X,p), t_vec, ...
                          [phi0_test(mi); 0], opts);
    CD_test(mi) = compute_CDbar(t_out, rad2deg(X_out(:,1)), sails_c, 2.0);
    if mod(mi,10)==0; fprintf('%d/%d | CD=%.4f\n',mi,N_test,CD_test(mi)); end
end

figure;
scatter(rad2deg(phi0_test), CD_test, 40, 'b', 'filled');
xlabel('$\theta_0$ (deg)','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',12,'Box','on');
fprintf('Std = %.6f\n', std(CD_test));
fprintf('CV  = %.4f%%\n', std(CD_test)/mean(CD_test)*100);

%% Local functions
function CD_bar=compute_CDbar(t_out,phi_deg_t,sails_c,d_theta_deg)
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
    CD_bar=sum(Cd_bins.*pdf_SAM)*d_theta_deg;
end
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
    rho=rho_arr(6)*p.rho_scale;
    q=0.5*rho*p.v_circ^2;
    sails=set_sails(); tau_total=0; phi_deg=rad2deg(phi);
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