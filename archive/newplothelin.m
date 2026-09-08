load('drag_sail_attitude_results.mat')

mu = 3.986004418e14; Re = 6.3781e6;
year_ep=2025; month_ep=6; day_ep=21; hour_utc=10;
f107a=150; f107=150; ap=15;
inc=deg2rad(98); RAAN=deg2rad(90);

clrs_alt   = lines(3);
alt_labels = {'350 km', '450 km', '650 km'};

fig = figure('Position',[60 60 700 450]);
hold on;

for alt_i = 1:3
    alt    = alt_vec(alt_i) * 1e3;
    a_orb  = Re + alt;
    T_orb  = 2*pi*sqrt(a_orb^3/mu);
    n_orb  = 2*pi/T_orb;
    t_out  = ref_results(alt_i,1).t;
    doy0   = day(datetime(year_ep,month_ep,day_ep),'dayofyear');

    rho_t = zeros(size(t_out));
    for k = 1:length(t_out)
        t  = t_out(k);
        nu = mod(n_orb*t, 2*pi);
        r_pf = a_orb*[cos(nu);sin(nu);0];
        cO=cos(RAAN); sO=sin(RAAN); ci=cos(inc); si=sin(inc);
        R = [cO,-sO*ci,sO*si; sO,cO*ci,-cO*si; 0,si,ci];
        r_eci = R*r_pf;
        utc_dt  = datetime(year_ep,month_ep,day_ep,hour_utc,0,0)+seconds(t);
        utc_vec = [utc_dt.Year,utc_dt.Month,utc_dt.Day, ...
                   utc_dt.Hour,utc_dt.Minute,utc_dt.Second];
        lla   = eci2lla(r_eci', utc_vec);
        doy   = doy0 + t/86400;
        utc_s = mod(hour_utc*3600+t, 86400);
        [~,rho_arr] = atmosnrlmsise00(lla(3),lla(1),lla(2), ...
                                      year_ep,doy,utc_s,f107a,f107,ap);
        rho_t(k) = rho_arr(6);
    end

    plot(t_out/60, rho_t, 'Color',clrs_alt(alt_i,:), 'LineWidth',1.5, ...
         'DisplayName',alt_labels{alt_i});
end

xlabel('Time (min)','Interpreter','latex');
ylabel('$\rho$ (kg m$^{-3}$)','Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex', ...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
savefig(fig,'plot_density_vs_time.fig');
exportgraphics(fig,'plot_density_vs_time.png','Resolution',300);
fprintf('Saved\n');