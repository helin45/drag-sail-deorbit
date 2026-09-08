%% Orbital lifetime vs altitude plot
clear; clc;

mu  = 3.986004418e14;
Re  = 6.3781e6;
m   = 5.8;
year_ep=2025; month_ep=6; day_ep=21; hour_utc=10;
doy0 = day(datetime(year_ep,month_ep,day_ep),'dayofyear');

load('drag_sail_attitude_results.mat')

h_vec_km = 200:10:800;

%% Solar activity scenarios (Option A)
%  f107a = 81-day average, f107 = daily, ap = geomagnetic index.
%  Kept as physically consistent trios (activity rises together).
solar(1) = struct('name','Solar min',  'f107a', 70, 'f107', 70, 'ap',  5);
solar(2) = struct('name','Solar mean', 'f107a',150, 'f107',150, 'ap', 15);
solar(3) = struct('name','Solar max',  'f107a',220, 'f107',220, 'ap', 30);
n_solar   = numel(solar);
mean_idx  = 2;   % scenario used for the existing beta / alpha plots

%% Build orbit-averaged density lookup, per solar scenario
fprintf('Computing orbit-averaged density (%d solar scenarios)...\n', n_solar);
rho_lut_all = zeros(n_solar, length(h_vec_km));
n_pts   = 50;
inc     = deg2rad(98);
RAAN    = deg2rad(90);

for s = 1:n_solar
    fprintf('\n--- %s (F10.7=%d, ap=%d) ---\n', ...
            solar(s).name, solar(s).f107, solar(s).ap);
    for k = 1:length(h_vec_km)
        alt_m = h_vec_km(k)*1e3;
        a_orb = Re + alt_m;
        T_orb = 2*pi*sqrt(a_orb^3/mu);
        n_orb = 2*pi/T_orb;

        rho_sum = 0;
        for j = 1:n_pts
            t  = (j-1)/n_pts * T_orb;
            nu = mod(n_orb*t, 2*pi);

            r_pf = a_orb*[cos(nu); sin(nu); 0];
            cO=cos(RAAN); sO=sin(RAAN); ci=cos(inc); si=sin(inc);
            R=[cO,-sO*ci,sO*si; sO,cO*ci,-cO*si; 0,si,ci];
            r_eci = R*r_pf;

            utc_dt  = datetime(year_ep,month_ep,day_ep,hour_utc,0,0)+seconds(t);
            utc_vec = [utc_dt.Year,utc_dt.Month,utc_dt.Day,...
                       utc_dt.Hour,utc_dt.Minute,utc_dt.Second];
            lla = eci2lla(r_eci', utc_vec);

            doy   = doy0 + t/86400;
            utc_s = mod(hour_utc*3600+t, 86400);

            [~, nd] = atmosnrlmsise00(lla(3),lla(1),lla(2),...
                                       year_ep,doy,utc_s,...
                                       solar(s).f107a, solar(s).f107, solar(s).ap);
            rho_sum = rho_sum + nd(6);   % total mass density (kg/m^3)
        end
        rho_lut_all(s,k) = rho_sum / n_pts;
    end
    fprintf('  done (%d km: %.3e -> %d km: %.3e kg/m3)\n', ...
            h_vec_km(1), rho_lut_all(s,1), h_vec_km(end), rho_lut_all(s,end));
end
rho_lut = rho_lut_all(mean_idx,:);   % existing plots below use the mean scenario
fprintf('Density profiles done\n');

%% ISS validation (uses mean scenario)
beta_ISS   = 100;
h_range_v  = (120e3:1e3:400e3);
rho_v      = interp1(h_vec_km*1e3, rho_lut, h_range_v, 'pchip');
a_v        = Re + h_range_v;
v_v        = sqrt(mu ./ a_v);
T_ISS      = trapz(h_range_v, beta_ISS./(rho_v.*v_v.*a_v)) / (365.25*24*3600);
fprintf('\nISS validation (beta=100, h0=400km): %.2f years\n', T_ISS);
fprintf('Expected: 1-2 years\n\n');

%% Beta values from actual results (proper ballistic coeff: m/(Cd*A_ref), kg/m^2)
A_ref_phi = arrayfun(@(pp) ref_geom.A_ref(ref_geom.phi_deg==pp), phi_vec);  % [1 x n_phi]
Cd_slice  = squeeze(Cd_bar_3D(2,:,:));        % [n_accom x n_phi] at 450 km
beta_all  = m ./ (Cd_slice .* A_ref_phi);     % broadcasts A_ref per phi
beta_all  = beta_all(:);
beta_min = min(beta_all); beta_max = max(beta_all); step = 0.05;
beta_selected = [];
current = beta_min;
while current <= beta_max + step
    [~,idx] = min(abs(beta_all - current));
    val = beta_all(idx);
    if isempty(beta_selected) || abs(val - beta_selected(end)) > step/2
        beta_selected(end+1) = val;
    end
    current = current + step;
end
beta_vals   = unique(round(beta_selected,3));
beta_labels = arrayfun(@(x) sprintf('$\\beta = %.3f$ kg m$^{-2}$',x), beta_vals,'UniformOutput',false);

%% Compute decay times Plot 1
fprintf('Computing decay times for Plot 1...\n');
T_decay_plot = zeros(length(h_vec_km), length(beta_vals));
for hi = 1:length(h_vec_km)
    h0_m = h_vec_km(hi)*1e3;
    for bi = 1:length(beta_vals)
        beta    = beta_vals(bi);
        h_range = (120e3:1e3:h0_m);
        if length(h_range) < 2; continue; end
        rho_i = interp1(h_vec_km*1e3, rho_lut, h_range, 'pchip');
        a_r   = Re + h_range;
        v_r   = sqrt(mu ./ a_r);
        T_s   = trapz(h_range, beta./(rho_i.*v_r.*a_r));
        T_decay_plot(hi,bi) = T_s / (365.25*24*3600);
    end
end
fprintf('Range: %.2f to %.1f years\n', min(T_decay_plot(:)), max(T_decay_plot(:)));

%% Plot 1
clrs = parula(length(beta_vals));
fig1 = figure('Position',[60 60 750 550]);
hold on;
for bi = 1:length(beta_vals)
    plot(T_decay_plot(:,bi), h_vec_km, ...
         'Color',clrs(bi,:),'LineWidth',1.8,...
         'DisplayName',beta_labels{bi});
end
set(gca,'XScale','log');
valid_T = T_decay_plot(T_decay_plot > 0);
xmin = 10^floor(log10(min(valid_T)));
xmax = 10^ceil(log10(max(valid_T)));
xlim([xmin xmax]); ylim([200 800]);
decades = log10(xmin):log10(xmax);
tick_labels = cell(size(decades));
for k = 1:length(decades)
    d = decades(k);
    if d >= 0
        tick_labels{k} = sprintf('%d', 10^d);
    else
        tick_labels{k} = sprintf('%.*f', -d, 10^d);
    end
end
xticks(10.^decades);
xticklabels(tick_labels);
text_x = 10^(log10(xmin) + 0.7*(log10(xmax)-log10(xmin)));
xline(5,'r--','LineWidth',2,'DisplayName','5-year guideline');
yline(350,'k:','LineWidth',1,'HandleVisibility','off');
yline(450,'k:','LineWidth',1,'HandleVisibility','off');
yline(650,'k:','LineWidth',1,'HandleVisibility','off');
text(text_x,355,'350 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(text_x,455,'450 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(text_x,655,'650 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
xlabel('Estimated Orbital Decay Time (years)','Interpreter','latex');
ylabel('Initial Orbital Altitude $h_0$ (km)','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex',...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
grid on;
% Mark where each beta curve crosses the 5-year guideline (if within altitude range)
hcol = h_vec_km(:);
for bi = 1:length(beta_vals)
    Tb = T_decay_plot(:,bi); good = Tb > 0;
    if any(good) && min(Tb(good)) <= 5 && max(Tb) >= 5
        [Tu,iu] = unique(Tb(good)); hg = hcol(good);
        h_cross = interp1(Tu, hg(iu), 5, 'pchip');
        plot(5, h_cross,'o','MarkerSize',6,'MarkerFaceColor',clrs(bi,:), ...
             'MarkerEdgeColor','k','LineWidth',0.75,'HandleVisibility','off');
        text(5*1.15, h_cross, sprintf('%.0f km', h_cross), ...
             'Interpreter','latex','FontSize',7,'VerticalAlignment','middle');
    end
end
lg = legend('Location','northwest','FontSize',8);
set(lg,'Interpreter','latex');
savefig(fig1,'plot_lifetime_vs_altitude.fig');
exportgraphics(fig1,'plot_lifetime_vs_altitude.png','Resolution',300);
fprintf('Saved plot 1\n');

%% Beta from accommodation at phi=75
ref_si   = find(phi_vec==75);
ref_ac   = find(accom_vals==0.75);
A_ref_75 = ref_geom.A_ref(ref_si);
accom_labels2 = arrayfun(@(x) sprintf('$\\alpha = %.2f$',x), accom_vals,'UniformOutput',false);
beta_from_accom = zeros(1,length(accom_vals));
for ac = 1:length(accom_vals)
    Cd_bar_75 = mean(squeeze(Cd_bar_3D(:,ac,ref_si)));
    beta_from_accom(ac) = m / (Cd_bar_75 * A_ref_75);
end

%% Compute decay times Plot 2
fprintf('Computing decay times for Plot 2...\n');
T_decay_accom = zeros(length(h_vec_km), length(accom_vals));
for hi = 1:length(h_vec_km)
    h0_m = h_vec_km(hi)*1e3;
    for ac = 1:length(accom_vals)
        beta    = beta_from_accom(ac);
        h_range = (120e3:1e3:h0_m);
        if length(h_range) < 2; continue; end
        rho_i = interp1(h_vec_km*1e3, rho_lut, h_range, 'pchip');
        a_r   = Re + h_range;
        v_r   = sqrt(mu ./ a_r);
        T_s   = trapz(h_range, beta./(rho_i.*v_r.*a_r));
        T_decay_accom(hi,ac) = T_s / (365.25*24*3600);
    end
end

%% Plot 2
clrs2 = flip(parula(length(accom_vals)));
fig2 = figure('Position',[60 60 750 550]);
hold on;
for ac = 1:length(accom_vals)
    plot(T_decay_accom(:,ac), h_vec_km, ...
         'Color',clrs2(ac,:),'LineWidth',1.8,...
         'DisplayName',accom_labels2{ac});
end
set(gca,'XScale','log');
xlim([0.1 300]); ylim([200 800]);
xticks([0.1 0.5 1 2 5 10 25 50 100 200 300]);
xticklabels({'0.1','0.5','1','2','5','10','25','50','100','200','300'});
xline(5,'r--','LineWidth',2,'DisplayName','5-year guideline');
yline(350,'k:','LineWidth',1,'HandleVisibility','off');
yline(450,'k:','LineWidth',1,'HandleVisibility','off');
yline(650,'k:','LineWidth',1,'HandleVisibility','off');
text(200,355,'350 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(200,455,'450 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(200,655,'650 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
title('$\phi = 75^\circ$, varying $\alpha$','Interpreter','latex');
xlabel('Estimated Orbital Decay Time (years)','Interpreter','latex');
ylabel('Initial Orbital Altitude $h_0$ (km)','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex',...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
grid on;
lg2 = legend('Location','northwest','FontSize',8);
set(lg2,'Interpreter','latex');
savefig(fig2,'plot_lifetime_vs_altitude_accom.fig');
exportgraphics(fig2,'plot_lifetime_vs_altitude_accom.png','Resolution',300);
fprintf('Saved plot 2\n');

%% =========================================================================
%  PLOT 3 - Solar-activity sensitivity (Option A)
%  Representative sail: phi=75 deg, alpha=0.75. min / mean / max solar,
%  with a shaded band spanning the solar-cycle spread.
% =========================================================================
Cd_rep   = mean(squeeze(Cd_bar_3D(:,ref_ac,ref_si)));   % phi=75, alpha=0.75
beta_rep = m / (Cd_rep * A_ref_75);
fprintf('\nPlot 3: representative beta = %.3f kg/m2 (phi=75, alpha=0.75)\n', beta_rep);

T_decay_rep = zeros(length(h_vec_km), n_solar);
for s = 1:n_solar
    rho_s = rho_lut_all(s,:);
    for hi = 1:length(h_vec_km)
        h0_m    = h_vec_km(hi)*1e3;
        h_range = (120e3:1e3:h0_m);
        if length(h_range) < 2; continue; end
        rho_i = interp1(h_vec_km*1e3, rho_s, h_range, 'pchip');
        a_r   = Re + h_range;
        v_r   = sqrt(mu ./ a_r);
        T_decay_rep(hi,s) = trapz(h_range, beta_rep./(rho_i.*v_r.*a_r)) ...
                            / (365.25*24*3600);
    end
end

% Solar min = lowest density = longest lifetime (rightmost curve)
% Solar max = highest density = shortest lifetime (leftmost curve)
[~, quiet_idx]  = min([solar.f107]);
[~, active_idx] = max([solar.f107]);
solar_clrs = [0.15 0.30 0.75;   % min  (blue)
              0.20 0.55 0.25;   % mean (green)
              0.80 0.35 0.10];  % max  (orange)

fig3 = figure('Position',[60 60 750 550]);
hold on;

% Shaded solar-cycle band (between active and quiet extremes)
valid = T_decay_rep(:,quiet_idx) > 0 & T_decay_rep(:,active_idx) > 0;
hb = h_vec_km(valid).';
patch_x = [T_decay_rep(valid,active_idx); flipud(T_decay_rep(valid,quiet_idx))];
patch_y = [hb; flipud(hb)];
fill(patch_x, patch_y, [0.5 0.5 0.5], 'FaceAlpha',0.15, ...
     'EdgeColor','none','HandleVisibility','off');

for s = 1:n_solar
    plot(T_decay_rep(:,s), h_vec_km, ...
         'Color',solar_clrs(s,:),'LineWidth',2.0, ...
         'DisplayName',sprintf('%s ($F_{10.7}=%d$)', solar(s).name, solar(s).f107));
end
set(gca,'XScale','log');
valid_T3 = T_decay_rep(T_decay_rep > 0);
xmin3 = 10^floor(log10(min(valid_T3)));
xmax3 = 10^ceil(log10(max([valid_T3(:); 5])));   % ensure 5-yr line is visible
xlim([xmin3 xmax3]); ylim([200 800]);
decades3 = log10(xmin3):log10(xmax3);
tl3 = cell(size(decades3));
for k = 1:length(decades3)
    d = decades3(k);
    if d >= 0, tl3{k} = sprintf('%d',10^d); else, tl3{k} = sprintf('%.*f',-d,10^d); end
end
xticks(10.^decades3); xticklabels(tl3);
tx3 = 10^(log10(xmin3) + 0.55*(log10(xmax3)-log10(xmin3)));
xline(5,'r--','LineWidth',2,'DisplayName','5-year guideline');
yline(350,'k:','LineWidth',1,'HandleVisibility','off');
yline(450,'k:','LineWidth',1,'HandleVisibility','off');
yline(650,'k:','LineWidth',1,'HandleVisibility','off');
text(tx3,355,'350 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(tx3,455,'450 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
text(tx3,655,'650 km','Interpreter','latex','FontSize',9,'HorizontalAlignment','right');
title(sprintf('$\\phi = 75^\\circ$, $\\alpha = 0.75$ ($\\beta = %.3f$ kg m$^{-2}$)', beta_rep), ...
      'Interpreter','latex');
xlabel('Estimated Orbital Decay Time (years)','Interpreter','latex');
ylabel('Initial Orbital Altitude $h_0$ (km)','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex',...
    'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
grid on;
% Mark where each solar curve crosses the 5-year guideline (if within altitude range)
hcol = h_vec_km(:);
for s = 1:n_solar
    Ts = T_decay_rep(:,s); good = Ts > 0;
    if any(good) && min(Ts(good)) <= 5 && max(Ts) >= 5
        [Tu,iu] = unique(Ts(good)); hg = hcol(good);
        h_cross = interp1(Tu, hg(iu), 5, 'pchip');
        plot(5, h_cross,'o','MarkerSize',7,'MarkerFaceColor',solar_clrs(s,:), ...
             'MarkerEdgeColor','k','LineWidth',0.75,'HandleVisibility','off');
        text(5*1.15, h_cross, sprintf('%.0f km', h_cross), ...
             'Interpreter','latex','FontSize',8,'VerticalAlignment','middle');
    end
end
lg3 = legend('Location','northwest','FontSize',9);
set(lg3,'Interpreter','latex');
savefig(fig3,'plot_lifetime_vs_altitude_solar.fig');
exportgraphics(fig3,'plot_lifetime_vs_altitude_solar.png','Resolution',300);
fprintf('Saved plot 3\n');