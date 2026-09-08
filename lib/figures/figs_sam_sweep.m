%% =========================================================================
%  ATTITUDE SAM SWEEP — ALL PLOTS (corrected)
%  Uses Cd_bar_fixed | data-driven ylims | no figure titles
% =========================================================================
clear; clc; close all;

load('attitude_SAM_sweep_FINAL.mat');   % gets Cd_bar_fixed, ref_results
load('longrun_SAM_sweep_FINAL.mat');    % gets Cd_SAM_long, Cd_random_long, etc.

% Override the old single-orbit Cd_bar_fixed with the 15-orbit SAM results
Cd_bar_fixed   = Cd_SAM_long;
Cd_random_3D   = Cd_random_long;
Cd_notumble_3D = Cd_notumble_long;

% Recompute ylim from new data
cd_min  = floor(min(Cd_bar_fixed(:)) * 20) / 20;
cd_max  = ceil( max(Cd_bar_fixed(:)) * 20) / 20;
cd_ylim = [cd_min - 0.05, cd_max + 0.05];

%% NOTE for Plot A5 (decay time colourmap):
%  load('decay_results.mat','decay_years_3D') before running

alt_vec    = [350, 450, 650];
phi_vec    = 45:5:85;
accom_vals = 0.50:0.05:1.00;
accom_vec  = {'0p50','0p55','0p60','0p65','0p70','0p75', ...
              '0p80','0p85','0p90','0p95','1p00'};

n_alts  = length(alt_vec);
n_phi   = length(phi_vec);
n_accom = length(accom_vals);

ref_ac_idx  = find(strcmp(accom_vec,'0p75'));
ref_phi_idx = find(phi_vec == 85);

cmap_accom = parula(n_accom);
cmap_phi   = turbo(n_phi);
alt_colors = [0.15 0.45 0.85;
              0.10 0.72 0.40;
              0.85 0.25 0.10];

cd_min  = floor(min(Cd_bar_fixed(:)) * 20) / 20;
cd_max  = ceil( max(Cd_bar_fixed(:)) * 20) / 20;
cd_ylim = [cd_min - 0.05, cd_max + 0.05];

mu = 3.986004418e14;
Re = 6.3781e6;


%% =========================================================================
%  PRE-COMPUTE: NRLMSISE density along orbit
% =========================================================================

year_ep = 2025; month_ep = 1; day_ep = 15; hour_utc = 0;
f107a = 150; f107 = 150; ap = 15;
n_orbit_pts = 360;

density_data(n_alts) = struct('rho',[],'nu_deg',[]);

fprintf('Computing NRLMSISE density profiles...\n');
for alt_i = 1:n_alts
    alt   = alt_vec(alt_i);
    a_orb = Re + alt*1e3;
    T_orb = 2*pi*sqrt(a_orb^3/mu);
    t_vec = linspace(0, T_orb, n_orbit_pts);
    rho_vec = zeros(1, n_orbit_pts);
    for j = 1:n_orbit_pts
        t     = t_vec(j);
        nu    = mod(2*pi/T_orb*t, 2*pi);
        r_eci = a_orb * [cos(nu); sin(nu); 0];
        utc_dt  = datetime(year_ep,month_ep,day_ep,hour_utc,0,0) + seconds(t);
        utc_vec = [utc_dt.Year utc_dt.Month utc_dt.Day ...
                   utc_dt.Hour utc_dt.Minute utc_dt.Second];
        lla   = eci2lla(r_eci', utc_vec);
        doy   = day(utc_dt,'dayofyear');
        utc_s = hour(utc_dt)*3600 + minute(utc_dt)*60 + second(utc_dt);
        [~,rho_arr] = atmosnrlmsise00(lla(3),lla(1),lla(2), ...
                      year_ep,doy,utc_s,f107a,f107,ap);
        rho_vec(j) = rho_arr(6);
    end
    density_data(alt_i).rho    = rho_vec;
    density_data(alt_i).nu_deg = rad2deg(linspace(0,2*pi,n_orbit_pts));
end
fprintf('Done.\n');

%% =========================================================================
%  SMOOTH INTERPOLATION GRID for colour plots
% =========================================================================

% Replace the PRE-COMPUTE smooth grid section with this:

% Replace the smooth grid pre-compute section with this:

n_interp   = 200;
accom_fine = linspace(accom_vals(1), accom_vals(end), n_interp);
phi_fine   = linspace(phi_vec(1),    phi_vec(end),    n_interp);
[AC_fine, PHI_fine] = meshgrid(accom_fine, phi_fine);

% Flat query points for griddata
[AC_pts, PHI_pts] = meshgrid(accom_vals, phi_vec);
ac_flat  = AC_pts(:);
phi_flat = PHI_pts(:);

%% =========================================================================
%  FIGURE 1: Cd_bar vs phi
% =========================================================================

figure('Position',[100 100 1400 450]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i); hold on; box on; grid on;
    for ac = 1:n_accom
        plot(phi_vec, squeeze(Cd_bar_fixed(alt_i,ac,:)), '-o', ...
             'Color', cmap_accom(ac,:), 'LineWidth', 1.2, 'MarkerSize', 5)
    end
    xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    ylabel('$\bar{C}_D$','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([43 87]); ylim(cd_ylim); xticks(phi_vec)
    if alt_i == n_alts
        legend(arrayfun(@(a) sprintf('$\\alpha = %.2f$',a), accom_vals, ...
               'UniformOutput',false),'Interpreter','latex', ...
               'Location','northwestoutside','FontSize',7)
    end
end

%% =========================================================================
%  FIGURE 2: Cd_bar vs accom
% =========================================================================

figure('Position',[100 100 1400 450]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i); hold on; box on; grid on;
    for si = 1:n_phi
        plot(accom_vals, squeeze(Cd_bar_fixed(alt_i,:,si)), '-o', ...
             'Color', cmap_phi(si,:), 'LineWidth', 1.2, 'MarkerSize', 5)
    end
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('$\bar{C}_D$','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([0.48 1.02]); ylim(cd_ylim); xticks(accom_vals)
    if alt_i == n_alts
        legend(arrayfun(@(p) sprintf('$\\phi = %d^\\circ$',p), phi_vec, ...
               'UniformOutput',false),'Interpreter','latex', ...
               'Location','northeastoutside','FontSize',7)
    end
end

%% =========================================================================
%  FIGURE 3: Pitch transient — coloured by phi (cmap_phi), panel = altitude
% =========================================================================

figure('Position',[100 100 700 750]);
for alt_i = 1:n_alts
    subplot(n_alts,1,alt_i); hold on; box on; grid on;
    all_theta_alt = [];
    for si = 1:n_phi
        t         = ref_results(alt_i, ref_ac_idx, si).t;
        phi_deg_t = wrapTo180(ref_results(alt_i, ref_ac_idx, si).phi_deg);
        plot(t/60, phi_deg_t, 'Color', cmap_phi(si,:), 'LineWidth', 0.8)
        all_theta_alt = [all_theta_alt; phi_deg_t];
    end
    yline(0,'k--','LineWidth',0.5)
    xlabel('Time (min)','Interpreter','latex')
    ylabel('$\theta$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([0 ref_results(alt_i,ref_ac_idx,1).T_orb_min])
    ylim([floor(min(all_theta_alt)*2)/2-1, ceil(max(all_theta_alt)*2)/2+1])
    if alt_i == 1
        legend(arrayfun(@(p) sprintf('$\\phi=%d^\\circ$',p), phi_vec, ...
               'UniformOutput',false),'Interpreter','latex', ...
               'NumColumns',3,'Location','northeast','FontSize',6)
    end
end

%% =========================================================================
%  FIGURE 4: 5-min zoom — coloured by ALTITUDE, panel = altitude + phi
% =========================================================================

phi_subset  = [45, 65, 85];
phi_sub_idx = arrayfun(@(p) find(phi_vec==p), phi_subset);

figure('Position',[100 100 900 850]);
sp = 1;
for si = 1:length(phi_subset)
    for alt_i = 1:n_alts
        subplot(length(phi_subset),n_alts,sp)
        hold on; box on; grid on;
        t         = ref_results(alt_i, ref_ac_idx, phi_sub_idx(si)).t;
        phi_deg_t = ref_results(alt_i, ref_ac_idx, phi_sub_idx(si)).phi_deg;
        idx5      = t <= 5*60;
        theta5    = wrapTo180(phi_deg_t(idx5));
        plot(t(idx5)/60, theta5, 'Color', alt_colors(alt_i,:), 'LineWidth', 1.2)
        yline(0,'k--','LineWidth',0.5)
        xlabel('Time (min)','Interpreter','latex')
        ylabel('$\theta$ (deg)','Interpreter','latex')
        title(sprintf('$h=%d$ km, $\\phi=%d^\\circ$', ...
              alt_vec(alt_i), phi_subset(si)),'Interpreter','latex')
        ylim([floor(min(theta5)*2)/2-0.5, ceil(max(theta5)*2)/2+0.5])
        sp = sp+1;
    end
end

%% =========================================================================
%  FIGURE 5: SAM PDF phi=85, alpha=0.75 — panel = altitude
% =========================================================================

figure('Position',[100 100 1200 380]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)
    pdf_vals = ref_results(alt_i, ref_ac_idx, ref_phi_idx).pdf_SAM;
    centres  = ref_results(alt_i, ref_ac_idx, ref_phi_idx).centres;
    bar(centres, pdf_vals, 'FaceColor', cmap_phi(ref_phi_idx,:), 'EdgeColor','none')
    xlabel('$\theta$ (deg)','Interpreter','latex')
    ylabel('PDF (deg$^{-1}$)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([-180 180]); xticks(-180:60:180)
    ylim([0, ceil(max(pdf_vals)*100)/100 + 0.001])
    grid on; box on
end

%% =========================================================================
%  FIGURE 6: Altitude sensitivity — panel = alpha value
% =========================================================================

accom_sel     = [0.50, 0.75, 1.00];
accom_sel_idx = arrayfun(@(a) find(abs(accom_vals-a)<1e-9), accom_sel);

figure('Position',[100 100 1200 380]);
for k = 1:length(accom_sel)
    subplot(1,3,k); hold on; box on; grid on;
    for alt_i = 1:n_alts
        plot(phi_vec, squeeze(Cd_bar_fixed(alt_i,accom_sel_idx(k),:)), '-o', ...
             'Color', alt_colors(alt_i,:), 'LineWidth', 1.4, 'MarkerSize', 5)
    end
    xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    ylabel('$\bar{C}_D$','Interpreter','latex')
    title(sprintf('$\\alpha = %.2f$', accom_sel(k)),'Interpreter','latex')
    xlim([43 87]); ylim(cd_ylim); xticks(phi_vec)
    legend(arrayfun(@(h) sprintf('$h = %d$ km',h), alt_vec, ...
           'UniformOutput',false),'Interpreter','latex','Location','northwest')
end

%% =========================================================================
%  NEW FIGURE A1: NRLMSISE density — all altitudes overlaid, one plot
% =========================================================================

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
for alt_i = 1:n_alts
    plot(density_data(alt_i).nu_deg, density_data(alt_i).rho, ...
         'Color', alt_colors(alt_i,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('$h = %d$ km', alt_vec(alt_i)))
end
xlabel('True Anomaly (deg)','Interpreter','latex')
ylabel('$\rho$ (kg m$^{-3}$)','Interpreter','latex')
xlim([0 360]); xticks(0:90:360)
legend('Interpreter','latex','Location','northeast')

%% =========================================================================
%  NEW FIGURE A2: PDFs for all phi, all altitudes — separate panels
%  n_alts rows x n_phi cols
% =========================================================================

figure('Position',[100 100 1800 500]);
sp = 1;
for alt_i = 1:n_alts
    for si = 1:n_phi
        subplot(n_alts, n_phi, sp)
        pdf_vals = ref_results(alt_i, ref_ac_idx, si).pdf_SAM;
        centres  = ref_results(alt_i, ref_ac_idx, si).centres;
        bar(centres, pdf_vals, 'FaceColor', cmap_phi(si,:), 'EdgeColor','none')
        xlim([-180 180]); xticks(-180:90:180)
        ylim([0, 0.012])
        if alt_i == n_alts
            xlabel('$\theta$ (deg)','Interpreter','latex','FontSize',7)
        end
        if si == 1
            ylabel(sprintf('%d km',alt_vec(alt_i)),'Interpreter','latex','FontSize',10)
        end
        title(sprintf('$\\phi=%d^\\circ$',phi_vec(si)),'Interpreter','latex','FontSize',7)
        set(gca,'FontSize',6)
        sp = sp + 1;
    end
end

%% =========================================================================
%  NEW FIGURE A3: PDFs for 3 accom coeffs — separate panels per alt and accom
%  3 rows (accom) x 3 cols (alt), fixed phi=85
% =========================================================================

figure('Position',[100 100 900 750]);
sp = 1;
for k = 1:length(accom_sel)
    for alt_i = 1:n_alts
        subplot(length(accom_sel), n_alts, sp)
        ac_idx   = accom_sel_idx(k);
        pdf_vals = ref_results(alt_i, ac_idx, ref_phi_idx).pdf_SAM;
        centres  = ref_results(alt_i, ac_idx, ref_phi_idx).centres;
        bar(centres, pdf_vals, 'FaceColor', alt_colors(alt_i,:), 'EdgeColor','none')
        xlim([-180 180]); xticks(-180:90:180)
        ylim([0, 0.012])
        if k == length(accom_sel)
            xlabel('$\theta$ (deg)','Interpreter','latex')
        end
        if alt_i == 1
            ylabel(sprintf('$\\alpha=%.2f$', accom_sel(k)),'Interpreter','latex')
        end
        title(sprintf('$h=%d$ km', alt_vec(alt_i)),'Interpreter','latex')
        sp = sp + 1;
    end
end

%% =========================================================================
%  NEW FIGURE A4: Smooth colourmap — beta_bar = 1/Cd_bar
% =========================================================================

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    beta_raw  = (1 ./ squeeze(Cd_bar_fixed(alt_i,:,:)))';
beta_fine = griddata(ac_flat, phi_flat, beta_raw(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, beta_fine, 40, 'LineColor','none')
    colormap(gca, turbo)
    cb = colorbar;
    cb.Label.String     = '$\bar{\beta}$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% =========================================================================
%  NEW FIGURE A5: Smooth colourmap — orbital decay time
% =========================================================================

if exist('decay_years_3D','var')

    figure('Position',[100 100 1400 420]);
    for alt_i = 1:n_alts
        subplot(1,3,alt_i)

        decay_raw  = squeeze(decay_years_3D(alt_i,:,:))';
    decay_fine = griddata(ac_flat, phi_flat, decay_raw(:), AC_fine, PHI_fine, 'cubic');

        contourf(accom_fine, phi_fine, decay_fine, 40, 'LineColor','none')
        colormap(gca, turbo)
        cb = colorbar;
        cb.Label.String      = 'Decay time (years)';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize    = 11;

        xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
        ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
        title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
        yticks(phi_vec); xticks(accom_vals); xtickangle(45)
    end

else
    fprintf('NOTE: decay_years_3D not found — skipping Figure A5.\n')
end

%% =========================================================================
%  NEW FIGURE A6: Bar chart — SAM vs random tumbling vs no tumbling
%  Fixed phi=75, alpha=0.75
% =========================================================================

phi_bar_idx = find(phi_vec == 75);
ac_bar_idx  = ref_ac_idx;

Cd_SAM_bar      = zeros(1,n_alts);
Cd_random_bar   = zeros(1,n_alts);
Cd_notumble_bar = zeros(1,n_alts);
for alt_i = 1:n_alts
    Cd_SAM_bar(alt_i)      = Cd_bar_fixed(alt_i, ac_bar_idx, phi_bar_idx);
    Cd_random_bar(alt_i)   = Cd_random_3D(alt_i, ac_bar_idx, phi_bar_idx);
    Cd_notumble_bar(alt_i) = Cd_notumble_3D(alt_i, ac_bar_idx, phi_bar_idx);
end

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
x = 1:n_alts; w = 0.25;
bar(x-w, Cd_SAM_bar,      w, 'FaceColor',[0.20 0.55 0.85],'DisplayName','SAM')
bar(x,   Cd_random_bar,   w, 'FaceColor',[0.10 0.75 0.40],'DisplayName','Random tumbling')
bar(x+w, Cd_notumble_bar, w, 'FaceColor',[0.85 0.30 0.10], ...
    'DisplayName','No tumbling ($\theta=0$)')
xticks(x)
xticklabels(arrayfun(@(h) sprintf('%d km',h), alt_vec,'UniformOutput',false))
ylabel('$C_D$','Interpreter','latex')
legend('Interpreter','latex','Location','northeast')
ylim([0, max([Cd_SAM_bar, Cd_random_bar, Cd_notumble_bar])*1.15])

%% =========================================================================
%  NEW FIGURE A7: Smooth colourmap — % difference SAM vs random tumbling
% =========================================================================

pct_diff_3D = 100 * (Cd_bar_fixed - Cd_random_3D) ./ Cd_random_3D;
pct_abs_max = ceil(max(abs(pct_diff_3D(:))));
clim_sym    = [-pct_abs_max, pct_abs_max];

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    pct_raw  = squeeze(pct_diff_3D(alt_i,:,:))';
pct_fine = griddata(ac_flat, phi_flat, pct_raw(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, pct_fine, 40, ['' ...
        'LineColor'],'none')
    colormap(gca, turbo)   % if you have aerospace toolbox
    clim(clim_sym)
    cb = colorbar;
    cb.Label.String      = '$\Delta \bar{C}_D$ (\%)';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

fprintf('\nAll figures generated.\n')

%% --- FIX: diverging colormap ---
n_c = 256;
cmap_div = [linspace(0.1,1,n_c/2)', linspace(0.1,1,n_c/2)', ones(n_c/2,1); ...
            ones(n_c/2,1), linspace(1,0.1,n_c/2)', linspace(1,0.1,n_c/2)'];

%% --- PLOT B1: Line graph Cd SAM vs random vs notumble ---
%  fixed alpha=0.75, all altitudes, varying phi

ac_fixed  = ref_ac_idx;   % alpha = 0.75
ls        = {'-o','-s','-^'};

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
for alt_i = 1:n_alts
    plot(phi_vec, squeeze(Cd_bar_fixed(alt_i, ac_fixed, :)), ...
         ls{alt_i}, 'Color', alt_colors(alt_i,:), 'LineWidth', 1.4, ...
         'MarkerSize', 6, 'DisplayName', sprintf('SAM $h=%d$ km', alt_vec(alt_i)))
end
for alt_i = 1:n_alts
    plot(phi_vec, squeeze(Cd_random_3D(alt_i, ac_fixed, :)), ...
         '--', 'Color', alt_colors(alt_i,:), 'LineWidth', 1.0, ...
         'DisplayName', sprintf('Random $h=%d$ km', alt_vec(alt_i)))
end
for alt_i = 1:n_alts
    plot(phi_vec, squeeze(Cd_notumble_3D(alt_i, ac_fixed, :)), ...
         ':', 'Color', alt_colors(alt_i,:), 'LineWidth', 1.0, ...
         'DisplayName', sprintf('No tumbling $h=%d$ km', alt_vec(alt_i)))
end
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
ylabel('$C_D$','Interpreter','latex')
xlim([43 87]); xticks(phi_vec)
legend('Interpreter','latex','Location','northwest','FontSize',7,'NumColumns',3)
title('$\alpha = 0.75$','Interpreter','latex')

%% --- PLOT B2: % difference SAM vs random — x=phi, y=accom, one panel per alt ---

pct_diff_3D = 100 * (Cd_bar_fixed - Cd_random_3D) ./ Cd_random_3D;

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    pct_raw  = squeeze(pct_diff_3D(alt_i,:,:))';  % [n_phi x n_accom]
    pct_fine = griddata(ac_flat, phi_flat, pct_raw(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, pct_fine, 40, 'LineColor','none')
    colormap(gca, cmap_div)
    pct_max = ceil(max(abs(pct_diff_3D(:))));
    clim([-pct_max pct_max])
    cb = colorbar;
    cb.Label.String      = '$(\bar{C}_{D,\mathrm{SAM}} - \bar{C}_{D,\mathrm{rand}}) / \bar{C}_{D,\mathrm{rand}}$ (\%)';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 9;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% --- PLOT B3: % difference SAM vs notumble colourmap ---

pct_diff_notumble = 100 * (Cd_bar_fixed - Cd_notumble_3D) ./ Cd_notumble_3D;

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    pct_raw  = squeeze(pct_diff_notumble(alt_i,:,:))';
    pct_fine = griddata(ac_flat, phi_flat, pct_raw(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, pct_fine, 40, 'LineColor','none')
    colormap(gca, cmap_div)
    pct_max2 = ceil(max(abs(pct_diff_notumble(:))));
    clim([-pct_max2 pct_max2])
    cb = colorbar;
    cb.Label.String      = '$(\bar{C}_{D,\mathrm{SAM}} - C_{D,0}) / C_{D,0}$ (\%)';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 9;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% --- PLOT B4: Cd_notumble colourmap (ram-facing performance) ---

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    raw_data = squeeze(Cd_notumble_3D(alt_i,:,:))';
    fine_data = griddata(ac_flat, phi_flat, raw_data(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, fine_data, 40, 'LineColor','none')
    colormap(gca, turbo)
    cb = colorbar;
    cb.Label.String      = '$C_{D}(\theta=0^\circ)$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% --- PLOT B5: Spider/radar chart — Cd SAM vs rand vs notumble for all phi ---
%  Fixed alpha=0.75, 450 km

alt_i_sel = 2;   % 450 km
ac_sel    = ref_ac_idx;

angles = deg2rad(phi_vec);
angles_closed = [angles, angles(1)];

SAM_vals      = squeeze(Cd_bar_fixed(alt_i_sel, ac_sel, :))';
rand_vals     = squeeze(Cd_random_3D(alt_i_sel, ac_sel, :))';
notumble_vals = squeeze(Cd_notumble_3D(alt_i_sel, ac_sel, :))';

SAM_c      = [SAM_vals, SAM_vals(1)];
rand_c     = [rand_vals, rand_vals(1)];
notumble_c = [notumble_vals, notumble_vals(1)];

figure('Position',[100 100 500 500]);
polarplot(angles_closed, SAM_c,      '-o', 'LineWidth', 1.5, 'DisplayName','SAM')
hold on
polarplot(angles_closed, rand_c,     '--s','LineWidth', 1.2, 'DisplayName','Random tumbling')
polarplot(angles_closed, notumble_c, ':^', 'LineWidth', 1.2, 'DisplayName','No tumbling')
legend('Interpreter','latex','Location','southoutside','NumColumns',3)
title('$h=450$ km, $\alpha=0.75$','Interpreter','latex')

%% --- PLOT B6: Absolute drag penalty (notumble - SAM) colourmap ---

drag_penalty = Cd_notumble_3D - Cd_bar_fixed;   % how much drag is lost to tumbling

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    pen_raw  = squeeze(drag_penalty(alt_i,:,:))';
    pen_fine = griddata(ac_flat, phi_flat, pen_raw(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, pen_fine, 40, 'LineColor','none')
    colormap(gca, turbo)
    cb = colorbar;
    cb.Label.String      = '$C_{D,0} - \bar{C}_{D,\mathrm{SAM}}$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% --- PLOT B7: Cd_bar_fixed colourmap (SAM result) ---

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    raw_data  = squeeze(Cd_bar_fixed(alt_i,:,:))';
    fine_data = griddata(ac_flat, phi_flat, raw_data(:), AC_fine, PHI_fine, 'cubic');

    contourf(accom_fine, phi_fine, fine_data, 40, 'LineColor','none')
    colormap(gca, turbo)
    cb = colorbar;
    cb.Label.String      = '$\bar{C}_{D,\mathrm{SAM}}$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

fprintf('\nExtra figures generated.\n')

%% =========================================================================
%  EXTENDED ANALYSIS PLOTS
%  Run after loading longrun_SAM_sweep_FINAL.mat and attitude_SAM_sweep_FINAL.mat
%  Requires: Cd_bar_fixed, Cd_random_3D, Cd_notumble_3D, alt_vec, phi_vec,
%            accom_vals, ac_flat, phi_flat, AC_fine, PHI_fine, alt_colors,
%            cmap_accom, cmap_phi, cmap_div from previous plot script
% =========================================================================

%% =========================================================================
%  C1: Tumbling penalty as fraction of notumble Cd
%  Line plot, x=phi, y=(Cd_notumble-Cd_SAM)/Cd_notumble, coloured by alpha
%  Fixed altitude = 450 km
% =========================================================================

alt_i_sel = 2;   % 450 km

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
for ac = 1:n_accom
    penalty = (squeeze(Cd_notumble_3D(alt_i_sel,ac,:)) - ...
               squeeze(Cd_bar_fixed(alt_i_sel,ac,:))) ./ ...
               squeeze(Cd_notumble_3D(alt_i_sel,ac,:)) * 100;
    plot(phi_vec, penalty, '-o', 'Color', cmap_accom(ac,:), ...
         'LineWidth', 1.2, 'MarkerSize', 5, ...
         'DisplayName', sprintf('$\\alpha=%.2f$', accom_vals(ac)))
end
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
ylabel('Tumbling Penalty (\%)','Interpreter','latex')
title('$h = 450$ km','Interpreter','latex')
xlim([43 87]); xticks(phi_vec)
legend('Interpreter','latex','Location','eastoutside','FontSize',7)

%% =========================================================================
%  C2: Cd_SAM / Cd_notumble ratio colourmap
%  Fraction of ideal drag recovered under tumbling
% =========================================================================

ratio_3D = Cd_bar_fixed ./ Cd_notumble_3D;

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)
    raw_data  = squeeze(ratio_3D(alt_i,:,:))';
    fine_data = griddata(ac_flat, phi_flat, raw_data(:), AC_fine, PHI_fine, 'cubic');
    contourf(accom_fine, phi_fine, fine_data, 40, 'LineColor','none')
    colormap(gca, parula)
    cb = colorbar;
    cb.Label.String      = '$\bar{C}_{D,\mathrm{SAM}} / C_{D,0}$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize    = 11;
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% =========================================================================
%  C3: Scatter plot Cd_SAM vs Cd_random for all (phi, alpha) combinations
%  One panel per altitude — deviation from diagonal = attitude dynamics matter
% =========================================================================

figure('Position',[100 100 1200 380]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)
    hold on; box on; grid on;

    sam_flat  = squeeze(Cd_bar_fixed(alt_i,:,:));
    rand_flat = squeeze(Cd_random_3D(alt_i,:,:));

    scatter(rand_flat(:), sam_flat(:), 30, cmap_accom(repmat(1:n_accom,n_phi,1)',:), 'filled')

    lims = [min([sam_flat(:);rand_flat(:)])*0.98, max([sam_flat(:);rand_flat(:)])*1.02];
    plot(lims, lims, 'k--', 'LineWidth', 1)   % diagonal
    xlim(lims); ylim(lims)
    xlabel('$\bar{C}_{D,\mathrm{rand}}$','Interpreter','latex')
    ylabel('$\bar{C}_{D,\mathrm{SAM}}$','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    axis square
end

%% =========================================================================
%  C4: Optimal phi as function of alpha for each altitude
%  argmax of Cd_SAM along phi dimension
% =========================================================================

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
for alt_i = 1:n_alts
    opt_phi = zeros(1, n_accom);
    for ac = 1:n_accom
        [~, idx] = max(squeeze(Cd_bar_fixed(alt_i, ac, :)));
        opt_phi(ac) = phi_vec(idx);
    end
    plot(accom_vals, opt_phi, '-o', 'Color', alt_colors(alt_i,:), ...
         'LineWidth', 1.4, 'MarkerSize', 6, ...
         'DisplayName', sprintf('$h=%d$ km', alt_vec(alt_i)))
end
xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
ylabel('Optimal $\phi$ (deg)','Interpreter','latex')
yticks(phi_vec); ylim([43 87])
legend('Interpreter','latex','Location','best')
grid on

%% =========================================================================
%  C5: Cd_notumble vs phi — theoretical maximum drag curve
%  All alpha, all altitudes
% =========================================================================

figure('Position',[100 100 1400 450]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i); hold on; box on; grid on;
    for ac = 1:n_accom
        plot(phi_vec, squeeze(Cd_notumble_3D(alt_i,ac,:)), '-o', ...
             'Color', cmap_accom(ac,:), 'LineWidth', 1.2, 'MarkerSize', 4)
    end
    xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    ylabel('$C_D(\theta=0^\circ)$','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([43 87]); xticks(phi_vec)
    if alt_i == n_alts
        legend(arrayfun(@(a) sprintf('$\\alpha=%.2f$',a), accom_vals, ...
               'UniformOutput',false),'Interpreter','latex', ...
               'Location','northwestoutside','FontSize',6)
    end
end

%% =========================================================================
%  C6: Feasibility map — binary pass/fail for Cd_SAM > threshold
%  Threshold = 1.5 (adjust as needed)
% =========================================================================

Cd_threshold = 1.5;

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    pass_fail = double(squeeze(Cd_bar_fixed(alt_i,:,:)) > Cd_threshold)';
    pass_fine = griddata(ac_flat, phi_flat, pass_fail(:), AC_fine, PHI_fine, 'linear');
    pass_fine = round(pass_fine);

    imagesc(accom_fine, phi_fine, pass_fine)
    axis xy
    colormap(gca, [0.85 0.25 0.10; 0.10 0.72 0.40])
    cb = colorbar;
    cb.Ticks = [0.25 0.75];
    cb.TickLabels = {'Fail','Pass'};
    cb.Label.String = sprintf('$\\bar{{C}}_D > %.1f$', Cd_threshold);
    cb.Label.Interpreter = 'latex';

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals); xtickangle(45)
end

%% =========================================================================
%  C7: 3D surface plot of Cd_SAM(phi, alpha) for each altitude
% =========================================================================

[AC_grid, PHI_grid] = meshgrid(accom_vals, phi_vec);

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)
    surf(AC_grid, PHI_grid, squeeze(Cd_bar_fixed(alt_i,:,:))', ...
         'EdgeColor','none','FaceAlpha',0.9)
    colormap(gca, turbo)
    xlabel('$\alpha$','Interpreter','latex')
    ylabel('$\phi$ (deg)','Interpreter','latex')
    zlabel('$\bar{C}_{D,\mathrm{SAM}}$','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    view(35, 30)
    cb = colorbar;
    cb.Label.String = '$\bar{C}_{D,\mathrm{SAM}}$';
    cb.Label.Interpreter = 'latex';
end

%% =========================================================================
%  C8: Grouped box plot — distribution of Cd_SAM across alpha for each phi
%  One figure per altitude
% =========================================================================

figure('Position',[100 100 1200 380]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)
    data = squeeze(Cd_bar_fixed(alt_i,:,:));   % [n_accom x n_phi]
    boxplot(data, 'Labels', arrayfun(@(p) sprintf('%d',p), phi_vec, ...
            'UniformOutput',false))
    xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    ylabel('$\bar{C}_{D,\mathrm{SAM}}$','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    grid on
end

%% =========================================================================
%  C9: Rank order heatmap of Cd_SAM across all 297 cases
%  Rows = altitude, cols = phi, colour = mean rank across alpha
% =========================================================================

all_vals = Cd_bar_fixed(:);
[~, sort_idx] = sort(all_vals);
ranks = zeros(size(all_vals));
ranks(sort_idx) = 1:numel(all_vals);
ranks_3D = reshape(ranks, n_alts, n_accom, n_phi);
mean_rank = squeeze(mean(ranks_3D, 2));   % [n_alts x n_phi]

figure('Position',[100 100 700 300]);
imagesc(phi_vec, 1:n_alts, mean_rank)
axis xy
colormap(turbo)
cb = colorbar;
cb.Label.String = 'Mean rank (higher = better $\bar{C}_D$)';
cb.Label.Interpreter = 'latex';
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
yticks(1:n_alts)
yticklabels(arrayfun(@(h) sprintf('%d km',h), alt_vec,'UniformOutput',false))
xticks(phi_vec)

%% =========================================================================
%  C10: Tornado sensitivity chart
%  Variance in Cd_SAM attributed to phi, alpha, altitude
% =========================================================================

var_phi  = var(mean(mean(Cd_bar_fixed, 1), 2), 0, 3);   % variance along phi
var_accom = var(mean(mean(Cd_bar_fixed, 1), 3), 0, 2);  % variance along accom
var_alt  = var(mean(mean(Cd_bar_fixed, 2), 3), 0, 1);   % variance along alt

contributions = [var_phi, var_accom, var_alt];
contributions = contributions / sum(contributions) * 100;
labels = {'$\phi$','$\alpha$','$h$'};

figure('Position',[100 100 500 300]);
barh(contributions, 'FaceColor', [0.20 0.55 0.85])
yticklabels(labels)
set(gca,'TickLabelInterpreter','latex')
xlabel('Percentage contribution to variance in $\bar{C}_D$ (\%)','Interpreter','latex')
grid on; box on
xlim([0 100])

%% =========================================================================
%  C11: Pareto front — Cd_SAM vs tumbling penalty
%  One point per (phi, alpha, alt), coloured by altitude
% =========================================================================

figure('Position',[100 100 600 450]);
hold on; box on; grid on;
for alt_i = 1:n_alts
    sam_v  = squeeze(Cd_bar_fixed(alt_i,:,:));
    pen_v  = squeeze(Cd_notumble_3D(alt_i,:,:)) - sam_v;
    scatter(sam_v(:), pen_v(:), 40, ...
            repmat(alt_colors(alt_i,:), numel(sam_v), 1), 'filled', ...
            'DisplayName', sprintf('$h=%d$ km', alt_vec(alt_i)))
end
xlabel('$\bar{C}_{D,\mathrm{SAM}}$','Interpreter','latex')
ylabel('Tumbling drag penalty $C_{D,0} - \bar{C}_{D,\mathrm{SAM}}$','Interpreter','latex')
legend('Interpreter','latex','Location','northwest')

%% =========================================================================
%  C12: Cd_SAM sensitivity — error bar style across altitudes
%  For each phi, show mean and std of Cd_SAM across alpha
% =========================================================================

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
offsets = [-0.3, 0, 0.3];
for alt_i = 1:n_alts
    mu_cd  = mean(squeeze(Cd_bar_fixed(alt_i,:,:)), 1);   % mean across alpha
    std_cd = std(squeeze(Cd_bar_fixed(alt_i,:,:)), 0, 1); % std across alpha
    errorbar(phi_vec + offsets(alt_i), mu_cd, std_cd, 'o', ...
             'Color', alt_colors(alt_i,:), 'LineWidth', 1.2, ...
             'MarkerSize', 5, 'CapSize', 5, ...
             'DisplayName', sprintf('$h=%d$ km', alt_vec(alt_i)))
end
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
ylabel('$\bar{C}_{D,\mathrm{SAM}}$ (mean $\pm$ std across $\alpha$)','Interpreter','latex')
xlim([43 87]); xticks(phi_vec)
legend('Interpreter','latex','Location','northwest')

%% =========================================================================
%  C13: Alpha threshold contour
%  For each phi and altitude, find alpha where Cd_SAM = 2.0
% =========================================================================

target_Cd = 2.0;

figure('Position',[100 100 600 420]);
hold on; box on; grid on;
for alt_i = 1:n_alts
    alpha_thresh = NaN(1, n_phi);
    for si = 1:n_phi
        cd_vec = squeeze(Cd_bar_fixed(alt_i, :, si));
        if max(cd_vec) >= target_Cd && min(cd_vec) <= target_Cd
            alpha_thresh(si) = interp1(cd_vec, accom_vals, target_Cd, 'pchip');
        end
    end
    plot(phi_vec, alpha_thresh, '-o', 'Color', alt_colors(alt_i,:), ...
         'LineWidth', 1.4, 'MarkerSize', 6, ...
         'DisplayName', sprintf('$h=%d$ km', alt_vec(alt_i)))
end
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
ylabel('$\alpha$ at which $\bar{C}_D = 2.0$','Interpreter','latex')
xlim([43 87]); xticks(phi_vec)
ylim([0.48 1.02])
legend('Interpreter','latex','Location','best')
yline(0.5,'k:'); yline(1.0,'k:')

fprintf('\nAll extended figures generated.\n')

