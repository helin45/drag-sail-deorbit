%% =========================================================================
%  ATTITUDE SAM SWEEP — ALL PLOTS (original + 7 new)
%  Uses Cd_bar_fixed | data-driven ylims | no figure titles
% =========================================================================
clear; clc; close all;

load('attitude_SAM_sweep_FINAL.mat');

%% NOTE for Plot 5 (decay time colourmap):
%  Load your decay lifetime results here. Expected variable:
%  decay_years_3D [n_alts x n_accom x n_phi]
%  e.g. load('decay_results.mat','decay_years_3D');

alt_vec    = [350, 450, 650];
phi_vec    = 45:5:85;
accom_vals = 0.50:0.05:1.00;
accom_vec  = {'0p50','0p55','0p60','0p65','0p70','0p75', ...
              '0p80','0p85','0p90','0p95','1p00'};

n_alts  = length(alt_vec);
n_phi   = length(phi_vec);
n_accom = length(accom_vals);

ref_ac_idx  = find(strcmp(accom_vec,'0p75'));
ref_phi_idx = find(phi_vec == 85);   % reference phi for PDF comparison plots

cmap_accom = parula(n_accom);
cmap_phi   = turbo(n_phi);
alt_colors = [0.15 0.45 0.85;   % 350 km — blue
              0.10 0.72 0.40;   % 450 km — green
              0.85 0.25 0.10];  % 650 km — red

% Global Cd ylim
cd_min  = floor(min(Cd_bar_fixed(:)) * 20) / 20;
cd_max  = ceil( max(Cd_bar_fixed(:)) * 20) / 20;
cd_ylim = [cd_min - 0.05, cd_max + 0.05];

mu = 3.986004418e14;
Re = 6.3781e6;

%% =========================================================================
%  PRE-COMPUTE: Cd_random and Cd_notumble from adbsat files
%  (needed for plots 6 and 7)
% =========================================================================

adbsat_base = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';

Cd_random_3D   = zeros(n_alts, n_accom, n_phi);
Cd_notumble_3D = zeros(n_alts, n_accom, n_phi);
alpha_vec_full = (-179.5:1:179.5)';

fprintf('Pre-computing Cd_random and Cd_notumble...\n');
for alt_i = 1:n_alts
    for ac = 1:n_accom
        for si = 1:n_phi
            fname = fullfile(adbsat_base, sprintf('%dkm', alt_vec(alt_i)), ...
                    sprintf('%ddeg_CLL_accom_%s.mat', phi_vec(si), accom_vec{ac}));
            raw    = load(fname);
            Cd_raw = -raw.aedb.aero.Cf_wX(:);
            Cd_raw(isnan(Cd_raw)) = 0;
            Cd_random_3D(alt_i, ac, si)   = mean(Cd_raw);
            Cd_notumble_3D(alt_i, ac, si) = interp1(alpha_vec_full, Cd_raw, 0, 'pchip');
        end
    end
end
fprintf('Done.\n');

%% =========================================================================
%  PRE-COMPUTE: NRLMSISE density along orbit (for plot A1)
% =========================================================================

year_ep = 2025; month_ep = 1; day_ep = 15; hour_utc = 0;
f107a = 150; f107 = 150; ap = 15;
n_orbit_pts = 360;

density_data(n_alts) = struct('t',[],'rho',[],'T_orb',[],'nu_deg',[]);

fprintf('Computing NRLMSISE density profiles...\n');
for alt_i = 1:n_alts
    alt   = alt_vec(alt_i);
    a_orb = Re + alt*1e3;
    T_orb = 2*pi*sqrt(a_orb^3/mu);
    t_vec = linspace(0, T_orb, n_orbit_pts);
    rho_vec = zeros(1, n_orbit_pts);

    for j = 1:n_orbit_pts
        t      = t_vec(j);
        nu     = mod(2*pi/T_orb * t, 2*pi);
        r_eci  = a_orb * [cos(nu); sin(nu); 0];

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

    density_data(alt_i).t      = t_vec;
    density_data(alt_i).rho    = rho_vec;
    density_data(alt_i).T_orb  = T_orb;
    density_data(alt_i).nu_deg = rad2deg(linspace(0,2*pi,n_orbit_pts));
end
fprintf('Done.\n');

%% =========================================================================
%  FIGURE 1: Cd_bar vs phi — panel title = altitude
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
               'Location','northwest','FontSize',7)
    end
end

%% =========================================================================
%  FIGURE 2: Cd_bar vs accom — panel title = altitude
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
               'Location','northeast','FontSize',7)
    end
end

%% =========================================================================
%  FIGURE 3: Pitch transient — coloured by ALTITUDE
% =========================================================================

figure('Position',[100 100 700 750]);
for alt_i = 1:n_alts
    subplot(n_alts,1,alt_i); hold on; box on; grid on;
    all_theta_alt = [];
    for si = 1:n_phi
        t         = ref_results(alt_i, ref_ac_idx, si).t;
        phi_deg_t = wrapTo180(ref_results(alt_i, ref_ac_idx, si).phi_deg);
        plot(t/60, phi_deg_t, 'Color', alt_colors(alt_i,:), 'LineWidth', 0.8)
        all_theta_alt = [all_theta_alt; phi_deg_t];
    end
    yline(0,'k--','LineWidth',0.5)
    xlabel('Time (min)','Interpreter','latex')
    ylabel('$\theta$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([0 ref_results(alt_i,ref_ac_idx,1).T_orb_min])
    th_min = floor(min(all_theta_alt)*2)/2 - 1;
    th_max = ceil( max(all_theta_alt)*2)/2 + 1;
    ylim([th_min th_max])
    if alt_i == 1
        legend(arrayfun(@(p) sprintf('$\\phi=%d^\\circ$',p), phi_vec, ...
               'UniformOutput',false),'Interpreter','latex', ...
               'NumColumns',3,'Location','northeast','FontSize',6)
    end
end

%% =========================================================================
%  FIGURE 4: 5-min zoom — panel title = altitude + phi
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
        th_min = floor(min(theta5)*2)/2 - 0.5;
        th_max = ceil( max(theta5)*2)/2 + 0.5;
        ylim([th_min th_max])
        sp = sp+1;
    end
end

%% =========================================================================
%  FIGURE 5: SAM PDF phi=85, alpha=0.75 — panel title = altitude
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
%  FIGURE 6: Altitude sensitivity — panel title = alpha value
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
%  NEW FIGURE A1: NRLMSISE density along orbit — panel title = altitude
% =========================================================================

figure('Position',[100 100 1200 380]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)
    plot(density_data(alt_i).nu_deg, density_data(alt_i).rho, ...
         'Color', alt_colors(alt_i,:), 'LineWidth', 1.5)
    xlabel('True Anomaly (deg)','Interpreter','latex')
    ylabel('$\rho$ (kg m$^{-3}$)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([0 360]); xticks(0:90:360)
    grid on; box on
end

%% =========================================================================
%  NEW FIGURE A2: PDFs for all phi, all altitudes
%  Colour varies by phi (cmap_phi), one panel per altitude
% =========================================================================

figure('Position',[100 100 1200 380]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i); hold on; box on; grid on;
    for si = 1:n_phi
        pdf_vals = ref_results(alt_i, ref_ac_idx, si).pdf_SAM;
        centres  = ref_results(alt_i, ref_ac_idx, si).centres;
        plot(centres, pdf_vals, 'Color', cmap_phi(si,:), 'LineWidth', 1.2)
    end
    xlabel('$\theta$ (deg)','Interpreter','latex')
    ylabel('PDF (deg$^{-1}$)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    xlim([-180 180]); xticks(-180:60:180); grid on
    if alt_i == n_alts
        legend(arrayfun(@(p) sprintf('$\\phi=%d^\\circ$',p), phi_vec, ...
               'UniformOutput',false),'Interpreter','latex', ...
               'Location','north','FontSize',6,'NumColumns',3)
    end
end

%% =========================================================================
%  NEW FIGURE A3: PDFs for 3 accom coeffs, all altitudes
%  Fixed phi=85, colour by altitude, one panel per accom
% =========================================================================

figure('Position',[100 100 1200 380]);
for k = 1:length(accom_sel)
    subplot(1,3,k); hold on; box on; grid on;
    for alt_i = 1:n_alts
        ac_idx   = accom_sel_idx(k);
        pdf_vals = ref_results(alt_i, ac_idx, ref_phi_idx).pdf_SAM;
        centres  = ref_results(alt_i, ac_idx, ref_phi_idx).centres;
        plot(centres, pdf_vals, 'Color', alt_colors(alt_i,:), 'LineWidth', 1.4)
    end
    xlabel('$\theta$ (deg)','Interpreter','latex')
    ylabel('PDF (deg$^{-1}$)','Interpreter','latex')
    title(sprintf('$\\alpha = %.2f$', accom_sel(k)),'Interpreter','latex')
    xlim([-180 180]); xticks(-180:60:180); grid on
    if k == length(accom_sel)
        legend(arrayfun(@(h) sprintf('$h=%d$ km',h), alt_vec, ...
               'UniformOutput',false),'Interpreter','latex','Location','north')
    end
end

%% =========================================================================
%  NEW FIGURE A4: Colourmap — beta_bar = 1/Cd_bar (ballistic coeff proxy)
%  x = accom, y = phi, one panel per altitude
% =========================================================================

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    beta_bar = 1 ./ squeeze(Cd_bar_fixed(alt_i,:,:));  % [n_accom x n_phi]

    imagesc(accom_vals, phi_vec, beta_bar')
    axis xy; colormap(gca, turbo)
    cb = colorbar;
    cb.Label.String = '$\bar{\beta}$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals)
    xtickangle(45)
end

%% =========================================================================
%  NEW FIGURE A5: Colourmap — orbital decay time (years)
%  Load decay_years_3D [n_alts x n_accom x n_phi] before running
% =========================================================================

if exist('decay_years_3D','var')

    figure('Position',[100 100 1400 420]);
    for alt_i = 1:n_alts
        subplot(1,3,alt_i)

        decay_slice = squeeze(decay_years_3D(alt_i,:,:));  % [n_accom x n_phi]

        imagesc(accom_vals, phi_vec, decay_slice')
        axis xy; colormap(gca, turbo)
        cb = colorbar;
        cb.Label.String = 'Decay time (years)';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 11;

        xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
        ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
        title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
        yticks(phi_vec); xticks(accom_vals)
        xtickangle(45)
    end

else
    fprintf('NOTE: decay_years_3D not found — skipping Figure A5.\n')
    fprintf('      Load your decay results mat file and re-run.\n')
end

%% =========================================================================
%  NEW FIGURE A6: Bar chart — SAM vs random tumbling vs no tumbling
%  Fixed phi=75, alpha=0.75, one group of bars per altitude
% =========================================================================

phi_bar_idx = find(phi_vec == 75);
ac_bar_idx  = ref_ac_idx;   % alpha=0.75

Cd_SAM_bar      = zeros(1, n_alts);
Cd_random_bar   = zeros(1, n_alts);
Cd_notumble_bar = zeros(1, n_alts);

for alt_i = 1:n_alts
    Cd_SAM_bar(alt_i)      = Cd_bar_fixed(alt_i, ac_bar_idx, phi_bar_idx);
    Cd_random_bar(alt_i)   = Cd_random_3D(alt_i, ac_bar_idx, phi_bar_idx);
    Cd_notumble_bar(alt_i) = Cd_notumble_3D(alt_i, ac_bar_idx, phi_bar_idx);
end

figure('Position',[100 100 600 420]);
hold on; box on; grid on;

x      = 1:n_alts;
w      = 0.25;
bar(x - w,   Cd_SAM_bar,      w, 'FaceColor', [0.20 0.55 0.85], 'DisplayName','SAM')
bar(x,       Cd_random_bar,   w, 'FaceColor', [0.10 0.75 0.40], 'DisplayName','Random tumbling')
bar(x + w,   Cd_notumble_bar, w, 'FaceColor', [0.85 0.30 0.10], 'DisplayName','No tumbling ($\theta=0$)')

xticks(x)
xticklabels(arrayfun(@(h) sprintf('%d km',h), alt_vec, 'UniformOutput',false))
ylabel('$C_D$','Interpreter','latex')
legend('Interpreter','latex','Location','northeast')
ylim([0, max([Cd_SAM_bar, Cd_random_bar, Cd_notumble_bar])*1.15])

%% =========================================================================
%  NEW FIGURE A7: Colourmap — % difference SAM vs random tumbling
%  x = accom, y = phi, one panel per altitude
% =========================================================================

pct_diff_3D = 100 * (Cd_bar_fixed - Cd_random_3D) ./ Cd_random_3D;  % [n_alts x n_accom x n_phi]

pct_abs_max = ceil(max(abs(pct_diff_3D(:))));
clim_sym    = [-pct_abs_max, pct_abs_max];

figure('Position',[100 100 1400 420]);
for alt_i = 1:n_alts
    subplot(1,3,alt_i)

    pct_slice = squeeze(pct_diff_3D(alt_i,:,:));   % [n_accom x n_phi]

    imagesc(accom_vals, phi_vec, pct_slice')
    axis xy; colormap(gca, 'coolwarm')
    clim(clim_sym)
    cb = colorbar;
    cb.Label.String = '$\Delta \bar{C}_D$ (\%)';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 11;

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex')
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex')
    title(sprintf('$h = %d$ km', alt_vec(alt_i)),'Interpreter','latex')
    yticks(phi_vec); xticks(accom_vals)
    xtickangle(45)
end

fprintf('\nAll figures generated.\n')
