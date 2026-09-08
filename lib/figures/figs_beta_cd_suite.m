load('drag_sail_attitude_results.mat')

B_3D = 1 ./ Cd_bar_3D;

%% Helpers
fmt = @(ax) set(ax, 'Box','on', 'LineWidth',1.5, ...
                'GridLineStyle','-', 'GridAlpha',0.3, ...
                'MinorGridLineStyle','-', 'MinorGridAlpha',0.1, ...
                'XMinorGrid','on', 'YMinorGrid','on', ...
                'FontSize', 11, 'TickLabelInterpreter','latex');

function save_fig(fig, name)
    savefig(fig, [name '.fig']);
    exportgraphics(fig, [name '.png'], 'Resolution', 300);
    fprintf('Saved: %s\n', name);
end

function cb = styled_colorbar(label_str)
    cb = colorbar;
    cb.Label.String         = label_str;
    cb.Label.Interpreter    = 'latex';
    cb.Label.FontSize       = 12;
    cb.TickLabelInterpreter = 'latex';
end

function smooth_surf(x, y, Z, clim_vals)
% Upsample data 10x for smooth colour transitions
    xi = linspace(min(x), max(x), 200);
    yi = linspace(min(y), max(y), 200);
    [Xi, Yi] = meshgrid(xi, yi);
    [X,  Y ] = meshgrid(x, y);
    Zi = interp2(X, Y, Z', Xi, Yi, 'spline');
    surf(Xi, Yi, zeros(size(Zi)), Zi, 'EdgeColor','none');
    view(2);
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    if nargin > 3
        clim(clim_vals);
    end
end

function smooth_contour(x, y, Z, n_levels, clim_vals)
% Upsample data 10x for smooth contour transitions
    xi = linspace(min(x), max(x), 200);
    yi = linspace(min(y), max(y), 200);
    [Xi, Yi] = meshgrid(xi, yi);
    [X,  Y ] = meshgrid(x, y);
    Zi = interp2(X, Y, Z', Xi, Yi, 'spline');
    contourf(xi, yi, Zi, n_levels, 'LineColor','none');
    if nargin > 4
        clim(clim_vals);
    end
end

clrs_sail  = turbo(length(phi_vec));
clrs_accom = parula(length(accom_vals));
clrs_alt   = lines(length(alt_vec));

alt_labels   = {'350 km', '450 km', '650 km'};
alt_names    = {'350km','450km','650km'};
phi_labels   = arrayfun(@(x) sprintf('$\\phi = %d^\\circ$', x), phi_vec, 'UniformOutput', false);
accom_labels = arrayfun(@(x) sprintf('$\\alpha = %.2f$', x), accom_vals, 'UniformOutput', false);
ref_ac = find(accom_vals == 0.75);
ref_si = find(phi_vec == 75);

Cd_lim = [min(Cd_bar_3D(:)) max(Cd_bar_3D(:))];
B_lim  = [min(B_3D(:))      max(B_3D(:))];

%% =========================================================
%  PLOT 1: Cd_bar vs phi per accom per altitude
% =========================================================
fig1 = figure('Position',[60 60 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for ac = 1:length(accom_vals)
        plot(phi_vec, squeeze(Cd_bar_3D(alt_i,ac,:)), '-o', ...
             'Color',clrs_accom(ac,:),'LineWidth',1.5,'MarkerSize',5, ...
             'DisplayName',accom_labels{ac});
        hold on;
    end
    xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    ylabel('$\bar{C}_D$','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i == 3
        lg = legend('Location','northwest','NumColumns',1,'FontSize',8);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig1,'plot01_Cd_vs_phi_per_accom');

%% =========================================================
%  PLOT 2: Cd_bar vs accom per sail per altitude
% =========================================================
fig2 = figure('Position',[60 500 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for si = 1:length(phi_vec)
        plot(accom_vals, squeeze(Cd_bar_3D(alt_i,:,si)), '-o', ...
             'Color',clrs_sail(si,:),'LineWidth',1.5,'MarkerSize',5, ...
             'DisplayName',phi_labels{si});
        hold on;
    end
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('$\bar{C}_D$','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i == 3
        lg = legend('Location','northeast','NumColumns',1,'FontSize',8);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig2,'plot02_Cd_vs_accom_per_sail');

%% =========================================================
%  PLOT 3: Cd_bar vs altitude per sail (accom=0.75)
% =========================================================
fig3 = figure('Position',[60 960 700 450]);
for si = 1:length(phi_vec)
    plot(alt_vec, squeeze(Cd_bar_3D(:,ref_ac,si)), '-o', ...
         'Color',clrs_sail(si,:),'LineWidth',1.5,'MarkerSize',6, ...
         'DisplayName',phi_labels{si});
    hold on;
end
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2);
set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig3,'plot03_Cd_vs_altitude_per_sail_accom075');

%% =========================================================
%  PLOT 4: Cd_bar vs altitude per accom (phi=75deg)
% =========================================================
fig4 = figure('Position',[800 60 700 450]);
for ac = 1:length(accom_vals)
    plot(alt_vec, squeeze(Cd_bar_3D(:,ac,ref_si)), '-o', ...
         'Color',clrs_accom(ac,:),'LineWidth',1.5,'MarkerSize',6, ...
         'DisplayName',accom_labels{ac});
    hold on;
end
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('$\bar{C}_D$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2);
set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig4,'plot04_Cd_vs_altitude_per_accom_phi75');

%% =========================================================
%  PLOT 5: Smooth heatmap Cd — phi vs accom per altitude
% =========================================================
fig5 = figure('Position',[60 60 1400 420]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    smooth_surf(accom_vals, phi_vec, squeeze(Cd_bar_3D(alt_i,:,:)), Cd_lim);
    colormap(parula);
    styled_colorbar('$\bar{C}_D$');
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
end
save_fig(fig5,'plot05_heatmap_Cd_phi_vs_accom');

%% =========================================================
%  PLOT 6: Smooth 3D surface Cd per altitude
% =========================================================
[ACCOM, PHI] = meshgrid(accom_vals, phi_vec);
accom_fine   = linspace(min(accom_vals), max(accom_vals), 100);
phi_fine     = linspace(min(phi_vec), max(phi_vec), 100);
[Af, Pf]     = meshgrid(accom_fine, phi_fine);

for alt_i = 1:3
    fig = figure('Position',[100+alt_i*30 100+alt_i*30 700 500]);
    Zf = interp2(ACCOM, PHI, squeeze(Cd_bar_3D(alt_i,:,:))', Af, Pf, 'spline');
    surf(Af, Pf, Zf, 'EdgeColor','none');
    colormap(parula);
    styled_colorbar('$\bar{C}_D$');
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    zlabel('$\bar{C}_D$','Interpreter','latex');
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
    view(45,30);
    save_fig(fig, sprintf('plot06_3Dsurf_Cd_%s',alt_names{alt_i}));
end

%% =========================================================
%  PLOT 7: Smooth 3D surface Cd — phi vs altitude (accom=0.75)
% =========================================================
[ALT_GRID, PHI_GRID] = meshgrid(alt_vec, phi_vec);
alt_fine  = linspace(min(alt_vec), max(alt_vec), 100);
[Alf, Pf2] = meshgrid(alt_fine, phi_fine);
Zf2 = interp2(ALT_GRID, PHI_GRID, squeeze(Cd_bar_3D(:,ref_ac,:))', Alf, Pf2, 'spline');

fig7 = figure('Position',[200 200 700 500]);
surf(Alf, Pf2, Zf2, 'EdgeColor','none');
colormap(parula);
styled_colorbar('$\bar{C}_D$');
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
zlabel('$\bar{C}_D$','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
view(45,30);
save_fig(fig7,'plot07_3Dsurf_Cd_phi_vs_altitude_accom075');

%% =========================================================
%  PLOT 8: Smooth contour Cd — phi vs accom per altitude
% =========================================================
fig8 = figure('Position',[60 60 1400 420]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    smooth_contour(accom_vals, phi_vec, squeeze(Cd_bar_3D(alt_i,:,:)), 20, Cd_lim);
    colormap(parula);
    styled_colorbar('$\bar{C}_D$');
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
end
save_fig(fig8,'plot08_contour_Cd_phi_vs_accom');

%% =========================================================
%  PLOT 9: Smooth contour Cd — phi vs altitude (accom=0.75)
% =========================================================
fig9 = figure('Position',[800 500 700 450]);
smooth_contour(alt_vec, phi_vec, squeeze(Cd_bar_3D(:,ref_ac,:)), 20);
colormap(parula);
styled_colorbar('$\bar{C}_D$');
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
save_fig(fig9,'plot09_contour_Cd_phi_vs_altitude_accom075');

%% =========================================================
%  PLOT 10: Sensitivity Cd to accom
% =========================================================
fig10 = figure('Position',[60 500 700 450]);
for alt_i = 1:3
    sens = zeros(1,length(phi_vec));
    for si = 1:length(phi_vec)
        cd_r = squeeze(Cd_bar_3D(alt_i,:,si));
        sens(si) = 100*(max(cd_r)-min(cd_r))/min(cd_r);
    end
    plot(phi_vec, sens, '-o','Color',clrs_alt(alt_i,:), ...
         'LineWidth',1.5,'MarkerSize',6,'DisplayName',alt_labels{alt_i});
    hold on;
end
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
ylabel('$(\Delta\bar{C}_D\,/\,\bar{C}_{D,\min})$ (\%)','Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig10,'plot10_sensitivity_Cd_accom');

%% =========================================================
%  PLOT 11: Sensitivity Cd to geometry
% =========================================================
fig11 = figure('Position',[800 500 700 450]);
for alt_i = 1:3
    geo_sens = zeros(1,length(accom_vals));
    for ac = 1:length(accom_vals)
        cd_r = squeeze(Cd_bar_3D(alt_i,ac,:));
        geo_sens(ac) = 100*(max(cd_r)-min(cd_r))/min(cd_r);
    end
    plot(accom_vals, geo_sens, '-o','Color',clrs_alt(alt_i,:), ...
         'LineWidth',1.5,'MarkerSize',6,'DisplayName',alt_labels{alt_i});
    hold on;
end
xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
ylabel('$(\Delta\bar{C}_D\,/\,\bar{C}_{D,\min})$ (\%)','Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig11,'plot11_sensitivity_Cd_geometry');

%% =========================================================
%  PLOT 12: SAM PDFs
% =========================================================
fig12 = figure('Position',[60 60 1400 800]);
for alt_i = 1:3
    for si = 1:length(phi_vec)
        subplot(3,length(phi_vec),(alt_i-1)*length(phi_vec)+si);
        bar(ref_results(alt_i,si).centres, ref_results(alt_i,si).pdf_SAM, 1, ...
            'FaceColor',clrs_sail(si,:),'EdgeColor','none','FaceAlpha',0.8);
        xlim([-20 20]);
        ttl = sprintf('%s, $\\phi=%d^\\circ$', alt_labels{alt_i}, phi_vec(si));
        title(ttl,'Interpreter','latex','FontSize',7);
        if si==1; ylabel('$p(\theta)$ (deg$^{-1}$)','Interpreter','latex'); end
        if alt_i==3; xlabel('$\theta$ (deg)','Interpreter','latex'); end
        set(gca,'Box','on','LineWidth',1.0,'FontSize',8,'TickLabelInterpreter','latex');
    end
end
save_fig(fig12,'plot12_SAM_PDFs_all');

%% =========================================================
%  PLOT 13: Pitch transient per sail (accom=0.75)
% =========================================================
fig13 = figure('Position',[60 60 1200 800]);
for alt_i = 1:3
    subplot(3,1,alt_i);
    for si = 1:length(phi_vec)
        plot(ref_results(alt_i,si).t/60, ref_results(alt_i,si).phi_deg, ...
             'Color',clrs_sail(si,:),'LineWidth',1.0,'DisplayName',phi_labels{si});
        hold on;
    end
    yline(0,'k--','HandleVisibility','off');
    ylabel('$\theta$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i==3; xlabel('Time (min)','Interpreter','latex'); end
    if alt_i==1
        lg = legend('NumColumns',3,'FontSize',8,'Location','northeast');
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig13,'plot13_pitch_transient_per_sail');

%% =========================================================
%  PLOT 14: Pitch transient per altitude — fixed phi=75deg
%  (uses ref_results which is accom=0.75 only)
% =========================================================
fig14 = figure('Position',[60 60 700 450]);
for alt_i = 1:3
    plot(ref_results(alt_i,ref_si).t/60, ref_results(alt_i,ref_si).phi_deg, ...
         'Color',clrs_alt(alt_i,:),'LineWidth',1.5,'DisplayName',alt_labels{alt_i});
    hold on;
end
yline(0,'k--','HandleVisibility','off');
xlabel('Time (min)','Interpreter','latex');
ylabel('$\theta$ (deg)','Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig14,'plot14_pitch_transient_all_altitudes_phi75');

%% =========================================================
%  PLOT 15: Zoomed pitch transient — first 5 min
% =========================================================
sel_sails  = [1 5 9];
sail_deg   = [45 65 85];
fig15 = figure('Position',[60 60 1200 800]);
sp = 0;
for si_k = 1:3
    si = sel_sails(si_k);
    for alt_i = 1:3
        sp = sp+1;
        subplot(3,3,sp);
        t  = ref_results(alt_i,si).t/60;
        ph = ref_results(alt_i,si).phi_deg;
        plot(t(t<=5), ph(t<=5), 'Color',clrs_alt(alt_i,:),'LineWidth',1.5);
        yline(0,'k--');
        xlabel('Time (min)','Interpreter','latex');
        ylabel('$\theta$ (deg)','Interpreter','latex');
        yline(0,'k--');
        ttl = sprintf('%s, $\\phi = %d^\\circ$', alt_labels{alt_i}, sail_deg(si_k));
        title(ttl,'Interpreter','latex','FontSize',9);
        fmt(gca);
    end
end
save_fig(fig15,'plot15_pitch_transient_zoomed_5min');

%% =========================================================
%  BALLISTIC COEFFICIENT PLOTS (beta_bar = 1/Cd_bar)
% =========================================================

% B1: beta vs phi per accom per altitude
fig_B1 = figure('Position',[60 60 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for ac = 1:length(accom_vals)
        plot(phi_vec, squeeze(B_3D(alt_i,ac,:)), '-o', ...
             'Color',clrs_accom(ac,:),'LineWidth',1.5,'MarkerSize',5, ...
             'DisplayName',accom_labels{ac});
        hold on;
    end
    xlabel('$\theta$ (deg)', 'Interpreter','latex');
    ylabel('$p(\theta)$ (deg$^{-1}$)', 'Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i==3
        lg = legend('Location','northeast','NumColumns',1,'FontSize',8);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig_B1,'plot_B01_beta_vs_phi_per_accom');

% B2: beta vs accom per sail per altitude
fig_B2 = figure('Position',[60 500 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for si = 1:length(phi_vec)
        plot(accom_vals, squeeze(B_3D(alt_i,:,si)), '-o', ...
             'Color',clrs_sail(si,:),'LineWidth',1.5,'MarkerSize',5, ...
             'DisplayName',phi_labels{si});
        hold on;
    end
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('$\bar{\beta}$','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i==3
        lg = legend('Location','northwest','NumColumns',1,'FontSize',8);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig_B2,'plot_B02_beta_vs_accom_per_sail');

% B3: beta vs altitude per sail (accom=0.75)
fig_B3 = figure('Position',[60 960 700 450]);
for si = 1:length(phi_vec)
    plot(alt_vec, squeeze(B_3D(:,ref_ac,si)), '-o', ...
         'Color',clrs_sail(si,:),'LineWidth',1.5,'MarkerSize',6, ...
         'DisplayName',phi_labels{si});
    hold on;
end
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('$\bar{\beta}$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2);
set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig_B3,'plot_B03_beta_vs_altitude_per_sail_accom075');

% B4: beta vs altitude per accom (phi=75deg)
fig_B4 = figure('Position',[800 60 700 450]);
for ac = 1:length(accom_vals)
    plot(alt_vec, squeeze(B_3D(:,ac,ref_si)), '-o', ...
         'Color',clrs_accom(ac,:),'LineWidth',1.5,'MarkerSize',6, ...
         'DisplayName',accom_labels{ac});
    hold on;
end
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('$\bar{\beta}$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2);
set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig_B4,'plot_B04_beta_vs_altitude_per_accom_phi75');

% B5: Smooth heatmap beta per altitude
fig_B5 = figure('Position',[60 60 1400 420]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    smooth_surf(accom_vals, phi_vec, squeeze(B_3D(alt_i,:,:)), B_lim);
    colormap(parula);
    styled_colorbar('$\bar{\beta}$');
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
end
save_fig(fig_B5,'plot_B05_heatmap_beta_phi_vs_accom');

% B6: Smooth 3D surface beta per altitude
for alt_i = 1:3
    fig = figure('Position',[100+alt_i*30 100+alt_i*30 700 500]);
    Zf_B = interp2(ACCOM, PHI, squeeze(B_3D(alt_i,:,:))', Af, Pf, 'spline');
    surf(Af, Pf, Zf_B, 'EdgeColor','none');
    colormap(parula);
    styled_colorbar('$\bar{\beta}$');
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    zlabel('$\bar{\beta}$','Interpreter','latex');
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
    view(45,30);
    save_fig(fig, sprintf('plot_B06_3Dsurf_beta_%s',alt_names{alt_i}));
end

% B7: Smooth 3D surface beta — phi vs altitude (accom=0.75)
fig_B7 = figure('Position',[200 200 700 500]);
Zf_B2 = interp2(ALT_GRID, PHI_GRID, squeeze(B_3D(:,ref_ac,:))', Alf, Pf2, 'spline');
surf(Alf, Pf2, Zf_B2, 'EdgeColor','none');
colormap(parula);
styled_colorbar('$\bar{\beta}$');
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
zlabel('$\bar{\beta}$','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
view(45,30);
save_fig(fig_B7,'plot_B07_3Dsurf_beta_phi_vs_altitude_accom075');

% B8: Smooth contour beta per altitude
fig_B8 = figure('Position',[60 60 1400 420]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    smooth_contour(accom_vals, phi_vec, squeeze(B_3D(alt_i,:,:)), 20, B_lim);
    colormap(parula);
    styled_colorbar('$\bar{\beta}$');
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
end
save_fig(fig_B8,'plot_B08_contour_beta_phi_vs_accom');

% B9: Smooth contour beta phi vs altitude (accom=0.75)
fig_B9 = figure('Position',[800 500 700 450]);
smooth_contour(alt_vec, phi_vec, squeeze(B_3D(:,ref_ac,:)), 20);
colormap(parula);
styled_colorbar('$\bar{\beta}$');
xlabel('Orbital Altitude (km)','Interpreter','latex');
ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
save_fig(fig_B9,'plot_B09_contour_beta_phi_vs_altitude_accom075');

% B10: Sensitivity beta to accom
fig_B10 = figure('Position',[60 500 700 450]);
for alt_i = 1:3
    sens = zeros(1,length(phi_vec));
    for si = 1:length(phi_vec)
        b_r = squeeze(B_3D(alt_i,:,si));
        sens(si) = 100*(max(b_r)-min(b_r))/min(b_r);
    end
    plot(phi_vec, sens, '-o','Color',clrs_alt(alt_i,:), ...
         'LineWidth',1.5,'MarkerSize',6,'DisplayName',alt_labels{alt_i});
    hold on;
end
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
ylabel('$(\Delta\bar{\beta}\,/\,\bar{\beta}_{\min})$ (\%)','Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig_B10,'plot_B10_sensitivity_beta_accom');

% B11: Sensitivity beta to geometry
fig_B11 = figure('Position',[800 500 700 450]);
for alt_i = 1:3
    geo_sens = zeros(1,length(accom_vals));
    for ac = 1:length(accom_vals)
        b_r = squeeze(B_3D(alt_i,ac,:));
        geo_sens(ac) = 100*(max(b_r)-min(b_r))/min(b_r);
    end
    plot(accom_vals, geo_sens, '-o','Color',clrs_alt(alt_i,:), ...
         'LineWidth',1.5,'MarkerSize',6,'DisplayName',alt_labels{alt_i});
    hold on;
end
xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
ylabel('$(\Delta\bar{\beta}\,/\,\bar{\beta}_{\min})$ (\%)','Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig_B11,'plot_B11_sensitivity_beta_geometry');

