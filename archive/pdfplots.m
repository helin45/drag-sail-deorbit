load('drag_sail_attitude_results.mat')

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

clrs_alt  = lines(3);
clrs_sail = turbo(length(phi_vec));
alt_labels = {'350 km', '450 km', '650 km'};
phi_labels = arrayfun(@(x) sprintf('$\\phi = %d^\\circ$', x), phi_vec, 'UniformOutput', false);
ref_ac = find(accom_vals == 0.75);

%% =========================================================
%  PDF PLOT 1: Scatter PDF — all altitudes, one subplot per sail angle
%  Fixed accom = 0.75
% =========================================================
fig_P1 = figure('Position',[60 60 1400 900]);
for si = 1:length(phi_vec)
    subplot(3,3,si);
    for alt_i = 1:3
        centres = ref_results(alt_i,si).centres;
        pdf     = ref_results(alt_i,si).pdf_SAM;
        scatter(centres, pdf, 18, 'filled', ...
                'MarkerFaceColor', clrs_alt(alt_i,:), ...
                'MarkerFaceAlpha', 0.7, ...
                'DisplayName', alt_labels{alt_i});
        hold on;
    end
    xlabel('$\phi$ (deg)', 'Interpreter','latex');
    ylabel('$p(\phi)$ (deg$^{-1}$)', 'Interpreter','latex');
    ttl = sprintf('$\\phi = %d^\\circ$', phi_vec(si));
    title(ttl, 'Interpreter','latex');
    if si == length(phi_vec)
        lg = legend('Location','best','FontSize',8);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig_P1, 'plot_PDF01_scatter_all_altitudes_per_sail');

%% =========================================================
%  PDF PLOT 2: Scatter PDF zoomed -10 to 10 deg, log scale
%  Fixed accom = 0.75, one subplot per sail angle
% =========================================================
fig_P2 = figure('Position',[60 60 1400 900]);
for si = 1:length(phi_vec)
    subplot(3,3,si);
    for alt_i = 1:3
        centres = ref_results(alt_i,si).centres;
        pdf     = ref_results(alt_i,si).pdf_SAM;
        mask    = centres >= -10 & centres <= 10 & pdf > 0;
        scatter(centres(mask), pdf(mask), 25, 'filled', ...
                'MarkerFaceColor', clrs_alt(alt_i,:), ...
                'MarkerFaceAlpha', 0.8, ...
                'DisplayName', alt_labels{alt_i});
        hold on;
    end
    set(gca, 'YScale','log');
    xlim([-10 10]);
    xlabel('$\phi$ (deg)', 'Interpreter','latex');
    ylabel('$p(\phi)$ (deg$^{-1}$)', 'Interpreter','latex');
    ttl = sprintf('$\\phi = %d^\\circ$', phi_vec(si));
    title(ttl, 'Interpreter','latex');
    if si == length(phi_vec)
        lg = legend('Location','best','FontSize',8);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig_P2, 'plot_PDF02_scatter_zoomed_log_per_sail');

%% =========================================================
%  PDF PLOT 3: Scatter PDF — all sail angles, one subplot per altitude
%  Full range, fixed accom = 0.75
% =========================================================
fig_P3 = figure('Position',[60 60 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for si = 1:length(phi_vec)
        centres = ref_results(alt_i,si).centres;
        pdf     = ref_results(alt_i,si).pdf_SAM;
        scatter(centres, pdf, 12, 'filled', ...
                'MarkerFaceColor', clrs_sail(si,:), ...
                'MarkerFaceAlpha', 0.6, ...
                'DisplayName', phi_labels{si});
        hold on;
    end
    xlabel('$\phi$ (deg)', 'Interpreter','latex');
    ylabel('$p(\phi)$ (deg$^{-1}$)', 'Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i == 3
        lg = legend('Location','best','FontSize',7,'NumColumns',1);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig_P3, 'plot_PDF03_scatter_all_sails_per_altitude');

%% =========================================================
%  PDF PLOT 4: Scatter PDF zoomed -10 to 10 deg, log scale
%  All sail angles, one subplot per altitude
% =========================================================
fig_P4 = figure('Position',[60 500 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    for si = 1:length(phi_vec)
        centres = ref_results(alt_i,si).centres;
        pdf     = ref_results(alt_i,si).pdf_SAM;
        mask    = centres >= -10 & centres <= 10 & pdf > 0;
        scatter(centres(mask), pdf(mask), 20, 'filled', ...
                'MarkerFaceColor', clrs_sail(si,:), ...
                'MarkerFaceAlpha', 0.8, ...
                'DisplayName', phi_labels{si});
        hold on;
    end
    set(gca, 'YScale','log');
    xlim([-10 10]);
    xlabel('$\phi$ (deg)', 'Interpreter','latex');
    ylabel('$p(\phi)$ (deg$^{-1}$)', 'Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i == 3
        lg = legend('Location','best','FontSize',7,'NumColumns',1);
        set(lg,'Interpreter','latex');
    end
    fmt(gca);
end
save_fig(fig_P4, 'plot_PDF04_scatter_zoomed_log_all_sails_per_altitude');

%% =========================================================
%  PDF PLOT 5: Overlaid scatter — all altitudes + all sails
%  Full range, coloured by altitude, one big plot
% =========================================================
fig_P5 = figure('Position',[60 60 900 500]);
for alt_i = 1:3
    for si = 1:length(phi_vec)
        centres = ref_results(alt_i,si).centres;
        pdf     = ref_results(alt_i,si).pdf_SAM;
        dn = sprintf('_dummy_%d_%d', alt_i, si);
        if si == 1
            sc = scatter(centres, pdf, 8, 'filled', ...
                    'MarkerFaceColor', clrs_alt(alt_i,:), ...
                    'MarkerFaceAlpha', 0.4, ...
                    'DisplayName', alt_labels{alt_i});
        else
            scatter(centres, pdf, 8, 'filled', ...
                    'MarkerFaceColor', clrs_alt(alt_i,:), ...
                    'MarkerFaceAlpha', 0.4, ...
                    'HandleVisibility','off');
        end
        hold on;
    end
end
xlabel('$\phi$ (deg)', 'Interpreter','latex');
ylabel('$p(\phi)$ (deg$^{-1}$)', 'Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig_P5, 'plot_PDF05_scatter_all_combined');

%% =========================================================
%  PDF PLOT 6: Zoomed -10 to 10, log scale
%  All altitudes + all sails overlaid, coloured by altitude
% =========================================================
fig_P6 = figure('Position',[60 560 900 500]);
for alt_i = 1:3
    for si = 1:length(phi_vec)
        centres = ref_results(alt_i,si).centres;
        pdf     = ref_results(alt_i,si).pdf_SAM;
        mask    = centres >= -10 & centres <= 10 & pdf > 0;
        if si == 1
            scatter(centres(mask), pdf(mask), 18, 'filled', ...
                    'MarkerFaceColor', clrs_alt(alt_i,:), ...
                    'MarkerFaceAlpha', 0.6, ...
                    'DisplayName', alt_labels{alt_i});
        else
            scatter(centres(mask), pdf(mask), 18, 'filled', ...
                    'MarkerFaceColor', clrs_alt(alt_i,:), ...
                    'MarkerFaceAlpha', 0.6, ...
                    'HandleVisibility','off');
        end
        hold on;
    end
end
set(gca,'YScale','log');
xlim([-10 10]);
xlabel('$\phi$ (deg)', 'Interpreter','latex');
ylabel('$p(\phi)$ (deg$^{-1}$)', 'Interpreter','latex');
lg = legend('Location','best'); set(lg,'Interpreter','latex');
fmt(gca);
save_fig(fig_P6, 'plot_PDF06_scatter_zoomed_log_all_combined');