load('drag_sail_attitude_results.mat')

adbsat_base = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';
accom_vec   = {'0p50','0p55','0p60','0p65','0p70','0p75','0p80','0p85','0p90','0p95','1p00'};

% Compute uniform tumbling Cd_bar (simple mean of Cd over all AoA)
% Cd_uniform(alt_i, ac, si) = mean(Cd) over -179.5:179.5
Cd_uniform = zeros(size(Cd_bar_3D));

for alt_i = 1:3
    alt     = alt_vec(alt_i);
    alt_dir = fullfile(adbsat_base, sprintf('%dkm', alt));
    for ac = 1:length(accom_vec)
        for si = 1:length(phi_vec)
            fname = fullfile(alt_dir, sprintf('%ddeg_CLL_accom_%s.mat', ...
                    phi_vec(si), accom_vec{ac}));
            raw   = load(fname);
            Cd    = -raw.aedb.aero.Cf_wX(:);
            Cd(isnan(Cd)) = 0;
            Cd_uniform(alt_i, ac, si) = mean(Cd);
        end
    end
    fprintf('Done altitude %dkm\n', alt);
end

% Percentage difference: (SAM - uniform) / uniform * 100
pct_diff = (Cd_bar_3D - Cd_uniform) ./ Cd_uniform * 100;

%% Setup
fmt_ax = @(ax) set(ax,'Box','on','LineWidth',1.5,'FontSize',11, ...
    'TickLabelInterpreter','latex','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.3,'MinorGridAlpha',0.1);

clrs_sail  = turbo(length(phi_vec));
clrs_accom = parula(length(accom_vec));
clrs_alt   = lines(3);
alt_labels = {'350 km','450 km','650 km'};
alt_names  = {'350km','450km','650km'};
phi_labels = arrayfun(@(x) sprintf('$\\phi=%d^\\circ$',x), phi_vec, 'UniformOutput',false);
accom_labels = arrayfun(@(x) sprintf('$\\alpha=%.2f$',x), accom_vals, 'UniformOutput',false);
ref_ac = find(accom_vals==0.75);
ref_si = find(phi_vec==75);

function cmap = redblue_cmap(n)
    half = floor(n/2);
    r1 = linspace(0.2,1,half)'; g1 = linspace(0.2,1,half)'; b1 = ones(half,1);
    r2 = ones(half,1); g2 = linspace(1,0.2,half)'; b2 = linspace(1,0.2,half)';
    cmap = [r1,g1,b1; r2,g2,b2];
end



%% Plot 1: SAM vs Uniform — Cd_bar vs phi, fixed accom=0.75
%  solid = SAM, dashed = uniform, one subplot per altitude
fig1 = figure('Position',[60 60 1200 400]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    plot(phi_vec, squeeze(Cd_bar_3D(alt_i,ref_ac,:)), '-o', ...
         'Color',[0.2 0.4 0.8],'LineWidth',1.8,'MarkerSize',6, ...
         'DisplayName','SAM');
    hold on;
    plot(phi_vec, squeeze(Cd_uniform(alt_i,ref_ac,:)), '--s', ...
         'Color',[0.8 0.2 0.2],'LineWidth',1.8,'MarkerSize',6, ...
         'DisplayName','Uniform tumbling');
    xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    ylabel('$\bar{C}_D$','Interpreter','latex');
    title(alt_labels{alt_i});
    if alt_i==3
        lg=legend('Location','northwest'); set(lg,'Interpreter','latex');
    end
    fmt_ax(gca);
end
savefig(fig1,'plot_SAM_vs_uniform_Cd_phi.fig');
exportgraphics(fig1,'plot_SAM_vs_uniform_Cd_phi.png','Resolution',300);
fprintf('Saved plot 1\n');

%% Plot 2: % difference vs phi — one line per altitude, fixed accom=0.75
fig2 = figure('Position',[60 500 700 450]);
for alt_i = 1:3
    plot(phi_vec, squeeze(pct_diff(alt_i,ref_ac,:)), '-o', ...
         'Color',clrs_alt(alt_i,:),'LineWidth',1.5,'MarkerSize',6, ...
         'DisplayName',alt_labels{alt_i});
    hold on;
end
yline(0,'k--','HandleVisibility','off');
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
ylabel('$(\bar{C}_{D,\mathrm{SAM}} - \bar{C}_{D,\mathrm{uniform}}) / \bar{C}_{D,\mathrm{uniform}}$ (\%)','Interpreter','latex');
lg=legend('Location','best'); set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig2,'plot_pctdiff_vs_phi.fig');
exportgraphics(fig2,'plot_pctdiff_vs_phi.png','Resolution',300);
fprintf('Saved plot 2\n');

%% Plot 3: % difference vs accom — one line per sail angle, fixed altitude=450km
fig3 = figure('Position',[800 60 700 450]);
for si = 1:length(phi_vec)
    plot(accom_vals, squeeze(pct_diff(2,:,si)), '-o', ...
         'Color',clrs_sail(si,:),'LineWidth',1.5,'MarkerSize',5, ...
         'DisplayName',phi_labels{si});
    hold on;
end
yline(0,'k--','HandleVisibility','off');
xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
ylabel('$(\bar{C}_{D,\mathrm{SAM}} - \bar{C}_{D,\mathrm{uniform}}) / \bar{C}_{D,\mathrm{uniform}}$ (\%)','Interpreter','latex');
lg=legend('Location','best','NumColumns',2); set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig3,'plot_pctdiff_vs_accom_450km.fig');
exportgraphics(fig3,'plot_pctdiff_vs_accom_450km.png','Resolution',300);
fprintf('Saved plot 3\n');

%% Plot 4: Heatmap of % difference — phi vs accom, per altitude
function smooth_surf(x,y,Z,clim_vals)
    xi=linspace(min(x),max(x),200); yi=linspace(min(y),max(y),200);
    [Xi,Yi]=meshgrid(xi,yi); [X,Y]=meshgrid(x,y);
    Zi=interp2(X,Y,Z',Xi,Yi,'spline');
    surf(Xi,Yi,zeros(size(Zi)),Zi,'EdgeColor','none'); view(2);
    xlim([min(x) max(x)]); ylim([min(y) max(y)]);
    if nargin>3; clim(clim_vals); end
end

pct_lim = [min(pct_diff(:)) max(pct_diff(:))];
fig4 = figure('Position',[60 60 1400 420]);
for alt_i = 1:3
    subplot(1,3,alt_i);
    smooth_surf(accom_vals, phi_vec, squeeze(pct_diff(alt_i,:,:)), pct_lim);
    colormap(redblue_cmap(256));   % diverging colourmap centred at 0
    cb = colorbar;
    cb.Label.String = '$(\bar{C}_{D,\mathrm{SAM}} - \bar{C}_{D,\mathrm{uniform}}) / \bar{C}_{D,\mathrm{uniform}}$ (\%)';
    cb.Label.Interpreter = 'latex'; cb.Label.FontSize=10;
    cb.TickLabelInterpreter = 'latex';
    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11,'TickLabelInterpreter','latex');
end
savefig(fig4,'plot_pctdiff_heatmap.fig');
exportgraphics(fig4,'plot_pctdiff_heatmap.png','Resolution',300);
fprintf('Saved plot 4\n');

%% Plot 5: SAM vs Uniform — absolute values side by side bar chart
%  Fixed phi=75deg, accom=0.75, all altitudes
fig5 = figure('Position',[60 60 600 450]);
x = categorical({'350 km','450 km','650 km'});
x = reordercats(x,{'350 km','450 km','650 km'});
SAM_vals     = squeeze(Cd_bar_3D(:,ref_ac,ref_si));
uniform_vals = squeeze(Cd_uniform(:,ref_ac,ref_si));
b = bar(x, [SAM_vals, uniform_vals], 0.6);
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.8 0.2 0.2];
ylabel('$\bar{C}_D$','Interpreter','latex');
lg = legend({'SAM','Uniform tumbling'},'Location','best');
set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig5,'plot_SAM_vs_uniform_bar.fig');
exportgraphics(fig5,'plot_SAM_vs_uniform_bar.png','Resolution',300);
fprintf('Saved plot 5\n');

