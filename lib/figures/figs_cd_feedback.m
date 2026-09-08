load('drag_sail_attitude_results.mat')

adbsat_base = fullfile(project_root(), 'adbsat_processed');
alt = 450; accom = '0p75';
clrs_sail = turbo(length(phi_vec));
phi_labels = arrayfun(@(x) sprintf('$\\phi = %d^\\circ$',x), phi_vec,'UniformOutput',false);

fmt_ax = @(ax) set(ax,'Box','on','LineWidth',1.5,'FontSize',11,...
    'TickLabelInterpreter','latex','XMinorGrid','on','YMinorGrid','on',...
    'GridAlpha',0.3,'MinorGridAlpha',0.1);

%% Plot 1: Full -180 to 180 from ADBSat (actual values)
fig1 = figure('Position',[60 60 700 450]);
hold on;
for si = 1:length(phi_vec)
    fname = fullfile(adbsat_base, sprintf('%dkm',alt), ...
            sprintf('%ddeg_CLL_accom_%s.mat', phi_vec(si), accom));
    raw  = load(fname);
    Cd   = -raw.aedb.aero.Cf_wX(:);
   nans = isnan(Cd);
    Cd(nans) = interp1(find(~nans), Cd(~nans), find(nans), 'pchip');
    theta_vec = (-179.5:1:179.5)';
    plot(theta_vec, Cd, 'Color',clrs_sail(si,:), 'LineWidth',1.5,...
         'DisplayName',phi_labels{si});
end
xlabel('Pitch Angle $\theta$ (deg)','Interpreter','latex');
ylabel('$C_D$','Interpreter','latex');
xlim([-180 180]); xticks(-180:45:180);
lg = legend('NumColumns',2,'Location','north','FontSize',8);
set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig1,'plot_ADM01_Cd_vs_theta_full.fig');
exportgraphics(fig1,'plot_ADM01_Cd_vs_theta_full.png','Resolution',300);
fprintf('Saved plot 1 (full range)\n');

%% Plot 2: 0 to 180 reflected to give -180 to 0
fig2 = figure('Position',[60 500 700 450]);
hold on;
for si = 1:length(phi_vec)
    fname = fullfile(adbsat_base, sprintf('%dkm',alt), ...
            sprintf('%ddeg_CLL_accom_%s.mat', phi_vec(si), accom));
    raw  = load(fname);
    Cd   = -raw.aedb.aero.Cf_wX(:);
   nans = isnan(Cd);
    Cd(nans) = interp1(find(~nans), Cd(~nans), find(nans), 'pchip');
    theta_vec = (-179.5:1:179.5)';

    mask_pos  = theta_vec >= 0;
    theta_pos = theta_vec(mask_pos);
    Cd_pos    = Cd(mask_pos);

    theta_neg = -flipud(theta_pos);
    Cd_neg    = flipud(Cd_pos);

    theta_full = [theta_neg; theta_pos];
    Cd_full    = [Cd_neg;    Cd_pos];

    plot(theta_full, Cd_full, 'Color',clrs_sail(si,:), 'LineWidth',1.5,...
         'DisplayName',phi_labels{si});
end
xlabel('Pitch Angle $\theta$ (deg)','Interpreter','latex');
ylabel('$C_D$','Interpreter','latex');
xlim([-180 180]); xticks(-180:45:180);
lg = legend('NumColumns',2,'Location','north','FontSize',8);
set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig2,'plot_ADM01_Cd_vs_theta_reflected.fig');
exportgraphics(fig2,'plot_ADM01_Cd_vs_theta_reflected.png','Resolution',300);
fprintf('Saved plot 2 (reflected)\n');