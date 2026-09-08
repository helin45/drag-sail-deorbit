load('drag_sail_attitude_results.mat')
adbsat_base = '/Users/helintaha/Library/CloudStorage/OneDrive-TheUniversityofManchester/Dissertation/MASTERCODE/adbsat_processed';
clrs        = turbo(length(phi_vec));
clrs_accom  = parula(length(accom_vals));
accom_vec   = {'0p50','0p55','0p60','0p65','0p70','0p75','0p80','0p85','0p90','0p95','1p00'};

fmt_ax = @(ax) set(ax,'Box','on','LineWidth',1.5,'FontSize',11, ...
    'TickLabelInterpreter','latex','XMinorGrid','on','YMinorGrid','on', ...
    'GridAlpha',0.3,'MinorGridAlpha',0.1);

%% ADM01: Cd vs theta per sail angle (accom=0.75, 450km)
fig = figure('Position',[60 60 700 450]);
for si = 1:length(phi_vec)
    fname = fullfile(adbsat_base,'450km',sprintf('%ddeg_CLL_accom_0p75.mat',phi_vec(si)));
    raw   = load(fname);
    aoa   = raw.aedb.aoa;
    Cd    = -raw.aedb.aero.Cf_wX;
    plot(rad2deg(aoa), Cd, 'Color',clrs(si,:), 'LineWidth',1.5, ...
         'DisplayName',sprintf('$\\phi = %d^\\circ$',phi_vec(si)));
    hold on;
end
xlim([-90 90]);
xlabel('Pitch Angle $\theta$ (deg)','Interpreter','latex');
ylabel('$C_D$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2); set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig,'plot_ADM01_Cd_vs_theta.fig');
exportgraphics(fig,'plot_ADM01_Cd_vs_theta.png','Resolution',300);
fprintf('Saved ADM01\n');

%% ADM02: Cm vs theta per sail angle (accom=0.75, 450km)
fig = figure('Position',[60 60 700 450]);
for si = 1:length(phi_vec)
    fname = fullfile(adbsat_base,'450km',sprintf('%ddeg_CLL_accom_0p75.mat',phi_vec(si)));
    raw   = load(fname);
    aoa   = raw.aedb.aoa;
    Cm    = raw.aedb.aero.Cm_BY;
    plot(rad2deg(aoa), Cm, 'Color',clrs(si,:), 'LineWidth',1.5, ...
         'DisplayName',sprintf('$\\phi = %d^\\circ$',phi_vec(si)));
    hold on;
end
xline(0,'k--','HandleVisibility','off');
yline(0,'k--','HandleVisibility','off');
xlim([-90 90]);
xlabel('Pitch Angle $\theta$ (deg)','Interpreter','latex');
ylabel('$C_m$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2); set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig,'plot_ADM02_Cm_vs_theta.fig');
exportgraphics(fig,'plot_ADM02_Cm_vs_theta.png','Resolution',300);
fprintf('Saved ADM02\n');

%% ADM03: Cd at theta=0 vs phi per accom (450km)
fig = figure('Position',[60 60 700 450]);
for ac = 1:length(accom_vals)
    Cd_at_zero = zeros(1,length(phi_vec));
    for si = 1:length(phi_vec)
        fname = fullfile(adbsat_base,'450km', ...
                sprintf('%ddeg_CLL_accom_%s.mat',phi_vec(si),accom_vec{ac}));
        raw   = load(fname);
        aoa   = raw.aedb.aoa;
        Cd    = -raw.aedb.aero.Cf_wX;
        [~,idx] = min(abs(rad2deg(aoa)));
        Cd_at_zero(si) = Cd(idx);
    end
    plot(phi_vec, Cd_at_zero, '-o', 'Color',clrs_accom(ac,:), 'LineWidth',1.5, ...
         'MarkerSize',5,'DisplayName',sprintf('$\\alpha=%.2f$',accom_vals(ac)));
    hold on;
end
xlim([-90 90]);
xlabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
ylabel('$C_D(\theta=0)$','Interpreter','latex');
lg = legend('Location','best','NumColumns',2); set(lg,'Interpreter','latex');
fmt_ax(gca);
savefig(fig,'plot_ADM03_Cd_at_zero_vs_phi.fig');
exportgraphics(fig,'plot_ADM03_Cd_at_zero_vs_phi.png','Resolution',300);
fprintf('Saved ADM03\n');