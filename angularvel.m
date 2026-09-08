%% Zoomed angular velocity — first 5 min, selected geometries, all altitudes
sel_sails  = [1 5 9];
sail_deg   = [45 65 85];

fig3 = figure('Position',[60 60 1200 800]);
sp = 0;
for si_k = 1:3
    si = sel_sails(si_k);
    for alt_i = 1:3
        sp = sp+1;
        subplot(3,3,sp);
        t     = ref_results(alt_i,si).t/60;
        theta = deg2rad(ref_results(alt_i,si).phi_deg);
        dtheta_dt_deg = rad2deg(gradient(theta, ref_results(alt_i,si).t));
        mask  = t <= 5;
        plot(t(mask), dtheta_dt_deg(mask), 'Color',clrs_alt(alt_i,:), 'LineWidth',1.5);
        yline(0,'k--');
        xlabel('Time (min)','Interpreter','latex');
        ylabel('$\dot{\theta}$ (deg s$^{-1}$)','Interpreter','latex');
        ttl = sprintf('%s, $\\phi = %d^\\circ$', alt_labels{alt_i}, sail_deg(si_k));
        title(ttl,'Interpreter','latex','FontSize',9);
        set(gca,'Box','on','LineWidth',1.5,'FontSize',9,'TickLabelInterpreter','latex', ...
            'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.3,'MinorGridAlpha',0.1);
    end
end
savefig(fig3,'plot_angular_velocity_zoomed_5min.fig');
exportgraphics(fig3,'plot_angular_velocity_zoomed_5min.png','Resolution',300);
fprintf('Saved zoomed angular velocity plot\n');