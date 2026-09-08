load('decay_times.mat')
load('sam_results.mat')

alt_labels = {'350 km','450 km','650 km'};
T_caps     = [2, 10, 220];

fig = figure('Position',[60 60 1400 450]);
for alt_i = 1:3
    subplot(1,3,alt_i);

    T_plot = squeeze(T_decay(alt_i,:,:));

    xi = linspace(min(accom_vals),max(accom_vals),300);
    yi = linspace(min(phi_vec),max(phi_vec),300);
    [Xi,Yi] = meshgrid(xi,yi);
    [X,Y]   = meshgrid(accom_vals,phi_vec);
    Zi      = interp2(X,Y,T_plot',Xi,Yi,'spline');
    Zi      = max(Zi,0);

    contourf(Xi,Yi,Zi,50,'LineColor','none');
    uistack(h_grey, 'top');
    colormap(gca,parula);
    cb = colorbar;
    cb.Label.String         = 'Estimated Decay Time (years)';
    cb.Label.Interpreter    = 'latex';
    cb.Label.FontSize       = 10;
    cb.TickLabelInterpreter = 'latex';
    clim([0 T_caps(alt_i)]);
    hold on;

    % Grey shading where T >= 5
    if max(Zi(:)) >= 5
        grey = double(Zi >= 5);
        grey(grey == 0) = NaN;
        h_grey = pcolor(Xi,Yi,grey);
        h_grey.FaceColor = [0.55 0.55 0.55];
        h_grey.FaceAlpha = 0.45;
        h_grey.EdgeColor = 'none';
    end

    % Red dashed 5-year contour only where it crosses
    if max(Zi(:)) >= 5 && min(Zi(:)) <= 5
        contour(Xi,Yi,Zi,[5 5],'r--','LineWidth',2.5);
    end

    xlabel('Accommodation Coefficient $\alpha$','Interpreter','latex');
    ylabel('Apex Half-Angle $\phi$ (deg)','Interpreter','latex');
    title(alt_labels{alt_i});
    set(gca,'Box','on','LineWidth',1.5,'FontSize',11, ...
        'TickLabelInterpreter','latex','Layer','top');
end

savefig(fig,'plot_decay_time_integrated.fig');
exportgraphics(fig,'plot_decay_time_integrated.png','Resolution',300);
fprintf('Saved\n');