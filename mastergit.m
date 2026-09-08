function results = run_drag_sail_sweep(config)

%% ===================== DEFAULT CONFIG =====================
if nargin == 0
    config = struct();
end

config = apply_defaults(config);

%% ===================== CONSTANTS =====================
mu = 3.986004418e14;
Re = 6.3781e6;

load(config.ref_geometry_file, 'ref_geom');
load(config.inertia_file, 'inertiaDB');

%% ===================== PREALLOC =====================
n_alts  = numel(config.alt_vec);
n_accom = numel(config.accom_vec);
n_sails = numel(config.phi_vec);

Cd_bar_3D = zeros(n_alts, n_accom, n_sails);
ref_accom_idx = find(strcmp(config.accom_vec, config.ref_accom));
ref_results = struct();

%% ===================== MAIN SWEEP =====================
for alt_i = 1:n_alts
    alt = config.alt_vec(alt_i);
    alt_dir = fullfile(config.adbsat_base, sprintf('%dkm', alt));

    for ac = 1:n_accom
        accom_str = config.accom_vec{ac};

        for si = 1:n_sails
            phi = config.phi_vec(si);

            fname = fullfile(alt_dir, ...
                sprintf('%ddeg_CLL_accom_%s.mat', phi, accom_str));

            assert(isfile(fname), 'Missing file: %s', fname);

            raw  = load(fname);
            aero = raw.aedb.aero;

            Cm = fillmissing(aero.Cm_BY(:), 'constant', 0);
            Cd = fillmissing(-aero.Cf_wX(:), 'constant', 0);

            I_idx = find(inertiaDB.angles_deg == phi);
            g_idx = find(ref_geom.phi_deg == phi);

            sails = struct( ...
                'alpha', (-179.5:1:179.5)', ...
                'Cm_pitch', Cm, ...
                'Cd', Cd, ...
                'Iyy', inertiaDB.I(2,2,I_idx), ...
                'A_ref', ref_geom.A_ref(g_idx), ...
                'L_ref', ref_geom.L_ref(g_idx) ...
            );

            %% ORBIT
            a_orb = Re + alt*1e3;
            T_orb = 2*pi * sqrt(a_orb^3 / mu);
            v_circ = sqrt(mu / a_orb);

            t_vec = linspace(0, T_orb, config.n_steps)';

            %% PARAM STRUCT
            p = build_param_struct(config, a_orb, T_orb, v_circ, sails);

            %% ODE
            opts = odeset('RelTol', config.ode_reltol, ...
                          'AbsTol', config.ode_abstol);

            [t_out, X_out] = ode45(@(t,X) attitude_ode(t, X, p), ...
                                   t_vec, [config.phi0; config.dphi0], opts);

            phi_deg_t = rad2deg(X_out(:,1));

            %% SAM PDF
            [pdf_SAM, centres, dt_k] = compute_sam_pdf( ...
                t_out, phi_deg_t, config.d_theta_deg);

            %% Cd_bar
            Cd_bar = compute_cd_bar(sails, centres, pdf_SAM, config.d_theta_deg);
            Cd_bar_3D(alt_i, ac, si) = Cd_bar;

            %% STORE REFERENCE
            if ac == ref_accom_idx
                ref_results(alt_i, si) = struct( ...
                    'alt_km', alt, ...
                    'phi_deg', phi, ...
                    'Cd_bar', Cd_bar, ...
                    't', t_out, ...
                    'phi_t', phi_deg_t, ...
                    'pdf', pdf_SAM ...
                );
            end
        end
    end
end

%% OUTPUT
results.Cd_bar_3D = Cd_bar_3D;
results.ref_results = ref_results;
results.config = config;

save(config.output_file, '-struct', 'results');

end

function config = apply_defaults(config)

defaults = struct( ...
    'adbsat_base', './adbsat_processed', ...
    'ref_geometry_file', 'ref_geometry.mat', ...
    'inertia_file', 'inertia_tensors.mat', ...
    'output_file', 'results.mat', ...
    'alt_vec', [350 450 650], ...
    'accom_vec', {'0p50','0p75','1p00'}, ...
    'phi_vec', 45:5:85, ...
    'ref_accom', '0p75', ...
    'n_steps', 3000, ...
    'ode_reltol', 1e-8, ...
    'ode_abstol', 1e-10, ...
    'phi0', deg2rad(5), ...
    'dphi0', 0, ...
    'd_theta_deg', 2, ...
    'year', 2025, ...
    'month', 6, ...
    'day', 21, ...
    'hour', 10, ...
    'f107', 150, ...
    'f107a', 150, ...
    'ap', 15 ...
);

fields = fieldnames(defaults);
for i = 1:numel(fields)
    if ~isfield(config, fields{i})
        config.(fields{i}) = defaults.(fields{i});
    end
end

end