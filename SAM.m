%% =========================================================================
%  SAM.m  --  Statistical Averaging Model
%
%  Stage 3 of the pipeline (ADM -> ROADM -> SAM).
%
%  Loads the ROADM pitch-angle time histories (roadm_results.mat), and for
%  every (phi, alpha, h) case:
%    1. builds the residence-time probability density p(theta) by binning
%       the dwell time of theta(t) over one bin width d_theta;
%    2. uses p(theta) as weights in a Riemann sum over the ADM drag-
%       coefficient profile C_D(theta) to get the orientation-averaged
%       drag coefficient Cd_bar;
%    3. forms the orientation-averaged ballistic coefficient
%       beta_bar = m / (Cd_bar * A_ref).
%
%  Output:  sam_results.mat
%      Cd_bar_3D   [n_alt x n_accom x n_phi]
%      beta_bar_3D [n_alt x n_accom x n_phi]
%      pdf_all     struct [n_alt x n_accom x n_phi] with .centres, .pdf_SAM
%      ref_results struct [n_alt x n_phi] for the reference accommodation
%                  slice (time history + PDF), consumed by make_figures
%      plus alt_vec, accom_vec, accom_vals, phi_vec, inertiaDB, ref_geom
%
%  Run ROADM first.
% =========================================================================
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'lib')));
cd(here);

%% ----- settings ------------------------------------------------------
d_theta_deg = 2.0;        % residence-time bin width (deg)
M_SAIL      = 5.8;        % spacecraft + sail mass (kg) for beta_bar
REF_ACCOM   = '0p75';     % accommodation slice stored in ref_results

if ~isfile('roadm_results.mat')
    error('roadm_results.mat not found -- run ROADM first');
end
S = load('roadm_results.mat');
roadm      = S.roadm;
sail_db    = S.sail_db;
alt_vec    = S.alt_vec;
accom_vec  = S.accom_vec;
accom_vals = S.accom_vals;
phi_vec    = S.phi_vec;
inertiaDB  = S.inertiaDB;
ref_geom   = S.ref_geom;

n_alts  = numel(alt_vec);
n_accom = numel(accom_vec);
n_phi   = numel(phi_vec);
ref_ac  = find(strcmp(accom_vec, REF_ACCOM), 1);
if isempty(ref_ac); ref_ac = find(abs(accom_vals - 0.75) < 1e-9, 1); end

edges   = (-180 : d_theta_deg : 180);
centres = edges(1:end-1) + d_theta_deg/2;
n_bins  = numel(centres);

Cd_bar_3D   = zeros(n_alts, n_accom, n_phi);
beta_bar_3D = zeros(n_alts, n_accom, n_phi);
pdf_all     = struct('centres', {}, 'pdf_SAM', {});
pdf_all(n_alts, n_accom, n_phi).centres = [];
ref_results = struct();

%% ----- averaging ---------------------------------------------------
for alt_i = 1:n_alts
    for ac = 1:n_accom
        for si = 1:n_phi
            t_out     = roadm(alt_i,ac,si).t;
            phi_deg_t = roadm(alt_i,ac,si).theta_deg;

            % residence-time PDF (same binning as the former MASTERCODE.m)
            dt_k = diff(t_out); dt_k(end+1) = dt_k(end); %#ok<AGROW>
            Delta_t = zeros(1, n_bins);
            for k = 1:numel(t_out)
                bi = find(edges <= phi_deg_t(k), 1, 'last');
                if ~isempty(bi) && bi <= n_bins
                    Delta_t(bi) = Delta_t(bi) + dt_k(k);
                end
            end
            T_total = sum(dt_k);
            pdf_SAM = Delta_t / (T_total * d_theta_deg);

            % orientation-averaged drag coefficient
            alpha_grid = sail_db(ac,si).alpha_deg;
            Cd_curve   = sail_db(ac,si).Cd;
            Cd_bins = zeros(1, n_bins);
            for bi = 1:n_bins
                Cd_bins(bi) = interp1(alpha_grid, Cd_curve, centres(bi), ...
                                      'pchip', 'extrap');
            end
            Cd_bar = sum(Cd_bins .* pdf_SAM) * d_theta_deg;

            A_ref = sail_db(ac,si).A_ref;
            Cd_bar_3D(alt_i,ac,si)   = Cd_bar;
            beta_bar_3D(alt_i,ac,si) = M_SAIL / (Cd_bar * A_ref);

            pdf_all(alt_i,ac,si).centres = centres;
            pdf_all(alt_i,ac,si).pdf_SAM = pdf_SAM;

            if ac == ref_ac
                ref_results(alt_i,si).label        = sprintf('%dkm_%ddeg', ...
                                                     alt_vec(alt_i), phi_vec(si));
                ref_results(alt_i,si).alt_km       = alt_vec(alt_i);
                ref_results(alt_i,si).phi_deg_sail = phi_vec(si);
                ref_results(alt_i,si).Iyy          = roadm(alt_i,ac,si).Iyy;
                ref_results(alt_i,si).t            = t_out;
                ref_results(alt_i,si).phi_deg      = phi_deg_t;
                ref_results(alt_i,si).centres      = centres;
                ref_results(alt_i,si).pdf_SAM      = pdf_SAM;
                ref_results(alt_i,si).Cd_bar       = Cd_bar;
                ref_results(alt_i,si).T_orb_min    = roadm(alt_i,ac,si).T_orb_s / 60;
            end

            fprintf('alt=%dkm | accom=%s | phi=%2ddeg | Cd_bar=%.4f | beta_bar=%.3f\n', ...
                    alt_vec(alt_i), accom_vec{ac}, phi_vec(si), Cd_bar, ...
                    beta_bar_3D(alt_i,ac,si));
        end
    end
end

%% ----- save ------------------------------------------------------
save('sam_results.mat', 'Cd_bar_3D', 'beta_bar_3D', 'pdf_all', 'ref_results', ...
     'alt_vec', 'accom_vec', 'accom_vals', 'phi_vec', 'inertiaDB', 'ref_geom', ...
     'd_theta_deg', 'M_SAIL', 'REF_ACCOM');

fprintf('\nSaved sam_results.mat\n');
fprintf('Cd_bar_3D   size %s   range [%.3f, %.3f]\n', mat2str(size(Cd_bar_3D)), ...
        min(Cd_bar_3D(:)), max(Cd_bar_3D(:)));
fprintf('beta_bar_3D range [%.3f, %.3f] kg/m^2\n', ...
        min(beta_bar_3D(:)), max(beta_bar_3D(:)));
fprintf('Next: make_figures\n');
