clear; close all; clc;

% ============================================================
%  KATSURAGI EQUATION 2 — Granular Penetration Velocity Fit
%
%  Reference: H. Katsuragi & D. J. Durian (2013)
%             "Unified force law for granular impact cratering"
%             Nature Physics
%
%  Physical Problem:
%    A sphere of mass m impacts a granular bed at velocity v0.
%    As it penetrates to depth z, it decelerates under gravity
%    and two granular resistive forces:
%      (1) inertial drag   ~ k·v²
%      (2) depth-dependent ~ k·d1·(2v/d1)·z  (quasi-static term)
%
%  Equation of motion integrated over depth gives Katsuragi Eq.2:
%
%    v²(z) = v₀²·exp(−2z/d₁)
%            − (k·d₁/m)·z
%            + [g·d₁ + k·d₁²/(2m)]·[1 − exp(−2z/d₁)]
%
%  Free fit parameters:
%    d₁  [m]    — characteristic inertial penetration depth
%    k   [kg/m] — inertial drag coefficient
%
%  Known/fixed constants (all in SI):
%    v₀  [m/s]  — surface entry velocity (taken directly from data)
%    m   [kg]   — projectile mass (derived from initial LAMMPS force)
%    g   [m/s²] — 9.8 (matches LAMMPS "fix bed all gravity 9.8")
%
%  LAMMPS Simulation Details (from log file):
%    Projectile type 3: diameter = 0.20 m, density = 500 kg/m³
%    → m = (4/3)π(0.10)³ × 500 = 2.094 kg  (matches force check)
%    Initial velocity set: vz = −3.00 m/s (downward)
%    Total simulation steps with thermo: 62119 rows
% ============================================================


% ============================================================
%  SECTION 1 — Load LAMMPS log file
% ============================================================
% The LAMMPS thermo_style in this simulation is:
%   Step   Atoms   v_vmaxx   v_vZZ   v_FZZ
%
%   col 1 → timestep number
%   col 2 → atom count = 200001 (fixed, used only as column check)
%   col 3 → v_vmaxx : |vz| of projectile  [m/s]
%   col 4 → v_vZZ   : z-position of projectile centroid  [m]
%   col 5 → v_FZZ   : net z-force on projectile  [N]
%
% Strategy: read every line, parse floating-point numbers.
%   Rows with exactly 5 numbers are thermo data rows.
%   Header lines, section titles, warnings → skipped automatically.
%
% All units in the LAMMPS file are already SI (units si in input).

filename = 'log (3).lammps';

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: "%s"\nCheck your working directory.', filename);
end

data = [];
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        nums = sscanf(line, '%f');
        if length(nums) == 5          % only thermo output rows
            data = [data; nums'];     %#ok<AGROW>
        end
    end
end
fclose(fid);

% Assign physical columns
v_raw = data(:, 3);   % projectile speed |vz|  [m/s]
z_pos = data(:, 4);   % projectile z-position  [m]
f_raw = data(:, 5);   % net z-force             [N]

fprintf('\n  File loaded successfully.\n');
fprintf('  Total thermo rows parsed : %d\n', size(data, 1));


% ============================================================
%  SECTION 2 — Compute projectile mass from first LAMMPS frame
% ============================================================
% At timestep 0 the projectile sits above the granular bed in free
% flight with no granular contact.  The only force acting is gravity:
%
%   F_z(t=0) = −m·g   →   m = |F_z(0)| / g
%
% LAMMPS uses g = 9.8 m/s² (fix bed all gravity 9.8 vector 0 0 -1).
% Using 9.8 here (not 9.81) is necessary for internal consistency.
%
% Independent cross-check from LAMMPS input script:
%   Sphere diameter d = 0.20 m  →  radius r = 0.10 m
%   Density ρ = 500 kg/m³
%   Volume V = (4/3)π r³ = (4/3)π(0.10)³ = 4.189 × 10⁻³ m³
%   m = ρ·V = 500 × 4.189e-3 = 2.094 kg  ✓  (matches force-derived value)

g = 9.8;                     % [m/s²] — gravitational acceleration in LAMMPS
m = abs(f_raw(1)) / g;       % [kg]   — projectile mass

fprintf('\n  === Projectile Properties ===\n');
fprintf('  m (from initial force F = -m·g) : %.6f kg\n', m);
fprintf('  F(t=0) = %.4f N  →  m = %.4f/%.1f = %.4f kg\n', ...
        abs(f_raw(1)), abs(f_raw(1)), g, m);


% ============================================================
%  SECTION 3 — Convert time-series to spatial profile
% ============================================================
% Katsuragi Eq.2 gives v as a function of depth z, not time.
% Sorting the time-series data by z-position maps the trajectory
% onto the spatial coordinate required by the model.
%
% Depth origin: penetration depth z is measured from the entry point
% (surface), NOT from absolute z-coordinate in LAMMPS box.
% We define: z_pen = z_pos − min(z_pos)  after sorting.
%   → z_pen = 0 at the deepest point (minimum LAMMPS z-position)
%   → z_pen = max at the surface (where the projectile enters)
% The depth direction is corrected in Section 4 below.

[z_sorted, sort_idx] = sort(z_pos);           % ascending z-position sort
v_sorted             = v_raw(sort_idx);        % velocities in same order

z     = z_sorted - min(z_sorted);             % penetration depth [m]
v_use = v_sorted;

fprintf('\n  z-position range  : %.5f m  to  %.5f m\n', ...
        min(z_pos), max(z_pos));
fprintf('  Penetration depth : %.5f m  (= %.2f mm)\n', ...
        max(z), max(z)*1000);


% ============================================================
%  SECTION 4 — Velocity trend diagnostic and direction correction
% ============================================================
% Physical requirement (Katsuragi model):
%   v must DECREASE from v₀ (at z=0, surface) to ~0 (at z=zmax, deep).
%
% After ascending z-sort:
%   z(1)   = 0   corresponds to min(z_pos) = deepest → v ≈ 0
%   z(end) = max corresponds to max(z_pos) = surface → v ≈ v₀
%
% Therefore v_use(end) > v_use(1) → depth axis is inverted.
% Fix: flip z so that z=0 is at the surface (high-velocity end).
% This is purely a coordinate relabelling, not data manipulation.

fprintf('\n  === Velocity Trend Check ===\n');
fprintf('  v at z_raw = 0   (deepest)  : %.6e m/s\n', v_use(1));
fprintf('  v at z_raw = max (surface)  : %.6e m/s\n', v_use(end));

if v_use(end) > v_use(1)
    % Velocity increases with z_raw → surface is at the high-z end.
    % Reverse depth direction: new z=0 at surface (max z_pos).
    fprintf('  [FIX] Velocity increases with z. Reversing depth axis.\n');
    fprintf('        (Transforms z so surface is at z=0 — Katsuragi convention.)\n\n');

    z     = max(z) - z;     % mirror: high z_raw → z=0 (surface)
    z     = flipud(z);      % reorder vector elements to match
    v_use = flipud(v_use);  % reorder velocities consistently
else
    fprintf('  [OK] Velocity decreases with z_raw. No reversal needed.\n\n');
end


% ============================================================
%  SECTION 5 — Initial velocity v₀ (directly from LAMMPS data)
% ============================================================
% v₀ is taken as the FIRST element of v_use after direction correction.
% This is the LAMMPS-recorded projectile speed at the surface entry.
%
% Per physical requirement: NO averaging over multiple early points,
% NO smoothing window.  LAMMPS provides the exact initial condition at
% each timestep — using v_use(1) directly is the correct approach.
%
% Note: LAMMPS sets velocity exactly to 3.00 m/s at step 4533500.
%       The first few timesteps show a tiny increase (~3.0006 m/s)
%       due to free-fall above the granular surface before first contact.
%       This is physical and present in the raw data.

v0 = v_use(1);    % [m/s] — surface entry velocity, direct from data

fprintf('  === Fit Inputs ===\n');
fprintf('  v₀ = %.6e m/s  (v_use(1), direct from LAMMPS)\n', v0);
fprintf('  m  = %.6e kg\n', m);
fprintf('  g  = %.2f m/s²  (matches LAMMPS gravity fix)\n', g);
fprintf('  N  = %d  data points\n\n', length(z));


% ============================================================
%  SECTION 6 — Katsuragi Equation 2 (velocity form, raw v in m/s)
% ============================================================
% The original Eq.2 in Katsuragi (2013) is written as v²/v₀².
% Here we multiply through by v₀² to obtain v²(z), then take sqrt.
%
%   v²(z) = v₀²·exp(−2z/d₁)
%            − (k·d₁/m)·z
%            + [g·d₁ + k·d₁²/(2m)]·[1 − exp(−2z/d₁)]
%
%   v(z) = sqrt( max(v²(z), 0) )
%
%   Physical meaning of each term:
%     ① v₀²·exp(−2z/d₁)
%        Exponential decay from impact — pure inertial drag term.
%        d₁ sets the e-folding depth: smaller d₁ → faster deceleration.
%
%     ② −(k·d₁/m)·z
%        Linear depth term — the quasi-static weight of overburden
%        on the projectile.  Grows with depth, always decelerating.
%
%     ③ [g·d₁ + k·d₁²/(2m)]·[1 − exp(−2z/d₁)]
%        Correction term that ensures v→0 properly at large z,
%        incorporating both gravity and depth-dependent drag together.
%
%   The max(·, 0) guard prevents sqrt of negative numbers during
%   Levenberg-Marquardt iteration steps far from the solution.
%
% p(1) = d₁  [m]
% p(2) = k   [kg/m]

eq2_vel = @(p, z) sqrt(max( ...
    v0^2 .* exp(-2.*z ./ p(1)) ...
    - (p(2) .* p(1) ./ m) .* z ...
    + (g .* p(1) + p(2) .* p(1).^2 ./ (2.*m)) ...
      .* (1 - exp(-2.*z ./ p(1))), ...
    0));


% ============================================================
%  SECTION 7 — Fitting bounds (physically motivated)
% ============================================================
% d₁ bounds:
%   Lower (max(z)×1e-3): d₁ cannot be negligibly small — if d₁→0,
%     the exponential decays within the first micron, physically absurd.
%     max(z)×1e-3 = 0.1% of total depth is a conservative lower floor.
%   Upper (max(z)×10): d₁ larger than 10× total penetration means the
%     exponential barely decays over the data range (≈ linear regime).
%     Keeping this as UB allows the fit to reach quasi-linear behavior
%     if warranted by the data.
%
% k bounds:
%   Lower (0): drag coefficient is non-negative by definition.
%   Upper (1e15): effectively unconstrained — LAMMPS granular forces
%     can span many orders of magnitude depending on contact stiffness.
%
% Note: lb(1) uses the 1e-3 factor as specified; ub(1) = max(z)×10
% provides a physically wide but bounded search range for d₁.

lb = [max(z)*1e-3,   0    ];   % lower bounds: [d₁_min,  k_min]
ub = [max(z)*10,     1e15 ];   % upper bounds: [d₁_max,  k_max]


% ============================================================
%  SECTION 8 — Multi-start nonlinear least-squares fitting
% ============================================================
% lsqcurvefit minimises  Σ [v_data(i) − v_model(z_i ; d₁, k)]²
% using the Levenberg-Marquardt/trust-region algorithm.
%
% Multiple initial guesses are tested to reduce sensitivity to local
% minima.  The candidate with the lowest residual SSE is kept.
%
% Initial guesses span the physically plausible parameter space:
%   d₁_init ≈ 30% of total depth (reasonable first estimate)
%   k_init  ~ m·v₀²/z_max  (kinetic energy balance rough estimate)

d1_init = max(z) * 0.30;
k_init  = m * v0^2 / max(z);

starts = [...
    d1_init,       k_init;       % primary guess
    d1_init*0.50,  k_init*5;     % shorter d₁, stronger drag
    d1_init*2.00,  k_init*0.10;  % longer  d₁, weaker drag
    max(z)*0.10,   k_init*10;    % very shallow d₁
    max(z)*0.50,   k_init*50  ]; % moderate d₁, high drag

opts = optimoptions('lsqcurvefit', ...
    'Display',                'off',    ...
    'MaxFunctionEvaluations',  50000,   ...
    'FunctionTolerance',       1e-14,   ...
    'StepTolerance',           1e-14);

fprintf('  Running multi-start fitting (%d starts)...\n', size(starts,1));

best_sse = Inf;
p_best   = [d1_init, k_init];   % safety fallback

for i = 1:size(starts, 1)
    try
        p_try = lsqcurvefit(eq2_vel, starts(i,:), z, v_use, lb, ub, opts);
        sse   = sum((v_use - eq2_vel(p_try, z)).^2);
        fprintf('    Start %d → d₁ = %.3e m,  k = %.3e,  SSE = %.6e\n', ...
                i, p_try(1), p_try(2), sse);
        if sse < best_sse
            best_sse = sse;
            p_best   = p_try;
        end
    catch ME
        fprintf('    Start %d → FAILED (%s)\n', i, ME.message);
    end
end

d1_opt = p_best(1);
k_opt  = p_best(2);

fprintf('\n  Best solution: d₁ = %.6e m,  k = %.6e kg/m\n\n', d1_opt, k_opt);


% ============================================================
%  SECTION 9 — Goodness-of-fit statistics
% ============================================================
% R² = 1 − SS_res/SS_tot
%   SS_res : sum of squared residuals between data and model
%   SS_tot : total variance of data around its mean
%   R²→1   : model explains nearly all variance (excellent fit)
%   R²→0   : model does no better than a horizontal mean line
%
% RMSE [m/s]: root-mean-square of velocity residuals.
%   Gives an intuitive, physics-unit measure of fit quality.
%   RMSE < 0.01 m/s would be <0.3% of v₀ = very good for DEM data.

v_fit  = eq2_vel(p_best, z);

SS_res = sum((v_use - v_fit).^2);
SS_tot = sum((v_use - mean(v_use)).^2);
R2     = 1 - SS_res / SS_tot;
rmse   = sqrt(mean((v_use - v_fit).^2));


% ============================================================
%  SECTION 10 — Final console summary
% ============================================================
fprintf('=======================================================\n');
fprintf('  KATSURAGI Eq.2 FIT RESULTS  (SI units throughout)\n');
fprintf('=======================================================\n');
fprintf('  d₁   = %.6e  m\n',    d1_opt);
fprintf('  k    = %.6e  kg/m\n', k_opt);
fprintf('  m    = %.6e  kg\n',   m);
fprintf('  v₀   = %.6e  m/s\n',  v0);
fprintf('  g    = %.2f         m/s²\n', g);
fprintf('  R²   = %.6f\n',       R2);
fprintf('  RMSE = %.6e  m/s\n',  rmse);
fprintf('  N    = %d  data points\n', length(z));
fprintf('  z_max= %.6f  m  (= %.2f mm)\n', max(z), max(z)*1000);
fprintf('=======================================================\n\n');


% ============================================================
%  SECTION 11 — Publication-quality plot: v vs z
%               (Katsuragi 2013 paper style)
% ============================================================
% x-axis: penetration depth z [m]  — increases right (deeper)
% y-axis: velocity v [m/s]         — decreases from v₀ at z=0
%
% Experimental data → thin solid blue line (dense LAMMPS time-series)
% Katsuragi Eq.2 fit → thick solid red line (smooth theoretical curve)

z_fine   = linspace(0, max(z), 3000)';
v_theory = eq2_vel(p_best, z_fine);

fig = figure('Color',       'white',   ...
             'Units',       'inches',  ...
             'Position',    [1 1 8 5.5], ...
             'InvertHardcopy', 'off');

% --- Experimental / simulation data ---
plot(z, v_use, '-', ...
    'Color',       [0.15 0.35 0.85], ...
    'LineWidth',    1.2,              ...
    'DisplayName', 'DEM Simulation (LAMMPS)');
hold on;

% --- Katsuragi Eq.2 theoretical fit ---
plot(z_fine, v_theory, '-', ...
    'Color',       [0.85 0.12 0.12], ...
    'LineWidth',    2.8,              ...
    'DisplayName', sprintf('Katsuragi Eq.2   R^2 = %.4f', R2));

% --- Annotation box: fitted parameters ---
ann_lines = { ...
    sprintf('d_1  = %.4e  m',    d1_opt), ...
    sprintf('k    = %.4e  kg/m', k_opt),  ...
    sprintf('R^2  = %.4f',       R2),     ...
    sprintf('RMSE = %.4e  m/s',  rmse)    };

annotation('textbox', [0.54 0.54 0.33 0.24], ...
    'String',          ann_lines,            ...
    'FitBoxToText',    'on',                 ...
    'BackgroundColor', [0.98 0.98 0.92],     ...
    'EdgeColor',       [0.45 0.45 0.45],     ...
    'LineWidth',       1.2,                  ...
    'FontSize',        10,                   ...
    'FontName',        'Courier New',        ...
    'Color',           'black');

% --- Axis labels ---
xlabel('Penetration depth   z   (m)',  'FontSize', 14, 'FontWeight', 'normal');
ylabel('Velocity   v   (m/s)',         'FontSize', 14, 'FontWeight', 'normal');
title ('Projectile Velocity vs. Penetration Depth — Katsuragi (2013) Eq.2', ...
       'FontSize', 13, 'FontWeight', 'bold');

% --- Legend ---
lgd = legend('Location', 'northeast', ...
             'FontSize',   12,         ...
             'Box',        'on',       ...
             'Color',      'white',    ...
             'TextColor',  'black');

% --- Axis limits and grid ---
xlim([0,  max(z) * 1.02]);
ylim([0,  v0     * 1.12]);
grid on; box on;

set(gca, ...
    'FontSize',      12,      ...
    'LineWidth',     1.2,     ...
    'GridAlpha',     0.20,    ...
    'GridLineStyle', '--',    ...
    'Color',         'white', ...
    'XColor',        'black', ...
    'YColor',        'black', ...
    'TickDir',       'out',   ...
    'TickLength',    [0.012 0.012]);

% --- Save figure (300 dpi, white background) ---
exportgraphics(fig, 'katsuragi_eq2_velocity_fit.png', ...
    'Resolution',      300,     ...
    'BackgroundColor', 'white');

fprintf('Figure saved  →  katsuragi_eq2_velocity_fit.png\n');
fprintf('Total fitted data points  :  %d\n', length(z));
