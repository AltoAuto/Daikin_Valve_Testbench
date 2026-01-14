%% ========================================================================
%   EXV Harness: Mass Flow vs. Opening & Pressure Drop
%   ----------------------------------------------------------------------
%   Purpose:
%       Stand-alone validation harness for the expansion valve (EXV) model.
%       We check that mass flow behaves sensibly as a function of:
%
%           • Valve opening fraction (10% → 100%)
%           • Pressure drop across the valve (Δp)
%
%       Outputs:
%           • Mass-flow tables in kg/s and lb/hr (Daikin units)
%           • Plot: mdot vs. opening at two Δp values
%           • PNG saved to /outputs folder
%
%   Notes:
%       – This is a "shape check" of the EXV map:
%           mdot should increase with opening and Δp.
%       – Detailed SH control authority is checked later.
% =========================================================================

clear; clc;
format compact;

%% ------------------------------------------------------------------------
%  File / directory setup
% ------------------------------------------------------------------------
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
OUTDIR   = fullfile(thisDir, 'outputs');
if ~exist(OUTDIR, 'dir'); mkdir(OUTDIR); end

modelName = "exv_model";   % EXV harness model name


%% ------------------------------------------------------------------------
%  EXV operating conditions (high side / low side)
% ------------------------------------------------------------------------

% Upstream pressure (high side) [Pa]
p_high_set = 2.5e6;                 % 2.5 MPa ≈ 25 bar (typical condensing)

% Downstream pressures (low side) [Pa]
p_low_vec      = [1.5e6, 1.0e6];    % 1.5 MPa and 1.0 MPa
temp_low_vec   = [21, 7];           % [°C] downstream saturation/air-side conditions

% Convert to bar for reporting
p_high_bar = p_high_set / 1e5;      % [bar]
p_low_bar  = p_low_vec  / 1e5;      % [bar]
dP_bar     = (p_high_set - p_low_vec) / 1e5;   % Δp in bar


%% ------------------------------------------------------------------------
%  EXV geometry / opening definition
% ------------------------------------------------------------------------

A_max = 7.5e-6;                     % [m^2] max valve restricted area
A_min = 0.02 * A_max;               % [m^2] min effective area (~2% of max)

% Valve opening fraction (0–1 → 10–100%)
opening_vec = [0.10 0.30 0.50 0.70 1.00];
nOpen = numel(opening_vec);
nPlow = numel(p_low_vec);

% Storage for mass flow (SI units: kg/s)
mdot_res = zeros(nOpen, nPlow);


%% ------------------------------------------------------------------------
%  Harness sweep: opening × downstream pressure
% ------------------------------------------------------------------------
for i = 1:nOpen
    for j = 1:nPlow

        opening_cmd  = opening_vec(i);   % EXV command (0–1)
        p_low_set    = p_low_vec(j);     % low-side pressure [Pa]
        temp_low_set = temp_low_vec(j);  % low-side temperature [°C]

        % Model is expected to read:
        %   p_high_set, p_low_set, temp_low_set, opening_cmd, A_max, A_min
        simOut = sim(modelName);

        % Extract EXV mass flow (steady-state value)
        mdot_sig      = simOut.logsout.getElement("mdot_exv").Values;
        mdot_res(i,j) = mdot_sig.Data(end);    % [kg/s]

    end
end


%% ------------------------------------------------------------------------
%  Convert mass flow to Daikin units (lb/hr)
% ------------------------------------------------------------------------
kgps_to_lbphr   = 3600 * 2.20462;
mdot_lbphr_res  = mdot_res * kgps_to_lbphr;


%% ------------------------------------------------------------------------
%  Display results — SI units (kg/s)
% ------------------------------------------------------------------------
disp("=== EXV Mass Flow Results (kg/s) ===");
T_SI = array2table(mdot_res, ...
    'VariableNames', { ...
        sprintf('DeltaP_%.1fbar', dP_bar(1)), ...
        sprintf('DeltaP_%.1fbar', dP_bar(2)) }, ...
    'RowNames', compose('Opening_%d_percent', round(opening_vec*100)));
disp(T_SI);


%% ------------------------------------------------------------------------
%  Display results — Daikin units (lb/hr)
% ------------------------------------------------------------------------
disp("=== EXV Mass Flow Results (lb/hr) ===");
T_LB = array2table(mdot_lbphr_res, ...
    'VariableNames', { ...
        sprintf('DeltaP_%.1fbar', dP_bar(1)), ...
        sprintf('DeltaP_%.1fbar', dP_bar(2)) }, ...
    'RowNames', compose('Opening_%d_percent', round(opening_vec*100)));
disp(T_LB);


%% ------------------------------------------------------------------------
%  Plot: mass flow vs opening (lb/hr)
% ------------------------------------------------------------------------
fig = figure('Color','w');

plot(opening_vec*100, mdot_lbphr_res(:,1), '-o', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('\\Delta p = %.1f bar', dP_bar(1))); hold on;

plot(opening_vec*100, mdot_lbphr_res(:,2), '-o', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('\\Delta p = %.1f bar', dP_bar(2)));

xlabel('Valve Opening (%)', 'FontSize', 12);
ylabel('Mass Flow (lb/hr)', 'FontSize', 12);
title('EXV Mass Flow vs Opening (Phase 1 Harness)', 'FontSize', 14);
legend('Location','northwest');
grid on;

% Save figure
fn = fullfile(OUTDIR,'EXV_MassFlow_vs_Opening.png');
print(fig, fn, '-dpng','-r600');

fprintf('\nEXV mass-flow harness plot saved to:\n  %s\n', fn);
fprintf('Phase-1 EXV harness sweep complete.\n');
