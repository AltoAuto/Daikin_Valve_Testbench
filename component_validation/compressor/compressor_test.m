%% ========================================================================
%   Compressor Harness: Mass-Flow vs. Speed & Discharge Pressure
%   ----------------------------------------------------------------------
%   Purpose:
%       Stand-alone validation harness to check that the compressor model
%       produces reasonable refrigerant mass-flow (ṁ) trends under:
%
%           • Multiple shaft speeds (50%, 100%, 150% of nominal)
%           • Multiple discharge pressures (2.2 MPa, 2.8 MPa)
%
%       Output:
%           • ṁ table in kg/s and Daikin unit (lb/hr)
%           • Plot of ṁ vs speed (lb/hr)
%           • PNG saved to /outputs folder
%
%   Notes:
%       – This harness is NOT measuring compressor efficiency or power.
%       – This is strictly a sanity test: "Does the compressor map
%         produce sensible mass-flow values and trends?"
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

modelName = "compressor_model";


%% ------------------------------------------------------------------------
%  Nominal operating point (design)
% ------------------------------------------------------------------------
N_nom = 60;                 % Nominal compressor speed (model units)
speed_factors = [0.5 1.0 1.5];
N_vec = N_nom * speed_factors;    % Sweep: [50%, 100%, 150%]


%% ------------------------------------------------------------------------
%  Discharge pressure sweep (Pa)
% ------------------------------------------------------------------------
p_dis_vec = [2.2e6, 2.8e6];       % Light load → heavy load


%% ------------------------------------------------------------------------
%  Preallocate arrays
% ------------------------------------------------------------------------
mdot_res = zeros(numel(N_vec), numel(p_dis_vec));   % Mass flow (kg/s)


%% ------------------------------------------------------------------------
%  Run harness sweep
% ------------------------------------------------------------------------
for i = 1:numel(N_vec)
    for j = 1:numel(p_dis_vec)

        % Set speed & discharge pressure for this case
        N_cmd     = N_vec(i);
        p_dis_set = p_dis_vec(j);

        % Run compressor_only model
        simOut = sim(modelName);

        % Extract compressor mass flow (steady-state value)
        mdot_sig = simOut.logsout.getElement("mdot_comp").Values;
        mdot_res(i,j) = mdot_sig.Data(end);

    end
end


%% ------------------------------------------------------------------------
%  Convert to Daikin units (lb/hr)
% ------------------------------------------------------------------------
kgps_to_lbphr = 3600 * 2.20462;
mdot_lbphr_res = mdot_res * kgps_to_lbphr;


%% ------------------------------------------------------------------------
%  Display results — SI units
% ------------------------------------------------------------------------
disp("=== Compressor Mass Flow Results (kg/s) ===");
T_SI = array2table(mdot_res, ...
    'VariableNames', {'pDis_2_2MPa','pDis_2_8MPa'}, ...
    'RowNames', {'50pct_speed','100pct_speed','150pct_speed'});
disp(T_SI);


%% ------------------------------------------------------------------------
%  Display results — Daikin units
% ------------------------------------------------------------------------
disp("=== Compressor Mass Flow Results (lb/hr) – Daikin Units ===");
T_LB = array2table(mdot_lbphr_res, ...
    'VariableNames', {'pDis_2_2MPa','pDis_2_8MPa'}, ...
    'RowNames', {'50pct_speed','100pct_speed','150pct_speed'});
disp(T_LB);


%% ------------------------------------------------------------------------
%  Plot: ṁ vs speed (in Daikin units)
% ------------------------------------------------------------------------
fig = figure('Color','w');

plot(N_vec, mdot_lbphr_res(:,1), '-o', 'LineWidth', 1.8, ...
     'DisplayName','p_{dis} = 2.2 MPa'); hold on;

plot(N_vec, mdot_lbphr_res(:,2), '-o', 'LineWidth', 1.8, ...
     'DisplayName','p_{dis} = 2.8 MPa');

xlabel('Compressor Speed (model units)', 'FontSize', 12);
ylabel('Mass Flow (lb/hr)', 'FontSize', 12);
title('Compressor Mass Flow vs Speed (Daikin Units)', 'FontSize', 14);
legend('Location','northwest');
grid on;

% Save figure (600 dpi)
fn = fullfile(OUTDIR,'Compressor_MassFlow_vs_Speed_Daikin.png');
print(fig, fn, '-dpng','-r600');

fprintf("\nPlot saved to: %s\n", fn);
fprintf("Phase-1 compressor harness complete.\n");
