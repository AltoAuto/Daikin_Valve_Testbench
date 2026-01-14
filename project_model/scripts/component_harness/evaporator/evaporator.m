%% Phase 1 – Evaporator TL–2P Harness Sweep (Operating Map)
% -------------------------------------------------------------------------
% Purpose:
%   Stand-alone evaporator harness for validation (fixed best geometry).
%
%   Sweeps:
%     - Refrigerant mass flow: mdot_ref_set
%     - TL inlet temperature: T_TL_in_set
%
%   For each operating point, the script:
%     - Simulates "evaporator_model"
%     - Computes evaporator duty Q_evap from refrigerant Δh
%     - Records TL-side ΔT and refrigerant-side Δp
%     - Writes a combined CSV of all cases
%     - Saves a plot Q_evap vs mdot_ref for different T_TL_in
% -------------------------------------------------------------------------

clear; clc;

thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
OUTDIR   = fullfile(thisDir, 'outputs');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

modelName = "evaporator_model";

%% =============================
%  Global / nominal design (same as optimizer)
%% =============================

% Refrigerant nominal (EVAPORATOR DESIGN POINT)
mdot_ref_nom   = 0.01997;                % [kg/s] design refrigerant mass flow
Q_evap_nom     = 0.8 * 3.517e3;          % [W] 4 ton ≈ 14.07 kW
Q_evap_nom_kW  = Q_evap_nom / 1e3;     % [kW]

% Evaporator-side pressure & temperatures (cycle)
p_evap_nom         = 9.55e5;           % [Pa]
T_evap_sat_nom_C   = 6.25;                % [°C]
SH_comp_nom        = 4.3;                % [K]
T_evap_out_nom_C   = T_evap_sat_nom_C + SH_comp_nom;   % [°C]

% Water-side nominal (chilled water loop)
cp_water              = 4180;          % [J/kg/K]
T_TL_evap_in_nom_C    = 12;            % [°C]
T_TL_evap_out_nom_C   = 7;             % [°C]
dT_TL_target          = T_TL_evap_in_nom_C - T_TL_evap_out_nom_C;  % [K] = 5
rho_Cu                = 8960;          % [kg/m^3]

% Water-side nominal mass flow (for target ΔT)
mdot_TL_evap_nom = Q_evap_nom / (cp_water * dT_TL_target);   % [kg/s]

% Fouling (both sides)
evap_fouling_2P  = 0.1;
evap_fouling_TL  = 0.1;

%% =============================
%  BEST EVAPORATOR GEOMETRY (from optimizer)
%% =============================

evap_tube_D_o   = 0.0100;   % [m]
evap_tube_D_i   = 0.0100;   % [m]
evap_tube_L     = 11;       % [m]
evap_2p_N_tubes = 22;       % [-]

evap_Dh_TL      = 0.0100;   % [m]
mdot_TL_evap    = 0.6731;   % [kg/s]
evap_TL_N_tubes = 17;       % [-]

% Derived quantities (same formulas as optimizer)
A_wall_cross          = (pi/4)*(evap_tube_D_o^2 - evap_tube_D_i^2);  % [m^2]
evap_m_wall_total     = 81; % [kg]
evap_crossec_area_ref = evap_2p_N_tubes * (pi * evap_tube_D_i^2 / 4) + 0.0006; % [m^2]

fprintf("\nUsing best evaporator geometry from Phase-1 optimizer:\n");
fprintf("  D_o = %.4f m, D_i = %.4f m, L = %.2f m, N_2p = %d\n", ...
    evap_tube_D_o, evap_tube_D_i, evap_tube_L, evap_2p_N_tubes);
fprintf("  TL: Dh = %.4f m, N_tubes = %d, mdot_TL = %.3f kg/s\n", ...
    evap_Dh_TL, evap_TL_N_tubes, mdot_TL_evap);

%% =============================
%  Sweep setup
%% =============================

% Refrigerant mass flow sweep [kg/s]
mdot_ref_vec = mdot_ref_nom * [0.5, 1.0, 1.5];

% Evaporator TL inlet temperatures [°C] around nominal
T_TL_in_vec  = [T_TL_evap_in_nom_C - 2, ...
                T_TL_evap_in_nom_C, ...
                T_TL_evap_in_nom_C + 2];

n_mdot = numel(mdot_ref_vec);
n_Tin  = numel(T_TL_in_vec);

Q_res_kW   = zeros(n_mdot, n_Tin);  % [kW]
dT_TL_res  = zeros(n_mdot, n_Tin);  % [K]
dp_res     = zeros(n_mdot, n_Tin);  % [Pa]

%% =============================
%  Main sweep (simple for-loops + assignin)
%% =============================

for i = 1:n_mdot
    for j = 1:n_Tin

        mdot_ref_set = mdot_ref_vec(i);   % [kg/s]
        T_TL_in_set  = T_TL_in_vec(j);    % [°C]

        % --- 2P geometry variables (same names as optimizer) ---
        assignin('base','evap_tube_D_o',    evap_tube_D_o);
        assignin('base','evap_tube_D_i',    evap_tube_D_i);
        assignin('base','evap_tube_L',      evap_tube_L);
        assignin('base','evap_2p_N_tubes',  evap_2p_N_tubes);
        assignin('base','evap_m_wall_total',     evap_m_wall_total);
        assignin('base','evap_crossec_area_ref', evap_crossec_area_ref);

        % --- ALSO define the names your error message wants ---
        assignin('base','m_wall_total',     evap_m_wall_total);
        assignin('base','crossec_area_ref', evap_crossec_area_ref);

        % --- TL geometry variables ---
        assignin('base','evap_Dh_TL',      evap_Dh_TL);
        assignin('base','evap_TL_N_tubes', evap_TL_N_tubes);
        assignin('base','mdot_TL_evap',    mdot_TL_evap);
        assignin('base','mdot_TL_set',     mdot_TL_evap);

        % TL mass flow vector name used in Flow Rate Source (TL)
        assignin('base','mdot_TL_evap_vec', mdot_TL_evap);

        % --- Boundary conditions (refrigerant side) ---
        assignin('base','p_evap_nom',       p_evap_nom);
        assignin('base','T_evap_sat_nom_C', T_evap_sat_nom_C);
        assignin('base','T_evap_out_nom_C', T_evap_out_nom_C);

        % --- Sweep variables ---
        assignin('base','mdot_ref_set', mdot_ref_set);
        assignin('base','T_TL_in_set',  T_TL_in_set);

        % --- Fouling factors ---
        assignin('base','evap_fouling_TL',  evap_fouling_TL);
        assignin('base','evap_fouling_2P',  evap_fouling_2P);

        % Run the model
        simOut = sim(modelName, "ReturnWorkspaceOutputs","on");

        if ~isempty(simOut.ErrorMessage)
            fprintf("Case (i=%d,j=%d) failed: %s\n", i, j, simOut.ErrorMessage);
            Q_res_kW(i,j)  = NaN;
            dT_TL_res(i,j) = NaN;
            dp_res(i,j)    = NaN;
            continue;
        end

        logsout = simOut.logsout;

        % --- Refrigerant enthalpies [kJ/kg] ---
        h_in  = logsout.getElement("h_evap_in").Values.Data(end);
        h_out = logsout.getElement("h_evap_out").Values.Data(end);

        % Q_evap in kW: mdot_ref_set [kg/s] * Δh [kJ/kg]
        Q_res_kW(i,j) = mdot_ref_set * (h_out - h_in);

        % --- Refrigerant pressure drop [Pa] ---
        p_in  = logsout.getElement("p_evap_in").Values.Data(end);
        p_out = logsout.getElement("p_evap_out").Values.Data(end);
        dp_res(i,j) = p_in - p_out;

        % --- TL-side temperature drop [K] (warm in → cold out) ---
        T_in  = logsout.getElement("T_TL_in").Values.Data(end);
        T_out = logsout.getElement("T_TL_out").Values.Data(end);
        dT_TL_res(i,j) = T_in - T_out;
    end
end

%% =============================
%  Display results in tables
%% =============================

disp("Evaporator Heat Transfer Q_{evap} (kW):");
T_Q = array2table(Q_res_kW, ...
    'VariableNames', compose('T_TL_in_%dC', round(T_TL_in_vec)), ...
    'RowNames',      compose('mdot_ref_%.3f', mdot_ref_vec));
disp(T_Q);

disp("Evaporator Water ΔT_{TL} (K):");
T_dT = array2table(dT_TL_res, ...
    'VariableNames', compose('T_TL_in_%dC', round(T_TL_in_vec)), ...
    'RowNames',      compose('mdot_ref_%.3f', mdot_ref_vec));
disp(T_dT);

disp("Evaporator Refrigerant Δp_{evap} (Pa):");
T_dp = array2table(dp_res, ...
    'VariableNames', compose('T_TL_in_%dC', round(T_TL_in_vec)), ...
    'RowNames',      compose('mdot_ref_%.3f', mdot_ref_vec));
disp(T_dp);

%% =============================
%  Save combined CSV
%% =============================

[mdot_grid, Ttl_grid] = ndgrid(mdot_ref_vec, T_TL_in_vec);

RefrigerantMassFlow_kg_s = mdot_grid(:);
WaterInletTemp_C         = Ttl_grid(:);
EvaporatorCapacity_kW    = Q_res_kW(:);
WaterTempDrop_K          = dT_TL_res(:);
RefrigerantDeltaP_Pa     = dp_res(:);

T_all = table(RefrigerantMassFlow_kg_s, WaterInletTemp_C, ...
              EvaporatorCapacity_kW, WaterTempDrop_K, ...
              RefrigerantDeltaP_Pa);

outfile = fullfile(OUTDIR, 'Evaporator_SweepResults.csv');
writetable(T_all, outfile);

fprintf('CSV saved to: %s\n', outfile);

%% =============================
%  Plot Q_evap vs mdot_ref
%% =============================

fig = figure;
for j = 1:n_Tin
    plot(mdot_ref_vec, Q_res_kW(:,j), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T_TL,in = %.1f C', T_TL_in_vec(j)));
    hold on;
end

xlabel('Refrigerant mass flow mdot_ref (kg/s)');
ylabel('Evaporator duty Q_{evap} (kW)');
title('Evaporator TL–2P Harness – Q_{evap} vs mdot_{ref}');
legend('Location','northwest');
grid on;

saveas(fig, fullfile(OUTDIR, 'Phase1_evap_Q_vs_mdot.png'));

fprintf('Evaporator sweep results saved in %s\n', OUTDIR);
