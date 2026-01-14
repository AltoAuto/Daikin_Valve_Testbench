%% Phase 1 – Condenser TL–2P Harness Sweep (FULL SETUP)
% -------------------------------------------------------------------------
% Purpose:
%   Stand-alone condenser harness for validation.
%
%   Sweeps:
%     - Refrigerant mass flow: mdot_ref_set
%     - TL inlet temperature: T_TL_in_set
%
%   For each operating point, the script:
%     - Simulates "condenser_model"
%     - Computes condenser duty Q_cond from refrigerant Δh
%     - Records TL-side ΔT and refrigerant-side Δp
%     - Writes a combined CSV of all cases
%     - Saves a plot Q_cond vs mdot_ref for different T_TL_in
% -------------------------------------------------------------------------

clear; clc;

%% === Output folder setup ===
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
OUTDIR   = fullfile(thisDir, 'outputs');

if ~exist(OUTDIR, 'dir')
    mkdir(OUTDIR);
end

modelName = "condenser_model";

%% === Global design constants ===

% --- Geometry (refrigerant side) ---
cond_tube_D_o      = 0.010;     % [m] tube outer diameter
cond_tube_D_i      = 0.008;     % [m] tube inner diameter
cond_tube_L        = 9.5000;    % [m] tube length per tube
cond_2p_N_tubes    = 22;        % [-] number of parallel tubes
cond_wall_thickness = (cond_tube_D_o - cond_tube_D_i) / 2;  % [m]
rho_Cu             = 8960;      % [kg/m^3] copper density

% Wall mass calculation (annulus)
A_wall_cross = (pi/4) * (cond_tube_D_o^2 - cond_tube_D_i^2);      % [m^2]
m_wall_total = rho_Cu * cond_2p_N_tubes * A_wall_cross * cond_tube_L; % [kg]

% Flow area for refrigerant
crossec_area_ref = cond_2p_N_tubes * (pi * cond_tube_D_i^2 / 4);  % [m^2]

% --- Thermal liquid (water) side ---
cond_TL_N_tubes = 13;      % [-]
cond_wall_k     = 380;     % [W/m/K]
cond_Dh_TL      = 0.010;   % [m] TL hydraulic diameter

cond_fouling_2P = 0.0;     % clean coil, 2P side
cond_fouling_TL = 0.0;     % clean coil, TL side

% --- Nominal operating conditions ---
cp_water = 4180;           % [J/kg/K]
rho_water = 1000;          % [kg/m^3]
v_TL      = 2.0;           % [m/s] (info only)

Q_evap_nom   = 4 * 3.517e3;        % [W] 4 ton evaporator
mdot_ref_nom = 0.05;               % [kg/s] nominal refrigerant mass flow

P_comp_nom   = 2.22e3;             % [W] compressor power
Q_cond_nom   = Q_evap_nom + P_comp_nom;   % [W] condenser duty (~16.32 kW)

SH_comp_nom  = 11.1;        % [K] superheat (Note: not used, for me to remember)
SC_cond_nom  = 8.3;         % [K] subcooling

p_evap_nom   = 1.0e6;       % [Pa]
p_cond_nom   = 2.5e6;       % [Pa]

T_evap_sat_nom_C = 5;       % [°C] (Note: not used, for me to remember)
T_cond_sat_nom_C = 45;      % [°C]

T_comp_in_nom_C  = T_evap_sat_nom_C + SH_comp_nom;   % [°C] (Note: not used, for me to remember)
T_comp_out_nom_C = T_cond_sat_nom_C + 20;            % [°C] assume +20 K SH
T_cond_out_nom_C = T_cond_sat_nom_C - SC_cond_nom;   % [°C]

% --- Condenser water-side design ---
T_TL_cond_in_nom_C  = 30;   % [°C] design inlet
T_TL_cond_out_nom_C = 35;   % [°C] design outlet
dT_TL_cond_nom      = T_TL_cond_out_nom_C - T_TL_cond_in_nom_C; % [K] = 5

mdot_TL_cond = Q_cond_nom / (cp_water * dT_TL_cond_nom);  % [kg/s] water flow

%% === Sweep setup ===
% Refrigerant mass flow sweep [kg/s]
mdot_ref_vec = mdot_ref_nom * [0.5, 1.0, 1.5];

% Condenser water inlet temperatures [°C]
T_TL_in_vec  = [26, 30, 34];

n_mdot = numel(mdot_ref_vec);
n_Tin  = numel(T_TL_in_vec);

% Preallocate result matrices
Q_res_kW   = zeros(n_mdot, n_Tin);  % [kW] condenser duty
dT_TL_res  = zeros(n_mdot, n_Tin);  % [K]  water ΔT
dp_res     = zeros(n_mdot, n_Tin);  % [Pa] refrigerant Δp

%% === Main sweep ===
for i = 1:n_mdot
    for j = 1:n_Tin

        % Set harness inputs for this operating point
        mdot_ref_set = mdot_ref_vec(i);   % [kg/s] → 2P mass flow source
        T_TL_in_set  = T_TL_in_vec(j);    % [°C]  → TL inlet reservoir

        % Run the condenser-only model
        simOut = sim(modelName, "ReturnWorkspaceOutputs","on");

        % --- Refrigerant enthalpies [kJ/kg] ---
        h_cond_in_sig  = simOut.logsout.getElement("h_cond_in").Values;
        h_cond_out_sig = simOut.logsout.getElement("h_cond_out").Values;

        h_in  = h_cond_in_sig.Data(end);
        h_out = h_cond_out_sig.Data(end);

        % Q_cond in kW: mdot_ref_set [kg/s] * Δh [kJ/kg]
        Q_res_kW(i,j) = mdot_ref_set * (h_out - h_in);

        % --- Refrigerant pressure drop [Pa] ---
        p_in_sig  = simOut.logsout.getElement("p_cond_in").Values;
        p_out_sig = simOut.logsout.getElement("p_cond_out").Values;

        p_in  = p_in_sig.Data(end);
        p_out = p_out_sig.Data(end);
        dp_res(i,j) = p_in - p_out;

        % --- TL-side temperature rise [K] ---
        T_in_sig  = simOut.logsout.getElement("T_TL_in").Values;
        T_out_sig = simOut.logsout.getElement("T_TL_out").Values;

        T_in  = T_in_sig.Data(end);
        T_out = T_out_sig.Data(end);

        dT_TL_res(i,j) = T_out - T_in;
    end
end

%% === Display results in tables ===

disp("Condenser Heat Transfer Q_{cond} (kW):");
T_Q = array2table(Q_res_kW, ...
    'VariableNames', compose('T_TL_in_%dC', T_TL_in_vec), ...
    'RowNames',      compose('mdot_ref_%.3f', mdot_ref_vec));
disp(T_Q);

disp("Condenser Water ΔT_{TL} (K):");
T_dT = array2table(dT_TL_res, ...
    'VariableNames', compose('T_TL_in_%dC', T_TL_in_vec), ...
    'RowNames',      compose('mdot_ref_%.3f', mdot_ref_vec));
disp(T_dT);

disp("Condenser Refrigerant Δp_cond (Pa):");
T_dp = array2table(dp_res, ...
    'VariableNames', compose('T_TL_in_%dC', T_TL_in_vec), ...
    'RowNames',      compose('mdot_ref_%.3f', mdot_ref_vec));
disp(T_dp);

%% === Save combined CSV (all operating points) ===

% Build grid of operating points
[mdot_grid, Ttl_grid] = ndgrid(mdot_ref_vec, T_TL_in_vec);

% Flatten grids and results
RefrigerantMassFlow_kg_s = mdot_grid(:);
WaterInletTemp_C          = Ttl_grid(:);
CondenserCapacity_kW      = Q_res_kW(:);
WaterTempRise_K           = dT_TL_res(:);
RefrigerantDeltaP_Pa      = dp_res(:);

% Create table
T_all = table(RefrigerantMassFlow_kg_s, WaterInletTemp_C, ...
              CondenserCapacity_kW, WaterTempRise_K, ...
              RefrigerantDeltaP_Pa);

% Save to CSV
outfile = fullfile(OUTDIR, 'Condenser_SweepResults.csv');
writetable(T_all, outfile);

fprintf('CSV saved to: %s\n', outfile);

%% === Plot Q_cond vs mdot_ref for each T_TL_in ===

figure;
for j = 1:n_Tin
    plot(mdot_ref_vec, Q_res_kW(:,j), '-o', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T_{TL,in} = %d°C', T_TL_in_vec(j)));
    hold on;
end

xlabel('Refrigerant mass flow \dot{m}_{ref} (kg/s)');
ylabel('Condenser duty Q_{cond} (kW)');
title('Condenser TL–2P Harness – Q_{cond} vs \dot{m}_{ref}');
legend('Location','northwest');
grid on;

saveas(gcf, fullfile(OUTDIR, 'Phase1_cond_Q_vs_mdot.png'));

fprintf('Condenser results saved in %s\n', OUTDIR);
