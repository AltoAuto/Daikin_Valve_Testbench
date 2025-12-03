%% Phase 1 – Condenser TL–2P Geometry Optimizer (parsim version)
% -------------------------------------------------------------------------
% Purpose:
%   Stand-alone condenser harness for validation & sizing.
%
%   We fix the compressor design-point refrigerant mass flow (mdot_ref_nom)
%   and nominal entering water temperature (T_TL_cond_in_nom_C), then:
%     - Sweep condenser 2P geometry (D_o, D_i, L, N_tubes)
%     - Sweep condenser TL geometry (cond_Dh_TL, cond_TL_N_tubes, mdot_TL_cond)
%   For each geometry, we:
%     - Run "condenser_model" with mdot_ref_set = mdot_ref_nom
%     - Extract Q_cond, ΔT_TL, and Δp_cond
%     - Compute a cost vs target Q_cond and target ΔT_TL
%   The script:
%     - Reports the best geometry at the design point
%     - Exports ALL cases to a CSV with readable names and units
% -------------------------------------------------------------------------

clear; clc;

thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
OUTDIR   = fullfile(thisDir, 'outputs');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

modelName = "condenser_model";

%% =============================
%  Global / nominal design
%% =============================

% Refrigerant nominal (COMPRESSOR DESIGN POINT)
mdot_ref_nom     = 0.05;                % [kg/s]
Q_evap_nom       = 4 * 3.517e3;         % [W] 4 ton
P_comp_nom       = 2.22e3;              % [W]
Q_cond_nom       = Q_evap_nom + P_comp_nom;    % [W]
Q_target_kW      = Q_cond_nom / 1e3;           % [kW] ~16.29

% Water-side nominal
cp_water            = 4180;             % [J/kg/K]
T_TL_cond_in_nom_C  = 30;               % [°C] design entering water temp
T_TL_cond_out_nom_C = 35;               % [°C] design leaving temp
dT_TL_target        = T_TL_cond_out_nom_C - T_TL_cond_in_nom_C;  % [K] = 5
rho_Cu              = 8960;             % [kg/m^3] copper density

% Pressures
p_cond_nom       = 2.5e6;               % [Pa]

% Temperatures (cycle)
T_cond_sat_nom_C = 45;                  % [°C]
SH_comp_nom      = 11.1;                % [K]
SC_cond_nom      = 8.3;                 % [K]
T_comp_out_nom_C = T_cond_sat_nom_C + 20;           % [°C] assume +20 K SH
T_cond_out_nom_C = T_cond_sat_nom_C - SC_cond_nom;  % [°C]

% Fouling (both sides)
cond_fouling_2P = 0.0;
cond_fouling_TL = 0.0;

%% =============================
%       SWEEP DIMENSIONS
%% =============================

% Refrigerant-side geometry (2P)
cond_tube_D_o_vec   = [0.010 0.011];    % [m] OD sweep
cond_tube_D_i_vec   = 0.008;            % [m] ID fixed
cond_tube_L_vec     = [9.5 10.5];       % [m] tube length
cond_2p_N_tubes_vec = 22;               % [-] number of 2P tubes

% Water-side geometry (near nominal)
mdot_TL_cond_nom = Q_cond_nom / (cp_water * dT_TL_target);   % [kg/s]
mdot_TL_cond_vec = mdot_TL_cond_nom;                         % single value
cond_Dh_TL_vec   = 0.010;                                    % [m]
cond_TL_N_tubes_vec = [13 16];                               % [-]

% TL inlet temperature fixed at design value
T_TL_in_set = T_TL_cond_in_nom_C;                            % [°C]

%% =============================
%   Build parameter grid
%% =============================

[Do_g, Di_g, L_g, N2p_g, TLflow_g, DhTL_g, NTL_g] = ndgrid( ...
    cond_tube_D_o_vec, ...
    cond_tube_D_i_vec, ...
    cond_tube_L_vec, ...
    cond_2p_N_tubes_vec, ...
    mdot_TL_cond_vec, ...
    cond_Dh_TL_vec, ...
    cond_TL_N_tubes_vec );

Do_list      = Do_g(:);      % [m] OD
Di_list      = Di_g(:);      % [m] ID
L_list       = L_g(:);       % [m] length
N2p_list     = N2p_g(:);     % [-] 2P tubes
mdotTL_list  = TLflow_g(:);  % [kg/s] TL flow
DhTL_list    = DhTL_g(:);    % [m] TL Dh
NTL_list     = NTL_g(:);     % [-] TL tubes

nCases = length(Do_list);
fprintf("Total geometry cases = %d\n", nCases);

%% =============================
%     Build parsim inputs
%% =============================

in(nCases,1) = Simulink.SimulationInput(modelName);

for k = 1:nCases

    % Unpack geometry for case k
    D_o       = Do_list(k);
    D_i       = Di_list(k);
    L_t       = L_list(k);
    N_2p      = N2p_list(k);
    mdot_TL   = mdotTL_list(k);
    Dh_TL     = DhTL_list(k);
    N_TL_tube = NTL_list(k);

    % Wall mass calc (2P side)
    A_wall_cross = (pi/4)*(D_o^2 - D_i^2);                 % [m^2]
    m_wall_total = rho_Cu * N_2p * A_wall_cross * L_t;     % [kg]

    % Refrigerant flow area (2P side)
    crossec_area_ref = N_2p * (pi * D_i^2 / 4);            % [m^2]

    % Create SimulationInput object
    s = Simulink.SimulationInput(modelName);

    % --- 2P geometry variables ---
    s = s.setVariable('cond_tube_D_o',    D_o);
    s = s.setVariable('cond_tube_D_i',    D_i);
    s = s.setVariable('cond_tube_L',      L_t);
    s = s.setVariable('cond_2p_N_tubes',  N_2p);
    s = s.setVariable('m_wall_total',     m_wall_total);
    s = s.setVariable('crossec_area_ref', crossec_area_ref);

    % --- TL geometry variables ---
    s = s.setVariable('cond_Dh_TL',      Dh_TL);
    s = s.setVariable('cond_TL_N_tubes', N_TL_tube);
    s = s.setVariable('mdot_TL_cond',    mdot_TL);
    s = s.setVariable('mdot_TL_set',     mdot_TL);

    % --- Boundary conditions (refrigerant side) ---
    s = s.setVariable('p_cond_nom',        p_cond_nom);
    s = s.setVariable('T_comp_out_nom_C',  T_comp_out_nom_C);
    s = s.setVariable('T_cond_out_nom_C',  T_cond_out_nom_C);

    % --- Compressor design-point (FIXED) ---
    s = s.setVariable('mdot_ref_set', mdot_ref_nom);

    % --- Water-side boundary (inlet temp FIXED at design) ---
    s = s.setVariable('T_TL_in_set', T_TL_in_set);

    % --- Fouling factors ---
    s = s.setVariable('cond_fouling_TL',  cond_fouling_TL);
    s = s.setVariable('cond_fouling_2P',  cond_fouling_2P);

    in(k) = s;
end

%% =============================
%           Run PARSIM
%% =============================

simOutArray = parsim(in, ...
    'UseFastRestart',"on", ...
    'ShowProgress',"on");

%% =============================
%        Evaluate results
%% =============================

bestCost = inf;
bestGeom = struct();
bestRes  = struct();

% Preallocate arrays for ALL CASES (for CSV export)
Refrigerant_TubeOD_m         = zeros(nCases,1);
Refrigerant_TubeID_m         = zeros(nCases,1);
Refrigerant_TubeLength_m     = zeros(nCases,1);
Refrigerant_NumParallelTubes = zeros(nCases,1);
Refrigerant_WallMass_kg      = zeros(nCases,1);
Refrigerant_FlowArea_m2      = zeros(nCases,1);

Water_NumTubes               = zeros(nCases,1);
Water_HydraulicDiameter_m    = zeros(nCases,1);
Water_MassFlow_kg_s          = zeros(nCases,1);
Water_InletTemp_C            = zeros(nCases,1);

Result_CondenserCapacity_kW  = zeros(nCases,1);
Result_WaterTempRise_K       = zeros(nCases,1);
Result_RefrigerantDeltaP_Pa  = zeros(nCases,1);
Result_Cost_J                = zeros(nCases,1);

for k = 1:nCases

    so = simOutArray(k);

    % Store input parameters for this case
    D_o       = Do_list(k);
    D_i       = Di_list(k);
    L_t       = L_list(k);
    N_2p      = N2p_list(k);
    mdot_TL   = mdotTL_list(k);
    Dh_TL     = DhTL_list(k);
    N_TL_tube = NTL_list(k);

    A_wall_cross = (pi/4)*(D_o^2 - D_i^2);
    m_wall_total = rho_Cu * N_2p * A_wall_cross * L_t;
    crossec_area_ref = N_2p * (pi * D_i^2 / 4);

    Refrigerant_TubeOD_m(k)         = D_o;
    Refrigerant_TubeID_m(k)         = D_i;
    Refrigerant_TubeLength_m(k)     = L_t;
    Refrigerant_NumParallelTubes(k) = N_2p;
    Refrigerant_WallMass_kg(k)      = m_wall_total;
    Refrigerant_FlowArea_m2(k)      = crossec_area_ref;

    Water_NumTubes(k)               = N_TL_tube;
    Water_HydraulicDiameter_m(k)    = Dh_TL;
    Water_MassFlow_kg_s(k)          = mdot_TL;
    Water_InletTemp_C(k)            = T_TL_in_set;

    if ~isempty(so.ErrorMessage)
        fprintf("Case %d failed: %s\n", k, so.ErrorMessage);

        % Mark result fields as NaN for this failed case
        Result_CondenserCapacity_kW(k) = NaN;
        Result_WaterTempRise_K(k)      = NaN;
        Result_RefrigerantDeltaP_Pa(k) = NaN;
        Result_Cost_J(k)               = NaN;
        continue;
    end

    % --- Read refrigerant-side enthalpies (kJ/kg) ---
    h_in  = so.logsout.getElement("h_cond_in").Values.Data(end);
    h_out = so.logsout.getElement("h_cond_out").Values.Data(end);

    % Q_cond in kW: mdot_ref_nom * Δh
    Q_kW = mdot_ref_nom * (h_out - h_in);

    % --- Water-side temperatures (°C) ---
    T_in  = so.logsout.getElement("T_TL_in").Values.Data(end);
    T_out = so.logsout.getElement("T_TL_out").Values.Data(end);
    dT_TL = T_out - T_in;   % [K]

    % --- Refrigerant pressure drop (Pa) ---
    p_in  = so.logsout.getElement("p_cond_in").Values.Data(end);
    p_out = so.logsout.getElement("p_cond_out").Values.Data(end);
    dp    = p_in - p_out;   % [Pa]

    % --- Cost function vs design targets ---
    err_Q  = Q_kW + Q_target_kW;   % [kW]  
    err_dT = dT_TL - dT_TL_target; % [K]

    J = 3*(err_Q)^2 + 1*(err_dT)^2;

    % Store results for this case
    Result_CondenserCapacity_kW(k) = Q_kW;
    Result_WaterTempRise_K(k)      = dT_TL;
    Result_RefrigerantDeltaP_Pa(k) = dp;
    Result_Cost_J(k)               = J;

    % Track best case
    if J < bestCost
        bestCost = J;

        bestGeom.cond_tube_D_o     = D_o;
        bestGeom.cond_tube_D_i     = D_i;
        bestGeom.cond_tube_L       = L_t;
        bestGeom.cond_2p_N_tubes   = N_2p;
        bestGeom.cond_Dh_TL        = Dh_TL;
        bestGeom.mdot_TL_cond      = mdot_TL;
        bestGeom.cond_TL_N_tubes   = N_TL_tube;

        bestRes.Q_kW  = Q_kW;
        bestRes.dT_TL = dT_TL;
        bestRes.dp    = dp;
    end
end

%% =============================
%             REPORT
%% =============================

fprintf("\n===== BEST CONDENSER GEOMETRY FOUND (Phase 1) =====\n");
disp(bestGeom);

fprintf("\nPerformance @ mdot_ref = %.3f kg/s and T_TL,in = %.1f °C:\n", ...
    mdot_ref_nom, T_TL_cond_in_nom_C);
fprintf("Q_cond(model) = %.2f kW (target %.2f kW)\n", bestRes.Q_kW, Q_target_kW);
fprintf("ΔT_TL(model)  = %.2f K  (target %.2f K)\n", bestRes.dT_TL, dT_TL_target);
fprintf("Δp_ref(model) = %.0f Pa\n", bestRes.dp);
fprintf("Cost J        = %.3e\n", bestCost);

%% =============================
%     EXPORT ALL CASES TO CSV
%% =============================

bestCase = struct();

% ---- Global / design parameters ----
bestCase.Design_RefrigerantMassFlow_kg_s      = mdot_ref_nom;
bestCase.Design_TargetCondenserCapacity_kW    = Q_target_kW;

bestCase.Design_CondensingPressure_Pa         = p_cond_nom;

bestCase.Design_CondSatTemp_C                 = T_cond_sat_nom_C;

bestCase.Design_CompressorOutletTemp_C        = T_comp_out_nom_C;
bestCase.Design_CondenserOutletTemp_C         = T_cond_out_nom_C;

bestCase.Design_WaterInletTemp_C              = T_TL_cond_in_nom_C;
bestCase.Design_WaterOutletTemp_C             = T_TL_cond_out_nom_C;
bestCase.Design_WaterTempRise_K               = dT_TL_target;

bestCase.Design_WaterCp_J_kgK                 = cp_water;
bestCase.Design_CopperDensity_kg_m3           = rho_Cu;

bestCase.Design_FoulingFactor_2P              = cond_fouling_2P;
bestCase.Design_FoulingFactor_TL              = cond_fouling_TL;

% ---- Best-case 2P condenser geometry ----
bestCase.TwoPhase_TubeOuterDiameter_m         = bestGeom.cond_tube_D_o;
bestCase.TwoPhase_TubeInnerDiameter_m         = bestGeom.cond_tube_D_i;
bestCase.TwoPhase_TubeLength_m                = bestGeom.cond_tube_L;
bestCase.TwoPhase_NumParallelTubes            = bestGeom.cond_2p_N_tubes;

% ---- Best-case TL-side geometry & operating condition ----
bestCase.TL_NumTubes                          = bestGeom.cond_TL_N_tubes;
bestCase.TL_HydraulicDiameter_m               = bestGeom.cond_Dh_TL;
bestCase.TL_MassFlow_kg_s                     = bestGeom.mdot_TL_cond;
bestCase.TL_InletTemp_C                       = T_TL_cond_in_nom_C;

% ---- Simulation results for best case ----
bestCase.Result_CondenserCapacity_kW          = bestRes.Q_kW;
bestCase.Result_WaterTempRise_K               = bestRes.dT_TL;
bestCase.Result_RefrigerantDeltaP_Pa          = bestRes.dp;
bestCase.Result_CostFunction_J                = bestCost;

% Convert struct -> 1-row table and write CSV
BestCaseTable = struct2table(bestCase);

csvName = fullfile(OUTDIR, 'Condenser_BestCase_AllParams.csv');
writetable(BestCaseTable, csvName);

fprintf('\nBest-case geometry + all condenser parameters saved to:\n  %s\n', csvName);