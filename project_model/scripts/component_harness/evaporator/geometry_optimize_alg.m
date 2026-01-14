%% Phase 1 – Evaporator TL–2P Geometry Optimizer (parsim version)
% -------------------------------------------------------------------------
% Purpose:
%   Stand-alone evaporator harness for validation & sizing.
%
%   We fix the compressor design-point refrigerant mass flow (mdot_ref_nom)
%   and nominal entering water temperature (T_TL_evap_in_nom_C), then:
%     - Sweep evaporator 2P geometry (D_o, D_i, L, N_tubes)
%     - Sweep evaporator TL geometry (evap_Dh_TL, evap_TL_N_tubes, mdot_TL_evap)
%   For each geometry, we:
%     - Run "evaporator_model" with mdot_ref_set = mdot_ref_nom
%     - Extract Q_evap, ΔT_TL, and Δp_evap
%     - Compute a cost vs target Q_evap and target ΔT_TL
%   The script:
%     - Reports the best geometry at the design point
%     - Exports ALL cases to a CSV with readable names and units
% -------------------------------------------------------------------------

clear; clc;

thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
OUTDIR   = fullfile(thisDir, 'outputs');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

modelName = "evaporator_model";   

%% =============================
%  Global / nominal design
% note: temperture coming in of the evap assumption is loose, lower than
% commented temps do not hit the design point delta T become negtaive, evap do not absorbe heat,
% looking in to corss flow instead of counter flow, use flow perpendicular instead maybe a good
% choice. Question to answer, since the temp are not ajustable comming in
% and out, what would you change to hit the design point Q? fouling factor
% some how played a huge role. The question to ask is how would you assume
% the temperature coming in and out, the data i found is through 2p full
% loop, ton is still low
%% =============================

% Refrigerant nominal (EVAPORATOR DESIGN POINT)
%mdot_ref_nom   = 0.172;                 % [kg/s] 10tons
mdot_ref_nom   = 0.01997;               % [kg/s] design refrigerant mass flow
%mdot_ref_nom   = 0.255;                  % [kg/s] design refrigerant mass flow
Q_evap_nom     = 0.8 * 3.517e3;           % [W] 4 ton ≈ 14.07 kW
Q_target_kW    = Q_evap_nom / 1e3;       % [kW]

% Evaporator-side pressure & temperatures (cycle)
p_evap_nom         = 0.972e6;            % [Pa] ~9.5 bar(a) at 5°C sat for R-32
T_evap_sat_nom_C   = 6.25;                % [°C] evaporating saturation temperature 6.15
SH_comp_nom        = 4.3;                % [K] desired superheat at evaporator outlet
T_evap_out_nom_C   = T_evap_sat_nom_C + SH_comp_nom;   % [°C] ≈ 13°C

% Water-side nominal (chilled water loop)
cp_water              = 4180;          % [J/kg/K]
T_TL_evap_in_nom_C    = 12;            % [°C] warm return from load
T_TL_evap_out_nom_C   = 7;             % [°C] chilled supply
dT_TL_target          = T_TL_evap_in_nom_C - T_TL_evap_out_nom_C;  % [K] = 5
rho_Cu                = 8960;          % [kg/m^3] copper density

% Water-side nominal mass flow (for target ΔT)
mdot_TL_evap_nom = 0.1;   % [kg/s] 0.8tons / 3 -> 15.1 tons

%mdot_TL_evap_nom = (Q_evap_nom) / (cp_water * dT_TL_target);   % [kg/s]
disp(mdot_TL_evap_nom);
% Fouling (both sides)
evap_fouling_2P  =0.1;
evap_fouling_TL  =0.1;

%% =============================
%       SWEEP DIMENSIONS
%% =============================

% Refrigerant-side geometry (2P)
evap_tube_D_o_vec   = 0.012;   % [m] OD sweep
evap_tube_D_i_vec   = 0.01;           % [m] ID fixed
evap_tube_L_vec     = 12;      % [m] tube length
evap_2p_N_tubes_vec = 22;              % [-] number of 2P tubes

% Water-side geometry (near nominal)
mdot_TL_evap_vec    = mdot_TL_evap_nom;    % single value for now
evap_Dh_TL_vec      = 0.010;               % [m]
evap_TL_N_tubes_vec = 18;             % [-]

% TL inlet temperature fixed at design value
T_TL_in_set = T_TL_evap_in_nom_C;          % [°C]

%% =============================
%   Build parameter grid
%% =============================

[Do_g, Di_g, L_g, N2p_g, TLflow_g, DhTL_g, NTL_g] = ndgrid( ...
    evap_tube_D_o_vec, ...
    evap_tube_D_i_vec, ...
    evap_tube_L_vec, ...
    evap_2p_N_tubes_vec, ...
    mdot_TL_evap_vec, ...
    evap_Dh_TL_vec, ...
    evap_TL_N_tubes_vec );

Do_list      = Do_g(:);      % [m] OD
Di_list      = Di_g(:);      % [m] ID
L_list       = L_g(:);       % [m] length
N2p_list     = N2p_g(:);     % [-] 2P tubes
mdotTL_list  = TLflow_g(:);  % [kg/s] TL flow
DhTL_list    = DhTL_g(:);    % [m] TL Dh
NTL_list     = NTL_g(:);     % [-] TL tubes

nCases = length(Do_list);
fprintf("Total evaporator geometry cases = %d\n", nCases);

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
    A_wall_cross   = (pi/4)*(D_o^2 - D_i^2);              % [m^2]
    m_wall_total   = rho_Cu * N_2p * A_wall_cross * L_t;  % [kg]

    % Refrigerant flow area (2P side)
    crossec_area_ref = N_2p * (pi * D_i^2 / 4) + 0.0006;           % [m^2]

    % Create SimulationInput object
    s = Simulink.SimulationInput(modelName);

    % --- 2P geometry variables (match these to your evap model workspace vars) ---
    s = s.setVariable('evap_tube_D_o',    D_o);
    s = s.setVariable('evap_tube_D_i',    D_i);
    s = s.setVariable('evap_tube_L',      L_t);
    s = s.setVariable('evap_2p_N_tubes',  N_2p);
    s = s.setVariable('evap_m_wall_total',     m_wall_total);
    s = s.setVariable('evap_crossec_area_ref', crossec_area_ref);

    % --- TL geometry variables ---
    s = s.setVariable('evap_Dh_TL',      Dh_TL);
    s = s.setVariable('evap_TL_N_tubes', N_TL_tube);
    s = s.setVariable('mdot_TL_evap',    mdot_TL);
    s = s.setVariable('mdot_TL_set',     mdot_TL);

    % --- Boundary conditions (refrigerant side) ---
    s = s.setVariable('p_evap_nom',       p_evap_nom);
    s = s.setVariable('T_evap_sat_nom_C', T_evap_sat_nom_C);
    s = s.setVariable('T_evap_out_nom_C', T_evap_out_nom_C);

    % --- Refrigerant design-point mass flow (FIXED) ---
    s = s.setVariable('mdot_ref_set', mdot_ref_nom);

    % --- Water-side boundary (inlet temp FIXED at design) ---
    s = s.setVariable('T_TL_in_set', T_TL_in_set);

    % --- Fouling factors ---
    s = s.setVariable('evap_fouling_TL',  evap_fouling_TL);
    s = s.setVariable('evap_fouling_2P',  evap_fouling_2P);

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

Result_EvaporatorCapacity_kW = zeros(nCases,1);
Result_WaterTempDrop_K       = zeros(nCases,1);
Result_RefrigerantDeltaP_Pa  = zeros(nCases,1);
Result_RefrigerantDeltaH_evap  =zeros(nCases,1);

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

    A_wall_cross     = (pi/4)*(D_o^2 - D_i^2);
    m_wall_total     = rho_Cu * N_2p * A_wall_cross * L_t;
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
        Result_EvaporatorCapacity_kW(k) = NaN;
        Result_WaterTempDrop_K(k)       = NaN;
        Result_RefrigerantDeltaP_Pa(k)  = NaN;
        Result_RefrigerantDeltaH_evap(k)  = NaN;
        Result_Cost_J(k)                = NaN;
        continue;
    end

    % --- Read refrigerant-side enthalpies (kJ/kg) ---
    % NOTE: these signal names must match your evaporator_model logsout
    h_in  = so.logsout.getElement("h_evap_in").Values.Data(end);
    h_out = so.logsout.getElement("h_evap_out").Values.Data(end);

    % Q_evap in kW: mdot_ref_nom * Δh  (h_out > h_in for evaporation)
    Q_kW = mdot_ref_nom * (h_out - h_in);

    % --- Water-side temperatures (°C) ---
    T_in  = so.logsout.getElement("T_TL_in").Values.Data(end);
    T_out = so.logsout.getElement("T_TL_out").Values.Data(end);
    dT_TL = T_in - T_out;   % [K] warm in → cold out

    % --- Refrigerant pressure drop (Pa) ---
    p_in  = so.logsout.getElement("p_evap_in").Values.Data(end);
    p_out = so.logsout.getElement("p_evap_out").Values.Data(end);
    dp    = p_in - p_out;   % [Pa]
    dh    = h_out - h_in;

    % --- Cost function vs design targets ---
    % FIXED BUG: use difference, not sum
    err_Q  = Q_kW - Q_target_kW;     % [kW]
    err_dT = dT_TL - dT_TL_target;   % [K]

    J = 3*(err_Q)^2 + 1*(err_dT)^2;

    % Store results for this case
    Result_EvaporatorCapacity_kW(k) = Q_kW;
    Result_WaterTempDrop_K(k)       = dT_TL;
    Result_RefrigerantDeltaP_Pa(k)  = dp;
    Result_RefrigerantDeltaH_evap(k)  = dh;

    Result_Cost_J(k)                = J;

    % Track best case
    if J < bestCost
        bestCost = J;

        bestGeom.evap_tube_D_o     = D_o;
        bestGeom.evap_tube_D_i     = D_i;
        bestGeom.evap_tube_L       = L_t;
        bestGeom.evap_2p_N_tubes   = N_2p;
        bestGeom.evap_Dh_TL        = Dh_TL;
        bestGeom.mdot_TL_evap      = mdot_TL;
        bestGeom.evap_TL_N_tubes   = N_TL_tube;
        bestGeom.CS = N_TL_tube * (pi *Dh_TL^2 / 4);
        bestGeom.CA = crossec_area_ref;
        bestRes.Q_kW  = Q_kW;
        bestRes.dT_TL = dT_TL;
        bestRes.dp    = dp;
        bestRes.dh    = dh;

    end
end

%% =============================
%             REPORT
%% =============================

fprintf("\n===== BEST EVAPORATOR GEOMETRY FOUND (Phase 1) =====\n");
disp(bestGeom);

fprintf("\nPerformance @ mdot_ref = %.3f kg/s and T_TL,in = %.1f °C:\n", ...
    mdot_ref_nom, T_TL_evap_in_nom_C);
fprintf("Q_evap(model) = %.2f kW (target %.2f kW)\n", bestRes.Q_kW, Q_target_kW);
fprintf("ΔT_TL(model)  = %.2f K  (target %.2f K)\n", bestRes.dT_TL, dT_TL_target);
fprintf("Δp_ref(model) = %.0f Pa\n", bestRes.dp);
fprintf("Δh_ref(model) = %.0f \n", bestRes.dh);
fprintf("Cost J        = %.3e\n", bestCost);


%% =============================
%     EXPORT ALL CASES TO CSV
%% =============================

bestCase = struct();

% ---- Global / design parameters ----
bestCase.Design_RefrigerantMassFlow_kg_s      = mdot_ref_nom;
bestCase.Design_TargetEvaporatorCapacity_kW   = Q_target_kW;

bestCase.Design_EvaporatingPressure_Pa        = p_evap_nom;
bestCase.Design_EvapSatTemp_C                 = T_evap_sat_nom_C;
bestCase.Design_EvaporatorOutletTemp_C        = T_evap_out_nom_C;

bestCase.Design_WaterInletTemp_C              = T_TL_evap_in_nom_C;
bestCase.Design_WaterOutletTemp_C             = T_TL_evap_out_nom_C;
bestCase.Design_WaterTempDrop_K               = dT_TL_target;

bestCase.Design_WaterCp_J_kgK                 = cp_water;
bestCase.Design_CopperDensity_kg_m3           = rho_Cu;

bestCase.Design_FoulingFactor_2P              = evap_fouling_2P;
bestCase.Design_FoulingFactor_TL              = evap_fouling_TL;

% ---- Best-case 2P evaporator geometry ----
bestCase.TwoPhase_TubeOuterDiameter_m         = bestGeom.evap_tube_D_o;
bestCase.TwoPhase_TubeInnerDiameter_m         = bestGeom.evap_tube_D_i;
bestCase.TwoPhase_TubeLength_m                = bestGeom.evap_tube_L;
bestCase.TwoPhase_NumParallelTubes            = bestGeom.evap_2p_N_tubes;

% ---- Best-case TL-side geometry & operating condition ----
bestCase.TL_NumTubes                          = bestGeom.evap_TL_N_tubes;
bestCase.TL_HydraulicDiameter_m               = bestGeom.evap_Dh_TL;
bestCase.TL_MassFlow_kg_s                     = bestGeom.mdot_TL_evap;
bestCase.TL_InletTemp_C                       = T_TL_evap_in_nom_C;

% ---- Simulation results for best case ----
bestCase.Result_EvaporatorCapacity_kW         = bestRes.Q_kW;
bestCase.Result_WaterTempDrop_K               = bestRes.dT_TL;
bestCase.Result_RefrigerantDeltaP_Pa          = bestRes.dp;
bestCase.Result_CostFunction_J                = bestCost;

% Convert struct -> 1-row table and write CSV
BestCaseTable = struct2table(bestCase);

csvName = fullfile(OUTDIR, 'Evaporator_BestCase_AllParams.csv');
writetable(BestCaseTable, csvName);

fprintf('\nBest-case geometry + all evaporator parameters saved to:\n  %s\n', csvName);
