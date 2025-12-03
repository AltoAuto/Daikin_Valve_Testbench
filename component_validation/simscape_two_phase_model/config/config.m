%% Project startup (auto-runs when the .prj opens)
% config/config.m

clc;  fprintf('Loading Valve Test System\n');

% Load project path
projRoot = fileparts(mfilename('fullpath')); 
projRoot = fileparts(projRoot);
addpath(fullfile(projRoot,'models'));
addpath(fullfile(projRoot,'scripts'));

%% set global varible
cfg = struct();
cfg.t_stop   = 60;          % [s] total simulation times 
cfg.dt_local = 1e-3;        % Simscape local solver sample time

cfg.comp.rpm_nom = 0.045;

% CFG for EXV valve demo (manual sweep + PI control)
% Units noted per field. 

%% -------------------- Physics Coil (Preparation) -----------------------
% Toggle: 1 = use geometry-based coil variant; 0 = keep system-level coil
cfg.coil.use_physics     = 1;

% Geometry (SI)
cfg.coil.D_i_m           = 12e-3;    % [m] inner diameter (e.g., 1/4" tube)
cfg.coil.L_tube_m        = 10.0;        % [m] length per tube
cfg.coil.N_tubes         = 20;          % [-] tubes in parallel
cfg.coil.rough_eps_m     = 1.5e-5;     % [m] roughness (drawn copper ~1e-5..3e-5)

% Air-side boundary for demo (replace with moist-air later)
cfg.coil.T_air_K         = 297.15;     % [K] ~24 °C

% Property wiring mode for HTC/Δp MATLAB Functions
%  "from_bus" -> signals rho, mu, k, Pr provided by property block at coil inlet
%  "fixed"    -> use these fallback constants (demo only)
cfg.coil.props_mode      = "from_bus";
cfg.coil.rho_ref_kgm3 = 30;     % Refrigerant density inside coil [kg/m^3]
cfg.coil.mu_Pas       = 1e-5;  % Dynamic viscosity of refrigerant [Pa·s]
cfg.coil.k_WmK        = 0.014;    % Thermal conductivity of refrigerant [W/(m·K)]
cfg.coil.Pr           = 1;     % Prandtl number (dimensionless)

% Solver guard rails for demo functions (avoid laminar corner cases)
cfg.coil.Re_min          = 3000;       % [-] clamp Re >= this in HTC/Δp functions

% Test sweep knobs (used by run scripts)
cfg.test.mdot_min_kgps   = 0.02;       % [kg/s]
cfg.test.mdot_max_kgps   = 0.12;       % [kg/s]
cfg.test.mdot_N          = 6;          % points in sweep

%% -------------------- Mode Selection --------------------
% 'manual'(0)  : open-loop step/sweep on EXV opening
% 'pi' (1)     : closed-loop SH control using PI
cfg.mode = 0;

%% -------------------- EXV Valve Model --------------------
% Mapping mode: 'area_sqrt' | 'area_table' | 'cv_table'
cfg.exv.mode = "area_sqrt";

% Actuator & limits
cfg.exv.u_min           = 0.05;    % [-]  lower opening limit
cfg.exv.u_max           = 0.95;    % [-]  upper opening limit
cfg.exv.rate_lim_per_s  = 0.03;    % [1/s] max slew of opening command
cfg.exv.tau_s           = 0.10;    % [s]   first-order actuator lag

% Effective area mapping (used if mode = 'area_sqrt' or 'area_table')
% A(u) = Amin + (Amax - Amin)*sqrt(u)   (good proxy since m_dot ~ A*sqrt(Δp))
%1.26e-5 2.4e-6  2.12e-5
cfg.exv.Amax_m2         =  2.4e-6;                    % [m^2] fully open eff. area
cfg.exv.Amin_m2         = 0.00005 * cfg.exv.Amax_m2;    % [m^2] leakage / never zero

% Optional explicit table (used if mode = 'area_table'); keep monotone
cfg.exv.area_map.u      = [0, 1];                    % [-]
cfg.exv.area_map.A_m2   = [cfg.exv.Amin_m2, cfg.exv.Amax_m2];

% Optional Cv table stub (used if mode = 'cv_table'); replace with real data
cfg.exv.cv_map.u        = [0, 1];                    % [-]
cfg.exv.cv_map.Cv       = [0, 1];                    % [US gal/min)/(psi^0.5)] (placeholder)

%% -------------------- Manual Sweep (Open-Loop Demo) --------------------
cfg.demo.step.u1        = 0.33;    % [-] initial opening
cfg.demo.step.u2        = 0.37;    % [-] final opening
cfg.demo.step.t_s       = 20;      % [s]  step time

%% -------------------- Superheat Control (Closed-Loop) ------------------
% Setpoint & ramp
cfg.sh.sp_K             = 4.6;     % [K] nominal SH setpoint
cfg.sh.ramp_start_K     = 12.0;    % [K] initial SP for soft-start
cfg.sh.ramp_time_s      = 15.0;    % [s] ramp from ramp_start_K -> sp_K
cfg.sh.lpf_tau_s        = 0.20;    % [s] optional SH low-pass

% PI controller (no D). Output is opening command in [u_min, u_max].
cfg.pi.Kp               = 0.12;    % [-/K]
cfg.pi.Ki               = 0.008;   % [1/s*K]
cfg.pi.Kaw              = 0.3 / cfg.pi.Kp;  % back-calc gain (heuristic)

% Safety clamps for SH signal (controller input)
cfg.sh.clip_min_K       = 0.0;     % [K]
cfg.sh.clip_max_K       = 30.0;    % [K]

%% -------------------- Housekeeping ----------------------
% Everything below here is convenience/meta 
cfg.meta.name           = "EXV Demo – R32";
cfg.meta.version        = "v1.0";
cfg.meta.notes          = "Area sqrt mapping; tune Kp/Ki on rig.";

% assign to base workspace
assignin('base','cfg',cfg);

% Buslog
if ~evalin('base','exist(''BusLogs'',''var'')') 
    BusLogs = Simulink.Bus;
    names = {'P1','P2','T_suct','P_suct','mdot','Opening','SH'};
    for i=1:numel(names)
        el = Simulink.BusElement;
        el.Name = names{i};
        el.DataType = 'double';
        el.Dimensions = 1;
        el.SampleTime = -1;
        el.Complexity = 'real';
        BusLogs.Elements(end+1) = el;
    end
    assignin('base','BusLogs',BusLogs);
end

fprintf('CFG loaded');
