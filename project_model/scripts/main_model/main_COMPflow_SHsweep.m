%% TXV Omega Sweep – Steady-State SH vs mdot Map with SH-Band Operating Regions
clc;
run(fullfile('config','config.m'));
cfg = evalin('base','cfg');

model = 'two_phase_main_loop';
load_system(model);

%% =========================
% Sweep setup
%% =========================
omega_min = 0.5e3;
omega_max = 2.6e4;
nCases    = 20;

omega_vec = linspace(omega_min, omega_max, nCases);

simIn = repmat(Simulink.SimulationInput(model), nCases, 1);

for k = 1:nCases
    omega_cmd = omega_vec(k);

    simIn(k) = simIn(k) ...
        .setVariable('omega_in', omega_cmd) ...
        .setModelParameter( ...
            'StopTime','2000', ...
            'SignalLogging','on', ...
            'SignalLoggingName','logsout', ...
            'SaveFormat','Dataset');
end

simOut = parsim(simIn, ...
    'ShowProgress','on', ...
    'TransferBaseWorkspaceVariables','on');

%% =========================
% Operating region definitions
%% =========================
% SH bands (what you want shaded)
SH_good_lo    = 4;    % K
SH_good_hi    = 10;   % K
SH_sketch_lo  = 1;    % K
SH_sketch_hi  = 4;    % K
SH_danger_hi  = 1;    % K  (<1 is danger)

% "quality/stability" flags (optional but useful)
SH_std_max  = 1.0;    % K (hunting)
SH_pp_max   = 5.0;    % K (hunting)
mdot_min    = 0.008;  % kg/s (starvation)
N           = 50;     % last-N samples for averaging

%% =========================
% Results table
%% =========================
results = table('Size',[nCases 11], ...
  'VariableTypes',{'double','double','double','double','double','string','logical','logical','logical','string','double'}, ...
  'VariableNames',{'omega_cmd','mdot_mean','SH_mean','SH_std','SH_pp','SH_band','flagFloodback','flagStarved','flagHunting','flagSummary','caseIndex'});

for k = 1:nCases
    results.caseIndex(k) = k;
    results.omega_cmd(k) = omega_vec(k);

    % Grab logs
    logsout = simOut(k).logsout;

    % NOTE: element(14) is fragile; prefer a name when you can.
    md_vec = logsout.getElement(14).Values.Data;
    SH_vec = logsout.getElement('SH').Values.Data;

    nmd = min(N, numel(md_vec));
    nSH = min(N, numel(SH_vec));

    tailmd = md_vec(end-nmd+1:end);
    tailSH = SH_vec(end-nSH+1:end);

    md_mean = mean(tailmd);
    SH_mean = mean(tailSH);
    SH_std  = std(tailSH);
    SH_pp   = max(tailSH) - min(tailSH);

    results.mdot_mean(k) = md_mean;
    results.SH_mean(k)   = SH_mean;
    results.SH_std(k)    = SH_std;
    results.SH_pp(k)     = SH_pp;

    % --- SH band classification (this replaces your 0/1) ---
    if SH_mean < SH_danger_hi
        band = "Danger (<1K)";
    elseif SH_mean >= SH_sketch_lo && SH_mean < SH_good_lo
        band = "Sketchy (1–4K)";
    elseif SH_mean >= SH_good_lo && SH_mean <= SH_good_hi
        band = "Good (4–10K)";
    else
        band = "High SH (>10K)";
    end
    results.SH_band(k) = band;

    % --- additional flags (optional) ---
    flagFloodback = (SH_mean <= 0.5);  % stronger than <1K band (0.5K margin)
    flagStarved   = (md_mean < mdot_min);
    flagHunting   = (SH_std > SH_std_max) || (SH_pp > SH_pp_max);

    results.flagFloodback(k) = flagFloodback;
    results.flagStarved(k)   = flagStarved;
    results.flagHunting(k)   = flagHunting;

    % human-readable summary
    summary = "OK";
    if flagHunting,   summary = summary + "|Hunting";   end
    if flagStarved,   summary = summary + "|Starved";   end
    if flagFloodback, summary = summary + "|Floodback"; end
    results.flagSummary(k) = summary;
end

writetable(results,'outputs/logs/TXV_omega_sweep_SHbands.csv');

%% =========================
% Plot: SH bands shaded + points colored by band
%% =========================
x = results.mdot_mean * 7936.64;
y = results.SH_mean;

xmin = min(x); xmax = max(x);
pad  = 0.05*(xmax-xmin + eps);
xmin = xmin - pad; xmax = xmax + pad;

figure; hold on;

% Shaded regions (y-bands)
patch([xmin xmax xmax xmin], [0 0 SH_danger_hi SH_danger_hi], [1 0.85 0.85], ...
    'EdgeColor','none', 'FaceAlpha',0.35);   % Danger

patch([xmin xmax xmax xmin], [SH_sketch_lo SH_sketch_lo SH_sketch_hi SH_sketch_hi], [1 0.95 0.80], ...
    'EdgeColor','none', 'FaceAlpha',0.35);   % Sketchy

patch([xmin xmax xmax xmin], [SH_good_lo SH_good_lo SH_good_hi SH_good_hi], [0.85 1 0.85], ...
    'EdgeColor','none', 'FaceAlpha',0.35);   % Good

% Scatter, colored by band
gscatter(x, y, results.SH_band);
xlabel('Mass Flow [lb/hr]');
ylabel('Superheat [K]');
title('TXV Operating Regions – SH Bands (Good / Sketchy / Danger)');
xlim([xmin xmax]);
ylim([0 max(12, max(y)+1)]);
grid on;

% Band edges
yline(SH_danger_hi,'--','1 K');
yline(SH_good_lo,'--','4 K');
yline(SH_good_hi,'--','10 K');

saveas(gcf,'outputs/plots/TXV_mdot_SH_SHbands.png');

%% Optional: second plot highlighting "bad dynamics" flags
% (shows if a point is Good SH but still hunting/starved)
figure; hold on;
scatter(x, y, 60, 'filled');
idx = results.flagHunting | results.flagStarved | results.flagFloodback;
scatter(x(idx), y(idx), 120, 'o'); % circles flagged points
xlabel('Mass Flow [kg/s]');
ylabel('Superheat [K]');
title('TXV SH Map – circled points have Hunting/Starved/Floodback flags');
grid on;
saveas(gcf,'outputs/plots/TXV_mdot_SH_flags.png');
