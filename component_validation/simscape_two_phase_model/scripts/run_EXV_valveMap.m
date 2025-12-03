%% Valve map mini-sweep  (Phase-2 characterization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: run_EXV_valveMap.m
% PURPOSE:
%   Perform a mini-sweep (Phase-2 characterization) to generate the valve
%   map: mass flow vs. opening and pressure drop across the EXV.
% OUTPUTS:
%   - /outputs/plots/EXV_valve_map.png
%   - /outputs/logs/valve_map.csv
% NOTES:
%   Compressor speed fixed; u_cmd stepped through test points.
%   Results are steady-state averages 
% AUTHOR: Aiden W.  |  Date: <10/29/2025>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
run(fullfile('config','config.m'));
cfg = evalin('base','cfg');

% Sweep definition
u_sweep = linspace(0.0003, 0.95, 20);
nCases  = numel(u_sweep);

%----------------------------------------------------------------------
% Build SimulationInput array for parsim (one SimInput per test point)
%----------------------------------------------------------------------
model = 'two_phase_main_loop';
load_system(model);

% Preallocate exactly nCases, clean
in = repmat(Simulink.SimulationInput(model), nCases, 1);

for k = 1:nCases
    in(k) = in(k).setVariable('u_cmd', u_sweep(k));
    in(k) = in(k).setModelParameter('StopTime','4000');
end

%----------------------------------------------------------------------
% Run simulations in parallel
%----------------------------------------------------------------------
simOut = parsim(in, ...
    'ShowProgress','on', ...
    'UseFastRestart','on');

%----------------------------------------------------------------------
% Post-process results
%----------------------------------------------------------------------
results = table('Size',[nCases 4], ...
    'VariableTypes', {'double','double','double','double'}, ...
    'VariableNames', {'u','A_m2','mdot_kgps','dP_Pa'});

N = 100;   % Last-N samples for averaging
nLast = @(v) min(N, numel(v));

for k = 1:nCases
    logsout = simOut(k).logsout;

    % extract signals
    md_vec   = logsout.getElement(14).Values.Data;   % md (kg/s)
    Pin_kPa  = logsout.getElement(11).Values.Data;   % upstream (kPa)
    Pout_kPa = logsout.getElement(3).Values.Data;    % downstream (kPa)

    % Last-N averaging
    idx_m  = nLast(md_vec);
    idx_Pi = nLast(Pin_kPa);
    idx_Po = nLast(Pout_kPa);

    mdot = mean(md_vec(end-idx_m+1:end));
    Pin  = mean(Pin_kPa(end-idx_Pi+1:end)) * 1e3;  % Pa
    Pout = mean(Pout_kPa(end-idx_Po+1:end)) * 1e3; % Pa
    dP   = Pin - Pout;

    u_cmd = u_sweep(k);
    A_m2  = cfg.exv.Amin_m2 + (cfg.exv.Amax_m2 - cfg.exv.Amin_m2)*sqrt(u_cmd);
    % Store row k
    results.u(k)        = u_cmd;
    results.A_m2(k)     = A_m2;
    results.mdot_kgps(k)= mdot;
    results.dP_Pa(k)    = dP;
end

% save results
writetable(results, fullfile('outputs','logs','valve_map.csv'));

%----------------------------------------------------------------------
% Poster-friendly plot (15 in x 6 in @ 300 dpi)
%----------------------------------------------------------------------

% Convert to intuitive units
mdot_lbhr = results.mdot_kgps * 3600 * 2.20462;   % kg/s -> lb/hr
dP_bar    = results.dP_Pa   / 1e2;               % Pa   -> bar

% Compute valve opening [%] from actual min/max area
A_min = cfg.exv.Amin_m2;
A_max = cfg.exv.Amax_m2;
open_pct = 100 * (results.A_m2 - A_min) / (A_max - A_min);

% Create figure
fig = figure('Units','inches','Position',[1 1 15 6],'PaperPositionMode','auto');
set(fig,'Color','w');   % white background

% Global styling
ax = axes(fig);
ax.FontName  = 'Arial';
ax.FontSize  = 16;
ax.LineWidth = 1.2;
hold(ax,'on');

yyaxis left
p1 = plot(open_pct, mdot_lbhr, '-', ...
    'LineWidth', 2, ...
    'MarkerSize', 6, ...
    'DisplayName','Mass flow');

ylabel('Mass flow [lb/hr]','FontWeight','bold');

yl_low = yline(80, ...
    'LineStyle','--', ...
    'LineWidth',3, ...                 % thicker for poster
    'Color',[0 0 0], ...               % solid black for visibility
    'Label','80 lb/hr (min controllable)', ...
    'FontSize',16, ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','bottom');


% Upper bound (1700 lb/hr)
yl_high = yline(1700, ...
    'LineStyle','--', ...
    'LineWidth',3, ...
    'Color',[0 0 0], ...
    'Label','1700 lb/hr (max controllable)', ...
    'FontSize',16, ...
    'LabelHorizontalAlignment','left', ...
    'LabelVerticalAlignment','top');

ylim([0, 1.1*max(mdot_lbhr)]);

yyaxis right
p2 = plot(open_pct, dP_bar, '--', ...
    'LineWidth', 2.0, ...
    'DisplayName','\DeltaP across EXV');

ylabel('\DeltaP [bar]','FontWeight','bold');
ylim([0, 1.2*max(dP_bar)]);

xlabel('Valve Opening [%]','FontWeight','bold');

title('EXV Flow Characterization — Mass Flow & \DeltaP vs Valve Opening', ...
    'FontWeight','bold');

box on;

% Legend (small, clean)
legend([p1 p2], 'Location','northwest','Box','on');

% Save high-quality PNG (15 x 6 inches)
print(fig, fullfile('outputs','plots','EXV_valve_map.png'), '-dpng', '-r300');
