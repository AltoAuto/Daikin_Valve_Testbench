%% main_EXV_openloop.m  —  Phase-2 demo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT: main_EXV_openloop.m
% PURPOSE:
%   Demonstrate dynamic EXV opening step (open-loop) and corresponding
%   superheat (SH) response in the two-phase refrigerant loop.
% OUTPUTS:
%   - /outputs/plots/EXV_SH_vs_Opening.png
%   - /outputs/logs/exv_openloop.mat, .csv
% AUTHOR: Aiden W.  |  Date: <10/29/2025>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run(fullfile('config','config.m'));
logsout = sim('two_phase_main_loop').logsout;       % run model

% Plot & save core results
t = logsout.getElement('SH').Values.Time;
SH = logsout.getElement('SH').Values.Data;
u  = logsout.getElement('A_m2').Values.Data;

A_min = cfg.exv.Amin_m2;
A_max = cfg.exv.Amax_m2;
open_pct = 100 * u  / A_max ;

fig = figure('Units','inches','Position',[1 1 15 6],'PaperPositionMode','auto');
set(fig,'Color','w');

ax = axes(fig);
ax.FontName  = 'Arial';
ax.FontSize  = 16;
ax.LineWidth = 1.2;
hold(ax,'on');

yyaxis left
plot(t, SH, 'LineWidth', 2.0);
ylabel('Superheat [K]','FontWeight','bold');

yyaxis right
plot(t, open_pct, '--','LineWidth', 2.0);
ylabel('Valve Opening [%]','FontWeight','bold');

xlabel('Time [s]','FontWeight','bold');

box on;
title('EXV Superheat Control — Step Response');

% Save high-quality PNG
outPath = fullfile('outputs','plots','EXV_SH_vs_Opening.png');

%% ===== Export logs to MAT + CSV =====
outLogsDir  = fullfile('outputs','logs');
outPlotsDir = fullfile('outputs','plots');
if ~exist(outLogsDir,'dir');  mkdir(outLogsDir);  end
if ~exist(outPlotsDir,'dir'); mkdir(outPlotsDir); end

getSig = @(name) logsout.getElement(name).Values;

% Base time vector = SH time
t = logsout.getElement('SH').Values.Time(:);

% SH on base time
SH = logsout.getElement('SH').Values.Data(:);

% EXV area A_m2 (interpolate to SH time if needed)
A_ts = logsout.getElement('A_m2').Values;
A_t  = A_ts.Time(:);
A    = A_ts.Data(:);
A_on_t = interp1(A_t, A, t, 'linear', 'extrap');

A_min = cfg.exv.Amin_m2;
A_max = cfg.exv.Amax_m2;

% Percent opening (use stroke-based % if Amin not zero)
open_pct = 100 * (A_on_t - A_min) ./ max(eps, (A_max - A_min));
open_pct = max(0, min(100, open_pct));  % clamp

% Optional: add mdot / pressures if you logged them (edit names as needed)
mdot_on_t = nan(size(t));
Psuct_on_t = nan(size(t));
Pdis_on_t  = nan(size(t));
names = logsout.getElementNames;
disp(names)
try
    md_ts = getSig(14);  % <-- change to your exact signal name
    mdot_on_t = interp1(md_ts.Time(:), md_ts.Data(:), t, 'linear', 'extrap');
end
try
    ps_ts = getSig(11);      % <-- change name if needed
    Psuct_on_t = interp1(ps_ts.Time(:), ps_ts.Data(:), t, 'linear', 'extrap');
end
try
    pd_ts = getSig(3);       % <-- change name if needed
    Pdis_on_t  = interp1(pd_ts.Time(:), pd_ts.Data(:), t, 'linear', 'extrap');
end

% Build table and write CSV
T = table(t, SH, A_on_t, open_pct, mdot_on_t, Psuct_on_t, Pdis_on_t, ...
    'VariableNames', {'t_s','SH_K','A_m2','open_pct','mdot_kgps','Psuct_MPa','Pdis_MPa'});

csvPath = fullfile(outLogsDir, 'exv_openloop.csv');
matPath = fullfile(outLogsDir, 'exv_openloop.mat');

writetable(T, csvPath);
save(matPath, 'T', 'cfg');

fprintf("Saved:\n  %s\n  %s\n", csvPath, matPath);

Npts = 20;
idx = round(linspace(1, numel(t), Npts));

t20  = t(idx);
SH20 = SH(idx);
A20  = A_on_t(idx);
open20 = open_pct(idx);
mdot_on_t20 =mdot_on_t(idx);
Psuct_on_t20=Psuct_on_t(idx);
Pdis_on_t20=Pdis_on_t(idx);

T20 = table(t20, SH20, A20, open20,mdot_on_t20,Psuct_on_t20,Pdis_on_t20, ...
    'VariableNames', {'t_s','SH_K','A_m2','open_pct','mdot_kgps','Psuct_MPa','Pdis_MPa'});

writetable(T20, fullfile('outputs','logs','exv_openloop_20pts.csv'));
