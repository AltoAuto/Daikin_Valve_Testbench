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

% Save high-quality PNG
outPath = fullfile('outputs','plots','EXV_SH_vs_Opening.png');
print(fig, outPath, '-dpng', '-r300');