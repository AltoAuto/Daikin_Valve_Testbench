%% run_compare_EXV_TXV.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:
%   Run EXV and TXV simulations and overlay key signals (SH, mdot)
%   on poster-quality figures (16 x 5 in).
%
% MODELS:
%   two_phase_main_loop_EXV        -> EXV model
%   Copy_of_two_phase_main_loop_EXV -> TXV model (copy with TXV logic)
%
% OUTPUTS:
%   - /outputs/plots/EXV_vs_TXV_SH.png
%   - /outputs/plots/EXV_vs_TXV_mdot.png
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

%% Load config and models
run(fullfile('config','config.m'));

model_exv = 'D:\Daikin\simscape_two_phase_model\models\two_phase_main_loop.slx';
model_txv = 'D:\Daikin\simscape_two_phase_model\models\Copy_of_two_phase_main_loop.slx';

stopTime = '2000';  % adjust if needed

%% Run simulations
simOut_EXV = sim(model_exv, 'StopTime', stopTime);
simOut_TXV = sim(model_txv, 'StopTime', stopTime);

logs_EXV = simOut_EXV.logsout;
logs_TXV = simOut_TXV.logsout;

%% --- Extract signals: Superheat ---
SH_exv_sig = logs_EXV.getElement('SH').Values;
SH_txv_sig = logs_TXV.getElement('SH').Values;

t_exv  = SH_exv_sig.Time;
SH_exv = SH_exv_sig.Data;

t_txv  = SH_txv_sig.Time;
SH_txv = SH_txv_sig.Data;

%% --- Extract signals: Mass flow (assumes element 14 is mdot in both models) ---
md_exv_sig = logs_EXV.getElement(14).Values;
md_txv_sig = logs_TXV.getElement(14).Values;

t_m_exv  = md_exv_sig.Time;
mdot_exv = md_exv_sig.Data;    % kg/s

t_m_txv  = md_txv_sig.Time;
mdot_txv = md_txv_sig.Data;    % kg/s

%% FIGURE 1: EXV vs TXV Superheat (16 x 5 in, 300 dpi)
fig1 = figure('Units','inches','Position',[1 1 16 5],'PaperPositionMode','auto');
set(fig1,'Color','w');

ax1 = axes(fig1);
ax1.FontName  = 'Arial';
ax1.FontSize  = 16;
ax1.LineWidth = 1.2;
hold(ax1,'on');

p1 = plot(t_exv, SH_exv, 'LineWidth', 2.2, 'Color',[0 0.45 0.74]);   % EXV (blue)
p2 = plot(t_txv, SH_txv, 'LineWidth', 2.2, 'Color',[0.85 0.33 0.1]); % TXV (orange)

xlabel('Time [s]','FontWeight','bold');
ylabel('Superheat [K]','FontWeight','bold');
yl = yline(5.3, '--', 'SH Setpoint = 5.3 K', 'LineWidth', 2, 'Color', [0.2 0.2 0.2], 'FontSize', 14, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
title('EXV vs TXV — Superheat Response','FontWeight','bold');

grid on;
box on;
legend([p1 p2], {'EXV','TXV'}, 'Location','best', 'Box','off');

outPath1 = fullfile('outputs','plots','EXV_vs_TXV_SH.png');
print(fig1, outPath1, '-dpng', '-r300');

%% FIGURE 2: EXV vs TXV Mass Flow (16 x 5 in, 300 dpi)
fig2 = figure('Units','inches','Position',[1 1 16 5],'PaperPositionMode','auto');
set(fig2,'Color','w');

ax2 = axes(fig2);
ax2.FontName  = 'Arial';
ax2.FontSize  = 16;
ax2.LineWidth = 1.2;
hold(ax2,'on');

p3 = plot(t_m_exv, mdot_exv*3600*2.20462, 'LineWidth', 2.2, 'Color',[0 0.45 0.74]);   % EXV
p4 = plot(t_m_txv, mdot_txv*3600*2.20462, 'LineWidth', 2.2, 'Color',[0.85 0.33 0.1]); % TXV

xlabel('Time [s]','FontWeight','bold');
ylabel('Mass flow [lb/hr]','FontWeight','bold');
title('EXV vs TXV — Mass Flow Through Evaporator','FontWeight','bold');

box on;
legend([p3 p4], {'EXV','TXV'}, 'Location','best', 'Box','off');

outPath2 = fullfile('outputs','plots','EXV_vs_TXV_mdot.png');
print(fig2, outPath2, '-dpng', '-r300');

disp('Saved:');
disp(outPath1);
disp(outPath2);
