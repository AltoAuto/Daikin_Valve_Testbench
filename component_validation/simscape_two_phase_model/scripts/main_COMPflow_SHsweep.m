%% TXV RPM Sweep – Steady-State SH vs mdot Map (Parallel)
clc;
run(fullfile('config','config.m'));
cfg = evalin('base','cfg');

model = 'two_phase_main_loop';

% Make sure model is loaded
load_system(model);

% RPM sweep (fraction of nominal)
rpm_frac = linspace(0.3, 1.5, 10);   % 30% → 150%

nCases  = numel(rpm_frac);
simIn   = repmat(Simulink.SimulationInput(model), nCases, 1);

for k = 1:nCases
    frac  = rpm_frac(k);
    omega = frac * 250;  % your nominal speed scaling (rad/s or whatever you're using)

    simIn(k) = simIn(k) ...
        .setVariable('omega_in', omega) ...          % sets variable in model workspace for that worker
        .setModelParameter('StopTime','2000');         % 30 s sim
end

% Run all simulations in parallel
simOut = parsim(simIn, ...
    'ShowProgress','on', ...
    'TransferBaseWorkspaceVariables','on');   % sends cfg and other stuff to workers

% Preallocate results table
results = table([],[],[],[],[],[], ...
  'VariableNames',{'rpm_frac','omega_rad','mdot_kgps','SH_K','Psuct_Pa','Pevap_Pa'});

N = 50;    % last-N samples for averaging

for k = 1:nCases
    logsout = simOut(k).logsout;

    % Adjust these names to exactly match your signals:
    md_vec   = logsout.getElement(14).Values.Data;
    SH_vec = logsout.getElement('SH').Values.Data;

    n = @(v) min(N,numel(v));
    mdot = mean(md_vec(end-n(md_vec)+1:end));
    SH = mean(SH_vec(end-n(SH_vec)+1:end));

    frac  = rpm_frac(k);
    omega = frac * 250;

    results = [results; {frac, omega, mdot, SH, NaN, NaN}];
end

writetable(results,'outputs/logs/TXV_RPM_sweep.csv');

%% quick plot: mdot vs SH
figure;
% if you want lb/hr on x-axis:
mdot_lbhr = results.mdot_kgps;
scatter(mdot_lbhr, results.SH_K, 50, 'filled');
xlabel('Mass Flow [kg/s]');
ylabel('Superheat [K]');
title('TXV Steady-State Performance – ṁ vs SH (RPM sweep)');
grid on;
saveas(gcf,'outputs/plots/TXV_mdot_SH_map.png');
