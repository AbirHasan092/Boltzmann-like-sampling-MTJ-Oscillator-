%% LC_Oscillator_Parameter_Set_Generator.m
% Generate statistically varied parameter sets for the nonlinear LC
% state-space oscillator used in Code_oscillator_BSN_LC_state_space_Heun.m.
%
% Primitive quantities are sampled independently:
%   L            inductance [H]
%   C            capacitance [F]
%   R            parallel loss resistance [ohm]
%   S            nonlinear-current scale [A]
%   GnMagnitude  magnitude of the negative small-signal conductance [S]
%
% Derived quantities are recalculated for every device:
%   Gn     = -GnMagnitude
%   omega0 = 1/sqrt(L*C)
%   f0_Hz  = omega0/(2*pi)
%
% A sampled set is accepted only if |Gn| > 1/R, so the nonlinear element
% can overcome tank loss and sustain oscillation. Run this script before
% LC_Oscillator_Ensemble_Heun_and_Initialization_Sensitivity.m.

clear; clc; close all;

%% User settings
Nsweep = 50;
parameter_seed = 24680;

% Independent relative one-standard-deviation spreads. These are study
% assumptions; replace them with measured process statistics when known.
rel_sigma.L = 0.1;
rel_sigma.C = 0.1;
rel_sigma.R = 0.1;
rel_sigma.S = 0.10;
rel_sigma.GnMagnitude = 0.10;

% Require a finite negative-conductance margin above the tank loss.
minimum_negative_conductance_margin = 1.05;

script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir,'lc_oscillator_parameter_sets.mat');

%% Baseline from Lai-Roychowdhury LC state-space model, Eq. (21)
baseline = struct();
baseline.reference = ['Y. Lai and J. Roychowdhury, Capturing oscillator ' ...
    'injection locking via nonlinear phase-domain macromodels, ' ...
    'IEEE Transactions on Microwave Theory and Techniques.'];
baseline.L = 2.5e-8/(2*pi);       % [H]
baseline.C = 4e-11/(2*pi);       % [F]
baseline.R = 100;                % [ohm]
baseline.S = 1/baseline.R;       % [A]
baseline.GnMagnitude = 10/baseline.R; % [S]
baseline.Gn = -baseline.GnMagnitude;  % [S]
baseline.omega0 = 1/sqrt(baseline.L*baseline.C); % [rad/s]
baseline.f0_Hz = baseline.omega0/(2*pi);         % [Hz]
baseline.T0 = 1/baseline.f0_Hz;                  % [s]
baseline.frequency_offset_Hz = 0;
baseline.frequency_offset_fraction = 0;
baseline.negative_conductance_margin = ...
    baseline.GnMagnitude*baseline.R;

%% Generate independent Gaussian parameter sets
if Nsweep<1 || Nsweep~=round(Nsweep)
    error('Nsweep must be a positive integer.');
end
if minimum_negative_conductance_margin<=1
    error('minimum_negative_conductance_margin must exceed one.');
end

rng(parameter_seed,'twister');
parameter_sets = repmat(baseline,Nsweep,1);

for set_index = 1:Nsweep
    device = baseline;
    accepted = false;
    attempt = 0;

    while ~accepted
        attempt = attempt+1;
        if attempt>10000
            error('Could not generate a self-oscillating parameter set.');
        end

        device.L = positive_gaussian(baseline.L,rel_sigma.L);
        device.C = positive_gaussian(baseline.C,rel_sigma.C);
        device.R = positive_gaussian(baseline.R,rel_sigma.R);
        device.S = positive_gaussian(baseline.S,rel_sigma.S);
        device.GnMagnitude = positive_gaussian( ...
            baseline.GnMagnitude,rel_sigma.GnMagnitude);

        accepted = device.GnMagnitude > ...
            minimum_negative_conductance_margin/device.R;
    end

    device.Gn = -device.GnMagnitude;
    device.omega0 = 1/sqrt(device.L*device.C);
    device.f0_Hz = device.omega0/(2*pi);
    device.T0 = 1/device.f0_Hz;
    device.frequency_offset_Hz = device.f0_Hz-baseline.f0_Hz;
    device.frequency_offset_fraction = ...
        device.frequency_offset_Hz/baseline.f0_Hz;
    device.negative_conductance_margin = ...
        device.GnMagnitude*device.R;

    parameter_sets(set_index) = device;
end

sampled_parameter_fields = {'L','C','R','S','GnMagnitude'};
sampled_parameter_names = {'L','C','R','S','|G_n|'};
sampled_parameter_units = {'H','F','ohm','A','S'};

save(output_file,'parameter_sets','baseline','Nsweep','parameter_seed', ...
    'rel_sigma','sampled_parameter_fields','sampled_parameter_names', ...
    'sampled_parameter_units','minimum_negative_conductance_margin');

fprintf('Saved %d LC parameter sets to:\n%s\n',Nsweep,output_file);
fprintf('Baseline: L = %.6f nH, C = %.6f pF, R = %.3f ohm\n', ...
    baseline.L*1e9,baseline.C*1e12,baseline.R);
fprintf('Baseline: S = %.6g A, Gn = %.6g S, f0 = %.9f GHz\n', ...
    baseline.S,baseline.Gn,baseline.f0_Hz/1e9);
fprintf('Generated f0 range = %.6f to %.6f GHz\n', ...
    min([parameter_sets.f0_Hz])/1e9,max([parameter_sets.f0_Hz])/1e9);

%% Variation overview
parameter_values = [ ...
    [parameter_sets.L].', ...
    [parameter_sets.C].', ...
    [parameter_sets.R].', ...
    [parameter_sets.S].', ...
    [parameter_sets.GnMagnitude].'];
baseline_values = [baseline.L baseline.C baseline.R baseline.S ...
    baseline.GnMagnitude];
percent_deviation = 100.*(parameter_values-baseline_values) ...
    ./baseline_values;

figure('Color','w','Name','LC device parameter variation', ...
    'NumberTitle','off');
subplot(2,1,1);
boxplot(percent_deviation,'Labels',sampled_parameter_names,'Whisker',1.5);
ylabel('Deviation from baseline (%)');
title('Independently sampled primitive LC parameters');
set(gca,'FontSize',12,'LineWidth',1.3,'TickDir','in');
box on;

subplot(2,1,2);
histogram([parameter_sets.f0_Hz]/1e9,max(8,round(sqrt(Nsweep))), ...
    'FaceColor',[0.18 0.52 0.82]);
xline(baseline.f0_Hz/1e9,'r--','Baseline','LineWidth',1.5);
xlabel('Derived free-running frequency (GHz)');
ylabel('Device count');
title('Frequency variation produced by L and C');
set(gca,'FontSize',12,'LineWidth',1.3,'TickDir','in');
box on;

%% Local utility
function value = positive_gaussian(mean_value,relative_sigma)
% Draw N(mean,(relative_sigma*mean)^2), rejecting nonphysical values.
value = mean_value+abs(mean_value)*relative_sigma*randn;
while ~isfinite(value) || value<=0
    value = mean_value+abs(mean_value)*relative_sigma*randn;
end
end
