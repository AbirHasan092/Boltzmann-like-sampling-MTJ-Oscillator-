%% Ensemble_simulation_and_IC_sensitivity_parametric.m
% Full nonlinear LC state-space Monte Carlo experiment with:
%   1) statistically varied LC devices,
%   2) stochastic-Heun integration of [inductor current; tank voltage],
%   3) a constant-amplitude FHI lock at the SHI-basin midpoint, removed
%      before the gamma-dependent FHI and parametric SHI sampling pulse,
%   4) independent SHI-only two-state calibration; fixed axis per device,
%      with every measured phase and binary outcome retained, and
%   5) FHI-phase sensitivity on one representative device.
% This is a full-state LC implementation of the paper's single-oscillator
% sampling sequence, not its reduced phase equation or coupled-OIM study.
% The measured probability, zero-input bias, and symmetry diagnostics must
% be checked before interpreting a tanh fit as Boltzmann-like sampling.
%
% Run LC_Oscillator_Parameter_Set_Generator.m first. Results are
% checkpointed after every device set and FHI-phase sensitivity test.

function Ensemble_simulation_and_IC_sensitivity_parametric( ...
    reference_frequency_mode,results_filename_override,fhi_start_time_s,config_overrides)
% Optional function arguments select the frequency-reference experiment.
% Calling this function without arguments selects
% measured-frequency behavior and a NEW parametric-model output filename.
if nargin<1 || isempty(reference_frequency_mode)
    reference_frequency_mode = 'measured_limit_cycle';
end
if nargin<2 || isempty(results_filename_override)
    results_filename_override = '';
end
if nargin<3 || isempty(fhi_start_time_s)
    if strcmp(reference_frequency_mode,'measured_limit_cycle')
        fhi_start_time_s = 100e-9;
    else
        fhi_start_time_s = 0;
    end
end
if nargin<4, config_overrides = struct(); end
clc; close all;

%% Files and run controls
script_dir = fileparts(mfilename('fullpath'));
parameter_file = fullfile(script_dir,'lc_oscillator_parameter_sets.mat');
if isempty(results_filename_override)
    results_file = fullfile(script_dir, ...
        '175n_lc_oscillator_2FHI_paramSHI.mat');
else
    results_file = fullfile(script_dir,results_filename_override);
end

% Keep all additive-SHI datasets intact. Parametric-SHI datasets receive a
% distinct suffix and cannot be resumed as if they came from the old model.
[result_dir,result_stem,result_ext] = fileparts(results_file);
if ~endsWith(result_stem,'_parametric_two_state')
    result_stem = [result_stem '_parametric_two_state'];
end
results_file = fullfile(result_dir,[result_stem result_ext]);

requested_number_of_parameter_sets = 50;
requested_parameter_set_indices = []; % Empty means all requested sets.
representative_parameter_set = 1;
resume_existing_results = true;
base_seed = 120;

%% Probability experiment
config = struct();
config.gamma_values = linspace(-1.5,1.5,51);
config.An = 0.3;                  % normalized noise multiplier
config.j_max = 300;              % trajectories per gamma point
config.As = 1;                   % normalized SHI peak envelope

% Locking FHI is constant across gamma; synaptic FHI alone is swept.
config.IlockFHI_A = 0.001;        % [A], nominal lock current; requires longer lock
config.IfhiUnit = 0.005;        % [A] per unit gamma, synaptic FHI only
% NEW PHYSICAL MODEL: time-varying conductance, NOT additive SHI current.
config.GshiUnit = 0.002;          % [S]; starting value, validated per device
config.shiCoupling = 'parametric_conductance_v1';
config.requireTwoStateValidation = true;
config.calibrationInitialPhases_deg = 0:15:345;
config.calibrationSettleCycles = 500;
config.calibrationClusterTolerance_deg = 3;
config.calibrationAntipodalTolerance_deg = 3;
config.calibrationDriftTolerance_deg = 1;
config.minimumFundamentalFraction = 0.05;
config.maximumInvalidPhaseFraction = 0.01;
% The preparation FHI ALWAYS turns off before sampling. This option only
% removes the synaptic FHI after SHI reaches full strength for readout.
config.removeFhiBeforeReadout = true;
config.lockFhiRampOffTime_s = 20e-9;
config.fhiRampOffTime_s = 20e-9; % synaptic FHI removal, after SHI is full
% The preparation voltage phase is measured relative to the calibrated
% bit-1 SHI state. A deterministic protocol-specific basin-boundary test
% refines the nominal +90-deg midpoint for each varied device.
% The synaptic FHI is in quadrature with the preparation FHI, as in the
% paper's phase model, rather than a second copy of the lock waveform.
config.phiFHI_deg = 90;
config.synapticPhaseOffset_deg = -90;
config.lockPhaseDefinition = 'measured_SHI_basin_midpoint_v1';
config.injectionProtocol = 'paper_FHI0_off_before_synaptic_FHI_and_SHI_v1';
config.slopeConvention = 'signed_synaptic_FHI_quadrature_v1';
config.phaseReferenceConvention = 'fixed_zero_phase';
config.bitAssignmentConvention = ...
    'independently_calibrated_antipodal_bit1_zero_fhi90_v2';
config.bit0Phase_deg = 180;      % angle in CALIBRATED coordinates
config.bit1Phase_deg = 0;
config.bitBoundaryRule = 'abs(wrapped_phase)<=90; exact ties are bit 1';
config.boundaryWarningWidth_deg = 5; % diagnostic only; no samples removed
config.referenceFrequencyMode = reference_frequency_mode;
config.phiSHI_deg = 0;
config.noiseCurrentBase = 175e-9; % [A/sqrt(Hz)] per unit An
config.fhiStartTime_s = fhi_start_time_s; % [s], FHI/reference start time t1

% Sequential injection schedule. The oscillator evolves autonomously until
% t1=fhiStartTime_s, then receives only the locking FHI. That current
% ramps fully to zero before the synaptic FHI and SHI ramps begin together.
% After SHI reaches full strength, the synaptic FHI may be removed for an
% SHI-only readout. Cycle counts refer to the selected reference frequency.
config.fhiLockCycles = 250;
config.synapticFhiRampTime_s = 2e-9; % fast relative to the SHI ramp
config.shiRampTime_s = 20e-9;     % [s], user-controlled SHI turn-on ramp
config.combinedSettleCycles = 100;
config.phaseWindowCycles = 20;
config.lockPhaseWindowCycles = 20;
config.lockPhaseTolerance_deg = 5;
config.lockPhaseDriftTolerance_deg = 2;
config.lockPhaseMaxIterations = 5;
config.lockInitialPhases_deg = 0:45:315;
config.lockMinimumAmplitudeFraction = 0.1;
config.autoCenterNeutralPhase = true;
config.neutralBisectionIterations = 9;

% Numerical resolution. The selected reference period defines the injection
% waveforms and integration time step.
config.stepsPerPeriod = 500;
config.startupCycles = 100;
config.saveEverySteps = 20;
config.initialVoltage = 1e-3;      % [V], initial tank-voltage seed

% Nonlinear free-running-frequency measurement used only by the strict
% no-detuning mode. Calibration is deterministic and precedes every
% stochastic probability curve.
config.frequencyMinSettleCycles = 100;
config.frequencyMaxSettleCycles = 500;
config.frequencyConvergencePeriods = 25;
config.frequencyConvergenceTolerance = 1e-3;

% Common oscillator starting phase for every production and sensitivity run.
config.parameter_initial_phase_deg = 0;

% Initialization-phase sensitivity on the representative device: sweep the
% measured preparation phase from 0 to 180 deg relative to the bit-1 SHI
% state. The synaptic FHI remains in quadrature with its drive. The initial state and the
% 0-deg reference remain fixed. The same random
% numbers are used for every row when the flag is true, isolating FHI-phase
% sensitivity from finite-noise differences.
config.phiFHI_sensitivity_deg =0:10:180;
config.representativeGammaTargets = [-0.5 -0.25 0 0.25 0.5];
config.use_common_noise_for_phiFHI_sensitivity = true;

% Set false after debugging if representative waveforms are unnecessary.
config.store_representative_trace = true;
config.representativeGamma = 0.48; % target gamma for the saved waveform

% Optional small-run overrides; structural/physical conventions stay explicit.
if ~isstruct(config_overrides) || ~isscalar(config_overrides)
    error('config_overrides must be a scalar struct.');
end
overrideNames = fieldnames(config_overrides);
for q = 1:numel(overrideNames)
    name = overrideNames{q};
    if ~isfield(config,name)
        error('Unknown config override: %s.',name);
    end
    config.(name) = config_overrides.(name);
end
% A changed gamma sweep should not leave the saved-trace or sensitivity
% targets outside its range unless the caller explicitly chose them.
if isfield(config_overrides,'gamma_values')
    if ~isfield(config_overrides,'representativeGamma')
        config.representativeGamma = config.gamma_values(end);
    end
    if ~isfield(config_overrides,'representativeGammaTargets') && ...
            numel(config.gamma_values)>=5
        targetIndices = round(linspace(1,numel(config.gamma_values),5));
        config.representativeGammaTargets = config.gamma_values(targetIndices);
    end
end

%% Load and validate parameter sets
if ~isfile(parameter_file)
    error(['Missing %s. Run LC_Oscillator_Parameter_Set_Generator.m ' ...
        'first, with both scripts in the same folder.'],parameter_file);
end

source = load(parameter_file);
required_parameter_fields = {'parameter_sets','baseline','rel_sigma', ...
    'sampled_parameter_fields','sampled_parameter_names', ...
    'sampled_parameter_units','parameter_seed'};
for field_index = 1:numel(required_parameter_fields)
    if ~isfield(source,required_parameter_fields{field_index})
        error('Parameter file is missing %s.', ...
            required_parameter_fields{field_index});
    end
end

parameter_sets = source.parameter_sets;
baseline = source.baseline;
rel_sigma = source.rel_sigma;
sampled_parameter_fields = source.sampled_parameter_fields;
sampled_parameter_names = source.sampled_parameter_names;
sampled_parameter_units = source.sampled_parameter_units;
parameter_seed = source.parameter_seed;

number_of_parameter_sets = min(requested_number_of_parameter_sets, ...
    numel(parameter_sets));
if isempty(requested_parameter_set_indices)
    parameter_set_indices = 1:number_of_parameter_sets;
else
    parameter_set_indices = requested_parameter_set_indices(:).';
end
if any(parameter_set_indices<1) || ...
        any(parameter_set_indices>number_of_parameter_sets) || ...
        any(parameter_set_indices~=round(parameter_set_indices))
    error('requested_parameter_set_indices contains an invalid index.');
end
if representative_parameter_set<1 || ...
        representative_parameter_set>number_of_parameter_sets
    error('representative_parameter_set is outside the requested ensemble.');
end

validate_config(config);
config.parameter_set_indices = parameter_set_indices;
config.number_of_parameter_sets = number_of_parameter_sets;
config.representative_parameter_set = representative_parameter_set;
config.base_seed = base_seed;
config.parameter_seed = parameter_seed;
config.baseline_signature = [baseline.L baseline.C baseline.R ...
    baseline.S baseline.Gn baseline.omega0];
config.parameter_signature = [[parameter_sets(1:number_of_parameter_sets).L].' ...
    [parameter_sets(1:number_of_parameter_sets).C].' ...
    [parameter_sets(1:number_of_parameter_sets).R].' ...
    [parameter_sets(1:number_of_parameter_sets).S].' ...
    [parameter_sets(1:number_of_parameter_sets).Gn].'];

gamma_values = config.gamma_values(:).';
Ngamma = numel(gamma_values);
NphiFHI = numel(config.phiFHI_sensitivity_deg);

%% Allocate restartable results
Pbit_all = nan(number_of_parameter_sets,Ngamma);
phase_axis_all_deg = nan(number_of_parameter_sets,1);
reference_frequency_all_Hz = nan(number_of_parameter_sets,1);
completed_parameter_sets = false(number_of_parameter_sets,1);

Pbit_phiFHI = nan(NphiFHI,Ngamma);
phase_axis_phiFHI_deg = nan(NphiFHI,1);
completed_phiFHI_conditions = false(NphiFHI,1);

representative_trace = empty_trace();
representative_phase_samples_deg = [];
phase_audit = struct('device',{cell(number_of_parameter_sets,1)}, ...
    'phiFHI',{cell(NphiFHI,1)});

if resume_existing_results && isfile(results_file)
    old = load(results_file);
    if compatible_checkpoint(old,config,number_of_parameter_sets, ...
            Ngamma,NphiFHI)
        Pbit_all = old.Pbit_all;
        phase_audit = old.phase_audit;
        phase_axis_all_deg = old.phase_axis_all_deg;
        reference_frequency_all_Hz = old.reference_frequency_all_Hz;
        completed_parameter_sets = old.completed_parameter_sets;
        Pbit_phiFHI = old.Pbit_phiFHI;
        phase_axis_phiFHI_deg = old.phase_axis_phiFHI_deg;
        completed_phiFHI_conditions = old.completed_phiFHI_conditions;
        if isfield(old,'representative_trace')
            representative_trace = old.representative_trace;
        end
        if isfield(old,'representative_phase_samples_deg')
            representative_phase_samples_deg = ...
                old.representative_phase_samples_deg;
        end
        fprintf('Resuming compatible results from:\n%s\n',results_file);
    else
        warning('Existing results are incompatible; starting a new run.');
    end
end

%% Device-parameter ensemble
fprintf('\n============================================================\n');
fprintf(' LC DEVICE VARIATION: FULL STATE-SPACE STOCHASTIC HEUN\n');
fprintf('============================================================\n');
if strcmp(config.referenceFrequencyMode,'measured_limit_cycle')
    fprintf(['Injection reference       = measured nonlinear limit-cycle ' ...
        'frequency (no detuning)\n']);
else
    fprintf('Fixed injection reference = %.9f GHz\n',baseline.f0_Hz/1e9);
end
fprintf('Parameter sets requested  = %d\n',number_of_parameter_sets);
fprintf('Gamma points              = %d\n',Ngamma);
fprintf('Trajectories per point    = %d\n',config.j_max);
fprintf('FHI phases tested         = %d\n',NphiFHI);
fprintf('FHI/reference start t1    = %.3f ns\n', ...
    config.fhiStartTime_s*1e9);
fprintf('FHI-only locking          = %d reference cycles\n', ...
    config.fhiLockCycles);
fprintf('Preparation FHI ramp-off  = %.3f ns (before sampling)\n', ...
    config.lockFhiRampOffTime_s*1e9);
fprintf('Synaptic FHI turn-on ramp = %.3f ns\n', ...
    config.synapticFhiRampTime_s*1e9);
fprintf('SHI turn-on ramp          = %.3f ns (synaptic FHI only)\n', ...
    config.shiRampTime_s*1e9);
fprintf('Readout settling          = %d reference cycles\n', ...
    config.combinedSettleCycles);
fprintf('Phase samples/trajectory  = 1\n');
fprintf('Parametric SHI amplitude  = %.3f mS (NOT mA)\n',config.GshiUnit*1e3);
fprintf('Lock/synaptic FHI scale   = %.3f / %.3f mA per gamma\n', ...
    config.IlockFHI_A*1e3,config.IfhiUnit*1e3);
fprintf('Nominal midpoint offset  = %.1f deg from bit 1\n',config.phiFHI_deg);
fprintf('Synaptic/lock drive offset = %+.1f deg\n', ...
    config.synapticPhaseOffset_deg);
fprintf('FHI drive phase           = calibrated for each device\n');
fprintf('Probability convention    = bit 1 near calibrated SHI-only 0 deg\n');
fprintf('FHI removed for readout   = %d\n',config.removeFhiBeforeReadout);
fprintf('============================================================\n');

for list_index = 1:numel(parameter_set_indices)
    set_index = parameter_set_indices(list_index);
    if completed_parameter_sets(set_index)
        fprintf('Parameter set %d already complete; skipping.\n',set_index);
        continue
    end

    device = parameter_sets(set_index);
    fprintf('\n========== LC parameter set %d of %d ==========\n', ...
        set_index,number_of_parameter_sets);
    fprintf(['L=%.6f nH, C=%.6f pF, R=%.3f ohm, ' ...
        'S=%.5g A, Gn=%.5g S, f0=%.6f GHz\n'], ...
        device.L*1e9,device.C*1e12,device.R,device.S,device.Gn, ...
        device.f0_Hz/1e9);
    if strcmp(config.referenceFrequencyMode,'measured_limit_cycle')
        fprintf(['Analytical tank frequency = %.6f GHz; nonlinear frequency ' ...
            'will be measured before injection.\n'],device.f0_Hz/1e9);
    else
        fprintf('Detuning from fixed reference = %+.4f %%\n', ...
            100*(device.f0_Hz/baseline.f0_Hz-1));
    end

    store_trace = config.store_representative_trace && ...
        set_index==representative_parameter_set;
    curve_seed = base_seed+100000*(set_index-1);
    [Pbit_all(set_index,:),phase_axis_all_deg(set_index), ...
        phase_samples_deg,candidate_trace, ...
        reference_frequency_all_Hz(set_index),curve_audit] = run_probability_curve( ...
        device,baseline,gamma_values, ...
        config.parameter_initial_phase_deg,config,curve_seed,store_trace);

    phase_audit.device{set_index} = curve_audit;
    if store_trace
        representative_trace = candidate_trace;
        representative_trace.parameter_set_index = set_index;
        representative_phase_samples_deg = phase_samples_deg;
    end
    completed_parameter_sets(set_index) = true;

    save(results_file,'Pbit_all','phase_axis_all_deg', ...
        'completed_parameter_sets','reference_frequency_all_Hz', ...
        'Pbit_phiFHI', ...
        'phase_axis_phiFHI_deg','completed_phiFHI_conditions', ...
        'gamma_values','parameter_sets','baseline','rel_sigma', ...
        'sampled_parameter_fields','sampled_parameter_names', ...
        'sampled_parameter_units','parameter_seed', ...
        'number_of_parameter_sets','parameter_set_indices', ...
        'representative_parameter_set','representative_trace', ...
        'representative_phase_samples_deg','phase_audit','config','-v7.3');
    fprintf('Checkpoint saved after parameter set %d.\n',set_index);
end

%% FHI-phase sensitivity on the representative device
if ~completed_parameter_sets(representative_parameter_set)
    error(['The representative parameter set is incomplete. Include it in ' ...
        'requested_parameter_set_indices before the sensitivity test.']);
end

device = parameter_sets(representative_parameter_set);
for phiFHI_index = 1:NphiFHI
    if completed_phiFHI_conditions(phiFHI_index)
        fprintf('FHI phase %.1f deg already complete; skipping.\n', ...
            config.phiFHI_sensitivity_deg(phiFHI_index));
        continue
    end

    if config.use_common_noise_for_phiFHI_sensitivity
        phiFHI_seed = base_seed+700000;
    else
        phiFHI_seed = base_seed+700000+100000*(phiFHI_index-1);
    end

    sensitivity_config = config;
    sensitivity_config.precomputedCalibration = ...
        phase_audit.device{representative_parameter_set}.two_state_calibration;
    sensitivity_config.autoCenterNeutralPhase = false;
    sensitivity_config.phiFHI_deg = ...
        config.phiFHI_sensitivity_deg(phiFHI_index);

    fprintf('\nFHI-phase sensitivity test %d/%d: phiFHI = %.1f deg\n', ...
        phiFHI_index,NphiFHI,sensitivity_config.phiFHI_deg);
    [Pbit_phiFHI(phiFHI_index,:), ...
        phase_axis_phiFHI_deg(phiFHI_index),~,~,~,curve_audit] = ...
        run_probability_curve(device,baseline,gamma_values, ...
        config.parameter_initial_phase_deg,sensitivity_config, ...
        phiFHI_seed,false);
    phase_audit.phiFHI{phiFHI_index} = curve_audit;
    completed_phiFHI_conditions(phiFHI_index) = true;

    save(results_file,'Pbit_all','phase_axis_all_deg', ...
        'completed_parameter_sets','reference_frequency_all_Hz', ...
        'Pbit_phiFHI', ...
        'phase_axis_phiFHI_deg','completed_phiFHI_conditions', ...
        'gamma_values','parameter_sets','baseline','rel_sigma', ...
        'sampled_parameter_fields','sampled_parameter_names', ...
        'sampled_parameter_units','parameter_seed', ...
        'number_of_parameter_sets','parameter_set_indices', ...
        'representative_parameter_set','representative_trace', ...
        'representative_phase_samples_deg','phase_audit','config','-v7.3');
end

%% Final metadata and compact diagnostic plots
nBits = 1;
nProbabilitySamples = config.j_max;
integration_method = ...
    ['stochastic Heun on full nonlinear LC states [i;v], with ' ...
    'multiplicative parametric conductance SHI'];
shi_coupling_equation = ...
    'I_SHI(t)=GshiUnit*a(t)*sin(2*omega_ref*(t-t1)+phiSHI)*v(t)';
phase_definition = ...
    ['windowed phase(v)-phase(selected zero-phase reference), ' ...
    'with a calibrated SHI-only bit axis and deterministic neutral ' ...
    'separatrix for each device'];
neutral_centering_method = ...
    'deterministic zero-input FHI0-release/SHI-ramp basin bisection';
sampling_sequence = ...
    ['autonomous evolution to t1 -> first FHI locks voltage to the ' ...
    'calibrated SHI-basin midpoint -> first FHI ramps fully off -> ' ...
    'quadrature synaptic FHI and parametric SHI begin together, with ' ...
    'synaptic FHI rising faster -> optional synaptic FHI removal after ' ...
    'SHI reaches full strength -> phase readout'];
if strcmp(config.referenceFrequencyMode,'measured_limit_cycle')
    variation_interpretation = ['The reference, FHI, SHI, phase kernel, ' ...
        'and integration period use each device''s measured nonlinear ' ...
        'limit-cycle frequency; frequency detuning is excluded.'];
else
    variation_interpretation = ['L and C variation produces physical ' ...
        'detuning because every device sees the same nominal reference.'];
end

save(results_file,'Pbit_all','phase_axis_all_deg', ...
    'completed_parameter_sets','reference_frequency_all_Hz', ...
    'Pbit_phiFHI', ...
    'phase_axis_phiFHI_deg','completed_phiFHI_conditions', ...
    'gamma_values','parameter_sets','baseline','rel_sigma', ...
    'sampled_parameter_fields','sampled_parameter_names', ...
    'sampled_parameter_units','parameter_seed', ...
    'number_of_parameter_sets','parameter_set_indices', ...
    'representative_parameter_set','representative_trace', ...
    'representative_phase_samples_deg','phase_audit','config','nBits', ...
    'nProbabilitySamples','integration_method','phase_definition', ...
    'sampling_sequence','variation_interpretation', ...
    'neutral_centering_method', ...
    'shi_coupling_equation','-v7.3');

complete_set_numbers = find(completed_parameter_sets & ...
    all(isfinite(Pbit_all),2));
if ~isempty(complete_set_numbers)
    figure('Color','w','Name','LC device-variation probability curves', ...
        'NumberTitle','off');
    hold on;
    colors = lines(numel(complete_set_numbers));
    for row = 1:numel(complete_set_numbers)
        set_index = complete_set_numbers(row);
        plot(gamma_values,Pbit_all(set_index,:),'o-', ...
            'Color',colors(row,:),'LineWidth',1.4,'MarkerSize',4, ...
            'DisplayName',sprintf('Set %d',set_index));
    end
    yline(0.5,'k:','P(1)=0.5');
    xlabel('\gamma'); ylabel('p = P(bit=1)');
    title('Varied nonlinear LC state-space devices');
    legend('Location','eastoutside');
    ylim([0 1]); grid on; box on;
end

if all(completed_phiFHI_conditions)
    figure('Color','w','Name','LC FHI-phase sensitivity', ...
        'NumberTitle','off');
    hold on;
    targets = config.representativeGammaTargets;
    selectedIndices = zeros(size(targets));
    for targetIndex = 1:numel(targets)
        [~,selectedIndices(targetIndex)] = min(abs( ...
            gamma_values-targets(targetIndex)));
    end
    colors = lines(numel(targets));
    for targetIndex = 1:numel(targets)
        gammaIndex = selectedIndices(targetIndex);
        plot(config.phiFHI_sensitivity_deg, ...
            Pbit_phiFHI(:,gammaIndex),'o-', ...
            'Color',colors(targetIndex,:),'LineWidth',1.4, ...
            'MarkerSize',4,'DisplayName',sprintf( ...
            '\\gamma=%+.3g',gamma_values(gammaIndex)));
    end
    yline(0.5,'k:');
    neutralAudit = phase_audit.device{representative_parameter_set};
    if isstruct(neutralAudit) && ...
            isfield(neutralAudit,'neutral_calibration') && ...
            neutralAudit.neutral_calibration.performed
        xline(neutralAudit.neutral_calibration.offset_deg,'k--', ...
            'Calibrated neutral phase','LabelVerticalAlignment','bottom');
    end
    xlabel('Measured initialization phase relative to bit 1 (deg)');
    ylabel('p = P(bit=1)');
    title(sprintf('Initialization-phase sensitivity, device set %d', ...
        representative_parameter_set));
    legend('Location','best');
    if numel(config.phiFHI_sensitivity_deg)>1
        xlim([min(config.phiFHI_sensitivity_deg) ...
            max(config.phiFHI_sensitivity_deg)]);
    end
    ylim([0 1]); grid on; box on;
end

if completed_parameter_sets(representative_parameter_set) && ...
        ~isempty(phase_audit.device{representative_parameter_set})
    calibrationPlot = phase_audit.device{representative_parameter_set} ...
        .two_state_calibration;
    figure('Color','w','Name','Independent SHI-only two-state check', ...
        'NumberTitle','off');
    plot(calibrationPlot.initial_phases_deg, ...
        calibrationPlot.final_phases_deg,'ko','MarkerFaceColor',[0.2 0.5 0.9]);
    hold on;
    yline(calibrationPlot.cluster_centers_deg(1),'b--','Bit 0');
    yline(calibrationPlot.cluster_centers_deg(2),'r--','Bit 1');
    xlabel('Initial free-running phase (degree)');
    ylabel('Final SHI-only phase (degree)');
    title(sprintf(['Validated parametric-SHI states: separation %.3f^\\circ, ' ...
        'G_{SHI}=%.3f mS'],calibrationPlot.cluster_separation_deg, ...
        1e3*calibrationPlot.GshiUnit_S));
    xlim([min(calibrationPlot.initial_phases_deg) ...
        max(calibrationPlot.initial_phases_deg)]);
    ylim([-180 180]); yticks(-180:45:180); grid on; box on;
end

fprintf('\nSaved LC ensemble and FHI-phase sensitivity results to:\n%s\n', ...
    results_file);

end

%% Local simulation functions
function [Pbit,phaseAxisDeg,phaseSamplesDeg,trace,referenceFrequencyHz, ...
    audit] = ...
    run_probability_curve(device,baseline,gammaValues,initialPhaseDeg, ...
    config,curveSeed,storeTrace)

if strcmp(config.referenceFrequencyMode,'measured_limit_cycle')
    [referenceOmega,referencePeriod,initialCurrent,initialVoltage, ...
        limitCycleInfo] = characterize_limit_cycle(device,initialPhaseDeg, ...
        config);
    fprintf(['  measured nonlinear reference = %.9f GHz ' ...
        '(period CV %.3g, amplitude CV %.3g)\n'], ...
        referenceOmega/(2*pi*1e9),limitCycleInfo.period_cv, ...
        limitCycleInfo.amplitude_cv);
else
    referenceOmega = baseline.omega0;
    referencePeriod = baseline.T0;
    [initialCurrent,initialVoltage] = initial_state_on_limit_cycle( ...
        device,initialPhaseDeg,config.stepsPerPeriod,config.startupCycles, ...
        config.initialVoltage);
    limitCycleInfo = struct('period_cv',nan,'amplitude_cv',nan, ...
        'settling_cycles',config.startupCycles);
end
referenceFrequencyHz = referenceOmega/(2*pi);
if isfield(config,'precomputedCalibration')
    calibration = config.precomputedCalibration;
else
    calibration = validate_parametric_states(device,referenceOmega, ...
        referencePeriod,config);
end
if ~calibration.passed
    error('SHI-only calibration failed: %s',calibration.failure_reason);
end
phaseAxisDeg = calibration.bit1_axis_deg;
fprintf('  validated SHI-only bit-1/bit-0 states: %.3f / %.3f deg; min A1 %.4g V\n', ...
    phaseAxisDeg,wrap_to_180(phaseAxisDeg+180), ...
    min(calibration.fundamental_amplitude_V));
if ~isfield(limitCycleInfo,'orbit_current_A')
    [~,~,~,~,lockOrbitInfo] = characterize_limit_cycle( ...
        device,initialPhaseDeg,config);
else
    lockOrbitInfo = limitCycleInfo;
end
nominalMidpointPhaseDeg = wrap_to_180(phaseAxisDeg+90);
nominalLockCalibration = calibrate_fhi_lock(device,referenceOmega, ...
    referencePeriod,lockOrbitInfo,nominalMidpointPhaseDeg,config);
% Keep the synaptic drive fixed relative to the SHI bit axis, even when
% the preparation phase is swept or its true separatrix is not 90 deg.
synapticDrivePhaseDeg = wrap_to_180( ...
    nominalLockCalibration.drive_phase_deg+config.synapticPhaseOffset_deg);
if config.autoCenterNeutralPhase
    neutralCalibration = find_neutral_separatrix(device,referenceOmega, ...
        referencePeriod,initialCurrent,initialVoltage,phaseAxisDeg, ...
        nominalLockCalibration.drive_phase_deg,config);
    preparationOffsetDeg = neutralCalibration.offset_deg;
else
    preparationOffsetDeg = config.phiFHI_deg;
    neutralCalibration = struct('performed',false, ...
        'offset_deg',preparationOffsetDeg,'bracket_deg',[], ...
        'method','requested_phase_offset_without_separatrix_search');
end
lockTargetPhaseDeg = wrap_to_180(phaseAxisDeg+preparationOffsetDeg);
if abs(wrap_to_180(lockTargetPhaseDeg-nominalMidpointPhaseDeg))<1e-12
    lockCalibration = nominalLockCalibration;
else
    lockCalibration = calibrate_fhi_lock(device,referenceOmega, ...
        referencePeriod,lockOrbitInfo,lockTargetPhaseDeg,config);
end
fhiDrivePhaseDeg = lockCalibration.drive_phase_deg;
fprintf(['  FHI lock target/measured/drive = %.2f / %.2f / %.2f deg; ' ...
    'minimum A1 %.4g V\n'],lockTargetPhaseDeg, ...
    lockCalibration.measured_phase_deg,fhiDrivePhaseDeg, ...
    min(lockCalibration.fundamental_amplitude_V));
fprintf(['  SHI bit-1 axis / preparation offset / synaptic drive = ' ...
    '%.2f / %.2f / %.2f deg\n'], ...
    phaseAxisDeg,preparationOffsetDeg,synapticDrivePhaseDeg);
dt = referencePeriod/config.stepsPerPeriod;
fhiStartTime = config.fhiStartTime_s;
lockOnlyEndTime = fhiStartTime+config.fhiLockCycles*referencePeriod;
lockRampEndTime = lockOnlyEndTime+config.lockFhiRampOffTime_s;
synapticStartTime = lockRampEndTime;
synapticRampEndTime = synapticStartTime+config.synapticFhiRampTime_s;
shiRampEndTime = synapticStartTime+config.shiRampTime_s;
measurementTime = shiRampEndTime ...
    +double(config.removeFhiBeforeReadout)*config.fhiRampOffTime_s ...
    +config.combinedSettleCycles*referencePeriod;
phaseWindowStartTime = measurementTime ...
    -config.phaseWindowCycles*referencePeriod;
lockWindowStartTime = lockOnlyEndTime ...
    -config.lockPhaseWindowCycles*referencePeriod;
physicalStopTime = measurementTime;
nSteps = ceil(physicalStopTime/dt);
nBits = 1;

Ngamma = numel(gammaValues);
stateSize = [Ngamma config.j_max];
iState = initialCurrent*ones(stateSize);
vState = initialVoltage*ones(stateSize);
sigmaVoltage = config.An*config.noiseCurrentBase/device.C;

zOscillator = complex(zeros(Ngamma,config.j_max,nBits));
zReference = complex(zeros(1,nBits));
windowWeightSum = 0;
zLockOscillator = complex(zeros(stateSize));
zLockReference = complex(0);

trace = empty_trace();
[initialLockEnvelope,initialSynapticEnvelope,initialShiEnvelope] = ...
    sequential_injection_envelopes(0,fhiStartTime,lockOnlyEndTime, ...
    lockRampEndTime,synapticStartTime,synapticRampEndTime, ...
    shiRampEndTime,config);
[~,zeroGammaIndex] = min(abs(gammaValues));
[~,representativeGammaIndex] = min(abs( ...
    gammaValues-config.representativeGamma));
if storeTrace
    nSavedMaximum = ceil(nSteps/config.saveEverySteps)+2;
    trace.time = zeros(nSavedMaximum,1);
    trace.normalized_time = zeros(nSavedMaximum,1);
    trace.voltage = zeros(nSavedMaximum,1);
    trace.current = zeros(nSavedMaximum,1);
    trace.fhi_current = zeros(nSavedMaximum,1);
    trace.fhi_lock_current = zeros(nSavedMaximum,1);
    trace.fhi_synaptic_current = zeros(nSavedMaximum,1);
    trace.shi_current = zeros(nSavedMaximum,1);
    trace.shi_conductance = zeros(nSavedMaximum,1);
    trace.fhi_envelope = zeros(nSavedMaximum,1);
    trace.synaptic_fhi_envelope = zeros(nSavedMaximum,1);
    trace.shi_envelope = zeros(nSavedMaximum,1);
    trace.save_count = 1;
    trace.time(1) = 0;
    trace.normalized_time(1) = 0;
    trace.voltage(1) = initialVoltage;
    trace.current(1) = initialCurrent;
    trace.fhi_envelope(1) = initialLockEnvelope;
    trace.synaptic_fhi_envelope(1) = initialSynapticEnvelope;
    trace.shi_envelope(1) = initialShiEnvelope;
    trace.gamma = gammaValues(representativeGammaIndex);
    trace.initial_phase_deg = initialPhaseDeg;
    trace.phiFHI_deg = config.phiFHI_deg;
    trace.preparation_offset_deg = preparationOffsetDeg;
    trace.lock_target_phase_deg = lockTargetPhaseDeg;
    trace.fhi_drive_phase_deg = fhiDrivePhaseDeg;
    trace.synaptic_drive_phase_deg = synapticDrivePhaseDeg;
    trace.lock_fhi_off_time_s = lockRampEndTime;
    trace.sampling_start_time_s = synapticStartTime;
    trace.shi_full_time_s = shiRampEndTime;
    trace.lock_calibration = lockCalibration;
    trace.neutral_calibration = neutralCalibration;
    trace.reference_frequency_Hz = referenceOmega/(2*pi);
    trace.reference_omega_rad_s = referenceOmega;
    trace.reference_start_time_s = fhiStartTime;
    trace.limit_cycle_period_cv = limitCycleInfo.period_cv;
    trace.limit_cycle_amplitude_cv = limitCycleInfo.amplitude_cv;
    trace.limit_cycle_settling_cycles = limitCycleInfo.settling_cycles;
end

rng(curveSeed,'twister');
simulationTimer = tic;
progressInterval = max(1,floor(nSteps/10));

for step = 1:nSteps
    tNow = (step-1)*dt;
    tNext = min(step*dt,physicalStopTime);
    dtStep = tNext-tNow;
    referenceCyclesNext = tNext/referencePeriod;

    % FHI0 prepares the calibrated SHI-basin midpoint. It is fully off
    % before the quadrature synaptic FHI and parametric SHI begin sampling.
    referenceTimeNow = tNow-fhiStartTime;
    referenceTimeNext = tNext-fhiStartTime;
    fhiReferenceNext = sin(referenceOmega*referenceTimeNext);
    lockInjectionNow = sin(referenceOmega*referenceTimeNow+ ...
        deg2rad(fhiDrivePhaseDeg));
    lockInjectionNext = sin(referenceOmega*referenceTimeNext+ ...
        deg2rad(fhiDrivePhaseDeg));
    synapticInjectionNow = sin(referenceOmega*referenceTimeNow+ ...
        deg2rad(synapticDrivePhaseDeg));
    synapticInjectionNext = sin(referenceOmega*referenceTimeNext+ ...
        deg2rad(synapticDrivePhaseDeg));

    [lockEnvelopeNow,synapticEnvelopeNow,shiEnvelopeNow] = ...
        sequential_injection_envelopes(tNow,fhiStartTime, ...
        lockOnlyEndTime,lockRampEndTime,synapticStartTime, ...
        synapticRampEndTime,shiRampEndTime,config);
    [lockEnvelopeNext,synapticEnvelopeNext,shiEnvelopeNext] = ...
        sequential_injection_envelopes(tNext,fhiStartTime, ...
        lockOnlyEndTime,lockRampEndTime,synapticStartTime, ...
        synapticRampEndTime,shiRampEndTime,config);

    lockFhiCurrentNow = lockEnvelopeNow*config.IlockFHI_A*lockInjectionNow;
    lockFhiCurrentNext = lockEnvelopeNext*config.IlockFHI_A*lockInjectionNext;
    synapticFhiCurrentNow = synapticEnvelopeNow*config.IfhiUnit ...
        *gammaValues(:)*synapticInjectionNow;
    synapticFhiCurrentNext = synapticEnvelopeNext*config.IfhiUnit ...
        *gammaValues(:)*synapticInjectionNext;
    fhiCurrentNow = lockFhiCurrentNow+synapticFhiCurrentNow;
    fhiCurrentNext = lockFhiCurrentNext+synapticFhiCurrentNext;
    shiConductanceNow = config.GshiUnit*shiEnvelopeNow* ...
        sin(2*referenceOmega*referenceTimeNow+deg2rad(config.phiSHI_deg));
    shiConductanceNext = config.GshiUnit*shiEnvelopeNext* ...
        sin(2*referenceOmega*referenceTimeNext+deg2rad(config.phiSHI_deg));

    externalCurrentNow = fhiCurrentNow;
    externalCurrentNext = fhiCurrentNext;
    dW = sqrt(dtStep)*randn(stateSize);
    [iState,vState] = stochastic_heun_lc_step(iState,vState, ...
        externalCurrentNow,externalCurrentNext,sigmaVoltage,dW,device, ...
        dtStep,shiConductanceNow,shiConductanceNext);

    if tNext>=lockWindowStartTime && tNext<=lockOnlyEndTime
        lockPosition = (tNext-lockWindowStartTime) ...
            /(lockOnlyEndTime-lockWindowStartTime);
        lockWeight = sin(pi*lockPosition)^2;
        lockKernel = exp(-1i*referenceOmega*tNext);
        zLockOscillator = zLockOscillator+lockWeight*vState*lockKernel;
        zLockReference = zLockReference ...
            +lockWeight*fhiReferenceNext*lockKernel;
    end

    % Measure after SHI has settled and, by default, FHI has been removed.
    % Store raw phase against the zero-phase reference; classify only after
    % subtracting the independently calibrated, gamma-independent axis.
    if tNext>=phaseWindowStartTime && tNext<=measurementTime
        windowPosition = (tNext-phaseWindowStartTime) ...
            /(measurementTime-phaseWindowStartTime);
        windowWeight = sin(pi*windowPosition)^2;
        phaseKernel = exp(-1i*referenceOmega*tNext);
        windowWeightSum = windowWeightSum+windowWeight;
        zOscillator(:,:,1) = zOscillator(:,:,1) ...
            +windowWeight*vState*phaseKernel;
        zReference(1) = zReference(1) ...
            +windowWeight*fhiReferenceNext*phaseKernel;
    end

    if storeTrace && (mod(step,config.saveEverySteps)==0 || step==nSteps)
        trace.save_count = trace.save_count+1;
        q = trace.save_count;
        representativeVoltage = vState(representativeGammaIndex,1);
        trace.time(q) = tNext;
        trace.normalized_time(q) = referenceCyclesNext;
        trace.voltage(q) = representativeVoltage;
        trace.current(q) = iState(representativeGammaIndex,1);
        trace.fhi_current(q) = fhiCurrentNext(representativeGammaIndex);
        trace.fhi_lock_current(q) = lockFhiCurrentNext;
        trace.fhi_synaptic_current(q) = ...
            synapticFhiCurrentNext(representativeGammaIndex);
        trace.shi_current(q) = shiConductanceNext*representativeVoltage;
        trace.shi_conductance(q) = shiConductanceNext;
        trace.fhi_envelope(q) = lockEnvelopeNext;
        trace.synaptic_fhi_envelope(q) = synapticEnvelopeNext;
        trace.shi_envelope(q) = shiEnvelopeNext;
    end

    if mod(step,5000)==0 && ...
            any(~isfinite(iState(:)) | ~isfinite(vState(:)))
        error('Non-finite LC state at reference cycle %.6g.', ...
            referenceCyclesNext);
    end
    if mod(step,progressInterval)==0 || step==nSteps
        fprintf('  state-space integration: %3.0f%% complete\n', ...
            100*step/nSteps);
    end
end

referencePhasors = reshape(zReference,1,1,nBits);
phaseSamplesDeg = wrap_to_180(rad2deg(angle( ...
    zOscillator.*conj(referencePhasors))));
if abs(zLockReference)==0 || any(abs(zLockOscillator(:))==0)
    error('Undefined FHI-only lock phasor before the synaptic input begins.');
end
lockedPhaseSamplesDeg = wrap_to_180(rad2deg(angle( ...
    zLockOscillator.*conj(zLockReference))));
lockedPhasorMean = mean(exp(1i*deg2rad(lockedPhaseSamplesDeg(:))));
lockedMeanPhaseDeg = wrap_to_180(rad2deg(angle(lockedPhasorMean)));
lockedMeanErrorDeg = wrap_to_180( ...
    lockedMeanPhaseDeg-lockTargetPhaseDeg);

% Phases are measured against the prescribed zero-phase sinusoid.
% Only the independent SHI-only calibration axis enters classification.
% No gamma-dependent fit, observed noisy axis or label reversal is used.
if any(~isfinite(phaseSamplesDeg(:))) || any(abs(zOscillator(:))==0) || ...
        any(abs(zReference(:))==0)
    error('Undefined measurement phasor: cannot assign a physical phase bit.');
end
fundamentalAmplitude = 2*abs(zOscillator)/windowWeightSum;
minimumAmplitude = config.minimumFundamentalFraction ...
    *calibration.free_running_amplitude_V;
lowAmplitude = fundamentalAmplitude<minimumAmplitude;
if any(mean(lowAmplitude,2)>config.maximumInvalidPhaseFraction)
    error(['Too many low-amplitude readouts for a physical phase bit. ' ...
        'Check noise, readout window and SHI strength; samples were NOT removed.']);
end
calibratedPhase = wrap_to_180(phaseSamplesDeg-phaseAxisDeg);
bitsGammaTrialBit = calibrated_phase_bits(calibratedPhase);
Pbit = mean(reshape(bitsGammaTrialBit,Ngamma,[]),2).';

% The observed doubled-phase axis is modulo 180 deg and diagnostic only.
% A small concentration means this axis does not describe tight clusters.
doubledMean = mean(exp(1i*2*deg2rad(phaseSamplesDeg)),2);
audit = struct();
audit.raw_relative_phase_deg = phaseSamplesDeg;
audit.calibrated_relative_phase_deg = calibratedPhase;
audit.fundamental_amplitude_V = fundamentalAmplitude;
audit.low_amplitude = lowAmplitude;
audit.minimum_amplitude_V = minimumAmplitude;
audit.two_state_calibration = calibration;
audit.lock_calibration = lockCalibration;
audit.neutral_calibration = neutralCalibration;
audit.locked_phase_samples_deg = lockedPhaseSamplesDeg;
audit.locked_mean_phase_deg = lockedMeanPhaseDeg;
audit.locked_mean_phase_error_deg = lockedMeanErrorDeg;
audit.locked_phase_concentration = abs(lockedPhasorMean);
audit.readout_fhi_removed = config.removeFhiBeforeReadout;
audit.binary_assignment = bitsGammaTrialBit;
audit.fixed_axis_deg = phaseAxisDeg;
audit.observed_axis_deg = 0.5*rad2deg(angle(doubledMean));
audit.observed_axis_deg(abs(doubledMean)<1e-12) = NaN;
audit.doubled_phase_concentration = abs(doubledMean);
audit.near_boundary = abs(abs(calibratedPhase)-90) ...
    <=config.boundaryWarningWidth_deg;
audit.gamma_values = gammaValues;
audit.phiFHI_deg = config.phiFHI_deg; % offset from calibrated bit-1 axis
audit.preparation_offset_deg = preparationOffsetDeg;
audit.lock_target_phase_deg = lockTargetPhaseDeg;
audit.fhi_drive_phase_deg = fhiDrivePhaseDeg;
audit.synaptic_drive_phase_deg = synapticDrivePhaseDeg;
audit.neutral_offset_error_deg = wrap_to_180( ...
    lockCalibration.measured_phase_deg-phaseAxisDeg-preparationOffsetDeg);
audit.lock_current_A = config.IlockFHI_A;
audit.synaptic_current_per_gamma_A = config.IfhiUnit;
audit.lock_fhi_off_before_sampling = ...
    lockRampEndTime<=synapticStartTime;
audit.lock_fhi_off_time_s = lockRampEndTime;
audit.sampling_start_time_s = synapticStartTime;
audit.shi_full_time_s = shiRampEndTime;
audit.nearest_zero_gamma = gammaValues(zeroGammaIndex);
if abs(audit.nearest_zero_gamma)<1e-12
    audit.zero_input_bias = Pbit(zeroGammaIndex)-0.5;
else
    audit.zero_input_bias = NaN;
end
mirrorError = nan(size(gammaValues));
for gammaIndex = 1:Ngamma
    [separation,mirrorIndex] = min(abs(gammaValues+gammaValues(gammaIndex)));
    if separation<1e-10
        mirrorError(gammaIndex) = abs( ...
            Pbit(gammaIndex)+Pbit(mirrorIndex)-1);
    end
end
audit.mirror_symmetry_mae = mean(mirrorError,'omitnan');
audit.reference_frequency_Hz = referenceFrequencyHz;
audit.assignment_convention = config.bitAssignmentConvention;
audit.slope_convention = config.slopeConvention;
audit.endpoint_probability_difference = Pbit(end)-Pbit(1);
audit.boundary_warning_width_deg = config.boundaryWarningWidth_deg;

if storeTrace
    fields_to_trim = {'time','normalized_time','voltage','current', ...
        'fhi_current','fhi_lock_current','fhi_synaptic_current', ...
        'shi_current','shi_conductance', ...
        'fhi_envelope','synaptic_fhi_envelope','shi_envelope'};
    for field_index = 1:numel(fields_to_trim)
        field_name = fields_to_trim{field_index};
        trace.(field_name) = trace.(field_name)(1:trace.save_count);
    end
    trace.stored = true;
    trace.runtime_seconds = toc(simulationTimer);
    trace.phase_axis_deg = phaseAxisDeg;
end

fprintf('  P(1) at gamma min/zero/max = %.3f / %.3f / %.3f\n', ...
    Pbit(1),Pbit(zeroGammaIndex),Pbit(end));
fprintf('  zero-input bias / mirror-symmetry MAE = %+.3f / %.3f\n', ...
    audit.zero_input_bias,audit.mirror_symmetry_mae);
fprintf('  stochastic lock mean phase error = %+.2f deg\n', ...
    lockedMeanErrorDeg);
if isfinite(audit.zero_input_bias) && ...
        abs(audit.zero_input_bias)>3*sqrt(0.25/config.j_max)
    warning(['The zero-input probability differs from 0.5 by more ' ...
        'than three nominal binomial standard errors. Inspect the ' ...
        'neutrality and phase diagnostics before claiming unbiased sampling.']);
end
end

function lock = calibrate_fhi_lock( ...
        device,omega,period,orbitInfo,targetPhaseDeg,config)
% Tune the common applied FHI phase until the measured voltage fundamental
% reaches the requested phase. Check eight autonomous starting phases so
% the result is a lock, rather than a lucky single initial condition.
drivePhaseDeg = wrap_to_180(targetPhaseDeg);
lock = struct('passed',false,'target_phase_deg',targetPhaseDeg, ...
    'drive_phase_deg',nan,'measured_phase_deg',nan, ...
    'phase_error_deg',nan,'maximum_phase_spread_deg',nan, ...
    'maximum_window_drift_deg',nan,'fundamental_amplitude_V',[], ...
    'initial_phases_deg',config.lockInitialPhases_deg, ...
    'iterations',0,'current_A',config.IlockFHI_A);
for iteration = 1:config.lockPhaseMaxIterations
    [finalPhase,previousPhase,amplitude] = simulate_fhi_only_lock( ...
        device,omega,period,orbitInfo,drivePhaseDeg,config);
    meanPhase = wrap_to_180(rad2deg(angle(mean( ...
        exp(1i*deg2rad(finalPhase))))));
    phaseError = wrap_to_180(meanPhase-targetPhaseDeg);
    spread = max(abs(wrap_to_180(finalPhase-meanPhase)));
    drift = max(abs(wrap_to_180(finalPhase-previousPhase)));
    amplitudeOk = all(amplitude>=config.lockMinimumAmplitudeFraction ...
        *orbitInfo.free_running_amplitude_V);
    lock.drive_phase_deg = drivePhaseDeg;
    lock.measured_phase_deg = meanPhase;
    lock.phase_error_deg = phaseError;
    lock.maximum_phase_spread_deg = spread;
    lock.maximum_window_drift_deg = drift;
    lock.fundamental_amplitude_V = amplitude;
    lock.iterations = iteration;
    lock.passed = amplitudeOk && ...
        abs(phaseError)<=config.lockPhaseTolerance_deg && ...
        spread<=config.lockPhaseTolerance_deg && ...
        drift<=config.lockPhaseDriftTolerance_deg;
    if lock.passed
        return
    end
    drivePhaseDeg = wrap_to_180(drivePhaseDeg-phaseError);
end
error(['The %.3f-mA first FHI did not lock the voltage to %.1f deg ' ...
    '(measured %.2f deg, drive %.2f deg, spread %.2f deg, ' ...
    'window drift %.2f deg). Adjust IlockFHI_A or fhiLockCycles.'], ...
    1e3*config.IlockFHI_A,targetPhaseDeg, ...
    lock.measured_phase_deg,lock.drive_phase_deg, ...
    lock.maximum_phase_spread_deg,lock.maximum_window_drift_deg);
end

function neutral = find_neutral_separatrix( ...
        device,omega,period,iInitial,vInitial,phaseAxisDeg, ...
        nominalDrivePhaseDeg,config)
% Locate the unstable zero-synaptic-input basin boundary under the actual
% FHI0 release and SHI-ramp protocol. This is deterministic: it does not
% fit or force the stochastic P(1|gamma=0) value.
lowerOffset = 0;
upperOffset = 180;
lowerBit = deterministic_zero_input_bit(lowerOffset);
upperBit = deterministic_zero_input_bit(upperOffset);
if lowerBit==upperBit
    error(['The two endpoints of the preparation-phase sweep reach the ' ...
        'same SHI basin. No neutral separatrix can be certified.']);
end
for iteration = 1:config.neutralBisectionIterations
    midpoint = 0.5*(lowerOffset+upperOffset);
    midpointBit = deterministic_zero_input_bit(midpoint);
    if midpointBit==lowerBit
        lowerOffset = midpoint;
    else
        upperOffset = midpoint;
    end
end
neutral = struct('performed',true, ...
    'offset_deg',0.5*(lowerOffset+upperOffset), ...
    'bracket_deg',[lowerOffset upperOffset], ...
    'endpoint_bits',[lowerBit upperBit], ...
    'iterations',config.neutralBisectionIterations, ...
    'method','deterministic_zero_input_SHI_ramp_basin_boundary_v1');
fprintf('  deterministic neutral offset = %.3f deg (bracket %.3f to %.3f)\n', ...
    neutral.offset_deg,lowerOffset,upperOffset);

    function bit = deterministic_zero_input_bit(offsetDeg)
        drivePhaseDeg = wrap_to_180( ...
            nominalDrivePhaseDeg+offsetDeg-90);
        dt = period/config.stepsPerPeriod;
        lockOnlyEnd = config.fhiLockCycles*period;
        lockRampEnd = lockOnlyEnd+config.lockFhiRampOffTime_s;
        synapticStart = lockRampEnd;
        synapticRampEnd = synapticStart+config.synapticFhiRampTime_s;
        shiRampEnd = synapticStart+config.shiRampTime_s;
        stopTime = shiRampEnd ...
            +double(config.removeFhiBeforeReadout)*config.fhiRampOffTime_s ...
            +config.combinedSettleCycles*period;
        windowStart = stopTime-config.phaseWindowCycles*period;
        i = iInitial;
        v = vInitial;
        zOscillator = complex(0);
        zReference = complex(0);
        for step = 1:ceil(stopTime/dt)
            t = (step-1)*dt;
            tn = min(step*dt,stopTime);
            h = tn-t;
            [lockNow,~,shiNow] = sequential_injection_envelopes( ...
                t,0,lockOnlyEnd,lockRampEnd,synapticStart, ...
                synapticRampEnd,shiRampEnd,config);
            [lockNext,~,shiNext] = sequential_injection_envelopes( ...
                tn,0,lockOnlyEnd,lockRampEnd,synapticStart, ...
                synapticRampEnd,shiRampEnd,config);
            currentNow = config.IlockFHI_A*lockNow ...
                *sin(omega*t+deg2rad(drivePhaseDeg));
            currentNext = config.IlockFHI_A*lockNext ...
                *sin(omega*tn+deg2rad(drivePhaseDeg));
            conductanceNow = config.GshiUnit*shiNow ...
                *sin(2*omega*t+deg2rad(config.phiSHI_deg));
            conductanceNext = config.GshiUnit*shiNext ...
                *sin(2*omega*tn+deg2rad(config.phiSHI_deg));
            [i,v] = stochastic_heun_lc_step(i,v,currentNow, ...
                currentNext,0,0,device,h,conductanceNow,conductanceNext);
            if ~isfinite(i) || ~isfinite(v)
                error('Non-finite state in neutral-separatrix calibration.');
            end
            if tn>=windowStart
                u = (tn-windowStart)/(stopTime-windowStart);
                weight = sin(pi*u)^2;
                kernel = exp(-1i*omega*tn);
                zOscillator = zOscillator+weight*v*kernel;
                zReference = zReference+weight*sin(omega*tn)*kernel;
            end
        end
        if abs(zOscillator)==0 || abs(zReference)==0
            error('Undefined readout phase in neutral calibration.');
        end
        relative = wrap_to_180(rad2deg(angle( ...
            zOscillator*conj(zReference)))-phaseAxisDeg);
        bit = calibrated_phase_bits(relative);
    end
end

function [finalPhase,previousPhase,amplitude] = simulate_fhi_only_lock( ...
        device,omega,period,orbitInfo,drivePhaseDeg,config)
orbitPhase = rad2deg(atan2(omega*device.L*orbitInfo.orbit_current_A, ...
    orbitInfo.orbit_voltage_V));
initialPhases = config.lockInitialPhases_deg(:).';
iState = zeros(size(initialPhases));
vState = zeros(size(initialPhases));
for k = 1:numel(initialPhases)
    [~,index] = min(abs(wrap_to_180(orbitPhase-initialPhases(k))));
    iState(k) = orbitInfo.orbit_current_A(index);
    vState(k) = orbitInfo.orbit_voltage_V(index);
end
dt = period/config.stepsPerPeriod;
stopTime = config.fhiLockCycles*period;
windowDuration = config.lockPhaseWindowCycles*period;
windowStart = [stopTime-2*windowDuration stopTime-windowDuration];
zOscillator = complex(zeros(2,numel(initialPhases)));
zReference = complex(zeros(2,1));
weightSum = zeros(2,1);
for step = 1:ceil(stopTime/dt)
    t = (step-1)*dt;
    tNext = min(step*dt,stopTime);
    h = tNext-t;
    currentNow = config.IlockFHI_A*sin( ...
        omega*t+deg2rad(drivePhaseDeg));
    currentNext = config.IlockFHI_A*sin( ...
        omega*tNext+deg2rad(drivePhaseDeg));
    [iState,vState] = stochastic_heun_lc_step(iState,vState, ...
        currentNow,currentNext,0,0,device,h,0,0);
    if any(~isfinite(iState) | ~isfinite(vState))
        error('Non-finite state during the first-FHI lock calibration.');
    end
    for window = 1:2
        u = (tNext-windowStart(window))/windowDuration;
        if u>=0 && u<=1
            weight = sin(pi*u)^2;
            kernel = exp(-1i*omega*tNext);
            zOscillator(window,:) = zOscillator(window,:) ...
                +weight*vState*kernel;
            zReference(window) = zReference(window) ...
                +weight*sin(omega*tNext)*kernel;
            weightSum(window) = weightSum(window)+weight;
        end
    end
end
if any(abs(zReference)==0) || any(abs(zOscillator(:))==0)
    error('Undefined phasor during the first-FHI lock calibration.');
end
relativePhase = wrap_to_180(rad2deg(angle( ...
    zOscillator.*conj(zReference))));
previousPhase = relativePhase(1,:);
finalPhase = relativePhase(2,:);
amplitude = 2*abs(zOscillator(2,:))/weightSum(2);
end

function [measuredOmega,measuredPeriod,iInitial,vInitial,info] = ...
    characterize_limit_cycle(device,targetPhaseDeg,config)
% Settle the autonomous nonlinear oscillator, verify period and amplitude
% convergence, and measure its free-running frequency from positive-going
% zero crossings of v(t). No analytical LC-frequency assumption is used in
% the returned reference frequency.

provisionalPeriod = device.T0;
dtCalibration = provisionalPeriod/config.stepsPerPeriod;
minimumSettleTime = ...
    config.frequencyMinSettleCycles*provisionalPeriod;
maximumSteps = ceil(config.frequencyMaxSettleCycles ...
    *provisionalPeriod/dtCalibration);
convergencePeriods = config.frequencyConvergencePeriods;
tolerance = config.frequencyConvergenceTolerance;

maximumRecordedCycles = 2*config.frequencyMaxSettleCycles+10;
crossingTimes = nan(maximumRecordedCycles,1);
cycleAmplitudes = nan(maximumRecordedCycles,1);
crossingCount = 0;
amplitudeCount = 0;

iState = 0;
vState = config.initialVoltage;
previousVoltage = vState;
previousTime = 0;
cycleMaximum = vState;
cycleMinimum = vState;
converged = false;
periodCv = inf;
amplitudeCv = inf;

for step = 1:maximumSteps
    currentTime = step*dtCalibration;
    [iNext,vNext] = rk4_lc_step( ...
        iState,vState,0,device,dtCalibration);
    if ~isfinite(iNext) || ~isfinite(vNext)
        error('Non-finite state during autonomous frequency calibration.');
    end

    cycleMaximum = max(cycleMaximum,vNext);
    cycleMinimum = min(cycleMinimum,vNext);

    if previousVoltage<0 && vNext>=0
        crossingFraction = -previousVoltage/(vNext-previousVoltage);
        crossingTime = previousTime+crossingFraction*dtCalibration;
        crossingCount = crossingCount+1;
        crossingTimes(crossingCount) = crossingTime;

        if crossingCount>1
            amplitudeCount = amplitudeCount+1;
            cycleAmplitudes(amplitudeCount) = ...
                0.5*(cycleMaximum-cycleMinimum);
        end
        cycleMaximum = vNext;
        cycleMinimum = vNext;

        if currentTime>=minimumSettleTime && ...
                crossingCount>=convergencePeriods+1 && ...
                amplitudeCount>=convergencePeriods
            recentPeriods = diff(crossingTimes( ...
                crossingCount-convergencePeriods:crossingCount));
            recentAmplitudes = cycleAmplitudes( ...
                amplitudeCount-convergencePeriods+1:amplitudeCount);
            meanPeriod = mean(recentPeriods);
            meanAmplitude = mean(abs(recentAmplitudes));
            periodCv = std(recentPeriods,0)/meanPeriod;
            amplitudeCv = std(recentAmplitudes,0)/meanAmplitude;
            if meanPeriod>0 && meanAmplitude>0 && ...
                    periodCv<=tolerance && amplitudeCv<=tolerance
                iState = iNext;
                vState = vNext;
                converged = true;
                break
            end
        end
    end

    iState = iNext;
    vState = vNext;
    previousVoltage = vState;
    previousTime = currentTime;
end

if ~converged
    error(['The autonomous oscillator did not meet the limit-cycle ' ...
        'frequency/amplitude convergence tolerance within %d nominal ' ...
        'cycles. Last period CV = %.4g; amplitude CV = %.4g.'], ...
        config.frequencyMaxSettleCycles,periodCv,amplitudeCv);
end

measuredPeriod = mean(recentPeriods);
measuredOmega = 2*pi/measuredPeriod;

% Record one measured period and select the state closest to the requested
% oscillator phase using the measured angular frequency.
dtOrbit = measuredPeriod/config.stepsPerPeriod;
iOrbit = zeros(config.stepsPerPeriod,1);
vOrbit = zeros(config.stepsPerPeriod,1);
for orbitIndex = 1:config.stepsPerPeriod
    [iState,vState] = rk4_lc_step(iState,vState,0,device,dtOrbit);
    iOrbit(orbitIndex) = iState;
    vOrbit(orbitIndex) = vState;
end

statePhase = atan2(measuredOmega*device.L*iOrbit,vOrbit);
targetPhase = deg2rad(targetPhaseDeg);
wrappedDistance = abs(atan2(sin(statePhase-targetPhase), ...
    cos(statePhase-targetPhase)));
[~,closestIndex] = min(wrappedDistance);
iInitial = iOrbit(closestIndex);
vInitial = vOrbit(closestIndex);

info = struct();
info.orbit_current_A = iOrbit;
info.orbit_voltage_V = vOrbit;
info.free_running_amplitude_V = 2*abs(sum(vOrbit.* ...
    exp(-1i*2*pi*(1:numel(vOrbit)).'/numel(vOrbit))))/numel(vOrbit);
info.period_cv = periodCv;
info.amplitude_cv = amplitudeCv;
info.settling_cycles = currentTime/provisionalPeriod;
info.crossings_used = convergencePeriods+1;
end

function [iInitial,vInitial] = initial_state_on_limit_cycle( ...
    device,targetPhaseDeg,stepsPerPeriod,startupCycles,initialVoltage)
% Settle the autonomous device, store its final period, and select the state
% closest to the requested physical phase. For an ideal LC orbit,
% theta=atan2(omega0*L*i,v), because v=A*cos(theta) and
% i=A/(omega0*L)*sin(theta).

dtStartup = device.T0/stepsPerPeriod;
nStartupSteps = startupCycles*stepsPerPeriod;
iState = 0;
vState = initialVoltage;
iOrbit = zeros(stepsPerPeriod,1);
vOrbit = zeros(stepsPerPeriod,1);

for step = 1:nStartupSteps
    [iState,vState] = rk4_lc_step(iState,vState,0,device,dtStartup);
    if step>nStartupSteps-stepsPerPeriod
        orbitIndex = step-(nStartupSteps-stepsPerPeriod);
        iOrbit(orbitIndex) = iState;
        vOrbit(orbitIndex) = vState;
    end
end

if any(~isfinite(iOrbit)) || any(~isfinite(vOrbit)) || ...
        max(abs(vOrbit))<1e-6
    error('Device did not reach a finite autonomous limit cycle.');
end

statePhase = atan2(device.omega0*device.L*iOrbit,vOrbit);
targetPhase = deg2rad(targetPhaseDeg);
wrappedDistance = abs(atan2(sin(statePhase-targetPhase), ...
    cos(statePhase-targetPhase)));
[~,closestIndex] = min(wrappedDistance);
iInitial = iOrbit(closestIndex);
vInitial = vOrbit(closestIndex);
end

function [iNext,vNext] = stochastic_heun_lc_step( ...
    i,v,IextNow,IextNext,sigmaVoltage,dW,device,dt,GNow,GNext)
[diNow,dvNow] = lc_state_derivative(i,v,IextNow,device,GNow);
iPredictor = i+diNow*dt;
vPredictor = v+dvNow*dt+sigmaVoltage*dW;
[diPredictor,dvPredictor] = lc_state_derivative( ...
    iPredictor,vPredictor,IextNext,device,GNext);
iNext = i+0.5*(diNow+diPredictor)*dt;
vNext = v+0.5*(dvNow+dvPredictor)*dt+sigmaVoltage*dW;
end

function [iNext,vNext] = rk4_lc_step(i,v,Iext,device,dt)
[k1i,k1v] = lc_state_derivative(i,v,Iext,device);
[k2i,k2v] = lc_state_derivative( ...
    i+0.5*dt*k1i,v+0.5*dt*k1v,Iext,device);
[k3i,k3v] = lc_state_derivative( ...
    i+0.5*dt*k2i,v+0.5*dt*k2v,Iext,device);
[k4i,k4v] = lc_state_derivative(i+dt*k3i,v+dt*k3v,Iext,device);
iNext = i+(dt/6)*(k1i+2*k2i+2*k3i+k4i);
vNext = v+(dt/6)*(k1v+2*k2v+2*k3v+k4v);
end

function [di_dt,dv_dt] = lc_state_derivative(i,v,Iext,device,Gshi)
if nargin<5, Gshi = 0; end
di_dt = v/device.L;
nonlinearCurrent = device.S*tanh((device.Gn/device.S).*v);
dv_dt = (-i-v/device.R-nonlinearCurrent+Iext+Gshi.*v)/device.C;
end

function [lockEnvelope,synapticEnvelope,shiEnvelope] = ...
    sequential_injection_envelopes(time,fhiStartTime,lockOnlyEndTime, ...
    lockRampEndTime,synapticStartTime,synapticRampEndTime, ...
    shiRampEndTime,config)
% Paper sequence: FHI0 prepares the neutral phase with SHI=0; FHI0 is
% completely removed before synaptic FHI and the low-to-high SHI ramp.
if time<fhiStartTime
    lockEnvelope = 0;
    synapticEnvelope = 0;
    shiEnvelope = 0;
    return
end

if time<=lockOnlyEndTime
    lockEnvelope = 1;
elseif time<lockRampEndTime
    fraction = (time-lockOnlyEndTime)/config.lockFhiRampOffTime_s;
    lockEnvelope = 0.5*(1+cos(pi*fraction));
else
    lockEnvelope = 0;
end
if time<=synapticStartTime
    synapticEnvelope = 0;
elseif time<synapticRampEndTime
    rampFraction = (time-synapticStartTime) ...
        /(synapticRampEndTime-synapticStartTime);
    synapticEnvelope = 0.5*(1-cos(pi*rampFraction));
else
    synapticEnvelope = 1;
    if config.removeFhiBeforeReadout && time>shiRampEndTime
        fraction = min(1,(time-shiRampEndTime)/config.fhiRampOffTime_s);
        synapticEnvelope = 0.5*(1+cos(pi*fraction));
    end
end
if time<=synapticStartTime
    shiEnvelope = 0;
elseif time<shiRampEndTime
    rampFraction = (time-synapticStartTime) ...
        /(shiRampEndTime-synapticStartTime);
    smoothFraction = 0.5*(1-cos(pi*rampFraction));
    shiEnvelope = config.As*smoothFraction;
else
    shiEnvelope = config.As;
end
end

function bits = calibrated_phase_bits(relativePhaseDeg)
% The independently calibrated 0-deg basin is bit 1; the antipodal
% 180-deg basin is bit 0. Exact +/-90-deg ties are assigned to bit 1.
if any(~isfinite(relativePhaseDeg(:)))
    error('Relative phase must be finite.');
end
wrappedPhaseDeg = wrap_to_180(relativePhaseDeg);
bits = abs(wrappedPhaseDeg)<=90;
end

function phaseWrappedDeg = wrap_to_180(phaseDeg)
phaseWrappedDeg = mod(phaseDeg+180,360)-180;
end

function trace = empty_trace()
trace = struct('stored',false,'time',[],'normalized_time',[], ...
    'voltage',[],'current',[],'fhi_current',[], ...
    'fhi_lock_current',[],'fhi_synaptic_current',[],'shi_current',[], ...
    'shi_conductance',[], ...
    'fhi_envelope',[],'synaptic_fhi_envelope',[], ...
    'shi_envelope',[],'save_count',0,'gamma',[], ...
    'initial_phase_deg',[],'phiFHI_deg',[], ...
    'preparation_offset_deg',[],'fhi_drive_phase_deg',[], ...
    'lock_target_phase_deg',[],'synaptic_drive_phase_deg',[], ...
    'lock_fhi_off_time_s',[],'sampling_start_time_s',[], ...
    'shi_full_time_s',[],'lock_calibration',[], ...
    'neutral_calibration',[],'runtime_seconds',[], ...
    'phase_axis_deg',[],'parameter_set_index',[], ...
    'reference_frequency_Hz',[],'reference_omega_rad_s',[], ...
    'reference_start_time_s',[],'limit_cycle_period_cv',[], ...
    'limit_cycle_amplitude_cv',[],'limit_cycle_settling_cycles',[]);
end

function validate_config(config)
validate_parametric_config(config);
if ~ischar(config.referenceFrequencyMode) || ...
        ~ismember(config.referenceFrequencyMode, ...
        {'fixed_nominal','measured_limit_cycle'})
    error(['config.referenceFrequencyMode must be ''fixed_nominal'' or ' ...
        '''measured_limit_cycle''.']);
end
if ~isscalar(config.fhiStartTime_s) || ...
        ~isfinite(config.fhiStartTime_s) || config.fhiStartTime_s<0
    error('config.fhiStartTime_s must be a finite nonnegative time.');
end
if config.j_max<1 || config.j_max~=round(config.j_max)
    error('config.j_max must be a positive integer.');
end
if config.stepsPerPeriod<40 || ...
        config.stepsPerPeriod~=round(config.stepsPerPeriod)
    error('config.stepsPerPeriod must be an integer of at least 40.');
end
if config.frequencyMinSettleCycles<2 || ...
        config.frequencyMinSettleCycles~=round( ...
        config.frequencyMinSettleCycles)
    error('config.frequencyMinSettleCycles must be an integer >= 2.');
end
if config.frequencyMaxSettleCycles<=config.frequencyMinSettleCycles || ...
        config.frequencyMaxSettleCycles~=round( ...
        config.frequencyMaxSettleCycles)
    error(['config.frequencyMaxSettleCycles must be an integer greater ' ...
        'than config.frequencyMinSettleCycles.']);
end
if config.frequencyConvergencePeriods<3 || ...
        config.frequencyConvergencePeriods~=round( ...
        config.frequencyConvergencePeriods)
    error('config.frequencyConvergencePeriods must be an integer >= 3.');
end
if ~isscalar(config.frequencyConvergenceTolerance) || ...
        ~isfinite(config.frequencyConvergenceTolerance) || ...
        config.frequencyConvergenceTolerance<=0
    error('config.frequencyConvergenceTolerance must be positive.');
end
if ~isscalar(config.initialVoltage) || ...
        ~isfinite(config.initialVoltage)
    error('config.initialVoltage must be a finite scalar voltage.');
end
if ~isscalar(config.fhiLockCycles) || ...
        ~isfinite(config.fhiLockCycles) || config.fhiLockCycles<1 || ...
        config.fhiLockCycles~=round(config.fhiLockCycles)
    error('config.fhiLockCycles must be a positive integer.');
end
if ~isscalar(config.synapticFhiRampTime_s) || ...
        ~isfinite(config.synapticFhiRampTime_s) || ...
        config.synapticFhiRampTime_s<=0
    error('config.synapticFhiRampTime_s must be positive and finite.');
end
if ~isscalar(config.shiRampTime_s) || ...
        ~isfinite(config.shiRampTime_s) || config.shiRampTime_s<=0
    error('config.shiRampTime_s must be a positive finite time.');
end
if config.synapticFhiRampTime_s>=config.shiRampTime_s
    error('The synaptic FHI must rise faster than the SHI sampling ramp.');
end
if ~isscalar(config.combinedSettleCycles) || ...
        ~isfinite(config.combinedSettleCycles) || ...
        config.combinedSettleCycles<1 || ...
        config.combinedSettleCycles~=round(config.combinedSettleCycles)
    error('config.combinedSettleCycles must be a positive integer.');
end
if ~isscalar(config.phaseWindowCycles) || ...
        ~isfinite(config.phaseWindowCycles) || ...
        config.phaseWindowCycles<=1 || ...
        config.phaseWindowCycles>config.combinedSettleCycles
    error(['config.phaseWindowCycles must exceed one and cannot exceed ' ...
        'config.combinedSettleCycles.']);
end
if any(~isfinite(config.gamma_values)) || isempty(config.gamma_values)
    error('config.gamma_values must contain finite values.');
end
if numel(config.gamma_values)<5 || any(diff(config.gamma_values)<=0)
    error('config.gamma_values must contain at least five increasing points.');
end
if ~isscalar(config.representativeGamma) || ...
        ~isfinite(config.representativeGamma) || ...
        config.representativeGamma<min(config.gamma_values) || ...
        config.representativeGamma>max(config.gamma_values)
    error(['config.representativeGamma must be a finite scalar within ' ...
        'the simulated config.gamma_values range.']);
end
if config.An<0 || config.noiseCurrentBase<0
    error('Noise amplitudes cannot be negative.');
end
if ~isscalar(config.phiFHI_deg) || ~isfinite(config.phiFHI_deg)
    error('config.phiFHI_deg must be a finite scalar angle in degrees.');
end
if isempty(config.phiFHI_sensitivity_deg) || ...
        any(~isfinite(config.phiFHI_sensitivity_deg))
    error('config.phiFHI_sensitivity_deg must contain finite angles.');
end
if any(config.phiFHI_sensitivity_deg<0 | ...
        config.phiFHI_sensitivity_deg>180)
    error('Sensitivity target phases must be within 0 to 180 degrees.');
end
if numel(config.representativeGammaTargets)~=5 || ...
        any(~isfinite(config.representativeGammaTargets)) || ...
        any(config.representativeGammaTargets<min(config.gamma_values)) || ...
        any(config.representativeGammaTargets>max(config.gamma_values))
    error('Set five representative gamma targets within the gamma sweep.');
end
selectedIndices = zeros(size(config.representativeGammaTargets));
for targetIndex = 1:numel(selectedIndices)
    [~,selectedIndices(targetIndex)] = min(abs(config.gamma_values ...
        -config.representativeGammaTargets(targetIndex)));
end
if numel(unique(selectedIndices))~=5
    error('Gamma sweep does not resolve five distinct representative inputs.');
end
end

function compatible = compatible_checkpoint(old,config,Nset,Ngamma,NphiFHI)
required = {'config','phase_audit','Pbit_all','phase_axis_all_deg', ...
    'reference_frequency_all_Hz','completed_parameter_sets','Pbit_phiFHI', ...
    'phase_axis_phiFHI_deg','completed_phiFHI_conditions'};
compatible = true;
for field_index = 1:numel(required)
    compatible = compatible && isfield(old,required{field_index});
end
if ~compatible
    return
end
compatible = isstruct(old.phase_audit) && ...
    isfield(old.phase_audit,'device') && isfield(old.phase_audit,'phiFHI');
if ~compatible
    return
end
compatible = isequaln(old.config,config) && ...
    numel(old.phase_audit.device)==Nset && ...
    numel(old.phase_audit.phiFHI)==NphiFHI && ...
    isequal(size(old.Pbit_all),[Nset Ngamma]) && ...
    isequal(size(old.reference_frequency_all_Hz),[Nset 1]) && ...
    isequal(size(old.Pbit_phiFHI),[NphiFHI Ngamma]) && ...
    isequal(size(old.completed_parameter_sets),[Nset 1]) && ...
    isequal(size(old.completed_phiFHI_conditions),[NphiFHI 1]);
end



function calibration = validate_parametric_states(device,omega,period,config)
% Independent deterministic SHI-only control, never using gamma or noisy bits.
[freeOmega,~,~,~,orbit] = characterize_limit_cycle(device,0,config);
phases = config.calibrationInitialPhases_deg(:).';
orbitPhase = rad2deg(atan2(freeOmega*device.L*orbit.orbit_current_A, ...
    orbit.orbit_voltage_V));
i = zeros(size(phases)); v = i;
for q = 1:numel(phases)
    [~,idx] = min(abs(wrap_to_180(orbitPhase-phases(q))));
    i(q) = orbit.orbit_current_A(idx);
    v(q) = orbit.orbit_voltage_V(idx);
end
dt = period/config.stepsPerPeriod;
rampEnd = config.shiRampTime_s;
stopTime = rampEnd+config.calibrationSettleCycles*period;
windowDuration = config.phaseWindowCycles*period;
windowStart = [stopTime-2*windowDuration stopTime-windowDuration];
z = complex(zeros(2,numel(phases))); weightSum = zeros(2,1);
zReference = complex(zeros(2,1));
nSteps = ceil(stopTime/dt);
for step = 1:nSteps
    t = (step-1)*dt; tn = min(step*dt,stopTime); h = tn-t;
    a0 = 0.5*(1-cos(pi*min(1,t/rampEnd)));
    a1 = 0.5*(1-cos(pi*min(1,tn/rampEnd)));
    g0 = config.GshiUnit*config.As*a0*sin(2*omega*t+deg2rad(config.phiSHI_deg));
    g1 = config.GshiUnit*config.As*a1*sin(2*omega*tn+deg2rad(config.phiSHI_deg));
    [i,v] = stochastic_heun_lc_step(i,v,0,0,0,0,device,h,g0,g1);
    if any(~isfinite(i) | ~isfinite(v))
        error('Non-finite state in independent parametric-SHI calibration.');
    end
    for window = 1:2
        u = (tn-windowStart(window))/windowDuration;
        if u>=0 && u<=1
            weight = sin(pi*u)^2;
            kernel = exp(-1i*omega*tn);
            z(window,:) = z(window,:)+weight*v*kernel;
            zReference(window) = zReference(window)+weight*sin(omega*tn)*kernel;
            weightSum(window) = weightSum(window)+weight;
        end
    end
end
amplitude = 2*abs(z)./weightSum;
phase = wrap_to_180(rad2deg(angle(z.*conj(zReference))));
axis = 0.5*rad2deg(angle(mean(exp(2i*deg2rad(phase(2,:))))));
relative = wrap_to_180(phase(2,:)-axis);
assignment = abs(relative)<=90;
distance = min(abs(relative),180-abs(relative));
drift = abs(wrap_to_180(phase(2,:)-phase(1,:)));
counts = [sum(~assignment) sum(assignment)];
center = nan(1,2);
for bit = 0:1
    selected = phase(2,assignment==bit);
    if ~isempty(selected)
        center(bit+1) = rad2deg(angle(mean(exp(1i*deg2rad(selected)))));
    end
end
separation = abs(wrap_to_180(diff(center)));
failures = {};
if any(amplitude(:)<config.minimumFundamentalFraction*orbit.free_running_amplitude_V)
    failures{end+1} = 'fundamental amplitude too small';
end
if any(counts<2)
    failures{end+1} = 'both basins were not reached by at least two initial phases';
end
if max(distance)>config.calibrationClusterTolerance_deg
    failures{end+1} = 'final phases do not form two tight antipodal clusters';
end
if ~isfinite(separation) || abs(separation-180)>config.calibrationAntipodalTolerance_deg
    failures{end+1} = 'cluster centers are not antipodal';
end
if max(drift)>config.calibrationDriftTolerance_deg
    failures{end+1} = 'phase has not settled between consecutive readout windows';
end
calibration = struct();
calibration.passed = isempty(failures);
calibration.failure_reason = strjoin(failures,'; ');
calibration.initial_phases_deg = phases;
calibration.final_phases_deg = phase(2,:);
calibration.previous_window_phases_deg = phase(1,:);
calibration.fundamental_amplitude_V = amplitude(2,:);
calibration.free_running_amplitude_V = orbit.free_running_amplitude_V;
calibration.bit0_axis_deg = wrap_to_180(axis+180);
calibration.bit1_axis_deg = axis;
calibration.cluster_centers_deg = center;
calibration.cluster_counts = counts;
calibration.cluster_separation_deg = separation;
calibration.maximum_cluster_error_deg = max(distance);
calibration.maximum_window_drift_deg = max(drift);
calibration.reference_frequency_Hz = omega/(2*pi);
calibration.GshiUnit_S = config.GshiUnit;
calibration.test = 'deterministic_SHI_only_multiple_initial_phases_v1';
% This finite test supports, but cannot mathematically prove, global bistability.
end

function validate_parametric_config(config)
validateattributes(config.IlockFHI_A,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(config.IfhiUnit,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(config.GshiUnit,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.As,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.lockFhiRampOffTime_s,{'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(config.fhiRampOffTime_s,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(config.synapticPhaseOffset_deg,{'numeric'}, ...
    {'scalar','real','finite'});
validateattributes(config.lockPhaseWindowCycles,{'numeric'}, ...
    {'scalar','integer','positive'});
if config.fhiLockCycles<=2*config.lockPhaseWindowCycles
    error('fhiLockCycles must exceed two lock-phase measurement windows.');
end
validateattributes(config.lockPhaseTolerance_deg,{'numeric'}, ...
    {'scalar','real','finite','positive','<',45});
validateattributes(config.lockPhaseDriftTolerance_deg,{'numeric'}, ...
    {'scalar','real','finite','positive','<',45});
validateattributes(config.lockPhaseMaxIterations,{'numeric'}, ...
    {'scalar','integer','positive'});
if ~islogical(config.autoCenterNeutralPhase) || ...
        ~isscalar(config.autoCenterNeutralPhase)
    error('autoCenterNeutralPhase must be a logical scalar.');
end
validateattributes(config.neutralBisectionIterations,{'numeric'}, ...
    {'scalar','integer','>=',5,'<=',16});
validateattributes(config.lockMinimumAmplitudeFraction,{'numeric'}, ...
    {'scalar','real','positive','<',1});
validateattributes(config.lockInitialPhases_deg,{'numeric'}, ...
    {'vector','real','finite','nonempty'});
if numel(unique(mod(config.lockInitialPhases_deg,360)))<4
    error('Use at least four initial phases to validate FHI locking.');
end
validateattributes(config.calibrationSettleCycles,{'numeric'}, ...
    {'scalar','integer','>',2*config.phaseWindowCycles});
validateattributes(config.minimumFundamentalFraction,{'numeric'}, ...
    {'scalar','real','positive','<',1});
validateattributes(config.maximumInvalidPhaseFraction,{'numeric'}, ...
    {'scalar','real','>=',0,'<',1});
validateattributes(config.calibrationInitialPhases_deg,{'numeric'}, ...
    {'vector','real','finite','nonempty'});
if numel(unique(mod(config.calibrationInitialPhases_deg,360)))<8
    error('Use at least eight distinct initial phases for two-state validation.');
end
for name = {'calibrationClusterTolerance_deg', ...
        'calibrationAntipodalTolerance_deg','calibrationDriftTolerance_deg'}
    validateattributes(config.(name{1}),{'numeric'}, ...
        {'scalar','real','finite','positive','<',45});
end
if ~isequal(config.requireTwoStateValidation,true)
    error('Two-state validation cannot be disabled in this implementation.');
end
if ~islogical(config.removeFhiBeforeReadout) || ~isscalar(config.removeFhiBeforeReadout)
    error('removeFhiBeforeReadout must be a logical scalar.');
end
if ~strcmp(config.shiCoupling,'parametric_conductance_v1') || ...
        ~strcmp(config.lockPhaseDefinition, ...
        'measured_SHI_basin_midpoint_v1') || ...
        ~strcmp(config.injectionProtocol, ...
        'paper_FHI0_off_before_synaptic_FHI_and_SHI_v1') || ...
        ~strcmp(config.bitAssignmentConvention, ...
        'independently_calibrated_antipodal_bit1_zero_fhi90_v2') || ...
        ~strcmp(config.slopeConvention, ...
        'signed_synaptic_FHI_quadrature_v1') || ...
        abs(wrap_to_180(config.synapticPhaseOffset_deg+90))>1e-9 || ...
        config.bit0Phase_deg~=180 || config.bit1Phase_deg~=0 || ...
        ~strcmp(config.bitBoundaryRule, ...
        'abs(wrapped_phase)<=90; exact ties are bit 1')
    error('Incompatible model/classification convention override.');
end
end
