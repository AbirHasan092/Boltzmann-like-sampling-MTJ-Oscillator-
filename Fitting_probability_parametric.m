%% LC_Oscillator_Parametric_Tanh_Fit_and_Sensitivity_Plots.m
% Fit every completed LC device probability curve to
%
%   P(bit=1) = 0.5*[1+tanh(Isyn_mA/I0_mA+b)],
%
% where I0_mA is fitted directly and beta_eff_mA_inv = 1/I0_mA.
% Isyn_mA = 1000*results.config.IfhiUnit*gamma is the signed synaptic
% FHI amplitude in mA. Beta therefore has units mA^(-1), not dimensionless.
% The device-beta PDF is shown as a smooth KDE curve and colored density points.
%
% construct an ensemble-effective response, plot parametric-SHI phase,
% amplitude, waveform, and FHI-phase diagnostics,
% and estimate standardized sensitivity of the fitted response
% to L, C, R, S, and |Gn|. This script uses fminsearch and does not require
% the Optimization or Statistics Toolboxes.
%
% Run LC_Oscillator_Parameter_Set_Generator.m and then
% Ensemble_simulation_and_IC_sensitivity_parametric.m first.

clear; clc; close all;

%% Files and plot controls
script_dir = fileparts(mfilename('fullpath'));

% This fitter accepts only the FHI0-off-before-sampling parametric-SHI data.
results_file = fullfile(script_dir, ...
    'lc_oscillator_2FHI_paramSHI_parametric_two_state.mat');




%'lc_oscillator_paper_aligned_FHI_SHI_parametric_two_state.mat');
output_suffix = '_paper_aligned_FHI_SHI_parametric_I0_current_PDF';
expected_reference_frequency_mode = 'measured_limit_cycle';
fit_file = fullfile(script_dir, ...
    ['lc_oscillator_tanh_fit_results' output_suffix '.mat']);
fit_table_file = fullfile(script_dir, ...
    ['lc_oscillator_tanh_fit_table' output_suffix '.csv']);
sensitivity_table_file = fullfile(script_dir, ...
    ['lc_oscillator_standardized_sensitivity_table' ...
    output_suffix '.csv']);
combined_barplot_file = fullfile(script_dir, ...
    ['lc_oscillator_all_parameter_deviations' output_suffix '.png']);
beta_barplot_file = fullfile(script_dir, ...
    ['lc_oscillator_effective_beta_by_sweep' output_suffix '.png']);

smooth_point_count = 501;
save_combined_parameter_barplot = true;
save_beta_barplot = true;
save_beta_pdf_plot = true;
% Empty uses an automatic bandwidth for both the KDE curve and PDF points.
beta_pdf_bandwidth = [0.15]; % units mA^(-1), the units of beta_eff

% FHI-phase scatter plot. Empty selects gamma min, nearest zero, and gamma
% max. Otherwise enter any nonempty vector of desired gamma values.
% The nearest saved gamma is used; the legend shows its actual Isyn in mA.
representative_gamma_targets = []; % use the five targets saved by ensemble
% Empty uses config.representative_parameter_set. Select another saved set
% to inspect its raw phases and unmodified bits at every gamma.
audit_parameter_set = [];
save_phase_audit_plots = true;

%% Load completed curves
if ~isfile(results_file)
    error(['Missing %s. Run the LC ensemble and FHI-phase-sensitivity ' ...
        'script first.'],results_file);
end

results = load(results_file);
required_fields = {'Pbit_all','completed_parameter_sets','gamma_values', ...
    'parameter_sets','baseline','rel_sigma','sampled_parameter_fields', ...
    'sampled_parameter_names','config','phase_audit', ...
    'Pbit_phiFHI','completed_phiFHI_conditions'};
for field_index = 1:numel(required_fields)
    if ~isfield(results,required_fields{field_index})
        error('Results file is missing %s.',required_fields{field_index});
    end
end
if ~isfield(results.config,'phaseReferenceConvention') || ...
        ~strcmp(results.config.phaseReferenceConvention,'fixed_zero_phase')
    error(['Results use the previous FHI-reference convention. Rerun the ' ...
        'updated ensemble simulation before fitting or plotting.']);
end
if ~isfield(results.config,'referenceFrequencyMode') || ...
        ~strcmp(results.config.referenceFrequencyMode, ...
        expected_reference_frequency_mode)
    error(['The selected results_dataset does not match the reference-' ...
        'frequency mode stored in the results file.']);
end

if ~isfield(results.config,'bitAssignmentConvention') || ...
        ~strcmp(results.config.bitAssignmentConvention, ...
        'independently_calibrated_antipodal_bit1_zero_fhi90_v2') || ...
        ~isfield(results.config,'slopeConvention') || ...
        ~strcmp(results.config.slopeConvention, ...
        'signed_synaptic_FHI_quadrature_v1') || ...
        ~isfield(results.config,'injectionProtocol') || ...
        ~strcmp(results.config.injectionProtocol, ...
        'paper_FHI0_off_before_synaptic_FHI_and_SHI_v1') || ...
        ~isfield(results.config,'lockPhaseDefinition') || ...
        ~strcmp(results.config.lockPhaseDefinition, ...
        'measured_SHI_basin_midpoint_v1') || ...
        ~isfield(results.config,'phiFHI_deg') || ...
        results.config.phiFHI_deg~=90 || ...
        ~isfield(results.config,'bit0Phase_deg') || ...
        results.config.bit0Phase_deg~=180 || ...
        ~isfield(results.config,'bit1Phase_deg') || ...
        results.config.bit1Phase_deg~=0 || ...
        ~isfield(results.config,'shiCoupling') || ...
        ~strcmp(results.config.shiCoupling,'parametric_conductance_v1')
    error(['The selected file is not calibrated two-FHI/parametric-SHI ' ...
        'data. Rerun the updated ensemble.']);
end
% Verify every completed curve against its saved raw phases BEFORE fitting.
for audit_index = find(results.completed_parameter_sets(:)).'
    validate_phase_audit(results.phase_audit.device{audit_index}, ...
        results.Pbit_all(audit_index,:),results.gamma_values, ...
        results.config.j_max);
end
for audit_index = find(results.completed_phiFHI_conditions(:)).'
    validate_phase_audit(results.phase_audit.phiFHI{audit_index}, ...
        results.Pbit_phiFHI(audit_index,:),results.gamma_values, ...
        results.config.j_max);
end

complete_mask = results.completed_parameter_sets(:) & ...
    all(isfinite(results.Pbit_all),2);
parameter_set_numbers = find(complete_mask);
if isempty(parameter_set_numbers)
    error('No complete parameter-set probability curves were found.');
end

Pdata = results.Pbit_all(complete_mask,:);
gamma_data = results.gamma_values(:).';
if ~isfield(results.config,'IfhiUnit') || ...
        ~isscalar(results.config.IfhiUnit) || ...
        ~isfinite(results.config.IfhiUnit) || results.config.IfhiUnit<=0
    error('config.IfhiUnit must be a positive current in amperes per gamma.');
end
Isyn_norm_A = results.config.IfhiUnit;
Isyn_norm_mA = 1e3*Isyn_norm_A;
Isyn_data_mA = Isyn_norm_mA*gamma_data;
current_unit = 'mA';
beta_eff_unit = 'mA^(-1)';
fit_parameterization = '0.5*(1+tanh(Isyn_mA/I0_mA+b))';
parameter_sets_used = results.parameter_sets(parameter_set_numbers);
Nsweep = size(Pdata,1);
% One shared palette for each device's data, tanh fit, and beta bar.
% Index by the saved set number so colors remain consistent if a set is
% incomplete and is omitted from the plots.
sweep_palette = turbo(numel(results.completed_parameter_sets));
sweep_colors = sweep_palette(parameter_set_numbers,:);
if size(Pdata,2)~=numel(gamma_data)
    error('Pbit_all dimensions do not match gamma_values.');
end
if any(Pdata(:)<0 | Pdata(:)>1)
    error('Pbit_all contains probabilities outside [0,1].');
end

% Use the saved Bernoulli outcomes directly in the binomial likelihood.
bit_counts = zeros(Nsweep,numel(gamma_data));
for row = 1:Nsweep
    audit = results.phase_audit.device{parameter_set_numbers(row)};
    bit_counts(row,:) = sum(double(audit.binary_assignment),2).';
end
nTrials = results.config.j_max;

gamma_fit = linspace(min(gamma_data),max(gamma_data),smooth_point_count);
Isyn_fit_mA = Isyn_norm_mA*gamma_fit;

if isempty(representative_gamma_targets)
    representative_gamma_targets = ...
        results.config.representativeGammaTargets;
elseif ~isnumeric(representative_gamma_targets) || ...
        ~isreal(representative_gamma_targets) || ...
        ~isvector(representative_gamma_targets) || ...
        any(~isfinite(representative_gamma_targets(:)))
    error('representative_gamma_targets must be a finite real numeric vector.');
end
representative_gamma_targets = representative_gamma_targets(:).';
if any(representative_gamma_targets<min(gamma_data) | ...
        representative_gamma_targets>max(gamma_data))
    error('Requested gamma values must lie within the saved simulation range.');
end
NrepresentativeGamma = numel(representative_gamma_targets);

[~,representative_gamma_indices] = min(abs( ...
    gamma_data(:)-representative_gamma_targets(:).'),[],1);
if numel(unique(representative_gamma_indices))~=NrepresentativeGamma
    error(['The representative gamma targets must map to ' ...
        'different simulated gamma points.']);
end
representative_gamma_actual = gamma_data(representative_gamma_indices);
if any(abs(representative_gamma_actual-representative_gamma_targets) ...
        >1e-10*max(1,max(abs(gamma_data))))
    warning('Using nearest saved gamma values: %s', ...
        mat2str(representative_gamma_actual,6));
end


%% Fit each varied device independently in physical-current units
% I0_mA is the fitted current scale. A negative I0 retains a decreasing
% response under the existing bit convention; probabilities are never flipped.
I0_mA = nan(Nsweep,1);
beta_eff_mA_inv = nan(Nsweep,1);
bias_b = zeros(Nsweep,1);
Isyn_at_half_mA = nan(Nsweep,1);
gamma_at_half = nan(Nsweep,1);
RMSE = zeros(Nsweep,1);
R_squared = zeros(Nsweep,1);
negative_log_likelihood = zeros(Nsweep,1);
binomial_deviance = zeros(Nsweep,1);
fit_exitflag = zeros(Nsweep,1);
fitted_probability_curves = zeros(Nsweep,smooth_point_count);

fit_options = optimset('Display','off','MaxFunEvals',2e4, ...
    'MaxIter',2e4,'TolX',1e-10,'TolFun',1e-12);
current_scale_mA = max(max(abs(Isyn_data_mA)),eps);

for row = 1:Nsweep
    probability = Pdata(row,:);
    transformed = 2*probability-1;
    clipped = max(min(transformed,1-1e-6),-1+1e-6);
    linearized = atanh(clipped);
    nonsaturated = abs(transformed)<0.95;
    if nnz(nonsaturated)>=2
        initial_line = polyfit(Isyn_data_mA(nonsaturated), ...
            linearized(nonsaturated),1);
    else
        initial_line = polyfit(Isyn_data_mA,linearized,1);
    end
    if any(~isfinite(initial_line))
        initial_line = [1/current_scale_mA 0];
    end
    if abs(initial_line(1))>1e-12/current_scale_mA
        initial_I0_mA = 1/initial_line(1);
    else
        initial_I0_mA = current_scale_mA;
    end
    initial_parameters = [initial_I0_mA initial_line(2)];

    counts = bit_counts(row,:);
    objective = @(parameters) binomial_negative_log_likelihood( ...
        parameters,Isyn_data_mA,counts,nTrials);
    [fitted,negative_log_likelihood(row),fit_exitflag(row)] = ...
        fminsearch(objective,initial_parameters,fit_options);
    if any(~isfinite(fitted)) || fitted(1)==0
        error('Variation %d returned an invalid I0 fit.',parameter_set_numbers(row));
    end
    if fit_exitflag(row)<=0
        warning('CurrentFit:NotConverged', ...
            'I0 fit for variation %d did not converge.',parameter_set_numbers(row));
    end

    I0_mA(row) = fitted(1);
    beta_eff_mA_inv(row) = 1/I0_mA(row);
    bias_b(row) = fitted(2);
    Isyn_at_half_mA(row) = -bias_b(row)*I0_mA(row);
    gamma_at_half(row) = Isyn_at_half_mA(row)/Isyn_norm_mA;

    predicted_at_data = tanh_probability(fitted,Isyn_data_mA);
    residual = probability-predicted_at_data;
    RMSE(row) = sqrt(mean(residual.^2));
    total_variation = sum((probability-mean(probability)).^2);
    if total_variation>0
        R_squared(row) = 1-sum(residual.^2)/total_variation;
    else
        R_squared(row) = NaN;
    end
    fitted_probability_curves(row,:) = ...
        tanh_probability(fitted,Isyn_fit_mA);

    observed_probability = counts/nTrials;
    term_one = zeros(size(counts));
    term_zero = zeros(size(counts));
    positive = counts>0;
    negative = counts<nTrials;
    term_one(positive) = counts(positive).*log(observed_probability(positive));
    term_zero(negative) = (nTrials-counts(negative)).* ...
        log(1-observed_probability(negative));
    saturated_nll = -sum(term_one+term_zero);
    binomial_deviance(row) = max(0,2*( ...
        negative_log_likelihood(row)-saturated_nll));
end
if any(beta_eff_mA_inv<=0)
    warning(['%d fitted curve(s) have nonpositive beta under the saved ' ...
        'fixed-bit convention. Labels were not flipped; inspect the ' ...
        'calibrated phase and endpoint diagnostics.'], ...
        nnz(beta_eff_mA_inv<=0));
end

% Keep explicit SI conversions and the previous normalized-gamma slope.
% These are derived quantities, not additional fit parameters.
I0_A = 1e-3*I0_mA;
beta_eff_A_inv = 1./I0_A;
beta_gamma_inv = Isyn_norm_mA*beta_eff_mA_inv;

%% Effective tanh fit across varied devices, also fitting I0 directly
mean_fitted_probability = mean(fitted_probability_curves,1);
mean_probability_data = mean(Pdata,1);
if Nsweep>1
    standard_error_probability = std(Pdata,0,1)/sqrt(Nsweep);
else
    standard_error_probability = zeros(size(gamma_data));
end

effective_initial_beta = mean(beta_eff_mA_inv);
if abs(effective_initial_beta)>1e-12/current_scale_mA
    effective_initial_I0_mA = 1/effective_initial_beta;
else
    effective_initial_I0_mA = current_scale_mA;
end
effective_initial = [effective_initial_I0_mA mean(bias_b)];
effective_objective = @(parameters) probability_sum_squares( ...
    parameters,Isyn_data_mA,mean_probability_data);
[effective_fit,~,effective_fit_exitflag] = ...
    fminsearch(effective_objective,effective_initial,fit_options);
if any(~isfinite(effective_fit)) || effective_fit(1)==0
    error('The ensemble-effective I0 fit is invalid.');
end
if effective_fit_exitflag<=0
    warning('CurrentFit:EffectiveNotConverged', ...
        'The ensemble-effective I0 fit did not converge.');
end

effective_I0_mA = effective_fit(1);
effective_I0_A = 1e-3*effective_I0_mA;
effective_beta_eff_mA_inv = 1/effective_I0_mA;
effective_beta_eff_A_inv = 1/effective_I0_A;
effective_beta_gamma_inv = Isyn_norm_mA*effective_beta_eff_mA_inv;
effective_bias_b = effective_fit(2);
effective_Isyn_at_half_mA = -effective_bias_b*effective_I0_mA;
effective_gamma_at_half = effective_Isyn_at_half_mA/Isyn_norm_mA;
effective_probability_curve = tanh_probability(effective_fit,Isyn_fit_mA);
effective_RMSE = sqrt(mean((tanh_probability(effective_fit,Isyn_data_mA) ...
    -mean_probability_data).^2));

%% Index of dispersion for the effective-beta calculation
% The effective response is obtained from the ensemble of first-stage
% device fits. Therefore, quantify its beta dispersion across the finite
% device-level beta values that enter that ensemble. The standard index of
% dispersion is the sample variance divided by the sample mean:
%
%   D_beta = var(beta_s)/mean(beta_s).
%
% A single second-stage effective beta is a scalar and does not itself
% define a distribution.
finite_beta_values = beta_eff_mA_inv(isfinite(beta_eff_mA_inv));
if numel(finite_beta_values)>=2
    beta_dispersion_mean = mean(finite_beta_values);
    beta_dispersion_variance = var(finite_beta_values,0);
    beta_mean_tolerance = 1e-12*max(1,max(abs(finite_beta_values)));
    if abs(beta_dispersion_mean)>beta_mean_tolerance
        beta_index_of_dispersion = beta_dispersion_variance/ ...
            beta_dispersion_mean;
    else
        beta_index_of_dispersion = NaN;
    end
else
    beta_dispersion_mean = NaN;
    beta_dispersion_variance = NaN;
    beta_index_of_dispersion = NaN;
end

% Physical-current beta dispersion has units mA^(-1). Its coefficient of
% variation is independent of the choice of current units.
beta_sample_SD = std(finite_beta_values,0);
if isfinite(beta_dispersion_mean) && ...
        abs(beta_dispersion_mean)>beta_mean_tolerance
    beta_CV_percent = 100*beta_sample_SD/abs(beta_dispersion_mean);
else
    beta_CV_percent = NaN;
end
I0_mean_mA = mean(I0_mA);
I0_sample_SD_mA = std(I0_mA,0);

%% Offset spread relative to the simulated input range
% gamma_at_half is the fitted horizontal offset (P(bit=1)=0.5). Its mean
% can be close to zero, so normalize the across-device sample standard
% deviation by the full gamma sweep width, not by the mean offset.
gamma_sweep_width = max(gamma_data)-min(gamma_data);
finite_offsets = gamma_at_half(isfinite(gamma_at_half));
if numel(finite_offsets)>=2
    offset_spread_std = std(finite_offsets,0);
    offset_spread_range = max(finite_offsets)-min(finite_offsets);
else
    offset_spread_std = NaN;
    offset_spread_range = NaN;
end
if gamma_sweep_width>0
    normalized_offset_spread = offset_spread_std/gamma_sweep_width;
    normalized_offset_range = offset_spread_range/gamma_sweep_width;
else
    normalized_offset_spread = NaN;
    normalized_offset_range = NaN;
end

midpoint_current_mean_mA = mean(Isyn_at_half_mA);
midpoint_current_SD_mA = std(Isyn_at_half_mA,0);
Isyn_sweep_width_mA = Isyn_norm_mA*gamma_sweep_width;

%% Primitive parameter deviations
parameter_fields = results.sampled_parameter_fields;
parameter_names = results.sampled_parameter_names;
Nparameter = numel(parameter_fields);
parameter_values = zeros(Nsweep,Nparameter);
baseline_values = zeros(1,Nparameter);
relative_sigmas = zeros(1,Nparameter);

for parameter_index = 1:Nparameter
    field_name = parameter_fields{parameter_index};
    parameter_values(:,parameter_index) = ...
        [parameter_sets_used.(field_name)].';
    baseline_values(parameter_index) = results.baseline.(field_name);
    relative_sigmas(parameter_index) = results.rel_sigma.(field_name);
end

theoretical_sigma = abs(baseline_values).*relative_sigmas;
parameter_z_scores = (parameter_values-baseline_values)./theoretical_sigma;
parameter_percent_deviation = 100.*(parameter_values-baseline_values) ...
    ./baseline_values;

% Derived frequency is included as an explanatory diagnostic, but it is not
% added to the regression because it is determined by L and C.
frequency_GHz = [parameter_sets_used.f0_Hz].'/1e9;
frequency_detuning_percent = 100*(frequency_GHz ...
    /(results.baseline.f0_Hz/1e9)-1);

%% Standardized response sensitivity
% Regression coefficients estimate the change in each standardized fitted
% response caused by a +1 sigma change of a primitive device parameter,
% while holding the other sampled parameters fixed in the linear model.
response_matrix = [beta_eff_mA_inv gamma_at_half bias_b];
response_names = {'Current slope \beta_{eff} (mA^{-1})', ...
    'Half-probability \gamma_{50}', ...
    'Tanh bias b'};
standardized_responses = zeros(size(response_matrix));
for response_index = 1:size(response_matrix,2)
    standardized_responses(:,response_index) = ...
        standardize_finite(response_matrix(:,response_index));
end

regression_coefficients = nan(Nparameter,size(response_matrix,2));
regression_R_squared = nan(1,size(response_matrix,2));
design_matrix = [ones(Nsweep,1) parameter_z_scores];
ridge_strength = 1e-8;
penalty = diag([0 ones(1,Nparameter)]);

for response_index = 1:size(response_matrix,2)
    response = standardized_responses(:,response_index);
    valid = all(isfinite(design_matrix),2) & isfinite(response);
    if nnz(valid)>=Nparameter+1
        coefficients = (design_matrix(valid,:).'*design_matrix(valid,:) ...
            +ridge_strength*penalty) ...
            \(design_matrix(valid,:).'*response(valid));
        regression_coefficients(:,response_index) = coefficients(2:end);
        predicted = design_matrix(valid,:)*coefficients;
        residual_sum = sum((response(valid)-predicted).^2);
        total_sum = sum((response(valid)-mean(response(valid))).^2);
        if total_sum>0
            regression_R_squared(response_index) = 1-residual_sum/total_sum;
        end
    else
        warning(['Only %d complete sets are available. At least %d are ' ...
            'needed for the full five-parameter sensitivity regression.'], ...
            nnz(valid),Nparameter+1);
    end
end


%% Tables, report, and saved numerical results
ParameterSet = parameter_set_numbers;
fit_results = table(ParameterSet,I0_mA,beta_eff_mA_inv,I0_A,beta_eff_A_inv, ...
    beta_gamma_inv,bias_b,Isyn_at_half_mA,gamma_at_half, ...
    RMSE,R_squared,negative_log_likelihood,binomial_deviance,fit_exitflag, ...
    frequency_GHz,frequency_detuning_percent);
disp(fit_results);

effective_fit_results = table(effective_I0_mA,effective_beta_eff_mA_inv, ...
    effective_I0_A,effective_beta_eff_A_inv,effective_beta_gamma_inv, ...
    effective_bias_b,effective_Isyn_at_half_mA,effective_gamma_at_half, ...
    effective_RMSE,beta_dispersion_mean,beta_sample_SD,beta_CV_percent, ...
    beta_dispersion_variance,beta_index_of_dispersion, ...
    midpoint_current_SD_mA,offset_spread_std,normalized_offset_spread, ...
    offset_spread_range,normalized_offset_range,gamma_sweep_width,Nsweep);
disp('Effective current-scale tanh fit across completed LC parameter sets:');
disp(effective_fit_results);

fprintf('\nFit model: P(1)=0.5*[1+tanh(Isyn_mA/I0_mA+b)].\n');
fprintf('I0 is fitted in mA; beta_eff=1/I0 is reported in mA^(-1).\n');
fprintf('Saved current scale: Isyn = gamma * %.8g mA.\n',Isyn_norm_mA);
fprintf('Mean I0 = %.8g mA; sample SD = %.8g mA.\n', ...
    I0_mean_mA,I0_sample_SD_mA);
fprintf('Mean beta_eff = %.8g mA^(-1); sample SD = %.8g mA^(-1).\n', ...
    beta_dispersion_mean,beta_sample_SD);
fprintf('Beta coefficient of variation = %.8g percent.\n',beta_CV_percent);
fprintf('D_beta = variance/mean = %.8g mA^(-1) (unit dependent).\n', ...
    beta_index_of_dispersion);
fprintf('Ensemble-effective I0 = %.8g mA; beta_eff = %.8g mA^(-1).\n', ...
    effective_I0_mA,effective_beta_eff_mA_inv);
fprintf('\nMidpoint Isyn,50 = -b*I0.\n');
fprintf('Mean midpoint current = %.8g mA; sample SD = %.8g mA.\n', ...
    midpoint_current_mean_mA,midpoint_current_SD_mA);
fprintf('SD(gamma_50) = %.8g; SD/full sweep width = %.8g.\n', ...
    offset_spread_std,normalized_offset_spread);

Parameter = parameter_names(:);
BetaSensitivity = regression_coefficients(:,1);
Gamma50Sensitivity = regression_coefficients(:,2);
BiasSensitivity = regression_coefficients(:,3);
sensitivity_results = table(Parameter,BetaSensitivity, ...
    Gamma50Sensitivity,BiasSensitivity);
disp('Standardized multivariable sensitivity coefficients:');
disp(sensitivity_results);

writetable(fit_results,fit_table_file);
writetable(sensitivity_results,sensitivity_table_file);

save(fit_file,'fit_results','effective_fit_results', ...
    'I0_mA','I0_A','beta_eff_mA_inv','beta_eff_A_inv','beta_gamma_inv', ...
    'bias_b','Isyn_at_half_mA','gamma_at_half', ...
    'effective_I0_mA','effective_I0_A','effective_beta_eff_mA_inv', ...
    'effective_beta_eff_A_inv','effective_beta_gamma_inv', ...
    'effective_bias_b','effective_Isyn_at_half_mA','effective_gamma_at_half', ...
    'effective_fit_exitflag','Isyn_data_mA','Isyn_fit_mA', ...
    'Isyn_norm_A','Isyn_norm_mA','current_unit','beta_eff_unit', ...
    'fit_parameterization','I0_mean_mA','I0_sample_SD_mA', ...
    'beta_sample_SD','beta_CV_percent', ...
    'midpoint_current_mean_mA','midpoint_current_SD_mA', ...
    'Isyn_sweep_width_mA', ...
    'sensitivity_results','gamma_data','Pdata','bit_counts','nTrials', ...
    'gamma_fit','fitted_probability_curves','mean_fitted_probability', ...
    'mean_probability_data','standard_error_probability', ...
    'effective_probability_curve', ...
    'beta_index_of_dispersion','beta_dispersion_mean', ...
    'beta_dispersion_variance','finite_beta_values', ...
    'gamma_sweep_width','finite_offsets','offset_spread_std', ...
    'offset_spread_range','normalized_offset_spread','normalized_offset_range', ...
    'representative_gamma_targets','representative_gamma_indices', ...
    'representative_gamma_actual', ...
    'parameter_set_numbers','parameter_sets_used','sweep_colors', ...
    'parameter_fields','parameter_names','parameter_values','parameter_z_scores', ...
    'parameter_percent_deviation','baseline_values','relative_sigmas', ...
    'regression_coefficients','regression_R_squared', ...
    'standardized_responses','frequency_GHz','frequency_detuning_percent', ...
    'results','-v7.3');


%% Figure 1: individual probability data and current-scale tanh fits
figure('Color','w','Name','LC current-scale tanh fits','NumberTitle','off');
hold on;
fit_handles = gobjects(Nsweep,1);
for row = 1:Nsweep
    color = sweep_colors(row,:);
    scatter(Isyn_data_mA,Pdata(row,:),34,color,'filled', ...
        'MarkerEdgeColor',color,'HandleVisibility','off');
    fit_handles(row) = plot(Isyn_fit_mA,fitted_probability_curves(row,:), ...
        'Color',color,'LineWidth',1.8,'DisplayName',sprintf( ...
        'Variation %d: I_0=%.3g mA, \\beta_{eff}=%.3g mA^{-1}, b=%+.3g', ...
        parameter_set_numbers(row),I0_mA(row),beta_eff_mA_inv(row),bias_b(row)));
end
xlabel('Synaptic FHI current I_{syn} (mA)'); ylabel('p = P(bit=1)');
title('$P(1)=\frac{1}{2}[1+\tanh(I_{\mathrm{syn}}/I_0+b)]$', ...
    'Interpreter','latex');
legend(fit_handles,'Location','eastoutside','FontSize',8);
xlim([min(Isyn_data_mA) max(Isyn_data_mA)]); ylim([-0.05 1.05]);
set(gca,'FontSize',13,'LineWidth',1.5,'TickDir','in');
box on; grid off;

%% Figure 2: ensemble-effective current response and data-based SEM
figure('Color','w','Name','Effective LC current-scale tanh response', ...
    'NumberTitle','off');
hold on;
for row = 1:Nsweep
    plot(Isyn_fit_mA,fitted_probability_curves(row,:),'--', ...
        'Color',[0.76 0.76 0.76],'LineWidth',0.9,'HandleVisibility','off');
end
errorbar(Isyn_data_mA,mean_probability_data,standard_error_probability, ...
    'o','Color',[0.10 0.35 0.75],'MarkerFaceColor',[0.82 0.89 1.00], ...
    'LineWidth',1.3,'CapSize',7,'DisplayName','Observed mean \pm SEM');
plot(Isyn_fit_mA,effective_probability_curve,'b-','LineWidth',3, ...
    'DisplayName',sprintf( ...
    'Effective: I_0=%.4g mA, \\beta_{eff}=%.4g mA^{-1}, b=%+.4g', ...
    effective_I0_mA,effective_beta_eff_mA_inv,effective_bias_b));
xlabel('Synaptic FHI current I_{syn} (mA)'); ylabel('p = P(bit=1)');
title(sprintf(['Parametric-SHI ensemble response from %d varied devices; ' ...
    '%d trials/point'],Nsweep,nTrials));
legend('Location','best');
xlim([min(Isyn_data_mA) max(Isyn_data_mA)]); ylim([-0.05 1.05]);
set(gca,'FontSize',14,'LineWidth',1.6,'TickDir','in');
box on; grid off;

%% Figure 3: sensitivity to FHI phase
if isfield(results,'Pbit_phiFHI') && ...
        isfield(results,'completed_phiFHI_conditions') && ...
        isfield(results.config,'phiFHI_sensitivity_deg')
    phiFHI_mask = results.completed_phiFHI_conditions(:) & ...
        all(isfinite(results.Pbit_phiFHI),2);
    if any(phiFHI_mask)
        phiFHI_probability = results.Pbit_phiFHI(phiFHI_mask,:);
        phiFHI_values = results.config.phiFHI_sensitivity_deg(phiFHI_mask);
        %% Probability versus measured lock-phase target for three inputs
        % Remove gamma = +/-0.48 from this plot only. Keep the original
        % representative inputs for all other diagnostics and stored results.
        sensitivity_gamma_tolerance = ...
            1e-10*max(1,max(abs(gamma_data)));
        sensitivity_representative_indices = find( ...
            abs(abs(representative_gamma_actual)-0.48) ...
            >sensitivity_gamma_tolerance);
        NsensitivityGamma = numel(sensitivity_representative_indices);
        figure('Color','w','Name', ...
            'Probability versus FHI initialization phase', ...
            'NumberTitle','off');
        hold on;
        gamma_colors = lines(NrepresentativeGamma);
        gamma_handles = gobjects(NsensitivityGamma,1);
        for sensitivity_index = 1:NsensitivityGamma
            representative_index = ...
                sensitivity_representative_indices(sensitivity_index);
            gamma_index = representative_gamma_indices( ...
                representative_index);
            gamma_handles(sensitivity_index) = plot( ...
                phiFHI_values,phiFHI_probability(:,gamma_index),'o-', ...
                'Color',gamma_colors(representative_index,:), ...
                'MarkerFaceColor',gamma_colors(representative_index,:), ...
                'MarkerSize',5,'LineWidth',1.5, ...
                'DisplayName',sprintf('I_{syn} = %+.4g mA', ...
                Isyn_data_mA(gamma_index)));
        end
        xlabel('Initialization phase relative to bit 1 (deg)');
        ylabel('p = P(bit=1)');
        title(sprintf(['Initialization-phase probability for ' ...
            '%d synaptic inputs, parameter set %d'], ...
            NsensitivityGamma,results.config.representative_parameter_set));
        legend(gamma_handles,'Location','best');
        xticks(sort(unique(phiFHI_values)));
        if numel(phiFHI_values)>1
            phase_margin = max(2,0.05*(max(phiFHI_values)-min(phiFHI_values)));
            xlim([min(phiFHI_values)-phase_margin ...
                max(phiFHI_values)+phase_margin]);
        end
        ylim([-0.05 1.05]);
        set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','in');
        grid off; box on;
    end
end

%% Figure 4: fitted response parameters by device set
figure('Color','w','Name','LC fitted-response parameters', ...
    'NumberTitle','off');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile;
beta_summary_bars = bar(parameter_set_numbers,beta_eff_mA_inv, ...
    'FaceColor','flat','EdgeColor','none');
beta_summary_bars.CData = sweep_colors;
ylabel('\beta_{eff} (mA^{-1})'); title('Inverse fitted current scale'); box on;
nexttile;
bar(parameter_set_numbers,Isyn_at_half_mA,'FaceColor',[0.90 0.45 0.20]);
ylabel('I_{syn,50} (mA)'); title('Half-probability current'); box on;
nexttile;
bar(parameter_set_numbers,RMSE,'FaceColor',[0.35 0.70 0.40]);
xlabel('Parameter-set number'); ylabel('RMSE'); title('Fit error'); box on;

%% Dedicated bar plot: beta by variation, colored to match the tanh fits
figBetaBySweep = figure('Color','w', ...
    'Name','Effective beta by device variation', ...
    'NumberTitle','off','Units','pixels','Position',[80 80 560 420], ...
    'PaperUnits','inches','PaperPosition',[0 0 5.6 4.2], ...
    'PaperSize',[5.6 4.2],'PaperPositionMode','manual');
beta_axis = axes('Parent',figBetaBySweep, ...
    'Position',[0.13 0.25 0.84 0.70]);
sweep_positions = 1:Nsweep;
beta_bars = bar(beta_axis,sweep_positions,beta_eff_mA_inv,0.80, ...
    'FaceColor','flat','EdgeColor','none');
beta_bars.CData = sweep_colors;
set(beta_axis,'XTick',sweep_positions, ...
    'XTickLabel',compose('Variation %d',parameter_set_numbers), ...
    'XTickLabelRotation',90,'FontSize',8,'LineWidth',1.0, ...
    'TickDir','in','TickLength',[0.010 0.025],'Layer','top');
beta_axis.XAxis.FontSize = 6; % fit all 50 rotated variation labels
xlim(beta_axis,[0.25 Nsweep+0.75]);
beta_lower = min(0,min(beta_eff_mA_inv));
beta_upper = max(0,max(beta_eff_mA_inv));
beta_padding = 0.12*max(1,beta_upper-beta_lower);
if beta_lower<0
    beta_lower = beta_lower-beta_padding;
end
ylim(beta_axis,[beta_lower beta_upper+beta_padding]);
xlabel(beta_axis,'Device parameter set (variation number)','FontSize',9);
ylabel(beta_axis,'Effective \beta_{eff} (mA^{-1})','FontSize',9);
grid(beta_axis,'off'); box(beta_axis,'on');
if save_beta_barplot
    % Fixed paper dimensions at 100 dpi produce a 560-by-420-pixel PNG.
    print(figBetaBySweep,beta_barplot_file,'-dpng','-r100');
end


%% Probability-density data points and smooth KDE curve
% A finite set of beta samples has no unique pointwise empirical PDF.
% Estimate density on a smooth grid and at each observed device beta by KDE.
% The red curve and colored markers use the same KDE bandwidth. No separate
% Gaussian distribution is fitted, and no histogram is displayed.
if isempty(beta_eff_mA_inv) || any(~isfinite(beta_eff_mA_inv)) || Nsweep<2
    warning('BetaPDF:Undefined', ...
        'The beta PDF requires at least two finite device-level beta values.');
else
    [sorted_beta_eff,beta_pdf_order] = sort(beta_eff_mA_inv(:),'ascend');
    beta_pdf_sample_count = numel(sorted_beta_eff);
    beta_pdf_estimation_method = ...
        'Gaussian KDE evaluated on a smooth grid and at device beta values';

    if isempty(beta_pdf_bandwidth)
        % Silverman's rule with robust scale; no Statistics Toolbox required.
        beta_pdf_SD = std(sorted_beta_eff,0);
        quartile_positions = 1+(beta_pdf_sample_count-1)*[0.25 0.75];
        beta_pdf_quartiles = interp1(1:beta_pdf_sample_count, ...
            sorted_beta_eff,quartile_positions,'linear');
        beta_pdf_IQR_scale = diff(beta_pdf_quartiles)/1.34;
        if beta_pdf_IQR_scale>0
            beta_pdf_scale = min(beta_pdf_SD,beta_pdf_IQR_scale);
        else
            beta_pdf_scale = beta_pdf_SD;
        end
        if ~isfinite(beta_pdf_scale) || beta_pdf_scale<=0
            beta_pdf_scale = 0.01*max(1,abs(mean(sorted_beta_eff)));
            warning('BetaPDF:ZeroSpread', ...
                ['All beta values coincide. The displayed KDE uses a finite ' ...
                 'visualization bandwidth, not a measured spread.']);
        end
        beta_pdf_bandwidth_used = ...
            0.9*beta_pdf_scale*beta_pdf_sample_count^(-1/5);
    else
        validateattributes(beta_pdf_bandwidth,{'numeric'}, ...
            {'scalar','real','finite','positive'});
        beta_pdf_bandwidth_used = beta_pdf_bandwidth;
    end

    beta_pdf_limits = [min(sorted_beta_eff)-4*beta_pdf_bandwidth_used ...
        max(sorted_beta_eff)+4*beta_pdf_bandwidth_used];
    beta_pdf_x = linspace(beta_pdf_limits(1),beta_pdf_limits(2),501);
    beta_pdf_density = gaussian_kernel_density( ...
        sorted_beta_eff,beta_pdf_x,beta_pdf_bandwidth_used);
    beta_pdf_at_devices = gaussian_kernel_density( ...
        sorted_beta_eff,sorted_beta_eff,beta_pdf_bandwidth_used);

    figBetaPDF = figure('Color','w','Name', ...
        'Effective-beta PDF data and KDE curve','NumberTitle','off', ...
        'Units','pixels','Position',[90 90 560 420]);
    beta_pdf_axes = axes('Parent',figBetaPDF);
    hold(beta_pdf_axes,'on');
    beta_pdf_curve_handle = plot(beta_pdf_axes,beta_pdf_x, ...
        beta_pdf_density,'r-','LineWidth',2.4, ...
        'DisplayName','Kernel density estimate');
    beta_pdf_scatter_handle = scatter(beta_pdf_axes,sorted_beta_eff, ...
        beta_pdf_at_devices,50,sweep_colors(beta_pdf_order,:), ...
        'filled','MarkerEdgeColor','flat','LineWidth',0.5, ...
        'DisplayName','PDF data points (KDE)');
    legend(beta_pdf_axes,[beta_pdf_scatter_handle beta_pdf_curve_handle], ...
        'Location','best','Box','off');
    xlabel(beta_pdf_axes,'Effective \beta_{eff} (mA^{-1})', ...
        'FontSize',15,'Interpreter','tex');
    % A density with respect to beta [mA^-1] has reciprocal units [mA].
    ylabel(beta_pdf_axes,'Probability density (mA)','FontSize',15);
    set(beta_pdf_axes,'FontSize',13,'LineWidth',1.4,'TickDir','in', ...
        'XColor','k','YColor','k');
    xlim(beta_pdf_axes,beta_pdf_limits);
    beta_pdf_max = max([beta_pdf_at_devices(:);beta_pdf_density(:)]);
    ylim(beta_pdf_axes,[0 1.12*beta_pdf_max]);
    box(beta_pdf_axes,'on'); grid(beta_pdf_axes,'off');

    beta_pdf_plot_file = fullfile(script_dir, ...
        ['lc_oscillator_effective_beta_probability_density' output_suffix '.png']);
    if save_beta_pdf_plot
        drawnow;
        exportgraphics(figBetaPDF,beta_pdf_plot_file,'Resolution',300);
        savefig(figBetaPDF,fullfile(script_dir, ...
            ['lc_oscillator_effective_beta_probability_density' output_suffix '.fig']));
        fprintf('Saved effective-beta PDF plot to:\n%s\n',beta_pdf_plot_file);
    end
    save(fit_file,'sorted_beta_eff','beta_pdf_order','beta_pdf_x', ...
        'beta_pdf_at_devices','beta_pdf_bandwidth_used','beta_pdf_sample_count', ...
        'beta_pdf_estimation_method','beta_pdf_density','-append');
end

%% Figure 5: standardized device-parameter sensitivity
figure('Color','w','Name','LC standardized parameter sensitivity', ...
    'NumberTitle','off');
bar(regression_coefficients,'grouped');
yline(0,'k-');
set(gca,'XTick',1:Nparameter,'XTickLabel',parameter_names, ...
    'FontSize',12,'LineWidth',1.3,'TickDir','in');
ylabel('Standardized response change per +1\sigma parameter change');
title(sprintf(['Multivariable linear sensitivity; R^2 = ' ...
    '[%.2f, %.2f, %.2f]'],regression_R_squared));
legend(response_names,'Location','best');
grid off; box on;

%% Figure 6: accumulated device-parameter z-score bar plot
figAllParameterBars = figure('Color','w', ...
    'Name','All LC parameter deviations','NumberTitle','off');

barplot_columns = min(5,Nsweep);
barplot_rows = ceil(Nsweep/barplot_columns);
barplot_layout = tiledlayout(barplot_rows,barplot_columns, ...
    'TileSpacing','compact','Padding','compact');
parameter_colors = turbo(Nparameter);

maximum_z = max(abs(parameter_z_scores(:)));
if ~isfinite(maximum_z) || maximum_z==0
    maximum_z = 1;
end
z_limit = max(2,ceil(maximum_z));
z_limits = [-z_limit z_limit];
z_ticks = unique([-z_limit 0 z_limit]);

first_bar_axis = gobjects(1);
for row = 1:Nsweep
    bar_axis = nexttile(barplot_layout);
    if row==1
        first_bar_axis = bar_axis;
    end
    bar_handle = bar(bar_axis,1:Nparameter, ...
        parameter_z_scores(row,:),0.78, ...
        'FaceColor','flat','EdgeColor','none');
    bar_handle.CData = parameter_colors;
    hold(bar_axis,'on');
    yline(bar_axis,0,'k-','LineWidth',0.7,'HandleVisibility','off');
    hold(bar_axis,'off');
    set(bar_axis,'XTick',[],'YTick',z_ticks,'YLim',z_limits, ...
        'XLim',[0.5 Nparameter+0.5],'FontSize',11, ...
        'LineWidth',1.1,'TickDir','in');
    text(bar_axis,0.04,0.82,sprintf('Variation %d', ...
        parameter_set_numbers(row)), ...
        'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','middle','FontSize',8,'FontWeight','bold', ...
        'BackgroundColor','w','Margin',0.5,'Clipping','on');
    grid(bar_axis,'off');
    box(bar_axis,'on');
end

z_score_label = ylabel(barplot_layout,'Z-score','FontSize',16);
z_score_label.FontWeight = 'bold';

legend_handles = gobjects(Nparameter,1);
hold(first_bar_axis,'on');
for parameter_index = 1:Nparameter
    legend_handles(parameter_index) = plot(first_bar_axis,nan,nan,'s', ...
        'MarkerSize',7,'MarkerFaceColor',parameter_colors(parameter_index,:), ...
        'MarkerEdgeColor','none','LineStyle','none');
end
hold(first_bar_axis,'off');
parameter_legend = legend(first_bar_axis,legend_handles,parameter_names, ...
    'Orientation','horizontal','NumColumns',Nparameter, ...
    'FontSize',9,'Box','off');
parameter_legend.Layout.Tile = 'south';

if save_combined_parameter_barplot
    saveas(figAllParameterBars,combined_barplot_file);
end

%% Figure 7: representative tank-voltage waveform, when available
if isfield(results,'representative_trace') && ...
        isfield(results.representative_trace,'stored') && ...
        results.representative_trace.stored
    trace = results.representative_trace;
    if ~isfield(trace,'voltage') || isempty(trace.voltage)
        warning(['The saved representative trace does not contain v(t). ' ...
            'Rerun the updated ensemble simulation to regenerate it.']);
    else
        % The phase reference is fixed at zero phase. The injected FHI
        % current in the separate current figure contains the calibrated
        % drive-phase shift, which differs from the voltage-phase target.
        if isfield(trace,'reference_omega_rad_s') && ...
                isscalar(trace.reference_omega_rad_s) && ...
                isfinite(trace.reference_omega_rad_s) && ...
                trace.reference_omega_rad_s>0
            plotted_reference_omega = trace.reference_omega_rad_s;
            plotted_reference_frequency = trace.reference_frequency_Hz;
        else
            plotted_reference_omega = results.baseline.omega0;
            plotted_reference_frequency = results.baseline.f0_Hz;
        end
        if isfield(trace,'reference_start_time_s') && ...
                isscalar(trace.reference_start_time_s) && ...
                isfinite(trace.reference_start_time_s)
            plotted_reference_start_time = trace.reference_start_time_s;
        else
            plotted_reference_start_time = 0;
        end
        reference_time = trace.time-plotted_reference_start_time;
        reference_is_active = reference_time>=0;
        fhi_reference_signal = zeros(size(trace.time));
        fhi_reference_signal(reference_is_active) = sin( ...
            plotted_reference_omega*reference_time(reference_is_active));
        measurement_window_start = trace.time(end) ...
            -results.config.phaseWindowCycles/plotted_reference_frequency;
        measurement_mask = trace.time>=measurement_window_start;
        if nnz(measurement_mask)>=2
            reference_plot_amplitude = 0.5*( ...
                max(trace.voltage(measurement_mask)) ...
                -min(trace.voltage(measurement_mask)));
        else
            reference_plot_amplitude = 0.5*( ...
                max(trace.voltage)-min(trace.voltage));
        end
        if ~isfinite(reference_plot_amplitude) || ...
                reference_plot_amplitude<=0
            reference_plot_amplitude = 1;
        end
        fhi_reference_voltage = ...
            reference_plot_amplitude*fhi_reference_signal;

        %% Representative tank-voltage/reference figure
        figure('Color','w','Name', ...
            'Representative v(t) and FHI reference','NumberTitle','off');
        plot(trace.normalized_time,trace.voltage,'LineWidth',1.3);
        hold on;
        plot(trace.normalized_time,fhi_reference_voltage,'k:', ...
            'LineWidth',1.5);
        xlabel('Reference cycles');
        ylabel('Tank voltage v(t) (V)');
        legend('v(t)','FHI reference (scaled)','Location','best');
        title(sprintf('Parameter set %d, \\gamma=%+.3g', ...
            trace.parameter_set_index,trace.gamma));
        xlim([trace.normalized_time(1) trace.normalized_time(end)]);
        grid off; box on;

        %% Separate locking FHI, synaptic FHI and SHI currents
        figure('Color','w','Name', ...
            'Representative locking FHI, synaptic FHI and SHI currents', ...
            'NumberTitle','off');
        plot(trace.normalized_time,trace.fhi_lock_current*1e3, ...
            'LineWidth',1.3);
        hold on;
        plot(trace.normalized_time,trace.fhi_synaptic_current*1e3, ...
            'LineWidth',1.3);
        plot(trace.normalized_time,trace.shi_current*1e3,'--', ...
            'LineWidth',1.3);
        xlabel('Reference cycles');
        ylabel('Current (mA)');
        legend('Lock FHI','Synaptic FHI','SHI','Location','best');
        xlim([trace.normalized_time(1) trace.normalized_time(end)]);
        grid off; box on;
        title(sprintf(['Parameter set %d, \\gamma=%+.3g; ' ...
            'FHI is removed before readout'], ...
            trace.parameter_set_index,trace.gamma));
    end
end

%% Audit figures: calibrated phases, axes, amplitudes, and outcomes
if isempty(audit_parameter_set)
    audit_parameter_set = results.config.representative_parameter_set;
end
if ~isscalar(audit_parameter_set) || ...
        ~ismember(audit_parameter_set,parameter_set_numbers)
    error('audit_parameter_set must be a completed parameter-set number.');
end
audit_data = results.phase_audit.device{audit_parameter_set};
raw_phase = squeeze(audit_data.raw_relative_phase_deg);
calibrated_phase = squeeze(audit_data.calibrated_relative_phase_deg);
raw_bits = squeeze(audit_data.binary_assignment);
fundamental_amplitude = squeeze(audit_data.fundamental_amplitude_V);
low_amplitude = squeeze(audit_data.low_amplitude);
audit_prefix = fullfile(script_dir,sprintf( ...
    'lc_parametric_phase_audit_set_%d%s', ...
    audit_parameter_set,output_suffix));
bit_colors = [0.10 0.45 0.80; 0.90 0.35 0.10];

% Classification is performed only after subtracting the independently
% calibrated SHI-only bit-1 axis. Every saved trajectory remains visible.
figCalibratedPhase = figure('Color','w', ...
    'Name','Calibrated parametric-SHI phases', ...
    'NumberTitle','off','Position',[80 80 1050 520]);
hold on;
gamma_matrix = repmat(gamma_data(:),1,size(calibrated_phase,2));
for bit_value = 0:1
    mask = raw_bits==bit_value;
    scatter(gamma_matrix(mask),calibrated_phase(mask),12, ...
        bit_colors(bit_value+1,:),'filled','MarkerFaceAlpha',0.35, ...
        'DisplayName',sprintf('Saved bit %d',bit_value));
end
draw_phase_guides();
xlabel('\gamma'); ylabel('Calibrated relative phase (deg)');
x_margin = 0.04*max(1,max(gamma_data)-min(gamma_data));
xlim([min(gamma_data)-x_margin max(gamma_data)+x_margin]);
title(sprintf(['Device %d: parametric-SHI phase after subtracting ' ...
    'the %.2f-deg calibration axis'], ...
    audit_parameter_set,audit_data.fixed_axis_deg));
legend('Location','eastoutside'); grid off; box on;

% Matched columns retain trajectory identity and expose the actual bit.
figCalibratedBits = figure('Color','w', ...
    'Name','Calibrated phase and binary assignment', ...
    'NumberTitle','off','Position', ...
    [100 100 min(1900,max(700,380*NrepresentativeGamma)) 720]);
bit_layout = tiledlayout(2,NrepresentativeGamma, ...
    'TileSpacing','compact','Padding','compact');
for column = 1:NrepresentativeGamma
    gamma_index = representative_gamma_indices(column);
    phases = calibrated_phase(gamma_index,:);
    bits = raw_bits(gamma_index,:);
    trajectory = 1:numel(bits);
    nexttile(bit_layout,column); hold on;
    scatter(trajectory,phases,18,bit_colors(double(bits)+1,:),'filled');
    draw_phase_guides();
    xlim([0.5 numel(bits)+0.5]);
    xlabel('Trajectory'); ylabel('Calibrated phase (deg)');
    title(sprintf('\\gamma=%+.3g; calibrated axis=0 deg', ...
        gamma_data(gamma_index))); grid off; box on;
    nexttile(bit_layout,NrepresentativeGamma+column); hold on;
    scatter(trajectory,double(bits),18, ...
        bit_colors(double(bits)+1,:),'filled');
    ylim([-0.15 1.15]); yticks([0 1]);
    xlim([0.5 numel(bits)+0.5]);
    xlabel('Trajectory'); ylabel('Saved binary assignment');
    title(sprintf('P(1)=%d/%d=%.4f',sum(bits),numel(bits),mean(bits)));
    grid off; box on;
end
title(bit_layout,sprintf(['Device %d: independently calibrated ' ...
    'antipodal classification'],audit_parameter_set));

% Physical axes are reported in the original zero-phase reference frame.
figPhaseAxis = figure('Color','w', ...
    'Name','Parametric-SHI phase-axis diagnostics', ...
    'NumberTitle','off','Position',[120 80 1000 900]);
axis_layout = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
plot(gamma_data,audit_data.fixed_axis_deg*ones(size(gamma_data)), ...
    'k-','LineWidth',1.8,'DisplayName','SHI-only calibrated bit-1 axis');
plot(gamma_data,audit_data.observed_axis_deg,'o','LineStyle','none', ...
    'DisplayName','Observed doubled-phase axis (diagnostic only)');
ylim([-95 95]); yticks(-90:45:90);
ylabel('Raw axis (deg, modulo 180)'); legend('Location','best');
grid off; box on;
nexttile;
plot(gamma_data,audit_data.doubled_phase_concentration,'o-');
ylim([0 1.05]); ylabel('Axial concentration');
title('Near 1: tight antipodal clustering; near 0: weakly defined axis');
grid off; box on;
nexttile;
plot(gamma_data,mean(audit_data.near_boundary,2),'o-');
ylim([0 1.05]); ylabel('Boundary fraction');
title(sprintf('Within %.1f deg of calibrated +/-90-deg boundaries', ...
    audit_data.boundary_warning_width_deg));
grid off; box on;
nexttile; hold on;
plot(gamma_data,mean(fundamental_amplitude,2),'o-', ...
    'DisplayName','Mean fundamental amplitude');
yline(audit_data.minimum_amplitude_V,'r--', ...
    'DisplayName','Minimum accepted amplitude');
xlabel('\gamma'); ylabel('Amplitude (V)');
title(sprintf('Low-amplitude fraction: max %.3g', ...
    max(mean(low_amplitude,2))));
legend('Location','best'); grid off; box on;
title(axis_layout,sprintf('Device %d: parametric-SHI readout diagnostics', ...
    audit_parameter_set));

% Show the deterministic multi-initial-condition test that defines the
% device-specific antipodal axis before any noisy probability curve.
calibration = audit_data.two_state_calibration;
figCalibration = figure('Color','w', ...
    'Name','Parametric-SHI two-state calibration', ...
    'NumberTitle','off','Position',[140 140 900 520]);
scatter(calibration.initial_phases_deg,calibration.final_phases_deg,70, ...
    calibration.initial_phases_deg,'filled','MarkerEdgeColor','k');
hold on;
yline(calibration.bit0_axis_deg,'k-','LineWidth',1.6, ...
    'DisplayName','Calibrated bit-0 axis');
yline(calibration.bit1_axis_deg,'k--','LineWidth',1.6, ...
    'DisplayName','Calibrated bit-1 axis');
xlabel('Initial oscillator phase (deg)');
ylabel('Final SHI-locked phase (deg)');
title(sprintf(['Deterministic SHI-only calibration: separation %.2f deg, ' ...
    'maximum drift %.2f deg'],calibration.cluster_separation_deg, ...
    calibration.maximum_window_drift_deg));
legend('Location','best'); ylim([-190 190]); yticks(-180:90:180);
grid off; box on;

% A row per trajectory permits direct reproduction of the plotted P(1).
[GammaIndex,Trajectory] = ndgrid(1:numel(gamma_data), ...
    1:size(calibrated_phase,2));
Gamma = gamma_matrix(:);
RawRelativePhase_deg = raw_phase(:);
CalibratedRelativePhase_deg = calibrated_phase(:);
FundamentalAmplitude_V = fundamental_amplitude(:);
LowAmplitude = logical(low_amplitude(:));
Bit = double(raw_bits(:));
NearBoundary = logical(audit_data.near_boundary(:));
FixedAxis_deg = audit_data.fixed_axis_deg*ones(numel(Bit),1);
phase_audit_table = table(GammaIndex(:),Gamma,Trajectory(:), ...
    RawRelativePhase_deg,CalibratedRelativePhase_deg, ...
    FundamentalAmplitude_V,LowAmplitude,FixedAxis_deg,Bit,NearBoundary);
writetable(phase_audit_table,[audit_prefix '_trajectories.csv']);

phase_axis_table = table(gamma_data(:), ...
    audit_data.fixed_axis_deg*ones(numel(gamma_data),1), ...
    audit_data.observed_axis_deg(:), ...
    audit_data.doubled_phase_concentration(:), ...
    mean(audit_data.near_boundary,2), ...
    mean(fundamental_amplitude,2),mean(low_amplitude,2), ...
    mean(raw_bits,2), ...
    'VariableNames',{'Gamma','CalibrationAxis_deg','ObservedAxis_deg', ...
    'DoubledPhaseConcentration','BoundaryFraction', ...
    'MeanFundamentalAmplitude_V','LowAmplitudeFraction','Pbit1'});
writetable(phase_axis_table,[audit_prefix '_axes.csv']);
if save_phase_audit_plots
    audit_figures = [figCalibratedPhase figCalibratedBits ...
        figPhaseAxis figCalibration];
    audit_names = {'_calibrated_phases','_binary_assignment', ...
        '_phase_axis_and_amplitude','_two_state_calibration'};
    for audit_index = 1:numel(audit_figures)
        exportgraphics(audit_figures(audit_index), ...
            [audit_prefix audit_names{audit_index} '.png'],'Resolution',200);
        savefig(audit_figures(audit_index), ...
            [audit_prefix audit_names{audit_index} '.fig']);
    end
end
save(fit_file,'audit_parameter_set','phase_audit_table', ...
    'phase_axis_table','-append');

fprintf('\nSaved fit results to:\n%s\n',fit_file);
fprintf('Saved fit table to:\n%s\n',fit_table_file);
fprintf('Saved standardized sensitivity table to:\n%s\n', ...
    sensitivity_table_file);
if save_beta_barplot
    fprintf('Saved labeled effective-beta bar plot to:\n%s\n', ...
        beta_barplot_file);
end

%% Local utilities

function probability = tanh_probability(parameters,Isyn_mA)
% parameters = [I0_mA, b]; I0 is the fitted physical-current scale.
if numel(parameters)~=2 || any(~isfinite(parameters)) || parameters(1)==0
    probability = nan(size(Isyn_mA));
    return
end
probability = 0.5*(1+tanh(Isyn_mA./parameters(1)+parameters(2)));
end

function nll = binomial_negative_log_likelihood( ...
        parameters,Isyn_mA,counts,nTrials)
probability = tanh_probability(parameters,Isyn_mA);
if any(~isfinite(probability))
    nll = realmax('double');
    return
end
probability = max(min(probability,1-1e-12),1e-12);
nll = -sum(counts.*log(probability) ...
    +(nTrials-counts).*log(1-probability));
if ~isfinite(nll)
    nll = realmax('double');
end
end

function sse = probability_sum_squares(parameters,Isyn_mA,observed)
probability = tanh_probability(parameters,Isyn_mA);
if any(~isfinite(probability))
    sse = realmax('double');
    return
end
sse = sum((probability-observed).^2);
end

function density = gaussian_kernel_density(samples,evaluation_points,bandwidth)
% Normalized Gaussian KDE; works without ksdensity or Statistics Toolbox.
standardized_distance = ...
    (evaluation_points(:)-samples(:).')/bandwidth;
density = mean(exp(-0.5*standardized_distance.^2),2) ...
    /(sqrt(2*pi)*bandwidth);
density = reshape(density,size(evaluation_points));
end

function standardized = standardize_finite(values)
standardized = nan(size(values));
valid = isfinite(values);
if nnz(valid)<2
    return
end
scale = std(values(valid),0);
if scale<=0 || ~isfinite(scale)
    standardized(valid) = 0;
else
    standardized(valid) = (values(valid)-mean(values(valid)))/scale;
end
end

function validate_phase_audit(audit,probability,gamma_values,nTrials)
required = {'raw_relative_phase_deg','calibrated_relative_phase_deg', ...
    'fundamental_amplitude_V','low_amplitude','minimum_amplitude_V', ...
    'two_state_calibration','lock_calibration','neutral_calibration', ...
    'preparation_offset_deg','lock_target_phase_deg', ...
    'neutral_offset_error_deg', ...
    'locked_mean_phase_deg','locked_mean_phase_error_deg', ...
    'locked_phase_concentration', ...
    'synaptic_drive_phase_deg','lock_fhi_off_before_sampling', ...
    'lock_fhi_off_time_s','sampling_start_time_s','shi_full_time_s', ...
    'zero_input_bias','mirror_symmetry_mae', ...
    'readout_fhi_removed','binary_assignment', ...
    'fixed_axis_deg','observed_axis_deg','doubled_phase_concentration', ...
    'near_boundary','gamma_values','assignment_convention', ...
    'slope_convention','endpoint_probability_difference', ...
    'boundary_warning_width_deg'};
if ~isstruct(audit) || ~all(isfield(audit,required))
    error(['Missing calibrated parametric-SHI audit fields for a ' ...
        'completed curve; rerun the updated ensemble.']);
end
raw_phase = squeeze(audit.raw_relative_phase_deg);
calibrated_phase = squeeze(audit.calibrated_relative_phase_deg);
amplitude = squeeze(audit.fundamental_amplitude_V);
low_amplitude = squeeze(audit.low_amplitude);
bits = squeeze(audit.binary_assignment);
expected_size = [numel(gamma_values) nTrials];
calibration = audit.two_state_calibration;
lock_calibration = audit.lock_calibration;
neutral_calibration = audit.neutral_calibration;
if ~strcmp(audit.assignment_convention, ...
        'independently_calibrated_antipodal_bit1_zero_fhi90_v2') || ...
        ~strcmp(audit.slope_convention, ...
        'signed_synaptic_FHI_quadrature_v1') || ...
        ~isscalar(audit.fixed_axis_deg) || ...
        ~isfinite(audit.fixed_axis_deg) || ...
        ~isscalar(audit.readout_fhi_removed) || ...
        ~logical(audit.readout_fhi_removed) || ...
        ~isstruct(lock_calibration) || ...
        ~isfield(lock_calibration,'passed') || ...
        ~logical(lock_calibration.passed) || ...
        ~isstruct(neutral_calibration) || ...
        ~isfield(neutral_calibration,'offset_deg') || ...
        ~isfinite(audit.preparation_offset_deg) || ...
        ~logical(audit.lock_fhi_off_before_sampling) || ...
        audit.lock_fhi_off_time_s>audit.sampling_start_time_s || ...
        audit.sampling_start_time_s>=audit.shi_full_time_s || ...
        ~isfinite(audit.lock_target_phase_deg) || ...
        ~isfinite(audit.synaptic_drive_phase_deg) || ...
        ~isfinite(audit.neutral_offset_error_deg) || ...
        ~isfinite(audit.locked_mean_phase_deg) || ...
        ~isfinite(audit.locked_mean_phase_error_deg) || ...
        ~isfinite(audit.locked_phase_concentration) || ...
        ~isequal(audit.gamma_values(:),gamma_values(:)) || ...
        ~isequal(size(raw_phase),expected_size) || ...
        ~isequal(size(calibrated_phase),expected_size) || ...
        ~isequal(size(amplitude),expected_size) || ...
        ~isequal(size(low_amplitude),expected_size) || ...
        ~isequal(size(bits),expected_size) || ...
        any(~isfinite(raw_phase(:))) || ...
        any(raw_phase(:)<-180 | raw_phase(:)>=180) || ...
        any(~isfinite(calibrated_phase(:))) || ...
        any(calibrated_phase(:)<-180 | calibrated_phase(:)>=180) || ...
        any(~isfinite(amplitude(:)) | amplitude(:)<0) || ...
        any(~ismember(bits(:),[0 1]))
    error('Invalid calibrated parametric-SHI audit metadata or dimensions.');
end
if ~isstruct(calibration) || ~isfield(calibration,'passed') || ...
        ~logical(calibration.passed) || ...
        ~isfield(calibration,'bit1_axis_deg') || ...
        abs(wrap_to_180(audit.fixed_axis_deg- ...
        calibration.bit1_axis_deg))>1e-10
    error('The independent deterministic two-state calibration is invalid.');
end
expected_calibrated = wrap_to_180(raw_phase-audit.fixed_axis_deg);
if any(abs(wrap_to_180(calibrated_phase(:) ...
        -expected_calibrated(:)))>1e-9)
    error('Saved calibrated phases do not match the saved calibration axis.');
end
expected_bits = abs(calibrated_phase)<=90;
if ~isequal(logical(bits),expected_bits)
    error('Saved bits differ from the calibrated antipodal phase rule.');
end
if ~isequal(logical(low_amplitude), ...
        amplitude<audit.minimum_amplitude_V)
    error('Saved low-amplitude flags are inconsistent with the threshold.');
end
expected_near_boundary = abs(abs(calibrated_phase)-90) ...
    <=audit.boundary_warning_width_deg;
if ~isequal(logical(audit.near_boundary),expected_near_boundary)
    error('Saved phase-boundary flags are inconsistent.');
end
recomputed_probability = mean(double(bits),2);
if any(~isfinite(probability(:))) || ...
        any(abs(recomputed_probability-probability(:))>1e-12)
    error('Saved P(bit=1) does not equal the mean of the saved binary outcomes.');
end

function angle_deg = wrap_to_180(angle_deg)
angle_deg = mod(angle_deg+180,360)-180;
end
end

function draw_phase_guides()
yline(0,'k-','HandleVisibility','off');
yline(180,'k-','HandleVisibility','off');
yline(-180,'k-','HandleVisibility','off');
yline(90,'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
yline(-90,'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
ylim([-190 190]); yticks(-180:90:180);
end

