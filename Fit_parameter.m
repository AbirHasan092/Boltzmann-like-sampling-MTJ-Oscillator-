%% fit_probability_tanh_gamma.m
% Fit MTJ probabilities using the signed Oersted-ring charge current Ic:
%
%    P(+z) = 0.5*(1+tanh(Ic_mA/I0_mA+b)).
%
% I0_mA is fitted directly; beta_eff_mA_inv=1/I0_mA has units mA^(-1).
% The default bias is fixed at +2, as requested. Set fit_bias_b=true to fit b
% instead, as in the LC oscillator fitter. A fixed +2 enforces P(Ic=0)=0.982.
% The spin current Is is NOT the fitted synaptic input.
% Plot a colored beta bar chart and a Gaussian-kernel PDF (KDE), not a CDF.
% Existing parameter, spin-current, and initialization plots are retained.
% Uses fminsearch; KDE does not require the Statistics Toolbox.

clear; clc; close all;

%% User settings
script_dir = fileparts(mfilename('fullpath'));
results_file = fullfile(script_dir,'INI_zonly_parameter_variation_results_0.18A.mat');
smooth_point_count = 1000;

% Nominal Oersted-ring current used to normalize the synaptic input:
%
%                Ic = gamma * nominal_Ic_mA.
nominal_Ic_mA = 10;               % [mA]; +/-10 mA maps to gamma = +/-1
% Fit only I0 for the explicitly requested +2 offset. Set true to fit b.
fit_bias_b = true;
fixed_bias_b = 2;
beta_pdf_bandwidth = [0.03]; % empty: automatic Silverman bandwidth, mA^(-1)
fit_file = fullfile(script_dir,'tanh_fit_results_I0_KDE.mat');
fit_table_file = fullfile(script_dir,'tanh_fit_table_I0_KDE.csv');
output_suffix = '_I0_current_KDE';
save_beta_pdf_plot = true;

% USER: set this to the folder where the parameter bar plots are saved.
% An absolute or relative folder path may be used.
barplot_output_folder = fullfile(script_dir,'INI2_parameter_zscore_barplots');

% Used only when variation_flags.only_vary_Is_direction is true. Each requested
% gamma value produces one probability curve in both the theta and phi
% panels. Values are dimensionless and must lie inside the simulated range.
% These targets are retained from the input script; gamma=Ic/nominal_Ic_mA.
representative_synaptic_gamma = [-0.5875 0 0.5875];
number_of_direction_probability_curves = 3;

% Number of subplot columns in the consolidated all-sweeps bar-plot figure.
concatenated_barplot_columns = 5;

%% Load probability points
if ~isfile(results_file)
    error(['Missing %s. Run ' ...
        'SOT_Elliptic_Ensemble_simulation_and_sensitivity_analysis.m ' ...
        'before running this fitting script.'],results_file);
end

S = load(results_file,'Ic_values','Pplus_all','parameter_sets', ...
    'variation_flags','trace','time','representative_sweep', ...
    'Ic_trace_actual');
if ~isfield(S,'Ic_values') || ~isfield(S,'Pplus_all') || ...
        ~isfield(S,'parameter_sets') || ~isfield(S,'variation_flags') || ...
        ~isfield(S,'trace') || ~isfield(S,'time') || ...
        ~isfield(S,'representative_sweep') || ...
        ~isfield(S,'Ic_trace_actual')
    error(['%s must contain Ic_values, Pplus_all, parameter_sets, and ' ...
        'variation_flags, plus the representative trace and time data. ' ...
        'Run the updated generator and ensemble scripts again.'], ...
        results_file);
end

Ic_mA = S.Ic_values(:).'*1e3;
if ~isscalar(nominal_Ic_mA) || ~isfinite(nominal_Ic_mA) || ...
        nominal_Ic_mA<=0
    error('nominal_Ic_mA must be a positive finite scalar.');
end

% Ic_mA is the physical Oersted current used in fitting and probability plots.
% gamma remains available for initialization diagnostics and compatibility.
gamma_values = Ic_mA/nominal_Ic_mA;
validateattributes(fit_bias_b,{'logical'},{'scalar'});
validateattributes(fixed_bias_b,{'numeric'},{'scalar','real','finite'});
current_unit = 'mA';
beta_eff_unit = 'mA^(-1)';
if fit_bias_b
    fit_parameterization = '0.5*(1+tanh(Ic_mA/I0_mA+b))';
else
    fit_parameterization = sprintf( ...
        '0.5*(1+tanh(Ic_mA/I0_mA%+.8g))',fixed_bias_b);
    fprintf(['Fixed bias b=%g: P(Ic=0)=%.8g for every fitted curve.\n' ...
        'This is imposed by the model, not inferred from the data.\n'], ...
        fixed_bias_b,0.5*(1+tanh(fixed_bias_b)));
end
Pdata = S.Pplus_all;

if isvector(Pdata)
    Pdata = reshape(Pdata,1,[]);
end
if size(Pdata,2)~=numel(Ic_mA) && size(Pdata,1)==numel(Ic_mA)
    Pdata = Pdata.';
end
if size(Pdata,2)~=numel(Ic_mA)
    error('Pplus_all dimensions do not match the number of Ic_values.');
end
if any(~isfinite(Pdata),'all') || any(Pdata<0,'all') || any(Pdata>1,'all')
    error('Pplus_all must contain finite probabilities between 0 and 1.');
end

Nsweep = size(Pdata,1);
if numel(S.parameter_sets)<Nsweep
    error('parameter_sets contains fewer entries than Pplus_all.');
end
parameter_sets_used = S.parameter_sets(1:Nsweep);
variation_flags = S.variation_flags;
representative_trace = S.trace;
representative_time = S.time(:).';
representative_sweep = S.representative_sweep;
representative_Ic_trace_actual = S.Ic_trace_actual;

if ~isscalar(representative_sweep) || ~isfinite(representative_sweep) || ...
        representative_sweep~=floor(representative_sweep) || ...
        representative_sweep<1 || representative_sweep>Nsweep
    error('representative_sweep must identify one of the fitted sweeps.');
end
required_trace_fields = {'m','Is','stored'};
for trace_field_index = 1:numel(required_trace_fields)
    if ~isfield(representative_trace, ...
            required_trace_fields{trace_field_index})
        error('Representative trace is missing %s.', ...
            required_trace_fields{trace_field_index});
    end
end
if ~representative_trace.stored || ...
        size(representative_trace.m,1)~=3 || ...
        size(representative_trace.Is,1)~=3 || ...
        size(representative_trace.m,2)~=numel(representative_time) || ...
        size(representative_trace.Is,2)~= ...
        numel(representative_time) || ...
        any(~isfinite(representative_time)) || ...
        any(~isfinite(representative_trace.m),'all') || ...
        any(~isfinite(representative_trace.Is),'all')
    error(['The representative trace or its time vector is invalid. ' ...
        'Rerun the updated ensemble script.']);
end
if ~isscalar(representative_Ic_trace_actual) || ...
        ~isfinite(representative_Ic_trace_actual)
    error('Ic_trace_actual must be a finite scalar.');
end
representative_gamma_trace_actual = ...
    representative_Ic_trace_actual*1e3/nominal_Ic_mA;
required_flag_names = {'vary_Is_direction_with_device_variation', ...
    'vary_Ic_direction_with_device_variation','only_vary_Is_direction'};
for flag_index = 1:numel(required_flag_names)
    if ~isfield(variation_flags,required_flag_names{flag_index})
        error('variation_flags is missing %s.', ...
            required_flag_names{flag_index});
    end
end
gamma_fit = linspace(min(gamma_values),max(gamma_values), ...
    smooth_point_count);


Ic_fit_mA = nominal_Ic_mA*gamma_fit;
I0_mA = nan(Nsweep,1);
beta_eff_mA_inv = nan(Nsweep,1);
bias_b = zeros(Nsweep,1);
RMSE = zeros(Nsweep,1);
R_squared = zeros(Nsweep,1);
fit_exitflag = zeros(Nsweep,1);

%% Fit each parameter variation directly in Oersted-current units
fit_options = optimset('Display','off','MaxIter',10000, ...
    'MaxFunEvals',20000,'TolX',1e-12,'TolFun',1e-12);
current_scale_mA = max(max(abs(Ic_mA)),eps);

for s = 1:Nsweep
    P = Pdata(s,:);
    transformed_probability = 2*P-1;
    clipped = max(min(transformed_probability,0.999),-0.999);
    linearized = atanh(clipped);
    nonsaturated = P>0.02 & P<0.98;
    if nnz(nonsaturated)>=2
        initial_line = polyfit(Ic_mA(nonsaturated),linearized(nonsaturated),1);
    else
        initial_line = polyfit(Ic_mA,linearized,1);
    end
    if any(~isfinite(initial_line))
        initial_line = [1/current_scale_mA 0];
    end
    if abs(initial_line(1))>1e-12/current_scale_mA
        initial_I0_mA = 1/initial_line(1);
    else
        initial_I0_mA = current_scale_mA;
    end
    if fit_bias_b
        initial_parameters = [initial_I0_mA initial_line(2)];
    else
        initial_parameters = initial_I0_mA;
    end
    objective = @(parameters) 4*current_tanh_sse( ...
        parameters,Ic_mA,P,fit_bias_b,fixed_bias_b);
    [fitted,~,fit_exitflag(s)] = ...
        fminsearch(objective,initial_parameters,fit_options);
    if any(~isfinite(fitted)) || fitted(1)==0
        error('Variation %d returned an invalid I0 fit.',s);
    end
    if fit_exitflag(s)<=0
        warning('CurrentFit:NotConverged', ...
            'I0 fit for variation %d did not converge.',s);
    end
    I0_mA(s) = fitted(1);
    beta_eff_mA_inv(s) = 1/I0_mA(s);
    if fit_bias_b
        bias_b(s) = fitted(2);
    else
        bias_b(s) = fixed_bias_b;
    end

    Ppredicted = current_tanh_probability([I0_mA(s) bias_b(s)],Ic_mA);
    residual = P-Ppredicted;
    RMSE(s) = sqrt(mean(residual.^2));
    total_variation = sum((P-mean(P)).^2);
    if total_variation>eps
        R_squared(s) = 1-sum(residual.^2)/total_variation;
    else
        R_squared(s) = NaN;
    end
end

% Preserve signed slopes and fixed binary labels, including decreasing fits.
I0_A = 1e-3*I0_mA;
beta_eff_A_inv = 1./I0_A;
% Explicit derived compatibility quantities; beta_gamma is dimensionless.
beta_per_mA = beta_eff_mA_inv;
beta_per_A = beta_eff_A_inv;
beta_gamma = nominal_Ic_mA*beta_eff_mA_inv;

%% Fit one effective current-scale tanh to all variation-level tanh curves
% Equal weights for all variations; retain the original second-stage method.
fitted_sweep_curves = zeros(Nsweep,smooth_point_count);
for s = 1:Nsweep
    fitted_sweep_curves(s,:) = current_tanh_probability( ...
        [I0_mA(s) bias_b(s)],Ic_fit_mA);
end
mean_fitted_curve = mean(fitted_sweep_curves,1);
if Nsweep>1
    standard_error_curve = std(fitted_sweep_curves,0,1)/sqrt(Nsweep);
else
    standard_error_curve = zeros(size(mean_fitted_curve));
end
standard_error_lower = max(0,mean_fitted_curve-standard_error_curve);
standard_error_upper = min(1,mean_fitted_curve+standard_error_curve);

initial_mean_beta = mean(beta_eff_mA_inv);
if abs(initial_mean_beta)>1e-12/current_scale_mA
    effective_initial_I0_mA = 1/initial_mean_beta;
else
    effective_initial_I0_mA = current_scale_mA;
end
if fit_bias_b
    effective_initial = [effective_initial_I0_mA mean(bias_b)];
else
    effective_initial = effective_initial_I0_mA;
end
effective_objective = @(parameters) current_tanh_sse( ...
    parameters,Ic_fit_mA,mean_fitted_curve,fit_bias_b,fixed_bias_b);
[effective_parameters,~,effective_fit_exitflag] = fminsearch( ...
    effective_objective,effective_initial,fit_options);
if any(~isfinite(effective_parameters)) || effective_parameters(1)==0
    error('The ensemble-effective I0 fit is invalid.');
end
if effective_fit_exitflag<=0
    warning('CurrentFit:EffectiveNotConverged', ...
        'The ensemble-effective I0 fit did not converge.');
end
effective_I0_mA = effective_parameters(1);
effective_I0_A = 1e-3*effective_I0_mA;
effective_beta_eff_mA_inv = 1/effective_I0_mA;
effective_beta_eff_A_inv = 1/effective_I0_A;
effective_beta_per_mA = effective_beta_eff_mA_inv;
effective_beta_per_A = effective_beta_eff_A_inv;
effective_beta_gamma = nominal_Ic_mA*effective_beta_eff_mA_inv;
if fit_bias_b
    effective_bias_b = effective_parameters(2);
else
    effective_bias_b = fixed_bias_b;
end
effective_probability_curve = current_tanh_probability( ...
    [effective_I0_mA effective_bias_b],Ic_fit_mA);

% Z-score normalization for comparing the first-stage fitted beta and
% bias distributions on the same dimensionless scale.
beta_standard_deviation = std(beta_eff_mA_inv,0);
if isfinite(beta_standard_deviation) && beta_standard_deviation>0
    beta_z_score = (beta_eff_mA_inv-mean(beta_eff_mA_inv))/ ...
        beta_standard_deviation;
else
    beta_z_score = zeros(size(beta_eff_mA_inv));
end

bias_standard_deviation = std(bias_b,0);
if isfinite(bias_standard_deviation) && bias_standard_deviation>0
    bias_z_score = (bias_b-mean(bias_b))/bias_standard_deviation;
else
    bias_z_score = zeros(size(bias_b));
end

%% Index of dispersion of the sweep-level fitted parameters
% The standard index of dispersion is the sample variance divided by the
% sample mean. It is evaluated across the first-stage sweep fits because a
% single second-stage effective fit does not itself define a distribution.
beta_dispersion_mean = mean(beta_eff_mA_inv);
beta_dispersion_variance = var(beta_eff_mA_inv,0);
bias_dispersion_mean = mean(bias_b);
bias_dispersion_variance = var(bias_b,0);

beta_mean_tolerance = 1e-12*max(1,max(abs(beta_eff_mA_inv)));
bias_mean_tolerance = 1e-12*max(1,max(abs(bias_b)));

if abs(beta_dispersion_mean)>beta_mean_tolerance
    beta_index_of_dispersion = beta_dispersion_variance/ ...
        beta_dispersion_mean;
else
    beta_index_of_dispersion = NaN;
end

if abs(bias_dispersion_mean)>bias_mean_tolerance
    bias_index_of_dispersion = bias_dispersion_variance/ ...
        bias_dispersion_mean;
else
    bias_index_of_dispersion = NaN;
end

beta_sample_SD = std(beta_eff_mA_inv,0);
if abs(beta_dispersion_mean)>beta_mean_tolerance
    beta_CV_percent = 100*beta_sample_SD/abs(beta_dispersion_mean);
else
    beta_CV_percent = NaN;
end
I0_mean_mA = mean(I0_mA);
I0_sample_SD_mA = std(I0_mA,0);

%% Normalized input-offset spread across sweep fits
% At P(+z)=0.5, Ic_50=-b*I0 and gamma_50=Ic_50/nominal_Ic_mA.
% With fixed b, midpoint variation follows I0 variation by construction.
% Its sample standard deviation is the current-offset standard deviation
% divided by nominal_Ic_mA. This measures horizontal offsets rather than
% the spread of the tanh intercept b. Do not divide by the mean offset,
% which may be zero for a symmetric distribution.
input_offset_gamma = nan(Nsweep,1);
offset_slope_tolerance = 1e-12;
valid_offset_fit = isfinite(I0_mA) & isfinite(beta_eff_mA_inv) & ...
    isfinite(bias_b) & abs(beta_eff_mA_inv)>offset_slope_tolerance;
input_offset_gamma(valid_offset_fit) = ...
    -bias_b(valid_offset_fit).*I0_mA(valid_offset_fit)/nominal_Ic_mA;
valid_offset_fit = valid_offset_fit & isfinite(input_offset_gamma);
input_offset_gamma(~valid_offset_fit) = NaN;
input_offset_mA = nominal_Ic_mA*input_offset_gamma;
normalized_offset_mean = NaN;
normalized_offset_spread = NaN;
% Do not silently discard invalid sweeps from a full-ensemble statistic.
if all(valid_offset_fit) && Nsweep>=2
    normalized_offset_mean = mean(input_offset_gamma);
    normalized_offset_spread = std(input_offset_gamma,0);
else
    warning('OffsetSpread:Undefined', ...
        ['Normalized offset spread is undefined: %d of %d sweeps have ' ...
        'finite midpoints and nonzero slopes. At least two sweeps, ' ...
        'with valid midpoints for every sweep, are required.'], ...
        nnz(valid_offset_fit),Nsweep);
end
offset_spread_mA = nominal_Ic_mA*normalized_offset_spread;
normalized_offset_spread_percent = 100*normalized_offset_spread;

%% Independently varied device and direction parameters for each sweep
% Derived quantities (V, Ns, Nx, Ny, and Nz) are intentionally excluded.
% Their values follow from Ms and the sampled ellipse dimensions and are not
% independent process variations in the Debashis parameter generator.
mtj_parameter_names = {'M_s','\alpha','T','D_x','D_y','t_F','H_k', ...
    'R_{ring}','\Delta I_{s,x}','\Delta I_{s,y}','\Delta I_{s,z}', ...
    '\Delta H_{Oe,x}','\Delta H_{Oe,y}','\Delta H_{Oe,z}'};
mtj_parameter_values = zeros(Nsweep,numel(mtj_parameter_names));

for s = 1:Nsweep
    p = parameter_sets_used(s);
    Is_offset = reshape(p.Is_offset,1,3);
    Oe_offset = reshape(p.Oe_offset,1,3);
    mtj_parameter_values(s,:) = [p.Ms p.alpha p.T ...
        p.diameter_x p.diameter_y p.tF p.Hk_Oe p.R_ring ...
        Is_offset Oe_offset];
end

mtj_parameter_means = mean(mtj_parameter_values,1);
mtj_parameter_standard_deviations = std(mtj_parameter_values,0,1);
mtj_parameter_z_scores = zeros(size(mtj_parameter_values));
for j = 1:size(mtj_parameter_values,2)
    sigma_j = mtj_parameter_standard_deviations(j);
    if isfinite(sigma_j) && sigma_j>0
        mtj_parameter_z_scores(:,j) = ...
            (mtj_parameter_values(:,j)-mtj_parameter_means(j))/sigma_j;
    end
end

%% Is-direction angles and selected-current probability data
% These arrays remain empty unless the Is-direction-only flag is enabled. The
% angles are calculated from the normalized direction actually used by
% the ensemble: Is_hat = (Is_ideal + offset)/|Is_ideal + offset|.
Is_initial_theta_deg = [];
Is_initial_phi_deg = [];
selected_synaptic_gamma_requested = [];
selected_synaptic_gamma_actual = [];
selected_synaptic_gamma_indices = [];
selected_direction_probability = [];

if variation_flags.only_vary_Is_direction
    if ~isscalar(number_of_direction_probability_curves) || ...
            ~isfinite(number_of_direction_probability_curves) || ...
            number_of_direction_probability_curves<1 || ...
            number_of_direction_probability_curves~= ...
            floor(number_of_direction_probability_curves)
        error(['number_of_direction_probability_curves must be a ' ...
            'positive integer.']);
    end
    if number_of_direction_probability_curves> ...
            numel(representative_synaptic_gamma)
        error(['The requested curve count exceeds the number of values ' ...
            'in representative_synaptic_gamma.']);
    end

    selected_synaptic_gamma_requested = ...
        representative_synaptic_gamma( ...
        1:number_of_direction_probability_curves);
    if any(~isfinite(selected_synaptic_gamma_requested)) || ...
            any(selected_synaptic_gamma_requested<min(gamma_values)) || ...
            any(selected_synaptic_gamma_requested>max(gamma_values))
        error(['Every representative synaptic gamma value must be ' ...
            'finite and inside [%.6g, %.6g].'], ...
            min(gamma_values),max(gamma_values));
    end

    selected_synaptic_gamma_indices = zeros(1, ...
        number_of_direction_probability_curves);
    selected_synaptic_gamma_actual = zeros(1, ...
        number_of_direction_probability_curves);
    for current_index = 1:number_of_direction_probability_curves
        [~,selected_synaptic_gamma_indices(current_index)] = min(abs( ...
            gamma_values-selected_synaptic_gamma_requested(current_index)));
        selected_synaptic_gamma_actual(current_index) = ...
            gamma_values(selected_synaptic_gamma_indices(current_index));
    end
    if numel(unique(selected_synaptic_gamma_indices))< ...
            number_of_direction_probability_curves
        error(['Two or more requested gamma values map to the same ' ...
            'simulated point. Choose values farther apart.']);
    end

    Is_initial_theta_deg = zeros(Nsweep,1);
    Is_initial_phi_deg = zeros(Nsweep,1);
    for s = 1:Nsweep
        varied_Is_direction = ...
            parameter_sets_used(s).Is_direction_ideal(:)+ ...
            parameter_sets_used(s).Is_offset(:);
        varied_direction_norm = norm(varied_Is_direction);
        if ~isfinite(varied_direction_norm) || varied_direction_norm<1e-12
            error('Sweep %d has an invalid Is direction.',s);
        end
        varied_Is_direction = varied_Is_direction/varied_direction_norm;
        Is_initial_theta_deg(s) = acosd(max(min( ...
            varied_Is_direction(3),1),-1));
        Is_initial_phi_deg(s) = atan2d( ...
            varied_Is_direction(2),varied_Is_direction(1));
    end
    selected_direction_probability = ...
        Pdata(:,selected_synaptic_gamma_indices);

    for current_index = 1:number_of_direction_probability_curves
        fprintf(['Direction plot gamma %d: requested %+.6g, ' ...
            'using simulated %+.6g (Ic=%+.6g mA).\n'],current_index, ...
            selected_synaptic_gamma_requested(current_index), ...
            selected_synaptic_gamma_actual(current_index), ...
            selected_synaptic_gamma_actual(current_index)*nominal_Ic_mA);
    end
end

effective_residual = mean_fitted_curve-effective_probability_curve;
effective_RMSE = sqrt(mean(effective_residual.^2));
effective_total_variation = sum((mean_fitted_curve- ...
    mean(mean_fitted_curve)).^2);
if effective_total_variation>eps
    effective_R_squared = 1-sum(effective_residual.^2)/ ...
        effective_total_variation;
else
    effective_R_squared = NaN;
end


%% Report and save fitted parameters
fprintf('\nPhysical fit input: Oersted-ring charge current Ic in mA.\n');
fprintf('Model: %s; beta_eff=1/I0 has units mA^(-1).\n',fit_parameterization);
fprintf('Gamma remains Ic/(%.6g mA) for initialization diagnostics.\n',nominal_Ic_mA);
Variation = (1:Nsweep).';
Sweep = Variation; % preserved index name for the existing plotting code
fit_results = table(Variation,I0_mA,beta_eff_mA_inv,I0_A,beta_eff_A_inv, ...
    beta_gamma,bias_b,input_offset_gamma,input_offset_mA,RMSE,R_squared, ...
    fit_exitflag);
disp(fit_results);

effective_fit_results = table(effective_I0_mA,effective_beta_eff_mA_inv, ...
    effective_I0_A,effective_beta_eff_A_inv,effective_beta_gamma, ...
    effective_bias_b,effective_RMSE,effective_R_squared,effective_fit_exitflag);
disp('Effective current-scale tanh fit to all variation-level tanh curves:');
disp(effective_fit_results);

fprintf('\nMean I0 = %.8g mA; sample SD = %.8g mA.\n', ...
    I0_mean_mA,I0_sample_SD_mA);
fprintf('Mean beta_eff = %.8g mA^(-1); sample SD = %.8g mA^(-1).\n', ...
    beta_dispersion_mean,beta_sample_SD);
fprintf('Beta coefficient of variation = %.8g percent.\n',beta_CV_percent);
fprintf('D_beta = variance/mean = %.8g mA^(-1) (unit dependent).\n', ...
    beta_index_of_dispersion);
if fit_bias_b
    fprintf('Fitted bias b: mean = %.8g; sample SD = %.8g.\n', ...
        mean(bias_b),std(bias_b,0));
else
    fprintf(['Bias b is FIXED at %g, not fitted. Its zero spread is a ' ...
        'model constraint, not evidence of device uniformity.\n'],fixed_bias_b);
end

fprintf('\nMidpoint Ic_50=-b*I0; gamma_50=Ic_50/Ic_norm.\n');
fprintf('Valid midpoint fits: %d/%d.\n',nnz(valid_offset_fit),Nsweep);
fprintf('SD(gamma_50) = %.8g (%.8g percent of nominal current).\n', ...
    normalized_offset_spread,normalized_offset_spread_percent);
fprintf('Midpoint-current sample SD = %.8g mA.\n',offset_spread_mA);
writetable(fit_results,fit_table_file);

save(fit_file,'fit_results', ...
    'I0_mA','I0_A','beta_eff_mA_inv','beta_eff_A_inv', ...
    'beta_gamma','beta_per_mA','beta_per_A','bias_b','fit_exitflag', ...
    'effective_I0_mA','effective_I0_A','effective_beta_eff_mA_inv', ...
    'effective_beta_eff_A_inv','effective_beta_gamma','effective_bias_b', ...
    'effective_fit_exitflag','I0_mean_mA','I0_sample_SD_mA', ...
    'beta_sample_SD','beta_CV_percent','Ic_fit_mA', ...
    'fit_bias_b','fixed_bias_b','fit_parameterization', ...
    'current_unit','beta_eff_unit', ...
    'input_offset_gamma','input_offset_mA','valid_offset_fit', ...
    'offset_slope_tolerance','normalized_offset_mean', ...
    'normalized_offset_spread','normalized_offset_spread_percent', ...
    'offset_spread_mA', ...
    'effective_fit_results','nominal_Ic_mA','gamma_values','gamma_fit', ...
    'Ic_mA','Pdata', ...
    'fitted_sweep_curves','mean_fitted_curve', ...
    'standard_error_curve','standard_error_lower', ...
    'standard_error_upper','effective_probability_curve', ...
    'beta_z_score','bias_z_score','beta_index_of_dispersion', ...
    'bias_index_of_dispersion','beta_dispersion_mean', ...
    'beta_dispersion_variance','bias_dispersion_mean', ...
    'bias_dispersion_variance','mtj_parameter_names', ...
    'mtj_parameter_values','mtj_parameter_means', ...
    'mtj_parameter_standard_deviations','mtj_parameter_z_scores', ...
    'variation_flags','Is_initial_theta_deg','Is_initial_phi_deg', ...
    'representative_synaptic_gamma', ...
    'number_of_direction_probability_curves', ...
    'selected_synaptic_gamma_requested', ...
    'selected_synaptic_gamma_actual','selected_synaptic_gamma_indices', ...
    'selected_direction_probability','representative_trace', ...
    'representative_time','representative_sweep', ...
    'representative_Ic_trace_actual','representative_gamma_trace_actual');


%% Plot probability points and fitted curves versus Oersted current
figFit = figure('Color','w','Name','Current-scale MTJ tanh fits', ...
    'NumberTitle','off');
clf(figFit); hold on;
colors = turbo(Nsweep);
fit_handles = gobjects(Nsweep,1);
for s = 1:Nsweep
    c = colors(s,:);
    scatter(Ic_mA,Pdata(s,:),48,c,'filled', ...
        'MarkerEdgeColor',c,'HandleVisibility','off');
    Pfit = current_tanh_probability([I0_mA(s) bias_b(s)],Ic_fit_mA);
    fit_handles(s) = plot(Ic_fit_mA,Pfit,'-','Color',c,'LineWidth',2.4, ...
        'DisplayName',sprintf( ...
        'Variation %d: I_0=%.4g mA, \\beta_{eff}=%.4g mA^{-1}, b=%+.4g', ...
        s,I0_mA(s),beta_eff_mA_inv(s),bias_b(s)));
end
xlabel('Oersted-ring current I_c (mA)','FontSize',18,'Interpreter','tex');
ylabel('P(+z)','FontSize',18);
if fit_bias_b
    title('$P(+z)=\frac{1}{2}[1+\tanh(I_c/I_0+b)]$', ...
        'Interpreter','latex');
else
    title(sprintf('$P(+z)=\\frac{1}{2}[1+\\tanh(I_c/I_0%+.8g)]$', ...
        fixed_bias_b),'Interpreter','latex');
end
legend(fit_handles,'Location','best','FontSize',11,'Interpreter','tex');
set(gca,'FontSize',16,'LineWidth',1.8,'TickDir','in');
xlim([min(Ic_mA) max(Ic_mA)]);
xticks(linspace(min(Ic_mA),max(Ic_mA),5));
ylim([-0.05 1.05]); box on; grid off;

%% Plot the effective tanh fitted from all variation-level tanh curves
figEffective = figure('Color','w','Name','Effective MTJ current response', ...
    'NumberTitle','off');
clf(figEffective); hold on;
individual_handle = plot(Ic_fit_mA,fitted_sweep_curves(1,:),'--', ...
    'Color',[0.72 0.72 0.72],'LineWidth',1.2, ...
    'DisplayName','Variation fits');
for s = 2:Nsweep
    plot(Ic_fit_mA,fitted_sweep_curves(s,:),'--', ...
        'Color',[0.72 0.72 0.72],'LineWidth',1.2,'HandleVisibility','off');
end
effective_handle = plot(Ic_fit_mA,effective_probability_curve,'r-', ...
    'LineWidth',3.0,'DisplayName',sprintf( ...
    'Effective: I_0=%.4g mA, \\beta_{eff}=%.4g mA^{-1}, b=%+.4g', ...
    effective_I0_mA,effective_beta_eff_mA_inv,effective_bias_b));
xlabel('Oersted-ring current I_c (mA)','FontSize',18,'Interpreter','tex');
ylabel('P(+z)','FontSize',18);
legend([individual_handle effective_handle], ...
    'Location','best','FontSize',11,'Interpreter','tex');
set(gca,'FontSize',16,'LineWidth',1.8,'TickDir','in');
xlim([min(Ic_mA) max(Ic_mA)]);
xticks(linspace(min(Ic_mA),max(Ic_mA),5));
ylim([-0.05 1.05]); box on; grid off;

%% Plot the pointwise standard error in a separate figure
figSEM = figure('Color','w','Name','SEM of variation tanh curves', ...
    'NumberTitle','off');
clf(figSEM);
plot(Ic_fit_mA,standard_error_curve,'b-','LineWidth',2.5);
xlabel('Oersted-ring current I_c (mA)','FontSize',18,'Interpreter','tex');
ylabel('SEM[P(+z)]','FontSize',18);
set(gca,'FontSize',16,'LineWidth',1.8,'TickDir','in');
xlim([min(Ic_mA) max(Ic_mA)]);
xticks(linspace(min(Ic_mA),max(Ic_mA),5));
ylim([0 max(standard_error_curve)*1.10+eps]);
box on; grid off;

%% Box plots of the first-stage sweep fit parameters
figFitParameters = figure('Color','w','Name', ...
    'Distribution of variation tanh parameters','NumberTitle','off');
clf(figFitParameters);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
boxplot(beta_z_score,'Labels',{''},'Whisker',1.5, ...
    'Symbol','r+');
ylabel('Z-score');
title('Standardized variation \beta_{eff}', ...
    'Interpreter','tex','FontSize',12);
set(gca,'XTick',1,'XTickLabel',{'\beta_{eff}'}, ...
    'TickLabelInterpreter','tex','FontSize',15, ...
    'LineWidth',1.5,'TickDir','in');
box on; grid off;

nexttile;
boxplot(bias_z_score,'Labels',{''},'Whisker',1.5,'Symbol','r+');
ylabel('Z-score');
if fit_bias_b
    title('Standardized variation bias b','FontSize',12);
else
    title(sprintf('Fixed bias b=%g (not fitted)',fixed_bias_b),'FontSize',12);
end
set(gca,'XTick',1,'XTickLabel',{'Bias b'}, ...
    'TickLabelInterpreter','tex','FontSize',15, ...
    'LineWidth',1.5,'TickDir','in');
box on; grid off;

%% Parameter z-score bars: one subplot per sweep, four sweeps per figure
% All active parameters are kept together in each sweep subplot. Figures
% are added automatically in pages of four sweeps using a 2-by-2 layout.
if variation_flags.only_vary_Is_direction
    % Only the three direct spin-current direction offsets are varied.
    active_plot_indices = 9:11;
else
    % Independently sampled primitive device parameters are always active.
    active_plot_indices = 1:8;
    if variation_flags.vary_Is_direction_with_device_variation
        active_plot_indices = [active_plot_indices 9:11];
    end
    if variation_flags.vary_Ic_direction_with_device_variation
        active_plot_indices = [active_plot_indices 12:14];
    end
end

active_parameter_names = mtj_parameter_names(active_plot_indices);
parameter_palette = turbo(numel(active_plot_indices));

% Create the user-designated output folder before saving the figures.
if ~isfolder(barplot_output_folder)
    [folder_created,folder_message] = mkdir(barplot_output_folder);
    if ~folder_created
        error('Could not create bar-plot folder %s: %s', ...
            barplot_output_folder,folder_message);
    end
end

maximum_z = max(abs(mtj_parameter_z_scores(:,active_plot_indices)), ...
    [],'all');
if ~isfinite(maximum_z) || maximum_z==0
    maximum_z = 1;
end
z_limits = 1.10*[-maximum_z maximum_z];

maximum_sweeps_per_deviation_figure = 4;
number_of_deviation_figures = ceil( ...
    Nsweep/maximum_sweeps_per_deviation_figure);
figMTJDeviations = gobjects(number_of_deviation_figures,1);

for deviation_figure_index = 1:number_of_deviation_figures
    first_sweep_in_figure = ...
        (deviation_figure_index-1)*maximum_sweeps_per_deviation_figure+1;
    last_sweep_in_figure = min( ...
        deviation_figure_index*maximum_sweeps_per_deviation_figure,Nsweep);

    figMTJDeviations(deviation_figure_index) = figure('Color','w', ...
        'Name',sprintf('Parameter deviations - variations %d to %d', ...
        first_sweep_in_figure,last_sweep_in_figure), ...
        'NumberTitle','off');
    clf(figMTJDeviations(deviation_figure_index));
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for sweep_index = first_sweep_in_figure:last_sweep_in_figure
        deviation_axes = nexttile;
        deviation_bars = bar(deviation_axes, ...
            1:numel(active_plot_indices), ...
            mtj_parameter_z_scores(sweep_index,active_plot_indices), ...
            'FaceColor','flat','EdgeColor','flat');
        deviation_bars.CData = parameter_palette;
        hold(deviation_axes,'on');
        yline(deviation_axes,0,'k-','LineWidth',1.0, ...
            'HandleVisibility','off');

        % ylabel(deviation_axes,'Z-score','FontSize',11);
        set(deviation_axes,'XTick',1:numel(active_plot_indices), ...
            'XTickLabel',active_parameter_names, ...
            'TickLabelInterpreter','tex','FontSize',12, ...
            'LineWidth',1.2,'TickDir','in');
        title(deviation_axes,sprintf('Variation %d',sweep_index), ...
            'FontSize',10,'Color','k');
        xtickangle(deviation_axes,90);
        xlim(deviation_axes,[0.25 numel(active_plot_indices)+0.75]);
        ylim(deviation_axes,z_limits);
        box(deviation_axes,'on');
        grid(deviation_axes,'off');
    end

    % Save each four-sweep page as a publication-quality image and as an
    % editable MATLAB figure in the user-designated folder.
    drawnow;
    barplot_base_name = sprintf( ...
        'parameter_zscores_sweeps_%03d_to_%03d', ...
        first_sweep_in_figure,last_sweep_in_figure);
    exportgraphics(figMTJDeviations(deviation_figure_index), ...
        fullfile(barplot_output_folder,[barplot_base_name '.png']), ...
        'Resolution',300);
    savefig(figMTJDeviations(deviation_figure_index), ...
        fullfile(barplot_output_folder,[barplot_base_name '.fig']));
end

%% Consolidated parameter z-score bars: all sweeps in one large figure
% Sweep numbers are placed inside the individual axes. Only the first
% subplot column retains y tick labels, and the tiled layout supplies one
% common y-axis label to maximize plotting space.
if ~isscalar(concatenated_barplot_columns) || ...
        ~isfinite(concatenated_barplot_columns) || ...
        concatenated_barplot_columns<1 || ...
        concatenated_barplot_columns~=floor(concatenated_barplot_columns)
    error('concatenated_barplot_columns must be a positive integer.');
end

concatenated_barplot_rows = ceil(Nsweep/concatenated_barplot_columns);
concatenated_figure_width = 1900;
concatenated_figure_height = max(1000,230*concatenated_barplot_rows);
two_sigma_z = 2;
concatenated_z_extent = max(maximum_z,two_sigma_z);
concatenated_z_limits = 1.10*[-concatenated_z_extent ...
    concatenated_z_extent];

figAllSweepDeviations = figure('Color','w','Name', ...
    sprintf('Parameter deviations - all %d variations',Nsweep), ...
    'NumberTitle','off','Units','pixels','Position', ...
    [40 40 concatenated_figure_width concatenated_figure_height]);
clf(figAllSweepDeviations);
all_sweeps_layout = tiledlayout(figAllSweepDeviations, ...
    concatenated_barplot_rows,concatenated_barplot_columns, ...
    'TileSpacing','compact','Padding','compact');
ylabel(all_sweeps_layout,'Z-score','FontSize',16,'FontWeight','bold');
parameter_legend_handles = gobjects(numel(active_plot_indices),1);
first_concatenated_axes = gobjects(1);

for sweep_index = 1:Nsweep
    concatenated_axes = nexttile(all_sweeps_layout);
    concatenated_bars = bar(concatenated_axes, ...
        1:numel(active_plot_indices), ...
        mtj_parameter_z_scores(sweep_index,active_plot_indices), ...
        'FaceColor','flat','EdgeColor','flat');
    concatenated_bars.CData = parameter_palette;
    hold(concatenated_axes,'on');

    % Dummy square markers from the first axes provide one legend entry for
    % every parameter while using exactly the corresponding bar color.
    if sweep_index==1
        first_concatenated_axes = concatenated_axes;
        for parameter_index = 1:numel(active_plot_indices)
            parameter_legend_handles(parameter_index) = plot( ...
                concatenated_axes,nan,nan,'s','LineStyle','none', ...
                'MarkerSize',12,'MarkerFaceColor', ...
                parameter_palette(parameter_index,:), ...
                'MarkerEdgeColor',parameter_palette(parameter_index,:));
        end
    end

    % Light-grey dashed reference limits at two standard deviations.
    yline(concatenated_axes,two_sigma_z,'--','Color', ...
        [0.78 0.78 0.78],'LineWidth',0.9,'HandleVisibility','off');
    yline(concatenated_axes,-two_sigma_z,'--','Color', ...
        [0.78 0.78 0.78],'LineWidth',0.9,'HandleVisibility','off');
    yline(concatenated_axes,0,'k-','LineWidth',0.8, ...
        'HandleVisibility','off');

    set(concatenated_axes,'XTick',[], ...
        'TickLabelInterpreter','tex','FontSize',11, ...
        'LineWidth',0.9,'TickDir','in');
    xlim(concatenated_axes,[0.25 numel(active_plot_indices)+0.75]);
    ylim(concatenated_axes,concatenated_z_limits);

    % Keep numeric y ticks only in the first column. All subplots still use
    % exactly the same z-score limits and tick positions.
    subplot_column = mod(sweep_index-1,concatenated_barplot_columns)+1;
    if subplot_column~=1
        set(concatenated_axes,'YTickLabel',[]);
    end

    text(concatenated_axes,0.03,0.94,sprintf('Variation %d',sweep_index), ...
        'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top','FontSize',9,'FontWeight','bold', ...
        'Color','k', ...
        'BackgroundColor','w','Margin',1,'Clipping','on');
    box(concatenated_axes,'on');
    grid(concatenated_axes,'off');
end

% One shared color-coded parameter legend replaces all subplot x labels.
parameter_legend = legend(first_concatenated_axes, ...
    parameter_legend_handles,active_parameter_names, ...
    'Orientation','horizontal', ...
    'NumColumns',min(7,numel(active_plot_indices)), ...
    'Interpreter','tex','FontSize',10,'Box','off');
parameter_legend.Layout.Tile = 'south';
parameter_legend.ItemTokenSize = [14 10];

drawnow;
exportgraphics(figAllSweepDeviations,fullfile(barplot_output_folder, ...
    'parameter_zscores_all_sweeps.png'),'Resolution',300);
savefig(figAllSweepDeviations,fullfile(barplot_output_folder, ...
    'parameter_zscores_all_sweeps.fig'));

%% Reproduce ensemble Figure 3: hard-axis spin-current components
representative_parameter_set = ...
    parameter_sets_used(representative_sweep);
representative_time_ns = representative_time*1e9;
representative_Is_magnitude = vecnorm(representative_trace.Is,2,1);

figRepresentativeSpinCurrent = figure('Color','w','Name', ...
    'Figure 3 - Hard-axis spin current','NumberTitle','off');
clf(figRepresentativeSpinCurrent);
representative_Is_component_axes = axes( ...
    figRepresentativeSpinCurrent);
representative_Is_component_handles = gobjects(3,1);
representative_Is_component_handles(1) = plot( ...
    representative_Is_component_axes,representative_time_ns, ...
    representative_trace.Is(1,:),'LineWidth',1.4);
hold(representative_Is_component_axes,'on');
representative_Is_component_handles(2) = plot( ...
    representative_Is_component_axes,representative_time_ns, ...
    representative_trace.Is(2,:),'LineWidth',1.4);
representative_Is_component_handles(3) = plot( ...
    representative_Is_component_axes,representative_time_ns, ...
    representative_trace.Is(3,:),'LineWidth',1.4);
xline(representative_Is_component_axes, ...
    representative_parameter_set.tIs_on*1e9,'--', ...
    'HandleVisibility','off');
xline(representative_Is_component_axes, ...
    representative_parameter_set.tIs_off*1e9,'--', ...
    'HandleVisibility','off');
xlabel(representative_Is_component_axes,'Time [ns]');
ylabel(representative_Is_component_axes,'I_s components [A]');
legend(representative_Is_component_axes, ...
    representative_Is_component_handles, ...
    {'I_{s,x}','I_{s,y}','I_{s,z}'},'Location','best', ...
    'Interpreter','tex');
title(representative_Is_component_axes, ...
    'Hard-axis spin-current initialization pulse');
xlim(representative_Is_component_axes, ...
    [representative_time_ns(1) representative_time_ns(end)]);
set(representative_Is_component_axes,'FontSize',13, ...
    'LineWidth',1.4,'TickDir','in');
box(representative_Is_component_axes,'on');
grid(representative_Is_component_axes,'on');

%% Reproduce ensemble Figure 4 from the saved representative trajectory

figRepresentativeMagnetizationIs = figure('Color','w','Name', ...
    'Representative magnetization and direct spin current', ...
    'NumberTitle','off');
clf(figRepresentativeMagnetizationIs);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

representative_m_axes = nexttile;
representative_m_handles = gobjects(3,1);
representative_m_handles(1) = plot(representative_m_axes, ...
    representative_time_ns,representative_trace.m(1,:), ...
    'LineWidth',1.1);
hold(representative_m_axes,'on');
representative_m_handles(2) = plot(representative_m_axes, ...
    representative_time_ns,representative_trace.m(2,:), ...
    'LineWidth',1.1);
representative_m_handles(3) = plot(representative_m_axes, ...
    representative_time_ns,representative_trace.m(3,:), ...
    'LineWidth',1.4);
xline(representative_m_axes,representative_parameter_set.tIs_on*1e9, ...
    '--','HandleVisibility','off');
xline(representative_m_axes,representative_parameter_set.tIs_off*1e9, ...
    '--','HandleVisibility','off');
ylabel(representative_m_axes,'m_i');
ylim(representative_m_axes,[-1.05 1.05]);
xlim(representative_m_axes,[representative_time_ns(1) ...
    representative_time_ns(end)]);
legend(representative_m_axes,representative_m_handles, ...
    {'m_x','m_y','m_z'},'Location','best','Interpreter','tex');
title(representative_m_axes,sprintf( ...
    'Representative Variation %d, \\gamma=%+.4g', ...
    representative_sweep,representative_gamma_trace_actual), ...
    'Interpreter','tex');
set(representative_m_axes,'FontSize',13,'LineWidth',1.4,'TickDir','in');
box(representative_m_axes,'on'); grid(representative_m_axes,'on');

representative_Is_axes = nexttile;
plot(representative_Is_axes,representative_time_ns, ...
    representative_Is_magnitude,'LineWidth',1.5);
hold(representative_Is_axes,'on');
xline(representative_Is_axes,representative_parameter_set.tIs_on*1e9, ...
    '--','HandleVisibility','off');
xline(representative_Is_axes,representative_parameter_set.tIs_off*1e9, ...
    '--','HandleVisibility','off');
xlabel(representative_Is_axes,'Time [ns]');
ylabel(representative_Is_axes,'|I_s| [A]');
xlim(representative_Is_axes,[representative_time_ns(1) ...
    representative_time_ns(end)]);
set(representative_Is_axes,'FontSize',13,'LineWidth',1.4, ...
    'TickDir','in');
box(representative_Is_axes,'on'); grid(representative_Is_axes,'on');

%% Is-only probability versus varied direct-spin-current direction angles
if variation_flags.only_vary_Is_direction
    % The split generator puts z-only offsets in the first 25 variations
    % and x/y-only offsets in the last 25 (26-50 for a 50-variation file).
    % Select groups BEFORE sorting angles, and retain original variation IDs
    % so each plotted probability stays paired with its own direction.
    if Nsweep<50
        error(['The split Is-direction plots require at least 50 variations ' ...
            'to select two nonoverlapping groups of 25.']);
    end
    theta_variation_indices = (1:25).';
    phi_variation_indices = (Nsweep-24:Nsweep).';
    direction_curve_colors = turbo( ...
        number_of_direction_probability_curves);
    % Use the actual simulated Oersted currents, not normalized gamma.
    current_labels = arrayfun(@(value) sprintf('I_c = %+.4g mA',value), ...
        Ic_mA(selected_synaptic_gamma_indices),'UniformOutput',false);

    figIsThetaProbability = figure('Color','w','Name', ...
        'Probability versus varied Is polar angle','NumberTitle','off');
    clf(figIsThetaProbability);
    theta_axes = axes(figIsThetaProbability);
    hold(theta_axes,'on');
    theta_handles = gobjects(number_of_direction_probability_curves,1);
    [sorted_theta,theta_group_order] = sort( ...
        Is_initial_theta_deg(theta_variation_indices));
    theta_order = theta_variation_indices(theta_group_order);
    for current_index = 1:number_of_direction_probability_curves
        theta_handles(current_index) = scatter(theta_axes,sorted_theta, ...
            selected_direction_probability(theta_order,current_index), ...
            55,direction_curve_colors(current_index,:),'filled', ...
            'MarkerEdgeColor',direction_curve_colors(current_index,:));
    end
    xlabel(theta_axes,'Spin-current polar angle \theta_{I_s} [deg]');
    ylabel(theta_axes,'P(m_z=+1)');
    title(theta_axes,'Probability versus I_s polar angle (z-only variations)');
    ylim(theta_axes,[-0.05 1.05]);
    set(theta_axes,'FontSize',13,'LineWidth',1.4,'TickDir','in');
    box(theta_axes,'on'); grid(theta_axes,'off');
    legend(theta_axes,theta_handles,current_labels,'Location','best', ...
        'FontSize',10,'Interpreter','tex');

    figIsPhiProbability = figure('Color','w','Name', ...
        'Probability versus varied Is azimuth','NumberTitle','off');
    clf(figIsPhiProbability);
    phi_axes = axes(figIsPhiProbability);
    hold(phi_axes,'on');
    phi_handles = gobjects(number_of_direction_probability_curves,1);
    [sorted_phi,phi_group_order] = sort( ...
        Is_initial_phi_deg(phi_variation_indices));
    phi_order = phi_variation_indices(phi_group_order);
    for current_index = 1:number_of_direction_probability_curves
        phi_handles(current_index) = scatter(phi_axes,sorted_phi, ...
            selected_direction_probability(phi_order,current_index), ...
            55,direction_curve_colors(current_index,:),'filled', ...
            'MarkerEdgeColor',direction_curve_colors(current_index,:));
    end
    xlabel(phi_axes,'Spin-current azimuth \phi_{I_s} [deg]');
    ylabel(phi_axes,'P(m_z=+1)');
    title(phi_axes,'Probability versus I_s azimuth (x/y-only variations)');
    ylim(phi_axes,[-0.05 1.05]);
    set(phi_axes,'FontSize',13,'LineWidth',1.4,'TickDir','in');
    box(phi_axes,'on'); grid(phi_axes,'off');
    legend(phi_axes,phi_handles,current_labels,'Location','best', ...
        'FontSize',10,'Interpreter','tex');
end

%% Effective beta for each device parameter set
% One bar per first-stage sweep fit, in the original parameter-set order.
% These are physical-current beta_eff values in mA^(-1), not z-scores or
% the single second-stage ensemble slope. Nsweep follows the data,
% so a 50-sweep results file produces exactly 50 bars.
figBetaByParameterSet = figure('Color','w','Name', ...
    'Effective beta by device parameter set','NumberTitle','off', ...
    'Units','pixels','Position',[60 80 max(1000,min(1900,30*Nsweep)) 550]);
beta_parameter_axes = axes(figBetaByParameterSet);
beta_parameter_bars = bar(beta_parameter_axes,Sweep,beta_eff_mA_inv, ...
    0.8,'FaceColor','flat','EdgeColor','none');
beta_parameter_bars.CData = colors;
xlabel(beta_parameter_axes,'Device parameter set (variation number)', ...
    'FontSize',18);
ylabel(beta_parameter_axes,'Effective \beta_{eff} (mA^{-1})', ...
    'FontSize',18,'Interpreter','tex');
set(beta_parameter_axes,'XTick',1:Nsweep, ...
    'XTickLabel',arrayfun(@(s) sprintf('Variation %d',s), ...
        1:Nsweep,'UniformOutput',false), ...
    'TickLabelInterpreter','none','XColor','k', ...
    'FontSize',12,'LineWidth',1.4,'TickDir','in');
if Nsweep>10
    xtickangle(beta_parameter_axes,90);
end
xlim(beta_parameter_axes,[0.25 Nsweep+0.75]);
box(beta_parameter_axes,'on');
grid(beta_parameter_axes,'off');

% Use the same user-designated folder as the parameter-deviation plots.
drawnow;
exportgraphics(figBetaByParameterSet,fullfile(barplot_output_folder, ...
    'effective_beta_current_by_parameter_set.png'),'Resolution',300);
savefig(figBetaByParameterSet,fullfile(barplot_output_folder, ...
    'effective_beta_current_by_parameter_set.fig'));

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
        beta_pdf_at_devices,50,colors(beta_pdf_order,:), ...
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

    beta_pdf_plot_file = fullfile(barplot_output_folder, ...
        ['mtj_effective_beta_probability_density' output_suffix '.png']);
    if save_beta_pdf_plot
        drawnow;
        exportgraphics(figBetaPDF,beta_pdf_plot_file,'Resolution',300);
        savefig(figBetaPDF,fullfile(barplot_output_folder, ...
            ['mtj_effective_beta_probability_density' output_suffix '.fig']));
        fprintf('Saved effective-beta PDF plot to:\n%s\n',beta_pdf_plot_file);
    end
    save(fit_file,'sorted_beta_eff','beta_pdf_order','beta_pdf_x', ...
        'beta_pdf_at_devices','beta_pdf_bandwidth_used','beta_pdf_sample_count', ...
        'beta_pdf_estimation_method','beta_pdf_density','-append');
end


fprintf('\nSaved current-scale MTJ fits to:\n%s\n',fit_file);
fprintf('Saved I0 and beta fit table to:\n%s\n',fit_table_file);

%% Local utilities
function probability = current_tanh_probability(parameters,Ic_mA)
% parameters=[I0_mA,b]; Ic is the signed Oersted current, not spin current Is.
if numel(parameters)~=2 || any(~isfinite(parameters)) || parameters(1)==0
    probability = nan(size(Ic_mA));
    return
end
probability = 0.5*(1+tanh(Ic_mA./parameters(1)+parameters(2)));
end

function sse = current_tanh_sse(parameters,Ic_mA,observed,fit_bias_b,fixed_bias_b)
if fit_bias_b
    if numel(parameters)~=2
        sse = realmax('double')/4;
        return
    end
    model_parameters = parameters;
else
    if numel(parameters)~=1
        sse = realmax('double')/4;
        return
    end
    model_parameters = [parameters(1) fixed_bias_b];
end
probability = current_tanh_probability(model_parameters,Ic_mA);
if any(~isfinite(probability))
    sse = realmax('double')/4;
else
    sse = sum((probability-observed).^2);
end
end

function density = gaussian_kernel_density(samples,evaluation_points,bandwidth)
% Normalized Gaussian KDE; works without ksdensity or Statistics Toolbox.
standardized_distance = ...
    (evaluation_points(:)-samples(:).')/bandwidth;
density = mean(exp(-0.5*standardized_distance.^2),2) ...
    /(sqrt(2*pi)*bandwidth);
density = reshape(density,size(evaluation_points));
end

