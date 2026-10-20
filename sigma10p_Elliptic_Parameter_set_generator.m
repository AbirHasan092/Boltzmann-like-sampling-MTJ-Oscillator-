%% SOT_10percent_Elliptic_Parameter_set_generator.m
% Independent Gaussian primitive parameters for the Debashis/Ostwal
% hard-axis spin-current macrospin model.
%
% Derived quantities are never sampled independently:
%   V        <- diameter_x, diameter_y, tF
%   Ns       <- Ms and V, using Ns=Ms*V/muB
%   Nx,Ny,Nz <- diameter_x, diameter_y, tF
%
% Nturns, z_offset, and the Is pulse timing are fixed controls.

clear; clc;

%% User settings
Nsweep = 50;
parameter_seed = 25680;
% Three variation-mode flags:
% 1) Add Is-direction variation to the normal device-parameter variation.
vary_Is_direction_with_device_variation = false;
% 2) Add Ic/Oersted-direction variation to device-parameter variation.
vary_Ic_direction_with_device_variation = true;
% 3) Vary only the Is direction; hold every other parameter at baseline.
%    First half: z offsets only (variations 1-25 when Nsweep=50).
%    Second half: x/y offsets only (variations 26-50 when Nsweep=50).
%    This mode takes precedence over the first two flags.
only_vary_Is_direction = true;

if only_vary_Is_direction && mod(Nsweep,2)~=0
    error('Nsweep must be even to split Is-only variations equally.');
end

script_dir = fileparts(mfilename('fullpath'));
output_file = fullfile(script_dir, ...
    'INI_parameter_sets_elliptic.mat');

% Independent relative 1-sigma spreads for primitive sampled quantities.
rel_sigma.Ms         = 0.10;
rel_sigma.alpha      = 0.10;
rel_sigma.temperature = 0.10;
rel_sigma.diameter   = 0.10;
rel_sigma.thickness  = 0.10;
rel_sigma.Hk         = 0.10;
rel_sigma.ring_radius = 0.10;

% Independent zero-mean Gaussian offsets for the +x spin-current hard axis
% and the +z Oersted-ring field direction.
offset_sigma.Is = [0.05 0.05 0.05];
offset_sigma.Oe = [0.05 0.05 0.05];

%% Debashis/Ostwal baseline means
baseline = struct();
baseline.reference = ['Ostwal, Debashis et al., Scientific Reports 8, ' ...
    '16689 (2018), doi:10.1038/s41598-018-34996-2'];
baseline.Ms_emu_cc = 1000;                % [emu/cm^3]
baseline.Ms        = 1000*baseline.Ms_emu_cc; % [A/m]; 1 emu/cm^3=1000 A/m
baseline.gamma     = 1.760859e11;         % [rad/(s T)]
baseline.diameter_x = 3.0e-7;             % [m], ellipse major axis
baseline.diameter_y = 1.0e-7;             % [m], ellipse minor axis
baseline.tF        = 1.0e-9;              % [m]
baseline.V         = pi*(baseline.diameter_x/2)*(baseline.diameter_y/2)* ...
    baseline.tF;                           % [m^3]
muB = 9.2740100783e-24;                    % [A m^2]
baseline.Ns        = baseline.Ms*baseline.V/muB;
[baseline.Nx,baseline.Ny,baseline.Nz] = ellipsoid_demag_factors( ...
    baseline.diameter_x,baseline.diameter_y,baseline.tF);
baseline.Hk_Oe     = 200;                 % [Oe], Debashis 2020 simulation value
baseline.Hk_mT     = baseline.Hk_Oe/10;   % [mT], 200 Oe = 20 mT
baseline.alpha     = 0.1;
baseline.T         = 300;                 % [K]
baseline.tIs_on    = 15e-9;               % [s]
baseline.tIs_off   = 30e-9;               % [s]
baseline.Is_offset = [0 0 0];

% The paper reports approximately 12 Oe at 10 mA and a 120-nm SiO2
% isolation layer. R_ring follows H=I/(2R) for the ideal loop centre.
baseline.Oe_reference_current = 10e-3;    % [A]
baseline.Oe_reference_field_Oe = 12;      % [Oe]
H_reference = baseline.Oe_reference_field_Oe*1000/(4*pi); % [A/m]
baseline.R_ring = baseline.Oe_reference_current/(2*H_reference); % [m]
baseline.Nturns     = 1;
baseline.z_offset   = 0;              % [m], fixed; not varied
baseline.Oe_offset  = [0 0 0];

% Ideal spin-current direction: +x, the hard axis of a z-easy-axis magnet.
baseline.Is_direction_ideal = [1 0 0];

rng(parameter_seed,'twister');
parameter_sets = repmat(baseline,Nsweep,1);

for s = 1:Nsweep
    p = baseline;
    if ~only_vary_Is_direction
        p.Ms       = positive_gaussian(baseline.Ms,rel_sigma.Ms);
        p.Ms_emu_cc = p.Ms/1000;
        p.alpha    = positive_gaussian(baseline.alpha,rel_sigma.alpha);
        p.T        = positive_gaussian(baseline.T,rel_sigma.temperature);
        p.diameter_x = positive_gaussian( ...
            baseline.diameter_x,rel_sigma.diameter);
        p.diameter_y = positive_gaussian( ...
            baseline.diameter_y,rel_sigma.diameter);
        p.tF       = positive_gaussian(baseline.tF,rel_sigma.thickness);
        p.Hk_Oe    = positive_gaussian(baseline.Hk_Oe,rel_sigma.Hk);
        p.Hk_mT    = p.Hk_Oe/10;
        p.R_ring   = positive_gaussian( ...
            baseline.R_ring,rel_sigma.ring_radius);

        % Recompute every derived quantity from the sampled primitives.
        p.V  = pi*(p.diameter_x/2)*(p.diameter_y/2)*p.tF;
        p.Ns = p.Ms*p.V/muB;
        [p.Nx,p.Ny,p.Nz] = ellipsoid_demag_factors( ...
            p.diameter_x,p.diameter_y,p.tF);
    end

    if only_vary_Is_direction
        if s<=Nsweep/2
            % Vary the z offset only; x and y offsets stay exactly zero.
            p.Is_offset = [0 0 offset_sigma.Is(3)*randn(1,1)];
        else
            % Vary x and y independently; the z offset stays exactly zero.
            p.Is_offset = [offset_sigma.Is(1:2).*randn(1,2) 0];
        end
    elseif vary_Is_direction_with_device_variation
        % Normal device-variation mode retains all three direction offsets.
        p.Is_offset = offset_sigma.Is.*randn(1,3);
    else
        p.Is_offset = [0 0 0];
    end

    if ~only_vary_Is_direction && ...
            vary_Ic_direction_with_device_variation
        p.Oe_offset = offset_sigma.Oe.*randn(1,3);
    else
        p.Oe_offset = [0 0 0];
    end
    parameter_sets(s) = p;
end

variation_flags = struct( ...
    'vary_Is_direction_with_device_variation', ...
        vary_Is_direction_with_device_variation, ...
    'vary_Ic_direction_with_device_variation', ...
        vary_Ic_direction_with_device_variation, ...
    'only_vary_Is_direction',only_vary_Is_direction);

save(output_file,'parameter_sets','baseline','Nsweep','parameter_seed', ...
    'rel_sigma','offset_sigma', ...
    'vary_Is_direction_with_device_variation', ...
    'vary_Ic_direction_with_device_variation', ...
    'only_vary_Is_direction','variation_flags');

fprintf('Saved %d parameter sets to %s\n',Nsweep,output_file);
fprintf('Baseline volume = %.6e m^3\n',baseline.V);
fprintf('Ns definition: Ms*V/muB (dimensionless spin count)\n');
fprintf('Baseline Ns = %.6e\n',baseline.Ns);
fprintf('Baseline Hk = %.1f mT = %.1f Oe\n', ...
    baseline.Hk_mT,baseline.Hk_Oe);
fprintf('Fixed ring z offset = %.1f nm\n',baseline.z_offset*1e9);
if only_vary_Is_direction
    fprintf('Variation mode: Is direction only; all device values fixed.\n');
    fprintf('  Variations 1-%d: z offset only; x/y offsets are zero.\n', ...
        Nsweep/2);
    fprintf('  Variations %d-%d: x/y offsets only; z offset is zero.\n', ...
        Nsweep/2+1,Nsweep);
else
    fprintf('Variation mode: device-parameter variation enabled.\n');
    fprintf('Is direction variation with device variation: %s\n', ...
        on_off(vary_Is_direction_with_device_variation));
    fprintf('Ic/Oersted direction variation with device variation: %s\n', ...
        on_off(vary_Ic_direction_with_device_variation));
end

%% Parameter-variation box plots
% Show only independently sampled primitive parameters. V, Ns, and the
% demagnetization factors are deliberately excluded because they are
% consequences of the sampled geometry and Ms.
names = {'Ms','alpha','T','diameter x','diameter y','tF','Hk','R ring'};
fields = {'Ms','alpha','T','diameter_x','diameter_y','tF', ...
    'Hk_Oe','R_ring'};
X = zeros(Nsweep,numel(fields));
for j = 1:numel(fields)
    values = arrayfun(@(p) p.(fields{j}),parameter_sets).';
    mean_value = baseline.(fields{j});
    X(:,j) = 100*(values-mean_value)/abs(mean_value);
end

O = zeros(Nsweep,6);
for s = 1:Nsweep
    O(s,:) = [parameter_sets(s).Is_offset parameter_sets(s).Oe_offset];
end

figBox = figure('Color','w','Name','Device parameter variation', ...
    'NumberTitle','off');
clf(figBox);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
boxplot(X,'Labels',names,'Whisker',1.5);
ylabel('Variation from mean [%]');
title('Parameter variation (box = 25th to 75th percentile)');
set(gca,'FontSize',11,'LineWidth',1.2,'TickLabelInterpreter','none');
xtickangle(35); box on; grid off;
nexttile;
boxplot(O,'Labels',{'I_s offset x','I_s offset y','I_s offset z', ...
    'Oe offset x','Oe offset y','Oe offset z'}, ...
    'Whisker',1.5);
ylabel('Direction offset');
title('Spin-current and Oersted-field direction-offset variation');
set(gca,'FontSize',11,'LineWidth',1.2,'TickLabelInterpreter','none');
box on; grid off;

function x = positive_gaussian(mu,relative_sigma)
x = mu+abs(mu)*relative_sigma*randn;
while x <= 0
    x = mu+abs(mu)*relative_sigma*randn;
end
end

function [Nx,Ny,Nz] = ellipsoid_demag_factors( ...
    diameter_x,diameter_y,thickness)
% Demagnetization factors for an ellipsoid with semi-axes a,b,c:
% Ni=(abc/2)*integral_0^Inf ds/[(s+ai^2)*Delta(s)].
semi_axes = [diameter_x diameter_y thickness]/2;
if any(~isfinite(semi_axes)) || any(semi_axes<=0)
    error('All ellipsoid dimensions must be positive and finite.');
end

% Scaling by the smallest semi-axis makes the integral dimensionless.
scaled_axes = semi_axes/min(semi_axes);
axis_squared = scaled_axes.^2;
prefactor = prod(scaled_axes)/2;
x_min = log(min(axis_squared))-35;
x_max = log(max(axis_squared))+35;
N = zeros(1,3);
for axis_index = 1:3
    integrand = @(x) demag_log_integrand( ...
        x,axis_squared,axis_index,prefactor);
    N(axis_index) = integral(integrand,x_min,x_max, ...
        'RelTol',1e-9,'AbsTol',1e-12);
end
N = max(N,0);
N = N/sum(N);
Nx=N(1); Ny=N(2); Nz=N(3);
end

function value = demag_log_integrand(x,axis_squared,axis_index,prefactor)
s = exp(x);
delta = sqrt((s+axis_squared(1)).*(s+axis_squared(2)).* ...
    (s+axis_squared(3)));
value = prefactor.*s./((s+axis_squared(axis_index)).*delta);
end

function text_value = on_off(flag_value)
if flag_value
    text_value = 'ENABLED';
else
    text_value = 'DISABLED';
end
end
