clc
clearvars
close all

% Gaussian-potential BSN sampler integrated with stochastic Heun.

Energy = 1;

a = 1;
b = 0.5;
c = 0.25;

if Energy == 1
    P = linspace(-0.9, 0.9, 101);
    E = -a.*(exp(-(P-b).^2/(2*c^2)) + exp(-(P+b).^2/(2*c^2)));
    Compiled = [P' E']; %#ok<NASGU>

    figure(2)
    plot(P, E, 'LineWidth', 3, 'Color', [1, 0.5, 0])
    ax = gca;
    ax.FontSize = 36;
    ax.LineWidth = 2.5;
    axis square
    box on
    yticks([-0.1 -0.05 0])
end

An = [0.1 0.2 0.3 0.4];
As = 1;
gamma_count = 31;
gamma = linspace(-1.5, 1.5, gamma_count);

tstop = 16;
tstep = 2e-4;
j_max = 1000;
T = 4;
tau = T;
N = tstop/(2*T);
initial = 0;
rngSeed = 1;
rng(rngSeed, 'twister')

nSteps = round(tstop/tstep);
time = (0:nSteps)'*tstep;
nBits = round(tstop/(2*T));
sampleSteps = round(((2*(1:nBits)-1)*T)/tstep);
sampleLookup = zeros(1, nSteps);
sampleLookup(sampleSteps) = 1:nBits;

p = zeros(gamma_count, numel(An));
sampled_bits = zeros(nBits, j_max, gamma_count, numel(An));
representative_trajectories = zeros(nSteps+1, gamma_count, numel(An));
accumulated_global_error_by_gamma = zeros(gamma_count, numel(An));
input_trace = driveTrace(time, As, T);
K = input_trace;

for z = 1:numel(An)
    for k = 1:gamma_count
        xState = initial*ones(1, j_max);
        bits = zeros(nBits, j_max);
        accumulatedTrajectoryError = zeros(1, j_max);
        representative_trajectories(1, k, z) = xState(1);

        for step = 1:nSteps
            tNow = time(step);
            tNext = time(step+1);
            driveNow = driveValue(tNow, As, T);
            driveNext = driveValue(tNext, As, T);
            dW = sqrt(tstep)*randn(1, j_max);

            % The original simulation passes the drive as the Gaussian-center
            % parameter; the energy-only value b is preserved above.
            driftNow = driftGaussian(xState, gamma(k), a, driveNow, c);
            predictor = xState + driftNow*tstep + An(z)*dW;
            driftPredictor = driftGaussian(predictor, gamma(k), a, driveNext, c);
            xState = xState + 0.5*(driftNow + driftPredictor)*tstep ...
                + An(z)*dW;

            secondDerivative = abs(driftSecondDerivativeGaussian( ...
                xState, a, driveNext, c));
            localTruncationError = 0.5.*secondDerivative.*tstep.^2;
            accumulatedTrajectoryError = accumulatedTrajectoryError ...
                + localTruncationError;

            representative_trajectories(step+1, k, z) = xState(1);
            sampleSlot = sampleLookup(step);
            if sampleSlot > 0
                bits(sampleSlot, :) = 0.5*(1 + sign(xState));
            end
        end

        accumulated_global_error_by_gamma(k, z) = ...
            max(accumulatedTrajectoryError);
        sampled_bits(:, :, k, z) = bits;
        p(k, z) = mean(bits(:));
        fprintf('Gaussian Heun: noise %d/%d, gamma %d/%d complete.\n', ...
            z, numel(An), k, gamma_count)
    end
end

accumulated_global_error = max(accumulated_global_error_by_gamma, [], 1);
local_error_definition = '0.5*|d^2 f(t,x)/dx^2|*tstep^2 per step';
global_error_definition = ...
    ['sum of local error along each trajectory; maximum over trajectories ' ...
     'and gamma for each noise level'];

figure(1)
for z = 1:numel(An)
    plot(gamma, p(:, z), 'LineWidth', 1.5)
    hold on
    scatter(gamma, p(:, z), 'filled')
end
axis square
xlabel('\gamma')
ylabel('Probability of bit 1')
legend(compose('An = %.1f', An), 'Location', 'best')

scriptPath = mfilename('fullpath');
outputDir = fileparts(scriptPath);
dataFile = fullfile(outputDir, 'Gaussian_BSN_Heun_data.mat');
integration_method = 'stochastic Heun predictor-corrector';
save(dataFile, 'p', 'gamma', 'An', 'As', 'a', 'b', 'c', 'Energy', ...
    'time', 'input_trace', 'K', 'sampled_bits', ...
    'representative_trajectories', 'tstop', 'tstep', 'j_max', 'T', ...
    'tau', 'N', 'initial', 'rngSeed', 'accumulated_global_error', ...
    'accumulated_global_error_by_gamma', 'local_error_definition', ...
    'global_error_definition', ...
    'integration_method', '-v7.3')
fprintf('Saved simulation data to %s\n', dataFile)

function value = driftGaussian(x, gammaValue, a, center, width)
value = -a.*(((x-center).*exp(-(x-center).^2/(2*width^2))/(width^2)) ...
    + ((x+center).*exp(-(x+center).^2/(2*width^2))/(width^2))) ...
    + gammaValue;
end

function value = driftSecondDerivativeGaussian(x, a, center, width)
left = x - center;
right = x + center;
leftTerm = exp(-left.^2./(2.*width.^2)) .* ...
    (left.^3./width.^6 - 3.*left./width.^4);
rightTerm = exp(-right.^2./(2.*width.^2)) .* ...
    (right.^3./width.^6 - 3.*right./width.^4);
value = -a.*(leftTerm + rightTerm);
end

function value = driveValue(t, amplitude, T)
phase = mod(t, 2*T);
if phase < T
    value = amplitude*(1 - exp(-phase/1e-2));
elseif phase > T
    value = amplitude*exp(-phase/1e-2);
else
    value = amplitude/2;
end
end

function values = driveTrace(t, amplitude, T)
phase = mod(t, 2*T);
values = zeros(size(t));
firstHalf = phase < T;
secondHalf = phase > T;
boundary = ~(firstHalf | secondHalf);
values(firstHalf) = amplitude*(1 - exp(-phase(firstHalf)/1e-2));
values(secondHalf) = amplitude*exp(-phase(secondHalf)/1e-2);
values(boundary) = amplitude/2;
end
