clc
clearvars
close all

% General A3 bistable sampler integrated with the stochastic Heun method.
% The SDE is dX = f(t,X)dt + An*dW. Because the noise is additive, the
% same Wiener increment is used in the predictor and corrector stages.

Energy = 0;

a = 1;
b = 2;
c = 0;

% Calculate the energy data regardless of whether its figure is displayed,
% so it is always available in A3_Heun_data.mat for later plotting.
P = linspace(-0.9, 0.9, 101);
E = -(a/2).*P.^2 + (b/4).*P.^4 + (c/6).*P.^6;
A3_energy = [P' E'];

if Energy == 1
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
K = input_trace; % Legacy-friendly name used by the Euler-Maruyama script.

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

            driftNow = driftA3(xState, gamma(k), driveNow, b, c);
            predictor = xState + driftNow*tstep + An(z)*dW;
            driftPredictor = driftA3(predictor, gamma(k), driveNext, b, c);
            xState = xState + 0.5*(driftNow + driftPredictor)*tstep ...
                + An(z)*dW;

            secondDerivative = abs(driftSecondDerivativeA3(xState, b, c));
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
        fprintf('A3 Heun: noise %d/%d, gamma %d/%d complete.\n', ...
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
dataFile = fullfile(outputDir, 'A3_Heun_data.mat');
integration_method = 'stochastic Heun predictor-corrector';
save(dataFile, 'p', 'gamma', 'An', 'As', 'a', 'b', 'c', 'Energy', ...
    'P', 'E', 'A3_energy', 'time', 'input_trace', 'K', 'sampled_bits', ...
    'representative_trajectories', 'tstop', 'tstep', 'j_max', 'T', ...
    'tau', 'N', 'initial', 'rngSeed', 'accumulated_global_error', ...
    'accumulated_global_error_by_gamma', 'local_error_definition', ...
    'global_error_definition', ...
    'integration_method', '-v7.3')
fprintf('Saved simulation data to %s\n', dataFile)

function value = driftA3(x, gammaValue, drive, b, c)
value = drive.*x - b.*x.^3 - c.*x.^5 + gammaValue;
end

function value = driftSecondDerivativeA3(x, b, c)
value = -6.*b.*x - 20.*c.*x.^3;
end

function value = driveValue(t, amplitude, T)
phase = mod(t, 2*T);
if phase < T
    value = amplitude*(1 - exp(-phase/1e-2));
elseif phase > T
    value = amplitude*exp(-phase/1e-2);
else
    % MATLAB heaviside(0) is 1/2, matching the original definition.
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
