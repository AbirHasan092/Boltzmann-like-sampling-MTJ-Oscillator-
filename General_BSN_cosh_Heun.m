clc
clearvars
close all

% Hyperbolic-cosh BSN sampler integrated with stochastic Heun.

Energy = 0;

a = 5;
b = 1;
c = 2.6;

if Energy == 1
    P = linspace(-1, 1, 101);
    E = a.*(cosh(b.*P) - 1) - c.*P.^2;
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

            % The original simulation passes the time-dependent drive as the
            % coefficient of 2*x; the energy-only value c is preserved above.
            driftNow = driftCosh(xState, gamma(k), a, b, driveNow);
            predictor = xState + driftNow*tstep + An(z)*dW;
            driftPredictor = driftCosh(predictor, gamma(k), a, b, driveNext);
            xState = xState + 0.5*(driftNow + driftPredictor)*tstep ...
                + An(z)*dW;

            secondDerivative = abs(driftSecondDerivativeCosh(xState, a, b));
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
        fprintf('Cosh Heun: noise %d/%d, gamma %d/%d complete.\n', ...
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
dataFile = fullfile(outputDir, 'Cosh_BSN_Heun_data.mat');
integration_method = 'stochastic Heun predictor-corrector';
save(dataFile, 'p', 'gamma', 'An', 'As', 'a', 'b', 'c', 'Energy', ...
    'time', 'input_trace', 'K', 'sampled_bits', ...
    'representative_trajectories', 'tstop', 'tstep', 'j_max', 'T', ...
    'tau', 'N', 'initial', 'rngSeed', 'accumulated_global_error', ...
    'accumulated_global_error_by_gamma', 'local_error_definition', ...
    'global_error_definition', ...
    'integration_method', '-v7.3')
fprintf('Saved simulation data to %s\n', dataFile)

function value = driftCosh(x, gammaValue, a, b, drive)
value = 2.*drive.*x - a.*b.*sinh(b.*x) + gammaValue;
end

function value = driftSecondDerivativeCosh(x, a, b)
value = -a.*b.^3.*sinh(b.*x);
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
