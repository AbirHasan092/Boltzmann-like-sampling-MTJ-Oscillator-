clc
clearvars
close all

% Fit all five Heun simulation results to
%
%     P(gamma) = 0.5*[1 + tanh(beta_eff*gamma)]
%
% One effective beta is fitted independently for every noise value. This
% script requires only core MATLAB functions; no Curve Fitting Toolbox is
% needed. Place this script beside the five Heun .mat files, or change
% dataDir below to the folder containing those files. For each case, the
% script also creates a bar plot of accumulated global-error proxy versus
% noise. The error data must be exported by the updated Heun simulations. If a
% .mat file contains recognizable energy-landscape data, an energy figure
% is created for that case as well. Probability plots include +/-1 binomial
% standard-error bars, computed from the measured probabilities using 2000
% independent trials per point. Five separate SE-versus-gamma figures and
% one combined SE figure are also saved. These quantify finite-trial sampling
% uncertainty, separately from the accumulated numerical-error proxy.
% Four additional probability-only figures group Oscillator/A3 and
% Even-order/Gaussian/Cosh, each with and without +/-1 SE bars.
% Probability data use small square markers and thin capped error bars,
% both matched to the corresponding noise-dependent fitted-curve color.
% SE values are not magnified or multiplied by two.

scriptPath = mfilename('fullpath');
dataDir = fileparts(scriptPath);
outputDir = dataDir;

% Independent binary trials per gamma/noise operating point.
numberOfTrials = 2000;
maximumStandardError = 1/(2*sqrt(numberOfTrials));

caseNames = {
    'A3'
    'Oscillator BSN'
    'Cosh BSN'
    'Even higher-order BSN'
    'Gaussian BSN'
    };

dataFiles = {
    'A3_Heun_data.mat'
    'Oscillator_BSN_Heun_data.mat'
    'Cosh_BSN_Heun_data.mat'
    'Even_order_BSN_Heun_data.mat'
    'Gaussian_BSN_Heun_data.mat'
    };

figurePrefixes = {
    'A3'
    'Oscillator_BSN'
    'Cosh_BSN'
    'Even_order_BSN'
    'Gaussian_BSN'
    };

panelLabels = {'(a)', '(b)', '(c)', '(d)', '(e)'};
fitEquation = 'P(gamma) = 0.5*(1 + tanh(beta_eff*gamma))';
fitResults = struct([]);

for caseIndex = 1:numel(dataFiles)
    dataPath = fullfile(dataDir, dataFiles{caseIndex});
    if ~isfile(dataPath)
        error('TanhFit:MissingDataFile', ...
            ['Cannot find %s. Run the corresponding Heun simulation first, ' ...
             'or change dataDir near the top of this script.'], dataPath)
    end

    % Load the complete file so optional energy variables can be detected.
    simulationData = load(dataPath);
    requiredVariables = {'p', 'gamma', 'An', 'accumulated_global_error'};
    for variableIndex = 1:numel(requiredVariables)
        if ~isfield(simulationData, requiredVariables{variableIndex})
            error('TanhFit:MissingVariable', ...
                '%s does not contain the variable %s.', dataPath, ...
                requiredVariables{variableIndex})
        end
    end

    gammaValues = simulationData.gamma(:);
    probability = simulationData.p;
    noiseValues = simulationData.An(:)';
    accumulatedGlobalError = simulationData.accumulated_global_error(:)';

    % Accept either p(gamma,noise) or its transpose.
    if size(probability, 1) ~= numel(gammaValues)
        if size(probability, 2) == numel(gammaValues)
            probability = probability';
        else
            error('TanhFit:InvalidProbabilitySize', ...
                'The dimensions of p in %s do not match gamma.', dataPath)
        end
    end

    numberOfCurves = size(probability, 2);
    if numel(noiseValues) ~= numberOfCurves
        error('TanhFit:InvalidNoiseSize', ...
            'The number of An values in %s does not match the columns of p.', ...
            dataPath)
    end
    if numel(accumulatedGlobalError) ~= numberOfCurves
        error('TanhFit:InvalidGlobalErrorSize', ...
            ['The number of accumulated_global_error values in %s does not ' ...
             'match the columns of p.'], dataPath)
    end

    if ~isreal(probability) || any(isinf(probability(:))) || ...
            any(probability(:) < 0 | probability(:) > 1)
        error('TanhFit:InvalidProbability', ...
            'Probabilities in %s must lie in [0,1], or be NaN for missing data.', ...
            dataPath)
    end
    % Plug-in binomial SE of each measured probability, not the fitted curve.
    % NaN entries remain missing. At p=0 or p=1 this estimate is zero;
    % it is not a confidence interval or proof of zero underlying uncertainty.
    standardError = sqrt(probability.*(1 - probability)/numberOfTrials);

    % Blue/green/yellow/red ordering mirrors the supplied plot template:
    % lower noise is shown in blue and higher noise in red.
    curveColors = makeCurveColors(numberOfCurves);
    betaEffective = nan(1, numberOfCurves);
    gammaFine = linspace(min(gammaValues), max(gammaValues), 501)';

    figureHandle = figure(caseIndex);
    set(figureHandle, 'Color', 'w', 'Name', ...
        sprintf('%s tanh fits', caseNames{caseIndex}), 'NumberTitle', 'off')
    axesHandle = axes(figureHandle);
    hold(axesHandle, 'on')
    fitLineHandles = gobjects(numberOfCurves, 1);
    legendLabels = cell(numberOfCurves, 1);

    for noiseIndex = 1:numberOfCurves
        measuredProbability = probability(:, noiseIndex);
        valid = isfinite(gammaValues) & isfinite(measuredProbability);
        if nnz(valid) < 3
            warning('TanhFit:InsufficientData', ...
                'Skipping %s at An = %g because fewer than three points are valid.', ...
                caseNames{caseIndex}, noiseValues(noiseIndex))
            continue
        end

        betaEffective(noiseIndex) = fitEffectiveBeta( ...
            gammaValues(valid), measuredProbability(valid));
        fittedProbability = tanhProbability(gammaFine, ...
            betaEffective(noiseIndex));

        fitLineHandles(noiseIndex) = plot(axesHandle, gammaFine, ...
            fittedProbability, '-', 'Color', curveColors(noiseIndex, :), ...
            'LineWidth', 2.4);
        legendLabels{noiseIndex} = sprintf('%.2f', betaEffective(noiseIndex));
    end

    % Draw measured points and SE bars above every fitted curve.
    plotBinomialErrorData(axesHandle, gammaValues, probability, standardError)

    xlabel(axesHandle, '\gamma', 'FontSize', 19)
    ylabel(axesHandle, 'Probability', 'FontSize', 15)
    title(axesHandle, caseNames{caseIndex}, 'FontWeight', 'normal')
    xlim(axesHandle, [min(gammaValues), max(gammaValues)])
    ylim(axesHandle, [0, 1])
    yticks(axesHandle, [0 0.5 1])
    if min(gammaValues) <= -1 && max(gammaValues) >= 1
        xticks(axesHandle, [-1 0 1])
    end
    axis(axesHandle, 'square')
    box(axesHandle, 'on')
    grid(axesHandle, 'off')
    axesHandle.FontSize = 13;
    axesHandle.LineWidth = 1.8;
    axesHandle.TickDir = 'in';

    text(axesHandle, 0.03, 0.94, panelLabels{caseIndex}, ...
        'Units', 'normalized', 'FontSize', 16, ...
        'VerticalAlignment', 'top')

    validLegend = isgraphics(fitLineHandles);
    legendHandle = legend(axesHandle, fitLineHandles(validLegend), ...
        legendLabels(validLegend), 'Location', 'southeast', 'Box', 'off', ...
        'FontSize', 14);
    legendHandle.Title.String = '\beta_{eff}';
    legendHandle.Title.FontSize = 16;

    fitResults(caseIndex).caseName = caseNames{caseIndex};
    fitResults(caseIndex).dataFile = dataFiles{caseIndex};
    fitResults(caseIndex).noiseValues = noiseValues;
    fitResults(caseIndex).betaEffective = betaEffective;
    fitResults(caseIndex).numberOfTrials = numberOfTrials;
    fitResults(caseIndex).standardError = standardError;
    fitResults(caseIndex).standardErrorDefinition = ...
        'sqrt(p.*(1-p)/N), measured binary-outcome probabilities; +/-1 SE';
    fitResults(caseIndex).accumulatedGlobalError = accumulatedGlobalError;
    if isfield(simulationData, 'accumulated_global_error_by_gamma')
        fitResults(caseIndex).accumulatedGlobalErrorByGamma = ...
            simulationData.accumulated_global_error_by_gamma;
    else
        fitResults(caseIndex).accumulatedGlobalErrorByGamma = [];
    end
    if isfield(simulationData, 'local_error_definition')
        fitResults(caseIndex).localErrorDefinition = ...
            simulationData.local_error_definition;
    else
        fitResults(caseIndex).localErrorDefinition = ...
            '0.5*|d^2 f(t,x)/dx^2|*tstep^2 per step';
    end
    if isfield(simulationData, 'global_error_definition')
        fitResults(caseIndex).globalErrorDefinition = ...
            simulationData.global_error_definition;
    else
        fitResults(caseIndex).globalErrorDefinition = ...
            ['sum along each trajectory; maximum over trajectories and ' ...
             'gamma for each noise level'];
    end
    fitResults(caseIndex).gamma = gammaValues;
    fitResults(caseIndex).measuredProbability = probability;
    fitResults(caseIndex).fitGamma = gammaFine;
    fitResults(caseIndex).fittedProbability = ...
        0.5*(1 + tanh(gammaFine*betaEffective));

    fprintf('\n%s\n', caseNames{caseIndex})
    fprintf('An\t\tbeta_eff\tAccumulated global error\n')
    for noiseIndex = 1:numberOfCurves
        fprintf('%.4g\t\t%.6g\t\t%.6g\n', noiseValues(noiseIndex), ...
            betaEffective(noiseIndex), accumulatedGlobalError(noiseIndex))
    end

    figureBaseName = sprintf('%s_tanh_fits', figurePrefixes{caseIndex});
    exportgraphics(figureHandle, ...
        fullfile(outputDir, [figureBaseName '.png']), 'Resolution', 300)
    savefig(figureHandle, fullfile(outputDir, [figureBaseName '.fig']))

    % Separate standard-error curves for every noise level in this case.
    % IDs 21--25 do not overlap existing fits, errors, energy, or composites.
    standardErrorFigure = figure(4*numel(dataFiles) + caseIndex);
    set(standardErrorFigure, 'Color', 'w', 'Name', ...
        sprintf('%s binomial standard error (N = %d)', ...
        caseNames{caseIndex}, numberOfTrials), 'NumberTitle', 'off')
    standardErrorAxes = axes(standardErrorFigure);
    plotStandardErrorPanel(standardErrorAxes, fitResults(caseIndex), ...
        maximumStandardError)
    title(standardErrorAxes, sprintf('%s (N = %d)', ...
        caseNames{caseIndex}, numberOfTrials), 'FontWeight', 'normal')
    axis(standardErrorAxes, 'square')
    standardErrorBaseName = sprintf('%s_standard_error_vs_gamma', ...
        figurePrefixes{caseIndex});
    exportgraphics(standardErrorFigure, ...
        fullfile(outputDir, [standardErrorBaseName '.png']), 'Resolution', 300)
    savefig(standardErrorFigure, ...
        fullfile(outputDir, [standardErrorBaseName '.fig']))

    % Create a separate accumulated-global-error-versus-noise bar plot.
    errorFigureHandle = figure(numel(dataFiles) + caseIndex);
    set(errorFigureHandle, 'Color', 'w', 'Name', ...
        sprintf('%s accumulated global error versus noise', ...
        caseNames{caseIndex}), ...
        'NumberTitle', 'off')
    errorAxesHandle = axes(errorFigureHandle);
    hold(errorAxesHandle, 'on')

    % Use evenly spaced bar positions so the outside labels have sufficient
    % horizontal room. The tick labels still show the actual noise values.
    barPositions = 1:numberOfCurves;
    barHandles = gobjects(numberOfCurves, 1);
    for noiseIndex = 1:numberOfCurves
        faceColor = curveColors(noiseIndex, :);
        % A darkened version of each face color blends the formerly black
        % border into the corresponding bar color.
        edgeColor = 0.62*faceColor;
        barHandles(noiseIndex) = bar(errorAxesHandle, ...
            barPositions(noiseIndex), accumulatedGlobalError(noiseIndex), ...
            0.65, ...
            'FaceColor', faceColor, 'EdgeColor', edgeColor, ...
            'LineWidth', 1.4);
    end

    xlabel(errorAxesHandle, 'Noise value, A_n', 'FontSize', 15)
    ylabel(errorAxesHandle, 'Accumulated global error', 'FontSize', 15)
    title(errorAxesHandle, sprintf('%s accumulated global error', ...
        caseNames{caseIndex}), 'FontWeight', 'normal')
    xticks(errorAxesHandle, barPositions)
    xticklabels(errorAxesHandle, compose('%.3g', noiseValues))
    xlim(errorAxesHandle, [0.25, numberOfCurves + 0.75])

    finiteGlobalError = accumulatedGlobalError(isfinite(accumulatedGlobalError));
    if isempty(finiteGlobalError) || max(finiteGlobalError) <= 0
        errorUpperLimit = 1;
    else
        errorUpperLimit = 1.22*max(finiteGlobalError);
    end
    ylim(errorAxesHandle, [0 errorUpperLimit])
    axis(errorAxesHandle, 'square')
    box(errorAxesHandle, 'on')
    grid(errorAxesHandle, 'on')
    errorAxesHandle.GridAlpha = 0;
    errorAxesHandle.FontSize = 13;
    errorAxesHandle.LineWidth = 1.8;
    errorAxesHandle.TickDir = 'in';

    errorFigureBaseName = sprintf('%s_accumulated_global_error_vs_noise', ...
        figurePrefixes{caseIndex});
    exportgraphics(errorFigureHandle, ...
        fullfile(outputDir, [errorFigureBaseName '.png']), ...
        'Resolution', 300)
    savefig(errorFigureHandle, ...
        fullfile(outputDir, [errorFigureBaseName '.fig']))

    % Plot an energy landscape only when compatible data is present.
    [hasEnergyData, energyState, energyValues, energyXAxisLabel, ...
        energySource] = extractEnergyData(simulationData, caseIndex);
    fitResults(caseIndex).hasEnergyData = hasEnergyData;
    fitResults(caseIndex).energySource = energySource;
    fitResults(caseIndex).energyState = energyState;
    fitResults(caseIndex).energyValues = energyValues;

    if hasEnergyData
        energyFigureHandle = figure(2*numel(dataFiles) + caseIndex);
        set(energyFigureHandle, 'Color', 'w', 'Name', ...
            sprintf('%s energy landscape', caseNames{caseIndex}), ...
            'NumberTitle', 'off')
        energyAxesHandle = axes(energyFigureHandle);
        energyLineHandles = plot(energyAxesHandle, energyState, ...
            energyValues, 'LineWidth', 2.8);

        if size(energyValues, 2) == 1
            energyLineHandles.Color = [0.00 0.750 0.750];
        else
            energyColors = turbo(size(energyValues, 2));
            for energyCurveIndex = 1:size(energyValues, 2)
                energyLineHandles(energyCurveIndex).Color = ...
                    energyColors(energyCurveIndex, :);
            end
            legend(energyAxesHandle, ...
                compose('Energy curve %d', 1:size(energyValues, 2)), ...
                'Location', 'best', 'Box', 'off')
        end

        xlabel(energyAxesHandle, energyXAxisLabel, 'FontSize', 19)
        ylabel(energyAxesHandle, 'Energy', 'FontSize', 15)
        title(energyAxesHandle, sprintf('%s energy landscape', ...
            caseNames{caseIndex}), 'FontWeight', 'normal')
        axis(energyAxesHandle, 'square')
        box(energyAxesHandle, 'on')
        grid(energyAxesHandle, 'off')
        energyAxesHandle.FontSize = 13;
        energyAxesHandle.LineWidth = 1.8;
        energyAxesHandle.TickDir = 'in';

        energyFigureBaseName = sprintf('%s_energy_landscape', ...
            figurePrefixes{caseIndex});
        exportgraphics(energyFigureHandle, ...
            fullfile(outputDir, [energyFigureBaseName '.png']), ...
            'Resolution', 300)
        savefig(energyFigureHandle, ...
            fullfile(outputDir, [energyFigureBaseName '.fig']))
        fprintf('Energy plot created from %s.\n', energySource)
    else
        fprintf('No stored energy data found for %s; energy plot skipped.\n', ...
            caseNames{caseIndex})
    end
end

% Use one common vertical scale for every accumulated global-error plot.
allAccumulatedGlobalError = [fitResults.accumulatedGlobalError];
finiteAllGlobalError = ...
    allAccumulatedGlobalError(isfinite(allAccumulatedGlobalError));
if isempty(finiteAllGlobalError) || max(finiteAllGlobalError) <= 0
    commonGlobalErrorUpperLimit = 1;
else
    commonGlobalErrorUpperLimit = 1.22*max(finiteAllGlobalError);
end

% Update and re-export the five error figures at the common scale.
for caseIndex = 1:numel(dataFiles)
    errorFigureHandle = figure(numel(dataFiles) + caseIndex);
    errorAxesHandle = findobj(errorFigureHandle, 'Type', 'axes');
    if ~isempty(errorAxesHandle)
        ylim(errorAxesHandle(1), [0 commonGlobalErrorUpperLimit])
        yticks(errorAxesHandle(1), ...
            linspace(0, commonGlobalErrorUpperLimit, 3))
        errorFigureBaseName = sprintf( ...
            '%s_accumulated_global_error_vs_noise', ...
            figurePrefixes{caseIndex});
        exportgraphics(errorFigureHandle, ...
            fullfile(outputDir, [errorFigureBaseName '.png']), ...
            'Resolution', 300)
        savefig(errorFigureHandle, ...
            fullfile(outputDir, [errorFigureBaseName '.fig']))
    end
end

% Composite figure 1: Oscillator in the left column and A3 in the right.
summaryFigure1 = figure(3*numel(dataFiles) + 1);
set(summaryFigure1, 'Color', 'w', 'Name', ...
    'Oscillator and A3 summary', 'NumberTitle', 'off', ...
    'Position', [100 60 900 700])
% The outer layout controls vertical spacing; each nested row layout
% independently uses compact horizontal spacing.
layout1 = tiledlayout(summaryFigure1, 2, 1, ...
    'TileSpacing', 'loose', 'Padding', 'compact');
layout1Top = tiledlayout(layout1, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
layout1Top.Layout.Tile = 1;
layout1Middle = tiledlayout(layout1, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
layout1Middle.Layout.Tile = 2;

oscillatorEnergyAxes = nexttile(layout1Top, 1);
plotEnergyPanel(oscillatorEnergyAxes, fitResults(2))
title(oscillatorEnergyAxes, 'Oscillator', 'FontWeight', 'normal', ...
    'Color', [0 0.5 0.5])

a3EnergyAxes = nexttile(layout1Top, 2);
plotEnergyPanel(a3EnergyAxes, fitResults(1))
title(a3EnergyAxes, 'A3 normal form', 'FontWeight', 'normal', ...
    'Color', [0 0.5 0.5])

oscillatorProbabilityAxes = nexttile(layout1Middle, 1);
plotProbabilityPanel(oscillatorProbabilityAxes, fitResults(2), false)

a3ProbabilityAxes = nexttile(layout1Middle, 2);
plotProbabilityPanel(a3ProbabilityAxes, fitResults(1), false)

summaryFigure1BaseName = 'Oscillator_A3_composite';
exportgraphics(summaryFigure1, ...
    fullfile(outputDir, [summaryFigure1BaseName '.png']), 'Resolution', 300)
savefig(summaryFigure1, ...
    fullfile(outputDir, [summaryFigure1BaseName '.fig']))

% Composite figure 2: colored probability data and fits only, without SE bars.
% Requested column order: Even order, Gaussian, Cosh.
summaryCaseOrder = [4 5 3];
summaryFigure2 = figure(3*numel(dataFiles) + 2);
set(summaryFigure2, 'Color', 'w', 'Name', ...
    'Even order, Gaussian, and Cosh summary', 'NumberTitle', 'off', ...
    'Position', [80 100 900 350])
layout2 = tiledlayout(summaryFigure2, 1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for columnIndex = 1:numel(summaryCaseOrder)
    caseIndex = summaryCaseOrder(columnIndex);
    probabilityAxes = nexttile(layout2, columnIndex);
    plotProbabilityPanel(probabilityAxes, fitResults(caseIndex), false)
    title(probabilityAxes, caseNames{caseIndex}, 'FontWeight', 'normal', ...
        'Color', [0 0.5 0.5])
end

summaryFigure2BaseName = 'Even_Gaussian_Cosh_composite';
exportgraphics(summaryFigure2, ...
    fullfile(outputDir, [summaryFigure2BaseName '.png']), 'Resolution', 300)
savefig(summaryFigure2, ...
    fullfile(outputDir, [summaryFigure2BaseName '.fig']))

% Combined view of the five sampling-SE plots, all on the same y scale.
standardErrorSummary = figure(5*numel(dataFiles) + 1);
set(standardErrorSummary, 'Color', 'w', 'Name', ...
    'All five cases: binomial standard error', 'NumberTitle', 'off', ...
    'Position', [80 80 1200 720])
standardErrorLayout = tiledlayout(standardErrorSummary, 2, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for caseIndex = 1:numel(dataFiles)
    standardErrorAxes = nexttile(standardErrorLayout, caseIndex);
    plotStandardErrorPanel(standardErrorAxes, fitResults(caseIndex), ...
        maximumStandardError)
    title(standardErrorAxes, caseNames{caseIndex}, 'FontWeight', 'normal')
end
title(standardErrorLayout, sprintf( ...
    'Binomial standard error of measured probabilities (N = %d)', numberOfTrials))
exportgraphics(standardErrorSummary, ...
    fullfile(outputDir, 'All_Heun_standard_error_composite.png'), 'Resolution', 300)
savefig(standardErrorSummary, ...
    fullfile(outputDir, 'All_Heun_standard_error_composite.fig'))

% Two probability-only groups, each duplicated with and without SE bars.
% Figure IDs 27--30 avoid all existing figures (1--26).
probabilityCaseGroups = {[2 1], [4 5 3]};
probabilityGroupTitles = {'Oscillator and A3', ...
    'Even order, Gaussian, and Cosh'};
probabilityGroupFileNames = { ...
    'Oscillator_A3_probability', ...
    'Even_Gaussian_Cosh_probability'};
probabilityVersionNames = {'with binomial SE', 'without binomial SE'};
probabilityVersionSuffixes = {'_with_binomial_SE', '_without_binomial_SE'};
probabilityPanelTitles = {'A3 normal form', 'Oscillator', 'Cosh', ...
    'Even higher-order', 'Gaussian'};
for groupIndex = 1:numel(probabilityCaseGroups)
    groupCases = probabilityCaseGroups{groupIndex};
    numberOfPanels = numel(groupCases);
    for versionIndex = 1:2
        showErrorBars = versionIndex == 1;
        probabilityFigure = figure(5*numel(dataFiles) + 1 + groupIndex + ...
            (versionIndex - 1)*numel(probabilityCaseGroups));
        set(probabilityFigure, 'Color', 'w', 'Name', ...
            [probabilityGroupTitles{groupIndex} ' probability ' ...
            probabilityVersionNames{versionIndex}], ...
            'NumberTitle', 'off', 'Position', ...
            [100 100 450*numberOfPanels 460])
        probabilityLayout = tiledlayout(probabilityFigure, 1, numberOfPanels, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        for panelIndex = 1:numberOfPanels
            caseIndex = groupCases(panelIndex);
            probabilityAxes = nexttile(probabilityLayout, panelIndex);
            plotProbabilityPanel(probabilityAxes, fitResults(caseIndex), showErrorBars)
            xlabel(probabilityAxes, '\gamma', 'FontSize', 19)
            ylabel(probabilityAxes, 'Probability', 'FontSize', 18)
            title(probabilityAxes, probabilityPanelTitles{caseIndex}, ...
                'FontWeight', 'normal')
            axis(probabilityAxes, 'square')
        end
        % Omit the overall title on figures with binomial SE bars.
        % Keep individual panel titles and the no-error-bar figure titles.
        if ~showErrorBars
            probabilityFigureTitle = sprintf( ...
                'Measured probability and tanh fits; error bars omitted (N = %d)', ...
                numberOfTrials);
            title(probabilityLayout, probabilityFigureTitle, ...
                'FontSize', 13, 'FontWeight', 'normal')
        end
        drawnow
        probabilityFigureBaseName = [probabilityGroupFileNames{groupIndex} ...
            probabilityVersionSuffixes{versionIndex}];
        exportgraphics(probabilityFigure, ...
            fullfile(outputDir, [probabilityFigureBaseName '.png']), ...
            'Resolution', 300)
        savefig(probabilityFigure, ...
            fullfile(outputDir, [probabilityFigureBaseName '.fig']))
    end
end

resultsFile = fullfile(outputDir, 'All_Heun_tanh_fit_results.mat');
save(resultsFile, 'fitResults', 'fitEquation', 'caseNames', 'dataFiles', ...
    'commonGlobalErrorUpperLimit', 'numberOfTrials', 'maximumStandardError')
fprintf('\nSaved all tanh fitting results to %s\n', resultsFile)

function probability = tanhProbability(gamma, beta)
probability = 0.5*(1 + tanh(beta.*gamma));
end

function beta = fitEffectiveBeta(gamma, probability)
gamma = gamma(:);
probability = probability(:);

% Estimate a useful starting point from the central part of the curve.
central = abs(gamma) <= max(0.25, 0.25*max(abs(gamma)));
if nnz(central) >= 2
    lineEstimate = polyfit(gamma(central), probability(central), 1);
    betaInitial = max(2*lineEstimate(1), 1e-3);
else
    betaInitial = 1;
end

% Optimize log(beta), which enforces beta > 0 without a toolbox.
objective = @(logBeta) sum((tanhProbability(gamma, exp(logBeta)) ...
    - probability).^2);
options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-12, ...
    'MaxFunEvals', 2000, 'MaxIter', 2000);
logBeta = fminsearch(objective, log(betaInitial), options);
beta = exp(logBeta);
end

function colors = makeCurveColors(numberOfCurves)
if numberOfCurves == 4
    colors = [
        0.00 0.20 1.00
        0.00 0.65 0.35
        1.00 0.55 0.00
        1.00 0.05 0.05
        ];
else
    colors = turbo(numberOfCurves);
end
end

function valueRange = rangeWithFallback(values)
valueRange = max(values) - min(values);
if valueRange <= 0
    valueRange = max(abs(values));
end
if valueRange <= 0
    valueRange = 1;
end
end

function [available, state, energy, xLabel, source] = ...
        extractEnergyData(data, caseIndex)
% Recognize the energy-data formats used by the supplied BSN scripts.
available = false;
state = [];
energy = [];
source = '';

if caseIndex == 2
    xLabel = '\vartheta';
else
    xLabel = 'P';
end

if isfield(data, 'osc_energy')
    [state, energy, available] = splitEnergyMatrix(data.osc_energy);
    if available
        source = 'osc_energy';
        xLabel = '\vartheta';
    end
end
if ~available && isfield(data, 'A3_energy')
    [state, energy, available] = splitEnergyMatrix(data.A3_energy);
    if available
        source = 'A3_energy';
    end
end
if ~available && isfield(data, 'Compiled')
    [state, energy, available] = splitEnergyMatrix(data.Compiled);
    if available
        source = 'Compiled';
    end
end
if ~available && isfield(data, 'vartheta') && ...
        isfield(data, 'oscillatorEnergy')
    state = data.vartheta(:);
    energy = orientEnergyValues(data.oscillatorEnergy, numel(state));
    available = ~isempty(energy);
    if available
        source = 'vartheta and oscillatorEnergy';
        xLabel = '\vartheta';
    end
end
if ~available && isfield(data, 'epsilon') && isfield(data, 'E2')
    state = data.epsilon(:);
    energy = orientEnergyValues(data.E2, numel(state));
    available = ~isempty(energy);
    if available
        source = 'epsilon and E2';
        xLabel = '\vartheta';
    end
end
if ~available && isfield(data, 'P') && isfield(data, 'E')
    state = data.P(:);
    energy = orientEnergyValues(data.E, numel(state));
    available = ~isempty(energy);
    if available
        source = 'P and E';
    end
end
if ~available && isfield(data, 'energyPosition') && ...
        isfield(data, 'energyValue')
    state = data.energyPosition(:);
    energy = orientEnergyValues(data.energyValue, numel(state));
    available = ~isempty(energy);
    if available
        source = 'energyPosition and energyValue';
    end
end

if available
    finiteRows = isfinite(state) & any(isfinite(energy), 2);
    state = state(finiteRows);
    energy = energy(finiteRows, :);
    available = numel(state) >= 2 && size(energy, 1) == numel(state);
end

if ~available
    state = [];
    energy = [];
    source = '';
end
end

function [state, energy, available] = splitEnergyMatrix(values)
state = [];
energy = [];
available = false;
if ~isnumeric(values) || ~ismatrix(values) || size(values, 2) < 2
    return
end
state = values(:, 1);
energy = values(:, 2:end);
available = numel(state) >= 2;
end

function energy = orientEnergyValues(values, numberOfStates)
energy = [];
if ~isnumeric(values) || isempty(values)
    return
end
if size(values, 1) == numberOfStates
    energy = values;
elseif size(values, 2) == numberOfStates
    energy = values';
end
end

function plotProbabilityPanel(axesHandle, result, showErrorBars)
if nargin < 3
    showErrorBars = true;
end
numberOfCurves = numel(result.noiseValues);
curveColors = makeCurveColors(numberOfCurves);
hold(axesHandle, 'on')
fitLineHandles = gobjects(numberOfCurves, 1);
legendLabels = cell(numberOfCurves, 1);

for noiseIndex = 1:numberOfCurves
    measuredProbability = result.measuredProbability(:, noiseIndex);
    valid = isfinite(result.gamma) & isfinite(measuredProbability);
    fitLineHandles(noiseIndex) = plot(axesHandle, result.fitGamma, ...
        result.fittedProbability(:, noiseIndex), '-', ...
        'Color', curveColors(noiseIndex, :), 'LineWidth', 2.0);
    legendLabels{noiseIndex} = sprintf('%.2f', ...
        result.betaEffective(noiseIndex));
end

% Match measured-point and error-bar colors to the fitted noise-level curves.
plotBinomialErrorData(axesHandle, result.gamma, ...
    result.measuredProbability, result.standardError, showErrorBars)

xlim(axesHandle, [min(result.gamma), max(result.gamma)])
ylim(axesHandle, [0 1])
yticks(axesHandle, [0 0.5 1])
if min(result.gamma) <= -1 && max(result.gamma) >= 1
    xticks(axesHandle, [-1 0 1])
end
box(axesHandle, 'on')
grid(axesHandle, 'off')
axesHandle.FontSize = 18;
axesHandle.LineWidth = 1.8;
axesHandle.TickDir = 'in';

validLegend = isgraphics(fitLineHandles);
legendHandle = legend(axesHandle, fitLineHandles(validLegend), ...
    legendLabels(validLegend), 'Location', 'southeast', 'Box', 'off', ...
    'FontSize', 16);
legendHandle.Title.String = '\beta_{eff}';
legendHandle.Title.FontSize = 16;
end

function plotBinomialErrorData(axesHandle, gamma, probability, standardError, showErrorBars)
% Each measured probability has +/-1 plug-in binomial SE, in probability
% units. The paper's factor of two for average magnetization is not used.
% SE plots use curve-colored squares and bars; no-SE plots use colored circles.
% Measured points have no connecting line; fits are drawn separately.
if nargin < 5
    showErrorBars = true;
end
gamma = gamma(:);
curveColors = makeCurveColors(size(probability, 2));
hold(axesHandle, 'on')
for noiseIndex = 1:size(probability, 2)
    measuredProbability = probability(:, noiseIndex);
    measuredSE = standardError(:, noiseIndex);
    valid = isfinite(gamma) & isfinite(measuredProbability) & isfinite(measuredSE);
    if ~any(valid)
        continue
    end
    if showErrorBars
        errorbar(axesHandle, gamma(valid), measuredProbability(valid), ...
            measuredSE(valid), 'LineStyle', 'none', ...
            'Color', curveColors(noiseIndex, :), ...
            'LineWidth', 0.8, 'CapSize', 6, 'Marker', 's', ...
            'MarkerSize', 3.5, ...
            'MarkerFaceColor', curveColors(noiseIndex, :), ...
            'MarkerEdgeColor', curveColors(noiseIndex, :), ...
            'HandleVisibility', 'off')
    else
        scatter(axesHandle, gamma(valid), measuredProbability(valid), ...
            60, curveColors(noiseIndex, :), 'filled', ...
            'MarkerEdgeColor', curveColors(noiseIndex, :), ...
            'HandleVisibility', 'off')
    end
end
end

function plotStandardErrorPanel(axesHandle, result, maximumStandardError)
numberOfCurves = numel(result.noiseValues);
curveColors = makeCurveColors(numberOfCurves);
hold(axesHandle, 'on')
% Sort only the plotting order; stored SE rows retain the original gamma order.
[sortedGamma, order] = sort(result.gamma);
for noiseIndex = 1:numberOfCurves
    plot(axesHandle, sortedGamma, result.standardError(order, noiseIndex), ...
        '-o', 'Color', curveColors(noiseIndex, :), ...
        'MarkerFaceColor', curveColors(noiseIndex, :), ...
        'MarkerSize', 3, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('A_n = %.3g', result.noiseValues(noiseIndex)))
end
yline(axesHandle, maximumStandardError, '--', ...
    'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
    'DisplayName', sprintf('Max SE = %.4f', maximumStandardError))
xlabel(axesHandle, '\gamma', 'FontSize', 15)
ylabel(axesHandle, 'Standard error of probability', 'FontSize', 13)
xlim(axesHandle, [min(result.gamma), max(result.gamma)])
ylim(axesHandle, [0 1.12*maximumStandardError])
yticks(axesHandle, [0 0.005 0.010 0.015])
if min(result.gamma) <= -1 && max(result.gamma) >= 1
    xticks(axesHandle, [-1 0 1])
end
box(axesHandle, 'on')
grid(axesHandle, 'off')
axesHandle.FontSize = 12;
axesHandle.LineWidth = 1.8;
axesHandle.TickDir = 'in';
legend(axesHandle, 'show', 'Location', 'southoutside', ...
    'NumColumns', 2, 'Box', 'off', 'FontSize', 10)
end

function plotEnergyPanel(axesHandle, result)
if ~result.hasEnergyData
    text(axesHandle, 0.5, 0.5, 'Energy data unavailable', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 11)
    box(axesHandle, 'on')
    axesHandle.XTick = [];
    axesHandle.YTick = [];
    return
end

hold(axesHandle, 'on')
lineHandles = plot(axesHandle, result.energyState, result.energyValues, ...
    'LineWidth', 2.4);
if size(result.energyValues, 2) == 1
    lineHandles.Color = [0.00 0.750 0.750];
else
    energyColors = turbo(size(result.energyValues, 2));
    for curveIndex = 1:size(result.energyValues, 2)
        lineHandles(curveIndex).Color = energyColors(curveIndex, :);
    end
end
box(axesHandle, 'on')
grid(axesHandle, 'off')
axesHandle.FontSize = 18;
axesHandle.LineWidth = 1.8;
axesHandle.TickDir = 'in';
end

function plotGlobalErrorPanel(axesHandle, result, commonUpperLimit)
numberOfCurves = numel(result.noiseValues);
curveColors = makeCurveColors(numberOfCurves);
barPositions = 1:numberOfCurves;
hold(axesHandle, 'on')

for noiseIndex = 1:numberOfCurves
    faceColor = curveColors(noiseIndex, :);
    edgeColor = 0.62*faceColor;
    bar(axesHandle, barPositions(noiseIndex), ...
        result.accumulatedGlobalError(noiseIndex), ...
        0.65, 'FaceColor', faceColor, 'EdgeColor', edgeColor, ...
        'LineWidth', 1.2);
end

xticks(axesHandle, barPositions)
xticklabels(axesHandle, compose('%.3g', result.noiseValues))
xlim(axesHandle, [0.25 numberOfCurves + 0.75])
ylim(axesHandle, [0 commonUpperLimit])
yticks(axesHandle, linspace(0, commonUpperLimit, 3))
box(axesHandle, 'on')
grid(axesHandle, 'off')
axesHandle.FontSize = 18;
axesHandle.LineWidth = 1.8;
axesHandle.TickDir = 'in';

end



