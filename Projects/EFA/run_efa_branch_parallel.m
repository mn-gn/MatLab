%% EFA SVD analysis for branch/parallel Debye mock data
% Data folder:
%   C:\Users\ming0\Desktop\MATLAB\mock data\branch_parallel_debye_mock_qt_20260518
%
% This script runs the same front/back SVD EFA analysis for branch and
% parallel data, without q-spectrum front/back matching.

clear;
clc;
close all;

projectDir = fileparts(mfilename("fullpath"));
dataDir = fullfile(getenv("USERPROFILE"), ...
    "Desktop", "MATLAB", "mock data", "branch_parallel_debye_mock_qt_20260518");
outputRoot = fullfile(projectDir, "results_branch_parallel");

if ~isfolder(dataDir)
    error("Data folder does not exist: %s", dataDir);
end

if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

qFile = fullfile(dataDir, "q.csv");
tFile = fullfile(dataDir, "t.csv");
datasetFiles = struct( ...
    "name", {"branch", "parallel"}, ...
    "file", { ...
        fullfile(dataDir, "branch_data_q_by_t.csv"), ...
        fullfile(dataDir, "parallel_data_q_by_t.csv")});

q = readNumericMatrix(qFile);
t = readNumericMatrix(tFile);
q = q(:);
t = t(:);

if any(~isfinite(q)) || any(~isfinite(t))
    error("q or t contains NaN/Inf values.");
end
if any(t <= 0)
    error("The time axis contains zero or negative values, so it cannot be shown on a log scale.");
end

summaryRows = table();

for datasetIndex = 1:numel(datasetFiles)
    datasetName = string(datasetFiles(datasetIndex).name);
    signalFile = datasetFiles(datasetIndex).file;

    if ~isfile(signalFile)
        error("Missing signal file: %s", signalFile);
    end

    outputDir = fullfile(outputRoot, datasetName);
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    deltaS = readNumericMatrix(signalFile);
    deltaS = orientSignalMatrix(deltaS, numel(q), numel(t), signalFile);

    if any(~isfinite(deltaS(:)))
        error("%s contains NaN/Inf values.", signalFile);
    end

    localT = t;
    localDeltaS = deltaS;
    if any(diff(localT) <= 0)
        [localT, order] = sort(localT);
        localDeltaS = localDeltaS(:, order);
        warning("%s time values were not strictly increasing. Signal columns were sorted by time.", datasetName);
    end

    opts = struct();
    opts.maxComponents = min([8, size(localDeltaS, 1), size(localDeltaS, 2)]);
    opts.thresholdMode = "spectralGap";
    opts.gapSearchStartComponent = 2;
    opts.gapSearchEndComponent = opts.maxComponents;
    opts.relativeSingularValueThreshold = 1e-2;
    opts.energyCutoff = 0.995;

    efa = efa_front_back_svd(localDeltaS, localT, opts);

    countTable = table( ...
        (1:numel(localT)).', ...
        localT, ...
        efa.frontComponentCount(:), ...
        efa.backComponentCount(:), ...
        efa.frontEnergyRank(:), ...
        efa.backEnergyRank(:), ...
        'VariableNames', { ...
            'ColumnIndex', ...
            'Time', ...
            'FrontThresholdComponentCount', ...
            'BackThresholdComponentCount', ...
            'FrontEnergyRank', ...
            'BackEnergyRank'});

    transitionTable = makeTransitionTable(efa);

    save(fullfile(outputDir, "efa_results.mat"), ...
        "q", "localT", "localDeltaS", "opts", "efa", "countTable", "transitionTable");
    writetable(countTable, fullfile(outputDir, "efa_component_counts.csv"));
    writetable(efa.componentSummary, fullfile(outputDir, "efa_component_summary.csv"));
    writetable(transitionTable, fullfile(outputDir, "efa_transition_times.csv"));
    writetable(efa.fullSvdTable, fullfile(outputDir, "full_svd_summary.csv"));
    writematrix(efa.frontSingularValues, fullfile(outputDir, "front_singular_values.csv"));
    writematrix(efa.backSingularValues, fullfile(outputDir, "back_singular_values.csv"));

    plotSignalMap(q, localT, localDeltaS, datasetName, outputDir);
    plotEfaTraces(localT, efa, datasetName, outputDir);
    plotComponentCounts(localT, countTable, datasetName, outputDir);
    plotTransitionTimes(localT, efa, datasetName, outputDir);

    fprintf("\n%s full-data SVD summary:\n", datasetName);
    disp(efa.fullSvdTable(1:min(height(efa.fullSvdTable), opts.maxComponents), :));
    fprintf("%s threshold mode: %s, threshold: %.6g, estimated full-data rank: %d\n", ...
        datasetName, efa.thresholdMode, efa.threshold, efa.estimatedFullRank);
    fprintf("%s front/back EFA transition summary:\n", datasetName);
    disp(transitionTable);

    summaryRows = [summaryRows; makeDatasetSummaryRow(datasetName, efa, transitionTable)]; %#ok<AGROW>
end

writetable(summaryRows, fullfile(outputRoot, "branch_parallel_efa_summary.csv"));

fprintf("\nSaved branch/parallel EFA outputs to:\n  %s\n", outputRoot);
disp(summaryRows);

function matrix = orientSignalMatrix(matrix, nQ, nT, signalFile)
    matrixSize = size(matrix);
    if isequal(matrixSize, [nQ, nT])
        return;
    end

    if isequal(matrixSize, [nT, nQ])
        matrix = matrix.';
        return;
    end

    error( ...
        "Signal matrix size in %s is %d-by-%d, but expected %d-by-%d (q-by-time) or %d-by-%d.", ...
        signalFile, matrixSize(1), matrixSize(2), nQ, nT, nT, nQ);
end

function matrix = readNumericMatrix(filename)
    matrix = readmatrix(filename, "FileType", "text");
    spaceDelimited = readmatrix(filename, "FileType", "text", "Delimiter", " ");
    if numel(spaceDelimited) > numel(matrix)
        matrix = spaceDelimited;
    end
end

function transitionTable = makeTransitionTable(efa)
    summary = efa.componentSummary;
    rankToReport = min(efa.estimatedFullRank, height(summary));

    component = summary.Component(1:rankToReport);
    forwardAppearanceColumn = summary.FrontFirstColumn(1:rankToReport);
    forwardAppearanceTime = summary.FrontFirstTime(1:rankToReport);
    backwardDisappearanceStartColumn = summary.BackFirstStartColumn(1:rankToReport);
    backwardDisappearanceTime = summary.BackFirstTime(1:rankToReport);
    pairedIntervalIsOrdered = backwardDisappearanceTime >= forwardAppearanceTime;

    transitionTable = table( ...
        component, ...
        forwardAppearanceColumn, ...
        forwardAppearanceTime, ...
        backwardDisappearanceStartColumn, ...
        backwardDisappearanceTime, ...
        pairedIntervalIsOrdered, ...
        'VariableNames', { ...
            'Component', ...
            'ForwardAppearanceColumn', ...
            'ForwardAppearanceTime', ...
            'BackwardDisappearanceStartColumn', ...
            'BackwardDisappearanceTime', ...
            'PairedIntervalIsOrdered'});
end

function summaryRow = makeDatasetSummaryRow(datasetName, efa, transitionTable)
    firstSingularValues = efa.fullSingularValues(1:min(3, numel(efa.fullSingularValues)));
    while numel(firstSingularValues) < 3
        firstSingularValues(end + 1) = NaN; %#ok<AGROW>
    end

    summaryRow = table( ...
        datasetName, ...
        efa.estimatedFullRank, ...
        efa.threshold, ...
        firstSingularValues(1), ...
        firstSingularValues(2), ...
        firstSingularValues(3), ...
        height(transitionTable), ...
        'VariableNames', { ...
            'Dataset', ...
            'EstimatedFullRank', ...
            'Threshold', ...
            'SV1', ...
            'SV2', ...
            'SV3', ...
            'TransitionRows'});
end

function plotSignalMap(q, t, deltaS, datasetName, outputDir)
    fig = figure("Color", "w", "Name", datasetName + " q-time signal map");
    surf(t(:).', q(:), deltaS, "EdgeColor", "none");
    view(2);
    axis tight;
    set(gca, "XScale", "log", "Layer", "top");
    colormap(parula);
    colorbar;
    xlabel("Time");
    ylabel("q");
    title(datasetName + " q-time signal map");
    exportgraphics(fig, fullfile(outputDir, "signal_map_logtime.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "signal_map_logtime.fig"));
end

function plotEfaTraces(t, efa, datasetName, outputDir)
    fig = figure("Color", "w", "Name", datasetName + " forward/backward EFA singular values");
    layout = tiledlayout(fig, 2, 1, "TileSpacing", "compact", "Padding", "compact");
    colors = lines(efa.maxComponents);
    labels = compose("SV%d", 1:efa.maxComponents);
    thresholdLog = log10(efa.threshold);

    nexttile(layout);
    hold on;
    for componentIndex = 1:efa.maxComponents
        plot(t, efa.frontLog10SingularValues(componentIndex, :), "-o", ...
            "Color", colors(componentIndex, :), ...
            "LineWidth", 1.2, ...
            "MarkerSize", 3);
    end
    yline(thresholdLog, "--k", "Threshold");
    addTransitionLines(gca, efa.componentSummary.FrontFirstTime, efa.estimatedFullRank, "forward");
    set(gca, "XScale", "log");
    grid on;
    xlabel("Time");
    ylabel("log10 singular value");
    title(datasetName + " forward EFA: SVD of columns 1:k");
    legend(labels, "Location", "eastoutside");

    nexttile(layout);
    hold on;
    for componentIndex = 1:efa.maxComponents
        plot(t, efa.backLog10SingularValues(componentIndex, :), "-o", ...
            "Color", colors(componentIndex, :), ...
            "LineWidth", 1.2, ...
            "MarkerSize", 3);
    end
    yline(thresholdLog, "--k", "Threshold");
    addTransitionLines(gca, efa.componentSummary.BackFirstTime, efa.estimatedFullRank, "backward");
    set(gca, "XScale", "log");
    grid on;
    xlabel("Time");
    ylabel("log10 singular value");
    title(datasetName + " backward EFA: SVD of columns k:end");
    legend(labels, "Location", "eastoutside");

    exportgraphics(fig, fullfile(outputDir, "efa_front_back_singular_values.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_front_back_singular_values.fig"));
end

function addTransitionLines(ax, transitionTimes, estimatedRank, scanDirection)
    rankToPlot = min(estimatedRank, numel(transitionTimes));
    if rankToPlot < 1
        return;
    end

    transitionTimes = transitionTimes(1:rankToPlot);
    finiteMask = isfinite(transitionTimes);
    if ~any(finiteMask)
        return;
    end

    if scanDirection == "forward"
        lineColor = [0.00 0.45 0.25];
    else
        lineColor = [0.70 0.10 0.10];
    end

    yLimits = ylim(ax);
    for componentIndex = find(finiteMask(:)).'
        xline(ax, transitionTimes(componentIndex), ":", ...
            "Color", lineColor, ...
            "LineWidth", 0.9);
        text(ax, transitionTimes(componentIndex), yLimits(2), sprintf("  PC%d", componentIndex), ...
            "Color", lineColor, ...
            "FontSize", 8, ...
            "VerticalAlignment", "top", ...
            "Rotation", 90);
    end
end

function plotComponentCounts(t, countTable, datasetName, outputDir)
    fig = figure("Color", "w", "Name", datasetName + " EFA component counts");
    hold on;
    stairs(t, countTable.FrontThresholdComponentCount, "-o", "LineWidth", 1.3, "MarkerSize", 3);
    stairs(t, countTable.BackThresholdComponentCount, "-s", "LineWidth", 1.3, "MarkerSize", 3);
    stairs(t, countTable.FrontEnergyRank, "--", "LineWidth", 1.1);
    stairs(t, countTable.BackEnergyRank, "--", "LineWidth", 1.1);
    set(gca, "XScale", "log");
    grid on;
    xlabel("Time");
    ylabel("Component count");
    title(datasetName + " number of components found by forward/backward EFA");
    legend( ...
        "Forward threshold count", ...
        "Backward threshold count", ...
        "Forward energy rank", ...
        "Backward energy rank", ...
        "Location", "best");
    exportgraphics(fig, fullfile(outputDir, "efa_component_counts.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_component_counts.fig"));
end

function plotTransitionTimes(t, efa, datasetName, outputDir)
    summary = efa.componentSummary;
    rankToPlot = min(efa.estimatedFullRank, height(summary));
    if rankToPlot < 1
        return;
    end

    component = summary.Component(1:rankToPlot);
    forwardTime = summary.FrontFirstTime(1:rankToPlot);
    backwardTime = summary.BackFirstTime(1:rankToPlot);

    fig = figure("Color", "w", "Name", datasetName + " EFA transition times");
    hold on;

    orderedColor = [0.25 0.25 0.25];
    crossedColor = [0.60 0.60 0.60];
    for rowIndex = 1:rankToPlot
        if isfinite(forwardTime(rowIndex)) && isfinite(backwardTime(rowIndex))
            if backwardTime(rowIndex) >= forwardTime(rowIndex)
                lineStyle = "-";
                lineColor = orderedColor;
            else
                lineStyle = "--";
                lineColor = crossedColor;
            end
            plot([forwardTime(rowIndex), backwardTime(rowIndex)], ...
                [component(rowIndex), component(rowIndex)], ...
                lineStyle, ...
                "Color", lineColor, ...
                "LineWidth", 1.1, ...
                "HandleVisibility", "off");
        end
    end

    scatter(forwardTime, component, 60, [0.00 0.45 0.25], ">", "filled", ...
        "DisplayName", "Forward appearance");
    scatter(backwardTime, component, 60, [0.70 0.10 0.10], "<", "filled", ...
        "DisplayName", "Backward disappearance");

    set(gca, "XScale", "log", "YDir", "reverse", "YTick", component);
    xlim([min(t), max(t)]);
    ylim([0.5, rankToPlot + 0.5]);
    grid on;
    xlabel("Time");
    ylabel("SVD component index");
    title(datasetName + " EFA singular-value transition times");
    legend("Location", "best");

    exportgraphics(fig, fullfile(outputDir, "efa_transition_times.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_transition_times.fig"));
end
