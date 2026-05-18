%% EFA SVD analysis for q-by-time mock data
% Expected input files:
%   <MATLAB root>/mock data/q.csv
%   <MATLAB root>/mock data/t_81points.csv
%   <MATLAB root>/mock data/q_deltaS_noist_10%_81points.csv
%
% The signal matrix is treated as rows = q points, columns = time points.
% Forward EFA uses columns 1:k. Backward EFA uses columns k:end.

clear;
clc;
close all;

projectDir = fileparts(mfilename("fullpath"));
matlabRoot = fileparts(fileparts(projectDir));
desktopDataDir = fullfile(getenv("USERPROFILE"), "Desktop", "MATLAB", "mock data");
defaultDataDir = fullfile(matlabRoot, "mock data");
if isfolder(desktopDataDir)
    dataDir = desktopDataDir;
else
    dataDir = defaultDataDir;
end
outputDir = fullfile(projectDir, "results");

if ~isfolder(outputDir)
    mkdir(outputDir);
end

qFile = fullfile(dataDir, "q.csv");
tFile = fullfile(dataDir, "t_81points.csv");
signalCandidates = [
    fullfile(dataDir, "q_delta_S_noise_10%_81points.csv")
    fullfile(dataDir, "q_deltaS_noist_10%_81points.csv")
    ];
signalFile = firstExistingFile(signalCandidates);

if signalFile == ""
    signalFile = signalCandidates(1);
end

requiredFiles = [qFile, tFile, signalFile];
missingFiles = requiredFiles(~isfile(requiredFiles));
if ~isempty(missingFiles)
    message = "Missing mock data file(s):" + newline + strjoin("  " + missingFiles, newline) + ...
        newline + newline + "Place the CSV files in: " + dataDir;
    error("%s", message);
end

q = readNumericMatrix(qFile);
t = readNumericMatrix(tFile);
deltaS = readNumericMatrix(signalFile);

q = q(:);
t = t(:);
deltaS = orientSignalMatrix(deltaS, numel(q), numel(t), signalFile);

if any(~isfinite(q)) || any(~isfinite(t)) || any(~isfinite(deltaS(:)))
    error("Input data contains NaN or Inf values. Clean the CSV files before EFA.");
end

if any(t <= 0)
    error("The time axis contains zero or negative values, so it cannot be shown on a log scale.");
end

if any(diff(t) <= 0)
    [t, order] = sort(t);
    deltaS = deltaS(:, order);
    warning("Time values were not strictly increasing. Signal columns were sorted by time.");
end

opts = struct();
opts.maxComponents = min([8, size(deltaS, 1), size(deltaS, 2)]);
opts.thresholdMode = "spectralGap";
opts.gapSearchStartComponent = 2;
opts.gapSearchEndComponent = opts.maxComponents;
opts.relativeSingularValueThreshold = 1e-2;
opts.energyCutoff = 0.995;

efa = efa_front_back_svd(deltaS, t, opts);

countTable = table( ...
    (1:numel(t)).', ...
    t, ...
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

save(fullfile(outputDir, "efa_results.mat"), "q", "t", "deltaS", "opts", "efa", "countTable");
writetable(countTable, fullfile(outputDir, "efa_component_counts.csv"));
writetable(efa.componentSummary, fullfile(outputDir, "efa_component_summary.csv"));
writetable(makeTransitionTable(efa), fullfile(outputDir, "efa_transition_times.csv"));
writetable(efa.fullSvdTable, fullfile(outputDir, "full_svd_summary.csv"));
writematrix(efa.frontSingularValues, fullfile(outputDir, "front_singular_values.csv"));
writematrix(efa.backSingularValues, fullfile(outputDir, "back_singular_values.csv"));

fprintf("Full-data SVD summary:\n");
disp(efa.fullSvdTable(1:min(height(efa.fullSvdTable), opts.maxComponents), :));
fprintf("Threshold mode: %s, threshold: %.6g, estimated full-data rank: %d\n", ...
    efa.thresholdMode, efa.threshold, efa.estimatedFullRank);
fprintf("\nForward/backward EFA component summary:\n");
disp(efa.componentSummary);

plotSignalMap(q, t, deltaS, outputDir);
plotEfaTraces(t, efa, outputDir);
plotComponentCounts(t, countTable, outputDir);
plotTransitionTimes(t, efa, outputDir);
spectralMatch = analyzeFrontBackQSpectra(q, deltaS, efa, outputDir);

fprintf("\nFront/back q-spectrum similarity matches:\n");
disp(spectralMatch.matchTable);

fprintf("\nSaved EFA outputs to:\n  %s\n", outputDir);

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

function filename = firstExistingFile(candidates)
    filename = "";
    for candidateIndex = 1:numel(candidates)
        if isfile(candidates(candidateIndex))
            filename = candidates(candidateIndex);
            return;
        end
    end
end

function matrix = readNumericMatrix(filename)
    matrix = readmatrix(filename, "FileType", "text");
    spaceDelimited = readmatrix(filename, "FileType", "text", "Delimiter", " ");
    if numel(spaceDelimited) > numel(matrix)
        matrix = spaceDelimited;
    end
end

function plotSignalMap(q, t, deltaS, outputDir)
    fig = figure("Color", "w", "Name", "q-time signal map");
    surf(t(:).', q(:), deltaS, "EdgeColor", "none");
    view(2);
    axis tight;
    set(gca, "XScale", "log", "Layer", "top");
    colormap(parula);
    colorbar;
    xlabel("Time");
    ylabel("q");
    title("q-time signal map");
    exportgraphics(fig, fullfile(outputDir, "signal_map_logtime.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "signal_map_logtime.fig"));
end

function plotEfaTraces(t, efa, outputDir)
    fig = figure("Color", "w", "Name", "Forward/backward EFA singular values");
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
    title("Forward EFA: SVD of columns 1:k");
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
    title("Backward EFA: SVD of columns k:end");
    legend(labels, "Location", "eastoutside");

    exportgraphics(fig, fullfile(outputDir, "efa_front_back_singular_values.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_front_back_singular_values.fig"));
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

    axes(ax);
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

function plotComponentCounts(t, countTable, outputDir)
    fig = figure("Color", "w", "Name", "EFA component counts");
    hold on;
    stairs(t, countTable.FrontThresholdComponentCount, "-o", "LineWidth", 1.3, "MarkerSize", 3);
    stairs(t, countTable.BackThresholdComponentCount, "-s", "LineWidth", 1.3, "MarkerSize", 3);
    stairs(t, countTable.FrontEnergyRank, "--", "LineWidth", 1.1);
    stairs(t, countTable.BackEnergyRank, "--", "LineWidth", 1.1);
    set(gca, "XScale", "log");
    grid on;
    xlabel("Time");
    ylabel("Component count");
    title("Number of components found by forward/backward EFA");
    legend( ...
        "Forward threshold count", ...
        "Backward threshold count", ...
        "Forward energy rank", ...
        "Backward energy rank", ...
        "Location", "best");
    exportgraphics(fig, fullfile(outputDir, "efa_component_counts.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_component_counts.fig"));
end

function plotTransitionTimes(t, efa, outputDir)
    summary = efa.componentSummary;
    rankToPlot = min(efa.estimatedFullRank, height(summary));
    if rankToPlot < 1
        return;
    end

    component = summary.Component(1:rankToPlot);
    forwardTime = summary.FrontFirstTime(1:rankToPlot);
    backwardTime = summary.BackFirstTime(1:rankToPlot);

    fig = figure("Color", "w", "Name", "EFA transition times");
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
    title("EFA singular-value transition times");
    legend("Location", "best");

    exportgraphics(fig, fullfile(outputDir, "efa_transition_times.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_transition_times.fig"));
end

function spectralMatch = analyzeFrontBackQSpectra(q, deltaS, efa, outputDir)
    summary = efa.componentSummary;
    rankToAnalyze = min(efa.estimatedFullRank, height(summary));

    if rankToAnalyze < 1
        spectralMatch = struct();
        spectralMatch.pairwiseTable = table();
        spectralMatch.matchTable = table();
        return;
    end

    frontSpectra = nan(numel(q), rankToAnalyze);
    backSpectra = nan(numel(q), rankToAnalyze);

    for componentIndex = 1:rankToAnalyze
        frontColumn = summary.FrontFirstColumn(componentIndex);
        if isfinite(frontColumn)
            frontMatrix = deltaS(:, 1:frontColumn);
            [uFront, ~, ~] = svd(frontMatrix, "econ");
            if componentIndex <= size(uFront, 2)
                frontSpectra(:, componentIndex) = normalizeSpectrum(uFront(:, componentIndex));
            end
        end

        backColumn = summary.BackFirstStartColumn(componentIndex);
        if isfinite(backColumn)
            backMatrix = deltaS(:, backColumn:end);
            [uBack, ~, ~] = svd(backMatrix, "econ");
            if componentIndex <= size(uBack, 2)
                backSpectra(:, componentIndex) = normalizeSpectrum(uBack(:, componentIndex));
            end
        end
    end

    [signedCosine, absoluteCosine, signedPearson, absolutePearson] = compareSpectra(frontSpectra, backSpectra);
    [frontMatch, backMatch, matchedAbsCosine] = matchComponents(absoluteCosine);

    pairwiseTable = makePairwiseSimilarityTable( ...
        signedCosine, absoluteCosine, signedPearson, absolutePearson);
    matchTable = makeMatchTable( ...
        summary, frontMatch, backMatch, matchedAbsCosine, signedCosine, signedPearson);

    writetable(pairwiseTable, fullfile(outputDir, "efa_q_spectrum_pairwise_similarity.csv"));
    writetable(matchTable, fullfile(outputDir, "efa_q_spectrum_matches.csv"));
    writematrix(signedCosine, fullfile(outputDir, "efa_q_spectrum_signed_cosine_matrix.csv"));
    writematrix(absoluteCosine, fullfile(outputDir, "efa_q_spectrum_abs_cosine_matrix.csv"));

    plotSpectrumSimilarityHeatmap(absoluteCosine, outputDir);
    plotMatchedSpectra(q, frontSpectra, backSpectra, frontMatch, backMatch, signedCosine, outputDir);

    spectralMatch = struct();
    spectralMatch.frontSpectra = frontSpectra;
    spectralMatch.backSpectra = backSpectra;
    spectralMatch.signedCosine = signedCosine;
    spectralMatch.absoluteCosine = absoluteCosine;
    spectralMatch.signedPearson = signedPearson;
    spectralMatch.absolutePearson = absolutePearson;
    spectralMatch.pairwiseTable = pairwiseTable;
    spectralMatch.matchTable = matchTable;
end

function spectrum = normalizeSpectrum(spectrum)
    spectrum = spectrum(:);
    normValue = norm(spectrum);
    if normValue > 0
        spectrum = spectrum ./ normValue;
    end
end

function [signedCosine, absoluteCosine, signedPearson, absolutePearson] = compareSpectra(frontSpectra, backSpectra)
    nFront = size(frontSpectra, 2);
    nBack = size(backSpectra, 2);

    signedCosine = nan(nFront, nBack);
    absoluteCosine = nan(nFront, nBack);
    signedPearson = nan(nFront, nBack);
    absolutePearson = nan(nFront, nBack);

    for frontIndex = 1:nFront
        frontVector = frontSpectra(:, frontIndex);
        for backIndex = 1:nBack
            backVector = backSpectra(:, backIndex);
            validMask = isfinite(frontVector) & isfinite(backVector);
            if nnz(validMask) < 2
                continue;
            end

            frontValid = frontVector(validMask);
            backValid = backVector(validMask);
            signedCosine(frontIndex, backIndex) = dot(frontValid, backValid) ./ ...
                max(norm(frontValid) * norm(backValid), realmin("double"));
            absoluteCosine(frontIndex, backIndex) = abs(signedCosine(frontIndex, backIndex));

            centeredFront = frontValid - mean(frontValid);
            centeredBack = backValid - mean(backValid);
            signedPearson(frontIndex, backIndex) = dot(centeredFront, centeredBack) ./ ...
                max(norm(centeredFront) * norm(centeredBack), realmin("double"));
            absolutePearson(frontIndex, backIndex) = abs(signedPearson(frontIndex, backIndex));
        end
    end
end

function [frontMatch, backMatch, matchedScores] = matchComponents(scoreMatrix)
    nFront = size(scoreMatrix, 1);
    nBack = size(scoreMatrix, 2);
    nMatch = min(nFront, nBack);

    frontMatch = (1:nMatch).';
    backMatch = nan(nMatch, 1);
    matchedScores = nan(nMatch, 1);

    if nMatch == 0
        return;
    end

    if nMatch <= 8
        permutations = perms(1:nBack);
        bestScore = -inf;
        bestPermutation = permutations(1, :);

        for permutationIndex = 1:size(permutations, 1)
            candidate = permutations(permutationIndex, 1:nMatch);
            candidateScores = scoreMatrix(sub2ind(size(scoreMatrix), (1:nMatch).', candidate(:)));
            score = sum(candidateScores, "omitnan");
            if score > bestScore
                bestScore = score;
                bestPermutation = candidate;
            end
        end

        backMatch = bestPermutation(:);
    else
        remainingBack = true(1, nBack);
        for frontIndex = 1:nMatch
            candidateScores = scoreMatrix(frontIndex, :);
            candidateScores(~remainingBack) = -inf;
            [~, backIndex] = max(candidateScores);
            backMatch(frontIndex) = backIndex;
            remainingBack(backIndex) = false;
        end
    end

    matchedScores = scoreMatrix(sub2ind(size(scoreMatrix), frontMatch, backMatch));
end

function pairwiseTable = makePairwiseSimilarityTable( ...
    signedCosine, absoluteCosine, signedPearson, absolutePearson)

    nFront = size(signedCosine, 1);
    nBack = size(signedCosine, 2);
    rowCount = nFront * nBack;

    frontComponent = zeros(rowCount, 1);
    backComponent = zeros(rowCount, 1);
    signedCosineValues = zeros(rowCount, 1);
    absoluteCosineValues = zeros(rowCount, 1);
    signedPearsonValues = zeros(rowCount, 1);
    absolutePearsonValues = zeros(rowCount, 1);

    rowIndex = 0;
    for frontIndex = 1:nFront
        for backIndex = 1:nBack
            rowIndex = rowIndex + 1;
            frontComponent(rowIndex) = frontIndex;
            backComponent(rowIndex) = backIndex;
            signedCosineValues(rowIndex) = signedCosine(frontIndex, backIndex);
            absoluteCosineValues(rowIndex) = absoluteCosine(frontIndex, backIndex);
            signedPearsonValues(rowIndex) = signedPearson(frontIndex, backIndex);
            absolutePearsonValues(rowIndex) = absolutePearson(frontIndex, backIndex);
        end
    end

    pairwiseTable = table( ...
        frontComponent, ...
        backComponent, ...
        signedCosineValues, ...
        absoluteCosineValues, ...
        signedPearsonValues, ...
        absolutePearsonValues, ...
        'VariableNames', { ...
            'FrontComponent', ...
            'BackComponent', ...
            'SignedCosine', ...
            'AbsoluteCosine', ...
            'SignedPearson', ...
            'AbsolutePearson'});
end

function matchTable = makeMatchTable( ...
    summary, frontMatch, backMatch, matchedAbsCosine, signedCosine, signedPearson)

    nMatch = numel(frontMatch);
    frontAppearanceTime = summary.FrontFirstTime(frontMatch);
    backDisappearanceTime = summary.BackFirstTime(backMatch);
    signedCosineMatch = signedCosine(sub2ind(size(signedCosine), frontMatch, backMatch));
    signedPearsonMatch = signedPearson(sub2ind(size(signedPearson), frontMatch, backMatch));
    interpretation = strings(nMatch, 1);

    for matchIndex = 1:nMatch
        if matchedAbsCosine(matchIndex) >= 0.85
            interpretation(matchIndex) = "same q-spectrum";
        elseif matchedAbsCosine(matchIndex) >= 0.60
            interpretation(matchIndex) = "partly similar";
        else
            interpretation(matchIndex) = "different q-spectrum";
        end
    end

    matchTable = table( ...
        frontMatch, ...
        frontAppearanceTime, ...
        backMatch, ...
        backDisappearanceTime, ...
        signedCosineMatch, ...
        matchedAbsCosine, ...
        signedPearsonMatch, ...
        abs(signedPearsonMatch), ...
        interpretation, ...
        'VariableNames', { ...
            'FrontComponent', ...
            'FrontAppearanceTime', ...
            'MatchedBackComponent', ...
            'BackDisappearanceTime', ...
            'SignedCosine', ...
            'AbsoluteCosine', ...
            'SignedPearson', ...
            'AbsolutePearson', ...
            'Interpretation'});
end

function plotSpectrumSimilarityHeatmap(absoluteCosine, outputDir)
    fig = figure("Color", "w", "Name", "Front/back q-spectrum similarity");
    imagesc(absoluteCosine, [0, 1]);
    axis image;
    colormap(parula);
    colorbar;
    xlabel("Backward component");
    ylabel("Forward component");
    title("q-spectrum similarity between front and back EFA components");

    nFront = size(absoluteCosine, 1);
    nBack = size(absoluteCosine, 2);
    set(gca, "XTick", 1:nBack, "YTick", 1:nFront);
    for frontIndex = 1:nFront
        for backIndex = 1:nBack
            text(backIndex, frontIndex, sprintf("%.2f", absoluteCosine(frontIndex, backIndex)), ...
                "HorizontalAlignment", "center", ...
                "Color", "w", ...
                "FontWeight", "bold");
        end
    end

    exportgraphics(fig, fullfile(outputDir, "efa_q_spectrum_similarity_heatmap.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_q_spectrum_similarity_heatmap.fig"));
end

function plotMatchedSpectra(q, frontSpectra, backSpectra, frontMatch, backMatch, signedCosine, outputDir)
    nMatch = numel(frontMatch);
    if nMatch < 1
        return;
    end

    fig = figure("Color", "w", "Name", "Matched front/back q spectra");
    layout = tiledlayout(fig, nMatch, 1, "TileSpacing", "compact", "Padding", "compact");

    for matchIndex = 1:nMatch
        frontIndex = frontMatch(matchIndex);
        backIndex = backMatch(matchIndex);
        signCorrection = sign(signedCosine(frontIndex, backIndex));
        if signCorrection == 0
            signCorrection = 1;
        end

        nexttile(layout);
        hold on;
        plot(q, frontSpectra(:, frontIndex), "LineWidth", 1.3);
        plot(q, signCorrection .* backSpectra(:, backIndex), "--", "LineWidth", 1.3);
        grid on;
        xlabel("q");
        ylabel("Normalized U(q)");
        title(sprintf("Front PC%d vs Back PC%d, |cos| = %.3f", ...
            frontIndex, backIndex, abs(signedCosine(frontIndex, backIndex))));
        legend("Forward", "Backward sign-corrected", "Location", "best");
    end

    exportgraphics(fig, fullfile(outputDir, "efa_q_spectrum_matched_overlays.png"), "Resolution", 200);
    savefig(fig, fullfile(outputDir, "efa_q_spectrum_matched_overlays.fig"));
end
