function efa = efa_front_back_svd(signalMatrix, timeAxis, opts)
%EFA_FRONT_BACK_SVD Run forward/backward evolving factor analysis by SVD.
%
%   efa = efa_front_back_svd(signalMatrix, timeAxis, opts)
%
%   signalMatrix rows are q points and columns are time points.
%   Forward EFA uses signalMatrix(:, 1:k).
%   Backward EFA uses signalMatrix(:, k:end).
%
%   opts fields:
%       maxComponents                    number of singular-value traces
%       thresholdMode                    "spectralGap" or "relative"
%       gapSearchStartComponent          first gap index for spectral gap
%       gapSearchEndComponent            last gap index for spectral gap
%       relativeSingularValueThreshold   fraction of the largest full-data SV
%       energyCutoff                     cumulative SV^2 energy rank cutoff

    if nargin < 3
        opts = struct();
    end

    [nQ, nT] = size(signalMatrix);
    if numel(timeAxis) ~= nT
        error("timeAxis must have one value per signalMatrix column.");
    end

    maxComponents = getOption(opts, "maxComponents", min(nQ, nT));
    maxComponents = min(maxComponents, min(nQ, nT));

    thresholdMode = string(getOption(opts, "thresholdMode", "spectralGap"));
    gapSearchStartComponent = getOption(opts, "gapSearchStartComponent", 2);
    gapSearchEndComponent = getOption(opts, "gapSearchEndComponent", maxComponents);
    relativeThreshold = getOption(opts, "relativeSingularValueThreshold", 1e-2);
    energyCutoff = getOption(opts, "energyCutoff", 0.995);

    if gapSearchStartComponent < 1
        error("gapSearchStartComponent must be at least 1.");
    end
    if gapSearchEndComponent < gapSearchStartComponent
        error("gapSearchEndComponent must be greater than or equal to gapSearchStartComponent.");
    end
    if relativeThreshold <= 0
        error("relativeSingularValueThreshold must be positive.");
    end
    if energyCutoff <= 0 || energyCutoff > 1
        error("energyCutoff must be in the interval (0, 1].");
    end

    fullSingularValues = svd(signalMatrix, "econ");
    [threshold, estimatedFullRank] = estimateThreshold( ...
        fullSingularValues, ...
        thresholdMode, ...
        relativeThreshold, ...
        gapSearchStartComponent, ...
        gapSearchEndComponent);

    frontSingularValues = nan(maxComponents, nT);
    backSingularValues = nan(maxComponents, nT);
    frontComponentCount = zeros(nT, 1);
    backComponentCount = zeros(nT, 1);
    frontEnergyRank = zeros(nT, 1);
    backEnergyRank = zeros(nT, 1);

    for columnIndex = 1:nT
        frontS = svd(signalMatrix(:, 1:columnIndex), "econ");
        nFront = min(maxComponents, numel(frontS));
        frontSingularValues(1:nFront, columnIndex) = frontS(1:nFront);
        frontComponentCount(columnIndex) = sum(frontS >= threshold);
        frontEnergyRank(columnIndex) = rankByEnergy(frontS, energyCutoff);

        backS = svd(signalMatrix(:, columnIndex:end), "econ");
        nBack = min(maxComponents, numel(backS));
        backSingularValues(1:nBack, columnIndex) = backS(1:nBack);
        backComponentCount(columnIndex) = sum(backS >= threshold);
        backEnergyRank(columnIndex) = rankByEnergy(backS, energyCutoff);
    end

    variance = fullSingularValues.^2;
    if sum(variance) == 0
        explainedVariance = zeros(size(variance));
    else
        explainedVariance = variance ./ sum(variance);
    end
    cumulativeVariance = cumsum(explainedVariance);

    fullSvdTable = table( ...
        (1:numel(fullSingularValues)).', ...
        fullSingularValues, ...
        explainedVariance, ...
        cumulativeVariance, ...
        fullSingularValues >= threshold, ...
        'VariableNames', { ...
            'Component', ...
            'SingularValue', ...
            'ExplainedVariance', ...
            'CumulativeVariance', ...
            'AboveThreshold'});

    componentSummary = summarizeComponents(frontSingularValues, backSingularValues, timeAxis(:), threshold);

    efa = struct();
    efa.timeAxis = timeAxis(:);
    efa.maxComponents = maxComponents;
    efa.thresholdMode = thresholdMode;
    efa.estimatedFullRank = estimatedFullRank;
    efa.gapSearchStartComponent = gapSearchStartComponent;
    efa.gapSearchEndComponent = gapSearchEndComponent;
    efa.relativeSingularValueThreshold = relativeThreshold;
    efa.energyCutoff = energyCutoff;
    efa.threshold = threshold;
    efa.fullSingularValues = fullSingularValues;
    efa.fullSvdTable = fullSvdTable;
    efa.frontSingularValues = frontSingularValues;
    efa.backSingularValues = backSingularValues;
    efa.frontLog10SingularValues = log10WithNan(frontSingularValues);
    efa.backLog10SingularValues = log10WithNan(backSingularValues);
    efa.frontComponentCount = frontComponentCount;
    efa.backComponentCount = backComponentCount;
    efa.frontEnergyRank = frontEnergyRank;
    efa.backEnergyRank = backEnergyRank;
    efa.componentSummary = componentSummary;
end

function [threshold, estimatedRank] = estimateThreshold( ...
    singularValues, thresholdMode, relativeThreshold, gapSearchStartComponent, gapSearchEndComponent)

    if isempty(singularValues) || singularValues(1) == 0
        threshold = realmin("double");
        estimatedRank = 0;
        return;
    end

    switch lower(thresholdMode)
        case "relative"
            threshold = relativeThreshold * singularValues(1);
            estimatedRank = sum(singularValues >= threshold);

        case "spectralgap"
            if numel(singularValues) < 2
                threshold = relativeThreshold * singularValues(1);
                estimatedRank = 1;
                return;
            end

            gapRatios = singularValues(1:end-1) ./ max(singularValues(2:end), realmin("double"));
            searchStart = min(max(1, gapSearchStartComponent), numel(gapRatios));
            searchEnd = min(max(searchStart, gapSearchEndComponent), numel(gapRatios));
            [~, relativeIndex] = max(gapRatios(searchStart:searchEnd));
            estimatedRank = searchStart + relativeIndex - 1;
            threshold = sqrt(singularValues(estimatedRank) * singularValues(estimatedRank + 1));

        otherwise
            error('Unknown thresholdMode "%s". Use "spectralGap" or "relative".', thresholdMode);
    end

    threshold = max(threshold, realmin("double"));
end

function value = getOption(opts, name, defaultValue)
    if isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    else
        value = defaultValue;
    end
end

function rankValue = rankByEnergy(singularValues, energyCutoff)
    if isempty(singularValues) || all(singularValues == 0)
        rankValue = 0;
        return;
    end

    energy = singularValues.^2;
    cumulativeEnergy = cumsum(energy) ./ sum(energy);
    rankValue = find(cumulativeEnergy >= energyCutoff, 1, "first");
end

function values = log10WithNan(values)
    isFinite = isfinite(values);
    values(isFinite) = log10(max(values(isFinite), realmin("double")));
end

function summary = summarizeComponents(frontSingularValues, backSingularValues, timeAxis, threshold)
    maxComponents = size(frontSingularValues, 1);
    nT = size(frontSingularValues, 2);

    component = (1:maxComponents).';
    frontFirstColumn = nan(maxComponents, 1);
    frontFirstTime = nan(maxComponents, 1);
    frontColumnsUsed = nan(maxComponents, 1);
    backFirstStartColumn = nan(maxComponents, 1);
    backFirstTime = nan(maxComponents, 1);
    backColumnsUsed = nan(maxComponents, 1);

    for componentIndex = 1:maxComponents
        frontIndex = find(frontSingularValues(componentIndex, :) >= threshold, 1, "first");
        if ~isempty(frontIndex)
            frontFirstColumn(componentIndex) = frontIndex;
            frontFirstTime(componentIndex) = timeAxis(frontIndex);
            frontColumnsUsed(componentIndex) = frontIndex;
        end

        backIndex = find(backSingularValues(componentIndex, :) >= threshold, 1, "last");
        if ~isempty(backIndex)
            backFirstStartColumn(componentIndex) = backIndex;
            backFirstTime(componentIndex) = timeAxis(backIndex);
            backColumnsUsed(componentIndex) = nT - backIndex + 1;
        end
    end

    summary = table( ...
        component, ...
        frontFirstColumn, ...
        frontFirstTime, ...
        frontColumnsUsed, ...
        backFirstStartColumn, ...
        backFirstTime, ...
        backColumnsUsed, ...
        'VariableNames', { ...
            'Component', ...
            'FrontFirstColumn', ...
            'FrontFirstTime', ...
            'FrontColumnsUsed', ...
            'BackFirstStartColumn', ...
            'BackFirstTime', ...
            'BackColumnsUsed'});
end
