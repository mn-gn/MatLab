%% 05_validate_debye_branch_parallel_open_generator.m
% Build branch and parallel mock data with Debye scattering spectra using
% atomic form factors, then validate open-generator / VarPro separation.

if exist("cfg_override", "var")
    cfg_override_local = cfg_override;
else
    cfg_override_local = [];
end
clearvars -except cfg_override_local;
if isstruct(cfg_override_local)
    cfg_override = cfg_override_local;
end
clc; close all;

cfg = struct();
cfg.randomSeed = 20260517;
cfg.outDir = fullfile(pwd, "debye_branch_parallel_open_generator_results");
cfg.nQ = 190;
cfg.nT = 121;
cfg.nDelay = 82;
cfg.qRange = [0.55, 7.0];       % A^-1
cfg.tRange = [1e-14, 1e-6];     % s
cfg.delayRange = [1e-14, 1e-6]; % s
cfg.noiseRel = 0.004;
cfg.rateLower = 1e5;
cfg.rateUpper = 5e13;
cfg.numStarts = 80;
cfg.maxIter = 900;
cfg.zeroLagWeight = 8.0;
cfg.shortDelayWeight = 2.0;
cfg.shortDelayScale = 2e-9;
cfg.zeroLagSurfaceWeight = 80;
cfg.sparsityWeight = 0.010;
cfg.spectrumNormWeight = 0.75;
cfg.pruneRelativeToSource = 2e-3;

if exist("cfg_override", "var") && isstruct(cfg_override)
    overrideNames = fieldnames(cfg_override);
    for iOverride = 1:numel(overrideNames)
        cfg.(overrideNames{iOverride}) = cfg_override.(overrideNames{iOverride});
    end
end

rng(cfg.randomSeed);
if ~exist(cfg.outDir, "dir")
    mkdir(cfg.outDir);
end

q = linspace(cfg.qRange(1), cfg.qRange(2), cfg.nQ).';
t = logspace(log10(cfg.tRange(1)), log10(cfg.tRange(2)), cfg.nT).';
delay_s = [0, logspace(log10(cfg.delayRange(1)), log10(cfg.delayRange(2)), cfg.nDelay - 1)].';
uSurface = shifted_surface(t, delay_s);
[~, zeroLagIndex] = min(abs(delay_s));
wSurface = surface_weights(numel(t), numel(delay_s), zeroLagIndex, delay_s, cfg);

fprintf("Debye branch/parallel validation\n");
fprintf("q points %d, time points %d, delay points %d\n", numel(q), numel(t), numel(delay_s));

branchData = make_branch_debye_mock(q, t, delay_s, uSurface, cfg);
parallelData = make_parallel_debye_mock(q, t, delay_s, uSurface, cfg);

branchFit = fit_branch_open_generator(branchData.X, branchData.Xmat, wSurface, t, delay_s, ...
    uSurface, zeroLagIndex, branchData.S_norm, cfg);
branchEval = evaluate_fit(branchData, branchFit, t, zeroLagIndex, "branch");

parallelFit = fit_parallel_open_model(parallelData.X, parallelData.Xmat, wSurface, t, delay_s, ...
    uSurface, zeroLagIndex, parallelData.S_norm, cfg);
parallelEval = evaluate_fit(parallelData, parallelFit, t, zeroLagIndex, "parallel");

save(fullfile(cfg.outDir, "workspace_debye_branch_parallel_validation.mat"), ...
    "cfg", "q", "t", "delay_s", "branchData", "parallelData", ...
    "branchFit", "parallelFit", "branchEval", "parallelEval", "-v7.3");

plot_case_results(q, t, delay_s, branchData, branchFit, branchEval, "branch", cfg);
plot_case_results(q, t, delay_s, parallelData, parallelFit, parallelEval, "parallel", cfg);

write_summary(cfg, branchData, branchFit, branchEval, parallelData, parallelFit, parallelEval);

figs = findall(0, "Type", "figure");
for i = 1:numel(figs)
    fig = figs(i);
    name = regexprep(string(fig.Name), "[^\w-]", "_");
    exportgraphics(fig, fullfile(cfg.outDir, sprintf("%02d_%s.png", i, name)), "Resolution", 170);
end

fprintf("\nBRANCH_RESULT zero=%.6g tensor=%.6g meanSpecCorr=%.4f meanConcCorr=%.4f edgeF1=%.3f fTrue=[%.3f %.3f] fEst=[%.3f %.3f]\n", ...
    branchFit.zeroRelerr, branchFit.fullRelerr, mean(branchEval.specCorr), mean(branchEval.concCorr), ...
    branchEval.edgeF1, branchData.truth.branchFraction(1), branchData.truth.branchFraction(2), ...
    branchEval.branchFractionEst(1), branchEval.branchFractionEst(2));
fprintf("BRANCH_RATES_TRUE %s\n", mat2str(branchData.truth.rates, 5));
fprintf("BRANCH_RATES_EST  %s\n", mat2str(branchFit.ratesPruned, 5));

fprintf("\nPARALLEL_RESULT zero=%.6g tensor=%.6g meanSpecCorr=%.4f meanConcCorr=%.4f rateRelErr=%.3g weightAbsErr=%.3g\n", ...
    parallelFit.zeroRelerr, parallelFit.fullRelerr, mean(parallelEval.specCorr), mean(parallelEval.concCorr), ...
    mean(abs(parallelEval.rateEstMatched - parallelEval.rateTrueMatched) ./ parallelEval.rateTrueMatched), ...
    mean(abs(parallelEval.weightEstMatched - parallelEval.weightTrueMatched)));
fprintf("PARALLEL_RATES_TRUE %s\n", mat2str(parallelData.truth.rates, 5));
fprintf("PARALLEL_RATES_EST  %s\n", mat2str(parallelFit.rates, 5));
fprintf("PARALLEL_WEIGHTS_TRUE %s\n", mat2str(parallelData.truth.weights, 5));
fprintf("PARALLEL_WEIGHTS_EST  %s\n", mat2str(parallelFit.weights, 5));
fprintf("OUTDIR=%s\n", cfg.outDir);

function data = make_branch_debye_mock(q, t, delay_s, uSurface, cfg)
    [elements, groundXYZ] = base_structure();
    [S, xyzSet] = select_distinct_debye_spectra(q, elements, groundXYZ, 360, 1101);
    S = S .* [1.00, 0.82, 1.14];

    rates = [7.2e11, 3.0e11, 0, 0, 2.7e8, 4.7e7]; % [k12 k13 k1G k23 k2G k3G]
    K = branch_K(rates);
    Cfull = generator_concentrations(uSurface, K, [1; 0; 0; 0]);
    Csurface = Cfull(:, 1:3);
    XcleanMat = S * Csurface.';
    Xmat = add_relative_noise(XcleanMat, cfg.noiseRel);
    X = reshape(Xmat, numel(q), numel(t), numel(delay_s));
    zeroRows = 1:numel(t);
    data = struct();
    data.kind = "branch";
    data.X = X;
    data.Xmat = Xmat;
    data.XcleanMat = XcleanMat;
    data.S_true = S;
    data.S_norm = sqrt(sum(S.^2, 1));
    data.C_true_surface = Csurface;
    data.C_true_t = Csurface(zeroRows, :);
    data.truth = struct();
    data.truth.rates = rates;
    data.truth.branchFraction = rates(1:2) ./ sum(rates(1:2));
    data.truth.elements = elements;
    data.truth.groundXYZ = groundXYZ;
    data.truth.xyzSet = xyzSet;
end

function data = make_parallel_debye_mock(q, t, delay_s, uSurface, cfg)
    [elements, groundXYZ] = base_structure();
    [S, xyzSet] = select_distinct_debye_spectra(q, elements, groundXYZ, 360, 2202);
    S = S .* [1.10, 0.90, 0.78];

    rates = [8.5e11, 1.2e9, 5.2e7];
    weights = [0.52, 0.31, 0.17];
    Csurface = parallel_concentrations(uSurface, rates, weights);
    XcleanMat = S * Csurface.';
    Xmat = add_relative_noise(XcleanMat, cfg.noiseRel);
    X = reshape(Xmat, numel(q), numel(t), numel(delay_s));
    zeroRows = 1:numel(t);
    data = struct();
    data.kind = "parallel";
    data.X = X;
    data.Xmat = Xmat;
    data.XcleanMat = XcleanMat;
    data.S_true = S;
    data.S_norm = sqrt(sum(S.^2, 1));
    data.C_true_surface = Csurface;
    data.C_true_t = Csurface(zeroRows, :);
    data.truth = struct();
    data.truth.rates = rates;
    data.truth.weights = weights;
    data.truth.elements = elements;
    data.truth.groundXYZ = groundXYZ;
    data.truth.xyzSet = xyzSet;
end

function fit = fit_branch_open_generator(X, Xmat, wSurface, t, delay_s, uSurface, zeroLagIndex, sNormTarget, cfg)
    nT = numel(t);
    qtTarget = X(:, :, zeroLagIndex);
    normZero = norm(qtTarget(:));
    normW = sqrt(sum(sum((Xmat.^2) .* wSurface.')));
    starts = branch_starts(cfg);
    opts = optimset("Display", "off", "MaxIter", cfg.maxIter, ...
        "MaxFunEvals", cfg.maxIter * 50, "TolX", 1e-7, "TolFun", 1e-8);
    bestScore = inf;
    best = [];
    for i = 1:size(starts, 1)
        obj = @(theta) branch_objective(theta, X, Xmat, wSurface, normW, qtTarget, normZero, ...
            t, uSurface, zeroLagIndex, sNormTarget, cfg);
        [thetaOpt, score] = fminsearch(obj, starts(i, :), opts);
        [score, model] = branch_objective(thetaOpt, X, Xmat, wSurface, normW, qtTarget, normZero, ...
            t, uSurface, zeroLagIndex, sNormTarget, cfg);
        if score < bestScore
            bestScore = score;
            best = model;
            best.theta = thetaOpt;
            best.startIndex = i;
            fprintf("branch start %02d improved: score %.6g zero %.5g rates %s\n", ...
                i, score, model.zeroRelerr, mat2str(model.rates, 4));
        end
    end
    fit = best;
    fit.rateNames = ["1->2", "1->3", "1->G", "2->3", "2->G", "3->G"];
    fit.ratesPruned = prune_branch_rates(fit.rates, cfg);
    fit.rateTable = table(cellstr(fit.rateNames.'), fit.rates(:), fit.ratesPruned(:), ...
        'VariableNames', {'transition', 'rate_raw', 'rate_pruned'});
    fit.C_at_t = fit.C_surface((1:nT) + (zeroLagIndex - 1) * nT, :);
end

function [score, model] = branch_objective(theta, X, Xmat, wSurface, normW, qtTarget, normZero, ...
        t, uSurface, zeroLagIndex, sNormTarget, cfg)
    theta = theta(:).';
    lo = log10(cfg.rateLower);
    hi = log10(cfg.rateUpper);
    thetaClip = min(max(theta, lo), hi);
    boundPenalty = sum((theta - thetaClip).^2);
    rates = 10 .^ thetaClip;
    K = branch_K(rates);
    Cfull = generator_concentrations(uSurface, K, [1; 0; 0; 0]);
    C = Cfull(:, 1:3);
    [S, XhatMat] = solve_weighted_spectra(Xmat, C, wSurface);
    Xhat = reshape(XhatMat, size(X));
    qtRecon = Xhat(:, :, zeroLagIndex);
    weightedRelerr = sqrt(sum(sum(((Xmat - XhatMat).^2) .* wSurface.'))) / max(normW, eps);
    zeroRelerr = norm(qtTarget(:) - qtRecon(:)) / max(normZero, eps);
    fullRelerr = norm(X(:) - Xhat(:)) / max(norm(X(:)), eps);
    sourceOut = [sum(rates(1:3)), sum(rates(4:5)), rates(6)];
    srcIdx = [1 1 1 2 2 3];
    sparsePenalty = sum(log1p(rates ./ max(sourceOut(srcIdx), eps))) / numel(rates);
    directGroundPenalty = (rates(3) / max(sourceOut(1), eps))^2;
    sNorm = sqrt(sum(S.^2, 1));
    spectrumNormPenalty = mean(log(max(sNorm, eps) ./ max(sNormTarget, eps)).^2);
    score = weightedRelerr + cfg.zeroLagWeight * zeroRelerr + ...
        cfg.sparsityWeight * sparsePenalty + cfg.spectrumNormWeight * spectrumNormPenalty + ...
        0.05 * directGroundPenalty + 100 * boundPenalty;
    model = struct("score", score, "rates", rates, "K", K, "C_surface", C, "S_q", S, ...
        "Xhat", Xhat, "qt_recon", qtRecon, "weightedRelerr", weightedRelerr, ...
        "zeroRelerr", zeroRelerr, "fullRelerr", fullRelerr);
end

function starts = branch_starts(cfg)
    lo = log10(cfg.rateLower);
    hi = log10(cfg.rateUpper);
    base = [
        5e11 3e11 1e6 1e6 3e8 5e7
        7e11 2e11 1e6 1e6 1e8 1e8
        8e11 1e6 1e6 2e8 1e6 5e7
        4e11 4e11 1e6 3e7 2e8 6e7
        1e12 2e11 1e6 1e6 1e9 1e8
        ];
    starts = zeros(cfg.numStarts, 6);
    nBase = min(size(base, 1), cfg.numStarts);
    starts(1:nBase, :) = log10(base(1:nBase, :));
    for i = nBase+1:cfg.numStarts
        center = starts(randi(nBase), :);
        starts(i, :) = min(max(center + 1.2 * randn(1, 6), lo), hi);
    end
end

function fit = fit_parallel_open_model(X, Xmat, wSurface, t, delay_s, uSurface, zeroLagIndex, sNormTarget, cfg)
    nT = numel(t);
    qtTarget = X(:, :, zeroLagIndex);
    normZero = norm(qtTarget(:));
    normW = sqrt(sum(sum((Xmat.^2) .* wSurface.')));
    starts = parallel_starts(cfg);
    opts = optimset("Display", "off", "MaxIter", cfg.maxIter, ...
        "MaxFunEvals", cfg.maxIter * 45, "TolX", 1e-7, "TolFun", 1e-8);
    bestScore = inf;
    best = [];
    for i = 1:size(starts, 1)
        obj = @(theta) parallel_objective(theta, X, Xmat, wSurface, normW, qtTarget, normZero, ...
            uSurface, zeroLagIndex, sNormTarget, cfg);
        [thetaOpt, score] = fminsearch(obj, starts(i, :), opts);
        [score, model] = parallel_objective(thetaOpt, X, Xmat, wSurface, normW, qtTarget, normZero, ...
            uSurface, zeroLagIndex, sNormTarget, cfg);
        if score < bestScore
            bestScore = score;
            best = model;
            best.theta = thetaOpt;
            best.startIndex = i;
            fprintf("parallel start %02d improved: score %.6g zero %.5g rates %s weights %s\n", ...
                i, score, model.zeroRelerr, mat2str(model.rates, 4), mat2str(model.weights, 4));
        end
    end
    fit = best;
    fit.C_at_t = fit.C_surface((1:nT) + (zeroLagIndex - 1) * nT, :);
end

function [score, model] = parallel_objective(theta, X, Xmat, wSurface, normW, qtTarget, normZero, ...
        uSurface, zeroLagIndex, sNormTarget, cfg)
    theta = theta(:).';
    lo = log10(cfg.rateLower);
    hi = log10(cfg.rateUpper);
    thetaRates = min(max(theta(1:3), lo), hi);
    boundPenalty = sum((theta(1:3) - thetaRates).^2);
    rates = 10 .^ thetaRates;
    weights = softmax_local(theta(4:6));
    C = parallel_concentrations(uSurface, rates, weights);
    [S, XhatMat] = solve_weighted_spectra(Xmat, C, wSurface);
    Xhat = reshape(XhatMat, size(X));
    qtRecon = Xhat(:, :, zeroLagIndex);
    weightedRelerr = sqrt(sum(sum(((Xmat - XhatMat).^2) .* wSurface.'))) / max(normW, eps);
    zeroRelerr = norm(qtTarget(:) - qtRecon(:)) / max(normZero, eps);
    fullRelerr = norm(X(:) - Xhat(:)) / max(norm(X(:)), eps);
    sNorm = sqrt(sum(S.^2, 1));
    spectrumNormPenalty = mean(log(max(sNorm, eps) ./ max(sNormTarget, eps)).^2);
    score = weightedRelerr + cfg.zeroLagWeight * zeroRelerr + ...
        cfg.spectrumNormWeight * spectrumNormPenalty + 100 * boundPenalty;
    model = struct("score", score, "rates", rates, "weights", weights, "C_surface", C, "S_q", S, ...
        "Xhat", Xhat, "qt_recon", qtRecon, "weightedRelerr", weightedRelerr, ...
        "zeroRelerr", zeroRelerr, "fullRelerr", fullRelerr);
end

function starts = parallel_starts(cfg)
    lo = log10(cfg.rateLower);
    hi = log10(cfg.rateUpper);
    baseRates = [
        8e11 1e9 5e7
        1e12 5e8 8e7
        5e11 2e9 4e7
        2e11 8e8 1e8
        ];
    baseWeights = [
        0.50 0.30 0.20
        0.34 0.33 0.33
        0.20 0.50 0.30
        0.60 0.20 0.20
        ];
    starts = zeros(cfg.numStarts, 6);
    nBase = min(size(baseRates, 1), cfg.numStarts);
    for i = 1:nBase
        starts(i, 1:3) = log10(baseRates(i, :));
        starts(i, 4:6) = log(max(baseWeights(i, :), 1e-6));
    end
    for i = nBase+1:cfg.numStarts
        j = randi(nBase);
        starts(i, 1:3) = min(max(log10(baseRates(j, :)) + 1.2 * randn(1, 3), lo), hi);
        w = rand(1, 3);
        w = w ./ sum(w);
        starts(i, 4:6) = log(w) + 0.5 * randn(1, 3);
    end
end

function eval = evaluate_fit(data, fit, t, zeroLagIndex, kind)
    nT = numel(t);
    Cest = fit.C_surface((1:nT) + (zeroLagIndex - 1) * nT, :);
    Ctrue = data.C_true_t;
    Sest = fit.S_q;
    Strue = data.S_true;
    [perm, specCorr, concCorr] = match_components(Strue, Sest, Ctrue, Cest);
    eval = struct();
    eval.perm = perm;
    eval.specCorr = specCorr;
    eval.concCorr = concCorr;
    eval.S_est_matched = Sest(:, perm);
    eval.C_est_matched = Cest(:, perm);
    if kind == "branch"
        truthNonzero = data.truth.rates(:) > 0;
        estNonzero = fit.ratesPruned(:) > 0;
        tp = sum(truthNonzero & estNonzero);
        fp = sum(~truthNonzero & estNonzero);
        fn = sum(truthNonzero & ~estNonzero);
        eval.edgePrecision = tp / max(tp + fp, eps);
        eval.edgeRecall = tp / max(tp + fn, eps);
        eval.edgeF1 = 2 * eval.edgePrecision * eval.edgeRecall / max(eval.edgePrecision + eval.edgeRecall, eps);
        r = fit.ratesPruned;
        eval.branchFractionEst = r(1:2) ./ max(sum(r(1:2)), eps);
    else
        eval.edgeF1 = NaN;
        eval.rateEstMatched = fit.rates(perm);
        eval.weightEstMatched = fit.weights(perm);
        eval.rateTrueMatched = data.truth.rates;
        eval.weightTrueMatched = data.truth.weights;
    end
end

function [permBest, specBest, concBest] = match_components(Strue, Sest, Ctrue, Cest)
    P = perms(1:3);
    scoreBest = -inf;
    permBest = P(1, :);
    specBest = zeros(1, 3);
    concBest = zeros(1, 3);
    for i = 1:size(P, 1)
        p = P(i, :);
        sc = zeros(1, 3);
        cc = zeros(1, 3);
        for k = 1:3
            sc(k) = abs(corr_safe(Strue(:, k), Sest(:, p(k))));
            cc(k) = abs(corr_safe(Ctrue(:, k), Cest(:, p(k))));
        end
        score = mean(sc) + mean(cc);
        if score > scoreBest
            scoreBest = score;
            permBest = p;
            specBest = sc;
            concBest = cc;
        end
    end
end

function r = corr_safe(a, b)
    a = a(:) - mean(a(:));
    b = b(:) - mean(b(:));
    r = (a.' * b) / max(norm(a) * norm(b), eps);
end

function K = branch_K(r)
    k12 = r(1); k13 = r(2); k1g = r(3);
    k23 = r(4); k2g = r(5); k3g = r(6);
    K = [-(k12+k13+k1g), 0, 0, 0; ...
          k12, -(k23+k2g), 0, 0; ...
          k13, k23, -k3g, 0; ...
          k1g, k2g, k3g, 0];
end

function C = generator_concentrations(u, K, c0)
    u = u(:);
    [V, L] = eig(K);
    coeff = V \ c0;
    C = zeros(numel(u), numel(c0));
    for m = 1:numel(c0)
        lambda = real(L(m, m));
        basis = double(u >= 0) .* exp(lambda .* max(u, 0));
        amp = real(V(:, m) * coeff(m));
        C = C + basis * amp.';
    end
    C = max(real(C), 0);
    rowSum = sum(C, 2);
    valid = rowSum > 1 + 1e-8;
    C(valid, :) = C(valid, :) ./ rowSum(valid);
end

function C = parallel_concentrations(u, rates, weights)
    u = u(:);
    C = zeros(numel(u), 3);
    for k = 1:3
        C(:, k) = weights(k) .* double(u >= 0) .* exp(-rates(k) .* max(u, 0));
    end
end

function [S, XhatMat] = solve_weighted_spectra(Xmat, C, w)
    G = C.' * (C .* w);
    ridge = 1e-9 * max(trace(G) / max(size(C, 2), 1), eps);
    S = ((Xmat .* w.') * C) / (G + ridge * eye(size(C, 2)));
    XhatMat = S * C.';
end

function ratesPruned = prune_branch_rates(rates, cfg)
    ratesPruned = rates(:).';
    sources = [1 1 1 2 2 3];
    for src = unique(sources)
        idx = sources == src;
        m = max(ratesPruned(idx));
        if m > 0
            ratesPruned(idx & ratesPruned < cfg.pruneRelativeToSource * m) = 0;
        end
    end
end

function y = add_relative_noise(x, rel)
    sigma = rel * rms(x(:));
    y = x + sigma .* randn(size(x));
end

function u = shifted_surface(t, delay_s)
    nT = numel(t);
    nD = numel(delay_s);
    u = zeros(nT * nD, 1);
    for j = 1:nD
        idx = (1:nT) + (j - 1) * nT;
        u(idx) = t(:) + delay_s(j);
    end
end

function w = surface_weights(nT, nD, zeroLagIndex, delay_s, cfg)
    wDelay = ones(nD, 1) + cfg.shortDelayWeight ./ (1 + delay_s(:) ./ cfg.shortDelayScale);
    wDelay(zeroLagIndex) = wDelay(zeroLagIndex) * cfg.zeroLagSurfaceWeight;
    w = repelem(wDelay, nT);
end

function s = softmax_local(x)
    x = x(:).' - max(x(:));
    ex = exp(x);
    s = ex ./ sum(ex);
end

function [elements, groundXYZ] = base_structure()
    elements = ["C","C","N","O","C","S","C","N","O","C","C","O"];
    groundXYZ = [
        -1.55 -0.80  0.10
        -0.55 -0.10 -0.05
         0.48  0.05  0.12
         1.45  0.40 -0.08
        -0.20  1.10  0.20
         0.85  1.45 -0.12
        -0.92 -1.42  0.36
         0.22 -1.35 -0.28
         1.36 -1.02  0.22
        -1.82  0.82 -0.36
         1.92  1.20  0.24
        -0.18  2.14 -0.18
        ];
end

function [S, xyzSet] = select_distinct_debye_spectra(q, elements, groundXYZ, nCandidates, seed)
    state = rng;
    rng(seed);
    nAtoms = size(groundXYZ, 1);
    I0 = debye_intensity(q, elements, groundXYZ);
    SCand = zeros(numel(q), nCandidates);
    xyzCand = zeros(nAtoms, 3, nCandidates);
    for c = 1:nCandidates
        xyz = random_candidate_xyz(groundXYZ, c);
        s = debye_intensity(q, elements, xyz) - I0;
        if norm(s) < 1e-10
            s = randn(size(s));
        end
        SCand(:, c) = s ./ max(norm(s), eps);
        xyzCand(:, :, c) = xyz;
    end
    corrMat = abs(corr(SCand));
    corrMat(1:nCandidates+1:end) = inf;
    [~, flat] = min(corrMat(:));
    [i1, i2] = ind2sub(size(corrMat), flat);
    scores = inf(nCandidates, 1);
    for c = 1:nCandidates
        if c ~= i1 && c ~= i2
            scores(c) = max(abs(corrMat(c, [i1, i2])));
        end
    end
    [~, i3] = min(scores);
    idx = [i1, i2, i3];
    S = SCand(:, idx);
    xyzSet = xyzCand(:, :, idx);
    rng(state);
end

function xyz = random_candidate_xyz(groundXYZ, candidateIndex)
    nAtoms = size(groundXYZ, 1);
    mode = mod(candidateIndex, 5);
    xyz = groundXYZ;
    switch mode
        case 0
            scale = 0.45 + 2.6 * rand(1, 3);
            xyz = xyz .* scale;
            xyz(:, 3) = xyz(:, 3) + (0.5 + 1.4 * rand) * sin((1:nAtoms).' * (0.7 + 1.8 * rand));
        case 1
            direction = randn(1, 3); direction = direction ./ norm(direction);
            sep = 1.5 + 5.0 * rand;
            xyz(1:floor(nAtoms/2), :) = xyz(1:floor(nAtoms/2), :) - sep * direction;
            xyz(floor(nAtoms/2)+1:end, :) = xyz(floor(nAtoms/2)+1:end, :) + sep * direction;
            xyz = xyz + 0.4 * randn(size(xyz));
        case 2
            theta = 2 * pi * rand;
            Rz = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
            xyz = (xyz * Rz.') .* [0.55 + 2.2 * rand, 0.55 + 2.2 * rand, 0.55 + 2.2 * rand];
            xyz(1:2:end, :) = xyz(1:2:end, :) + (0.8 + 2.5 * rand) * randn(ceil(nAtoms/2), 3);
        case 3
            clusters = [randn(1, 3); randn(1, 3); randn(1, 3)];
            clusters = clusters ./ vecnorm(clusters, 2, 2);
            clusters = clusters .* (2.0 + 4.5 * rand(3, 1));
            for a = 1:nAtoms
                xyz(a, :) = 0.35 * groundXYZ(a, :) + clusters(mod(a-1, 3)+1, :) + 0.35 * randn(1, 3);
            end
        otherwise
            xyz(:, 1) = (0.35 + 3.0 * rand) * xyz(:, 1);
            xyz(:, 2) = (0.35 + 3.0 * rand) * xyz(:, 2);
            xyz(:, 3) = (0.35 + 3.0 * rand) * xyz(:, 3);
            xyz = xyz + linspace(-2.5, 2.5, nAtoms).' .* (randn(1, 3) * 0.8);
    end
end

function I = debye_intensity(q, elements, xyz)
    q = q(:);
    nQ = numel(q);
    nAtoms = size(xyz, 1);
    f = zeros(nQ, nAtoms);
    for i = 1:nAtoms
        f(:, i) = atomic_form_factor(q, elements(i));
    end
    rij = zeros(nAtoms, nAtoms);
    for i = 1:nAtoms
        for j = 1:nAtoms
            rij(i, j) = norm(xyz(i, :) - xyz(j, :));
        end
    end
    I = zeros(nQ, 1);
    for iq = 1:nQ
        val = 0;
        for i = 1:nAtoms
            for j = 1:nAtoms
                qr = q(iq) * rij(i, j);
                if rij(i, j) < 1e-12
                    sinc = 1;
                else
                    sinc = sin(qr) / qr;
                end
                val = val + f(iq, i) * f(iq, j) * sinc;
            end
        end
        I(iq) = val;
    end
end

function f = atomic_form_factor(q, element)
    s2 = (q(:) ./ (4 * pi)).^2;
    switch upper(string(element))
        case "C"
            a = [2.3100, 1.0200, 1.5886, 0.8650]; b = [20.8439, 10.2075, 0.5687, 51.6512]; c = 0.2156;
        case "N"
            a = [12.2126, 3.1322, 2.0125, 1.1663]; b = [0.0057, 9.8933, 28.9975, 0.5826]; c = -11.5290;
        case "O"
            a = [3.0485, 2.2868, 1.5463, 0.8670]; b = [13.2771, 5.7011, 0.3239, 32.9089]; c = 0.2508;
        case "S"
            a = [6.9053, 5.2034, 1.4379, 1.5863]; b = [1.4679, 22.2151, 0.2536, 56.1720]; c = 0.8669;
        otherwise
            error("Unknown element %s", element);
    end
    f = c + zeros(size(q(:)));
    for k = 1:4
        f = f + a(k) .* exp(-b(k) .* s2);
    end
end

function plot_case_results(q, t, delay_s, data, fit, eval, kind, cfg)
    k = char(kind);
    figure("Name", k + "_01_true_and_recovered_concentrations");
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    plot(t, data.C_true_t, "LineWidth", 1.5);
    set(gca, "XScale", "log"); grid on; ylim([-0.02, 1.05]);
    xlabel("Time (s)"); ylabel("Concentration");
    title(k + " true concentrations");
    legend(["species 1", "species 2", "species 3"], "Location", "best");
    nexttile;
    plot(t, eval.C_est_matched, "--", "LineWidth", 1.5);
    set(gca, "XScale", "log"); grid on; ylim([-0.02, 1.05]);
    xlabel("Time (s)"); ylabel("Concentration");
    title(sprintf("%s recovered, mean corr %.3f", k, mean(eval.concCorr)));
    legend(["matched 1", "matched 2", "matched 3"], "Location", "best");

    figure("Name", k + "_02_debye_spectra_recovery");
    tiledlayout(1, 3, "TileSpacing", "compact");
    for i = 1:3
        nexttile;
        st = data.S_true(:, i) ./ max(norm(data.S_true(:, i)), eps);
        se = eval.S_est_matched(:, i) ./ max(norm(eval.S_est_matched(:, i)), eps);
        if corr_safe(st, se) < 0
            se = -se;
        end
        plot(q, st, "k-", "LineWidth", 1.5); hold on;
        plot(q, se, "--", "LineWidth", 1.3); hold off;
        grid on; xlabel("q (A^{-1})"); ylabel("normalized Delta S(q)");
        title(sprintf("species %d corr %.3f", i, eval.specCorr(i)));
    end

    figure("Name", k + "_03_zero_lag_reconstruction");
    [~, zidx] = min(abs(delay_s));
    tiledlayout(1, 3, "TileSpacing", "compact");
    clim = max(abs(data.X(:, :, zidx)), [], "all");
    nexttile;
    imagesc(t, q, data.X(:, :, zidx)); set(gca, "XScale", "log", "YDir", "normal"); colorbar; caxis([-clim clim]);
    xlabel("Time (s)"); ylabel("q"); title(k + " noisy target");
    nexttile;
    imagesc(t, q, fit.qt_recon); set(gca, "XScale", "log", "YDir", "normal"); colorbar; caxis([-clim clim]);
    xlabel("Time (s)"); ylabel("q"); title(sprintf("fit zero err %.3g", fit.zeroRelerr));
    nexttile;
    imagesc(t, q, data.X(:, :, zidx) - fit.qt_recon); set(gca, "XScale", "log", "YDir", "normal"); colorbar; caxis([-0.15*clim 0.15*clim]);
    xlabel("Time (s)"); ylabel("q"); title("residual");

    if kind == "branch"
        figure("Name", k + "_04_open_generator_rates");
        bar([data.truth.rates(:), fit.ratesPruned(:)]);
        set(gca, "XTick", 1:6, "XTickLabel", ["1->2","1->3","1->G","2->3","2->G","3->G"]);
        xtickangle(30); grid on; ylabel("rate (s^{-1})");
        legend(["truth", "estimated/pruned"], "Location", "best");
        title(sprintf("branch edge F1 %.3f, f true %.2f/%.2f, f est %.2f/%.2f", ...
            eval.edgeF1, data.truth.branchFraction(1), data.truth.branchFraction(2), ...
            eval.branchFractionEst(1), eval.branchFractionEst(2)));
    else
        figure("Name", k + "_04_parallel_rates_and_weights");
        tiledlayout(1, 2, "TileSpacing", "compact");
        nexttile;
        bar([eval.rateTrueMatched(:), eval.rateEstMatched(:)]);
        set(gca, "YScale", "log"); grid on; ylabel("rate (s^{-1})");
        legend(["truth", "estimated"], "Location", "best");
        title("matched decay rates");
        nexttile;
        bar([eval.weightTrueMatched(:), eval.weightEstMatched(:)]);
        ylim([0 0.7]); grid on; ylabel("initial/channel weight");
        legend(["truth", "estimated"], "Location", "best");
        title("matched weights");
    end
end

function write_summary(cfg, branchData, branchFit, branchEval, parallelData, parallelFit, parallelEval)
    fid = fopen(fullfile(cfg.outDir, "summary.txt"), "w");
    if fid < 0
        return;
    end
    fprintf(fid, "Debye equation with atomic form factors validation\n");
    fprintf(fid, "noiseRel %.6g\n", cfg.noiseRel);
    fprintf(fid, "\n[branch]\n");
    fprintf(fid, "zero_lag_relerr %.12g\n", branchFit.zeroRelerr);
    fprintf(fid, "full_tensor_relerr %.12g\n", branchFit.fullRelerr);
    fprintf(fid, "weighted_tensor_relerr %.12g\n", branchFit.weightedRelerr);
    fprintf(fid, "spec_corr %s\n", mat2str(branchEval.specCorr, 8));
    fprintf(fid, "conc_corr %s\n", mat2str(branchEval.concCorr, 8));
    fprintf(fid, "edge_precision %.12g\n", branchEval.edgePrecision);
    fprintf(fid, "edge_recall %.12g\n", branchEval.edgeRecall);
    fprintf(fid, "edge_F1 %.12g\n", branchEval.edgeF1);
    fprintf(fid, "rates_true %s\n", mat2str(branchData.truth.rates, 12));
    fprintf(fid, "rates_est_pruned %s\n", mat2str(branchFit.ratesPruned, 12));
    fprintf(fid, "branch_fraction_true %s\n", mat2str(branchData.truth.branchFraction, 12));
    fprintf(fid, "branch_fraction_est %s\n", mat2str(branchEval.branchFractionEst, 12));
    fprintf(fid, "\n[parallel]\n");
    fprintf(fid, "zero_lag_relerr %.12g\n", parallelFit.zeroRelerr);
    fprintf(fid, "full_tensor_relerr %.12g\n", parallelFit.fullRelerr);
    fprintf(fid, "weighted_tensor_relerr %.12g\n", parallelFit.weightedRelerr);
    fprintf(fid, "spec_corr %s\n", mat2str(parallelEval.specCorr, 8));
    fprintf(fid, "conc_corr %s\n", mat2str(parallelEval.concCorr, 8));
    fprintf(fid, "rates_true %s\n", mat2str(parallelData.truth.rates, 12));
    fprintf(fid, "rates_est_raw %s\n", mat2str(parallelFit.rates, 12));
    fprintf(fid, "weights_true %s\n", mat2str(parallelData.truth.weights, 12));
    fprintf(fid, "weights_est_raw %s\n", mat2str(parallelFit.weights, 12));
    fprintf(fid, "rates_est_matched %s\n", mat2str(parallelEval.rateEstMatched, 12));
    fprintf(fid, "weights_est_matched %s\n", mat2str(parallelEval.weightEstMatched, 12));
    fclose(fid);
end
