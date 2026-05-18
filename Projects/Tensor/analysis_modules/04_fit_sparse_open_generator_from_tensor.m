%% 04_fit_sparse_open_generator_from_tensor.m
% Fit one sparse first-order acyclic generator inside the tensor-derived
% kinetic manifold.  This is not a branch-vs-sequential model sweep:
% all peak-ordered active/ground edges are present at once, and the
% unnecessary rates are suppressed by fit quality plus sparsity.

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
cfg.resultFile = fullfile(pwd, "tensor_mode_closure_rotation_results", "workspace_result.mat");
cfg.outDir = fullfile(pwd, "sparse_open_generator_results");
cfg.randomSeed = 101;
cfg.numStarts = 44;
cfg.maxIter = 650;
cfg.zeroLagWeight = 7.0;
cfg.shortDelayWeight = 3.0;
cfg.shortDelayScale_s = 3e-9;
cfg.zeroLagSurfaceWeight = 120;
cfg.sparsityWeight = 0.018;
cfg.initialGroundPenalty = 0.80;
cfg.earlyGroundWeight = 0.15;
cfg.peakOrderPenalty = 0.45;
cfg.minPeakGapLog10 = 1.0;
cfg.minPeakHeight = 0.25;
cfg.peakHeightPenalty = 0.20;
cfg.rateLower_s_inv = 1e3;
cfg.rateUpper_s_inv = 1e14;
cfg.pruneRelativeToSource = 2e-3;
cfg.useIrfConvolution = true;

if exist("cfg_override", "var") && isstruct(cfg_override)
    overrideNames = fieldnames(cfg_override);
    for iOverride = 1:numel(overrideNames)
        cfg.(overrideNames{iOverride}) = cfg_override.(overrideNames{iOverride});
    end
end

scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) > 0
    cd(scriptDir);
end

R = load(cfg.resultFile);
q = R.q(:);
t = R.t(:);
delay_s = R.delay_s(:);
X = R.X;
best_par_struct = R.best_par_struct;

[nQ, nT, nD] = size(X);
[~, zeroLagIndex] = min(abs(delay_s));
uSurface = surface_shifted_time_axis(t, delay_s, R.cfg.shiftMode);
wSurface = make_surface_fit_weights(nT, nD, zeroLagIndex, delay_s, cfg);
Xmat = reshape(X, nQ, nT*nD);
normXw = sqrt(sum(sum((Xmat.^2) .* wSurface.')));
qt_target = X(:, :, zeroLagIndex);
normD0 = norm(qt_target(:));

tau0 = sort(mean([R.tau_time(:), R.tau_delay(:)], 2), "ascend");
rate0 = 1 ./ tau0(:).';
fprintf("Tensor lifetime seed tau(s): %s\n", mat2str(tau0(:).', 8));

[thetaBest, modelBest] = fit_sparse_generator_varpro( ...
    X, Xmat, wSurface, normXw, qt_target, normD0, t, delay_s, uSurface, ...
    zeroLagIndex, rate0, best_par_struct, cfg);

if ~exist(cfg.outDir, "dir")
    mkdir(cfg.outDir);
end

model = modelBest; %#ok<NASGU>
theta = thetaBest; %#ok<NASGU>
save(fullfile(cfg.outDir, "workspace_sparse_generator_result.mat"), ...
    "cfg", "q", "t", "delay_s", "model", "theta", "best_par_struct", "-v7.3");

plot_sparse_generator_results(q, t, qt_target, modelBest, cfg);
figs = findall(0, "Type", "figure");
for i = 1:numel(figs)
    fig = figs(i);
    safe = regexprep(string(fig.Name), "[^\w가-힣-]", "_");
    exportgraphics(fig, fullfile(cfg.outDir, sprintf("%02d_%s.png", i, safe)), "Resolution", 170);
end

disp("Sparse generator rate table:");
disp(modelBest.rate_table);
fprintf("SPARSE_METRIC zero=%.12g tensor=%.12g weighted=%.12g overlap=%.12g earlyG=%.12g endG=%.12g score=%.12g\n", ...
    modelBest.zero_lag_relerr, modelBest.full_tensor_relerr, modelBest.weighted_tensor_relerr, ...
    modelBest.overlap, modelBest.early_ground_mean, modelBest.C_ground_at_t(end), modelBest.score);
fprintf("SPARSE_METRIC peaks=%s peakTimes=%s groundAtPeaks=%s\n", ...
    mat2str(max(modelBest.C_at_t, [], 1), 6), ...
    mat2str(modelBest.peak_time_s(:).', 6), ...
    mat2str(modelBest.C_ground_at_t(modelBest.peak_index).', 6));
fprintf("OUTDIR=%s\n", cfg.outDir);

function [thetaBest, best] = fit_sparse_generator_varpro(X, Xmat, wSurface, normXw, ...
        qt_target, normD0, t, delay_s, uSurface, zeroLagIndex, rate0, best_par_struct, cfg)
    rng(cfg.randomSeed);
    starts = make_generator_starts(rate0, cfg);
    opts = optimset("Display", "off", "MaxIter", cfg.maxIter, ...
        "MaxFunEvals", cfg.maxIter * 35, "TolX", 1e-7, "TolFun", 1e-8);

    bestScore = inf;
    thetaBest = starts(1, :);
    best = struct();
    for iStart = 1:size(starts, 1)
        obj = @(theta) sparse_generator_objective(theta, X, Xmat, wSurface, normXw, ...
            qt_target, normD0, t, delay_s, uSurface, zeroLagIndex, best_par_struct, cfg);
        [thetaOpt, score] = fminsearch(obj, starts(iStart, :), opts);
        [score, model] = sparse_generator_objective(thetaOpt, X, Xmat, wSurface, normXw, ...
            qt_target, normD0, t, delay_s, uSurface, zeroLagIndex, best_par_struct, cfg);
        if score < bestScore
            bestScore = score;
            thetaBest = thetaOpt;
            best = model;
            best.best_start = iStart;
            fprintf("start %02d improved: score %.8g, zero %.6g, tensor %.6g, rates %s\n", ...
                iStart, score, model.zero_lag_relerr, model.weighted_tensor_relerr, ...
                mat2str(model.rates_s_inv, 4));
        end
    end
end

function starts = make_generator_starts(rate0, cfg)
    lo = log10(cfg.rateLower_s_inv);
    hi = log10(cfg.rateUpper_s_inv);
    tiny = max(cfg.rateLower_s_inv * 10, min(rate0) * 1e-4);

    seq = [rate0(1), tiny, tiny, rate0(2), tiny, rate0(3)];
    branch = [0.55*rate0(1), 0.45*rate0(1), tiny, tiny, rate0(2), rate0(3)];
    mixed = [0.65*rate0(1), 0.12*rate0(1), tiny, rate0(2), tiny, rate0(3)];
    cascadeFastGround = [rate0(1), tiny, tiny, rate0(2), 0.05*rate0(2), rate0(3)];

    base = [seq; branch; mixed; cascadeFastGround];
    starts = zeros(cfg.numStarts, 6);
    for i = 1:min(size(base, 1), cfg.numStarts)
        starts(i, :) = log10(min(max(base(i, :), cfg.rateLower_s_inv), cfg.rateUpper_s_inv));
    end
    for i = size(base, 1)+1:cfg.numStarts
        jitter = 0.8 * randn(1, 6);
        if mod(i, 2) == 0
            seed = log10(seq);
        else
            seed = log10(branch);
        end
        starts(i, :) = min(max(seed + jitter, lo), hi);
    end
end

function [score, model] = sparse_generator_objective(theta, X, Xmat, wSurface, normXw, ...
        qt_target, normD0, t, delay_s, uSurface, zeroLagIndex, best_par_struct, cfg)
    theta = theta(:).';
    lo = log10(cfg.rateLower_s_inv);
    hi = log10(cfg.rateUpper_s_inv);
    thetaClip = min(max(theta, lo), hi);
    boundPenalty = sum((theta - thetaClip).^2);
    rates = 10 .^ thetaClip;
    K = generator_matrix_from_rates(rates);

    CfullSurface = generator_concentrations(uSurface, K, best_par_struct, cfg);
    Csurface = max(CfullSurface(:, 1:3), 0);
    [S_q, XhatMat] = solve_weighted_spectra(Xmat, Csurface, wSurface);
    Xhat = reshape(XhatMat, size(X));

    zeroRows = (1:numel(t)) + (zeroLagIndex - 1) * numel(t);
    C_at_t = Csurface(zeroRows, :);
    C_ground_at_t = max(CfullSurface(zeroRows, 4), 0);
    qt_recon = Xhat(:, :, zeroLagIndex);

    weightedRelerr = sqrt(sum(sum(((Xmat - XhatMat).^2) .* wSurface.'))) / max(normXw, eps);
    zeroRelerr = norm(qt_target(:) - qt_recon(:)) / max(normD0, eps);
    fullRelerr = norm(X(:) - Xhat(:)) / max(norm(X(:)), eps);

    [peakVal, peakIdx] = max(C_at_t, [], 1);
    peakTime = t(peakIdx);
    peakOrderPenalty = peak_order_penalty(peakTime, cfg);
    peakHeightPenalty = mean(max(cfg.minPeakHeight - peakVal, 0).^2);
    overlap = concentration_overlap_metric(C_at_t);

    sourceOut = [sum(rates(1:3)), sum(rates(4:5)), rates(6)];
    nonterminalGroundFrac = (rates(3) + rates(5)) / max(sum(rates), eps);
    initialGroundFrac = rates(3) / max(sourceOut(1), eps);
    sparsityPenalty = sum(log1p(rates ./ max(sourceOut([1, 1, 1, 2, 2, 3]), eps))) / numel(rates);
    earlyGroundMean = mean(C_ground_at_t(t <= max(peakTime)));

    score = weightedRelerr + cfg.zeroLagWeight * zeroRelerr ...
        + cfg.sparsityWeight * sparsityPenalty ...
        + cfg.initialGroundPenalty * initialGroundFrac.^2 ...
        + cfg.earlyGroundWeight * earlyGroundMean.^2 ...
        + cfg.peakOrderPenalty * peakOrderPenalty ...
        + cfg.peakHeightPenalty * peakHeightPenalty ...
        + 100 * boundPenalty;

    model = struct();
    model.score = score;
    model.K = K;
    model.rates_s_inv = rates;
    model.rate_names = ["1->2", "1->3", "1->ground", "2->3", "2->ground", "3->ground"];
    model.S_q = S_q;
    model.C_at_t = C_at_t;
    model.C_ground_at_t = C_ground_at_t;
    model.C_surface = Csurface;
    model.C_ground_surface = max(CfullSurface(:, 4), 0);
    model.qt_recon = qt_recon;
    model.Xhat = Xhat;
    model.weighted_tensor_relerr = weightedRelerr;
    model.full_tensor_relerr = fullRelerr;
    model.zero_lag_relerr = zeroRelerr;
    model.overlap = overlap;
    model.peak_index = peakIdx;
    model.peak_time_s = peakTime;
    model.peak_height = peakVal;
    model.initial_ground_fraction = initialGroundFrac;
    model.nonterminal_ground_fraction = nonterminalGroundFrac;
    model.early_ground_mean = earlyGroundMean;
    model.rate_table = make_rate_table(model.rate_names, rates, cfg);
end

function K = generator_matrix_from_rates(r)
    k12 = r(1); k13 = r(2); k1g = r(3);
    k23 = r(4); k2g = r(5); k3g = r(6);
    K = [-(k12+k13+k1g), 0, 0, 0; ...
          k12, -(k23+k2g), 0, 0; ...
          k13, k23, -k3g, 0; ...
          k1g, k2g, k3g, 0];
end

function C = generator_concentrations(u, K, par, cfg)
    u = u(:);
    c0 = [1; 0; 0; 0];
    [V, L] = eig(K);
    coeff = V \ c0;
    C = zeros(numel(u), 4);
    for m = 1:4
        lambda = real(L(m, m));
        if cfg.useIrfConvolution && lambda < -eps
            basis = conv_gauss_expdec_local(u, 1, 0, -1/lambda, par.t0, par.FWHM);
        elseif cfg.useIrfConvolution
            basis = 0.5 * erfc(-(u - par.t0) ./ (sqrt(2) * (par.FWHM / 2.3548)));
        else
            basis = double(u >= 0) .* exp(lambda * max(u, 0));
        end
        amp = real(V(:, m) * coeff(m));
        C = C + basis * amp.';
    end
    C = max(real(C), 0);
    rowSum = sum(C, 2);
    over = rowSum > 1 + 1e-8;
    C(over, :) = C(over, :) ./ rowSum(over);
end

function [S_q, XhatMat] = solve_weighted_spectra(Xmat, Csurface, wSurface)
    G0 = Csurface.' * (Csurface .* wSurface);
    ridge = 1e-9 * max(trace(G0) / max(size(Csurface, 2), 1), eps);
    S_q = ((Xmat .* wSurface.') * Csurface) / (G0 + ridge * eye(size(Csurface, 2)));
    XhatMat = S_q * Csurface.';
end

function penalty = peak_order_penalty(peakTime, cfg)
    x = log10(max(peakTime(:), realmin));
    reversePenalty = mean(max(-diff(x), 0).^2);
    gapPenalty = mean(max(cfg.minPeakGapLog10 - diff(x), 0).^2);
    penalty = reversePenalty + gapPenalty;
end

function T = make_rate_table(names, rates, cfg)
    pruned = rates(:);
    source = [1; 1; 1; 2; 2; 3];
    for src = unique(source).'
        idx = source == src;
        maxRate = max(pruned(idx));
        if maxRate > 0
            pruned(idx & pruned < cfg.pruneRelativeToSource * maxRate) = 0;
        end
    end
    [sortedRates, order] = sort(pruned, "descend");
    T = table(cellstr(names(order).'), sortedRates, ...
        'VariableNames', {'transition', 'rate_s_inv'});
end

function wSurface = make_surface_fit_weights(nT, nD, zeroLagIndex, delay_s, cfg)
    wDelay = ones(nD, 1);
    if cfg.shortDelayWeight > 0
        wDelay = wDelay + cfg.shortDelayWeight ./ (1 + delay_s(:) ./ cfg.shortDelayScale_s);
    end
    if cfg.zeroLagSurfaceWeight > 1
        wDelay(zeroLagIndex) = wDelay(zeroLagIndex) * cfg.zeroLagSurfaceWeight;
    end
    wSurface = repelem(wDelay, nT);
end

function uSurface = surface_shifted_time_axis(t, delay_s, shiftMode)
    nT = numel(t);
    nD = numel(delay_s);
    uSurface = zeros(nT*nD, 1);
    for k = 1:nD
        idx = (1:nT) + (k - 1) * nT;
        switch string(shiftMode)
            case "evaluate_later_time"
                uSurface(idx) = t(:) + delay_s(k);
            case "delay_signal"
                uSurface(idx) = t(:) - delay_s(k);
            otherwise
                error("Unknown shift mode: %s", shiftMode);
        end
    end
end

function overlap = concentration_overlap_metric(C)
    C = max(C, 0);
    rowSum = sum(C, 2);
    valid = rowSum > 1e-12;
    if ~any(valid)
        overlap = inf;
        return;
    end
    P = C(valid, :) ./ rowSum(valid);
    entropy = -sum(P .* log(P + eps), 2);
    overlap = mean(entropy) / log(size(C, 2));
end

function y = conv_gauss_expdec_local(t, A, A0, tau, t0, FWHM)
    t = t(:);
    sigma = FWHM / 2.3548;
    y = A0 * ones(numel(t), 1);
    inv_sqrt2_sigma = 1 / (sqrt(2) * sigma);
    gauss_env = exp(-((t - t0).^2) ./ (2 * sigma^2));
    Z = ((sigma^2)/tau - (t - t0)) * inv_sqrt2_sigma;
    comp = zeros(numel(t), 1);
    idx_pos = (Z >= 0);
    if any(idx_pos)
        comp(idx_pos) = gauss_env(idx_pos) .* erfcx(Z(idx_pos));
    end
    idx_neg = (Z < 0);
    if any(idx_neg)
        E = (sigma^2)/(2*tau^2) - (t(idx_neg) - t0)/tau;
        comp(idx_neg) = exp(E) .* erfc(Z(idx_neg));
    end
    y = y + 0.5 * A .* comp;
end

function plot_sparse_generator_results(q, t, qt_target, model, cfg)
    figure("Name", "01_sparse_generator_concentrations");
    plot(t, [model.C_at_t, model.C_ground_at_t], "LineWidth", 1.5);
    set(gca, "XScale", "log");
    grid on;
    xlabel("Time (s)");
    ylabel("Concentration");
    title("Sparse first-order generator concentrations");
    legend(["component 1", "component 2", "component 3", "silent ground"], "Location", "best");
    ylim([-0.02, 1.05]);

    figure("Name", "02_sparse_generator_q_spectra");
    plot(q, model.S_q, "LineWidth", 1.4);
    hold on; yline(0, "k--"); hold off;
    grid on;
    xlabel("q");
    ylabel("Species spectrum (a.u.)");
    title("Sparse generator active species spectra");
    legend(["component 1", "component 2", "component 3"], "Location", "best");

    figure("Name", "03_sparse_generator_zero_lag_reconstruction");
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    imagesc(t, q, qt_target);
    set(gca, "XScale", "log", "YDir", "normal");
    xlabel("Time (s)"); ylabel("q"); title("Zero-lag target"); colorbar;
    nexttile;
    imagesc(t, q, model.qt_recon);
    set(gca, "XScale", "log", "YDir", "normal");
    xlabel("Time (s)"); ylabel("q");
    title(sprintf("Sparse generator fit %.4f", 1 - model.zero_lag_relerr)); colorbar;

    figure("Name", "04_sparse_generator_rates");
    bar(model.rate_table.rate_s_inv);
    set(gca, "XTick", 1:height(model.rate_table), "XTickLabel", model.rate_table.transition);
    xtickangle(35);
    grid on;
    ylabel("rate (s^{-1})");
    title("Sparse acyclic generator rates");

    fid = fopen(fullfile(cfg.outDir, "summary.txt"), "w");
    if fid >= 0
        fprintf(fid, "zero_lag_relerr %.12g\n", model.zero_lag_relerr);
        fprintf(fid, "weighted_tensor_relerr %.12g\n", model.weighted_tensor_relerr);
        fprintf(fid, "full_tensor_relerr %.12g\n", model.full_tensor_relerr);
        fprintf(fid, "overlap %.12g\n", model.overlap);
        fprintf(fid, "rates %s\n", mat2str(model.rates_s_inv, 12));
        fclose(fid);
    end
end
