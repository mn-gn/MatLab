%% unique_tensor_modes_closure_rotation.m
% Shifted-time tensor를 먼저 CP/PARAFAC로 분해해 lifetime eigenmode를 얻고,
% 그 eigenmode 공간 안에서 nonnegative closure rotation으로 species profile을 찾는다.
%
% 핵심 흐름:
%   1) D(q,t)를 SVD-RSV fitting으로 smooth kinetic manifold로 만든다.
%   2) X(q,t,delay)=D_model(q,t+delay) 텐서를 만든다.
%   3) CP-ALS로 X ~= sum_r A_r(q) B_r(t) C_r(delay)를 분해한다.
%   4) B_r(t)C_r(delay)가 unique tensor eigenmode라고 보고,
%      그 span 안에서 c_active>=0, sum(c_active)<=1인 회전을 찾는다.
%   5) c_ground=1-sum(c_active)는 silent ground로만 계산하고 spectrum은 fit하지 않는다.

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

%% ===================== 설정 =====================
cfg = struct();
cfg.qFile    = "q.csv";
cfg.tFile    = "t_81points.csv";
cfg.dataFile = "q_delta_S_noise_10%_81points.csv";

cfg.numPC = 3;
cfg.numTau = 3;
cfg.numStart = 100;
cfg.fixSVDSign = true;
cfg.centerData = false;

cfg.useLogTimeGridForTensor = true;
cfg.numTensorTimePoints = [];
cfg.shiftMode = "evaluate_later_time";  % t + delay

delay_min_s = 1e-14;
delay_max_s = 100e-9;
num_log_delay = 100;
cfg.delay_s = [0, logspace(log10(delay_min_s), log10(delay_max_s), num_log_delay)];

cfg.cpRank = 3;
cfg.cpMaxIters = 300;
cfg.cpTol = 1e-8;
cfg.cpRandomSeed = 7;
cfg.cpVerbose = false;

cfg.rotationNumStarts = 180;
cfg.rotationMaxIter = 700;
cfg.rotationRandomSeed = 43;
cfg.rotationOverlapWeight = 0.25;
cfg.rotationNegativePenalty = 1e2;
cfg.rotationClosurePenalty = 1e2;
cfg.rotationConditionPenalty = 1e-3;
cfg.rotationZeroLagWeight = 35.0;
cfg.rotationZeroLagSolveWeight = 120;
cfg.rotationShortDelaySolveWeight = 3;
cfg.rotationShortDelayScale_s = 3e-9;
cfg.rotationPeakPenalty = 20;
cfg.rotationMinComponentPeak = 0.20;
cfg.rotationPeakTarget = 0.90;
cfg.rotationPurityPenalty = 12;
cfg.rotationMinPeakPurity = 0.58;
cfg.rotationEndpointPenalty = 10;
cfg.rotationActiveSumMonotonicPenalty = 18;
cfg.rotationEarlyGroundPenalty = 200;
cfg.rotationEarlyGroundMaxPenalty = 500;
cfg.rotationPeakGroundPenalty = 300;
cfg.rotationPeakSpreadPenalty = 35;
cfg.rotationMinPeakLogGap = 1.2;
cfg.rotationMinLastPeakLogFraction = 0.55;
cfg.rotationInitialStatePenalty = 2;
cfg.applyTerminalGroundProjection = true;
cfg.solveSpectraOnProjectedSurface = true;
cfg.rotationGeneratorPenalty = 0.20;
cfg.rotationGeneratorTimeWeightAlpha = 0.25;
cfg.rotationNonterminalGroundRatePenalty = 0.05;
cfg.topologyUsePeakOrderedAcyclic = true;
cfg.topologyTimeWeightAlpha = 0.5;
cfg.topologyMinRelativeSourceRate = 1e-3;

cfg.outDir = fullfile(pwd, "unique_tensor_modes_closure_rotation_results");

if exist("cfg_override", "var") && isstruct(cfg_override)
    overrideNames = fieldnames(cfg_override);
    for iOverride = 1:numel(overrideNames)
        cfg.(overrideNames{iOverride}) = cfg_override.(overrideNames{iOverride});
    end
end

%% ===================== 경로 / helper =====================
scriptDir = fileparts(mfilename("fullpath"));
if strlength(scriptDir) > 0
    cd(scriptDir);
end
if exist("fit_rsv", "file") ~= 2
    error("fit_rsv.m 파일이 현재 폴더 또는 MATLAB path에 없습니다.");
end
ensure_asinh_helpers_exist();

%% ===================== 데이터 / SVD / RSV fitting =====================
q = read_numeric_table(cfg.qFile);
t_raw = read_numeric_table(cfg.tFile);
D = read_numeric_table(cfg.dataFile);
q = q(:);
t_raw = t_raw(:);

if size(D, 1) == numel(q) && size(D, 2) == numel(t_raw)
    % OK
elseif size(D, 1) == numel(t_raw) && size(D, 2) == numel(q)
    D = D.';
else
    error("데이터 크기가 q/t와 맞지 않습니다.");
end

if cfg.centerData
    D_mean_q = mean(D, 2);
    D_work = D - D_mean_q;
else
    D_mean_q = zeros(size(D, 1), 1);
    D_work = D;
end

[U, S, V] = svd(D_work, "econ");
numPC = min([cfg.numPC, size(U, 2), size(V, 2)]);
U_pc = U(:, 1:numPC);
S_pc = S(1:numPC, 1:numPC);
V_pc = V(:, 1:numPC);

if cfg.fixSVDSign
    for r = 1:numPC
        [~, idxMax] = max(abs(V_pc(:, r)));
        if V_pc(idxMax, r) < 0
            U_pc(:, r) = -U_pc(:, r);
            V_pc(:, r) = -V_pc(:, r);
        end
    end
end

[best_par_struct, ~, ~] = fit_rsv(t_raw, V_pc, cfg.numTau, cfg.numStart);
close all;

if cfg.useLogTimeGridForTensor
    t_pos = t_raw(t_raw > 0);
    nTensorTime = cfg.numTensorTimePoints;
    if isempty(nTensorTime)
        nTensorTime = numel(t_raw);
    end
    t = logspace(log10(min(t_pos)), log10(max(t_pos)), nTensorTime).';
else
    t = t_raw;
end

fprintf("RSV fitted tau(s): ");
fprintf("%.6g ", sort(best_par_struct.tau(:)));
fprintf("\n");

%% ===================== shifted-time tensor 생성 =====================
nQ = numel(q);
nT = numel(t);
nDelay = numel(cfg.delay_s);
X = zeros(nQ, nT, nDelay);

for k = 1:nDelay
    t_eval = shifted_time_local(t, cfg.delay_s(k), cfg);
    V_shift = eval_fitted_rsv(best_par_struct, t_eval);
    X(:, :, k) = U_pc * S_pc * V_shift.' + D_mean_q;
end

fprintf("Shifted tensor: %d(q) x %d(time) x %d(delay)\n", size(X));

%% ===================== CP tensor decomposition =====================
cp = cp_als_signed_q_nonneg_td(X, cfg.cpRank, cfg);
[~, zeroLagIndex] = min(abs(cfg.delay_s));

% CP eigenmode surface Z_eig(:,r)=B_r(t)C_r(delay)
A_eig = cp.A_q .* reshape(cp.lambda, 1, []);
B_eig = cp.B_time;
C_eig = cp.C_delay;
Z_eig = mode_surface_matrix(B_eig, C_eig);
E0_eig = B_eig .* C_eig(zeroLagIndex, :);
qt_target = X(:, :, zeroLagIndex);
qt_cp = A_eig * E0_eig.';

tau_time = estimate_factor_tau(t, B_eig);
tau_delay = estimate_factor_tau(cfg.delay_s(:), C_eig);
fprintf("CP tensor relerr: %.6g, fit %.6f\n", cp.relerr, cp.fit);
fprintf("CP time-factor tau estimates(s): ");
fprintf("%.6g ", tau_time);
fprintf("\n");
fprintf("CP delay-factor tau estimates(s): ");
fprintf("%.6g ", tau_delay);
fprintf("\n");

%% ===================== closure-constrained species rotation =====================
cfg.tensor_time_s = t;
rot = closure_rotate_tensor_modes(X, Z_eig, E0_eig, cfg);
qt_rot = rot.S_q * rot.C_at_t.';

fprintf("Closure rotation tensor relerr: %.6g, zero-lag relerr: %.6g\n", ...
    rot.tensor_relerr, rot.zero_lag_relerr);
fprintf("Closure max error: %.6g\n", rot.closure_max_error);

topology = infer_sparse_active_ground_topology(t, rot.C_at_t, cfg);
disp("Inferred active-state transition rates from rotated concentration profiles:");
disp(topology.rate_table);

%% ===================== 결과 저장 / 그림 =====================
if ~exist(cfg.outDir, "dir")
    mkdir(cfg.outDir);
end

result = struct();
result.cfg = cfg;
result.q = q;
result.t = t;
result.delay_s = cfg.delay_s;
result.X = X;
result.best_par_struct = best_par_struct;
result.cp = cp;
result.A_eig = A_eig;
result.B_eig = B_eig;
result.C_eig = C_eig;
result.tau_time = tau_time;
result.tau_delay = tau_delay;
result.rotation = rot;
result.topology = topology;
save(fullfile(cfg.outDir, "workspace_result.mat"), "-struct", "result", "-v7.3");

plot_results(q, t, cfg.delay_s, X, qt_target, qt_cp, qt_rot, cp, A_eig, B_eig, C_eig, rot, topology);

figs = findall(0, "Type", "figure");
for i = 1:numel(figs)
    figure(figs(i));
    drawnow;
    name = get(figs(i), "Name");
    if strlength(string(name)) == 0
        name = sprintf("figure_%02d", i);
    end
    safe = regexprep(char(name), "[^a-zA-Z0-9_.-]", "_");
    exportgraphics(figs(i), fullfile(cfg.outDir, sprintf("%02d_%s.png", i, safe)), "Resolution", 160);
end

fprintf("OUTDIR=%s\n", cfg.outDir);

%% ===================== Local functions =====================
function cp = cp_als_signed_q_nonneg_td(X, R, cfg)
    rng(cfg.cpRandomSeed);
    [nQ, nT, nD] = size(X);
    Xnorm = norm(X(:));

    X1 = reshape(X, nQ, nT*nD);
    X2 = reshape(permute(X, [2 1 3]), nT, nQ*nD);
    X3 = reshape(permute(X, [3 1 2]), nD, nQ*nT);

    [A, B, C] = initialize_cp_factors(X1, X2, X3, R);
    fitHistory = zeros(cfg.cpMaxIters, 1);
    prevRelerr = inf;

    for iter = 1:cfg.cpMaxIters
        K = khatri_rao_local(C, B);
        A = X1 * K / (K.' * K + 1e-10 * eye(R));

        K = khatri_rao_local(C, A);
        B = max(X2 * K / (K.' * K + 1e-10 * eye(R)), 0);
        [A, B] = absorb_norm_into_left(A, B);

        K = khatri_rao_local(B, A);
        C = max(X3 * K / (K.' * K + 1e-10 * eye(R)), 0);
        [A, C] = absorb_norm_into_left(A, C);

        Xhat = reconstruct_cp_local(A, B, C);
        relerr = norm(X(:) - Xhat(:)) / max(Xnorm, eps);
        fitHistory(iter) = relerr;

        if cfg.cpVerbose && (iter == 1 || mod(iter, 25) == 0)
            fprintf("CP iter %3d: relerr %.6g\n", iter, relerr);
        end
        if abs(prevRelerr - relerr) < cfg.cpTol * max(prevRelerr, 1)
            fitHistory = fitHistory(1:iter);
            break;
        end
        prevRelerr = relerr;
    end

    lambda = vecnorm(A, 2, 1);
    for r = 1:R
        if lambda(r) > 0
            A(:, r) = A(:, r) / lambda(r);
        end
    end

    [~, order] = sort(estimate_factor_tau((1:nT).', B), "ascend");
    A = A(:, order);
    B = B(:, order);
    C = C(:, order);
    lambda = lambda(order);

    Xhat = reconstruct_cp_local(A .* reshape(lambda, 1, []), B, C);
    relerr = norm(X(:) - Xhat(:)) / max(Xnorm, eps);

    cp = struct();
    cp.A_q = A;
    cp.B_time = B;
    cp.C_delay = C;
    cp.lambda = lambda(:).';
    cp.relerr = relerr;
    cp.fit = 1 - relerr;
    cp.num_iters = numel(fitHistory);
    cp.fit_history = fitHistory;
    cp.Xhat = Xhat;
end

function [A, B, C] = initialize_cp_factors(X1, X2, X3, R)
    [U1, ~, ~] = svd(X1, "econ");
    [U2, ~, ~] = svd(X2, "econ");
    [U3, ~, ~] = svd(X3, "econ");
    A = U1(:, 1:R);
    B = abs(U2(:, 1:R));
    C = abs(U3(:, 1:R));
    B(B == 0) = eps;
    C(C == 0) = eps;
end

function [A, F] = absorb_norm_into_left(A, F)
    for r = 1:size(F, 2)
        s = norm(F(:, r));
        if s > 0
            F(:, r) = F(:, r) / s;
            A(:, r) = A(:, r) * s;
        end
    end
end

function K = khatri_rao_local(A, B)
    R = size(A, 2);
    K = zeros(size(A, 1)*size(B, 1), R);
    for r = 1:R
        K(:, r) = kron(A(:, r), B(:, r));
    end
end

function Xhat = reconstruct_cp_local(A, B, C)
    [nQ, R] = size(A);
    nT = size(B, 1);
    nD = size(C, 1);
    Xhat = zeros(nQ, nT, nD);
    for r = 1:R
        Xhat = Xhat + reshape(A(:, r), nQ, 1, 1) .* reshape(B(:, r), 1, nT, 1) .* reshape(C(:, r), 1, 1, nD);
    end
end

function Z = mode_surface_matrix(B, C)
    nT = size(B, 1);
    nD = size(C, 1);
    R = size(B, 2);
    Z = zeros(nT*nD, R);
    for k = 1:nD
        idx = (1:nT) + (k-1)*nT;
        Z(idx, :) = B .* C(k, :);
    end
end

function tau = estimate_factor_tau(t, F)
    t = t(:);
    R = size(F, 2);
    tau = nan(1, R);
    for r = 1:R
        y = F(:, r);
        y = y / max(y);
        good = isfinite(t) & isfinite(y) & t >= 0 & y > 0.05 & y < 0.98;
        if nnz(good) < 4
            good = isfinite(t) & isfinite(y) & t >= 0 & y > 0.02;
        end
        if nnz(good) < 2
            tau(r) = NaN;
            continue;
        end
        obj = @(p) sum((y(good) - exp(-t(good) ./ 10.^p)).^2);
        p0 = log10(max(median(t(good)), eps));
        tau(r) = 10.^fminsearch(obj, p0, optimset("Display", "off"));
    end
end

function rot = closure_rotate_tensor_modes(X, Z0, E0, cfg)
    rng(cfg.rotationRandomSeed);
    R = size(Z0, 2);
    [nQ, nT, nD] = size(X);
    Xmat = reshape(X, nQ, nT*nD);
    [~, zeroLagIndex] = min(abs(cfg.delay_s));
    wSurface = make_surface_fit_weights(nT, nD, zeroLagIndex, cfg);
    normX2 = sum(sum((Xmat.^2) .* wSurface.'));
    XZ0 = (Xmat .* wSurface.') * Z0;
    G0 = Z0.' * (Z0 .* wSurface);
    D0 = X(:, :, zeroLagIndex);
    normD02 = sum(D0(:).^2);

    starts = make_rotation_starts(Z0, E0, R, cfg);
    opts = optimset("Display", "off", "MaxIter", cfg.rotationMaxIter, ...
        "MaxFunEvals", cfg.rotationMaxIter * 30);

    bestScore = inf;
    bestT = eye(R);
    bestStats = struct();
    for iStart = 1:size(starts, 3)
        p0 = starts(:, :, iStart);
        obj = @(p) rotation_objective(reshape(p, R, R), Z0, E0, XZ0, G0, normX2, D0, normD02, cfg);
        [pOpt, score] = fminsearch(obj, p0(:), opts);
        if score < bestScore
            bestScore = score;
            bestT = reshape(pOpt, R, R);
            [~, bestStats] = rotation_objective(bestT, Z0, E0, XZ0, G0, normX2, D0, normD02, cfg);
        end
    end

    T = bestStats.T_scaled;
    Z = Z0 * T;
    C_at_t = E0 * T;
    C_ground = max(1 - sum(max(C_at_t, 0), 2), 0);
    G = T.' * G0 * T + 1e-10 * eye(R);
    S_q = (XZ0 * T) / G;
    XhatMat = S_q * Z.';
    Xhat = reshape(XhatMat, nQ, nT, nD);

    [~, order] = sort_component_profiles(C_at_t);
    C_at_t = C_at_t(:, order);
    Z = Z(:, order);
    S_q = S_q(:, order);
    T = T(:, order);
    XhatMat = S_q * Z.';
    Xhat = reshape(XhatMat, nQ, nT, nD);

    [~, zeroLagIndex] = min(abs(cfg.delay_s));
    qt_target = X(:, :, zeroLagIndex);
    qt_recon = Xhat(:, :, zeroLagIndex);
    C_raw_at_t = max(C_at_t, 0);
    C_ground_raw_at_t = max(1 - sum(C_raw_at_t, 2), 0);
    C_surface_raw = max(Z, 0);
    C_surface_use = C_surface_raw;
    C_ground_surface = max(1 - sum(C_surface_use, 2), 0);
    S_q_tensor = S_q;
    zero_lag_relerr_tensor = norm(qt_target(:) - qt_recon(:)) / max(norm(qt_target(:)), eps);
    raw_tensor_relerr = norm(X(:) - Xhat(:)) / max(norm(X(:)), eps);

    if isfield(cfg, "applyTerminalGroundProjection") && cfg.applyTerminalGroundProjection
        [C_use_at_t, C_ground_use_at_t, projectionInfo] = terminal_ground_projection(C_raw_at_t, cfg, cfg.tensor_time_s(:));

        if isfield(cfg, "solveSpectraOnProjectedSurface") && cfg.solveSpectraOnProjectedSurface
            uSurface = surface_shifted_time_axis(nT, nD, cfg);
            [C_surface_use, C_ground_surface, surfaceProjectionInfo] = terminal_ground_projection( ...
                C_surface_raw, cfg, uSurface, projectionInfo.terminal_time_s);
            Guse = C_surface_use.' * (C_surface_use .* wSurface) + 1e-10 * eye(R);
            S_q = ((Xmat .* wSurface.') * C_surface_use) / Guse;
            XhatMat = S_q * C_surface_use.';
            Xhat = reshape(XhatMat, nQ, nT, nD);
            zeroRows = (1:nT) + (zeroLagIndex - 1) * nT;
            C_use_at_t = C_surface_use(zeroRows, :);
            C_ground_use_at_t = C_ground_surface(zeroRows);
            qt_recon = Xhat(:, :, zeroLagIndex);
            projectionInfo.surface = surfaceProjectionInfo;
        else
            Guse = C_use_at_t.' * C_use_at_t + 1e-10 * eye(R);
            S_q = (qt_target * C_use_at_t) / Guse;
            qt_recon = S_q * C_use_at_t.';
        end
    else
        C_use_at_t = C_raw_at_t;
        C_ground_use_at_t = C_ground_raw_at_t;
        projectionInfo = struct("applied", false);
    end

    rot = struct();
    rot.T = T;
    rot.S_q = S_q;
    rot.S_q_tensor = S_q_tensor;
    rot.C_at_t = C_use_at_t;
    rot.C_ground_at_t = C_ground_use_at_t;
    rot.C_at_t_raw = C_raw_at_t;
    rot.C_ground_at_t_raw = C_ground_raw_at_t;
    rot.C_surface = C_surface_use;
    rot.C_surface_raw = C_surface_raw;
    rot.C_ground_surface = C_ground_surface;
    rot.Xhat = Xhat;
    rot.raw_tensor_relerr = raw_tensor_relerr;
    rot.tensor_relerr = norm(X(:) - Xhat(:)) / max(norm(X(:)), eps);
    rot.zero_lag_relerr = norm(qt_target(:) - qt_recon(:)) / max(norm(qt_target(:)), eps);
    rot.zero_lag_relerr_tensor = zero_lag_relerr_tensor;
    rot.closure_max_error = max(abs(sum(rot.C_at_t, 2) + rot.C_ground_at_t - 1));
    rot.overlap = concentration_overlap_metric(rot.C_at_t);
    rot.projection = projectionInfo;
    rot.objective = bestScore;
    rot.objective_stats = bestStats;
end

function [Cproj, Gproj, info] = terminal_ground_projection(C, cfg, rowTime, terminalTime)
    Cproj = max(C, 0);
    nRows = size(Cproj, 1);
    if nargin < 3 || isempty(rowTime)
        rowTime = (1:nRows).';
    else
        rowTime = rowTime(:);
    end
    [~, peakIdx] = max(Cproj, [], 1);
    if nargin < 4 || isempty(terminalTime) || ~isfinite(terminalTime)
        [terminalTime, iLastPeak] = max(rowTime(peakIdx));
        terminalIdx = peakIdx(iLastPeak);
    else
        [~, terminalIdx] = min(abs(rowTime - terminalTime));
    end
    activeSum = sum(Cproj, 2);
    preMask = rowTime <= terminalTime + eps(max(1, abs(terminalTime)));
    scaleRows = find(preMask & activeSum > eps);
    Cproj(scaleRows, :) = Cproj(scaleRows, :) ./ activeSum(scaleRows);
    Gproj = max(1 - sum(Cproj, 2), 0);
    Gproj(preMask) = 0;

    info = struct();
    info.applied = true;
    info.terminal_index = terminalIdx;
    if isfield(cfg, "tensor_time_s")
        info.terminal_time_s = cfg.tensor_time_s(terminalIdx);
    else
        info.terminal_time_s = NaN;
    end
    info.terminal_time_s = terminalTime;
    info.raw_ground_max_before_terminal = max(max(1 - activeSum(preMask), 0));
end

function uSurface = surface_shifted_time_axis(nT, nD, cfg)
    t = cfg.tensor_time_s(:);
    uSurface = zeros(nT * nD, 1);
    for k = 1:nD
        idx = (1:nT) + (k - 1) * nT;
        uSurface(idx) = shifted_time_local(t, cfg.delay_s(k), cfg);
    end
end

function starts = make_rotation_starts(Z0, E0, R, cfg)
    nStart = max(cfg.rotationNumStarts, 3);
    starts = zeros(R, R, nStart);
    starts(:, :, 1) = eye(R);
    starts(:, :, 2) = simplex_rotation_start(E0);
    starts(:, :, 3) = simplex_rotation_start(Z0);
    for i = 4:nStart
        [Q, ~] = qr(randn(R, R), 0);
        starts(:, :, i) = Q + 0.25 * randn(R, R);
    end
end

function T = simplex_rotation_start(Z)
    R = size(Z, 2);
    Z = max(Z, 0);
    idx = zeros(1, R);
    [~, idx(1)] = max(sum(Z.^2, 2));
    basis = Z(idx(1), :);
    for r = 2:R
        P = basis' / (basis*basis' + 1e-12) * basis;
        residual = Z - Z * P;
        [~, idx(r)] = max(sum(residual.^2, 2));
        basis = Z(idx(1:r), :);
    end
    T = pinv(Z(idx, :));
end

function [score, stats] = rotation_objective(T, Z0, E0, XZ0, G0, normX2, D0, normD02, cfg)
    Zraw = Z0 * T;
    Craw = E0 * T;
    scale = max(sum(max(Zraw, 0), 2));
    if ~isfinite(scale) || scale <= eps
        score = inf;
        stats = struct();
        return;
    end
    T = T / scale;
    Z = Z0 * T;
    C = E0 * T;
    Cpos = max(C, 0);
    if isfield(cfg, "applyTerminalGroundProjection") && cfg.applyTerminalGroundProjection
        [Cphys, ~, ~] = terminal_ground_projection(Cpos, cfg);
    else
        Cphys = Cpos;
    end

    G = T.' * G0 * T + 1e-10 * eye(size(T, 2));
    XZ = XZ0 * T;
    fitTerm = trace((XZ / G) * XZ.');
    rss = max(normX2 - fitTerm, 0);
    relerr = sqrt(rss / max(normX2, eps));

    G0lag = Cphys.' * Cphys + 1e-10 * eye(size(T, 2));
    D0C = D0 * Cphys;
    fitTerm0 = trace((D0C / G0lag) * D0C.');
    rss0 = max(normD02 - fitTerm0, 0);
    zeroLagRelerr = sqrt(rss0 / max(normD02, eps));

    negPenalty = mean(min(Z(:), 0).^2) + mean(min(C(:), 0).^2);
    sumZ = sum(max(Z, 0), 2);
    sumC = sum(Cphys, 2);
    closurePenalty = mean(max(sumZ - 1, 0).^2) + mean(max(sumC - 1, 0).^2) + (max(sumZ) - 1)^2;
    endpointPenalty = (sumC(1) - 1)^2 + sumC(end)^2;
    componentPeaks = max(Cphys, [], 1);
    peakPenalty = mean(max(cfg.rotationMinComponentPeak - componentPeaks, 0).^2) ...
        + 0.1 * mean(max(cfg.rotationPeakTarget - componentPeaks, 0).^2);
    peakPurity = component_peak_purity(Cphys);
    purityPenalty = mean(max(cfg.rotationMinPeakPurity - peakPurity, 0).^2);
    activeSumMonotonicPenalty = mean(max(diff(sumC), 0).^2);
    [~, peakIdx] = max(Cphys, [], 1);
    lastActivePeakIdx = max(peakIdx);
    ground = max(1 - sumC, 0);
    preTerminalGround = ground(1:lastActivePeakIdx);
    earlyGroundPenalty = mean(preTerminalGround.^2);
    earlyGroundMaxPenalty = max(preTerminalGround).^2;
    peakGroundPenalty = mean(ground(peakIdx).^2);
    peakSpreadPenalty = component_peak_spread_penalty(Cphys, cfg);
    initialStatePenalty = initial_single_state_penalty(Cphys, peakIdx);
    [generatorPenalty, generatorStats] = generator_consistency_penalty(Cphys, cfg);
    overlap = concentration_overlap_metric(Cphys);
    condPenalty = log10(cond(T) + 1)^2;

    score = relerr + cfg.rotationNegativePenalty * negPenalty ...
        + cfg.rotationZeroLagWeight * zeroLagRelerr ...
        + cfg.rotationClosurePenalty * closurePenalty ...
        + cfg.rotationEndpointPenalty * endpointPenalty ...
        + cfg.rotationPeakPenalty * peakPenalty ...
        + cfg.rotationPurityPenalty * purityPenalty ...
        + cfg.rotationActiveSumMonotonicPenalty * activeSumMonotonicPenalty ...
        + cfg.rotationEarlyGroundPenalty * earlyGroundPenalty ...
        + cfg.rotationEarlyGroundMaxPenalty * earlyGroundMaxPenalty ...
        + cfg.rotationPeakGroundPenalty * peakGroundPenalty ...
        + cfg.rotationPeakSpreadPenalty * peakSpreadPenalty ...
        + cfg.rotationInitialStatePenalty * initialStatePenalty ...
        + cfg.rotationGeneratorPenalty * generatorPenalty ...
        + cfg.rotationNonterminalGroundRatePenalty * generatorStats.nonterminalGroundFraction ...
        + cfg.rotationOverlapWeight * overlap ...
        + cfg.rotationConditionPenalty * condPenalty;

    stats = struct();
    stats.T_scaled = T;
    stats.relerr = relerr;
    stats.zeroLagRelerr = zeroLagRelerr;
    stats.negPenalty = negPenalty;
    stats.closurePenalty = closurePenalty;
    stats.endpointPenalty = endpointPenalty;
    stats.peakPenalty = peakPenalty;
    stats.componentPeaks = componentPeaks;
    stats.peakPurity = peakPurity;
    stats.purityPenalty = purityPenalty;
    stats.activeSumMonotonicPenalty = activeSumMonotonicPenalty;
    stats.earlyGroundPenalty = earlyGroundPenalty;
    stats.earlyGroundMaxPenalty = earlyGroundMaxPenalty;
    stats.peakGroundPenalty = peakGroundPenalty;
    stats.peakSpreadPenalty = peakSpreadPenalty;
    stats.initialStatePenalty = initialStatePenalty;
    stats.generatorPenalty = generatorPenalty;
    stats.generatorStats = generatorStats;
    stats.overlap = overlap;
    stats.condPenalty = condPenalty;
end

function [penalty, stats] = generator_consistency_penalty(C, cfg)
    C = max(C, 0);
    [nT, R] = size(C);
    if nT < 4 || R < 2
        penalty = 0;
        stats = struct("relerr", 0, "nonterminalGroundFraction", 0, "rates", []);
        return;
    end

    t = cfg.tensor_time_s(:);
    dCdt = zeros(size(C));
    for r = 1:R
        dCdt(:, r) = gradient(C(:, r), t);
    end

    [~, peakIdx] = max(C, [], 1);
    tWeight = max(t, min(t(t > 0))).^cfg.rotationGeneratorTimeWeightAlpha;
    A = [];
    sourceIndex = zeros(0, 1);
    isGroundEdge = false(0, 1);
    for src = 1:R
        for dst = 1:R
            if dst == src || peakIdx(dst) <= peakIdx(src)
                continue;
            end
            col = zeros(nT, R);
            col(:, src) = -tWeight .* C(:, src);
            col(:, dst) = tWeight .* C(:, src);
            A(:, end+1) = col(:); %#ok<AGROW>
            sourceIndex(end+1, 1) = src; %#ok<AGROW>
            isGroundEdge(end+1, 1) = false; %#ok<AGROW>
        end
        col = zeros(nT, R);
        col(:, src) = -tWeight .* C(:, src);
        A(:, end+1) = col(:); %#ok<AGROW>
        sourceIndex(end+1, 1) = src; %#ok<AGROW>
        isGroundEdge(end+1, 1) = true; %#ok<AGROW>
    end

    B = tWeight .* dCdt;
    b = B(:);
    if isempty(A) || norm(b) <= eps
        rates = zeros(size(A, 2), 1);
        relerr = 0;
    else
        rates = lsqnonneg(A, b, optimset("Display", "off", "MaxIter", 500));
        relerr = norm(A * rates - b) / max(norm(b), eps);
    end

    lastPeakSource = find(peakIdx == max(peakIdx), 1);
    nonterminalGround = isGroundEdge & sourceIndex ~= lastPeakSource;
    nonterminalGroundFraction = sum(rates(nonterminalGround)) / max(sum(rates), eps);

    penalty = relerr;
    stats = struct();
    stats.relerr = relerr;
    stats.nonterminalGroundFraction = nonterminalGroundFraction;
    stats.rates = rates;
end

function penalty = initial_single_state_penalty(C, peakIdx)
    [~, order] = sort(peakIdx, "ascend");
    firstComp = order(1);
    startFrac = C(1, :) ./ (sum(C(1, :)) + eps);
    other = setdiff(1:size(C, 2), firstComp);
    penalty = (1 - startFrac(firstComp)).^2 + sum(startFrac(other).^2);
end

function penalty = component_peak_spread_penalty(C, cfg)
    [~, peakIdx] = max(C, [], 1);
    if isfield(cfg, "tensor_time_s")
        x = log10(max(cfg.tensor_time_s(:), realmin));
    else
        x = (1:size(C, 1)).';
    end
    peakX = sort(x(peakIdx), "ascend");
    gapPenalty = 0;
    if numel(peakX) > 1
        gapPenalty = mean(max(cfg.rotationMinPeakLogGap - diff(peakX), 0).^2);
    end
    span = max(x) - min(x);
    minLastPeakX = min(x) + cfg.rotationMinLastPeakLogFraction * span;
    lastPeakPenalty = max(minLastPeakX - max(peakX), 0).^2;
    penalty = gapPenalty + lastPeakPenalty;
end

function wSurface = make_surface_fit_weights(nT, nD, zeroLagIndex, cfg)
    wDelay = ones(nD, 1);
    delay = cfg.delay_s(:);
    if isfield(cfg, "rotationShortDelaySolveWeight") && cfg.rotationShortDelaySolveWeight > 0
        scale = max(cfg.rotationShortDelayScale_s, eps);
        wDelay = wDelay + cfg.rotationShortDelaySolveWeight ./ (1 + delay ./ scale);
    end
    if isfield(cfg, "rotationZeroLagSolveWeight") && cfg.rotationZeroLagSolveWeight > 1
        wDelay(zeroLagIndex) = wDelay(zeroLagIndex) * cfg.rotationZeroLagSolveWeight;
    end
    wSurface = repelem(wDelay, nT);
end

function purity = component_peak_purity(C)
    R = size(C, 2);
    purity = zeros(1, R);
    rowSum = sum(C, 2) + eps;
    for r = 1:R
        [~, idx] = max(C(:, r));
        purity(r) = C(idx, r) / rowSum(idx);
    end
end

function [peakIdx, order] = sort_component_profiles(C)
    [~, peakIdx] = max(C, [], 1);
    [~, order] = sort(peakIdx, "ascend");
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

function topology = infer_sparse_active_ground_topology(t, C, cfg)
    C = max(C, 0);
    R = size(C, 2);
    dCdt = zeros(size(C));
    for r = 1:R
        dCdt(:, r) = gradient(C(:, r), t);
    end
    [~, peakIdx] = max(C, [], 1);
    tWeight = max(t(:), min(t(t > 0))).^cfg.topologyTimeWeightAlpha;

    rateNames = strings(0, 1);
    sourceIndex = zeros(0, 1);
    A = [];
    for src = 1:R
        for dst = 1:R
            if dst == src
                continue;
            end
            if cfg.topologyUsePeakOrderedAcyclic && peakIdx(dst) <= peakIdx(src)
                continue;
            end
            col = zeros(numel(t), R);
            col(:, src) = -tWeight .* C(:, src);
            col(:, dst) = tWeight .* C(:, src);
            A(:, end+1) = col(:); %#ok<AGROW>
            rateNames(end+1, 1) = sprintf("%d->%d", src, dst); %#ok<AGROW>
            sourceIndex(end+1, 1) = src; %#ok<AGROW>
        end
        col = zeros(numel(t), R);
        col(:, src) = -tWeight .* C(:, src);
        A(:, end+1) = col(:); %#ok<AGROW>
        rateNames(end+1, 1) = sprintf("%d->ground", src); %#ok<AGROW>
        sourceIndex(end+1, 1) = src; %#ok<AGROW>
    end
    b = (tWeight .* dCdt);
    b = b(:);
    rates = lsqnonneg(A, b, optimset("Display", "off", "MaxIter", 2000));
    for src = 1:R
        idx = sourceIndex == src;
        maxSourceRate = max(rates(idx));
        if maxSourceRate > 0
            rates(idx & rates < cfg.topologyMinRelativeSourceRate * maxSourceRate) = 0;
        end
    end
    [sortedRates, order] = sort(rates, "descend");

    topology = struct();
    topology.peak_index = peakIdx;
    topology.time_weight_alpha = cfg.topologyTimeWeightAlpha;
    topology.uses_peak_ordered_acyclic_edges = cfg.topologyUsePeakOrderedAcyclic;
    topology.min_relative_source_rate = cfg.topologyMinRelativeSourceRate;
    topology.rate_names = rateNames;
    topology.rates = rates;
    topology.sorted_names = rateNames(order);
    topology.sorted_rates = sortedRates;
    topology.rate_table = table(cellstr(rateNames(order)), sortedRates(:), ...
        'VariableNames', {'transition', 'rate_s_inv'});
end

function plot_results(q, t, delay_s, X, qt_target, qt_cp, qt_rot, cp, A_eig, B_eig, C_eig, rot, topology)
    figure("Name", "01_CP_tensor_fit_history");
    semilogy(cp.fit_history, "LineWidth", 1.4);
    grid on;
    xlabel("ALS iteration");
    ylabel("relative error");
    title(sprintf("CP tensor decomposition fit %.5f", cp.fit));

    figure("Name", "02_CP_lifetime_eigenmodes");
    tiledlayout(1, 3, "TileSpacing", "compact");
    nexttile;
    plot(t, B_eig, "LineWidth", 1.3);
    grid on; xlabel("Time (s)"); ylabel("B_r(t)");
    title("CP time factors");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile;
    plot(delay_s, C_eig, "LineWidth", 1.3);
    grid on; xlabel("Delay (s)"); ylabel("C_r(delay)");
    title("CP delay factors");
    if all(delay_s > 0), set(gca, "XScale", "log"); end
    nexttile;
    plot(q, A_eig, "LineWidth", 1.3);
    hold on; yline(0, "k--"); hold off;
    grid on; xlabel("q"); ylabel("Eigen spectrum");
    title("CP eigen spectra");

    figure("Name", "03_CP_zero_lag_reconstruction");
    clim_all = [min([qt_target(:); qt_cp(:)]), max([qt_target(:); qt_cp(:)])];
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile; imagesc(t, q, qt_target); axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q"); title("Zero-lag target");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile; imagesc(t, q, qt_cp); axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q"); title(sprintf("CP zero-lag fit %.4f", 1 - norm(qt_target(:)-qt_cp(:))/norm(qt_target(:))));
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "04_Closure_rotated_concentrations");
    plot(t, [rot.C_at_t, rot.C_ground_at_t], "LineWidth", 1.5);
    grid on;
    xlabel("Time (s)");
    ylabel("Concentration");
    title("Closure rotation from unique tensor eigenmodes");
    legend([compose("component %d", 1:size(rot.C_at_t, 2)), "silent ground"], "Location", "best");
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "05_Closure_rotated_q_spectra");
    plot(q, rot.S_q, "LineWidth", 1.5);
    hold on; yline(0, "k--"); hold off;
    grid on;
    xlabel("q");
    ylabel("Species spectrum (a.u.)");
    title("Closure-rotated active species spectra");
    legend(compose("component %d", 1:size(rot.S_q, 2)), "Location", "best");

    figure("Name", "06_Closure_zero_lag_reconstruction");
    clim_all = [min([qt_target(:); qt_rot(:)]), max([qt_target(:); qt_rot(:)])];
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile; imagesc(t, q, qt_target); axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q"); title("Zero-lag target");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile; imagesc(t, q, qt_rot); axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q"); title(sprintf("Closure rotation fit %.4f", 1 - rot.zero_lag_relerr));
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "07_Inferred_sparse_topology");
    nShow = min(9, numel(topology.sorted_rates));
    bar(topology.sorted_rates(1:nShow));
    grid on;
    set(gca, "XTick", 1:nShow, "XTickLabel", topology.sorted_names(1:nShow));
    xtickangle(35);
    ylabel("rate (s^{-1})");
    title("Sparse active-ground transition fit after tensor rotation");
end

function u = shifted_time_local(t, delay, cfg)
    switch cfg.shiftMode
        case "delay_signal"
            u = t(:) - delay;
        case "evaluate_later_time"
            u = t(:) + delay;
        otherwise
            error("알 수 없는 cfg.shiftMode: %s", cfg.shiftMode);
    end
end

function V_model = eval_fitted_rsv(parStruct, tQuery)
    tQuery = tQuery(:);
    A = parStruct.A;
    tau = parStruct.tau(:);
    A0 = parStruct.A0(:);
    FWHM = parStruct.FWHM;
    t0 = parStruct.t0;
    numRSV = size(A, 1);
    V_model = zeros(numel(tQuery), numRSV);
    for j = 1:numRSV
        V_model(:, j) = conv_gauss_expdec_local(tQuery, A(j, :), A0(j), tau, t0, FWHM);
    end
end

function y = conv_gauss_expdec_local(t, A, A0, tau, t0, FWHM)
    t = t(:);
    tau = tau(:);
    sigma = FWHM / 2.3548;
    y = A0 * ones(numel(t), 1);
    inv_sqrt2_sigma = 1 / (sqrt(2) * sigma);
    gauss_env = exp(-((t - t0).^2) ./ (2 * sigma^2));
    for k = 1:numel(tau)
        Z = ((sigma^2)/tau(k) - (t - t0)) * inv_sqrt2_sigma;
        comp = zeros(numel(t), 1);
        idx_pos = (Z >= 0);
        if any(idx_pos)
            comp(idx_pos) = gauss_env(idx_pos) .* erfcx(Z(idx_pos));
        end
        idx_neg = (Z < 0);
        if any(idx_neg)
            E = (sigma^2)/(2*tau(k)^2) - (t(idx_neg) - t0)/tau(k);
            comp(idx_neg) = exp(E) .* erfc(Z(idx_neg));
        end
        y = y + 0.5 * A(k) .* comp;
    end
end

function M = read_numeric_table(fileName)
    if exist(fileName, "file") ~= 2
        error("파일을 찾을 수 없습니다: %s", fileName);
    end
    txt = fileread(fileName);
    lines = regexp(txt, "\r\n|\n|\r", "split");
    rows = {};
    maxCols = 0;
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if isempty(line), continue; end
        line = regexprep(line, "[,;\t]", " ");
        vals = sscanf(line, "%f").';
        if isempty(vals), continue; end
        rows{end+1, 1} = vals; %#ok<AGROW>
        maxCols = max(maxCols, numel(vals));
    end
    M = nan(numel(rows), maxCols);
    for i = 1:numel(rows)
        M(i, 1:numel(rows{i})) = rows{i};
    end
    if any(isnan(M(:)))
        M = readmatrix(fileName);
    end
end

function ensure_asinh_helpers_exist()
    if exist("asinhspace", "file") ~= 2
        fid = fopen("asinhspace.m", "w");
        if fid >= 0
            fprintf(fid, "function x = asinhspace(a,b,n)\nif nargin<3, n=100; end\nx=sinh(linspace(asinh(a),asinh(b),n));\nend\n");
            fclose(fid);
        end
    end
    if exist("asinhplot", "file") ~= 2
        fid = fopen("asinhplot.m", "w");
        if fid >= 0
            fprintf(fid, "function h = asinhplot(x,y,varargin)\nh=plot(x,y,varargin{:});\nend\n");
            fclose(fid);
        end
    end
end
