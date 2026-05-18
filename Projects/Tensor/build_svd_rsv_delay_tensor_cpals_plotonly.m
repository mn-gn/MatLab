%% build_svd_rsv_delay_tensor_cpals_plotonly.m
% q / time / q-delta-S 데이터를 읽고,
% 1) SVD 수행
% 2) 주성분 RSV = V(:,1:numPC)를 fit_rsv로 피팅
% 3) 피팅된 RSV 모델로 시간 지연 텐서 D(q,t,delay) 생성
% 4) 시간축/tau축 factor에 비음수 제약을 둔 CP-ALS 수행
% 5) CP 결과를 다시 q-t 2차원 데이터로 재구성
% 6) 아래 그림만 표시
%    - q축 intermediate spectra
%    - t축 concentration profiles (nonnegative)
%    - 재구성된 q-t 2D map
%
% 주의:
%   - 결과 파일(MAT/CSV)은 저장하지 않습니다.
%   - q-mode factor(중간체 스펙트럼)는 비음수 제약을 주지 않으므로
%     음수 값을 가질 수 있습니다.
%   - fit_rsv.m 내부에서 생성되는 중간 figure는 모두 닫고,
%     최종 figure만 다시 그립니다.

clear; clc;

%% ===================== 사용자 설정 =====================
cfg = struct();

% 입력 파일명
cfg.qFile    = "q.csv";
cfg.tFile    = "t_81points.csv";
cfg.dataFile = "q_delta_S_noise_10%_81points.csv";

% SVD / RSV 피팅 설정
cfg.numPC    = 3;
cfg.numTau   = 3;
cfg.numStart = 100;

% 데이터 전처리
cfg.centerData = false;
cfg.fixSVDSign = true;
cfg.useLogTimeGridForTensor = true;
cfg.numTensorTimePoints = [];  % []이면 원본 time point 개수와 동일하게 logspace 재구성

% 시간 지연 목록 [초]
% [0s, 0.01 ps, 0.1 ps, 1 ps, 10 ps, 100 ps, 1000 ps, 10 ns, 100 ns]
delay_min_s = 1e-14;   % 0.01 ps
delay_max_s = 100e-9;  % 100 ns
num_log_delay = 100;   % one side log-spaced lag count
cfg.useSignedTauAxis = false;
cfg.outOfRangeWeight = 1;
cfg.zeroSignalBeforeTimeZero = false;
cfg.signalTimeZero_s = 0;

pos_delay_s = logspace(log10(delay_min_s), log10(delay_max_s), num_log_delay);
if cfg.useSignedTauAxis
    cfg.delay_s = [-fliplr(pos_delay_s), 0, pos_delay_s];
else
    cfg.delay_s = [0, pos_delay_s];
end

% delay label 자동 생성
cfg.delay_labels = make_delay_labels(cfg.delay_s);

% 시간 지연 방향
%   "delay_signal"       : D_delay(t) = D_model(t - delay)
%   "evaluate_later_time": D_delay(t) = D_model(t + delay)
cfg.shiftMode = "evaluate_later_time";

% CP-ALS 설정
cfg.doCPALS      = true;
cfg.decompMethod = "kinetic_tensor_scan";  % "kinetic_tensor_scan": fit the whole shifted-time tensor, "kinetic_scan": q-t model scan, "shift_mcr": blind shift MCR, "sequential_varpro": validation, "cp_als": separable CP
cfg.cpRank       = 3;
cfg.cpMaxIters   = 500;
cfg.cpTol        = 1e-8;
cfg.cpInit       = "svd";
cfg.cpRandomSeed = 1;
cfg.cpVerbose    = true;
cfg.cpNonnegativeModes = ["time", "tau"];
cfg.shiftMcrNumGrid = 500;
cfg.shiftMcrNumStarts = 20;
cfg.shiftMcrUseNNLS = false;
cfg.shiftMcrShapeConstraint = "unimodal_smooth";  % "none" or "unimodal_smooth"; does not impose a reaction topology.
cfg.shiftMcrSmoothWindow = 9;
cfg.sequentialFitRates = true;
cfg.sequentialNumStarts = 20;
cfg.sequentialTauLower_s = 1e-14;
cfg.sequentialTauUpper_s = 1e-6;
cfg.kineticScanModels = ["sequential", "branch"];  % closure-compatible species models. "parallel" is DAS-like, not a closed species model.
cfg.kineticScanNumStarts = 24;
cfg.kineticScanMaxIter = 400;
cfg.kineticScanTauLower_s = 1e-14;
cfg.kineticScanTauUpper_s = 1e-6;
cfg.kineticUseIRFConvolution = true;
cfg.kineticTensorSelectionMetric = "aic";  % "aic", "bic", "relerr", or "overlap"
cfg.kineticUseSilentGroundClosure = true;
cfg.kineticPlotSilentGround = true;

% 3차원 텐서를 2차원 q-t 데이터로 다시 줄이는 방법
%   "mean_delay" : delay 차원 평균
%   "sum_delay"  : delay 차원 합
%   "delay_index": 특정 delay slice 사용
cfg.qtCollapseMode = "mean_delay";
cfg.plotDelayIndex = 1;   % qtCollapseMode = "delay_index"일 때 사용

% 비교 그림을 같이 볼지 여부
cfg.showComparisonMap = true;

%% ===================== 경로 / helper 확인 =====================
scriptDir = fileparts(mfilename('fullpath'));
if strlength(scriptDir) > 0
    cd(scriptDir);
end

if exist("fit_rsv", "file") ~= 2
    error("fit_rsv.m 파일이 현재 폴더 또는 MATLAB path에 없습니다.");
end

ensure_asinh_helpers_exist();

%% ===================== 데이터 읽기 =====================
q = read_numeric_table(cfg.qFile);
t = read_numeric_table(cfg.tFile);
D = read_numeric_table(cfg.dataFile);

q = q(:);
t = t(:);

% D가 q x t 형태가 되도록 자동 확인/전치
if size(D,1) == numel(q) && size(D,2) == numel(t)
    % OK: D(q,t)
elseif size(D,1) == numel(t) && size(D,2) == numel(q)
    D = D.';
    warning("입력 데이터가 t x q 형태로 감지되어 q x t로 전치했습니다.");
else
    error("데이터 크기가 q/t와 맞지 않습니다. size(D)=[%d %d], length(q)=%d, length(t)=%d", ...
        size(D,1), size(D,2), numel(q), numel(t));
end

[nQ, nT] = size(D);
fprintf("Loaded data: D = %d(q) x %d(time)\n", nQ, nT);

t_data = t;
D_data = D; %#ok<NASGU>
nT_data = nT;

%% ===================== SVD =====================
if cfg.centerData
    D_mean_q = mean(D, 2);
    D_work = D - D_mean_q;
else
    D_mean_q = zeros(nQ, 1);
    D_work = D;
end

[U, S, V] = svd(D_work, "econ");
numPC = min([cfg.numPC, size(U,2), size(V,2)]);
U_pc = U(:, 1:numPC);
S_pc = S(1:numPC, 1:numPC);
V_pc = V(:, 1:numPC);   % time x numPC

% SVD 부호 정렬: RSV의 최대 절대값 지점이 양수
if cfg.fixSVDSign
    for k = 1:numPC
        [~, idxMax] = max(abs(V_pc(:, k)));
        if V_pc(idxMax, k) < 0
            U_pc(:, k) = -U_pc(:, k);
            V_pc(:, k) = -V_pc(:, k);
        end
    end
end

%% ===================== RSV 피팅 =====================
[best_par_struct, V_fit, fit_save_struct] = fit_rsv(t, V_pc, cfg.numTau, cfg.numStart); %#ok<NASGU>

% fit_rsv.m 내부 figure는 닫고, 아래에서 최종 그림만 다시 그림
close all;

if cfg.useLogTimeGridForTensor
    t_pos = t_data(t_data > 0);
    if isempty(t_pos)
        error("logspace tensor time grid를 만들려면 양수 time point가 필요합니다.");
    end

    nTensorTime = cfg.numTensorTimePoints;
    if isempty(nTensorTime)
        nTensorTime = nT_data;
    end

    t = logspace(log10(min(t_pos)), log10(max(t_pos)), nTensorTime).';
    nT = numel(t);
    V_fit_tensor_time = eval_fitted_rsv(best_par_struct, t);
    D_tensor_time = U_pc * S_pc * V_fit_tensor_time.';
    D_tensor_time = D_tensor_time + D_mean_q; %#ok<NASGU>

    fprintf("Tensor time axis reconstructed on logspace: %d points, %.4g s to %.4g s\n", ...
        nT, min(t), max(t));
else
    V_fit_tensor_time = V_fit; %#ok<NASGU>
    D_tensor_time = D; %#ok<NASGU>
end

%% ===================== 시간 지연 텐서 생성 =====================
nDelay = numel(cfg.delay_s);
delay_tensor = zeros(nQ, nT, nDelay);
V_fit_shifted = zeros(nT, numPC, nDelay);
valid_weight_time_tau = zeros(nT, nDelay);

for iDelay = 1:nDelay
    delay = cfg.delay_s(iDelay);
    t_eval = shifted_time_local(t, delay, cfg);

    if cfg.zeroSignalBeforeTimeZero
        preTimeZero = t_eval < cfg.signalTimeZero_s;
    else
        preTimeZero = false(size(t_eval));
    end

    isValidEval = true(size(t_eval));
    valid_weight_time_tau(:, iDelay) = cfg.outOfRangeWeight;
    valid_weight_time_tau(isValidEval, iDelay) = 1;

    t_eval_model = t_eval;
    t_eval_model(preTimeZero) = cfg.signalTimeZero_s;
    V_shift = eval_fitted_rsv(best_par_struct, t_eval_model);
    V_shift(preTimeZero, :) = 0;
    V_fit_shifted(:, :, iDelay) = V_shift;

    D_shift = U_pc * S_pc * V_shift.';
    D_shift = D_shift + D_mean_q;
    D_shift(:, preTimeZero) = 0;
    delay_tensor(:, :, iDelay) = D_shift;
end

fprintf("Delay tensor created: %d(q) x %d(time) x %d(delay)\n", size(delay_tensor));

%% ===================== 텐서 분해 =====================
switch lower(string(cfg.decompMethod))
    case "kinetic_tensor_scan"
        [kin_tensor_result, cp_reconstructed_tensor] = first_order_kinetic_tensor_scan( ...
            delay_tensor, valid_weight_time_tau, t, cfg.delay_s, best_par_struct, cfg);

        cp_result = kin_tensor_result;
        q_spectra = kin_tensor_result.best.S_q;
        t_conc = kin_tensor_result.best.C_at_t;
        collapse_weights = [];
        collapse_label = "first-order kinetic tensor scan"; %#ok<NASGU>
        [~, zeroLagIndex] = min(abs(cfg.delay_s));
        qt_target = delay_tensor(:, :, zeroLagIndex);
        qt_recon = kin_tensor_result.best.qt_recon;

        fprintf("Kinetic tensor scan done: best=%s, rank=%d, weighted relerr=%.6g, fit=%.6f\n", ...
            kin_tensor_result.best.name, kin_tensor_result.rank, kin_tensor_result.relerr, kin_tensor_result.fit);
        fprintf("Best tensor tau/rate parameters(s): ");
        fprintf("%.6g ", kin_tensor_result.best.tau_s);
        fprintf("\n");

        plot_kinetic_tensor_scan_results(q, t, qt_target, qt_recon, kin_tensor_result, cfg);

    case "kinetic_scan"
        [kin_result, cp_reconstructed_tensor] = first_order_kinetic_model_scan( ...
            D_tensor_time, t, U_pc, S_pc, best_par_struct, cfg);

        cp_result = kin_result;
        q_spectra = kin_result.best.S_q;
        t_conc = kin_result.best.C_time;
        qt_target = D_tensor_time;
        qt_recon = kin_result.best.qt_recon;
        collapse_weights = [];
        collapse_label = "first-order kinetic model scan"; %#ok<NASGU>

        fprintf("Kinetic scan done: best=%s, rank=%d, relerr=%.6g, fit=%.6f\n", ...
            kin_result.best.name, kin_result.rank, kin_result.relerr, kin_result.fit);
        fprintf("Best tau/rate parameters(s): ");
        fprintf("%.6g ", kin_result.best.tau_s);
        fprintf("\n");

        plot_kinetic_scan_results(q, t, qt_target, qt_recon, kin_result, cfg);

    case "sequential_varpro"
        [seq_result, cp_reconstructed_tensor] = sequential_varpro_causal_fit( ...
            delay_tensor, valid_weight_time_tau, t, cfg.delay_s, ...
            U_pc, S_pc, D_mean_q, best_par_struct, cfg);

        cp_result = seq_result;
        q_spectra = seq_result.S_q;
        t_conc = seq_result.C_at_t;
        qt_recon = seq_result.qt_recon;
        collapse_weights = [];
        collapse_label = "sequential kinetic model / zero lag";
        [~, zeroLagIndex] = min(abs(cfg.delay_s));
        qt_target = delay_tensor(:, :, zeroLagIndex);

        fprintf("Sequential varpro done: rank=%d, weighted relerr=%.6g, fit=%.6f\n", ...
            cp_result.rank, cp_result.relerr, cp_result.fit);
        fprintf("Sequential tau(s): ");
        fprintf("%.6g ", cp_result.tau_s);
        fprintf("\n");

        plot_sequential_results(q, t, q_spectra, t_conc, qt_target, qt_recon, seq_result, cfg);

    case "shift_mcr"
        [shift_result, cp_reconstructed_tensor] = shift_mcr_causal_als( ...
            delay_tensor, valid_weight_time_tau, t, cfg.delay_s, ...
            U_pc, S_pc, D_mean_q, best_par_struct, cfg);

        cp_result = shift_result;
        q_spectra = shift_result.S_q;
        t_conc = shift_result.C_at_t;
        qt_recon = shift_result.qt_recon;
        collapse_weights = [];
        collapse_label = "zero lag / original time";
        [~, zeroLagIndex] = min(abs(cfg.delay_s));
        qt_target = delay_tensor(:, :, zeroLagIndex);

        fprintf("Shift-MCR done: rank=%d, iterations=%d, weighted relerr=%.6g, fit=%.6f\n", ...
            cp_result.rank, cp_result.num_iters, cp_result.relerr, cp_result.fit);

        %% ===================== 최종 그림 =====================
        plot_shift_mcr_results(q, t, q_spectra, t_conc, qt_target, qt_recon, shift_result, cfg);

    case "cp_als"
        if ~cfg.doCPALS
            error("cfg.doCPALS = true 로 설정하세요.");
        end

        cp_opts = struct();
        cp_opts.maxIters   = cfg.cpMaxIters;
        cp_opts.tol        = cfg.cpTol;
        cp_opts.init       = cfg.cpInit;
        cp_opts.randomSeed = cfg.cpRandomSeed;
        cp_opts.verbose    = cfg.cpVerbose;
        cp_opts.nonnegativeModes = cfg.cpNonnegativeModes;
        cp_opts.nnlsOptions = optimset('Display', 'off', 'MaxIter', 2000);

        [cp_result, cp_reconstructed_tensor] = cp_als_time_nonneg(delay_tensor, valid_weight_time_tau, cfg.cpRank, cp_opts);

        fprintf("CP-ALS done: rank=%d, iterations=%d, final relative error=%.6g, fit=%.6f\n", ...
            cp_result.rank, cp_result.num_iters, cp_result.relerr, cp_result.fit);

        [qt_target, qt_recon, q_spectra, t_conc, collapse_weights, collapse_label] = ...
            collapse_cp_to_qt(delay_tensor, cp_result, cfg); %#ok<ASGLU>

        fprintf("Collapsed q-t matrix mode: %s\n", collapse_label);
        plot_final_results(q, t, q_spectra, t_conc, qt_target, qt_recon, cp_result, collapse_label, cfg);

    otherwise
        error("알 수 없는 cfg.decompMethod: %s", cfg.decompMethod);
end


%% ===================== Local functions =====================
function [res, Xhat] = first_order_kinetic_tensor_scan(X, W_time_tau, t, lag_axis, best_par_struct, cfg)
% Fit the full shifted-time tensor with first-order kinetic candidate models.
% X(q,t,delay) ~= S(q,:) * C(t +/- delay; kinetic parameters)'.
    modelNames = string(cfg.kineticScanModels);
    nModel = numel(modelNames);
    model_results = repmat(empty_kinetic_tensor_model_result(), 1, nModel);

    for iModel = 1:nModel
        model_results(iModel) = fit_first_order_tensor_candidate( ...
            X, W_time_tau, t, lag_axis, modelNames(iModel), best_par_struct, cfg);
    end

    score = kinetic_tensor_selection_scores(model_results, cfg);
    [~, bestIdx] = min(score);
    best = model_results(bestIdx);

    Xhat = best.tensor_recon;
    res = struct();
    res.method = "kinetic_tensor_scan";
    res.rank = size(best.C_at_t, 2);
    res.num_iters = best.num_iters;
    res.relerr = best.relerr;
    res.full_relerr = best.full_relerr;
    res.fit = 1 - best.relerr;
    res.matrix_fit = 1 - best.zero_lag_relerr;
    res.best_start = best.best_start;
    res.constraint = sprintf("first-order kinetic tensor scan; S_q signed; C(%s)>=0; selected by %s", ...
        shifted_time_label_local(cfg), cfg.kineticTensorSelectionMetric);
    res.best = best;
    res.models = model_results;
    res.selection_score = score;

    fprintf("Kinetic tensor scan summary:\n");
    for i = 1:numel(model_results)
        m = model_results(i);
        fprintf("  %-10s tensor_relerr=%.6g zero_relerr=%.6g fit=%.6f dAIC=%.4g overlap=%.4f closure_err=%.3g tau(s)=", ...
            m.name, m.relerr, m.zero_lag_relerr, m.fit, m.aic - min([model_results.aic]), m.overlap, m.closure_max_error);
        fprintf("%.4g ", m.tau_s);
        fprintf("\n");
    end
end

function empty = empty_kinetic_tensor_model_result()
    empty = struct( ...
        "name", "", ...
        "tau_s", [], ...
        "rate_s_inv", [], ...
        "S_q", [], ...
        "C_at_t", [], ...
        "C_ground_at_t", [], ...
        "C_surface", [], ...
        "C_ground_surface", [], ...
        "qt_recon", [], ...
        "tensor_recon", [], ...
        "relerr", inf, ...
        "full_relerr", inf, ...
        "zero_lag_relerr", inf, ...
        "fit", -inf, ...
        "rss", inf, ...
        "aic", inf, ...
        "bic", inf, ...
        "overlap", inf, ...
        "closure_max_error", NaN, ...
        "delay_relerr", [], ...
        "num_iters", 0, ...
        "best_start", 0);
end

function best = fit_first_order_tensor_candidate(X, W_time_tau, t, lag_axis, modelName, best_par_struct, cfg)
    nRate = kinetic_model_num_rates(modelName);
    lb = log10(cfg.kineticScanTauLower_s);
    ub = log10(cfg.kineticScanTauUpper_s);
    starts = make_kinetic_rate_starts(nRate, best_par_struct.tau(:), lb, ub, cfg);

    opts = optimset("Display", "off", ...
        "MaxIter", cfg.kineticScanMaxIter, ...
        "MaxFunEvals", cfg.kineticScanMaxIter * max(30, 12*nRate));

    bestScore = inf;
    best = empty_kinetic_tensor_model_result();

    for iStart = 1:size(starts, 1)
        obj = @(p) kinetic_tensor_candidate_score( ...
            X, W_time_tau, t, lag_axis, modelName, p, lb, ub, best_par_struct, cfg);
        [pOpt, score, ~, output] = fminsearch(obj, starts(iStart, :), opts);
        model = build_kinetic_tensor_candidate( ...
            X, W_time_tau, t, lag_axis, modelName, pOpt, lb, ub, best_par_struct, cfg);

        if score < bestScore
            bestScore = score;
            best = model;
            best.num_iters = output.iterations;
            best.best_start = iStart;
        end
    end
end

function score = kinetic_tensor_candidate_score(X, W_time_tau, t, lag_axis, modelName, p, lb, ub, best_par_struct, cfg)
    pClip = min(max(p(:).', lb), ub);
    boundPenalty = sum((p(:).' - pClip).^2);
    model = build_kinetic_tensor_candidate( ...
        X, W_time_tau, t, lag_axis, modelName, pClip, lb, ub, best_par_struct, cfg);

    scores = kinetic_tensor_selection_scores(model, cfg);
    score = scores(1) + 1e4 * boundPenalty;
end

function scores = kinetic_tensor_selection_scores(models, cfg)
    scores = zeros(1, numel(models));
    for i = 1:numel(models)
        switch lower(string(cfg.kineticTensorSelectionMetric))
            case "aic"
                scores(i) = models(i).aic;
            case "bic"
                scores(i) = models(i).bic;
            case "relerr"
                scores(i) = models(i).relerr;
            case "overlap"
                scores(i) = models(i).overlap;
            otherwise
                error("알 수 없는 cfg.kineticTensorSelectionMetric: %s", cfg.kineticTensorSelectionMetric);
        end
    end
end

function model = build_kinetic_tensor_candidate(X, W_time_tau, t, lag_axis, modelName, p, lb, ub, best_par_struct, cfg)
    p = min(max(p(:).', lb), ub);
    [Z, C_at_t, tau_s, rate_s_inv, Z_ground, C_ground_at_t] = kinetic_surface_matrix(t, lag_axis, modelName, p, best_par_struct, cfg);
    [S_q, XhatMat] = solve_spectra_for_kinetic_tensor(X, W_time_tau, Z);

    [nQ, nT, nLag] = size(X);
    Xhat = reshape(XhatMat, nQ, nT, nLag);

    [~, zeroLagIndex] = min(abs(lag_axis));
    qt_recon = Xhat(:, :, zeroLagIndex);
    qt_target = X(:, :, zeroLagIndex);

    [S_q, C_at_t, Z, Xhat, Z_ground, C_ground_at_t] = sort_kinetic_tensor_components( ...
        S_q, C_at_t, Z, Xhat, t, lag_axis, modelName, p, best_par_struct, cfg);
    qt_recon = Xhat(:, :, zeroLagIndex);
    XhatMat = reshape(Xhat, nQ, nT*nLag);

    w = reshape(W_time_tau, nT*nLag, 1);
    if isempty(w)
        w = ones(nT*nLag, 1);
    end
    rss = sum(sum(((reshape(X, nQ, nT*nLag) - XhatMat).^2) .* w.'));
    normX = sqrt(sum(sum((reshape(X, nQ, nT*nLag).^2) .* w.')));
    relerr = sqrt(rss) / max(normX, eps);
    full_relerr = norm(X(:) - Xhat(:)) / max(norm(X(:)), eps);
    zero_lag_relerr = norm(qt_target(:) - qt_recon(:)) / max(norm(qt_target(:)), eps);

    delay_relerr = zeros(1, nLag);
    for k = 1:nLag
        wk = W_time_tau(:, k).';
        Ek = (X(:, :, k) - Xhat(:, :, k)).^2;
        Nk = X(:, :, k).^2;
        delay_relerr(k) = sqrt(sum(sum(Ek .* wk))) / max(sqrt(sum(sum(Nk .* wk))), eps);
    end

    nObs = nQ * nnz(w > 0);
    nParam = numel(rate_s_inv) + nQ * size(C_at_t, 2);

    model = empty_kinetic_tensor_model_result();
    model.name = char(lower(string(modelName)));
    model.tau_s = tau_s(:).';
    model.rate_s_inv = rate_s_inv(:).';
    model.S_q = S_q;
    model.C_at_t = C_at_t;
    model.C_ground_at_t = C_ground_at_t;
    model.C_surface = Z;
    model.C_ground_surface = Z_ground;
    model.qt_recon = qt_recon;
    model.tensor_recon = Xhat;
    model.relerr = relerr;
    model.full_relerr = full_relerr;
    model.zero_lag_relerr = zero_lag_relerr;
    model.fit = 1 - relerr;
    model.rss = rss;
    model.aic = nObs * log(max(rss / max(nObs, 1), eps)) + 2 * nParam;
    model.bic = nObs * log(max(rss / max(nObs, 1), eps)) + log(max(nObs, 2)) * nParam;
    model.overlap = concentration_overlap_metric(C_at_t);
    if cfg.kineticUseSilentGroundClosure
        model.closure_max_error = max(abs(sum(C_at_t, 2) + C_ground_at_t - 1));
    end
    model.delay_relerr = delay_relerr;
end

function [Z, C_at_t, tau_s, rate_s_inv, Z_ground, C_ground_at_t] = kinetic_surface_matrix(t, lag_axis, modelName, p, best_par_struct, cfg)
    [C_at_t, tau_s, rate_s_inv, C_ground_at_t] = kinetic_concentrations_for_model(t, modelName, p, best_par_struct, cfg);

    nT = numel(t);
    nLag = numel(lag_axis);
    uAll = zeros(nT*nLag, 1);

    for k = 1:nLag
        u = shifted_time_local(t, lag_axis(k), cfg);
        idx = (1:nT) + (k-1)*nT;
        uAll(idx) = u;
    end

    [Z, ~, ~, Z_ground] = kinetic_concentrations_for_model(uAll, modelName, p, best_par_struct, cfg);
    if cfg.zeroSignalBeforeTimeZero
        Z(uAll < cfg.signalTimeZero_s, :) = 0;
        Z_ground(uAll < cfg.signalTimeZero_s, :) = 1;
    end
end

function [C, tau_s, rate_s_inv, C_ground] = kinetic_concentrations_for_model(tQuery, modelName, p, best_par_struct, cfg)
    tau_s = 10.^p(:).';
    modelName = lower(string(modelName));

    switch modelName
        case "parallel"
            if cfg.kineticUseSilentGroundClosure
                error("parallel model은 closure-compatible species model이 아닙니다. sequential 또는 branch를 사용하세요.");
            end
            tau_s = sort(tau_s, "ascend");
            C = kinetic_parallel_concentrations(tQuery, tau_s, best_par_struct, cfg);
            rate_s_inv = 1 ./ tau_s;
            C_ground = max(1 - sum(C, 2), 0);

        case "sequential"
            tau_s = sort(tau_s, "ascend");
            rate_s_inv = 1 ./ tau_s;
            K = [-rate_s_inv(1), 0, 0; ...
                  rate_s_inv(1), -rate_s_inv(2), 0; ...
                  0, rate_s_inv(2), -rate_s_inv(3)];
            c0 = [1; 0; 0];
            C = convolved_first_order_concentrations(tQuery, K, c0, best_par_struct, cfg);
            C_ground = silent_ground_from_active(C, cfg);

        case "branch"
            rate_s_inv = 1 ./ tau_s;
            K = [-(rate_s_inv(1) + rate_s_inv(2)), 0, 0; ...
                  rate_s_inv(1), -rate_s_inv(3), 0; ...
                  rate_s_inv(2), 0, -rate_s_inv(4)];
            c0 = [1; 0; 0];
            C = convolved_first_order_concentrations(tQuery, K, c0, best_par_struct, cfg);
            C_ground = silent_ground_from_active(C, cfg);

        otherwise
            error("알 수 없는 kinetic model: %s", modelName);
    end
end

function [S_q, XhatMat] = solve_spectra_for_kinetic_tensor(X, W_time_tau, Z)
    [nQ, nT, nLag] = size(X);
    Xmat = reshape(X, nQ, nT*nLag);

    if nargin < 2 || isempty(W_time_tau)
        w = ones(nT*nLag, 1);
    else
        w = reshape(W_time_tau, nT*nLag, 1);
    end

    G0 = Z.' * (Z .* w);
    ridge = 1e-8 * max(trace(G0) / max(size(Z, 2), 1), eps);
    G = G0 + ridge * eye(size(Z, 2));
    S_q = ((Xmat .* w.') * Z) / G;
    XhatMat = S_q * Z.';
end

function [S_q, C_at_t, Z, Xhat, Z_ground, C_ground_at_t] = sort_kinetic_tensor_components(S_q, C_at_t, Z, Xhat, t, lag_axis, modelName, p, best_par_struct, cfg)
    [~, peakIdx] = max(C_at_t, [], 1);
    [~, order] = sort(peakIdx, "ascend");
    S_q = S_q(:, order);
    C_at_t = C_at_t(:, order);
    Z = Z(:, order);

    [nQ, nT, nLag] = size(Xhat);
    XhatMat = S_q * Z.';
    Xhat = reshape(XhatMat, nQ, nT, nLag);

    % Recompute the shifted surface after sorting so plots and saved results
    % use the same component order as C_at_t.
    [Z0, ~, ~, ~, Z_ground, C_ground_at_t] = kinetic_surface_matrix(t, lag_axis, modelName, p, best_par_struct, cfg);
    Z = Z0(:, order);
end

function plot_kinetic_tensor_scan_results(q, t, qt_target, qt_recon, kin_tensor_result, cfg)
    best = kin_tensor_result.best;
    modelNames = string({kin_tensor_result.models.name});
    relerrs = [kin_tensor_result.models.relerr];
    aics = [kin_tensor_result.models.aic];

    figure("Name", "Kinetic tensor model scan scores");
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    bar(relerrs);
    grid on;
    set(gca, "XTick", 1:numel(modelNames), "XTickLabel", modelNames);
    ylabel("weighted tensor relative error");
    title("Full tensor fit error");
    nexttile;
    bar(aics - min(aics));
    grid on;
    set(gca, "XTick", 1:numel(modelNames), "XTickLabel", modelNames);
    ylabel("\DeltaAIC");
    title("Model selection penalty");

    figure("Name", "Kinetic tensor target vs zero-lag reconstruction");
    clim_all = [min([qt_target(:); qt_recon(:)]), max([qt_target(:); qt_recon(:)])];
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    imagesc(t, q, qt_target);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title("Tensor zero-lag target");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile;
    imagesc(t, q, qt_recon);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title(sprintf("%s tensor model (fit %.4f)", best.name, best.fit));
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Kinetic tensor concentration profiles");
    if cfg.kineticPlotSilentGround && ~isempty(best.C_ground_at_t)
        plot(t, [best.C_at_t, best.C_ground_at_t], "LineWidth", 1.5);
        legendLabels = [compose("component %d", 1:size(best.C_at_t, 2)), "silent ground"];
    else
        plot(t, best.C_at_t, "LineWidth", 1.5);
        legendLabels = compose("component %d", 1:size(best.C_at_t, 2));
    end
    grid on;
    xlabel("Time (s)");
    ylabel("Concentration");
    title(sprintf("Best tensor-constrained model: %s", best.name));
    legend(legendLabels, "Location", "best");
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Kinetic tensor q spectra");
    plot(q, best.S_q, "LineWidth", 1.5);
    hold on; yline(0, "k--", "LineWidth", 1.0); hold off;
    grid on;
    xlabel("q");
    ylabel("Species spectrum (a.u.)");
    title(sprintf("%s tensor-constrained spectra", best.name));
    legend(compose("component %d", 1:size(best.S_q, 2)), "Location", "best");

    figure("Name", "Kinetic tensor delay-wise residual");
    plot(cfg.delay_s, best.delay_relerr, "o-", "LineWidth", 1.3);
    grid on;
    xlabel("Delay (s)");
    ylabel("relative error per delay slice");
    title("Full tensor residual by delay");
    if all(cfg.delay_s > 0)
        set(gca, "XScale", "log");
    elseif any(cfg.delay_s > 0)
        set(gca, "XScale", "linear");
    end

    fprintf("Kinetic tensor selected model: %s\n", best.name);
    if strcmpi(best.name, "branch")
        fprintf("Branch parameter order: tau12, tau13, tau2->ground, tau3->ground.\n");
    elseif strcmpi(best.name, "sequential")
        fprintf("Sequential parameter order: tau1->2, tau2->3, tau3->ground.\n");
    end
end

function [res, Xhat] = first_order_kinetic_model_scan(Y, t, U_pc, S_pc, best_par_struct, cfg)
% Compare topology-neutral first-order candidate models by variable projection.
% Spectra are solved linearly for each concentration model; only kinetic rates
% are optimized nonlinearly.
    modelNames = string(cfg.kineticScanModels);
    nModel = numel(modelNames);
    model_results = repmat(empty_kinetic_model_result(), 1, nModel);

    for iModel = 1:nModel
        model_results(iModel) = fit_first_order_candidate(Y, t, modelNames(iModel), best_par_struct, cfg);
    end

    [~, bestIdx] = min([model_results.aic]);
    best = model_results(bestIdx);
    direct = direct_rsv_exponential_components(U_pc, S_pc, best_par_struct, t);

    Xhat = best.qt_recon;
    res = struct();
    res.method = "kinetic_scan";
    res.rank = size(best.C_time, 2);
    res.num_iters = best.num_iters;
    res.relerr = best.relerr;
    res.full_relerr = best.relerr;
    res.fit = 1 - best.relerr;
    res.matrix_fit = 1 - best.relerr;
    res.best_start = best.best_start;
    res.constraint = "first-order model scan; S_q signed; C(t)>=0; selected by AIC";
    res.best = best;
    res.models = model_results;
    res.direct_rsv_exp = direct;
end

function empty = empty_kinetic_model_result()
    empty = struct( ...
        "name", "", ...
        "tau_s", [], ...
        "rate_s_inv", [], ...
        "S_q", [], ...
        "C_time", [], ...
        "qt_recon", [], ...
        "relerr", inf, ...
        "fit", -inf, ...
        "rss", inf, ...
        "aic", inf, ...
        "bic", inf, ...
        "overlap", inf, ...
        "num_iters", 0, ...
        "best_start", 0);
end

function best = fit_first_order_candidate(Y, t, modelName, best_par_struct, cfg)
    nRate = kinetic_model_num_rates(modelName);
    lb = log10(cfg.kineticScanTauLower_s);
    ub = log10(cfg.kineticScanTauUpper_s);
    starts = make_kinetic_rate_starts(nRate, best_par_struct.tau(:), lb, ub, cfg);

    opts = optimset("Display", "off", ...
        "MaxIter", cfg.kineticScanMaxIter, ...
        "MaxFunEvals", cfg.kineticScanMaxIter * max(30, 12*nRate));

    bestScore = inf;
    best = empty_kinetic_model_result();

    for iStart = 1:size(starts, 1)
        obj = @(p) kinetic_candidate_score(Y, t, modelName, p, lb, ub, best_par_struct, cfg);
        [pOpt, score, exitflag, output] = fminsearch(obj, starts(iStart, :), opts); %#ok<ASGLU>
        model = build_kinetic_candidate(Y, t, modelName, pOpt, lb, ub, best_par_struct, cfg);

        if score < bestScore
            bestScore = score;
            best = model;
            best.num_iters = output.iterations;
            best.best_start = iStart;
        end
    end
end

function nRate = kinetic_model_num_rates(modelName)
    switch lower(string(modelName))
        case {"parallel", "sequential"}
            nRate = 3;
        case "branch"
            nRate = 4;
        otherwise
            error("알 수 없는 kinetic model: %s", modelName);
    end
end

function starts = make_kinetic_rate_starts(nRate, fittedTau, lb, ub, cfg)
    rng(cfg.cpRandomSeed + 307 + nRate);
    nStart = max(1, cfg.kineticScanNumStarts);
    starts = lb + (ub - lb) * rand(nStart, nRate);

    tau = sort(fittedTau(:).');
    if isempty(tau)
        tau = logspace(lb, ub, nRate);
    end
    while numel(tau) < nRate
        tau(end+1) = min(10^ub, tau(end) * 3); %#ok<AGROW>
    end

    if nRate == 4 && numel(fittedTau) >= 3
        baseTau = sort(fittedTau(:).');
        baseTau = [baseTau(1), sqrt(baseTau(1)*baseTau(2)), baseTau(2), baseTau(3)];
    else
        baseTau = tau(1:nRate);
    end
    starts(1, :) = min(max(log10(baseTau), lb), ub);

    for i = 2:min(nStart, 1+nRate)
        jitter = 0.35 * randn(1, nRate);
        starts(i, :) = min(max(starts(1, :) + jitter, lb), ub);
    end
end

function score = kinetic_candidate_score(Y, t, modelName, p, lb, ub, best_par_struct, cfg)
    pClip = min(max(p(:).', lb), ub);
    boundPenalty = sum((p(:).' - pClip).^2);
    model = build_kinetic_candidate(Y, t, modelName, pClip, lb, ub, best_par_struct, cfg);
    score = model.aic + 1e4 * boundPenalty;
end

function model = build_kinetic_candidate(Y, t, modelName, p, lb, ub, best_par_struct, cfg) %#ok<INUSD>
    p = min(max(p(:).', lb), ub);
    tau_s = 10.^p;
    modelName = lower(string(modelName));

    switch modelName
        case "parallel"
            tau_s = sort(tau_s, "ascend");
            C = kinetic_parallel_concentrations(t, tau_s, best_par_struct, cfg);
            rate_s_inv = 1 ./ tau_s;

        case "sequential"
            tau_s = sort(tau_s, "ascend");
            rate_s_inv = 1 ./ tau_s;
            K = [-rate_s_inv(1), 0, 0; ...
                  rate_s_inv(1), -rate_s_inv(2), 0; ...
                  0, rate_s_inv(2), -rate_s_inv(3)];
            c0 = [1; 0; 0];
            C = convolved_first_order_concentrations(t, K, c0, best_par_struct, cfg);

        case "branch"
            rate_s_inv = 1 ./ tau_s;
            K = [-(rate_s_inv(1) + rate_s_inv(2)), 0, 0; ...
                  rate_s_inv(1), -rate_s_inv(3), 0; ...
                  rate_s_inv(2), 0, -rate_s_inv(4)];
            c0 = [1; 0; 0];
            C = convolved_first_order_concentrations(t, K, c0, best_par_struct, cfg);

        otherwise
            error("알 수 없는 kinetic model: %s", modelName);
    end

    [S_q, C, qt_recon] = solve_spectra_for_kinetic_concentrations(Y, C);
    [S_q, C] = sort_kinetic_components_by_peak(S_q, C);
    qt_recon = S_q * C.';

    err = Y - qt_recon;
    rss = sum(err(:).^2);
    relerr = sqrt(rss) / max(norm(Y(:)), eps);
    nObs = numel(Y);
    nParam = numel(rate_s_inv) + size(Y, 1) * size(C, 2);

    model = empty_kinetic_model_result();
    model.name = char(modelName);
    model.tau_s = tau_s(:).';
    model.rate_s_inv = rate_s_inv(:).';
    model.S_q = S_q;
    model.C_time = C;
    model.qt_recon = qt_recon;
    model.relerr = relerr;
    model.fit = 1 - relerr;
    model.rss = rss;
    model.aic = nObs * log(max(rss / nObs, eps)) + 2 * nParam;
    model.bic = nObs * log(max(rss / nObs, eps)) + log(nObs) * nParam;
    model.overlap = concentration_overlap_metric(C);
end

function C = kinetic_parallel_concentrations(t, tau_s, best_par_struct, cfg)
    t = t(:);
    C = zeros(numel(t), numel(tau_s));
    for k = 1:numel(tau_s)
        if cfg.kineticUseIRFConvolution
            C(:, k) = conv_gauss_expdec_local(t, 1, 0, tau_s(k), best_par_struct.t0, best_par_struct.FWHM);
        else
            C(:, k) = exp(-max(t - best_par_struct.t0, 0) ./ tau_s(k));
        end
    end
    C = normalize_kinetic_concentrations(C);
end

function C = convolved_first_order_concentrations(t, K, c0, best_par_struct, cfg)
    t = t(:);
    R = numel(c0);
    C = zeros(numel(t), R);

    if cfg.kineticUseIRFConvolution
        [V, L] = eig(K);
        coeff = V \ c0;
        for m = 1:R
            lambda = real(L(m, m));
            if lambda >= -eps
                basis = double(t >= best_par_struct.t0);
            else
                tau = -1 / lambda;
                basis = conv_gauss_expdec_local(t, 1, 0, tau, best_par_struct.t0, best_par_struct.FWHM);
            end
            amp = real(V(:, m) * coeff(m));
            C = C + basis * amp.';
        end
    else
        for i = 1:numel(t)
            u = max(t(i) - best_par_struct.t0, 0);
            C(i, :) = (expm(K * u) * c0).';
        end
    end

    C(abs(C) < 1e-12) = 0;
    C = max(real(C), 0);
    if ~cfg.kineticUseSilentGroundClosure
        C = normalize_kinetic_concentrations(C);
    end
    C = enforce_active_closure(C, cfg);
end

function C = enforce_active_closure(C, cfg)
    if ~cfg.kineticUseSilentGroundClosure
        return;
    end

    C = max(real(C), 0);
    activeSum = sum(C, 2);
    over = activeSum > 1;
    if any(over)
        C(over, :) = C(over, :) ./ activeSum(over);
    end
end

function C_ground = silent_ground_from_active(C, cfg)
    if cfg.kineticUseSilentGroundClosure
        C_ground = max(1 - sum(C, 2), 0);
    else
        C_ground = zeros(size(C, 1), 1);
    end
end

function C = normalize_kinetic_concentrations(C)
    C = max(real(C), 0);
    for r = 1:size(C, 2)
        scale = max(C(:, r));
        if scale > 0
            C(:, r) = C(:, r) / scale;
        end
    end
end

function [S_q, C, qt_recon] = solve_spectra_for_kinetic_concentrations(Y, C)
    good = max(C, [], 1) > 1e-10;
    C = C(:, good);
    G = C.' * C + 1e-12 * eye(size(C, 2));
    S_q = (Y * C) / G;
    qt_recon = S_q * C.';
end

function [S_q, C] = sort_kinetic_components_by_peak(S_q, C)
    [~, idx] = max(C, [], 1);
    [~, order] = sort(idx, "ascend");
    C = C(:, order);
    S_q = S_q(:, order);
end

function overlap = concentration_overlap_metric(C)
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

function direct = direct_rsv_exponential_components(U_pc, S_pc, best_par_struct, t)
    tau = best_par_struct.tau(:).';
    basis = zeros(numel(t), numel(tau));
    for k = 1:numel(tau)
        basis(:, k) = conv_gauss_expdec_local(t, 1, 0, tau(k), best_par_struct.t0, best_par_struct.FWHM);
    end
    spectra = (U_pc * S_pc) * best_par_struct.A;
    offset_spectrum = (U_pc * S_pc) * best_par_struct.A0(:);
    direct = struct();
    direct.tau_s = tau;
    direct.basis = basis;
    direct.spectra = spectra;
    direct.offset_spectrum = offset_spectrum;
    direct.reconstruction = spectra * basis.' + offset_spectrum;
end

function plot_kinetic_scan_results(q, t, qt_target, qt_recon, kin_result, cfg)
    best = kin_result.best;
    modelNames = string({kin_result.models.name});
    relerrs = [kin_result.models.relerr];
    aics = [kin_result.models.aic];

    figure("Name", "Kinetic model scan scores");
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    bar(relerrs);
    grid on;
    set(gca, "XTick", 1:numel(modelNames), "XTickLabel", modelNames);
    ylabel("relative error");
    title("Variable-projection fit error");
    nexttile;
    bar(aics - min(aics));
    grid on;
    set(gca, "XTick", 1:numel(modelNames), "XTickLabel", modelNames);
    ylabel("\DeltaAIC");
    title("Model selection penalty");

    figure("Name", "Kinetic target vs reconstructed q-t maps");
    clim_all = [min([qt_target(:); qt_recon(:)]), max([qt_target(:); qt_recon(:)])];
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    imagesc(t, q, qt_target);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title("RSV-fitted q-t target");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile;
    imagesc(t, q, qt_recon);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title(sprintf("%s reconstruction (fit %.4f)", best.name, best.fit));
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Kinetic concentration profiles");
    plot(t, best.C_time, "LineWidth", 1.5);
    grid on;
    xlabel("Time (s)");
    ylabel("Concentration (normalized)");
    title(sprintf("Best first-order model: %s", best.name));
    legend(compose("component %d", 1:size(best.C_time, 2)), "Location", "best");
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Kinetic q spectra");
    plot(q, best.S_q, "LineWidth", 1.5);
    hold on; yline(0, "k--", "LineWidth", 1.0); hold off;
    grid on;
    xlabel("q");
    ylabel("Species spectrum (a.u.)");
    title(sprintf("%s model spectra", best.name));
    legend(compose("component %d", 1:size(best.S_q, 2)), "Location", "best");

    direct = kin_result.direct_rsv_exp;
    figure("Name", "Direct RSV exponential DAS");
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    plot(t, direct.basis, "LineWidth", 1.3);
    grid on;
    xlabel("Time (s)");
    ylabel("IRF-convolved exponential basis");
    title("Direct fitted exponential basis");
    legend(compose("\\tau %.3g ps", direct.tau_s * 1e12), "Location", "best");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile;
    plot(q, direct.spectra, "LineWidth", 1.3);
    hold on; yline(0, "k--", "LineWidth", 1.0); hold off;
    grid on;
    xlabel("q");
    ylabel("DAS amplitude (a.u.)");
    title("Direct lifetime-associated spectra");
    legend(compose("\\tau %.3g ps", direct.tau_s * 1e12), "Location", "best");

    fprintf("Kinetic model scan summary:\n");
    for i = 1:numel(kin_result.models)
        m = kin_result.models(i);
        fprintf("  %-10s relerr=%.6g fit=%.6f dAIC=%.4g overlap=%.4f tau(s)=", ...
            m.name, m.relerr, m.fit, m.aic - min(aics), m.overlap);
        fprintf("%.4g ", m.tau_s);
        fprintf("\n");
    end

    if strcmpi(best.name, "branch")
        fprintf("Branch parameter order: tau12, tau13, tau2->ground, tau3->ground.\n");
    elseif strcmpi(best.name, "sequential")
        fprintf("Sequential parameter order: tau1->2, tau2->3, tau3->ground.\n");
    end

    fprintf("Direct RSV fitted tau(s): ");
    fprintf("%.6g ", sort(direct.tau_s));
    fprintf("\n");

    if cfg.kineticUseIRFConvolution
        fprintf("Kinetic concentrations include Gaussian-IRF convolution from fit_rsv.\n");
    end
end

function M = read_numeric_table(fileName)
% 공백/콤마/탭이 섞인 numeric text 파일을 행 구조 유지하며 읽는다.
    if exist(fileName, "file") ~= 2
        error("파일을 찾을 수 없습니다: %s", fileName);
    end

    txt = fileread(fileName);
    lines = regexp(txt, "\r\n|\n|\r", "split");
    rows = {};
    maxCols = 0;

    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end
        line = regexprep(line, "[,;\t]", " ");
        vals = sscanf(line, "%f").';
        if isempty(vals)
            continue;
        end
        rows{end+1, 1} = vals; %#ok<AGROW>
        maxCols = max(maxCols, numel(vals));
    end

    if isempty(rows)
        error("숫자 데이터를 읽지 못했습니다: %s", fileName);
    end

    nRows = numel(rows);
    M = nan(nRows, maxCols);
    for i = 1:nRows
        vals = rows{i};
        M(i, 1:numel(vals)) = vals;
    end

    if any(isnan(M(:)))
        error("행마다 열 개수가 다릅니다. 파일 형식을 확인하세요: %s", fileName);
    end
end

function V_model = eval_fitted_rsv(parStruct, tQuery)
% fit_rsv가 반환한 best_par_struct를 이용해 임의 시간 tQuery에서 RSV 모델 평가.
    tQuery = tQuery(:);

    A    = parStruct.A;
    tau  = parStruct.tau(:);
    A0   = parStruct.A0(:);
    FWHM = parStruct.FWHM;
    t0   = parStruct.t0;

    numRSV = size(A, 1);
    V_model = zeros(numel(tQuery), numRSV);

    for j = 1:numRSV
        V_model(:, j) = conv_gauss_expdec_local(tQuery, A(j, :), A0(j), tau, t0, FWHM);
    end
end

function y = conv_gauss_expdec_local(t, A, A0, tau, t0, FWHM)
% fit_rsv.m 내부 모델과 동일한 Gaussian-IRF convolved multi-exponential model.
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

function [qt_target, qt_recon, q_spectra, t_conc, w_delay, labelText] = collapse_cp_to_qt(delay_tensor, cp, cfg)
% CP 결과를 다시 q-t 2차원 데이터로 줄인다.
% q_spectra: q축 intermediate spectra (음수 가능)
% t_conc   : t축 concentration profiles (비음수)
% qt_recon : q x t 재구성 행렬

    lambda_row = reshape(cp.lambda(:), 1, []);
    C = cp.C_delay;

    switch string(cfg.qtCollapseMode)
        case "mean_delay"
            w_delay = mean(C, 1);
            qt_target = mean(delay_tensor, 3);
            labelText = "mean over delay";

        case "sum_delay"
            w_delay = sum(C, 1);
            qt_target = sum(delay_tensor, 3);
            labelText = "sum over delay";

        case "delay_index"
            idx = cfg.plotDelayIndex;
            if idx < 1 || idx > size(C,1)
                error("cfg.plotDelayIndex가 delay 범위를 벗어났습니다.");
            end
            w_delay = C(idx, :);
            qt_target = delay_tensor(:, :, idx);
            labelText = sprintf("delay slice index %d", idx);

        otherwise
            error("알 수 없는 cfg.qtCollapseMode: %s", cfg.qtCollapseMode);
    end

    % 선택한 delay-collapse 가중치를 q-mode에 흡수
    q_spectra = cp.A_q .* reshape(lambda_row .* w_delay, 1, []);

    % 시간축 농도 profile은 nonnegative factor를 그대로 사용
    t_conc = cp.B_time;

    % q-t 2D 재구성
    qt_recon = q_spectra * t_conc.';
end

function plot_final_results(q, t, q_spectra, t_conc, qt_target, qt_recon, cp_result, collapse_label, cfg)
    nComp = size(q_spectra, 2);

    % 1) q축 intermediate spectra
    figure("Name", "Intermediate spectra on q-axis");
    plot(q, q_spectra, "LineWidth", 1.5);
    hold on;
    yline(0, "k--", "LineWidth", 1.0);
    hold off;
    grid on;
    xlabel("q");
    ylabel("Intermediate spectrum (a.u.)");
    title(sprintf("CP intermediate spectra on q (rank = %d)", cp_result.rank));
    legend(compose("component %d", 1:nComp), "Location", "best");

    % 2) t축 concentration profiles (nonnegative)
    figure("Name", "Concentration profiles on time-axis");
    plot(t, t_conc, "LineWidth", 1.5);
    grid on;
    xlabel("Time (s)");
    ylabel("Concentration profile (a.u.)");
    title("CP time-mode concentration profiles (nonnegative)");
    legend(compose("component %d", 1:nComp), "Location", "best");
    if all(t > 0)
        set(gca, "XScale", "log");
    end

    % 3) 재구성된 q-t 2D map
    figure("Name", "Reconstructed q-t map from CP");
    imagesc(t, q, qt_recon);
    axis xy;
    colorbar;
    xlabel("Time (s)");
    ylabel("q");
    title(sprintf("CP reconstructed q-t map (%s)", collapse_label));
    if all(t > 0)
        set(gca, "XScale", "log");
    end

    % 4) 원본 collapse와 비교
    if cfg.showComparisonMap
        clim_all = [min([qt_target(:); qt_recon(:)]), max([qt_target(:); qt_recon(:)])];
        figure("Name", "Target vs CP reconstructed q-t maps");
        tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");

        nexttile;
        imagesc(t, q, qt_target);
        axis xy;
        colorbar;
        caxis(clim_all);
        xlabel("Time (s)");
        ylabel("q");
        title(sprintf("Target q-t map (%s)", collapse_label));
        if all(t > 0)
            set(gca, "XScale", "log");
        end

        nexttile;
        imagesc(t, q, qt_recon);
        axis xy;
        colorbar;
        caxis(clim_all);
        xlabel("Time (s)");
        ylabel("q");
        title(sprintf("CP reconstructed q-t map (%s)", collapse_label));
        if all(t > 0)
            set(gca, "XScale", "log");
        end
    end
end

function [res, Xhat] = sequential_varpro_causal_fit(X, W_time_tau, t, lag_axis, U_pc, S_pc, D_mean_q, best_par_struct, cfg)
% Fit X(q,t,lag) with a causal first-order sequential model:
% A1 -> A2 -> ... -> AR -> ground, and X(q,t,lag)=S(q,:)C(t+lag)'.
    R = cfg.cpRank;
    u_axis = make_shift_mcr_time_grid(t, lag_axis, cfg);
    Y = eval_reconstructed_matrix_at_time(U_pc, S_pc, D_mean_q, best_par_struct, u_axis);

    tau0 = sort(best_par_struct.tau(:), "ascend");
    if numel(tau0) < R
        tau0 = [tau0; logspace(log10(min(t(t > 0))), log10(max(t)), R-numel(tau0)).'];
    end
    tau0 = tau0(1:R);
    tau0 = min(max(tau0, cfg.sequentialTauLower_s), cfg.sequentialTauUpper_s);

    if cfg.sequentialFitRates
        [tau_s, S_q, C_time, matrixRelerr, bestStart] = fit_sequential_rates_varpro(Y, u_axis, tau0, cfg);
    else
        tau_s = tau0(:);
        C_time = sequential_concentration_local(u_axis, tau_s);
        S_q = solve_spectra_for_concentration(Y, C_time);
        matrixRelerr = norm(Y(:) - reshape(S_q * C_time.', [], 1)) / norm(Y(:));
        bestStart = 1;
    end

    [S_q, C_time] = normalize_shift_mcr_columns(S_q, C_time);
    Xhat = reconstruct_shift_mcr_tensor(S_q, C_time, u_axis, t, lag_axis, cfg);
    normXw = weighted_tensor_norm_local(X, W_time_tau);
    tensorRelerr = weighted_tensor_relerr_local(X, Xhat, W_time_tau, normXw);
    fullRelerr = norm(X(:) - Xhat(:)) / norm(X(:));

    C_at_t = interpolate_shift_mcr_concentration(u_axis, C_time, t);
    qt_recon = S_q * C_at_t.';

    res = struct();
    res.rank = R;
    res.S_q = S_q;
    res.C_time = C_time;
    res.u_axis = u_axis;
    res.C_at_t = C_at_t;
    res.qt_recon = qt_recon;
    res.tau_s = tau_s(:);
    res.rate_s_inv = 1 ./ tau_s(:);
    res.fit = 1 - tensorRelerr;
    res.relerr = tensorRelerr;
    res.full_relerr = fullRelerr;
    res.matrix_fit = 1 - matrixRelerr;
    res.matrix_relerr = matrixRelerr;
    res.fit_history = [];
    res.relerr_history = [];
    res.num_iters = 0;
    res.best_start = bestStart;
    res.weight_valid_fraction = nnz(W_time_tau > 0) / numel(W_time_tau);
    res.constraint = sprintf("sequential varpro: A1->...->AR->ground; S_q signed; C_seq(%s)>=0; C(u<0)=0", ...
        shifted_time_label_local(cfg));
end

function [tauBest, SBest, CBest, relBest, bestStart] = fit_sequential_rates_varpro(Y, u_axis, tau0, cfg)
    R = numel(tau0);
    lb = log10(cfg.sequentialTauLower_s);
    ub = log10(cfg.sequentialTauUpper_s);
    p0 = log10(tau0(:));
    relBest = inf;
    tauBest = tau0(:);
    SBest = [];
    CBest = [];
    bestStart = 1;

    opts = optimset('Display', 'off', 'MaxIter', 500, 'MaxFunEvals', 5000, 'TolX', 1e-8, 'TolFun', 1e-8);
    for iStart = 1:cfg.sequentialNumStarts
        if iStart == 1
            pStart = p0;
        else
            pStart = p0 + 0.6 * randn(R, 1);
            pStart = min(max(pStart, lb), ub);
        end

        obj = @(p) sequential_varpro_objective(p, Y, u_axis, lb, ub);
        pFit = fminsearch(obj, pStart, opts);
        [relerr, S_q, C_time, tau_s] = sequential_varpro_objective(pFit, Y, u_axis, lb, ub);

        if relerr < relBest
            relBest = relerr;
            tauBest = tau_s;
            SBest = S_q;
            CBest = C_time;
            bestStart = iStart;
        end
    end
end

function [relerr, S_q, C_time, tau_s] = sequential_varpro_objective(p, Y, u_axis, lb, ub)
    p = p(:);
    penalty = sum(max(lb - p, 0).^2 + max(p - ub, 0).^2);
    p = min(max(p, lb), ub);
    tau_s = sort(10.^p, "ascend");
    C_time = sequential_concentration_local(u_axis, tau_s);
    S_q = solve_spectra_for_concentration(Y, C_time);
    Yhat = S_q * C_time.';
    relerr = norm(Y(:) - Yhat(:)) / norm(Y(:)) + 10 * penalty;
end

function S_q = solve_spectra_for_concentration(Y, C_time)
    G = C_time.' * C_time + 1e-12 * eye(size(C_time, 2));
    S_q = (Y * C_time) / G;
end

function C = sequential_concentration_local(u_axis, tau_s)
    u_axis = u_axis(:);
    tau_s = tau_s(:);
    R = numel(tau_s);
    k = 1 ./ tau_s;

    K = zeros(R, R);
    for r = 1:R
        K(r, r) = -k(r);
        if r < R
            K(r+1, r) = k(r);
        end
    end

    c0 = zeros(R, 1);
    c0(1) = 1;
    C = zeros(numel(u_axis), R);
    for i = 1:numel(u_axis)
        u = u_axis(i);
        if u < 0
            continue;
        end
        C(i, :) = (expm(K * u) * c0).';
    end
    C = max(C, 0);
end

function plot_sequential_results(q, t, q_spectra, t_conc, qt_target, qt_recon, seq_result, cfg)
    nComp = size(q_spectra, 2);
    uLabel = shifted_time_label_local(cfg);

    figure("Name", "Sequential concentration profiles");
    plot(seq_result.u_axis, seq_result.C_time, "LineWidth", 1.5);
    grid on;
    xlabel("Shifted model time u = " + uLabel + " (s)");
    ylabel("Concentration");
    title("First-order sequential concentration profiles");
    legend(compose("component %d", 1:nComp), "Location", "best");
    if all(seq_result.u_axis > 0), set(gca, "XScale", "log"); end

    figure("Name", "Sequential concentrations at measured time");
    plot(t, t_conc, "LineWidth", 1.5);
    grid on;
    xlabel("Measured time t (s)");
    ylabel("Concentration");
    title("Sequential C_r(t) on measured time grid");
    legend(compose("component %d", 1:nComp), "Location", "best");
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Sequential q spectra");
    plot(q, q_spectra, "LineWidth", 1.5);
    hold on; yline(0, "k--", "LineWidth", 1.0); hold off;
    grid on;
    xlabel("q");
    ylabel("Species spectrum (a.u.)");
    title(sprintf("Sequential species spectra, tau = [%s] ps", sprintf("%.3g ", seq_result.tau_s * 1e12)));
    legend(compose("component %d", 1:nComp), "Location", "best");

    figure("Name", "Sequential target vs reconstructed q-t maps");
    clim_all = [min([qt_target(:); qt_recon(:)]), max([qt_target(:); qt_recon(:)])];
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    imagesc(t, q, qt_target);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title("Target q-t map at zero lag");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile;
    imagesc(t, q, qt_recon);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title(sprintf("Sequential reconstruction (fit %.4f)", seq_result.fit));
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Sequential causal concentration surfaces");
    tiledlayout(nComp, 1, "TileSpacing", "compact");
    for r = 1:nComp
        surf_r = zeros(numel(t), numel(cfg.delay_s));
        for k = 1:numel(cfg.delay_s)
            u = shifted_time_local(t, cfg.delay_s(k), cfg);
            c_eval = interpolate_shift_mcr_concentration(seq_result.u_axis, seq_result.C_time(:, r), u);
            if cfg.zeroSignalBeforeTimeZero
                c_eval(u < cfg.signalTimeZero_s) = 0;
            end
            surf_r(:, k) = c_eval;
        end
        nexttile;
        imagesc(1:numel(cfg.delay_s), t, surf_r);
        axis xy; colorbar;
        xlabel("Signed lag index");
        ylabel("Time (s)");
        title(sprintf("component %d: C_r(%s)", r, uLabel));
        if all(t > 0), set(gca, "YScale", "log"); end
    end
end

function [res, Xhat] = shift_mcr_causal_als(X, W_time_tau, t, lag_axis, U_pc, S_pc, D_mean_q, best_par_struct, cfg)
% Decompose X(q,t,lag) with the same shifted-time convention used to build X.
% X(q,t,lag) ~= sum_r S_q(q,r) * C_r(u), where u is t-delay or t+delay.
    R = cfg.cpRank;
    u_axis = make_shift_mcr_time_grid(t, lag_axis, cfg);
    Y = eval_reconstructed_matrix_at_time(U_pc, S_pc, D_mean_q, best_par_struct, u_axis);

    nnlsOptions = optimset('Display', 'off', 'MaxIter', 2000);
    maxIters = cfg.cpMaxIters;
    tol = cfg.cpTol;
    nStarts = cfg.shiftMcrNumStarts;
    rng(cfg.cpRandomSeed);

    [nQ, nU] = size(Y);
    normY = norm(Y(:));
    bestRelerr = inf;
    best = struct();

    for iStart = 1:nStarts
        C_time = initialize_shift_mcr_concentration(Y, R, iStart);
        prevFit = -Inf;
        fitHistory = nan(maxIters, 1);
        relerrHistory = nan(maxIters, 1);

        for iter = 1:maxIters
            G = C_time.' * C_time + 1e-12 * eye(R);
            S_q = (Y * C_time) / G;

            if cfg.shiftMcrUseNNLS
                for iu = 1:nU
                    C_time(iu, :) = lsqnonneg(S_q, Y(:, iu), nnlsOptions).';
                end
            else
                Gc = S_q.' * S_q + 1e-12 * eye(R);
                C_time = (Y.' * S_q) / Gc;
                C_time = max(C_time, 0);
            end

            C_time = apply_shift_mcr_shape_constraint(C_time, cfg);

            [S_q, C_time] = normalize_shift_mcr_columns(S_q, C_time);

            Yhat = S_q * C_time.';
            relerr = norm(Y(:) - Yhat(:)) / normY;
            fit = 1 - relerr;
            fitHistory(iter) = fit;
            relerrHistory(iter) = relerr;

            if iter > 1 && abs(fit - prevFit) < tol
                fitHistory = fitHistory(1:iter);
                relerrHistory = relerrHistory(1:iter);
                break;
            end
            prevFit = fit;

            if iter == maxIters
                fitHistory = fitHistory(1:iter);
                relerrHistory = relerrHistory(1:iter);
            end
        end

        if relerr < bestRelerr
            bestRelerr = relerr;
            best.S_q = S_q;
            best.C_time = C_time;
            best.fitHistory = fitHistory;
            best.relerrHistory = relerrHistory;
            best.num_iters = numel(fitHistory);
            best.start = iStart;
        end
    end

    [best.S_q, best.C_time] = sort_shift_mcr_components(best.S_q, best.C_time);
    Xhat = reconstruct_shift_mcr_tensor(best.S_q, best.C_time, u_axis, t, lag_axis, cfg);

    normXw = weighted_tensor_norm_local(X, W_time_tau);
    tensorRelerr = weighted_tensor_relerr_local(X, Xhat, W_time_tau, normXw);
    fullRelerr = norm(X(:) - Xhat(:)) / norm(X(:));

    C_at_t = interpolate_shift_mcr_concentration(u_axis, best.C_time, t);
    qt_recon = best.S_q * C_at_t.';

    res = struct();
    res.rank = R;
    res.S_q = best.S_q;
    res.C_time = best.C_time;
    res.u_axis = u_axis;
    res.C_at_t = C_at_t;
    res.qt_recon = qt_recon;
    res.fit = 1 - tensorRelerr;
    res.relerr = tensorRelerr;
    res.full_relerr = fullRelerr;
    res.matrix_fit = 1 - bestRelerr;
    res.matrix_relerr = bestRelerr;
    res.fit_history = best.fitHistory;
    res.relerr_history = best.relerrHistory;
    res.num_iters = best.num_iters;
    res.best_start = best.start;
    res.weight_valid_fraction = nnz(W_time_tau > 0) / numel(W_time_tau);
    res.constraint = sprintf("shift MCR: A_q signed; C_r(%s) >= 0; zero-before-time0=%d", ...
        shifted_time_label_local(cfg), cfg.zeroSignalBeforeTimeZero);
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

function label = shifted_time_label_local(cfg)
    switch cfg.shiftMode
        case "delay_signal"
            label = "t-delay";
        case "evaluate_later_time"
            label = "t+delay";
        otherwise
            label = "shifted time";
    end
end

function u_axis = make_shift_mcr_time_grid(t, lag_axis, cfg)
    t = t(:);
    lag_axis = lag_axis(:).';

    uSamples = t;
    for k = 1:numel(lag_axis)
        uSamples = [uSamples; shifted_time_local(t, lag_axis(k), cfg)]; %#ok<AGROW>
    end
    uSamples = uSamples(isfinite(uSamples));

    if cfg.zeroSignalBeforeTimeZero
        uSamples = uSamples(uSamples >= cfg.signalTimeZero_s);
        uSamples = [uSamples; cfg.signalTimeZero_s]; %#ok<AGROW>
    end

    uMin = min(uSamples);
    uMax = max(uSamples);
    if uMin > 0
        uGrid = logspace(log10(uMin), log10(uMax), cfg.shiftMcrNumGrid).';
    elseif uMax < 0
        uGrid = linspace(uMin, uMax, cfg.shiftMcrNumGrid).';
    else
        scale = min(abs(uSamples(uSamples ~= 0)));
        if isempty(scale) || scale <= 0
            scale = max(abs([uMin, uMax, 1])) * 1e-3;
        end
        z = linspace(asinh(uMin / scale), asinh(uMax / scale), cfg.shiftMcrNumGrid).';
        uGrid = scale * sinh(z);
    end

    u_axis = unique([uSamples; uGrid], "sorted");
end

function Y = eval_reconstructed_matrix_at_time(U_pc, S_pc, D_mean_q, best_par_struct, u_axis)
    V_u = eval_fitted_rsv(best_par_struct, u_axis(:));
    Y = U_pc * S_pc * V_u.';
    Y = Y + D_mean_q;
end

function C_time = initialize_shift_mcr_concentration(Y, R, iStart)
    nU = size(Y, 2);
    if iStart == 1
        [~, ~, V0] = svd(Y, "econ");
        nHave = min([size(V0, 2), R]);
        C_time = zeros(nU, R);
        C_time(:, 1:nHave) = abs(V0(:, 1:nHave));
        if nHave < R
            C_time(:, nHave+1:R) = rand(nU, R-nHave);
        end
    else
        C_time = rand(nU, R);
    end

    for r = 1:R
        if norm(C_time(:, r)) == 0
            C_time(:, r) = rand(nU, 1);
        end
    end
end

function [S_q, C_time] = normalize_shift_mcr_columns(S_q, C_time)
    R = size(C_time, 2);
    for r = 1:R
        scale = max(C_time(:, r));
        if scale > 0
            C_time(:, r) = C_time(:, r) / scale;
            S_q(:, r) = S_q(:, r) * scale;
        end
    end
end

function C_time = apply_shift_mcr_shape_constraint(C_time, cfg)
    switch lower(string(cfg.shiftMcrShapeConstraint))
        case "none"
            C_time = max(C_time, 0);
        case "unimodal_smooth"
            win = max(1, round(cfg.shiftMcrSmoothWindow));
            for r = 1:size(C_time, 2)
                c = max(C_time(:, r), 0);
                if win > 1
                    c = movmean(c, win, "Endpoints", "shrink");
                end
                c = project_unimodal_local(c);
                if win > 1
                    c = movmean(c, win, "Endpoints", "shrink");
                    c = project_unimodal_local(c);
                end
                C_time(:, r) = max(c, 0);
            end
        otherwise
            error("알 수 없는 cfg.shiftMcrShapeConstraint: %s", cfg.shiftMcrShapeConstraint);
    end
end

function y = project_unimodal_local(y)
    y = max(y(:), 0);
    if all(y == 0)
        return;
    end
    [~, peakIdx] = max(y);
    left = isotonic_increasing_local(y(1:peakIdx));
    right = -isotonic_increasing_local(-y(peakIdx:end));
    peakVal = max(left(end), right(1));
    left(end) = peakVal;
    right(1) = peakVal;
    y = [left; right(2:end)];
    y = max(y, 0);
end

function yfit = isotonic_increasing_local(y)
    y = y(:);
    n = numel(y);
    level = zeros(n, 1);
    weight = zeros(n, 1);
    startIdx = zeros(n, 1);
    endIdx = zeros(n, 1);
    nBlocks = 0;

    for i = 1:n
        nBlocks = nBlocks + 1;
        level(nBlocks) = y(i);
        weight(nBlocks) = 1;
        startIdx(nBlocks) = i;
        endIdx(nBlocks) = i;

        while nBlocks > 1 && level(nBlocks-1) > level(nBlocks)
            newWeight = weight(nBlocks-1) + weight(nBlocks);
            newLevel = (weight(nBlocks-1)*level(nBlocks-1) + weight(nBlocks)*level(nBlocks)) / newWeight;
            level(nBlocks-1) = newLevel;
            weight(nBlocks-1) = newWeight;
            endIdx(nBlocks-1) = endIdx(nBlocks);
            nBlocks = nBlocks - 1;
        end
    end

    yfit = zeros(n, 1);
    for b = 1:nBlocks
        yfit(startIdx(b):endIdx(b)) = level(b);
    end
end

function [S_q, C_time] = sort_shift_mcr_components(S_q, C_time)
    [~, peakIdx] = max(C_time, [], 1);
    [~, order] = sort(peakIdx, "ascend");
    S_q = S_q(:, order);
    C_time = C_time(:, order);
end

function C_eval = interpolate_shift_mcr_concentration(u_axis, C_time, uQuery)
    C_eval = zeros(numel(uQuery), size(C_time, 2));
    valid = (uQuery(:) >= min(u_axis)) & (uQuery(:) <= max(u_axis));
    if any(valid)
        C_eval(valid, :) = interp1(u_axis, C_time, uQuery(valid), "pchip", 0);
        C_eval(valid, :) = max(C_eval(valid, :), 0);
    end
end

function Xhat = reconstruct_shift_mcr_tensor(S_q, C_time, u_axis, t, lag_axis, cfg)
    nQ = size(S_q, 1);
    nT = numel(t);
    nLag = numel(lag_axis);
    Xhat = zeros(nQ, nT, nLag);

    for k = 1:nLag
        u = shifted_time_local(t, lag_axis(k), cfg);
        C_eval = interpolate_shift_mcr_concentration(u_axis, C_time, u);
        if cfg.zeroSignalBeforeTimeZero
            C_eval(u < cfg.signalTimeZero_s, :) = 0;
        end
        Xhat(:, :, k) = S_q * C_eval.';
    end
end

function plot_shift_mcr_results(q, t, q_spectra, t_conc, qt_target, qt_recon, shift_result, cfg)
    nComp = size(q_spectra, 2);
    uLabel = shifted_time_label_local(cfg);

    figure("Name", "Shift-MCR q spectra");
    plot(q, q_spectra, "LineWidth", 1.5);
    hold on; yline(0, "k--", "LineWidth", 1.0); hold off;
    grid on;
    xlabel("q");
    ylabel("Species spectrum (a.u.)");
    title(sprintf("Shift-constrained species spectra (rank = %d)", shift_result.rank));
    legend(compose("component %d", 1:nComp), "Location", "best");

    figure("Name", "Shift-MCR concentration profiles");
    plot(shift_result.u_axis, shift_result.C_time, "LineWidth", 1.5);
    grid on;
    xlabel("Shifted model time u = " + uLabel + " (s)");
    ylabel("Concentration profile (a.u.)");
    title("Shared concentration profiles C_r(" + uLabel + ")");
    legend(compose("component %d", 1:nComp), "Location", "best");
    if all(shift_result.u_axis > 0)
        set(gca, "XScale", "log");
    end

    figure("Name", "Shift-MCR concentrations at measured time");
    plot(t, t_conc, "LineWidth", 1.5);
    grid on;
    xlabel("Measured time t (s)");
    ylabel("Concentration profile (a.u.)");
    title("C_r(t) on original measured time grid");
    legend(compose("component %d", 1:nComp), "Location", "best");
    if all(t > 0)
        set(gca, "XScale", "log");
    end

    figure("Name", "Shift-MCR target vs reconstructed q-t maps");
    clim_all = [min([qt_target(:); qt_recon(:)]), max([qt_target(:); qt_recon(:)])];
    tiledlayout(1, 2, "TileSpacing", "compact");
    nexttile;
    imagesc(t, q, qt_target);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title("Target q-t map at zero lag");
    if all(t > 0), set(gca, "XScale", "log"); end
    nexttile;
    imagesc(t, q, qt_recon);
    axis xy; colorbar; caxis(clim_all);
    xlabel("Time (s)"); ylabel("q");
    title(sprintf("Shift-MCR q-t reconstruction (fit %.4f)", shift_result.fit));
    if all(t > 0), set(gca, "XScale", "log"); end

    figure("Name", "Shift-MCR concentration surfaces");
    tiledlayout(nComp, 1, "TileSpacing", "compact");
    for r = 1:nComp
        surf_r = zeros(numel(t), numel(cfg.delay_s));
        for k = 1:numel(cfg.delay_s)
            u = shifted_time_local(t, cfg.delay_s(k), cfg);
            c_eval = interpolate_shift_mcr_concentration(shift_result.u_axis, shift_result.C_time(:, r), u);
            if cfg.zeroSignalBeforeTimeZero
                c_eval(u < cfg.signalTimeZero_s) = 0;
            end
            surf_r(:, k) = c_eval;
        end
        nexttile;
        imagesc(1:numel(cfg.delay_s), t, surf_r);
        axis xy; colorbar;
        xlabel("Signed lag index");
        ylabel("Time (s)");
        title(sprintf("component %d: C_r(%s)", r, uLabel));
        if all(t > 0), set(gca, "YScale", "log"); end
    end
end

function [cp, Xhat] = cp_als_time_nonneg(X, W_time_tau, rankR, opts)
% Weighted CP-ALS for X(q,time,tau) with optional nonnegative constraints.
% Model: X(i,j,k) ~= sum_r lambda(r) * A_q(i,r) * B_time(j,r) * C_tau(k,r)
% Default constraint: B_time(:,r) >= 0 and C_tau(:,r) >= 0. A_q is signed.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, "maxIters"),   opts.maxIters = 500; end
    if ~isfield(opts, "tol"),        opts.tol = 1e-8; end
    if ~isfield(opts, "init"),       opts.init = "svd"; end
    if ~isfield(opts, "randomSeed"), opts.randomSeed = 1; end
    if ~isfield(opts, "verbose"),    opts.verbose = true; end
    if ~isfield(opts, "nonnegativeModes"), opts.nonnegativeModes = ["time", "tau"]; end
    if ~isfield(opts, "nnlsOptions"), opts.nnlsOptions = optimset('Display', 'off', 'MaxIter', 2000); end
    nonnegModes = lower(string(opts.nonnegativeModes));
    useTimeNonneg = any(nonnegModes == "time" | nonnegModes == "all");
    useTauNonneg = any(nonnegModes == "tau" | nonnegModes == "delay" | nonnegModes == "all");

    X = double(X);
    [nQ, nT, nDelay] = size(X);
    R = rankR;
    if R < 1 || R ~= round(R)
        error("CP rank는 양의 정수여야 합니다.");
    end

    if nargin < 2 || isempty(W_time_tau)
        W_time_tau = ones(nT, nDelay);
    end
    W_time_tau = double(W_time_tau);
    if ~isequal(size(W_time_tau), [nT, nDelay])
        error("W_time_tau 크기는 [nTime, nTau]여야 합니다.");
    end
    W_time_tau(~isfinite(W_time_tau) | W_time_tau < 0) = 0;
    sqrtW_time_tau = sqrt(W_time_tau);

    rng(opts.randomSeed);

    X1 = reshape(X, nQ, nT*nDelay);
    X2 = reshape(permute(X, [2 1 3]), nT, nQ*nDelay);
    X3 = reshape(permute(X, [3 1 2]), nDelay, nQ*nT);

    [A, B, C] = initialize_cp_factors_local(X, R, opts.init);
    if useTimeNonneg
        B = abs(B);
    end
    if useTauNonneg
        C = abs(C);
    end
    for r = 1:R
        if useTimeNonneg && norm(B(:, r)) == 0
            B(:, r) = rand(nT, 1);
        end
        if useTauNonneg && norm(C(:, r)) == 0
            C(:, r) = rand(nDelay, 1);
        end
    end

    normXw = weighted_tensor_norm_local(X, W_time_tau);
    if normXw == 0
        error("가중치가 적용된 입력 텐서의 norm이 0입니다. CP 분해를 수행할 수 없습니다.");
    end

    fit_history = nan(opts.maxIters, 1);
    relerr_history = nan(opts.maxIters, 1);
    prevFit = -Inf;
    reg = 1e-12;

    for iter = 1:opts.maxIters
        % Mode 1: q factor A_q, weighted unconstrained LS
        Z = khatri_rao_local(C, B);
        w1 = reshape(sqrtW_time_tau, [], 1);
        Zw = Z .* w1;
        X1w = X1 .* reshape(w1, 1, []);
        G = Zw.' * Zw + reg * eye(R);
        A = (X1w * Zw) / G;

        % Mode 2: time factor B_time, row-wise weighted LS/NNLS
        Z = khatri_rao_local(C, A);
        for j = 1:nT
            w2 = repelem(sqrtW_time_tau(j, :).', nQ);
            if all(w2 == 0)
                B(j, :) = 0;
                continue;
            end
            Zw = Z .* w2;
            y = X2(j, :).' .* w2;
            if useTimeNonneg
                B(j, :) = lsqnonneg(Zw, y, opts.nnlsOptions).';
            else
                G = Zw.' * Zw + reg * eye(R);
                B(j, :) = (y.' * Zw) / G;
            end
        end

        % Mode 3: tau/lag factor C_tau, row-wise weighted LS/NNLS
        Z = khatri_rao_local(B, A);
        for k = 1:nDelay
            w3 = repelem(sqrtW_time_tau(:, k), nQ);
            if all(w3 == 0)
                C(k, :) = 0;
                continue;
            end
            Zw = Z .* w3;
            y = X3(k, :).' .* w3;
            if useTauNonneg
                C(k, :) = lsqnonneg(Zw, y, opts.nnlsOptions).';
            else
                G = Zw.' * Zw + reg * eye(R);
                C(k, :) = (y.' * Zw) / G;
            end
        end

        [A, B, C] = stabilize_cp_columns_local(A, B, C);

        Xhat_iter = cp_reconstruct_local(A, B, C, ones(R,1));
        relerr = weighted_tensor_relerr_local(X, Xhat_iter, W_time_tau, normXw);
        fit = 1 - relerr;

        fit_history(iter) = fit;
        relerr_history(iter) = relerr;

        if opts.verbose && (iter == 1 || mod(iter, 10) == 0 || iter == opts.maxIters)
            fprintf("CP-ALS iter %4d: weighted relerr = %.6g, fit = %.6f\n", iter, relerr, fit);
        end

        if iter > 1 && abs(fit - prevFit) < opts.tol
            fit_history = fit_history(1:iter);
            relerr_history = relerr_history(1:iter);
            break;
        end
        prevFit = fit;

        if iter == opts.maxIters
            fit_history = fit_history(1:iter);
            relerr_history = relerr_history(1:iter);
        end
    end

    lambda = zeros(R, 1);
    for r = 1:R
        nr = norm(A(:, r));
        if nr > 0
            lambda(r) = nr;
            A(:, r) = A(:, r) / nr;
        else
            lambda(r) = 0;
        end
    end

    Xhat = cp_reconstruct_local(A, B, C, lambda);
    final_relerr = weighted_tensor_relerr_local(X, Xhat, W_time_tau, normXw);
    final_fit = 1 - final_relerr;
    full_relerr = norm(X(:) - Xhat(:)) / norm(X(:));

    cp = struct();
    cp.rank = R;
    cp.lambda = lambda;
    cp.A_q = A;
    if useTimeNonneg
        cp.B_time = max(B, 0);
    else
        cp.B_time = B;
    end
    if useTauNonneg
        cp.C_delay = max(C, 0);
    else
        cp.C_delay = C;
    end
    cp.C_tau = cp.C_delay;
    cp.fit = final_fit;
    cp.relerr = final_relerr;
    cp.full_relerr = full_relerr;
    cp.fit_history = fit_history;
    cp.relerr_history = relerr_history;
    cp.num_iters = numel(fit_history);
    cp.nonnegativeModes = nonnegModes;
    cp.weight_valid_fraction = nnz(W_time_tau > 0) / numel(W_time_tau);
    cp.constraint = sprintf("weighted CP; A_q signed; B_time nonnegative=%d; C_tau nonnegative=%d", ...
        useTimeNonneg, useTauNonneg);
end
function nrm = weighted_tensor_norm_local(X, W_time_tau)
    nrm2 = 0;
    for k = 1:size(X, 3)
        Xk = X(:, :, k);
        wk = reshape(W_time_tau(:, k), 1, []);
        nrm2 = nrm2 + sum((Xk.^2) .* wk, "all");
    end
    nrm = sqrt(nrm2);
end

function relerr = weighted_tensor_relerr_local(X, Xhat, W_time_tau, normXw)
    diff = X - Xhat;
    relerr = weighted_tensor_norm_local(diff, W_time_tau) / normXw;
end
function [A, B, C] = initialize_cp_factors_local(X, R, initMethod)
    [nQ, nT, nDelay] = size(X);
    initMethod = string(initMethod);

    switch lower(initMethod)
        case "svd"
            X1 = reshape(X, nQ, nT*nDelay);
            X2 = reshape(permute(X, [2 1 3]), nT, nQ*nDelay);
            X3 = reshape(permute(X, [3 1 2]), nDelay, nQ*nT);

            A = leading_left_singular_vectors_local(X1, R);
            B = abs(leading_left_singular_vectors_local(X2, R));
            C = leading_left_singular_vectors_local(X3, R);

        case "random"
            A = randn(nQ, R);
            B = rand(nT, R);
            C = randn(nDelay, R);

        otherwise
            error("알 수 없는 CP 초기화 방법입니다: %s", initMethod);
    end

    A(~isfinite(A)) = 0;
    B(~isfinite(B)) = 0;
    C(~isfinite(C)) = 0;

    for r = 1:R
        if norm(A(:, r)) == 0, A(:, r) = randn(nQ, 1); end
        if norm(B(:, r)) == 0, B(:, r) = rand(nT, 1); end
        if norm(C(:, r)) == 0, C(:, r) = randn(nDelay, 1); end
    end
end

function Ulead = leading_left_singular_vectors_local(M, R)
    [Utmp, ~, ~] = svd(M, "econ");
    nRows = size(M, 1);
    nHave = min(size(Utmp, 2), R);
    Ulead = zeros(nRows, R);
    if nHave > 0
        Ulead(:, 1:nHave) = Utmp(:, 1:nHave);
    end
    if nHave < R
        Ulead(:, nHave+1:R) = randn(nRows, R-nHave);
    end
end

function Z = khatri_rao_local(A, B)
    if size(A, 2) ~= size(B, 2)
        error("Khatri-Rao 입력 행렬의 column 수가 같아야 합니다.");
    end
    R = size(A, 2);
    Z = zeros(size(A,1)*size(B,1), R);
    for r = 1:R
        Z(:, r) = kron(A(:, r), B(:, r));
    end
end

function [A, B, C] = stabilize_cp_columns_local(A, B, C)
    R = size(A, 2);
    for r = 1:R
        nb = norm(B(:, r));
        if nb > 0
            B(:, r) = B(:, r) / nb;
            A(:, r) = A(:, r) * nb;
        end

        nc = norm(C(:, r));
        if nc > 0
            C(:, r) = C(:, r) / nc;
            A(:, r) = A(:, r) * nc;
        end
    end
end

function Xhat = cp_reconstruct_local(A, B, C, lambda)
    [nQ, R] = size(A);
    nT = size(B, 1);
    nDelay = size(C, 1);
    Xhat = zeros(nQ, nT, nDelay);

    for r = 1:R
        Xhat = Xhat + lambda(r) * reshape(A(:, r), [nQ, 1, 1]) .* ...
                                reshape(B(:, r), [1, nT, 1]) .* ...
                                reshape(C(:, r), [1, 1, nDelay]);
    end
end

function ensure_asinh_helpers_exist()
% fit_rsv.m 내부 plotting에서 필요한 asinhspace/asinhplot helper 생성
    if exist("asinhspace", "file") ~= 2
        fid = fopen("asinhspace.m", "w");
        if fid >= 0
            fprintf(fid, [ ...
                "function x = asinhspace(a,b,n)\n" ...
                "if nargin < 3, n = 100; end\n" ...
                "if a > 0 && b > 0\n" ...
                "    x = logspace(log10(a), log10(b), n);\n" ...
                "else\n" ...
                "    s = max(abs([a b]))/1000;\n" ...
                "    if s == 0, s = 1; end\n" ...
                "    u = linspace(asinh(a/s), asinh(b/s), n);\n" ...
                "    x = s*sinh(u);\n" ...
                "end\n" ...
                "x = x(:);\n" ...
                "end\n"]);
            fclose(fid);
        end
    end

    if exist("asinhplot", "file") ~= 2
        fid = fopen("asinhplot.m", "w");
        if fid >= 0
            fprintf(fid, [ ...
                "function h = asinhplot(x,y,varargin)\n" ...
                "h = plot(x,y,varargin{:});\n" ...
                "if all(x(:) > 0)\n" ...
                "    try, set(gca, 'XScale', 'log'); catch, end\n" ...
                "end\n" ...
                "end\n"]);
            fclose(fid);
        end
    end
end

function labels = make_delay_labels(delay_s)
% 초 단위 delay_s를 보기 좋은 label로 자동 변환
% 예:
%   0          -> "0s"
%   1e-14      -> "0.01ps"
%   1e-12      -> "1ps"
%   1e-9       -> "1ns"
%   1e-6       -> "1us"

    delay_s = delay_s(:);
    labels = strings(size(delay_s));

    for i = 1:numel(delay_s)
        x = delay_s(i);

        if x == 0
            labels(i) = "0s";

        elseif abs(x) < 1e-12
            labels(i) = sprintf("%.4gfs", x / 1e-15);

        elseif abs(x) < 1e-9
            labels(i) = sprintf("%.4gps", x / 1e-12);

        elseif abs(x) < 1e-6
            labels(i) = sprintf("%.4gns", x / 1e-9);

        elseif abs(x) < 1e-3
            labels(i) = sprintf("%.4gus", x / 1e-6);

        elseif abs(x) < 1
            labels(i) = sprintf("%.4gms", x / 1e-3);

        else
            labels(i) = sprintf("%.4gs", x);
        end
    end

    labels = cellstr(labels);
end

