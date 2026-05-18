%% build_svd_rsv_delay_tensor.m
% q / time / q-delta-S 행렬을 읽고,
% 1) SVD 수행
% 2) 주성분 RSV = V(:,1:numPC)를 fit_rsv로 피팅
% 3) 피팅된 RSV 모델을 이용해 지정된 시간 지연별로 전체 행렬을 재생성
% 4) D(q,t,delay) 형태의 시간 지연 텐서로 저장
%
% 필요 파일:
%   - q.csv
%   - t_81points.csv
%   - q_delta_S_noise_10%_81points.csv
%   - fit_rsv.m
%
% 출력:
%   - svd_rsv_delay_tensor_result.mat
%   - output_delay_slices/*.csv        각 시간 지연별 재생성 행렬
%   - output_delay_slices/delay_block_q_by_timeDelay.csv
%
% NOTE
%   fit_rsv.m 내부에서 asinhspace/asinhplot을 호출합니다.
%   현재 경로에 해당 함수가 없으면, 이 스크립트가 단순 helper m-file을 자동 생성합니다.

clear; clc;

%% ===================== 사용자 설정 =====================
cfg = struct();

% 입력 파일명
cfg.qFile    = "q.csv";
cfg.tFile    = "t_81points.csv";
cfg.dataFile = "q_delta_S_noise_10%_81points.csv";

% SVD / RSV 피팅 설정
cfg.numPC    = 3;     % 사용할 주성분 개수. 필요하면 1,2,3,...으로 변경
cfg.numTau   = 3;     % fit_rsv의 exponential component 개수
cfg.numStart = 100;    % MultiStart 시작점 개수. 빠르게 테스트하려면 1~5로 줄이세요.
cfg.useFittedTauAxis = false; % true이면 fitted lifetime을 축으로 사용. 기본은 log-spaced tau axis.

% 데이터 전처리
cfg.centerData = false;   % true면 q별 평균을 빼고 SVD 후 재구성 시 다시 더함
cfg.fixSVDSign = false;    % true면 각 RSV의 최대 절대값 지점이 양수가 되도록 SVD 부호 정렬
cfg.useLogTimeGridForTensor = true;
cfg.numTensorTimePoints = [];  % []이면 원본 time point 개수와 동일하게 logspace 재구성

% tau 축 목록 [초]
% RSV fitting 결과는 log-spaced tau 축에서 비어 있는 시간점을 평가/보간하는 모델로 사용한다.
tau_min_s = 1e-14;       % 0.01 ps
tau_max_s = 100e-9;      % 100 ns
num_log_tau = 100;       % 0을 제외한 log-spaced tau 개수
cfg.includeZeroTauSlice = true;
if cfg.includeZeroTauSlice
    cfg.delay_s = [0, logspace(log10(tau_min_s), log10(tau_max_s), num_log_tau)];
else
    cfg.delay_s = logspace(log10(tau_min_s), log10(tau_max_s), num_log_tau);
end
cfg.delay_labels = cellstr(compose("tau_%.4gs", cfg.delay_s(:)));
if cfg.includeZeroTauSlice
    cfg.delay_labels{1} = "tau_0s";
end

% 시간 지연의 부호/방향
%   "delay_signal"       : 지연된 신호 D(t-delay)를 생성. 일반적인 time-delay shift.
%   "evaluate_later_time": 더 늦은 관측시간 D(t+delay)를 평가.
cfg.shiftMode = "evaluate_later_time";

% 저장 설정
cfg.outMatFile = "svd_rsv_tau_tensor_result.mat";
cfg.outCsvDir  = "output_tau_slices";
cfg.saveCsvSlices = true;

% 그림 설정
cfg.makePlots = true;

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
singVals = diag(S);
explained = singVals.^2 ./ sum(singVals.^2);
explained_cumsum = cumsum(explained);

numPC = min([cfg.numPC, size(U,2), size(V,2)]);
U_pc = U(:, 1:numPC);
S_pc = S(1:numPC, 1:numPC);
V_pc = V(:, 1:numPC);       % RSV: time x numPC

% SVD 부호는 임의이므로, RSV의 가장 큰 절대값 지점이 양수로 오도록 정렬
if cfg.fixSVDSign
    for k = 1:numPC
        [~, idxMax] = max(abs(V_pc(:, k)));
        if V_pc(idxMax, k) < 0
            U_pc(:, k) = -U_pc(:, k);
            V_pc(:, k) = -V_pc(:, k);
            U(:, k)    = -U(:, k);
            V(:, k)    = -V(:, k);
        end
    end
end

fprintf("Using numPC = %d. Cumulative explained variance = %.4f\n", ...
    numPC, explained_cumsum(numPC));

%% ===================== RSV 피팅 =====================
% fit_rsv 입력: t = num_t x 1, V = num_t x num_rsv
[best_par_struct, V_fit, fit_save_struct] = fit_rsv(t, V_pc, cfg.numTau, cfg.numStart);

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
    V_fit_tensor_time = eval_fitted_rsv(best_par_struct, t); %#ok<NASGU>
    D_tensor_time = U_pc * S_pc * V_fit_tensor_time.';
    D_tensor_time = D_tensor_time + D_mean_q; %#ok<NASGU>

    fprintf("Tensor time axis reconstructed on logspace: %d points, %.4g s to %.4g s\n", ...
        nT, min(t), max(t));
else
    V_fit_tensor_time = V_fit; %#ok<NASGU>
    D_tensor_time = D; %#ok<NASGU>
end

if cfg.useFittedTauAxis
    cfg.delay_s = best_par_struct.tau(:).';
    cfg.delay_labels = cellstr(compose("tau%d_%.4gs", (1:numel(cfg.delay_s)).', cfg.delay_s(:)));
end

tau_axis_s = cfg.delay_s(:);
tau_axis_labels = cfg.delay_labels(:);

%% ===================== 시간 지연 텐서 생성 =====================
nDelay = numel(tau_axis_s);
delay_tensor = zeros(nQ, nT, nDelay);
V_fit_shifted = zeros(nT, numPC, nDelay);
t_eval_all = zeros(nT, nDelay);

for iDelay = 1:nDelay
    tauShift = tau_axis_s(iDelay);

    switch cfg.shiftMode
        case "delay_signal"
            % 신호 자체를 delay만큼 늦춤: D_delay(t) = D_model(t - delay)
            t_eval = t - tauShift;
        case "evaluate_later_time"
            % 원래 모델을 더 늦은 시간에서 평가: D_delay(t) = D_model(t + delay)
            t_eval = t + tauShift;
        otherwise
            error("알 수 없는 cfg.shiftMode: %s", cfg.shiftMode);
    end

    t_eval_all(:, iDelay) = t_eval;

    % fitted RSV model을 임의 시간축에서 평가
    V_shift = eval_fitted_rsv(best_par_struct, t_eval);
    V_fit_shifted(:, :, iDelay) = V_shift;

    % SVD 저차원 모델로 전체 q x t 행렬 재생성
    D_shift = U_pc * S_pc * V_shift.';

    % centerData를 썼다면 평균 성분을 다시 더함
    D_shift = D_shift + D_mean_q;

    delay_tensor(:, :, iDelay) = D_shift;
end

% 2D block 형태도 함께 생성: q x [time(delay1), time(delay2), ...]
delay_block_q_by_timeDelay = reshape(delay_tensor, nQ, nT*nDelay);
tau_tensor = delay_tensor;
tau_block_q_by_timeTau = delay_block_q_by_timeDelay;

fprintf("Tau tensor created: %d(q) x %d(time) x %d(tau)\n", ...
    size(delay_tensor,1), size(delay_tensor,2), size(delay_tensor,3));

%% ===================== 저장 =====================
save(cfg.outMatFile, ...
    "cfg", "q", "t", "D", "D_work", "D_mean_q", ...
    "U", "S", "V", "singVals", "explained", "explained_cumsum", ...
    "numPC", "U_pc", "S_pc", "V_pc", ...
    "best_par_struct", "V_fit", "fit_save_struct", ...
    "tau_axis_s", "tau_axis_labels", "t_eval_all", "V_fit_shifted", ...
    "tau_tensor", "tau_block_q_by_timeTau", ...
    "delay_tensor", "delay_block_q_by_timeDelay", ...
    "-v7.3");

fprintf("Saved MAT result: %s\n", cfg.outMatFile);

if cfg.saveCsvSlices
    if ~exist(cfg.outCsvDir, "dir")
        mkdir(cfg.outCsvDir);
    end

    writematrix(q, fullfile(cfg.outCsvDir, "q.csv"));
    writematrix(t, fullfile(cfg.outCsvDir, "t_original_s.csv"));
    writematrix(cfg.delay_s(:), fullfile(cfg.outCsvDir, "delay_s.csv"));
    writematrix(tau_axis_s, fullfile(cfg.outCsvDir, "tau_axis_s.csv"));
    writematrix(delay_block_q_by_timeDelay, fullfile(cfg.outCsvDir, "delay_block_q_by_timeDelay.csv"));
    writematrix(tau_block_q_by_timeTau, fullfile(cfg.outCsvDir, "tau_block_q_by_timeTau.csv"));

    for iDelay = 1:nDelay
        safeLabel = regexprep(tau_axis_labels{iDelay}, "[^a-zA-Z0-9_.-]", "_");
        outName = sprintf("D_tau_%02d_%s.csv", iDelay, safeLabel);
        writematrix(delay_tensor(:, :, iDelay), fullfile(cfg.outCsvDir, outName));
    end

    fprintf("Saved CSV slices in folder: %s\n", cfg.outCsvDir);
end

%% ===================== 간단 플롯 =====================
if cfg.makePlots
    figure("Name", "SVD singular values");
    semilogy(singVals, "o-");
    grid on;
    xlabel("Component index");
    ylabel("Singular value");
    title("SVD singular values");

    figure("Name", "Explained variance");
    plot(100*explained_cumsum, "o-");
    grid on;
    xlabel("Component index");
    ylabel("Cumulative explained variance (%)");
    title("Cumulative explained variance");

    % 원래 RSV와 피팅 RSV 비교
    figure("Name", "RSV fit summary");
    tiledlayout(numPC, 1, "TileSpacing", "compact");
    for k = 1:numPC
        nexttile;
        plot(t, V_pc(:, k), "k.", "MarkerSize", 12); hold on;
        plot(t, V_fit(:, k), "r-", "LineWidth", 1.5);
        grid on;
        xlabel("Time (s)");
        ylabel(sprintf("RSV%d", k));
        legend("SVD RSV", "fit_rsv", "Location", "best");
    end

    % 첫 번째 q 지점에 대해 delay별 trace 예시
    figure("Name", "Example tau-axis traces at first q");
    hold on;
    for iDelay = 1:nDelay
        plot(t, squeeze(tau_tensor(1, :, iDelay)), "DisplayName", tau_axis_labels{iDelay});
    end
    hold off; grid on;
    xlabel("Original time grid t (s)");
    ylabel(sprintf("Reconstructed signal at q = %.4g", q(1)));
    legend("Location", "best");
    title("Tau-axis reconstructed traces");
end

%% ===================== Local functions =====================
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
        % 쉼표/세미콜론/탭을 공백으로 통일
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

function ensure_asinh_helpers_exist()
% fit_rsv.m의 plotting 부분에서 필요한 asinhspace/asinhplot이 없을 때
% 최소 동작 helper를 현재 폴더에 생성한다.
    if exist("asinhspace", "file") ~= 2
        fid = fopen("asinhspace.m", "w");
        if fid < 0
            warning("asinhspace.m helper를 생성하지 못했습니다.");
        else
            fprintf(fid, [ ...
                "function x = asinhspace(a,b,n)\n" ...
                "%% Minimal helper generated by build_svd_rsv_delay_tensor.m\n" ...
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
        if fid < 0
            warning("asinhplot.m helper를 생성하지 못했습니다.");
        else
            fprintf(fid, [ ...
                "function h = asinhplot(x,y,varargin)\n" ...
                "%% Minimal helper generated by build_svd_rsv_delay_tensor.m\n" ...
                "h = plot(x,y,varargin{:});\n" ...
                "if all(x(:) > 0)\n" ...
                "    try, set(gca, 'XScale', 'log'); catch, end\n" ...
                "end\n" ...
                "end\n"]);
            fclose(fid);
        end
    end
end


