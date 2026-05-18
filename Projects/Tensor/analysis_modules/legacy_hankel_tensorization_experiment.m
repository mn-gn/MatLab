%% 1. 데이터 불러오기 및 SVD 수행 (기존과 동일)
clc; clear; close all;
data = readmatrix('q_delta_S_noise_10%_81points.csv');
q = readmatrix('q.csv');
t = readmatrix('t_81points.csv');

[U, S, V] = svd(data, 'econ');

%% 2. fit_rsv 함수를 사용하여 주요 RSV 피팅
rsv_fit_range = 1:3;  % 피팅을 적용할 RSV 성분 범위 (예: 1~2번)
num_tau = 3;          % 다중 지수 감쇠 성분(tau)의 개수
num_start = 100;        % MultiStart 최적화 시작점 개수

% 선택한 RSV에 대해서만 피팅 함수 호출
[best_par_struct, V_fit_original, save_struct] = fit_rsv(t, V(:, rsv_fit_range), num_tau, num_start);

%% 3. 동시간 간격(Uniform Time Interval) 그리드 생성
num_uniform_points = 1000;
t_uniform = linspace(1e-13, 1e-8, num_uniform_points)';
V_uniform = zeros(length(t_uniform), size(V, 2));

% ==============================================================
%% 4. best_par_struct를 활용한 모델 기반 V 재구성
% ==============================================================
% fit_rsv.m 의 수식을 참조하여 파라미터 추출
sigma = best_par_struct.FWHM / 2.3548;
inv_sqrt2_sigma = 1 / (sqrt(2) * sigma);
t0 = best_par_struct.t0;
tau = best_par_struct.tau;

% 공통으로 사용되는 가우시안 포락선 (수치적 안정을 위해 사용)
gauss_env_uniform = exp(-((t_uniform - t0).^2) ./ (2 * sigma^2));

for i = 1:length(rsv_fit_range)
    rsv_idx = rsv_fit_range(i);
    A = best_par_struct.A(i, :);   % 해당 RSV의 진폭 벡터
    A0 = best_par_struct.A0(i);    % 해당 RSV의 오프셋
    
    % 오프셋으로 y_calc 초기화
    y_calc = A0 * ones(length(t_uniform), 1);
    
    % 각 tau 성분별로 감쇠 곡선 계산 및 합산
    for k = 1:length(tau)
        % 여오차 함수에 들어갈 인자 Z 계산
        Z = ((sigma^2)/tau(k) - (t_uniform - t0)) * inv_sqrt2_sigma;
        comp = zeros(length(t_uniform), 1);
        
        % 수치적 안정성을 위한 조건 분기 (fit_rsv.m 내부 논리 반영)
        idx_pos = (Z >= 0);
        idx_neg = (Z < 0);
        
        % 조건 1: Z >= 0 인 경우 (erfcx를 사용하여 지수항 오버플로우 방지)
        if any(idx_pos)
            comp(idx_pos) = (A(k)/2) * gauss_env_uniform(idx_pos) .* erfcx(Z(idx_pos));
        end
        
        % 조건 2: Z < 0 인 경우 (일반적인 exp와 erfc 사용)
        if any(idx_neg)
            term1 = (sigma^2)/(2*tau(k)^2) - (t_uniform(idx_neg) - t0)/tau(k);
            comp(idx_neg) = (A(k)/2) * exp(term1) .* erfc(Z(idx_neg));
        end
        
        y_calc = y_calc + comp;
    end
    
    % 계산된 모델 곡선을 동시간 간격 V 행렬에 할당
    V_uniform(:, rsv_idx) = y_calc;
end

%% 5. 나머지 RSV 포인트들에 대한 스플라인(Spline) 보간
rsv_spline_range = 4:6; % 스플라인을 적용할 RSV 성분 범위 (예: 3~5번)
for i = rsv_spline_range
    V_uniform(:, i) = spline(t, V(:, i), t_uniform);
end

%% 6. 전체 산란 데이터 재구성
selected_rsvs = [rsv_fit_range, rsv_spline_range];
U_selected = U(:, selected_rsvs);
S_selected = S(selected_rsvs, selected_rsvs);
V_uniform_selected = V_uniform(:, selected_rsvs);

% 동시간 간격으로 재구성된 최종 데이터 (A_new = U * S * V_new')
data_reconstructed = U_selected * S_selected * V_uniform_selected';

%% 7. 결과 비교 시각화
figure;
subplot(2, 1, 1);
plot(t, V(:, 1), 'ko', 'DisplayName', 'Original RSV 1'); hold on;
plot(t, V_fit_original(:, 1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fitted RSV 1');
set(gca, 'XScale', 'log');
title('Original RSV vs Fitted RSV (on Original Time)');
legend; xlabel('Time (t)'); ylabel('Amplitude');

subplot(2, 1, 2);
plot(t_uniform, V_uniform(:, 1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed RSV 1');
set(gca, 'XScale', 'linear');
title('Reconstructed RSV (on Uniform Time)');
legend; xlabel('Uniform Time'); ylabel('Amplitude');


% ==============================================================
%% 8. Hankel Tensorization (시간축 기준)
% ==============================================================

% 재구성된 데이터(data_reconstructed)의 크기 확인
% num_q: q 포인트의 개수, N: 동시간 간격 시간축의 전체 길이
[num_q, N] = size(data_reconstructed);

% 윈도우 길이 L 설정 (사용자가 분석 목적에 맞게 조절 가능)
% 일반적으로 시계열 분석에서는 N의 절반에 가까운 값을 많이 사용합니다.
L = floor(N / 2); 
K = N - L + 1;

% Hankel 행렬의 인덱스 생성
% 벡터화된 연산을 위해 [L x K] 크기의 인덱스 행렬을 만듭니다.
% 1열은 1~L, 2열은 2~L+1, ..., K열은 K~N의 인덱스를 갖습니다.
hankel_idx = (1:L)' + (0:K-1); 

% 3차원 Hankel Tensor 초기화 (크기: [num_q x L x K])
hankel_tensor = zeros(num_q, L, K);

% 각 q 포인트(행)에 대해 시계열 데이터를 Hankel 행렬로 변환하여 텐서에 할당
for i = 1:num_q
    % 1. i번째 q의 시계열 데이터를 1차원 벡터로 추출
    ts = data_reconstructed(i, :);
    
    % 2. 슬라이딩 윈도우 인덱싱을 통해 L x K 크기의 Hankel 행렬 생성
    H_matrix = ts(hankel_idx);
    
    % 3. 좌항의 3차원 슬라이스(1 x L x K) 크기에 맞게 형태 변환하여 할당
    hankel_tensor(i, :, :) = reshape(H_matrix, [1, L, K]);
end

% 변환 결과 출력
disp('Hankel Tensor의 크기 [num_q, L, K]:');
disp(size(hankel_tensor));

%% 9. 특정 q 포인트의 Hankel 행렬 시각화 예시
% 첫 번째 q 포인트(q index = 1)의 Hankel 행렬 성분을 시각화합니다.
figure;
imagesc(squeeze(hankel_tensor(1, :, :)));
title('Hankel Matrix for the 1st q-point');
xlabel('Window Index (K)');
ylabel('Time Delay Index (L)');
colorbar;

% ==============================================================
%% 10. CP-ALS Tensor Decomposition (혼합 제약 조건 적용)
% Mode 1 (q축): 비음수 제약 없음 (Unconstrained)
% Mode 2 (시간 지연 L), Mode 3 (시간 윈도우 K): 비음수 제약 (Non-negative)
% ==============================================================

% 텐서 차원 확인
[I, J, K] = size(hankel_tensor);

% 분석할 성분(Rank)의 개수 설정 (사용자가 원하는 성분 수로 조절하세요)
R = 3; 

% 1. 텐서의 Matricization (Unfolding)
% MATLAB 메모리 구조(Column-major)에 맞게 각 차원으로 텐서를 펼칩니다.
X1 = reshape(hankel_tensor, I, J * K);
X2 = reshape(permute(hankel_tensor, [2, 1, 3]), J, I * K);
X3 = reshape(permute(hankel_tensor, [3, 1, 2]), K, I * J);

% 2. 인자 행렬(Factor Matrices) 무작위 초기화
A = rand(I, R); % q축 (제약 없음)
B = rand(J, R); % 시간 지연축 (비음수)
C = rand(K, R); % 시간 윈도우축 (비음수)

% ALS 알고리즘 설정
max_iter = 100;
tol = 1e-4;
err_old = inf;

disp('Starting CP-ALS with mixed constraints...');

for iter = 1:max_iter
    
    % -----------------------------------------------------
    % Update A (Mode 1: q축) - 비음수 제약 없음 (Unconstrained)
    % A = X1 * pinv(KhatriRao(C, B)') 식을 효율적으로 계산합니다.
    % -----------------------------------------------------
    KR_CB = khatrirao_prod(C, B); % 크기: (K*J) x R
    V_A = (C'*C) .* (B'*B);       % 정규 방정식의 분모 부분
    A = (X1 * KR_CB) *pinv(V_A);       % 일반 최소제곱 업데이트
    
    % -----------------------------------------------------
    % Update B (Mode 2: 시간 지연축 L) - 비음수 제약 (Non-negative)
    % lsqnonneg 를 사용하여 행렬 B의 각 행에 대해 비음수 제약을 보장합니다.
    % -----------------------------------------------------
    KR_CA = khatrirao_prod(C, A); % 크기: (K*I) x R
    for j = 1:J
        B(j, :) = lsqnonneg(KR_CA, X2(j, :)')';
    end
    
    % -----------------------------------------------------
    % Update C (Mode 3: 시간 윈도우축 K) - 비음수 제약 (Non-negative)
    % -----------------------------------------------------
    KR_BA = khatrirao_prod(B, A); % 크기: (J*I) x R
    for k = 1:K
        C(k, :) = lsqnonneg(KR_BA, X3(k, :)')';
    end
    
    % 스케일 정규화 (B와 C의 크기가 무한히 커지는 것을 방지하고 A에 흡수)
    for r = 1:R
        norm_B = norm(B(:, r));
        if norm_B > 1e-12 % 0에 너무 가까운 경우 나누지 않음
            B(:, r) = B(:, r) / norm_B;
        else
            norm_B = 1; % 0일 경우 A에 0이 곱해지는 것을 방지
        end
        
        norm_C = norm(C(:, r));
        if norm_C > 1e-12
            C(:, r) = C(:, r) / norm_C;
        else
            norm_C = 1;
        end
        
        A(:, r) = A(:, r) * norm_B * norm_C;
    end
    
    % 오차(Error) 계산 및 수렴 판정
    KR_CB_new = khatrirao_prod(C, B);
    X1_hat = A * KR_CB_new';
    err = norm(X1 - X1_hat, 'fro') / norm(X1, 'fro');
    
    if abs(err_old - err) < tol
        fprintf('Converged at iteration %d with relative error: %.4f\n', iter, err);
        break;
    end
    err_old = err;
    
    if mod(iter, 10) == 0
        fprintf('Iteration %d, Error: %.4f\n', iter, err);
    end
end

% ==============================================================
%% 11. 추출된 성분 시각화
% ==============================================================
figure('Name', 'CP-ALS 3 Components', 'Position', [100, 100, 1200, 800]);

for r = 1:R
    % Mode 1: q축 스펙트럼 (제약 없음 -> 양수/음수 혼재)
    subplot(R, 3, (r-1)*3 + 1);
    plot(q,A(:, r), 'k-', 'LineWidth', 1.5);
    title(sprintf('Component %d - Mode 1 (q-axis)', r));
    ylabel('Amplitude');
    % 마지막 행에만 x축 라벨 추가
    if r == R, xlabel('q index'); end 
    grid on;
    
    % Mode 2: 시간 지연축 (비음수 제약)
    subplot(R, 3, (r-1)*3 + 2);
    plot(B(:, r), 'b-', 'LineWidth', 1.5);
    title(sprintf('Component %d - Mode 2 (Delay)', r));
    if r == R, xlabel('Delay Index (L)'); end
    grid on;
    
    % Mode 3: 시간 윈도우축 (비음수 제약)
    subplot(R, 3, (r-1)*3 + 3);
    plot(C(:, r), 'r-', 'LineWidth', 1.5);
    title(sprintf('Component %d - Mode 3 (Window)', r));
    if r == R, xlabel('Window Index (K)'); end
    grid on;
end



% ==============================================================
%% 12. De-Hankelization (대각 평균화) 및 실제 농도(Kinetics) 재구성
% ==============================================================
% Mode 2(B)와 Mode 3(C)를 조합하여 원래의 시간 길이 N에 대한 농도 궤적을 복원합니다.
% N = L + K - 1 이며, 이는 t_uniform의 길이와 정확히 일치합니다.

disp('Reconstructing actual kinetics via De-Hankelization...');

% 전체 시간 포인트 N에 대한 각 성분의 농도를 저장할 행렬 (크기: N x R)
kinetics_reconstructed = zeros(N, R);

for r = 1:R
    % r번째 성분의 Hankel 행렬 복원 (L x K 크기)
    H_r = B(:, r) * C(:, r)'; 
    
    % Anti-diagonal(반대각선) 평균화를 통한 1D 시계열(농도) 복원
    for n = 1:N
        % Hankel 행렬의 특성상 인덱스 l + k - 1 = n 인 성분들이 같은 시간대(n)를 의미합니다.
        % 따라서 k = n - l + 1 이 됩니다.
        l_min = max(1, n - K + 1);
        l_max = min(L, n);
        
        % 시간 n에 해당하는 대각 성분들의 인덱스 배열 생성
        l_idx = l_min:l_max;
        k_idx = n - l_idx + 1;
        
        % 행렬에서 해당 성분들을 추출하여 평균값 계산
        vals = zeros(1, length(l_idx));
        for idx = 1:length(l_idx)
            vals(idx) = H_r(l_idx(idx), k_idx(idx));
        end
        
        % 재구성된 농도 행렬에 할당
        kinetics_reconstructed(n, r) = mean(vals);
    end
end

% ==============================================================
%% 13. 실제 시간(t_uniform)에 따른 최종 주성분 시각화
% ==============================================================
% 기존의 Mode 2, Mode 3 개별 플롯 대신, 물리적으로 의미 있는 통합 결과를 보여줍니다.
figure('Name', 'Final Components vs Actual Time', 'Position', [150, 150, 1000, 450]);

for r = 1:R
    % Mode 1: 구조적 지문 (q축 산란 패턴)
    % A 행렬을 사용자의 q 벡터에 맞춰 그립니다.
    subplot(1, 2, 1);
    plot(q, A(:, r), 'LineWidth', 1.5, 'DisplayName', sprintf('Component %d', r));
    hold on;
    
    % 재구성된 농도 궤적 (시간축 동역학)
    % 복원된 농도를 실제 물리적 시간축인 t_uniform 에 맞춰 그립니다.
    subplot(1, 2, 2);
    plot(t_uniform, kinetics_reconstructed(:, r), 'LineWidth', 1.5, 'DisplayName', sprintf('Component %d', r));
    hold on;
end

% 좌측 그래프 (q-space) 장식
subplot(1, 2, 1);
title('Mode 1: Structural Fingerprint (q-axis)');
xlabel('q (\AA^{-1})'); ylabel('Scattering Amplitude (\Delta S)');
legend('Location', 'best'); grid on;

% 우측 그래프 (Time-space) 장식
subplot(1, 2, 2);
title('Reconstructed Kinetics vs Actual Time (t\_uniform)');
xlabel('Time (t\_uniform)'); ylabel('Concentration (a.u.)');
% 만약 초기 피코초 반응부터 수 나노초까지 한눈에 보시려면 'XScale'을 'log'로 변경하세요.
set(gca, 'XScale', 'linear'); 
legend('Location', 'best'); grid on;


%% ==============================================================
% 로컬 함수: Khatri-Rao Product 계산
% ==============================================================
function KR = khatrirao_prod(Mat1, Mat2)
    % 두 행렬의 열 단위 크로네커 곱(Kronecker product)을 수행합니다.
    R_cols = size(Mat1, 2);
    KR = zeros(size(Mat1, 1) * size(Mat2, 1), R_cols);
    for r = 1:R_cols
        KR(:, r) = kron(Mat1(:, r), Mat2(:, r));
    end
end