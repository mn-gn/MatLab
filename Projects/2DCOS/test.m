clc; clear; close all;


q_file = "C:\Users\ming0\Desktop\mock data\q.csv";
time_file = "C:\Users\ming0\Desktop\mock data\t_81points.csv";
data_file = "C:\Users\ming0\Desktop\mock data\q_delta_S_noise_10%_81points.csv";

q = readmatrix(q_file);
t = readmatrix(time_file);
A = readmatrix(data_file);
A = A./q;

figure
imagesc(q,t,A');
set(gca,'YScale','log');
colormap("turbo")
% clim([-1.6e-4, 1.6e4])

window = [1e-14, 1e-4]; %[10e-12, 100e-12];

% [Phi, Psi] = run_q_space_2dcos(q,t,A,window,3);
[Phi_r, Psi_r, r] = run_r_space_2dcos(q,t,A,window,6);

figure;
plot(r,diag(Phi_r));
% plot(q,diag(Phi));

function [Phi, Psi] = run_q_space_2dcos(q, t, A_full, t_window, num_pcs)
    % run_q_space_2dcos: 불균일 시간 간격(Log-scale 등) 대응 및 PCA 노이즈 제거
    %
    % [입력]
    % q: q 벡터 (M x 1), t: 시간 벡터 (1 x N), A_full: 데이터 행렬 (M x N)
    % t_window: [start, end] (전체는 []), num_pcs: 유지할 PC 개수 (없으면 [])

    %% 1. 데이터 검증 및 방향 통일
    q = q(:); t = t(:)'; 
    [num_q, num_t] = size(A_full);
    if length(q) ~= num_q || length(t) ~= num_t
        error('데이터 차원이 일치하지 않습니다. q:%d, t:%d, Matrix:%dx%d', length(q), length(t), num_q, num_t);
    end

    %% 2. 윈도우(Window) 설정
    if nargin < 4 || isempty(t_window)
        t_window = [min(t), max(t)];
    end
    idx = find(t >= t_window(1) & t <= t_window(2));
    if length(idx) < 3
        error('윈도우 내 데이터가 부족합니다 (최소 3개 필요).');
    end
    
    t_sel = t(idx);
    A_sel = A_full(:, idx);
    N = length(t_sel);
    fprintf('분석 구간: %e ~ %e (%d points)\n', min(t_sel), max(t_sel), N);

    %% 3. 동적 스펙트럼 및 PCA 노이즈 제거
    A_mean = mean(A_sel, 2);
    A_dyn = A_sel - A_mean;

    if nargin >= 5 && ~isempty(num_pcs) && num_pcs > 0
        disp(['PCA 노이즈 제거 중... (PC: ', num2str(num_pcs), ')']);
        [U, S, V] = svd(A_dyn, 'econ');
        num_pcs = min(num_pcs, N);
        A_dyn = U(:, 1:num_pcs) * S(1:num_pcs, 1:num_pcs) * V(:, 1:num_pcs)';
    end

    %% 4. 불균일 간격 대응 2DCOS 계산 (Unevenly Spaced 2DCOS)
    disp('2D 상관 행렬 계산 중 (물리적 시간 가중치 적용)...');
    
    % (1) 동시성 스펙트럼 (Synchronous)
    Phi = (A_dyn * A_dyn') / (N - 1);

    % (2) 비동시성 스펙트럼 (Asynchronous) 
    % 물리적 시간 간격을 반영한 Noda-Hilbert Weight Matrix 계산
    W = zeros(N, N);
    
    % 사다리꼴 공식 기반의 적분 가중치 (Quadrature weights)
    dt = zeros(1, N);
    dt(1) = (t_sel(2) - t_sel(1)) / 2;
    dt(N) = (t_sel(N) - t_sel(N-1)) / 2;
    for i = 2:N-1
        dt(i) = (t_sel(i+1) - t_sel(i-1)) / 2;
    end
    
    % Hilbert 변환 커널 생성
    for j = 1:N
        for k = 1:N
            if j ~= k
                % dt(j)와 dt(k)의 기하평균을 사용하여 완벽한 반대칭성(Antisymmetry) 보장
                W(j, k) = (1/pi) * (sqrt(dt(j) * dt(k)) / (t_sel(k) - t_sel(j)));
            end
        end
    end
        
    % 비동시성 계산 (정규화 계수 적용)
    Psi = (A_dyn * W * A_dyn') / (N - 1);

    %% 5. 시각화 (Smooth Heatmap)
    disp('그래프 생성 중...');
    cmap = custom_bwr_colormap();
    figure('Position', [100, 100, 1100, 450]);

    subplot(1, 2, 1);
    vmax_phi = max(abs(Phi(:)));
    contourf(q, q, Phi, 100, 'LineColor', 'none'); 
    colormap(gca, cmap); caxis([-vmax_phi, vmax_phi]); colorbar; axis square;
    title('\Phi (Synchronous)', 'FontSize', 13); xlabel('q_1'); ylabel('q_2');

    subplot(1, 2, 2);
    vmax_psi = max(abs(Psi(:)));
    contourf(q, q, Psi, 100, 'LineColor', 'none'); 
    colormap(gca, cmap); caxis([-vmax_psi, vmax_psi]); colorbar; axis square;
    title('\Psi (Asynchronous)', 'FontSize', 13); xlabel('q_1'); ylabel('q_2');

    figure;
    M = Psi.*Phi;
    vmax_M = max(abs(M(:)));
    contourf(q,q,M,100,'LineColor','none');
    colormap(gca,cmap);caxis([-vmax_M, vmax_M]);colorbar;axis square;


    disp('분석 완료!');
end

function [Phi_r, Psi_r, r] = run_r_space_2dcos(q, t, A_full, t_window, num_pcs)
    % run_r_space_2dcos: r-공간(실제 거리) 2DCOS 분석 함수
    % (Sine Fourier Transform + 불균일 시간 간격 대응 + PCA 노이즈 제거)
    %
    % [입력]
    % q: q 벡터 (M x 1)
    % t: 시간 벡터 (1 x N)
    % A_full: q-공간 산란 신호 행렬 (M x N)
    % t_window: [시작시간, 종료시간] (전체는 [])
    % num_pcs: 유지할 PC 개수 (노이즈 제거용, 없으면 [])
    %
    % [출력]
    % Phi_r: r-공간 동시성 행렬
    % Psi_r: r-공간 비동시성 행렬
    % r: 분석에 사용된 r (거리) 벡터

    %% 1. 데이터 검증 및 방향 통일
    q = q(:); t = t(:)'; 
    [num_q, num_t] = size(A_full);
    if length(q) ~= num_q || length(t) ~= num_t
        error('데이터 차원이 일치하지 않습니다. q:%d, t:%d, Matrix:%dx%d', length(q), length(t), num_q, num_t);
    end

    %% 2. 윈도우(Window) 설정
    if nargin < 4 || isempty(t_window)
        t_window = [min(t), max(t)];
    end
    idx = find(t >= t_window(1) & t <= t_window(2));
    if length(idx) < 3
        error('윈도우 내 데이터가 부족합니다 (최소 3개 필요).');
    end
    
    t_sel = t(idx);
    A_sel = A_full(:, idx);
    N = length(t_sel);
    fprintf('분석 시간 구간: %e ~ %e (%d points)\n', min(t_sel), max(t_sel), N);

    %% 3. 사인 푸리에 변환 (Sine Fourier Transform: q -> r)
    disp('q-공간 데이터를 r-공간으로 푸리에 변환 중입니다 (f_Au(q)^2 보정 포함)...');
    
    r = linspace(1.0, 10.0, 1000)'; 
    dq = gradient(q);
    alpha = 0.03; 
    damping = exp(-alpha * (q.^2));
    
    % ---------- [추가된 부분] 금(Au)의 원자 산란 인자 계산 ----------
    % Cromer-Mann 계수 (Au 기준)
    a = [37.3027, 14.9306, 10.3425, 2.01229];
    b = [1.00810, 6.52550, 16.5100, 76.9117];
    c = 14.3992;
    
    s = q / (4 * pi); % 산란 벡터 s 정의
    f_Au = c * ones(size(q));
    for i = 1:4
        f_Au = f_Au + a(i) * exp(-b(i) * (s.^2));
    end
    % -----------------------------------------------------------
    
    % 적분식 내부 계산 (f_Au.^2 로 나누는 과정 추가됨!)
    % q * (dS(q) / f_Au(q)^2) * exp(-alpha*q^2) * dq
    A_weighted = (A_sel ./ (f_Au.^2)) .* (q .* damping .* dq);       
    
    integral_part = sin(r * q') * A_weighted;         
    A_r = (r / (2 * pi^2)) .* integral_part;

    %% 4. 동적 스펙트럼 및 PCA 노이즈 제거 (r-공간 기준)
    A_mean = mean(A_r, 2);
    A_dyn = A_r - A_mean;

    if nargin >= 5 && ~isempty(num_pcs) && num_pcs > 0
        disp(['PCA 노이즈 제거 중... (PC: ', num2str(num_pcs), ')']);
        [U, S, V] = svd(A_dyn, 'econ');
        num_pcs = min(num_pcs, N);
        A_dyn = U(:, 1:num_pcs) * S(1:num_pcs, 1:num_pcs) * V(:, 1:num_pcs)';
    end

    %% 5. 불균일 간격 대응 2DCOS 계산 (Unevenly Spaced 2DCOS)
    disp('r-공간 2D 상관 행렬 계산 중...');
    
    % (1) 동시성 스펙트럼 (Synchronous)
    Phi_r = (A_dyn * A_dyn') / (N - 1);

    % (2) 비동시성 스펙트럼 (Asynchronous) 
    W = zeros(N, N);
    dt = zeros(1, N);
    dt(1) = (t_sel(2) - t_sel(1)) / 2;
    dt(N) = (t_sel(N) - t_sel(N-1)) / 2;
    for i = 2:N-1
        dt(i) = (t_sel(i+1) - t_sel(i-1)) / 2;
    end
    
    for j = 1:N
        for k = 1:N
            if j ~= k
                % dt(j)와 dt(k)의 기하평균을 사용하여 완벽한 반대칭성(Antisymmetry) 보장
                W(j, k) = (1/pi) * (sqrt(dt(j) * dt(k)) / (t_sel(k) - t_sel(j)));
            end
        end
    end
    
    Psi_r = (A_dyn * W * A_dyn') / (N - 1);

    %% 6. 시각화 (Smooth Heatmap)
    disp('그래프 생성 중...');
    cmap = custom_bwr_colormap(); % 별도의 colormap 함수가 있다고 가정
    
    % ---------- [새로 추가된 부분] r-공간 전체 변환 데이터 (A_r) 시각화 ----------
    figure('Name', 'Time-Resolved Difference RDF', 'Position', [50, 100, 600, 450]);
    % X축: 시간(t_sel), Y축: 거리(r), Z축: 방사 분포 함수(A_r)
    surf(r, t_sel, A_r', 'EdgeColor', 'none'); % 경계선 없이 매끄럽게 출력
    set(gca, 'YScale', 'log');                % 시간(t) 축을 log scale로 설정
    set(gca, 'Ydir', 'reverse');
    view(0, 90);                              % 위에서 정면으로 내려다보게 설정 (2D 맵 형태)
    
    colormap(gca, cmap);
    vmax_Ar = max(abs(A_r(:)));               % 양수, 음수 색상 균형 맞추기
    caxis([-vmax_Ar, vmax_Ar]);
    colorbar;
    axis tight;
    
    title('Time-Resolved r-space Data (r^2\DeltaS)', 'FontSize', 13);
    xlabel('r (\AA)');
    ylabel('Time (t)');
    % -------------------------------------------------------------------------
    
    % ---------- 기존 2DCOS 분석 결과 시각화 ----------
    figure('Name', '2D Correlation Spectroscopy (r-space)', 'Position', [700, 100, 1100, 450]);
    
    subplot(1, 2, 1);
    vmax_phi = max(abs(Phi_r(:)));
    contourf(r, r, Phi_r, 100, 'LineColor', 'none'); 
    colormap(gca, cmap); caxis([-vmax_phi, vmax_phi]); colorbar; axis square;
    title('\Phi (Synchronous in r-space)', 'FontSize', 13); 
    xlabel('r_1 (\AA)'); ylabel('r_2 (\AA)');
    
    subplot(1, 2, 2);
    vmax_psi = max(abs(Psi_r(:)));
    contourf(r, r, Psi_r, 100, 'LineColor', 'none'); 
    colormap(gca, cmap); caxis([-vmax_psi, vmax_psi]); colorbar; axis square;
    title('\Psi (Asynchronous in r-space)', 'FontSize', 13); 
    xlabel('r_1 (\AA)'); ylabel('r_2 (\AA)');
    
    disp('분석 완료!');
end


function cmap = custom_bwr_colormap()
    n = 256; c = ones(n, 3);
    c(1:n/2, 1:2) = repmat(linspace(0, 1, n/2)', 1, 2);
    c(n/2+1:end, 2:3) = repmat(linspace(1, 0, n/2)', 1, 2);
    cmap = c;
end