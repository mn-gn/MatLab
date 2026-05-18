clc; clear; close all;

%% [1] 환경 설정 및 모델 로드
% 학습된 모델과 정규화 기준값 로드
modelPath = 'N2N_Xray_Model.mat';
if ~exist(modelPath, 'file'), error('모델 파일이 없습니다.'); end
load(modelPath); % net, maxLogVal 로드 

% 파일 및 경로 설정
h5_file_path = "Z:\KEK_2512_h5\results\KEK_2512\I3_MeCN\run5\raw\images\image_pairs.h5"; 
on_dataset   = "/delays/delay_1.50e-10/on"; 
off_dataset  = "/delays/delay_1.50e-10/off"; 
poni_path    = "Z:\KEK_2512_h5\KEK_2512_1.poni";
mask_path    = "Z:\KEK_2512_h5\run1_neg_001-mask.edf";

%% [2] 전체 프레임 수 감지 및 네트워크 재조립
info_on  = h5info(h5_file_path, on_dataset);
info_off = h5info(h5_file_path, off_dataset);

n_frames = info_on.Dataspace.Size(1);
m_frames = info_off.Dataspace.Size(1);
[H, W]   = deal(info_on.Dataspace.Size(2), info_on.Dataspace.Size(3));

% 마스크 로드 및 Transpose 정렬 [cite: 2026-02-11]
mask = read_mask(mask_path)'; 

% 네트워크 재조립 (Inference용) [cite: 2026-02-10]
lgraph = layerGraph(net);
newInput = imageInputLayer([H W 1], 'Name', 'input', 'Normalization', 'none');
lgraph = replaceLayer(lgraph, 'input', newInput);
netFull = assembleNetwork(lgraph);

%% [3] ON/OFF 전수 복원 및 Raw 평균 계산
% 평균 저장용 변수들 (Raw vs Denoised)
avg_raw_on   = zeros(H, W, 'single');
avg_raw_off  = zeros(H, W, 'single');
avg_deno_on  = zeros(H, W, 'single');
avg_deno_off = zeros(H, W, 'single');

fprintf('1. ON 데이터셋 전수 처리 중 (%d frames)...\n', n_frames);
for i = 1:n_frames
    img = squeeze(single(h5read(h5_file_path, on_dataset, [i 1 1], [1 H W])))'; % Transpose [cite: 2026-02-11]
    
    % Raw 누적
    avg_raw_on = avg_raw_on + img / n_frames;
    
    % N2N 복원 및 누적 [cite: 2026-02-10, 2026-02-11]
    norm_img = log(img + 1) / maxLogVal;
    deno = exp(predict(netFull, reshape(norm_img, [H,W,1,1])) * maxLogVal) - 1;
    avg_deno_on = avg_deno_on + squeeze(deno) / n_frames;
end

fprintf('2. OFF 데이터셋 전수 처리 중 (%d frames)...\n', m_frames);
for j = 1:m_frames
    img = squeeze(single(h5read(h5_file_path, off_dataset, [j 1 1], [1 H W])))'; % Transpose [cite: 2026-02-11]
    
    % Raw 누적
    avg_raw_off = avg_raw_off + img / m_frames;
    
    % N2N 복원 및 누적
    norm_img = log(img + 1) / maxLogVal;
    deno = exp(predict(netFull, reshape(norm_img, [H,W,1,1])) * maxLogVal) - 1;
    avg_deno_off = avg_deno_off + squeeze(deno) / m_frames;
end

% 모든 평균 이미지에 마스크(mask==1) 적용 [cite: 2026-02-11]
avg_raw_on(mask == 1)   = 0;  avg_raw_off(mask == 1)  = 0;
avg_deno_on(mask == 1)  = 0;  avg_deno_off(mask == 1) = 0;

%% [4] 1D Azimuthal Integration (Raw & Denoised)
fprintf('3. Azimuthal Integration 수행 중...\n');
detectors = py.importlib.import_module('pyFAI.detectors');
azimuthal = py.importlib.import_module('pyFAI.azimuthalIntegrator');
ai = AI(poni_path, detectors, azimuthal);
num_q = 1024;
ai_opt = pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask);

% 4가지 1D Curve 추출
[c_raw_on, q] = ai_start(ai, avg_raw_on, ai_opt, num_q);
[c_raw_off, ~] = ai_start(ai, avg_raw_off, ai_opt, num_q);
[c_deno_on, ~] = ai_start(ai, avg_deno_on, ai_opt, num_q);
[c_deno_off, ~] = ai_start(ai, avg_deno_off, ai_opt, num_q);

%% [5] 정규화 및 차분($\Delta S$) 계산 (q = 1~6)
n_idx = (q >= 1) & (q <= 6);

% Raw 차분 계산
ds_raw = (c_raw_on / mean(c_raw_on(n_idx))) - (c_raw_off / mean(c_raw_off(n_idx)));

% N2N 차분 계산
ds_deno = (c_deno_on / mean(c_deno_on(n_idx))) - (c_deno_off / mean(c_deno_off(n_idx)));

%% [6] 시각화 및 비교 (동시 Plot)
figure('Name', 'Raw vs N2N Integration Comparison', 'Position', [50, 50, 1200, 800]);
t = tiledlayout(2, 1, 'TileSpacing', 'compact');

% --- Top: 1D Intensity Profiles ---
nexttile;
plot(q, c_raw_on, 'k:', 'LineWidth', 1, 'DisplayName', 'Raw ON'); hold on;
plot(q, c_raw_off, 'b:', 'LineWidth', 1, 'DisplayName', 'Raw OFF');
plot(q, c_deno_on, 'r-', 'LineWidth', 1.5, 'DisplayName', 'N2N Denoised ON');
plot(q, c_deno_off, 'g-', 'LineWidth', 1.5, 'DisplayName', 'N2N Denoised OFF');
grid on; ylabel('Intensity (a.u.)');
legend('Location', 'northeast', 'NumColumns', 2);
title('1D Intensity Profiles: Raw vs N2N Denoised');

% --- Bottom: Difference Signals (\DeltaS) ---
nexttile;
plot(q, ds_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Raw \DeltaS'); hold on;
plot(q, ds_deno, 'r-', 'LineWidth', 2, 'DisplayName', 'N2N Denoised \DeltaS');
grid on; xlabel('q (\AA^{-1})'); ylabel('\DeltaS (Normalized)');
legend('Location', 'northeast');
title('Comparison of Difference Signals (\DeltaS = ON - OFF)');

sgtitle('Noise2Noise Total Performance Evaluation: Final Comparison');

%% [보조 함수]
function ai = AI(poni_path, detectors, azimuthal)
    poni = read_poni(poni_path); 
    det_shape = py.tuple({int32(2048), int32(2048)});
    detector = detectors.Detector(poni.PixelSize1, poni.PixelSize2, pyargs('max_shape', det_shape));
    ai_params = pyargs('dist', poni.Distance, 'poni1', poni.Poni1, 'poni2', poni.Poni2, ...
        'rot1', poni.Rot1, 'rot2', poni.Rot2, 'rot3', poni.Rot3, ...
        'detector', detector, 'wavelength', poni.Wavelength);
    ai = azimuthal.AzimuthalIntegrator(ai_params);
end

function [curve, q] = ai_start(ai, Data, ai_option, num_q)
    res = ai.integrate1d(double(Data), int32(num_q), ai_option);
    curve = double(res{2}.tolist())';
    q = double(res{1}.tolist())';
end