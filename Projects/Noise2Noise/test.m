clc; clear; close all;

%% [1] 데이터 및 마스크 로드
h5_file_path = "Z:\KEK_2512_h5\results\KEK_2512\I3_MeOH\run4\raw\images\image_pairs.h5"; 
dataset_name = "/delays/delay_1.00e-09/on"; 
mask_path    = "Z:\KEK_2512_h5\run1_neg_001-mask.edf";

% 첫 번째 프레임 로드 [Frame=1, H=2048, W=2048]
fprintf('이미지 로드 중...\n');
img_raw = h5read(h5_file_path, dataset_name, [1 1 1], [1 2048 2048]);
img = squeeze(single(img_raw))'; % 2D 행렬로 변환

% 마스크 로드 (마스킹 영역 = 1)
mask = read_mask(mask_path);

%% [2] 시각화 및 축 확인
figure('Name', 'Axis Orientation Check', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.6]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact');

% --- 1. 원본 이미지 (Log Scale) ---
nexttile;
imagesc(log10(img + 1)); 
axis image; % 이미지 비율 유지 (축 확인의 핵심!)
colorbar;
title('1. Raw Image (Log10)');
xlabel('Width (Column)'); ylabel('Height (Row)');

% --- 2. 마스크 확인 ---
nexttile;
imagesc(mask); 
axis image;
title('2. Mask (1=Masked)');
xlabel('W'); ylabel('H');

% --- 3. 오버레이 확인 (이미지 위에 마스크 겹치기) ---
nexttile;
% 이미지를 배경으로 깔고, 마스크가 1인 곳을 붉은색으로 표시
imagesc(log10(img + 1)); hold on;
axis image;
h = imagesc(mask);
h.AlphaData = mask * 0.5; % 마스크 영역만 50% 투명도로 겹침
colormap(gca, [hot(256); 1 0 0]); % 기본은 hot, 마스크는 빨강
title('3. Overlay Check');
xlabel('W'); ylabel('H');

sgtitle('KEK Data Axis Alignment: Check if red areas match the gaps/beamstop');

fprintf('축 확인 가이드:\n');
fprintf('1. 이미지의 빔스탑 위치와 마스크의 붉은 영역이 일치합니까?\n');
fprintf('2. 가로/세로 비율(2048x2048)이 정사각형으로 보입니까?\n');