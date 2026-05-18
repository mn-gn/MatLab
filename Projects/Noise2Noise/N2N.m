clc; clear; close all;

%% [1] 데이터 로드 및 파라미터 설정
h5_file_path = "Z:\KEK_2512_h5\results\KEK_2512\I3_MeOH\run4\raw\images\image_pairs.h5"; 
dataset_name = "/delays/delay_1.00e-09/on"; 
mask_path    = "Z:\KEK_2512_h5\run1_neg_001-mask.edf";
model_save_name = 'N2N_Xray_Model.mat';

% 이미지 구조: [374(Frames) x 2048(H) x 2048(W)]
fprintf('1. H5 데이터 로드 중...\n');
raw_data = single(h5read(h5_file_path, dataset_name));
[num_total, H, W] = size(raw_data); 

%% [2] 전처리: Outlier 제거 (Bad Shot Removal)
fprintf('2. 전처리: Outlier 제거 중...\n');
mean_intensities = squeeze(mean(mean(raw_data, 2), 3));
threshold_low = mean(mean_intensities) - 2*std(mean_intensities);
valid_idx = find(mean_intensities > threshold_low);

clean_stack = raw_data(valid_idx, :, :);
num_clean = length(valid_idx);
fprintf('   -> 유효 프레임: %d / %d\n', num_clean, num_total);

%% [3] 전처리: Log-Scale 및 Normalization
fprintf('3. 전처리: Log 변환 및 정규화 중...\n');
data_log = log(clean_stack + 1);
maxLogVal = max(data_log(:)); % 이 값은 나중에 복원 시 꼭 필요합니다.
data_norm = data_log / maxLogVal;

%% [4] 마스킹(Mask=1)을 고려한 유효 패치 추출
fprintf('4. 마스크 필터링(Mask=1 제외) 및 패치 추출 중...\n');
mask = read_mask(mask_path); % EDF 마스크 로드 (마스킹 영역 = 1)

patchSize = 64;
numPatchesPerFrame = 40; % 데이터 양에 따라 조절하세요.
trainInput = [];
trainTarget = [];

% Noise2Noise Pair 구성 (Frame i -> i+1)
for i = 1:2:num_clean-1
    img1 = squeeze(data_norm(i, :, :))';
    img2 = squeeze(data_norm(i+1, :, :))';
    
    count = 0;
    while count < numPatchesPerFrame
        y = randi([1, H - patchSize + 1]);
        x = randi([1, W - patchSize + 1]);
        
        % 해당 영역의 마스크 추출
        maskPatch = mask(y:y+patchSize-1, x:x+patchSize-1);
        
        % [핵심] 마스킹 된 영역(1)이 패치 내에 하나라도 있으면 패스
        if any(maskPatch(:) == 1)
            continue; 
        end
        
        % 유효한 영역만 학습 데이터에 추가
        trainInput(:,:,1,end+1) = img1(y:y+patchSize-1, x:x+patchSize-1);
        trainTarget(:,:,1,end+1) = img2(y:y+patchSize-1, x:x+patchSize-1);
        count = count + 1;
    end
end
fprintf('   -> 총 %d개의 유효 패치로 학습을 준비합니다.\n', size(trainInput, 4));

%% [5] 네트워크 정의 및 학습
fprintf('5. 모델 학습 시작...\n');

layers = [
    imageInputLayer([patchSize patchSize 1], 'Name', 'input', 'Normalization', 'none')
    convolution2dLayer(3, 64, 'Padding', 'same')
    reluLayer()
    convolution2dLayer(3, 64, 'Padding', 'same')
    reluLayer()
    convolution2dLayer(3, 64, 'Padding', 'same')
    reluLayer()
    convolution2dLayer(3, 1, 'Padding', 'same')
    regressionLayer()
];

options = trainingOptions('adam', ...
    'MaxEpochs', 30, ...
    'MiniBatchSize', 128, ...
    'InitialLearnRate', 1e-3, ...
    'Plots', 'training-progress', ...
    'Shuffle', 'every-epoch', ...
    'Verbose', false);

net = trainNetwork(trainInput, trainTarget, layers, options);

%% [6] 모델 저장
fprintf('6. 모델 및 파라미터 저장 중...\n');
save(model_save_name, 'net', 'maxLogVal');
fprintf('학습 완료! 모델이 %s에 저장되었습니다.\n', model_save_name);