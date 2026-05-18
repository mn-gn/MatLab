% TRXL 데이터 SVD 분석 및 UMAP 입력 데이터 생성 스크립트 (RSV & SRSV 비교)
clear; clc;

% --- 1. 설정 및 데이터 로드 ---
input_file = 'C:\Users\ming0\Desktop\MATLAB\Reaction_path\TRXL_Processed_Data.mat';
output_dir = 'C:\Users\ming0\Desktop\MATLAB\Reaction_path'; % UMAP 데이터 저장 폴더

% 출력 폴더가 없으면 생성
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 이전 스크립트에서 저장한 변수들 (q_vector, time_vector, dS_matrix) 로드
load(input_file);
fprintf('데이터 로드 완료: %s\n', input_file);

% --- 2. SVD (특이값 분해) 수행 ---
% 공간(q)을 행으로, 시간(time)을 열로 맞추기 위해 행렬 전치
A = dS_matrix'; 

fprintf('SVD 분해를 수행 중입니다...\n');
[U, S, V] = svd(A, 'econ');

% 특이값 추출
s_values = diag(S);

% RSV (Right Singular Vectors) 및 SRSV (Scaled Right Singular Vectors) 계산
RSV = V;
SRSV = V * S;

% UMAP 분석에 활용할 상위 주성분(PC) 개수 설정
num_components = 5; 
umap_input_RSV = RSV(:, 1:num_components);
umap_input_SRSV = SRSV(:, 1:num_components);

% --- 3. UMAP용 데이터 파일 저장 ---
save_filepath_mat = fullfile(output_dir, 'UMAP_Input_Data.mat');
save_filepath_csv_rsv = fullfile(output_dir, 'UMAP_Input_RSV.csv');
save_filepath_csv_srsv = fullfile(output_dir, 'UMAP_Input_SRSV.csv');

% MAT 파일에 모두 저장
save(save_filepath_mat, 'time_vector', 'umap_input_RSV', 'umap_input_SRSV', 'num_components', 's_values');

% CSV 파일 각각 저장 (첫 번째 열은 시간 벡터)
writematrix([time_vector, umap_input_RSV], save_filepath_csv_rsv);
writematrix([time_vector, umap_input_SRSV], save_filepath_csv_srsv);

fprintf('\n✅ UMAP용 데이터 저장 완료:\n');
fprintf(' - MAT 통합 파일: %s\n', save_filepath_mat);
fprintf(' - CSV 파일 (RSV): %s\n', save_filepath_csv_rsv);
fprintf(' - CSV 파일 (SRSV): %s\n\n', save_filepath_csv_srsv);

