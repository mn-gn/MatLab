% TRXL 데이터 전처리 및 변수 추출 최적화 스크립트 (진행상황 출력 포함)
clear; clc;

% --- 사용자 설정 파라미터 ---
start_run = 1;
end_run = 100;
base_dir = "Z:\PAL_XFEL_2509\ssd2\ue_250907_FXL\scratch\resultsTRXL\I3m_MeCN"; 
save_filename = 'C:\Users\ming0\Desktop\MATLAB\Reaction_path\TRXL_Processed_Data.mat'; % 저장될 파일 이름
q_min = 1.0;
q_max = 7.0;

% --- 데이터 수집 초기화 ---
all_time_delays = [];
all_weights = [];
all_DS_matrix = [];
q_vector = []; 
idx = 1;

fprintf('데이터 수집을 시작합니다...\n');
fprintf('========================================\n');

for i = start_run:end_run
    run_name = sprintf('run%03d', i);
    target_dir = fullfile(base_dir, run_name, 'DiffAve');
    
    if ~exist(target_dir, 'dir')
        fprintf('[경고] %s 폴더가 존재하지 않아 건너뜁니다.\n', run_name);
        continue; 
    end
    
    files = dir(fullfile(target_dir, 'diff_av_*.txt'));
    
    for j = 1:length(files)
        filename = files(j).name;
        tokens = regexp(filename, 'diff_av_([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', 'tokens');
        if isempty(tokens), continue; end
        
        filepath = fullfile(target_dir, filename);
        
        % 1. 가중치(number_of_images) 추출
        fid = fopen(filepath, 'r');
        first_line = fgetl(fid);
        fclose(fid);
        
        weight_tokens = regexp(first_line, 'number_of_images:\s*(\d+)', 'tokens');
        weight_val = 1; % 기본값
        if ~isempty(weight_tokens)
            weight_val = str2double(weight_tokens{1}{1});
        end
        
        % 2. 데이터 로드 및 q 범위 필터링
        data = importdata(filepath);
        if isstruct(data)
            data_mat = data.data;
        else
            data_mat = data;
        end
        
        valid_idx = (data_mat(:, 1) >= q_min) & (data_mat(:, 1) <= q_max);
        
        if isempty(q_vector)
            q_vector = data_mat(valid_idx, 1);
        end
        
        filtered_ds = data_mat(valid_idx, 2);
        if length(filtered_ds) ~= length(q_vector), continue; end
        
        % 3. 임시 배열에 누적
        all_time_delays(idx, 1) = str2double(tokens{1}{1});
        all_weights(idx, 1) = weight_val;
        all_DS_matrix(idx, :) = filtered_ds';
        idx = idx + 1;
    end
    fprintf('%s 데이터 읽기 완료\n', run_name);
end

fprintf('----------------------------------------\n');
fprintf('총 수집된 데이터 파일 수: %d개\n', idx-1);

% --- 가중평균 및 정렬 ---
if ~isempty(all_time_delays)
    fprintf('동일한 Time Delay에 대해 가중평균을 계산 중입니다...\n');
    unique_delays = unique(all_time_delays);
    num_unique = length(unique_delays);

    time_vector = zeros(num_unique, 1);
    dS_matrix = zeros(num_unique, length(q_vector));

    for k = 1:num_unique
        t_idx = (all_time_delays == unique_delays(k));
        target_weights = all_weights(t_idx);
        
        % 동일 시간대 데이터의 가중평균 계산
        dS_matrix(k, :) = sum(all_DS_matrix(t_idx, :) .* target_weights, 1) / sum(target_weights);
        time_vector(k, 1) = unique_delays(k);
    end

    % 시간 순 정렬
    [time_vector, sort_idx] = sort(time_vector);
    dS_matrix = dS_matrix(sort_idx, :);

    % --- 데이터 저장 및 변수 삭제 (최적화) ---
    save(save_filename, 'q_vector', 'time_vector', 'dS_matrix');
    clearvars -except q_vector time_vector dS_matrix;
    
    fprintf('가중평균 완료. 최종 고유 Time Delay 수: %d개\n', length(time_vector));
    fprintf('========================================\n');
    disp('✅ 최종 데이터 변수 추출 완료: q_vector, time_vector, dS_matrix');
else
    fprintf('처리할 데이터가 없습니다.\n');
end