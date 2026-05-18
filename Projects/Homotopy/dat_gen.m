% 호모토피(Optimal Transport) 및 UMAP 분석을 위한 TRXL 데이터 전처리 스크립트
clear; clc;

% --- 사용자 설정 파라미터 ---
start_run = 1;
end_run = 40;
base_dir = "Z:\PAL_XFEL_2509\ssd2\ue_250907_FXL\scratch\resultsTRXL\Mn_MeOH"; 
q_min = 1.0;
q_max = 7.0;
save_dir = "C:\Users\ming0\Desktop\homotopy\Mn_MeOH";

if ~exist(save_dir, 'dir'), mkdir(save_dir); end

% 정규표현식: 숫자 부분만 추출
pattern = 'diff_av_([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';

% 모든 파일 데이터를 임시로 담을 배열 초기화
all_time_delays = [];
all_weights = [];
all_DS_matrix = [];
q_axis = []; % 호모토피 계산의 '기하학적 거리(Cost matrix)' 기준이 될 필수 축

idx = 1;
fprintf('데이터를 수집 중입니다...\n');

for i = start_run:end_run
    run_name = sprintf('run%03d', i);
    target_dir = fullfile(base_dir, run_name, 'DiffAve');
    
    if ~exist(target_dir, 'dir'), continue; end
    
    files = dir(fullfile(target_dir, 'diff_av_*.txt'));
    
    for j = 1:length(files)
        filename = files(j).name;
        tokens = regexp(filename, pattern, 'tokens');
        
        if isempty(tokens)
            continue;
        end
        
        delay_val = str2double(tokens{1}{1});
        
        try
            filepath = fullfile(target_dir, filename);
            
            % 💡 [수정] 파일의 첫 번째 줄을 읽어 number_of_images 추출
            fid = fopen(filepath, 'r');
            first_line = fgetl(fid);
            fclose(fid);
            
            weight_val = 1; % 기본 가중치 (추출 실패 시)
            weight_tokens = regexp(first_line, 'number_of_images:\s*(\d+)', 'tokens');
            if ~isempty(weight_tokens)
                weight_val = str2double(weight_tokens{1}{1});
            end
            
            data = importdata(filepath);
            
            % importdata 반환값 처리 (버전/파일형식 호환성)
            if isstruct(data)
                data_mat = data.data;
            else
                data_mat = data;
            end
            
            q_col = data_mat(:, 1);
            ds_col = data_mat(:, 2);
            
            valid_idx = (q_col >= q_min) & (q_col <= q_max);
            filtered_q = q_col(valid_idx);
            filtered_ds = ds_col(valid_idx);
            
            % 첫 번째 파일에서 공통 q-axis 추출
            if isempty(q_axis)
                q_axis = filtered_q(:)';
            end
            
            % q-axis 길이가 다르면 에러 방지를 위해 스킵
            if length(filtered_ds) ~= length(q_axis)
                fprintf('경고: q 범위 불일치 (Skip) - %s\n', filename);
                continue;
            end
            
            % 행렬에 전체 데이터 임시 누적 (이후 가중평균 처리를 위함)
            all_time_delays(idx, 1) = delay_val;
            all_weights(idx, 1) = weight_val;
            all_DS_matrix(idx, :) = filtered_ds(:)';
            
            idx = idx + 1;
        catch ME
            fprintf('파일 읽기 오류 무시: %s (%s)\n', filename, ME.message);
        end
    end
    fprintf('%s 읽기 완료\n', run_name);
end

fprintf('총 수집된 데이터 파일 수: %d개\n', idx-1);

% --- 💡 [수정] Time Delay 기준으로 가중평균 계산 ---
if ~isempty(all_time_delays)
    fprintf('동일한 타임 딜레이에 대해 이미지 수를 기반으로 가중평균을 계산합니다...\n');
    
    unique_delays = unique(all_time_delays);
    num_unique = length(unique_delays);
    
    final_time_delays = zeros(num_unique, 1);
    final_DS_matrix = zeros(num_unique, length(q_axis));
    
    for k = 1:num_unique
        current_delay = unique_delays(k);
        
        % 현재 타임 딜레이에 해당하는 데이터 인덱스 찾기
        target_idx = (all_time_delays == current_delay);
        
        target_data = all_DS_matrix(target_idx, :);
        target_weights = all_weights(target_idx);
        
        total_weight = sum(target_weights);
        
        if total_weight > 0
            % 가중평균 계산: sum(데이터 * 가중치) / 전체 가중치 합
            weighted_sum = sum(target_data .* target_weights, 1);
            final_DS_matrix(k, :) = weighted_sum / total_weight;
        else
            % 예외 처리: 가중치 합이 0인 경우 단순 평균
            final_DS_matrix(k, :) = mean(target_data, 1);
        end
        
        final_time_delays(k, 1) = current_delay;
    end
    
    % --- 데이터 정렬 (Time Delay 기준) ---
    [sorted_time_delays, sort_idx] = sort(final_time_delays);
    sorted_DS_matrix = final_DS_matrix(sort_idx, :);
    
    fprintf('가중평균 완료. 최종 산출된 고유 Time Delay 수: %d개\n', length(sorted_time_delays));
    
    % --- 🚀 파이썬 호모토피 학습을 위한 최적화 파일 내보내기 ---
    
    % 1. MAT 파일 저장 (Python scipy.io.loadmat 으로 즉시 로드 가능)
    % ※ 런 넘버(run_nums)는 여러 런이 혼합되었으므로 제거했습니다.
    save_path_mat = fullfile(save_dir, 'TRXL_Homotopy_Input.mat');
    save(save_path_mat, 'q_axis', 'sorted_DS_matrix', 'sorted_time_delays');
    fprintf('\n✅ Python 연동용 MAT 파일 저장 완료: %s\n', save_path_mat);
    
    % 2. CSV 파일 저장 (Pandas 로딩용 보조 데이터)
    writematrix([q_axis; sorted_DS_matrix], fullfile(save_dir, 'TRXL_DS_matrix.csv'));
    writematrix(sorted_time_delays, fullfile(save_dir, 'TRXL_time_delays.csv'));
    fprintf('✅ CSV 데이터셋 저장 완료!\n');
else
    fprintf('저장할 데이터가 없습니다.\n');
end