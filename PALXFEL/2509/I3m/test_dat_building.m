% ==========================================
% 1. 설정 및 경로 (MATLAB 원본 기준 복원)
% ==========================================
scan_path = 'Z:\PAL_XFEL_2509\ssd1\scan';
header = 'I3m_MeCN';
X_center = 12.699;
repetitionRate = 60;
UpsideDown = false;
auc_range = [1.5, 4.5];
q_target_range = [1, 7];

% ==========================================
% 2. 타임 딜레이 목록 정제 (부동소수점 오차 우회)
% ==========================================
disp('런 폴더 스캔 및 타임 딜레이 추출 중...');
run_folders = dir(fullfile(scan_path, sprintf('%s_*', header)));
all_delays = [];

for i = 1:length(run_folders)
    rf_path = fullfile(run_folders(i).folder, run_folders(i).name);
    sinfo = fullfile(rf_path, 'scanInfo', 'scanInfo.h5');
    if ~isfile(sinfo), continue; end
    
    try
        % 2D 스캔 제외 (m2 모터 존재 여부 확인)
        info = h5info(sinfo, '/run/scan00001/motor');
        motor_names = {info.Datasets.Name};
        if ismember('m2', motor_names), continue; end
        
        delays = h5read(sinfo, '/run/scan00001/motor/m1');
        all_delays = [all_delays; delays(:)];
    catch
        continue;
    end
end

% MATLAB의 강력한 허용 오차 기반 고유값 추출 (1e-4 오차 허용)
valid_delays = uniquetol(all_delays, 1e-4);

disp('==============================');
disp('  유효 타임 딜레이 목록');
disp('==============================');
for i = 1:length(valid_delays)
    fprintf('  [%2d] : %8.3f ps\n', i, valid_delays(i));
end

idx = input('\n분석할 번호를 입력하세요: ');
target_delay = valid_delays(idx);

% ==========================================
% 3. 데이터 로드 및 전처리
% ==========================================
disp(['데이터 추출 시작: ', num2str(target_delay), ' ps']);
all_ds = [];
q_vec = [];

for r_idx = 1:length(run_folders)
    rf_path = fullfile(run_folders(r_idx).folder, run_folders(r_idx).name);
    sinfo = fullfile(rf_path, 'scanInfo', 'scanInfo.h5');
    if ~isfile(sinfo), continue; end
    
    try
        t_delays = h5read(sinfo, '/run/scan00001/motor/m1');
    catch
        continue;
    end
    
    % 타겟 딜레이 매칭 (Tolerance 적용)
    matched_idx = find(abs(t_delays - target_delay) < 1e-4);
    if isempty(matched_idx), continue; end
    
    pulse_files = dir(fullfile(rf_path, 'pulseInfo', '*.h5'));
    if isempty(pulse_files), continue; end
    
    for m = 1:length(matched_idx)
        f_idx = matched_idx(m);
        if f_idx > length(pulse_files), continue; end
        
        p_path = fullfile(pulse_files(f_idx).folder, pulse_files(f_idx).name);
        h5name = pulse_files(f_idx).name;
        
        try
            % Shot 목록 읽기 및 PID 파싱
            info = h5info(p_path);
            shots = {info.Groups.Name}; 
            pids = zeros(1, length(shots));
            for s = 1:length(shots)
                parts = strsplit(shots{s}, '_');
                pids(s) = str2double(parts{end});
            end
            
            % 가속기 진단 데이터 로드
            data_paths = {'eh1qbpm1_totalsum', 'ohqbpm1_totalsum', 'ohqbpm2_totalsum', 'eh1oxc_pos'};
            diag_vals = zeros(length(data_paths), length(shots));
            for dp_idx = 1:length(data_paths)
                dp_file = fullfile(rf_path, data_paths{dp_idx}, h5name);
                for s = 1:length(shots)
                    diag_vals(dp_idx, s) = h5read(dp_file, [shots{s}]);
                end
            end
            
            % 아웃라이어 제거 (MATLAB 내장 isoutlier 활용)
            valid_mask = ~any(isoutlier(diag_vals, 'median', 2), 1);
            if ~any(valid_mask), continue; end
            
            pids_eff = pids(valid_mask);
            valid_shots = shots(valid_mask);
            
            % 강도 데이터(AI) 로드
            ai_file = fullfile(rf_path, 'eh1rayMXAI_int', h5name);
            ai_data_temp = h5read(ai_file, [valid_shots{1}]);
            ai_data = zeros(length(ai_data_temp), length(valid_shots));
            for s = 1:length(valid_shots)
                ai_data(:, s) = h5read(ai_file, [valid_shots{s}]);
            end
            
            % q_vector 생성 (최초 1회)
            if isempty(q_vec)
                tth_file = fullfile(rf_path, 'eh1rayMXAI_tth', h5name);
                tth_deg = h5read(tth_file, [valid_shots{1}]);
                tth = deg2rad(tth_deg);
                q_vec = 4 * pi * sin(tth / 2) / (12.3984 / X_center);
            end
            
            % ON/OFF 분리
            is_on = mod(pids_eff, 360 / (repetitionRate / 2)) == 0;
            if UpsideDown, is_on = ~is_on; end
            
            if sum(is_on) < 2 || sum(~is_on) < 2, continue; end
            
            % AUC 정규화
            q_mask_auc = (q_vec > auc_range(1)) & (q_vec < auc_range(2));
            auc_vals = sum(ai_data(q_mask_auc, :), 1);
            
            ai_on = ai_data(:, is_on) ./ auc_vals(is_on);
            ai_off = ai_data(:, ~is_on) ./ auc_vals(~is_on);
            
            % KNN 기반 샷 페어링 (MATLAB 내장 knnsearch 활용)
            auc_on = auc_vals(is_on)';
            auc_off = auc_vals(~is_on)';
            match_idx = knnsearch(auc_off, auc_on);
            
            ds_batch = ai_on - ai_off(:, match_idx);
            
            % 유효 q 범위 추출 및 누적
            q_mask_final = (q_vec >= q_target_range(1)) & (q_vec <= q_target_range(2));
            all_ds = [all_ds, ds_batch(q_mask_final, :)];
            
        catch
            continue;
        end
    end
end

if isempty(all_ds)
    disp('유효한 데이터가 없습니다.');
    return;
end

raw_ds = all_ds'; % (Shot x q_vec) 형태로 전치
q_vec_final = q_vec((q_vec >= q_target_range(1)) & (q_vec <= q_target_range(2)));

