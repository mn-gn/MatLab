clc; clear; close all;
%%
% DS(q), q, time delay, eh1qbpm1, ohqbpm1, ohqbpm2, 

scan_path = "Z:\PAL_XFEL_2509\ssd1\scan";
header = "I3m_MeCN";
X_center = 12.699;
repetitionRate = 60;
UpsideDown = false;
DIR = DIR_lst(scan_path,header);
auc_range = [1.5, 4.5];

q_target_range = [1, 7]; % 남겨둘 유효 q 범위 (데이터에 맞게 수정하세요)
save_dir = 'Z:\mngn_test';



%%
for i = 1: length(DIR)
    save_name = sprintf("I3m_Dataset%.3d.h5",i);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
export_filename = fullfile(save_dir, save_name);

if isfile(export_filename)
    delete(export_filename); % 덮어쓰기 방지 (새로 시작할 때마다 초기화)
end

is_h5_created = false; % HDF5 구조가 생성되었는지 확인하는 플래그
current_shot_idx = 1;  % 현재까지 파일에 쓰인 총 샷의 위치(인덱스) 추적
chunk_size = 1024;     % I/O 속도 향상을 위한 데이터 청크 단위
    run_path = fullfile(scan_path,DIR(i));
    scaninfo_file = fullfile(run_path,'scanInfo','scanInfo.h5');
%% 2D SCAN MODE 판별
    mode = '2dscan';
    try
        h5read(scaninfo_file,"/run/scan00001/motor/m2");
    catch
        mode = 'normal';
    end

    if mode == "2dscan"
        fprintf('run %d: skip the 2dscan mode \n',i)
        continue
    end
    
    clear scaninfo_file
%% TIME DELAY LOOP
timedelays = h5read(run_path+"\scanInfo\scanInfo.h5","/run/scan00001/motor/m1");
timedelays = sort(timedelays);

pulseinfo_file = fullfile(run_path,'pulseInfo','*.h5');
pulseinfo_file = dir(pulseinfo_file);

for j = 1: length(timedelays)
    timedelay = timedelays(j);                                              % ps
    h5name = pulseinfo_file(j).name;                                        % "001_001_00i.h5"
    pulseinfo_path = fullfile(pulseinfo_file(j).folder, h5name);
    pulseinfo = h5info(pulseinfo_path);

    pulsename = {pulseinfo.Datasets.Name}';

    pulseinfo = split({pulseinfo.Datasets.Name}','_');
    
    pid = (str2double(pulseinfo(:,2)));

    %% GET q VECTOR
    tth_path = fullfile(run_path,"eh1rayMXAI_tth",h5name);
    tth = deg2rad(h5read(tth_path,"/"+pulsename{1}));
    lambda = 12.3984/X_center;
    q_vector = 4*pi*sin(tth/2)/lambda;

    num_pulse = length(pulsename);
    idx_all = (1 : num_pulse)';
    
    %% OUTLIER (eh1qbpm1, ohqbpm1, ohqbpm2, eh1oxc_pos)

    I0_path = fullfile(run_path,"eh1qbpm1_totalsum",h5name);
    oh1_path = fullfile(run_path,"ohqbpm1_totalsum",h5name);
    oh2_path = fullfile(run_path,"ohqbpm2_totalsum",h5name);
    pos_path = fullfile(run_path,"eh1oxc_pos",h5name);

    I0 = zeros(1,num_pulse);
    oh1 = zeros(1,num_pulse);
    oh2 = zeros(1,num_pulse);
    pos = zeros(1,num_pulse);

    for k = 1 : num_pulse
        I0(k) = h5read(I0_path,"/"+pulsename{k});
        oh1(k) = h5read(oh1_path,"/"+pulsename{k});
        oh2(k) = h5read(oh2_path,"/"+pulsename{k});
        pos(k) = h5read(pos_path,"/"+pulsename{k});
    end

    idx_eff = idx_all(~isoutlier(I0)&~isoutlier(oh1)&~isoutlier(oh2)&~isoutlier(pos));
    num_eff = length(idx_eff);
    
    I0 = I0(idx_eff);
    oh1 = oh1(idx_eff);
    oh2 = oh2(idx_eff);
    pos = pos(idx_eff);

    pid = pid(idx_eff);

    filtered_idx = 1 : num_eff;

    %% GET AI
    AI_path = fullfile(run_path,"eh1rayMXAI_int",h5name);

    AI = zeros(length(q_vector), num_eff);
    for k = 1 : num_eff
        
        AI(:,k) = h5read(AI_path,"/"+pulsename{idx_eff(k)});
    end
    

    %% ON OFF indexing
    divisor = 360/(repetitionRate/2);
    remainders = mod(pid,divisor);
    is_on = (remainders == 0);
    if UpsideDown, is_on = ~is_on; end

    on_idx = filtered_idx(is_on);
    off_idx = filtered_idx(~is_on);



    %% AUC
    
    logical_roi = q_vector > auc_range(1) & q_vector < auc_range(2);
    q_roi = q_vector(logical_roi);
    dq = diff(q_roi);
    W = zeros(length(q_roi),1);
                
    if ~isempty(dq)
    W(1) = dq(1)/2;
    W(2:end-1) = (dq(1:end-1) + dq(2:end)) / 2;
    W(end) = dq(end) / 2;
    else
    W(1) = 1;
    end
                
    data_roi = AI(logical_roi,:);
    auc = data_roi'*W;
    %% SHOT PAIRING
    

    AI_on = AI(:,on_idx)./auc(on_idx)';
    AI_off = AI(:,off_idx)./auc(off_idx)';

    match_idx = knnsearch(auc(off_idx),auc(on_idx));
    matched_off = AI_off(:,match_idx);

    DS = AI_on - matched_off;
    % plot(q_vector,mean(DS,2));
    % 
    % hold on;
%% --- [추가] 2. q 범위 잘라내기 및 실시간 HDF5 이어쓰기 (Append) ---
    % 1) 타겟 q 범위에 해당하는 인덱스 찾기
    valid_q_idx = q_vector >= q_target_range(1) & q_vector <= q_target_range(2);
    DS_truncated = DS(valid_q_idx, :); % [q_len x num_shots]
    q_vector_truncated = q_vector(valid_q_idx);
    num_on_shots = size(DS_truncated, 2);

    if num_on_shots > 0
        q_len = length(q_vector_truncated);

        % 2) 최초 1회만 HDF5 데이터셋 구조 생성 (무한 차원 Inf 적용)
        if ~is_h5_created
            disp('최초 HDF5 뼈대 생성 중...');
            % q_vector는 고정 크기이므로 바로 생성 후 씁니다.
            h5create(export_filename, '/data/q_vector', [q_len, 1]);
            h5write(export_filename, '/data/q_vector', q_vector_truncated);

            % 나머지 데이터는 샷 개수가 늘어나야 하므로 Inf로 설정하고 ChunkSize를 줍니다.
            h5create(export_filename, '/data/signal', [q_len, Inf], 'ChunkSize', [q_len, chunk_size]);
            h5create(export_filename, '/metadata/timedelay', [1, Inf], 'ChunkSize', [1, chunk_size]);
            h5create(export_filename, '/metadata/run_id', [1, Inf], 'ChunkSize', [1, chunk_size]);
            h5create(export_filename, '/metadata/I0', [1, Inf], 'ChunkSize', [1, chunk_size]);
            h5create(export_filename, '/metadata/oh1', [1, Inf], 'ChunkSize', [1, chunk_size]);
            h5create(export_filename, '/metadata/oh2', [1, Inf], 'ChunkSize', [1, chunk_size]);
            h5create(export_filename, '/metadata/pos', [1, Inf], 'ChunkSize', [1, chunk_size]);
            
            is_h5_created = true; % 구조 생성 완료 플래그
        end

        % 3) 메타데이터 배열 준비 (가로 벡터 형태 1xN 으로 통일)
        delay_array  = repmat(timedelay, 1, num_on_shots);
        run_id_array = repmat(i, 1, num_on_shots); % i는 런 인덱스
        % reshape를 통해 확실하게 1xN 행벡터로 만들어 차원 에러 방지
        I0_array  = reshape(I0(on_idx), 1, []);
        oh1_array = reshape(oh1(on_idx), 1, []);
        oh2_array = reshape(oh2(on_idx), 1, []);
        pos_array = reshape(pos(on_idx), 1, []);

        % 4) HDF5 파일에 현재 루프의 데이터 이어쓰기 (시작 인덱스 지정)
        % h5write(파일명, 변수명, 쓸데이터, [행시작, 열시작], [행크기, 열크기])
        h5write(export_filename, '/data/signal', DS_truncated, [1, current_shot_idx], [q_len, num_on_shots]);
        h5write(export_filename, '/metadata/timedelay', delay_array, [1, current_shot_idx], [1, num_on_shots]);
        h5write(export_filename, '/metadata/run_id', run_id_array, [1, current_shot_idx], [1, num_on_shots]);
        h5write(export_filename, '/metadata/I0', I0_array, [1, current_shot_idx], [1, num_on_shots]);
        h5write(export_filename, '/metadata/oh1', oh1_array, [1, current_shot_idx], [1, num_on_shots]);
        h5write(export_filename, '/metadata/oh2', oh2_array, [1, current_shot_idx], [1, num_on_shots]);
        h5write(export_filename, '/metadata/pos', pos_array, [1, current_shot_idx], [1, num_on_shots]);

        % 5) 다음 이어쓰기를 위해 인덱스 업데이트
        current_shot_idx = current_shot_idx + num_on_shots;
    end
end
end
%% --- [추가] 3. 작업 종료 메시지 ---
% 이제 한 번에 병합하는 과정이 필요 없습니다. 이미 하드디스크에 완벽하게 쓰였습니다.
fprintf('\n데이터셋 실시간 저장 완료!\n');
fprintf('총 저장된 유효 샷 개수: %d 개\n', current_shot_idx - 1);
fprintf('저장 경로: %s\n', export_filename);

%%
function folder_list = DIR_lst(scan_path,header)
items = dir(fullfile(scan_path,header+'_*'));
is_folder = [items.isdir];
folders = items(is_folder);
folder_list = {folders.name}';
folder_list = string(folder_list);
end