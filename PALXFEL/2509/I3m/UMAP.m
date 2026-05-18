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

%%
for i = 1  %1 : length(DIR)
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

for j = 50  % : 15%: 200 %1 : length(timedelays)
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
    plot(q_vector,mean(DS,2));

    hold on;
    I0 = I0(on_idx);
    oh1 = oh1(on_idx);
    oh2 = oh2(on_idx);
    pos = pos(on_idx);
end












end



%%
function folder_list = DIR_lst(scan_path,header)
items = dir(fullfile(scan_path,header+'_*'));
is_folder = [items.isdir];
folders = items(is_folder);
folder_list = {folders.name}';
folder_list = string(folder_list);
end