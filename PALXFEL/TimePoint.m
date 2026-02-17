classdef TimePoint < handle

    properties
        index (1,1) int32           % timedelay index for the run
        DelayValue (1,1) double     % timedelay for the timepoint
        Unit (1,1) string = "ps"    % timdelay unit, default = "ps"
        q  (:,1) double             % q vecter
        AIfile (:,1) string         % Aximuthal integrated shot path(*.h5)
        shot (:,1) SingleShot       % shot information (TimeStamp, PID, AI_int, I0, pumpOn)
        paring_idx (:,1) int32      % Pump laser on and off paring index
        normfactor (:,1) double
    end
    
    %% methods
    methods
        
        function obj = TimePoint(index,timedelay,unit) % Constructor
            if nargin > 0
                obj.index = index;
                obj.DelayValue = timedelay;
                
            end
            if nargin > 2
                obj.Unit = unit;
            end
        end
        
        function obj = get_q(obj,config)
            h5path = sprintf("%s/eh1rayMXAI_tth",config.dir);
            namePattern = sprintf("%s/*.h5",h5path);
            d = dir(namePattern);
           
            filename = h5path + '/' + d(1).name;
            TS_PID = h5info(filename).Datasets(1).Name;
            tth = deg2rad(h5read(filename,'/'+string(TS_PID)));
            lambda = 12.398/config.XrayCenter;
            obj.q = 4*pi*sin(tth/2)/lambda;
            
        end

        function obj = AI_path(obj,config)
            AIfolder = sprintf("%s/eh1rayMXAI_int",config.dir);     
            
            d = dir(fullfile(AIfolder,sprintf("*%.3d.h5",obj.index)));
            obj.AIfile = string(fullfile({d.folder},{d.name})');
        end

        function obj = singleshot(obj, config)
            file_batches = cell(1,length(obj.AIfile));  % 파일별로 싱글샷을 뭉텡이로 묶기 위해 빈 셀 만듬
            divisor = 360/(config.repetitionRate/2);    % 후에 레이저 온오프 판단하기 위해서 디바이저로 pid를 나눌꺼임

            %% file loop(different scans for the timedelay)
            for pathidx = 1 : length(obj.AIfile)
                h5path = char(obj.AIfile(pathidx));

                try
                    h5 = h5info(h5path);
                    allnames = string({h5.Datasets.Name})';

                    num_shots = length(allnames);

                    if num_shots == 0
                        continue;
                    end

                    current_batch(num_shots,1) = SingleShot;

                    splitdata = split(allnames,'_');
                    TimeStamps = str2double(splitdata(:,1));
                    PIDs = str2double(splitdata(:,2));

                    remainders = mod(PIDs,divisor);
                    is_on = (remainders == 0);

                    if config.UpsideDown, is_on = ~is_on; end

                    for i = 1 : num_shots
                        subpath = "/" + allnames(i);
                        data = h5read(h5path, subpath);

                        current_batch(i) = SingleShot(TimeStamps(i),PIDs(i),data,is_on(i));
                    end

                    file_batches{pathidx} = current_batch;

                    clear current_batch;

                catch ME
                    fprintf('  Warning: 파일 읽기 실패 %s (%s)\n', currentPath, ME.message);
                end
            end

            if ~isempty(file_batches)
                obj.shot = vertcat(file_batches{:});
            else
                obj.shot = SingleShot.empty;
            end
        end
        
        function obj = get_I0(obj, config)
            % Safety Check: Return immediately if no shots are loaded
            if isempty(obj.shot), return; end

            % Define I0 folder path
            I0folder = fullfile(config.dir, 'eh1qbpm1_totalsum');
            filePattern = sprintf("*_%.3d.h5", obj.index);
            d = dir(fullfile(I0folder, filePattern));

            if isempty(d)
                fprintf('Warning: I0 file cannot be found (Index: %d)\n', obj.index);
                return;
            end

            % Initialize buffers for batch processing
            pid_batch = {};
            i0_batch = {};

            for i = 1 : length(d)
                filePath = fullfile(d(i).folder, d(i).name);

                try
                    h5 = h5info(filePath);
                    allnames = {h5.Datasets.Name}; % Cell array
                    num_data = length(allnames);

                    if num_data == 0, continue; end

                    temp_i0s = zeros(num_data, 1, 'double');
                    name_str = string(allnames)';
                    split_res = split(name_str, '_');
                    temp_pids = int64(str2double(split_res(:, 2)));

                    for j = 1 : num_data
                        temp_i0s(j) = h5read(filePath, "/" + allnames{j});
                    end

                    pid_batch{end+1} = temp_pids; 
                    i0_batch{end+1} = temp_i0s; 

                catch ME
                    fprintf('  Warning: I0 read error %s (%s)\n', d(i).name, ME.message);
                end

            end

            if isempty(pid_batch), return; end

            all_I0_PIDs = vertcat(pid_batch{:});
            all_I0_Vals = vertcat(i0_batch{:});
            my_PIDs = [obj.shot.PID];
            [Lia, Locb] = ismember(my_PIDs, all_I0_PIDs);

            if any(Lia)
                matched_vals = all_I0_Vals(Locb(Lia));
                val_cells = num2cell(matched_vals);
                [obj.shot(Lia).I0] = val_cells{:};
            end

            fprintf('  I0 Matching: %d / %d shots matched.\n', sum(Lia), length(obj.shot));
        end

        function obj = AUC_cal(obj, config)
            if isempty(obj.shot), return; end
            % 
            % logical_roi = find(obj.q>=config.AUC_range(1)&obj.q<=config.AUC_range(2));
            % 
            % q_roi = obj.q(logical_roi);
            % dq = diff(q_roi);
            % W = zeros(length(q_roi),1);
            % 
            % if ~isempty(dq)
            %     W(1) = dq(1)/2;
            %     W(2:end-1) = (dq(1:end-1) + dq(2:end)) / 2;
            %     W(end) = dq(end) / 2;
            % else
            %     W(1) = 1;
            % end

            all_data = [obj.shot.AI_int];
            % data_roi = all_data(logical_roi,:);
            % 
            % calculated_AUCs = data_roi' * W;
            % 
            % auc_cells = num2cell(calculated_AUCs);
            auc_cells = num2cell(obj.area_under_curve(all_data, obj.q, config.AUC_range));
            [obj.shot.AUC] = auc_cells{:};

            fprintf('  AUC calculated for %d shots (Vectorized).\n', length(obj.shot));
        end

        function obj = shot_paring(obj)
            if isempty(obj.shot), return; end

            obj.paring_idx = zeros(length(obj.shot),1);

            all_aucs = [obj.shot.AUC];
            all_laser = [obj.shot.pumpOn];

            on_idx = find(all_laser);
            off_idx = find(~all_laser);
            
            if isempty(on_idx) || isempty(off_idx)
                fprintf('  Warning: Cannot pair (Missing On or Off shots).\n');
                return;
            end

            on_vals = all_aucs(on_idx);
            off_vals = all_aucs(off_idx);

            matches_local = dsearchn(off_vals(:), on_vals(:));

            matched_idx = off_idx(matches_local);

            obj.paring_idx = matched_idx;
            fprintf('Raw Pairing done for %d shots.\n', length(on_idx));
        end
    
        function obj = normalize(obj, config)
            if isempty(obj.shot), return; end
            all_data = [obj.shot.AI_int];
            obj.normfactor = obj.area_under_curve(all_data, obj.q, config.norm_range);

        end
    end

    methods (Static, Access = private)
        function area_double = area_under_curve(data, q, roi_range)
                logical_roi = q>=roi_range(1) & q<=roi_range(2);
                q_roi = q(logical_roi);
                dq = diff(q_roi);
                W = zeros(length(q_roi),1);
                
                if ~isempty(dq)
                    W(1) = dq(1)/2;
                    W(2:end-1) = (dq(1:end-1) + dq(2:end)) / 2;
                    W(end) = dq(end) / 2;
                else
                    W(1) = 1;
                end
                
                data_roi = data(logical_roi,:);
                area_double = data_roi'*W;

            end
    end
    
    
    
    
    
end

    

