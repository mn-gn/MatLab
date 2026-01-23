function poni_struct = read_poni(filename)
    % READ_PONI : pyFAI PONI 파일을 읽어옵니다. (JSON 호환 버전)
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('파일을 열 수 없습니다: %s', filename);
    end

    poni_struct = struct();

    while ~feof(fid)
        tline = fgetl(fid);
        
        % 1. 주석(#)이나 빈 줄 무시
        if isempty(tline) || startsWith(strtrim(tline), '#')
            continue;
        end
        
        % 2. 첫 번째 콜론(:) 위치 찾기 (이게 핵심! ⭐)
        idx = strfind(tline, ':');
        if isempty(idx)
            continue;
        end
        first_colon = idx(1); % 무조건 맨 처음 콜론만 봅니다.
        
        % 3. 키(Key)와 값(Value) 분리
        key = strtrim(tline(1:first_colon-1));
        val_str = strtrim(tline(first_colon+1:end));
        
        % 4. 숫자인지 체크하고 변환
        val_num = str2double(val_str);
        
        if ~isnan(val_num)
            % 숫자면 바로 저장
            poni_struct.(key) = val_num;
        else
            % 숫자가 아니면 (Detector_config 같은 경우)
            % 혹시 JSON 데이터인가? (중괄호 {} 체크)
            if startsWith(val_str, '{') && endsWith(val_str, '}')
                try
                    % 최신 매트랩(R2016b+)이면 JSON 해석 시도
                    json_data = jsondecode(val_str);
                    poni_struct.(key) = json_data;
                    
                    % 편의를 위해 PixelSize도 밖으로 꺼내줌 (있으면)
                    if isfield(json_data, 'pixel1')
                        poni_struct.PixelSize1 = json_data.pixel1;
                        poni_struct.PixelSize2 = json_data.pixel2;
                    end
                catch
                    % 옛날 매트랩이면 그냥 문자열로 저장
                    poni_struct.(key) = val_str;
                end
            else
                % 그냥 문자열 저장 (Detector 이름 등)
                poni_struct.(key) = val_str;
            end
        end
    end
    
    fclose(fid);
end