function [mask, header] = read_mask(filename)
    % READ_MASK : .edf 형식의 마스크 파일을 읽어옵니다.
    % 사용법: mask = read_mask('Mask_run4.edf');
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('파일을 열 수 없습니다: %s', filename);
    end

    % 1. 헤더 읽기 ( "{" 부터 "}" 까지)
    header_text = '';
    header_end = false;
    
    while ~feof(fid)
        line = fgets(fid);
        header_text = [header_text, line];
        
        % 헤더 끝 표시('}')를 찾으면 멈춤
        if contains(line, '}')
            header_end = true;
            break; 
        end
    end
    
    if ~header_end
        fclose(fid);
        error('EDF 헤더의 끝(})을 찾을 수 없습니다.');
    end
    
    % 2. 헤더에서 정보 파싱 (Width, Height, Type)
    dim_1 = regexp(header_text, 'Dim_1\s*=\s*(\d+)', 'tokens');
    dim_2 = regexp(header_text, 'Dim_2\s*=\s*(\d+)', 'tokens');
    data_type_str = regexp(header_text, 'DataType\s*=\s*(\w+)', 'tokens');
    byte_order_str = regexp(header_text, 'ByteOrder\s*=\s*(\w+)', 'tokens');
    
    if isempty(dim_1) || isempty(dim_2)
        fclose(fid);
        error('이미지 크기 정보(Dim_1, Dim_2)를 찾을 수 없습니다.');
    end
    
    width = str2double(dim_1{1}{1});  % Dim_1 = Col (가로)
    height = str2double(dim_2{1}{1}); % Dim_2 = Row (세로)
    
    % 3. 데이터 타입 결정
    precision = 'uint16'; % 기본값
    if ~isempty(data_type_str)
        type_val = data_type_str{1}{1};
        switch type_val
            case 'UnsignedShort'
                precision = 'uint16';
            case 'SignedInteger'
                precision = 'int32';
            case 'UnsignedInteger'
                precision = 'uint32';
            case 'FloatValue'
                precision = 'single';
            case 'DoubleValue'
                precision = 'double';
            case 'UnsignedByte'
                precision = 'uint8';
        end
    end
    
    % 4. 바이트 순서 (Endian)
    machinefmt = 'ieee-le'; % 기본 리틀 엔디안
    if ~isempty(byte_order_str)
        if strcmpi(byte_order_str{1}{1}, 'HighByteFirst')
            machinefmt = 'ieee-be';
        end
    end

    % 5. 바이너리 데이터 읽기
    % (fgets가 끝난 지점부터 바로 데이터가 시작됨)
    % 주의: EDF 헤더 뒤에 줄바꿈(\n)이나 공백이 있을 수 있어 위치 조정이 필요할 수 있음
    
    % 현재 위치 저장
    pos = ftell(fid);
    
    % 데이터 읽기
    [data, count] = fread(fid, width * height, precision, machinefmt);
    
    % 만약 데이터가 덜 읽혔으면, 헤더 패딩(보통 1024바이트) 문제일 수 있음
    if count ~= width * height
        % 헤더가 보통 1024, 512 바이트 단위로 끊김. 강제로 이동해봄.
        fseek(fid, 1024, 'bof'); 
        [data, count] = fread(fid, width * height, precision, machinefmt);
        
        if count ~= width * height
             fclose(fid);
             error('데이터 크기가 맞지 않습니다. 헤더 파싱 오류 가능성.');
        end
    end
    
    fclose(fid);
    
    % 6. 형상 변환 (1D -> 2D)
    % EDF는 Row-major(C스타일)로 저장되므로, MATLAB(Column-major)으로 읽은 뒤 Transpose 필요
    mask = reshape(data, [width, height])'; 
end