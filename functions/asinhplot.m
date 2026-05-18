function h = asinhplot(x, y, varargin)
    % 기본 스케일
    S = 1e-15; 
    plot_args = {}; % plot에 전달할 인수들

    % 사용자가 입력한 추가 옵션(varargin) 분석
    if ~isempty(varargin)
        % 첫 번째 추가 인수가 숫자(스칼라)면 S로 취급
        if isnumeric(varargin{1}) && isscalar(varargin{1})
            S = varargin{1};
            plot_args = varargin(2:end); % 나머지는 플롯 옵션
        else
            % 숫자가 아니면(예: 'r-') 전부 플롯 옵션으로 취급
            plot_args = varargin;
        end
    end

    % 1. 데이터를 asinh 공간으로 변환하여 플롯 (옵션 전달)
    x_trans = asinh(x / S);
    h_temp = plot(x_trans, y, plot_args{:}); 
    grid on;

    % 2. 눈금(Ticks) 자동 설정
    x_min = min(x); 
    x_max = max(x);
    
    % 데이터 최대 절댓값 확인 (0일 경우를 대비해 최솟값을 S로 보정)
    max_val = max(abs([x_min, x_max]));
    if max_val == 0
        max_val = S;
    end
    
    % 10의 지수승 베이스 생성
    powers = floor(log10(S)):1:(ceil(log10(max_val)) + 1);
    base_ticks = [0, 10.^powers, -10.^powers];
    real_ticks = sort(unique(base_ticks));
    
    % --- 수정된 부분: 데이터 범위를 바깥에서 감싸는 10의 지수승 경계 찾기 ---
    % x_min보다 작거나 같은 tick 중 가장 큰 값
    tick_min = max(real_ticks(real_ticks <= x_min));
    % x_max보다 크거나 같은 tick 중 가장 작은 값
    tick_max = min(real_ticks(real_ticks >= x_max));
    
    % 예외 처리 (데이터가 없을 경우 등)
    if isempty(tick_min), tick_min = x_min; end
    if isempty(tick_max), tick_max = x_max; end
    
    % 예외 처리 (데이터가 정확히 특정 10의 지수승 값 하나만 있을 경우 xlim 에러 방지)
    if tick_min >= tick_max
        idx = find(real_ticks == tick_min, 1);
        if ~isempty(idx) && idx > 1
            tick_min = real_ticks(idx - 1);
        else
            tick_min = tick_min - S;
        end
        if ~isempty(idx) && idx < length(real_ticks)
            tick_max = real_ticks(idx + 1);
        else
            tick_max = tick_max + S;
        end
    end
    
    % 계산된 경계(tick_min, tick_max)에 해당하는 범위 내의 10의 지수승만 남김
    real_ticks = real_ticks(real_ticks >= tick_min & real_ticks <= tick_max);
    
    % x축 tick 적용 및 xlim을 외부 tick에 맞춰 확장
    xticks(asinh(real_ticks / S));
    xlim([asinh(tick_min / S), asinh(tick_max / S)]); % 축도 바깥 tick까지 뻗도록 설정
    % -------------------------------------------------------------------
    
    % 라벨 포맷팅
    labels = cell(size(real_ticks));
    for i = 1:length(real_ticks)
        if real_ticks(i) == 0
            labels{i} = '0';
        else
            labels{i} = sprintf('%.0e', real_ticks(i));
        end
    end
    xticklabels(labels);
    xlabel(sprintf('Linear within \\pm%g / Log scale outside', S));

    % 출력 인수가 요청되었을 때만 h 반환
    if nargout > 0
        h = h_temp;
    end
end