function y = math_sinc(x)
    % MATH_SINC: 숫자, 행렬, 심볼릭 모두 지원하는 sin(x)/x
    
    if isa(x, 'sym')
        % 1. 심볼릭 입력인 경우: 수식으로 반환하되 0에서 limit 처리 가능하게 함
        % 심볼릭에서는 x=0 대입 시 NaN이 뜨면 limit(y, x, 0)으로 해결하거나
        % 아래와 같이 piecewise를 사용하여 정의할 수 있습니다.
        y = piecewise(x == 0, 1, x ~= 0, sin(x)/x);
    else
        % 2. 수치형(숫자, 행렬) 입력인 경우
        y = ones(size(x));
        idx = (x ~= 0);
        y(idx) = sin(x(idx)) ./ x(idx);
    end
end