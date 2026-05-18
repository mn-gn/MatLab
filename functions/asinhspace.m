function y = asinhspace(x1, x2, n, S)
    % asinhspace_femto: 0 근처에서는 선형, 외부에서는 로그 간격으로 점 생성
    % x1: 시작값
    % x2: 끝값
    % n: 점의 개수
    % S: 선형-로그 전환 기준 스케일 (입력하지 않으면 기본값 1e-15)

    if nargin < 4
        S = 1e-15; % 기본 스케일: 1 펨토
    end

    % asinh 공간에서 균일하게 나눈 뒤 다시 sinh로 복원
    y = S * sinh(linspace(asinh(x1/S), asinh(x2/S), n));
end