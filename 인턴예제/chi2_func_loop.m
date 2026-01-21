function [f, S_fit_local] = chi2_func_loop(par,~)
    % 중요: 함수 안에서 데이터를 직접 읽거나 persistent 변수로 관리
    persistent q_data S_data
    if isempty(q_data)
        temp = readtable("Example#4_I3_deltaS.dat");
        q_data = temp{:, 1};
        S_data = temp{:, 2};
    end

    % 모델 계산
    S_fit_local = par(7) .* (Debye_I(q_data, [par(4) par(5) par(6)]) - Debye_I(q_data, [par(1) par(2) par(3)]));
    
    % Chi-square (벡터 연산)
    f = sum((S_data - S_fit_local).^2);
end