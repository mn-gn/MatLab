clc; clear; close all;
load("Mn_U_prime.mat");
for i =1 : 3
qIq = U_prime(:,i);

p = polyfit(q, qIq, 2); 
background_baseline = polyval(p, q);

% 순수하게 진동하는 분자 신호만 남깁니다! (이게 가장 중요합니다)
qIq_oscillations = qIq - background_baseline; 

% 윈도우 함수 적용 (노이즈 억제)
qmax = max(q);
window = sin(pi * q / qmax) ./ (pi * q / qmax);
qIq_ready = qIq_oscillations .* window;

% =========================================================
% 3. 푸리에 사인 변환 (직접 수치 적분)
% =========================================================
% 찾고자 하는 r의 범위를 촘촘하게 설정합니다. (예: 0.1 ~ 10 Angstrom)
r = linspace(0.1, 10, 1000); 
F_r = zeros(size(r)); % 결과를 저장할 빈 배열

% 각 r 값에 대해 적분 수행
for i = 1:length(r)
    % 피적분 함수: q * I(q) * sin(q * r)
    integrand = qIq_ready .* sin(q .* r(i));
    
    % trapz 함수를 이용하여 q에 대해 수치 적분
    F_r(i) = trapz(q, integrand);
end

% =========================================================
% 4. 결과 시각화
% =========================================================
% figure('Position', [100, 100, 800, 600]);



% 변환된 r-space 그래프
F_r_magnitude = F_r;
plot(r, F_r_magnitude, 'LineWidth', 1.5);
hold on
end
title('Fourier Sine Transform: F(r) (r-space)');
xlabel('r (A)');
ylabel('F(r)');
grid on;

% 시각적 보조선 (용매 신호와 분자 신호가 분리됨을 확인)
xline(1.0, '--k', 'Low-r (용매/배경 영역)', 'LabelOrientation', 'horizontal');

