clc; close all;

c = parcluster; 
fprintf('현재 프로파일 이름: %s\n', c.Profile);
fprintf('최대 허용 워커 수: %d\n', c.NumWorkers);
start_parallel(20);

data = readtable("Example#4_I3_deltaS.dat");
q_data = data.Var1;
S_data = data.Var2;

w = 1; %q_data; % 가중치

model = @(r, q) r(7).*(Debye_I(q,[r(4) r(5) r(6)]) - Debye_I(q, [r(1) r(2) r(3)]));
chi2_func = @(r) sum((w.*(S_data - model(r, q_data))).^2);

% r = [r_grd (1 2 theta) ...
%      r_exd (1 2 theta)...
%      scaling_factor]

r0 = [2.0 2.0 2.0 ...
      2.0 2.0 2.0 ...
      2];  % initial condition

lb = [1.0 1.0 0.0 ...
      1.0 1.0 0.0 ...
      0];  % lower bound

ub = [6.0 6.0 pi ...
      6.0 6.0 pi ...
      6.0];  % upper bound

option = optimoptions('fmincon', ...
                      'Algorithm','sqp');

problem = createOptimProblem('fmincon', ...
                             'objective', chi2_func, ...
                             'options', option, ...
                             'x0', r0, ...
                             'lb', lb, ...
                             'ub', ub);

% gs = GlobalSearch('Display','iter','NumTrialPoints',1000);
% [r_final, fval] = run(gs,problem);
ms = MultiStart('Display', 'iter', 'UseParallel',true);
[r_final, fval, exitflag, output, solutions] = run(ms, problem, 100);

r_final_sorted = r_final;

for i = 1 : 2
    r_final_sorted(3*i) = sqrt(r_final_sorted(3*i-1)^2 + r_final_sorted(3*i-2)^2 - 2 * r_final_sorted(3*i-1) * r_final_sorted(3*i-2) * cos(r_final_sorted(3*i)));
    % cond_sorted(3*i) = sqrt(cond_sorted(3*i-1)^2 + cond_sorted(3*i-2)^2 - 2*cond_sorted(3*i-1) * cond_sorted(3*i-2) * cos(cond_sorted(3*i)));
end

fprintf('  결과 \n g  : %.4f %.4f %.4f \n e1 : %.4f %.4f %.4f \n e2 : %.4f %.4f %.4f \n scaling factor %.4f \n',r_final_sorted);
% fprintf('  정답 \n g  : %.4f %.4f %.4f \n e1 : %.4f %.4f %.4f \n e2 : %.4f %.4f %.4f \n scaling factor %.4f \n',cond_sorted);
S_fit = model(r_final, q_data);
%% figure
figure();
hold on;
plot(q_data, S_data,'ko','DisplayName','Experimental Data');
plot(q_data, S_fit,'r-','lineWidth',2,'DisplayName','Best Fit');
xlabel('q'); ylabel('\Delta S');
legend;
grid on;
title(['Best Fit Results (min \chi^2: ', num2str(fval), ')']);

% A = zeros(50,1);
% 
% for i = 1 : 50
%     A(i) = solutions(i).Fval;
% end