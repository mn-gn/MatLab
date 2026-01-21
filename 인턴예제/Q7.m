clc; clear; close all;

start_parallel(20);
%% 

M = readmatrix('Example#7_I3_for_reaction_mechanism.dat');
T = M(1, 2:end)';
Q = M(2:end, 1);
DS = M(2:end, 2:end);

[U, S, V] = svd(DS);

rank = 3;

LSV = U(:,1:rank);
RSV = V(:,1:rank);

SLSV = LSV;
SRSV = RSV;

for i = 1 : rank
    SLSV(:,i) = SLSV(:,i) * S(i,i);
    SRSV(:,i) = SRSV(:,i) * S(i,i);
end
ExpDec2 = @(t_local, tau) (t_local(1) + t_local(2)*exp(-T*tau(1)) + t_local(3)*exp(-T*tau(2)));
model = @(t) [ExpDec2(t(1:3),t(10:11)),...
              ExpDec2(t(4:6),t(10:11)),...
              ExpDec2(t(7:9),t(10:11))] ;
obj = @(t) sum((SRSV - model(t)).^2,"all");

option = optimoptions('fmincon',...
                      'Algorithm','sqp');

problem = createOptimProblem('fmincon', ...
                             'objective', obj, ...
                             'options', option, ...
                             'x0', [1 1 1 1 1 1 1 1 1 1 1]);

ms = MultiStart('Display', 'iter', 'UseParallel',true);
[t_final, fval, exitflag, output, solutions] = run(ms, problem, 200);

time_k = sort(t_final(10:11),'descend');

fit_RSV = model(t_final);
%% 

figure();
ax1 = subplot(3, 1, 1);
semilogx(T, SRSV(:,1), 'ro', 'LineWidth', 0.2); hold on;
semilogx(T,fit_RSV(:,1),'b','Linewidth',1.2);
% grid on;


ax2 = subplot(3, 1, 2);
semilogx(T, SRSV(:,2), 'ro', 'LineWidth', 0.2); hold on;
semilogx(T,fit_RSV(:,2),'b','Linewidth',1.2);
% grid on;


ax3 = subplot(3, 1, 3);
semilogx(T, SRSV(:,3), 'ro', 'LineWidth', 0.2); hold on;
semilogx(T,fit_RSV(:,3),'b','Linewidth',1.2);
% grid on;

linkaxes([ax1, ax2, ax3], 'x');

%% 

y0 = [1; 0; 0];

ode_system = @(t,y) [
    -time_k(1)*y(1);
    time_k(1)*y(1)-time_k(2)*y(2);
    time_k(2)*y(2)];

[t, C_matrix] = ode45(ode_system, T, y0);

figure;
semilogx(T,C_matrix);
title('Concentration');
xlabel('time');
ylabel('concentration');

U_prime = DS/C_matrix';
DS_prime = U_prime * C_matrix';

Discrepancy = sum((DS - DS_prime).^2,"all");


%%
Q_large = Q(1:end);
QU_prime_large = U_prime(1:end,:) .* [Q_large,Q_large,Q_large];

aff_coeff = [20.1472 4.347 18.9949 0.3814 7.5138 27.766 2.2735 66.8776 4.0712];
f1 = aff(aff_coeff,Q_large);
f2 = f1.^2;



model2 = @(r) [Delta_S3([r(1:3) r(4:6) r(13)],Q_large,f2).*Q_large, ...
               Delta_S3([r(1:3) r(7:9) r(13)],Q_large,f2).*Q_large, ...
               Delta_S3([r(1:3) r(10:12) r(13)],Q_large,f2).*Q_large];

obj2 = @(r) sum((QU_prime_large - model2(r)).^2,"all");

problem2 = createOptimProblem('fmincon', ...
                             'objective', obj2, ...
                             'options', option, ...
                             'x0', [3.1 3.3 5.6 ...
                                    3.0 3.3 4.0 ...
                                    3.1 3.4 6.5 ...
                                    3.1 1000 1000 ...
                                    1], ...
                                    ...
                             'lb', [2 2 2 ...
                                    2 2 2 ...
                                    2 2 2 ...
                                    2 1000 1000 ...
                                    0], ...
                                    ...
                             'ub', [10 10 10 ...
                                    10 10 10 ...
                                    10 10 10 ...
                                    10 1000 1000 ...
                                    1]);


[r_final, fval2] = run(ms, problem2, 1000);
figure;
fit = model2(r_final);
plot(Q_large,fit); hold on;
plot(Q_large,QU_prime_large,'o');

figure();
ax1 = subplot(3, 1, 1);
plot(Q_large,QU_prime_large(:,1), 'ro', 'LineWidth', 0.2); hold on;
plot(Q_large,fit(:,1),'b','Linewidth',1.2);
% grid on;


ax2 = subplot(3, 1, 2);
plot(Q_large,QU_prime_large(:,2), 'ro', 'LineWidth', 0.2); hold on;
plot(Q_large,fit(:,2),'b','Linewidth',1.2);
% grid on;


ax3 = subplot(3, 1, 3);
plot(Q_large,QU_prime_large(:,3), 'ro', 'LineWidth', 0.2); hold on;
plot(Q_large,fit(:,3),'b','Linewidth',1.2);
% grid on;

linkaxes([ax1, ax2, ax3], 'x');