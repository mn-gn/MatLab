clc; clear; close all;

load("MnPEPC3_v2.mat");
[U, S, V] = svd(PEPC3);

y0 = [1; 0; 0];

tau = [317e-15, 7.2e-12, 1.28e-9];
time_k = 1./tau;

% Case 1
% f = 0.0996;
% ode_system = @(t,y) [
%     -f*time_k(1)*y(1) - (1-f)*time_k(2)*y(1);
%     f*time_k(1)*y(1) + (1-f)*time_k(2)*y(1) - time_k(3)*y(2);
%     time_k(3)*y(2)];

% Case 2
% f = 0.5418;
% ode_system = @(t,y) [
%     -time_k(1)*y(1);
%     time_k(1)*y(1) - f*time_k(2)*y(2) - (1-f)*time_k(3)*y(2);
%     f*time_k(2)*y(2) + (1-f)*time_k(3)*y(2)];

% Case 3
f = 0.0995;
ode_system = @(t,y) [
    -f*time_k(1)*y(1) - (1-f)*time_k(2)*y(1);
    f*time_k(1)*y(1) - time_k(3)*y(2);
    (1-f)*time_k(2)*y(1) + f*time_k(3)*y(2)];


time = logspace(-13,log10(3e-9),1000);
[~, C_matrix] = ode45(ode_system, time, y0);


figure;
plot(time,C_matrix,'LineWidth',3);
set(gca,'xscale','log');
legend('I1','I2','I3')
ylim([0 1]);
grid on

C = zeros(22,3);
t0 = 133e-15;
new_time = t - t0;
[~, C_temp] = ode45(ode_system, new_time(5:end), y0);

C(5:end,:) = C_temp;
C = C(:,1:3);
%% Get U' vector
P = C\V;
U_prime = U*S*P';

%% fitting
lb= 0;
ub = 1;
x0 = 0.5;

chi2 = @(f) (fit(PEPC3,t,f));

opt = optimoptions('lsqnonlin', ...
                   'Algorithm', 'trust-region-reflective', ...
                   'MaxIterations', 300, ...
                   'MaxFunctionEvaluations', 5e4, ...
                   'Display', 'off');

% --- Wrap into an optimization problem ---
problem = createOptimProblem('lsqnonlin', ...
                             'objective', chi2, ...
                             'options', opt, ...
                             'lb', lb, ...
                             'ub', ub, ...
                             'x0', x0);
ms = MultiStart('Display', 'iter', ...
                'UseParallel', true);

% --- Run MultiStart ---
[best_par, best_chi2, exitflag, output, solutions] = run(ms, problem, 100);
%%

function chi2 = fit(A,t,f)
    [U, S, V] = svd(A);
    y0 = [1; 0; 0];

    
    % tau = [317e-15, 7.2e-12, 1.28e-9];
    tau = [317e-15, 1.28e-9, 7.2e-12];
    time_k = 1./tau;
    % case1
    % ode_system = @(t,y) [
    % -f*time_k(1)*y(1) - (1-f)*time_k(2)*y(1);
    % f*time_k(1)*y(1) + (1-f)*time_k(2)*y(1) - time_k(3)*y(2);
    % time_k(3)*y(2)];
    

    % case2
    % ode_system = @(t,y) [
    % -time_k(1)*y(1);
    % time_k(1)*y(1) - f*time_k(2)*y(2) - (1-f)*time_k(3)*y(2);
    % f*time_k(2)*y(2) + (1-f)*time_k(3)*y(2)];

    % case3
    ode_system = @(t,y) [
    -f*time_k(1)*y(1) - (1-f)*time_k(2)*y(1);
    f*time_k(1)*y(1) - time_k(3)*y(2);
    (1-f)*time_k(2)*y(1) + f*time_k(3)*y(2)];

    
    % yeah
    C = zeros(22,3);
    t0 = 133e-15;
    new_time = t - t0;
    [~, C_temp] = ode15s(ode_system, new_time(5:end), y0);
    C(5:end,:) = C_temp;
    C = C(:,1:3);
    
    P = C\V;
    A_prime = U*S*P'*C';
    chi2 = norm(A-A_prime);
end