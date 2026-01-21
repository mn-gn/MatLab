clc; close all;




data = readtable("Example#4_I3_deltaS.dat");
q_data = data{:, 1}; 
S_data = data{:, 2};

num_of_par = 7;
par_num = (1:num_of_par)';

step_size = 0.01 * ones(num_of_par,1);
lb = ones(num_of_par,1);
ub = 6 * ones(num_of_par,1);
% initial_guess = 2 * ones(num_of_par,1);

for i = 1 : 2
    lb(3*i) = 0;
    ub(3*i) = pi;
end
lb(num_of_par) = 0; 


run_num = 50;

% best_chi2 = zeros(size(run_num));
% best_pars = zeros(num_of_par,1);
% best_S = zeros(size(S_data));

all_pars = zeros(7,50);
all_chi2 = zeros(50,1);
all_S_fit = zeros(80,50);

minuit_cmd = 'call fcn ; seek 100 ; minimize ; improve ; call fcn';
for i = 1 : run_num
    initial_guess = rand(num_of_par,1).*(ub-lb) + lb;
    lims = [par_num, step_size, lb, ub];

    data_mat = zeros(2); %% dummy matrix for function input (ignore)
    obj_handle = @(p) chi2_func_loop(p);

    [pars, errs, chi2] = fminuit('chi2_func_loop',initial_guess,data_mat, '-c', minuit_cmd, '-s',lims); %% Executing fitting loop.
    [~, S_fit_temp] = chi2_func_loop(pars);

    all_pars(:, i) = pars;
    all_chi2(i) = chi2;
    all_S_fit(:, i) = S_fit_temp;
end

[min_chi2, best_idx] = min(all_chi2);
best_S_fit = all_S_fit(:, best_idx);
best_par = all_pars(:,best_idx);

figure();
plot(q_data, S_data, 'ko', 'DisplayName', 'Experimental Data'); hold on;
plot(q_data, best_S_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Best Fit');
xlabel('q'); ylabel('\Delta S');
legend;
grid on;
title(['Best Fit Results (min \chi^2: ', num2str(min_chi2), ')']);
