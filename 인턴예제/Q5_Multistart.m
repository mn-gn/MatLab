clc; clear; close all;

c = parcluster; 
fprintf('현재 프로파일 이름: %s\n', c.Profile);
fprintf('최대 허용 워커 수: %d\n', c.NumWorkers);
start_parallel(20);

T = readtable("Example#5_Au3_deltaS_time_series.dat");
T = T(2:351,:);

q_data = T.Var1;
S_data_all = table2array(T(:,2:end));
[num_q, num_frames] = size(S_data_all); % 350 x 101


aff_coeff = readtable("atomic_form_factor_coefficients.dat",'VariableNamingRule','preserve').("Au");
f1 = aff(aff_coeff,q_data);
f2 = f1.^2;

chi2_func = @(p) total_chi2(p,q_data,S_data_all,f2)*1e-4;

% global parameter
r_g_0 = [2.0; 2.0; 2.0; 0.0005];
lb_g  = [1.0; 1.0; 0.0; 0.0001];
ub_g  = [8.0; 8.0; pi;  0.0010];

% local parameter
r_l_0 = [2.0; 2.0; 2.0];
lb_l  = [1.0; 1.0; 0.0];
ub_l  = [8.0; 8.0; pi];

% all constraints
p0 = [r_g_0; repmat(r_l_0, num_frames, 1)];
lb = [lb_g;  repmat(lb_l,  num_frames, 1)];
ub = [ub_g;  repmat(ub_l,  num_frames, 1)];


% option = optimoptions('fmincon', ...
%                       'Algorithm','sqp');

option = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'HessianApproximation', 'lbfgs', ...
    'MaxFunctionEvaluations', 2e5, ...
    'StepTolerance',1e-6);

problem = createOptimProblem('fmincon', ...
                             'objective', chi2_func, ...
                             'options', option, ...
                             'x0', p0, ...
                             'lb', lb, ...
                             'ub', ub);


ms = MultiStart('Display', 'iter', 'UseParallel',true);
[p_final, fval, exitflag, output, solutions] = run(ms, problem, 1);

final_r123 = p_final(1:3);
final_r7   = p_final(4);
final_locals = reshape(p_final(5:end), 3, num_frames);

fprintf('\n=== 피팅 완료 ===\n');
fprintf('Global r1, r2, r3: %.4f, %.4f, %.4f\n', final_r123);
fprintf('Global Scale r7: %.6f\n', final_r7);

S = zeros(size(S_data_all));

for i = 1 : 101
    r = [final_r123; final_locals(:,i); final_r7];
    S(:,i) = Delta_S3(r,q_data,f2);
    
end

asdf = final_locals;
for i = 1: 101

asdf(3,i) = sqrt(final_locals(1,i)^2 + final_locals(2,i)^2 - 2*final_locals(1,i)*final_locals(2,i)*cos(final_locals(3,i)));
asdf(:,i) = sort(asdf(:,i));

end

figure;
subplot(2,1,1);
hold on;
first_frame_r = [final_r123; final_locals(:,1); final_r7];
S_fit_1 = Delta_S3(q_data, first_frame_r, f2);
plot(q_data, S_mat(:,1), 'k.', 'DisplayName', 'Exp Frame 1');
plot(q_data, S_fit_1, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Best Fit');
legend; title('First Frame Fitting');

subplot(2,1,2);
plot(1:num_frames, final_locals(1,:), 'b-o', 'DisplayName', 'r4 (Local)');
hold on;
plot(1:num_frames, final_locals(2,:), 'r-s', 'DisplayName', 'r5 (Local)');
xlabel('Frame Number'); ylabel('Parameter Value');
legend; title('Local Parameter Evolution');


function total_err = total_chi2(p,q,S_all,f2)
g_123 = p(1:3);
g_7 = p(4);

l_par = reshape(p(5:end),3,[]);

num_f = size(S_all,2);
total_err = 0;
for j = 1:num_f
    rj = [g_123; l_par(:,j); g_7];

    S_model = Delta_S3(rj,q,f2);
    err = S_all(:,j) - S_model;
    total_err = total_err + sum(err.*err);
end
end

