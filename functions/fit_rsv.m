function [best_par_struct, V_fit, save_struct] = fit_rsv(t, V, num_tau, num_start)

%
% Fit RSVs with a Gaussian-IRF convolved multi-exponential decay model.
%
% INPUT
%    t         : time vector (row/column). Internally converted to a column vector.
%    V         : (num_t x num_rsv) matrix of experimental RSVs (each column is one RSV).
%    num_tau   : number of exponential components (default = 1)
%    num_start : number of MultiStart start points (default = 1)
%
% OUTPUT
%    best_par_struct       : struct containing best-fit parameters
%                   .A     : (num_rsv x num_tau) amplitudes
%                   .tau   : (num_tau x 1) time constants (shared across RSVs)
%                   .A0    : (num_rsv x 1) offsets
%                   .FWHM  : scalar IRF FWHM (shared)
%                   .t0    : scalar time-zero shift (shared)
%                   .chi2  : best chi2 returned by MultiStart
%    V_fit                 : (num_t x num_rsv) fitted curves evaluated on the original time grid
%    save_struct           : struct containing optimization outputs
%               .exitflag  : exitflag from MultiStart
%               .output    : output struct from MultiStart
%               .solutions : Solution objects from MultiStart
%
% NOTE
%   This function uses parallel evaluation for MultiStart if a parallel pool
%   is already running. If no pool exists, it runs in serial mode.
%   (e.g., start_parallel(...) can be called before fit_rsv to enable parallel.)

arguments
    t double
    V double
    num_tau   (1,1) double {mustBeInteger, mustBePositive} = 1
    num_start (1,1) double {mustBeInteger, mustBePositive} = 1
end

% Ensure t is a column vector
t = t(:);

% Dimensions
[~, num_rsv] = size(V);

% --- Parameter vector layout ---
% par = [A(:); tau(:); A0(:); FWHM; t0]
i1      = num_rsv * num_tau; % end index for A(:)
i2      = i1 + num_tau;      % end index for tau(:)
i3      = i2 + num_rsv;      % end index for A0(:)
num_par = i3 + 2;            % + [FWHM; t0]

% --- Bound constraints ---
lb = zeros(num_par, 1);
ub = zeros(num_par, 1);
typ = ones(num_par, 1);

% Amplitudes A
lb(1:i1)  = -1;
ub(1:i1)  = 1;
typ(1:i1) = 1;

% Time constants tau (log scale)
lb(i1+1:i2)  = log10(1e-13);
ub(i1+1:i2)  = log10(3.16e-1);
typ(i1+1:i2) = 1;

% Offsets A0
lb(i2+1:i3)  = -1;
ub(i2+1:i3)  = 1;
typ(i2+1:i3) = 1;

% IRF FWHM (log scale)
lb(i3+1)  = log10(50e-15);
ub(i3+1)  = log10(10e-12);
typ(i3+1) = 1;

% IRF t0
lb(i3+2)  = -1e-12;
ub(i3+2)  = 1e-12;
typ(i3+2) = 1e-15;

% --- Initial guess ---
x0 = (ub + lb) / 2;

% --- Residual function for lsqnonlin ---
chi = @(par) (V - model_rsv(par, t));

% --- Solver options ---
opt = optimoptions('lsqnonlin', ...
                   'Algorithm', 'trust-region-reflective', ...
                   'TypicalX', typ, ...
                   'MaxIterations', 300, ...
                   'MaxFunctionEvaluations', 5e4, ...
                   'Display', 'off');

% --- Wrap into an optimization problem ---
problem = createOptimProblem('lsqnonlin', ...
                             'objective', chi, ...
                             'options', opt, ...
                             'lb', lb, ...
                             'ub', ub, ...
                             'x0', x0);

% --- Enable parallel MultiStart only if a parallel pool is already running ---
if ~isempty(gcp('nocreate'))
    parallel_opt = true;
else
    parallel_opt = false;
end

% --- MultiStart configuration ---
ms = MultiStart('Display', 'iter', ...
                'UseParallel', parallel_opt);

% --- Run MultiStart ---
[best_par, best_chi2, exitflag, output, solutions] = run(ms, problem, num_start);

% --- Unpack best parameter vector ---
best_A    = reshape(best_par(1:i1), [num_rsv, num_tau]);
best_tau  = 10.^best_par(i1+1:i2);
best_A0   = best_par(i2+1:i3);
best_FWHM = 10.^best_par(i3+1);
best_t0   = best_par(i3+2);

% --- Package best-fit parameters into a struct ---
best_par_struct      = struct();
best_par_struct.A    = best_A;
best_par_struct.tau  = best_tau;
best_par_struct.A0   = best_A0;
best_par_struct.FWHM = best_FWHM;
best_par_struct.t0   = best_t0;
best_par_struct.chi2 = best_chi2;

% --- Evaluate fitted model ---
V_fit = model_rsv(best_par, t);

% --- Save optimization diagnostics ---
save_struct = struct();
save_struct.exitflag   = exitflag;
save_struct.output     = output;
save_struct.solutions  = solutions;

% --- Draw figures ---
t_dense = asinhspace(t(1), t(end), 10000);
t_dense = t_dense(:);
V_fit_dense = model_rsv(best_par, t_dense);
close all;
for i = 1:num_rsv
    figure();
    % set(gca,'xscale','log')
    hold on
    asinhplot(t, V(:, i), 'k.', ...
        "MarkerSize", 15)
    asinhplot(t_dense, V_fit_dense(:, i), 'r-', ...
        'LineWidth', 2)
    hold off
    title(sprintf('RSV%d Fitting', i))
    xlabel('Time (s)')
    ylabel(sprintf('RSV%d (a.u.)', i))
    legend('Exp', 'Fit')
end

    figure();
    % set(gca,'xscale','log')
    hold on
for i = 1:num_rsv
    asinhplot(t_dense, V_fit_dense(:, i), '-', ...
        'LineWidth', 2)
end
    hold off

    function V_model = model_rsv(par, t)

        %
        % Evaluate the global multi-exponential model for all RSVs.
        %
        % INPUT
        %   par : parameter vector [A(:); tau(:); A0(:); FWHM; t0]
        %
        % OUTPUT
        %   V_model : (num_t x num_rsv) model RSVs

        % Unpack parameter vector
        A    = reshape(par(1:i1), [num_rsv, num_tau]);
        tau  = 10.^par(i1+1:i2);
        A0   = par(i2+1:i3);
        FWHM = 10.^par(i3+1);
        t0   = par(i3+2);

        % Compute model for each RSV column
        V_model = zeros(length(t), num_rsv);
        for j = 1:num_rsv
            V_model(:, j) = conv_gauss_expdec(t, A(j, :), A0(j), tau, t0, FWHM);
        end
    end
end

function y = conv_gauss_expdec(t, A, A0, tau, t0, FWHM)

%
% Gaussian-IRF convolved sum of exponential decays.
%
% INPUT
%   t    : time vector (length = num_t). Row or column vector accepted.
%   A    : amplitude vector (length = n). Row or column vector accepted.
%   A0   : scalar offset
%   tau  : time-constant vector (length = n). Row or column vector accepted.
%   t0   : scalar time-zero shift
%   FWHM : scalar Gaussian IRF full width at half maximum
%
% OUTPUT
%   y    : column vector [num_t x 1]

% Enforce expected shapes
t   = t(:);
tau = tau(:);

num_tau = length(tau);

% Convert FWHM to sigma for a Gaussian
sigma = FWHM / 2.3548;

% Initialize with offset
y = A0 * ones(length(t), 1);

% Precompute constants for efficiency
inv_sqrt2_sigma = 1 / (sqrt(2) * sigma);
gauss_env = exp(-((t - t0).^2) ./ (2 * sigma^2));

% Closed-form expression using complementary error function
for k = 1:num_tau
    % 여오차 함수에 들어갈 인자 Z 계산
    Z = ((sigma^2)/tau(k) - (t - t0)) * inv_sqrt2_sigma;
    
    % 현재 tau에 대한 성분 배열 초기화
    comp = zeros(length(t), 1);
    
    % 조건 1: Z >= 0 (작은 tau로 인한 폭발 방지, erfcx 사용)
    idx_pos = (Z >= 0);
    if any(idx_pos)
        comp(idx_pos) = gauss_env(idx_pos) .* erfcx(Z(idx_pos));
    end
    
    % 조건 2: Z < 0 (큰 t에서 발생하는 0 * Inf = NaN 방지, 원래 수식 사용)
    idx_neg = (Z < 0);
    if any(idx_neg)
        % 원래의 지수 항 계산 (Z < 0 이면 이 값은 무조건 음수이므로 안전함)
        E = (sigma^2)/(2*tau(k)^2) - (t(idx_neg) - t0)/tau(k);
        comp(idx_neg) = exp(E) .* erfc(Z(idx_neg));
    end
    
    % 최종 결과에 현재 성분 추가
    y = y + 0.5 * A(k) .* comp;
end
end
