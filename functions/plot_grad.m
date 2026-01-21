function plot_grad(q, t, A, color_map)

%
% Plots multiple curves of A(q, t) with a gradient colormap, labeling by t.
%
% INPUT
% q         : vector (row OR column)
% t         : vector (row OR column)
% A         : (num_q x num_t) matrix
% color_map : (optional) 'jet', 'turbo', 'cool' (default = 'jet')
%
% NOTES
% • Any unsupported colormap name will fall back to 'jet'.

q = q(:);
t = t(:);

% ---- Default colormap ----
if nargin < 4 || isempty(color_map)
    color_map = 'jet';
end

% ---- Restrict to valid options ----
valid_cmaps = {'jet', 'turbo', 'cool'};
if ~ismember(lower(color_map), valid_cmaps)
    warning('Unsupported colormap "%s". Using ''jet'' instead.', color_map);
    color_map = 'jet';
end

% ---- Checks ----
if length(q) ~= size(A, 1)
    error('Length of q (%d) must match number of rows in A (%d).', length(q), size(A, 1));
end

if length(t) ~= size(A, 2)
    error('Length of t (%d) must match number of columns in A (%d).', length(t), size(A, 2));
end

[~, num_t] = size(A);

% ---- Generate colors ----
switch lower(color_map)
    case 'jet'
        colors = jet(num_t);
    case 'turbo'
        colors = turbo(num_t);
    case 'cool'
        colors = cool(num_t);
end

% ---- Plot curves ----
hold on
for i = 1:num_t
    plot(q, A(:, i), '-', 'Color', colors(i, :), 'LineWidth', 0.7);
end
hold off

% ---- Colormap & colorbar ----
colormap(colors);
cb = colorbar;
ylabel(cb, 't (s)')
cb.TickLength = 0;

% ---- Colorbar tick sampling ----
max_labels = 20;
if num_t == 1
    cb.Ticks      = 0.5;
    cb.TickLabels = {sprintf('%.1e', t)};
else
    if num_t <= max_labels
        idx = 1:num_t;
        pos = linspace(0, 1, length(idx));
    else
        idx = round(linspace(1, num_t, max_labels));
        pos = linspace(0, 1, length(idx));
    end
    cb.Ticks      = pos;
    cb.TickLabels = arrayfun(@(x) sprintf('%.2e', x), t(idx), 'UniformOutput', false);
end

xlim([min(q), max(q)]);
xlabel('q (Å^{-1})')
ylabel('\DeltaS(q) (a.u.)')
title('\DeltaS(q, t)')
end
