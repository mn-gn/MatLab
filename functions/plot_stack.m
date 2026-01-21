function plot_stack(q, t, A, offset)

%
% Stacked line plot of A(q, t) with vertical offset.
%
% INPUT
% q      : vector of q-values (row or column)
% t      : vector of t-labels (row or column)
% A      : (num_q x num_t) matrix
% offset : scalar vertical offset between traces (positive = upward, negative = downward)

q = q(:);
t = t(:);

if length(q) ~= size(A, 1)
    error('Length of q (%d) must match number of rows in A (%d).', length(q), size(A, 1));
end

if length(t) ~= size(A, 2)
    error('Length of t (%d) must match number of columns in A (%d).', length(t), size(A, 2));
end

[num_q, num_t] = size(A);

offsets = (0:num_t - 1) * offset;

hold on
for i = 1:num_t
    plot(q, zeros(num_q, 1) + offsets(i), 'k:', 'LineWidth', 0.5);
    plot(q, A(:, i) + offsets(i), 'k-', 'LineWidth', 0.7);
end

if offset > 0
    yticks(offsets);
    yticklabels(arrayfun(@(x) sprintf('%.2e', x), t, 'UniformOutput', false));
else
    yticks(fliplr(offsets));
    yticklabels(flip(arrayfun(@(x) sprintf('%.2e', x), t, 'UniformOutput', false)));
end

axis tight
set(gca, 'TickLength', [0, 0]);
xlabel('q (Å^{-1})')
ylabel('t (s)')
title('\DeltaS(q, t)')
hold off
end
