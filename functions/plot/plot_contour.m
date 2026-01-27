function plot_contour(q, t, A, color_scale, color_map)

%
% Pseudocolor contour (filled surface) of A(q, t) with flexible colorbar options.
%
% INPUT
% q           : vector (row or column)
% t           : vector (row or column)
% A           : (num_q x num_t) matrix
% color_scale : abs(c_min and c_max) OR [c_min, c_max], default = [min(A(:)), max(A(:))] (optional)
% color_map   : string (e.g., 'jet', 'turbo', 'gray', 'parula'), default = 'jet' (optional)

q = q(:);
t = t(:);

if length(q) ~= size(A, 1)
    error('Length of q (%d) must match number of rows in A (%d).', length(q), size(A, 1));
end

if length(t) ~= size(A, 2)
    error('Length of t (%d) must match number of columns in A (%d).', length(t), size(A, 2));
end

[num_q, num_t] = size(A);

% ---- Defaults ----
if nargin < 4 || isempty(color_scale)
    cmin = min(A(:));
    cmax = max(A(:));
elseif isscalar(color_scale)
    cmin = -abs(color_scale);
    cmax = abs(color_scale);
else
    cmin = color_scale(1);
    cmax = color_scale(2);
end

if nargin < 5 || isempty(color_map)
    color_map = 'jet';
end

% ---- Plot ----
surf(q, t, A.', 'EdgeColor', 'none', 'FaceColor', 'interp');
set(gca, 'YDir', 'reverse');
xlim([min(q), max(q)]);
ylim([min(t), max(t)]);
clim([cmin, cmax]);
axis tight
axis square
xlabel('q (Å^{-1})')
ylabel('t (s)')
title('\DeltaS(q, t)')
view(2); % top-down view

% ---- Colormap ----
switch lower(color_map)
    case 'jet'
        colormap(jet);
    case 'turbo'
        colormap(turbo);
    case 'gray'
        colormap(flipud(gray));
    case 'parula'
        colormap(parula);
    otherwise
        warning('Unknown color_map "%s". Using default "jet".', color_map);
        colormap(jet);
end

% ---- Colorbar ----
cb = colorbar;
ylabel(cb, '\DeltaS(q, t) (a.u.)')
end
