function start_parallel(num_worker, pool_type)

%
% Open and manage a parallel pool (threads or local).
%
% INPUT
%   num_worker :
%       0 (default)
%           - Auto mode.
%           - If a pool already exists with the same pool_type, keep it as-is.
%           - If no pool exists, start a pool using the maximum
%             number of workers allowed by the selected profile.
%       Positive integer N
%           - Explicit mode.
%           - If a pool exists with a different number of workers,
%             restart the pool using N workers (auto-reduced if needed).
%   pool_type :
%       "threads" (default)
%       "local"

arguments
    num_worker (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    pool_type  (1,1) string {mustBeMember(pool_type, ["threads","local"])} = "local"
end

p = gcp('nocreate');
auto_mode = (num_worker == 0);

% If a pool exists, decide whether to keep or restart
if ~isempty(p)

    % Detect current pool type without touching p.Cluster (threads-safe)
    if isa(p, 'parallel.ThreadPool')
        current_type = "threads";
    elseif isa(p, 'parallel.ProcessPool')
        current_type = "local";
    else
        current_type = "unknown";
    end

    % If pool type differs, always restart
    if current_type ~= pool_type
        delete(p);
        p = [];
    else
        % Same pool type
        % - Auto mode: keep as-is
        % - Explicit : restart only if different NumWorkers
        if auto_mode || p.NumWorkers == num_worker
            fprintf('Parallel pool already running (%d workers, %s).\n', ...
                    p.NumWorkers, pool_type);
            return
        else
            delete(p);
            p = [];
        end
    end
end

% Decide requested workers in auto mode
if auto_mode
    switch pool_type
        case "threads"
            num_worker = feature('numcores');
        case "local"
            c = parcluster('local');
            num_worker = c.NumWorkers;
    end
end

% Start pool (reduce if capped)
try
    parpool(pool_type, num_worker);
catch ME
    if contains(ME.message, 'Too many workers requested', 'IgnoreCase', true)
        tok = regexp(ME.message, 'maximum of\s+(\d+)\s+workers', 'tokens', 'once');
        if isempty(tok)
            rethrow(ME);
        end
        num_worker_requested = num_worker;
        num_worker = str2double(tok{1});
        fprintf(['Requested %d workers exceeds the maximum allowed.\n' ...
                 'Using %d workers instead (%s pool).\n'], ...
                 num_worker_requested, num_worker, pool_type);
        parpool(pool_type, num_worker);
    else
        rethrow(ME);
    end
end

fprintf('Parallel pool started (%d workers, %s).\n', num_worker, pool_type);
end
