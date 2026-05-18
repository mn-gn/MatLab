function start_parallel(num_worker, pool_type)

%
% Open and manage a parallel pool (local or threads).
%
% INPUT
%   num_worker :
%       0 (default)
%           - Auto mode (always maximum workers for the requested pool_type).
%       Positive integer N
%           - Explicit mode (restart and try N workers; auto-reduce if capped).
%   pool_type :
%       "local" (default)
%           - Process-based parallel pool on the local machine.
%       "threads"
%           - Thread-based parallel pool (shared memory, parfor only).
%
% NOTE
%   - Always restart the pool when this function is called.
%   - Auto mode (num_worker = 0) : always start with the maximum allowed workers.
%   - Explicit mode              : start with num_worker (auto-reduce if capped).

arguments
    num_worker (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    pool_type  (1,1) string {mustBeMember(pool_type, ["local","threads"])} = "local"
end

auto_mode = (num_worker == 0);

% Always shut down any existing pool first
p = gcp('nocreate');
if ~isempty(p)
    delete(p);
end

% Decide requested workers in auto mode (maximum for the selected pool_type)
if auto_mode
    switch pool_type
        case "local"
            c = parcluster('local');
            num_worker = c.NumWorkers;
        case "threads"
            num_worker = feature('numcores');
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
