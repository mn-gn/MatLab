function start_parallel(numWorkers)
    % START_PARALLEL: 병렬 풀을 개방함
    % 입력: numWorkers (선택 사항) - 사용할 코어 수
    
    % 1. 현재 병렬 풀 상태 확인
    p = gcp('nocreate');
    
    % 2. 사용 가능한 최대 코어 수 확인
    myCluster = parcluster('local');
    maxCores = myCluster.NumWorkers;
    
    % 3. 입력 인자가 없거나 최대치를 초과하면 최대 코어로 설정
    if nargin < 1 || isempty(numWorkers)
        numWorkers = maxCores;
    elseif numWorkers > maxCores
        fprintf('경고: 요청한 코어 수(%d)가 가용 최대치(%d)보다 많습니다. %d개로 조정합니다.\n', ...
                numWorkers, maxCores, maxCores);
        numWorkers = maxCores;
    end
    
    % 4. 병렬 풀 제어
    if isempty(p)
        fprintf('--- 병렬 컴퓨팅 활성화 (%d/%d 코어) ---\n', numWorkers, maxCores);
        parpool('local', numWorkers);
    else
        % 이미 풀이 열려있는데 코어 수가 다른 경우, 기존 풀 닫고 새로 열기
        if p.NumWorkers ~= numWorkers
            fprintf('기존 병렬 풀(%d 코어)을 닫고 새 풀(%d 코어)을 개방합니다...\n', ...
                    p.NumWorkers, numWorkers);
            delete(p);
            parpool('local', numWorkers);
        else
            fprintf('요청하신 %d개의 코어가 이미 작동 중입니다.\n', numWorkers);
        end
    end
end