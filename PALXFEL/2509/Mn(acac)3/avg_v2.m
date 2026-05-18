clc; clear; close all;


S = zeros(2000,22);
W = zeros(4,22);

for i = 3 : 6
    path = sprintf("Z:\\PAL_XFEL_2509\\ssd2\\ue_250907_FXL\\scratch\\resultsTRXL\\Mn_MeOH\\run00%d",i);
    avg_path = path + "\DiffAve\";
    DIR = dir(fullfile(avg_path,"*.txt"));
    pattern = '[+-]?\d+\.?\d*[eE][+-]?\d+';
    names = string({DIR.name});
    extractedStrings = regexp(names, pattern, 'match', 'once');
    timeValues = str2double(extractedStrings);
    [sortedTimes, sortIdx] = sort(timeValues);
    sortedFileNames = names(sortIdx);
    nshot_path = path + "\diffMat\nshots.dat";
    nshot = readmatrix(nshot_path);

    W(i-2,:) = nshot(:,2);

    for j = 1 : length(names)
        subpath = avg_path+sortedFileNames(j);
        A = readmatrix(subpath);
        q = A(:,1);
        dS = A(:,2);
        norm_f = sum(dS(q>=1 & q<=7));
        norm_S = dS/norm_f;
        w = nshot(j,2);

        S(:,j) = S(:,j) + norm_S*w;

    end
    
end

S = S./sum(W,1);
t = sortedTimes';
colNames = "T_" + string(t); 
finalTable = array2table(S, 'VariableNames', colNames);
finalTable = [table(q, 'VariableNames', "q"), finalTable];

T_q_cut = finalTable(q>=1 & q<=7,:);
Table = table2array(T_q_cut(:,2:end));
% [U s V] = svd(Table);
% 
% sV = diag(s)'.*V;

%% PEPC
heat_data = readmatrix("C:\Users\ming0\Documents\카카오톡 받은 파일\MeOH_comps.dat");

q_heat = heat_data(:,1);
heat_basis = heat_data(:,2:end);


 % for i = 1:size(heat_basis, 2)
 %        heat_basis(:,i) = heat_basis(:,i).*q_heat;
 % end

heat_basis_interp = interp1(q_heat,heat_basis,T_q_cut.q,"linear","extrap");

heat_basis_interp_smooth = smoothdata(heat_basis_interp,1,"movmean",10);

 % plot_range = 1:3;
 %    for i = plot_range
 %        figure
 %        hold on
 %        plot(q, heat_basis_interp(:, i), 'b-')
 %        plot(q, heat_basis_interp_smooth(:, i), 'r-')
 %        hold off
 %    end

 heat_basis_interp = heat_basis_interp_smooth;

 A_correct = zeros(size(Table,1),size(Table,2),size(heat_basis,2));
 heating = A_correct;

 num_heat_basis = size(heat_basis,2);
 for j = 1 : num_heat_basis
     for i = 1 : length(t)
         x = lsqr(heat_basis_interp(:,1:j),Table(:,i));
         arti = heat_basis_interp(:,1:j)*x;
         A_correct(:,i,j) = Table(:,i)-arti;
         heating(:,i,j) = arti;
     end
 end
%%
% residual = Table;
% 
% for j = 1 : num_heat_basis
%     for i = 1 : length(t)
%         % 1. 현재 단계의 단일 기저(basis)만 선택
%         current_basis = heat_basis_interp(:, j);
% 
%         % 2. 원본 데이터가 아닌, '이전 단계까지 제거되고 남은 잔차'에 대해 피팅
%         x = lsqr(current_basis, residual(:, i));
% 
%         % 3. j번째 기저에 의해 현재 단계에서 새롭게 추정된 아티팩트
%         current_arti = current_basis * x;
% 
%         % 4. 잔차 업데이트 (누적 제거의 핵심): 현재 잔차에서 방금 구한 아티팩트를 뺌
%         residual(:, i) = residual(:, i) - current_arti;
% 
%         % 5. 결과 저장
%         % A_correct에는 누적 제거되고 남은 최종 잔차를 저장
%         A_correct(:, i, j) = residual(:, i);
% 
%         % heating에는 지금까지 제거된 아티팩트의 총합을 누적하여 저장
%         if j == 1
%             heating(:, i, j) = current_arti;
%         else
%             heating(:, i, j) = heating(:, i, j-1) + current_arti;
%         end
%     end
% end
%%
  U = zeros(length(T_q_cut.q), length(t), num_heat_basis);
    S = zeros(length(t), length(t), num_heat_basis);
    V = zeros(length(t), length(t), num_heat_basis);
    
    for i = 1:num_heat_basis
        [U(:, :, i), S(:, :, i), V(:, :, i)] = svd(A_correct(:, :, i), "econ");
    end

 
    q = T_q_cut.q;

    % num_rsv_sel = 2;
    % figure
    % hold on
    % legends = strings(1, num_heat_basis + 1);
    % for i = 1:num_heat_basis
    %     plot(t, V(:, num_rsv_sel, i), '.-')
    %     legends(i) = sprintf('%d components eliminated', i);
    % end

    
    
    for num_heat_sel = 1 : 1
    
    scale = [-0.02, 0.02];
    
    figure('Color','white','pos',[10 10 800 800]);
    surf(q,t, q'.*Table','EdgeColor', 'None', 'facecolor', 'interp');
    set(gca,'Ydir','reverse','TickDir','out')
    colormap jet
    zlim(scale)
    clim(scale)
    view(2);
    xlim([q(1) q(end)]);
    ylim([t(1) t(end)]);
    xlabel('q')
    ylabel('t')
    title('Original q\DeltaS(q, t)')


    % figure('Color','white','pos',[10 10 800 800]);
    % surf(q,t, q.*heating(:, :, num_heat_sel)','EdgeColor', 'None', 'facecolor', 'interp');
    % set(gca,'Ydir','reverse','TickDir','out')
    % colormap jet
    % zlim(scale)
    % clim(scale)
    % view(2);
    % xlim([q(1) q(end)]);
    % ylim([t(1) t(end)]);
    % xlabel('q')
    % ylabel('t')
    % title('Heating q\DeltaS(q, t)')


    figure('Color','white','pos',[10 10 800 800]);
    surf(q,t, q'.*A_correct(:, :, num_heat_sel)','EdgeColor', 'None', 'facecolor', 'interp');
    set(gca,'Ydir','reverse','TickDir','out')
    colormap jet
    zlim(scale)
    clim(scale)
    view(2);
    xlim([q(1) q(end)]);
    ylim([t(1) t(end)]);
    xlabel('q')
    ylabel('t')
    title(['PEPC-treated q\DeltaS(q, t) heat sel :', num2str(num_heat_sel)]);

    end