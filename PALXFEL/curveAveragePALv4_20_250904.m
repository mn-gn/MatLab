%% Section 1 | Specification
parentFolderH5='/xfel/ffs/dat/scan';
% parentFolderH5='/xfel/ffs/dat/ue_250401_FXL/rawData';
parentFolderResults='/xfel/ffs/dat/ue_250918_FXL/scratch/resultsTRXL'; % Parent folder of results folders.

h5header = 'Pt_phenyl';


runNumber = []; %AutoRun

% I3m MeCN
normRange =[4.0 6.0];
SVDRange = [1.0 6.0];
AUCRange = [1.0 2.3];
autoAUCValue=0.01; % Manual AUC setting;


%%% operation Options
% XFEL options
XrayCenter = 14.983; % (KeV) 250905
repetitionRate = 60; % Hz(1/s)

% Saving options

% outputWriteOn=false; % Operon
outputWriteOn=true; % Operon
figureSave=true;

eachScanSave=false; % true or 1 means making (scan numer)_diff_av_(time-delays) file in scan_process of each run directory
diffsSave=false;

% pumpI0 options
pumpI0On=false;
pumpI0Correction=false;
correctionFunc=@(pumpI0) (7.5e6*pumpI0+0.0920);
% correctionFunc=@(pumpI0) ( (pumpI0-min(pumpI0)-(max(pumpI0)-min(pumpI0))*4) / (min(pumpI0)-(max(pumpI0)-min(pumpI0))*4) );

OXCcorrection = [8.8e-22 -3.1e-18 2e-15 6.5e-13]; %%250830 | 16:46
isOXCcorrection = false;
%isOXCcorrection = true;
tdOXCcorrect = 2.5e-12;
binOXCcorrect = 20e-15;

% Draw options
drawRange=[0.5 7.0]; % Draw q range (2Dscan diffAUC MAX)
plotOption='qSq';
timedelayDrawPoint = 3e-12;

mode = 'TRXL'; % 'TRXL','2Dscan' But autoModer will be operate
config = encapsulatorConfig(...
    parentFolderH5,parentFolderResults,h5header,mode,...
    normRange,SVDRange,AUCRange,autoAUCValue,XrayCenter,repetitionRate,...
    outputWriteOn,figureSave,eachScanSave,diffsSave,pumpI0On,pumpI0Correction,...
    correctionFunc,drawRange,plotOption,timedelayDrawPoint,...
    OXCcorrection,isOXCcorrection,tdOXCcorrect,binOXCcorrect...
    );
% config.forceLatestRun = false;
config.forceLatestRun = true;

% config=figureEncapsulator(config);
welcomeComment;

%% AutoRun
if isempty(runNumber)
    pauseTime = 3; % Wating time for new analysis
    % minimumH5File=2; %% Minimum chipal number for anlysis that run
    config.isAutoRun = true;
    autorunCurveAveragePAL(pauseTime,config);
end
%% ManualRun
for idx = runNumber
    config.isAutoRun = false;
    config.folderh5=sprintf('%s/%s_%.5d_DIR',parentFolderH5,h5header,idx);
    if ~exist('h5Struct')
        h5StructPre=struct;
    elseif exist('h5Struct')
        h5StructPre=h5Struct;
    end
    [h5Struct,diffStruct,config]=curveAveragePAL(h5StructPre,config);
end
%% MainFunction
function [h5Struct,diffStruct,config]=curveAveragePAL(h5StructPre,config)
% config.debug_numH5 = 10; % This is for debug, please config.debug_numH5 = []; for typical data
config.debug_numH5 = []; % This is for debug, please config.debug_numH5 = []; for typical data

% Additional optionsS
config.isDiffAveBin = false;
config.isDiffAve = true;
config.autoModerOn=true;
% config.verbose=false;
config.verbose=true;
config.numInt=4;
% config.numInt=5;
config.plotEachScan=false;
config.figResolution=150;
config.noUpsideDown=true;
config.scatterdROI = config.AUCRange;
config.I0cutoff = 1e-10; % eh1qbpm1 %%%%%%% in ci den te beam drop


config.channelROI{1} = {'eh1qbpm1_totalsum',1,[]}; % name | dataType : (1 = scalar) | memory
% config.channelROI{1} = {'ohqbpm2_totalsum',1,[]}; % name | dataType : (1 = scalar) | memory
config.channelROI{2} = {'eh1pdcom_channel0p',1,[]}; % name | dataType : (1 = scalar) | memory
config.channelROI{3} = {'eh1rayMXAI_int',2,[]}; % name | dataType : (2 = matrix) | memory
config.channelROI{4} = {'eh1rayMXAI_tth',0,[]}; % name | dataType : (2 = matrix) | memory
config.channelROI{5} = {'pulseInfo',0,[]}; % name | dataType : (1 = scalar) | memory
config.channelROI{6} = {'scanInfo',0,[]}; % name | dataType : (0 = scanInfo) | memory

if config.isOXCcorrection
    config.channelROI{7} = {'eh1oxc_pos',1,[]}; % name | dataType : (1 = scalar) | memory
end


% Decapsulator
isOXCcorrectionOrg = config.isOXCcorrection;

outputWriteOn=config.outputWriteOn;
eachScanSave=config.eachScanSave;
figureSave=config.figureSave;
diffsSave=config.diffsSave;

% Start time recorder
timeStart=sprintf('%s',datetime('now'));

headerNameCheck(config);
[
    config.q,...
    config.numH5,...
    h5Struct,...
    ]=readH5Struct(h5StructPre,config);
clearvars h5StructPre;

% [
%     config.numH5,...
%     h5Struct,...
%     ]=fixAppendingErrorH5Strcut(h5Struct,config);
[
    config.mode,...
    config.run,...
    config.posi1List,...
    config.posi2List,...
    config.posi3List,...
    config.numShotList,...
    config.numData,...
    config.td_org,...
    h5Struct,...
    ]=readH5ScanInfo(h5Struct,config);

switch config.mode
    case '2Dscan'
        config.isOXCcorrection = false;
        fprintf('  2Dscan mode | turn off the OXC correction \n');
end

if config.isOXCcorrection
    [
        h5Struct,...
        ] = getTimingJitterFromOXCpos(h5Struct,config);
end

[
    config.qIsAUC,...
    config.qIsNorm,...
    config.qIsSVD,...
    h5Struct,...
    ]=calAUCH5(h5Struct,config);
[
    h5Struct,...
    ]=getShotPair(h5Struct,config);
[
    h5Struct,...
    ]=getDiffAUC(h5Struct,config);
[
    h5Struct,...
    ]=goodSignalCheckChiStruct(h5Struct,config);
config.cutoffAUC = config.autoAUCValue;
[
    h5Struct,...
    ]=getDiffAveSTD(h5Struct,config);



%%%%%%%%%%%%%%%%%%%%%%%% Dr Strange %%%%%%%%%%%%%%%%%
[config.weAreEndGameNow,config.progress] = DrStrange(config.numData,length(h5Struct));

if ~isempty(config.debug_numH5)
    fprintf('  Bebug mode | We must end this game.\n')
    config.weAreEndGameNow = true;
end


if all(config.td_org > config.tdOXCcorrect) || (length(config.td_org) < 10)
    config.isOXCcorrection = false;
    fprintf('  turn off All timedelay is above the OXCcorrect ...\n');
end

if config.isOXCcorrection && config.isAutoRun
    if config.progress < 0.4
        config.isOXCcorrection = false;
        fprintf('  turn off the OXC correction ... | progress %.1f %% \n',config.progress*100);
    else
        fprintf('  OXC correction now operating !!! | progress %.1f %% \n',config.progress*100);
    end
end
fprintf('---------------------------------------------------------\n\n');
fprintf('  isOXCcorrection | %d\n\n',config.isOXCcorrection);
fprintf('---------------------------------------------------------\n');

% 
% if config.progress > 0.4 || ~config.isAutoRun
%     if config.isOXCcorrection
%         fprintf('  OXC correction now operating !!! | progress %.1f %% \n',config.progress*100);
%     end
%     else
%     fprintf('  turn off the OXC correction ... | progress %.1f %% \n',config.progress*100);
%     config.isOXCcorrection = false;
% end



%%%%%%%%%%%%%%%%%%%%%%%%%% Make diff Str.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch config.mode
    case '2Dscan'
        [
            diffStruct,...
            ]=setDiffStructWithMotor(h5Struct,config);
        [
            diffStruct,...
            ]=getFileIndex3DscanDiffStruct(diffStruct,config);
        [
            diffStruct,...
            ]=getAveSTDpumpI0(h5Struct,diffStruct,config);
        % [
        %     diffStruct,...
        %     numBadCurve,...
        %     ]=badDiffStructRemove(diffStruct,config);
    case 'TRXL'
        %%%%%%%%% no OXC correciton %%%%%%%%%%%
        [
            config.timedelay,...
            config.numTimedelay,...
            diffStruct,...
            ]=setTimedelayDiffStruct(h5Struct,config);
        [
            diffStruct,...
            ]=getFileIndexDiffStruct(h5Struct,diffStruct,config);
        [
            diffStruct,...
            ]=getAveSTDpumpI0(h5Struct,diffStruct,config);
        % [
        %     diffStruct,...
        %     numBadCurve,...
        %     ]=badDiffStructRemove(diffStruct,config);
        [
            config.eh1rayMXAI_intMean,...
            ]=geteh1rayMXAI_intMean(h5Struct);
        %%%%%% get diffBinStruct %%%%%%%%%
        if config.isOXCcorrection
            [
                diffBinStruct,...
                ] = getdiffBinStruct(h5Struct,diffStruct,config);
            % [
            %     diffBinStruct,...
            %     ] = removeBadPointDiffStruct(diffBinStruct);
            tdOXCcorrectOrg = config.tdOXCcorrect; 
            config.tdOXCcorrect = 1e12;
            [
                diffBinAllStruct,...
                ] = getdiffBinStruct(h5Struct,diffStruct,config);
            config.tdOXCcorrect = tdOXCcorrectOrg;
        end
        %%%%%% post process of diffStruct %%%%%%%%
        if config.isOXCcorrection
            [
                config.diffAllwoOXC,...
                config.timedelaywoOXC,...
                config.numTimedelaywoOXC,...
                ]=makeDiffAll(diffStruct);
            [
                config.UwoOXC,...
                config.SwoOXC,...
                config.VwoOXC,...
                ]=doSVD(config.diffAllwoOXC,config);
            [
                config.diffAll,...
                config.timedelay,...
                config.numTimedelay,...
                ]=makeDiffAll(diffBinStruct);
            [
                config.U,...
                config.S,...
                config.V,...
                ]=doSVD(config.diffAll,config);
        else
            [
                config.diffAll,...
                config.timedelay,...
                config.numTimedelay,...
                ]=makeDiffAll(diffStruct);
            [
                config.U,...
                config.S,...
                config.V,...
                ]=doSVD(config.diffAll,config);
        end
    otherwise
        fprintf('  "mode" List | "TRXL", "2Dscan"\n')
        error('       "Mode" name error. Must check mode name.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Draw Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotStartComment;
switch config.mode
    case '2Dscan'
        [
            config.scan2DFig,...
            ]=plot2Dscan(11111,diffStruct,h5Struct,config);
    case 'TRXL'
        [
            config.fig2d,...
            config.fig1d,...
            ]=drawDifferenceScatteringCurves(33333,config);
        [
            config.SVFig,...
            ]=plotSVD(66667,config);

        [
            config.SCUCVFig,...
            ] = plotSCUCV(99999,config);
        [
            config.OXCstaticFig,...
            ] = plotOXCstatic(11312,h5Struct);
        if config.isOXCcorrection
            [
                config.OXCcompFig,...
                ] = plotOXCComparison(11311,config);
        end
        if config.pumpI0On
            [
                config.pumpI0Fig,...
                ]=plotPumpI0(77777,diffStruct,config);
        end
        if config.plotEachScan
            [
                config.lastScanSVFig,...
                ]=drawLastScanSV(111,h5Struct,config);
        end



end
try
    [
        config.cutoffFig,...
        config.numGood,...
        config.R2,...
        config.I0Fig,...
        ]=plotCutoff(h5Struct,config);
catch
end
plotDoneComment;
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%% Write output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch config.mode
    case '2Dscan'
        DiffavSave=false;

        SVDSave=false;
        eachScanSave=false;
        diffsSave=false;
    case 'TRXL'
        DiffavSave=true;
        SVDSave=true;
end

% config.weAreEndGameNow = DrStrange(config.numData,length(h5Struct)-numBadCurve);
% 
% if ~isempty(config.debug_numH5)
%     fprintf('  Bebug mode | We must end this game.\n')
%     config.weAreEndGameNow = true;
% end

if outputWriteOn
    [
        config.folderResults,...
        config.folderDiffAve,...
        config.folderDiffAveCorrection,...
        config.folderEachScan,...
        config.folderDiffs,...
        config.folderValue,...
        config.folderImg,...
        config.folderSVD...
        ]=outputFolderAllocator(config);


    writeSpecification(config);
    try
        writeRsqaure(config);
    catch
    end

    if SVDSave
        writeSVDresults(config);
    end
    if eachScanSave
        writeEachScan(h5Struct,config);
    end
    if diffsSave
        writeDiffs(h5Struct,config);
    end
    % if config.isAutoRun && ~config.weAreEndGameNow
    %     fprintf('  Still not a endgame yet\n');
    % else
    %     fprintf('  We are the endgame now. Start saving\n');
    if DiffavSave
        if config.isOXCcorrection
            binVal = round(config.binOXCcorrect*1e15);

            if config.isDiffAveBin
                outputName = sprintf('DiffAve%d',binVal);
                writeOutputDiffav(diffBinAllStruct,outputName,config);
            end

            outputName = sprintf('diffMat%d',binVal);
            writeOutputDiffavMatrix(diffBinAllStruct,outputName,config);

            if config.isDiffAve
                writeOutputDiffav(diffStruct,'DiffAve',config);
            end
            writeOutputDiffavMatrix(diffStruct,'diffMat',config);
        else
            if config.isDiffAve
                writeOutputDiffav(diffStruct,'DiffAve',config);
            end
            writeOutputDiffavMatrix(diffStruct,'diffMat',config);
        end

    end

    if DiffavSave && config.pumpI0Correction
        if config.isOXCcorrection
            writeOutputDiffavCorrection(diffBinStruct,config);
        else
            writeOutputDiffavCorrection(diffStruct,config);
        end
    end

    writeStaticMean(h5Struct,config);
    writeValue(h5Struct,config);

    % imgSaveFlag = config.weAreEndGameNow;

    if figureSave
        if ~config.isAutoRun || config.weAreEndGameNow
            fprintf('  Save figures in progress ...\n');
            writeImg(config);
        else
            fprintf('  Skip saving figure\n');
        end
    end
    % end
end
config.isOXCcorrection = isOXCcorrectionOrg;
finishComment(timeStart);
end
%% SubFunction | Same order as in the main function
function bins = makeZeroCenteredBins(X, binSize)
% X: 1×N (row vector)
x = X(:);
xmin = min(x); xmax = max(x);

kmin = floor((xmin + binSize/2)/binSize);
kmax = ceil((xmax - binSize/2)/binSize);
centers = (kmin:kmax) * binSize;   % 0 포함
edges   = (centers(1)-binSize/2) : binSize : (centers(end)+binSize/2);

binIdx = discretize(x, edges, 'IncludedEdge','left'); % [left,right)
bins.edges   = edges;
bins.centers = centers(:);
bins.binIdx  = binIdx;   % 각 column이 속한 bin 번호
bins.nbins   = numel(centers);
end
function acc = initBinAccumulator(nbins, M)
acc.n    = zeros(nbins,1,'uint32');
acc.mean = zeros(nbins, M);
acc.M2   = zeros(nbins, M);
end

function acc = updateBinAccumulator(acc, k, colVec)
% colVec: M×1 column vector
if isnan(k), return; end
n0 = acc.n(k);
n1 = n0 + 1;
acc.n(k) = n1;

v = colVec(:)'; % 1×M row
delta  = v - acc.mean(k,:);
acc.mean(k,:) = acc.mean(k,:) + delta / double(n1);
delta2 = v - acc.mean(k,:);
acc.M2(k,:)   = acc.M2(k,:) + delta .* delta2;
end

function out = finalizeBins(bins, acc, binSize)
% FINALIZEBINS
%   - bins: struct from makeZeroCenteredBins
%       .centers  (nbins×1)
%       .binIdx   (N×1)
%   - acc: struct from init/updateBinAccumulator
%       .n     (nbins×1)   % uint32 sample counts per bin
%       .mean  (nbins×M)   % running means per bin (Welford)
%       .M2    (nbins×M)   % sum of squared deviations per bin (Welford)
%   - binSize: scalar (bin width)
%
% Returns
%   out.edges    : (nb_nonempty+1)×1  reconstructed edges (non-contiguous bins allowed)
%   out.centers  : nb_nonempty×1      bin centers (0 포함)
%   out.mean     : nb_nonempty×M      mean curves per bin
%   out.sem      : nb_nonempty×M      SEM = sqrt( variance / n )
%   out.count    : nb_nonempty×1      sample count per bin
%   out.binIdx   : N×1                original bin index per sample (unchanged)
%   out.nonEmptyMask : nbins×1        mask of bins kept

% 1) 기본 정보
nb      = numel(bins.centers);
counts  = double(acc.n(:));            % nb×1
nonEmpty = counts > 0;                  % keep only bins with samples
nEff     = counts(nonEmpty);            % nb_nonempty×1

% 2) 평균/분산 원천값을 먼저 잘라낸다 (연속 괄호 인덱싱 금지!)
mean_non = acc.mean(nonEmpty, :);       % nb_nonempty×M
M2_non   = acc.M2(nonEmpty, :);         % nb_nonempty×M

% 3) 분산/SEM 계산 (Welford)
%    unbiased variance: var = M2/(n-1), for n>=2; else NaN
varMat = nan(size(M2_non));             % nb_nonempty×M
ge2    = (nEff >= 2);                   % 논리마스크

if any(ge2)
    denom = (nEff(ge2) - 1);            % nb_ge2×1
    % R2016b+ : 암시적 확장
    varMat(ge2, :) = M2_non(ge2, :) ./ denom;
    % 구버전 호환이 필요하면 다음 라인 사용:
    % varMat(ge2, :) = bsxfun(@rdivide, M2_non(ge2, :), denom);
end

semMat = nan(size(varMat));
if any(ge2)
    semMat(ge2, :) = sqrt(varMat(ge2, :) ./ nEff(ge2));
    % 구버전:
    % semMat(ge2, :) = sqrt(bsxfun(@rdivide, varMat(ge2, :), nEff(ge2)));
end

% 4) 남는 centers와 edges 재구성
centers_kept = bins.centers(nonEmpty);  % nb_nonempty×1
% 비어있는 bin을 제거했으므로 edges도 해당 centers만을 기준으로 재작성
% (연속성이 보장되지는 않지만, 각 center의 좌/우 경계로 충분)
edges = [centers_kept - binSize/2; centers_kept(end) + binSize/2];

% 5) 결과 포장
out.edges        = edges(:);
out.centers      = centers_kept(:);
out.mean         = mean_non';
out.sem          = semMat';
out.count        = nEff(:);
out.binIdx       = bins.binIdx;         % 원본 bin 번호는 유지 (필요시 재매핑 가능)
out.nonEmptyMask = nonEmpty(:);
end









function config=encapsulatorConfig(...
    parentFolderH5,...
    parentFolderResults,...
    h5header,...
    mode,...
    normRange,...
    SVDRange,...
    AUCRange,...
    autoAUCValue,...
    XrayCenter,...
    repetitionRate,...
    outputWriteOn,...
    figureSave,...
    eachScanSave,...
    diffsSave,...
    pumpI0On,...
    pumpI0Correction,...
    correctionFunc,...
    drawRange,...
    plotOption,...
    timedelayDrawPoint,...
    OXCcorrection,...
    isOXCcorrection,...
    tdOXCcorrect,...
    binOXCcorrect...
    )
config.h5header=h5header;
config.parentFolderH5=parentFolderH5;
config.parentFolderResults=parentFolderResults;
config.mode=mode;
config.XrayCenter=XrayCenter;
config.repetitionRate=repetitionRate;
config.normRange=normRange;
config.SVDRange=SVDRange;
config.AUCRange=AUCRange;
config.drawRange=drawRange;
config.autoAUCValue=autoAUCValue;
config.outputWriteOn=outputWriteOn;
config.eachScanSave=eachScanSave;
config.pumpI0On=pumpI0On;
config.pumpI0Correction=pumpI0Correction;
config.correctionFunc=correctionFunc;
config.figureSave=figureSave;

config.diffsSave=diffsSave;
config.plotOption=plotOption;
config.timedelayDrawPoint=timedelayDrawPoint;

config.OXCcorrection = OXCcorrection;
config.isOXCcorrection = isOXCcorrection;
config.tdOXCcorrect = tdOXCcorrect;
config.binOXCcorrect = binOXCcorrect;

end
function welcomeComment
% PPiAAc;
fprintf('---------------------------------------------------------\n');
fprintf('  curveAveragePALv4.m\n\n')
fprintf('  Made by Jungmin Kim\n')
fprintf('  jkim9486@gmail.com\n\n')
fprintf('  Institute for Basic Science\n')
fprintf('  Center for Advanced Reaction Dynamics\n')
TodaysProverb;
fprintf('---------------------------------------------------------\n');
end
function headerNameCheck(config)
parentFolderH5=config.parentFolderH5;
h5header=config.h5header;
folderh5=config.folderh5;

dirName=sprintf('%s/%s*/scanInfo/*.h5',parentFolderH5,h5header);
dirList=dir(dirName);
if isempty(dirList)
    % fprintf('  H5 PATH\n  cd %s\n\n',folderh5);
    fprintf('  ERROR !!! ERROR !!! h5header Name IS WRONG !!! FIRE !!!\n');
    fprintf('  ERROR !!! ERROR !!! h5header Name IS WRONG !!! FIRE !!!\n');
    fprintf('  ERROR !!! ERROR !!! h5header Name IS WRONG !!! FIRE !!!\n');
    fprintf('  ERROR !!! ERROR !!! h5header Name IS WRONG !!! FIRE !!!\n');
    fprintf('  ERROR !!! ERROR !!! h5header Name IS WRONG !!! FIRE !!!\n');
    return;
else
    % fprintf('  H5 PATH\n  cd %s\n\n',folderh5);
end
end
function [q,numH5,h5Struct]=readH5Struct(h5StructPre,config)
%% Section 0 | Decapsulation
XrayCenter = config.XrayCenter;
h5header = config.h5header;
repetitionRate = config.repetitionRate;
folderh5 = config.folderh5;
verbose = config.verbose;
noUpsideDown = config.noUpsideDown;
debug_numH5 = config.debug_numH5;
channelROI = config.channelROI;
%% Section 1 | Read q
tic
% Read q H5 file
% checker =1;
while true
    try
        pathh5=sprintf('%s/eh1rayMXAI_tth',folderh5);
        namePattern=sprintf('%s/*.h5',pathh5);
        h5Struct=dir(namePattern);
        fileName=sprintf('%s/%s',pathh5,h5Struct(1).name);
        h5QFile=h5info(fileName);
        break;
    catch
        fprintf('  Wait for the first h5 files\n');
        pause(3);
    end
end


timestampPID=h5QFile.Datasets(1).Name;
timestampPIDChar=sprintf('/%s',timestampPID);
twoTheta=h5read(fileName,timestampPIDChar);
twoThetaRadian = twoTheta/180 * pi;
XrayWaveLength=12.398/XrayCenter;
q = 4*pi*sin(twoThetaRadian/2)/XrayWaveLength;
%% Section 2 | Read premative pulse Info
% pathh5=sprintf('%s/pulseInfo',folderh5);
pathh5=sprintf('%s/eh1rayMXAI_int',folderh5);
namePattern=sprintf('%s/*.h5',pathh5);
h5Struct=dir(namePattern);
if ~isempty(debug_numH5)
    if length(h5Struct) < debug_numH5
        debug_numH5 = length(h5Struct);
    end
    h5Struct = h5Struct(1:debug_numH5);
end
numH5 = length(h5Struct);
%% Section 3 | Inherit from h5StructPre
[h5PairList,h5ExistIdx,h5NoExistIdx]=getH5ProcessListWithSkip(h5Struct,h5StructPre);

% h5NoExistIdx = h5NoExistIdx(1:end-1);


for idx=1:length(h5ExistIdx)
    idx_h5 = h5ExistIdx(idx);
    idx_Pre = h5PairList(idx);
    h5Struct(idx_h5).PID=h5StructPre(idx_Pre).PID;
    h5Struct(idx_h5).timestampPIDString=h5StructPre(idx_Pre).timestampPIDString;
    h5Struct(idx_h5).laserOn=h5StructPre(idx_Pre).laserOn;
    h5Struct(idx_h5).position1=h5StructPre(idx_Pre).position1;
    h5Struct(idx_h5).position2=h5StructPre(idx_Pre).position2;
    h5Struct(idx_h5).position3=h5StructPre(idx_Pre).position3;
    h5Struct(idx_h5).numShot=h5StructPre(idx_Pre).numShot;
    h5Struct(idx_h5).run=h5StructPre(idx_Pre).run;

    for cdx = 1:length(channelROI)
        dataType = channelROI{cdx}{2};
        if dataType ~= 0
            channelName = channelROI{cdx}{1};
            h5Struct(idx_h5).(channelName) = h5StructPre(idx_Pre).(channelName);
        end
    end
end
clearvars h5StructPre;
%% Section 4 | Read data all
for idx=h5NoExistIdx
    %% Section 4-1 | Check data exitance
    name=h5Struct(idx).name;
    folder=h5Struct(idx).folder;
    fileName=sprintf('%s/%s',pathh5,name);

    if verbose
        fprintf('  Reading h5 | %s\n',name);
    end
    while 1
        try
            h5File=h5info(fileName);
            break;
        catch
            pauseTime=3;
            pause(pauseTime)
            fprintf('  %s Reading in progress. Pause %d\n',name,pauseTime);
        end
    end

    %% Section 4-2 | Allocation for shot by shot
    numShot = length(h5File.Datasets);
    timestampPIDString=strings(1,numShot); % test
    for cdx = 1:length(channelROI)
        if channelROI{cdx}{2} ~= 0 % 1 or 2 allocation empty memory
            channelROI{cdx}{3} = [];
        end
    end

    %% Section 4-3 | shot-by-shot reading
    for shotIdx=1:numShot
        try
            timestampPID=h5File.Datasets(shotIdx).Name;
        catch
            fprintf('  The Last timestampPID : %s\n',timestampPID);
            fprintf('  The numShot : %d\n',numShot);
            errorMessage=sprintf('Error on the %dth h5 file / %d shots.',idx,shotIdx);
            error(errorMessage); %#ok<*SPERR>
        end
        timestampPIDChar=sprintf('/%s',timestampPID);
        PID(shotIdx)=getPID(timestampPIDChar);
        textBuff=sprintf('/%s',timestampPID);
        timestampPIDString(shotIdx)=string(textBuff);
        for cdx = 1:length(channelROI)
            switch channelROI{cdx}{2}
                case {1,2}
                    filename = sprintf('%s/%s/%s',folderh5,channelROI{cdx}{1},name);
                    try
                        valNow = h5read(filename,timestampPIDChar);
                    catch
                        switch config.channelROI{cdx}{2}
                            case 1
                                valNow = 0;
                            case 2
                                valNow = zeros(2000,1);
                            case 3
                                valNow = [];
                        end
                    end
                    channelROI{cdx}{3} = [channelROI{cdx}{3} valNow];
            end
        end
    end
    %% Section 5 | Check upsideDown of signals
    PID = PID(1:numShot);
    if noUpsideDown
        laserOn=~mod(PID,360/(repetitionRate/2)); % Normal
    else
        laserOn=~~mod(PID,360/(repetitionRate/2)); % UpsideDown

        fprintf('\n       Detector UPSIDE DOWN MODE !!!\n\n');
        fprintf('--------------------------------------------------------------\n');
    end
    %% Section 6 | Complete the h5Struct
    [position1,position2,position3]=getPosition(name);
    run=getRun(folder,h5header);
    h5Struct(idx).PID=PID;
    h5Struct(idx).timestampPIDString=timestampPIDString;
    h5Struct(idx).laserOn=laserOn;
    h5Struct(idx).position1=position1;
    h5Struct(idx).position2=position2;
    h5Struct(idx).position3=position3;
    h5Struct(idx).numShot=numShot;
    h5Struct(idx).run=run;
    for cdx = 1:length(channelROI) % numChannel
        channelName = channelROI{cdx}{1};
        switch channelROI{cdx}{2}
            case {1,2}
                h5Struct(idx).(channelName)=channelROI{cdx}{3};
        end
    end
    % numShot = size(h5Struct(10).eh1rayMXAI_int,2);
    % disp('debug');
end
%% Section 7 | goodBye
numH5Process=length(h5NoExistIdx);
timeReading=toc;
dataReadingDoneComment(timeReading,numH5Process);
    function dataReadingDoneComment(timeReading,numH5Process)
        message0=sprintf('  Data reading done');
        message1=sprintf('  Total time spend reading data       | %.4fs',timeReading);
        message2=sprintf('  Average time spent reading one data | %.4fs',timeReading/numH5Process);

        fprintf('%s\n\n%s\n%s\n%s\n',message0,message1,message2);
        disp('--------------------------------------------------------------')
    end

end
function [numH5,h5Struct]=fixAppendingErrorH5Strcut(h5StructBroken,config)
% We encount the error on 2211 beamtime that combine two h5 file to one h5
% with a number of image number was double.
numH5Broken=config.numH5;

numShotList=nan(numH5Broken,1);

for idx=1:numH5Broken
    PID=h5StructBroken(idx).PID;
    numShot=length(PID);
    numShotList(idx)=numShot;
end

[~,badidx]=rmoutliers(numShotList);
goodidx=~badidx;
h5Struct=h5StructBroken(goodidx);
numH5=length(h5Struct);
if any(badidx)
    fprintf('  Bad h5 file(s) removed\n\n');
    fprintf('  Removed h5 idx | %d (ImageNumber %d) \n',[[find(badidx)] [numShotList(badidx)]]');
else
    fprintf('  There is no bad h5. GREAT!\n');
end
disp('--------------------------------------------------------------')

% Get numShotAvg
% numShotAvg=mean(numShotList);

% idxCounter=0;
% for idx=1:numH5
%     % Read inputs
%     name=h5StructBroken(idx).name;
%     folder=h5StructBroken(idx).folder;
%     filedate=h5StructBroken(idx).date;
%     eh1rayMXAI_int=h5StructBroken(idx).eh1rayMXAI_int;
%     PID=h5StructBroken(idx).PID;
%     timestampPIDString=h5StructBroken(idx).timestampPIDString;
%     laserOn=h5StructBroken(idx).laserOn;
%     position1=h5StructBroken(idx).position1;
%     position2=h5StructBroken(idx).position2;
%     position3=h5StructBroken(idx).position3;
%     run=h5StructBroken(idx).run;
%
%     numShot=numShotList(idx);
%
%     if numShot<=1.1*numShotAvg
%        idxCounter=idxCounter+1;
%        h5Struct(idxCounter)=writeH5Struct(name,folder,filedate,eh1rayMXAI_int, ...
%            PID,timestampPIDString,laserOn, ...
%             position1,position2,position3,run);
%     else
%         fprintf('  Too many shot Number\n');
%     end
%
%
%
% end

%     function h5Single=writeH5Struct(name,folder,filedate,eh1rayMXAI_int,PID, ...
%             timestampPIDString,laserOn, ...
%             position1,position2,position3,run)
%         h5Single.name=name;
%         h5Single.folder=folder;
%         h5Single.date=filedate;
%         h5Single.eh1rayMXAI_int=eh1rayMXAI_int;
%         h5Single.PID=PID;
%         h5Single.timestampPIDString=timestampPIDString;
%         h5Single.laserOn=laserOn;
%         h5Single.position1=position1;
%         h5Single.position2=position2;
%         h5Single.position3=position3;
%         h5Single.run=run;
%     end
end
function [mode,run,posi1List,posi2List,posi3List,numShotList,numData,td,h5Struct]=readH5ScanInfo(h5Struct,config)
mode=config.mode;
autoModerOn=config.autoModerOn;

folderh5=config.folderh5;
numH5=length(h5Struct);

pathh5=sprintf('%s/scanInfo',folderh5);
fileName=sprintf('%s/scanInfo.h5',pathh5);
fprintf('  Reading scanInfo in progress\n\n');

%%% autoModer %%%%%%%%%%%%%%%%%%%
switch mode
    case '2Dscan'
        try
            motor1ValueByPosition = h5read(fileName,'/run/scan00001/motor/m1');
            motor2ValueByPosition = h5read(fileName,'/run/scan00001/motor/m2');
            fprintf('  mode | 2Dscan\n');
        catch
            if autoModerOn
                timedelayByPosition = h5read(fileName,'/run/scan00001/motor/m1');
                fprintf('  autoModer changing mode | TRXL\n');
                mode='TRXL';
            else
                error('      Should check the mode name');
            end
        end
    case 'TRXL'
        try
            motor1ValueByPosition = h5read(fileName,'/run/scan00001/motor/m1');
            motor2ValueByPosition = h5read(fileName,'/run/scan00001/motor/m2');
            if autoModerOn
                fprintf('  autoModer changing mode | 2Dscan\n');
                mode='2Dscan';
            else
                timedelayByPosition = h5read(fileName,'/run/scan00001/motor/m1');
                fprintf('  2D TRXL (with scan) | TRXL\n');
            end
        catch
            timedelayByPosition = h5read(fileName,'/run/scan00001/motor/m1');
            fprintf('  mode | TRXL\n');
        end
    otherwise
        error('      Should check the mode name');
end

%%% Read scanInfo %%%%%%%%%%%%%%%%%
for idx=1:numH5
    run=h5Struct(idx).run;
    position1=h5Struct(idx).position1;
    position2=h5Struct(idx).position2;
    position3=h5Struct(idx).position3;

    switch mode
        case '2Dscan'
            h5Struct(idx).motor1=motor1ValueByPosition(position1);
            h5Struct(idx).motor2=motor2ValueByPosition(position2);
            h5Struct(idx).motor3=position3;
            name=sprintf('run%.3d_position1_%.3d_position2_%.3d.chi_pal',run,position1,position2);
        case 'TRXL'
            timeOn=timedelayByPosition(position1)*1e-12;
            scan=position2;
            h5Struct(idx).timedelay=timeOn;
            h5Struct(idx).scan=scan;
            h5Struct(idx).motor3=position3;
            name=sprintf('run%.3d_scan%.3d_position%.3d_delay_%.3f.chi_pal',run,scan,position1,timeOn*1e12);
    end

    h5Struct(idx).chiName=name;
    h5Struct(idx).run=run;
    h5Struct(idx).position1=position1;
    h5Struct(idx).position2=position2;
    h5Struct(idx).position3=position3;
end

[run,...
    posi1List,...
    posi2List,...
    posi3List,...
    numShotList,...
    ]=getListFromH5(h5Struct,config);

switch mode
    case '2Dscan'
        td = [];
    case 'TRXL'
        td = [];
        for idx = 1:length(h5Struct)
            tdNow = h5Struct(idx).timedelay;
            td = [td tdNow];
        end
end

%%% Make numMotor %%%

numData = count_numData;
% numData = count_numData(posi1List,posi2List,posi3List);

readH5ScanInfoDoneComment(run,posi1List,posi2List,numShotList,numData);
    function readH5ScanInfoDoneComment(run,posi1List,posi2List,numShotList,numData)
        message1=sprintf('  Run number                         | %d',run);
        message2=sprintf('  NumPos1 (the number of time delay) | %d',length(unique(posi1List)));
        message3=sprintf('  NumPos2 (the number of scan)       | %d',length(unique(posi2List)));
        message4=sprintf('  The number of whole position       | %d',numData);
        message5=sprintf('  AvgShotNumber(On+Off)              | %d',round(mean(numShotList,2)));

        fprintf('%s\n%s\n%s\n%s\n%s\n',message1,message2,message3,message4,message5);
        disp('--------------------------------------------------------------')
    end


% function numData = count_numData(posi1List,posi2List,posi3List)
%     pos1 = unique(posi1List);
%     pos2 = unique(posi2List);
%     pos3 = unique(posi3List);
%
%     num_pos1 = length(pos1);
%     num_pos2 = length(pos2);
%     num_pos3 = length(pos3);
%
%     numData = num_pos1*num_pos2*num_pos3;
% end
    function numData = count_numData
        switch mode
            case '2Dscan'
                numMotor1 = length(motor1ValueByPosition);
                numMotro2 = length(motor2ValueByPosition);
                numData = numMotor1*numMotro2;
            case 'TRXL'
                numData = length(timedelayByPosition);
        end
    end


end

function h5Struct = getTimingJitterFromOXCpos(h5Struct,config)
numH5 = config.numH5;
OXCcorrection = config.OXCcorrection;

for idx = 1:numH5
    h5Struct(idx).jitter = polyval(OXCcorrection,double(h5Struct(idx).eh1oxc_pos));
    % h5Struct(idx).jitter = rand(1,length(h5Struct(idx).laserOn))*100*1e-15;
    h5Struct(idx).tds = h5Struct(idx).timedelay + h5Struct(idx).jitter;

    onJitter = h5Struct(idx).jitter(h5Struct(idx).laserOn);

    [~,badidx] = rmoutliers(onJitter,'gesd');
    goodidx=~badidx;

    h5Struct(idx).OXCgood = goodidx;
end
end

% function [mode,modeConfig]=autoModer(config)
% % autoModerOn=config.autoModerOn;
%
%
%
% switch mode
%     case '2Dscan'
%     case
% end
%
% folderh5=config.folderh5;
% numH5=config.numH5;
% pathh5=sprintf('%s/scanInfo',folderh5);
% fileName=sprintf('%s/scanInfo.h5',pathh5);
% end
%
%
% motor1ValueByPosition = h5read(fileName,'/run/scan00001/motor/m1');
% motor2ValueByPosition = h5read(fileName,'/run/scan00001/motor/m2');
%
%
%
%
% %%% old %%%%%
% autoModerOn=false;
% modeOrg=config.mode;
% if autoModerOn
% posi2List=config.posi2List;
% numPos2=length(unique(posi2List));
% if numPos2>1
%     mode='2Dscan';
%     fprintf('  autoModer Operated\n');
%     fprintf('  NumPos2 | %d\n',numPos2);
%     fprintf('  Original mode | %s\n',modeOrg);
%     fprintf('  Modified mode | %s\n',mode);
% end
% modeConfig=mode;
% end
% mode=modeOrg;
% modeConfig=modeOrg;
% end
function [qIsAUC,qIsNorm,qIsSVD,h5Struct]=calAUCH5(h5Struct,config)
% Decapsulation
AUCRange=config.AUCRange;
normRange=config.normRange;
SVDRange=config.SVDRange;
q=config.q;
numH5=config.numH5;
I0cutoff = config.I0cutoff;
channelROI = config.channelROI;

[qIsAUC,qIsNorm,qIsSVD] = indexingQ(q,AUCRange,normRange,SVDRange);
for idx=1:numH5
    % [h5Struct(idx).AUC,h5Struct(idx).normfactor]=getAUCandNorm(h5Struct(idx), ...
    %     qIsAUC,qIsNorm,q);
    dataMatrix = h5Struct(idx).eh1rayMXAI_int;
    h5Struct(idx).AUC = AUC_Jkim(dataMatrix,q,AUCRange);
    h5Struct(idx).normfactor = AUC_Jkim(dataMatrix,q,normRange);
    h5Struct(idx).I0good = I0cutoff<h5Struct(idx).(channelROI{1}{1});
end
end
function h5Struct=getShotPair(h5Struct,config)
fprintf('  Pairing in progress\n');
numH5=config.numH5;
% I0cutoff = config.I0cutoff;
% channelROI = config.channelROI;

for idx=1:numH5
    laserOn = h5Struct(idx).laserOn;

    % I0 = h5Struct(idx).(channelROI{1}{1});
    I0good = h5Struct(idx).I0good;
    %     numShot=h5Struct(idx).numShot;
    AUC=h5Struct(idx).AUC;

    laserOff=~laserOn;
    try
    AUCOn=AUC(laserOn); %I0
    catch
        disp('')
    end
    AUCOff=AUC(laserOff&I0good);
    % AUCOff=AUC(laserOff&I0cutoff>I0);
    numShotOn=sum(laserOn);

    shotPairBuff=nan(1,numShotOn);
    for shotIdx=1:numShotOn
        
        %         refAUC=h5Struct(h5Struct(idx).refFile).AUC; %#ok<PFBNS>
        %         AUCDifference=abs(AUCOff-h5Struct(idx).AUC(shotIdx));
        if isempty(AUCOff)
            shotPairBuff(shotIdx) = shotIdx;
        else
            try

                AUCDifference=abs(AUCOff-AUCOn(shotIdx));
                AUCMinIndex=find(AUCDifference==min(AUCDifference));
                % try
                shotPairBuff(shotIdx)=AUCMinIndex(1);
            catch
                disp('debug')
            end
        end
        % catch
        %     % size(AUCMinIndex)
        %     % disp(shotIdx)
        %     % disp(numShotOn)
        %     fprintf('       There is no laser on data');
        %     fprintf('       Should check repeition rate');
        % end
    end
    shotPairCell{idx}=shotPairBuff;
end
for idx=1:numH5
    h5Struct(idx).shotPair=shotPairCell{idx};
end
fprintf('  Pairing Done\n');
disp('--------------------------------------------------------------')
end
function h5Struct=getDiffAUC(h5Struct,config)
%% Section 1 |
pumpI0Correction=config.pumpI0Correction;
q=config.q;
qIsAUC=config.qIsAUC;
numH5=config.numH5;
noUpsideDown=config.noUpsideDown;
correctionFunc=config.correctionFunc;
pumpI0On=config.pumpI0On;

fprintf('  Obtaining diffs in progress\n');
if pumpI0Correction
    fprintf('  pumpI0Correction On\n');
end
for idx=1:numH5
    laserOn=h5Struct(idx).laserOn;
    eh1rayMXAI_int=h5Struct(idx).eh1rayMXAI_int;
    normfactor=h5Struct(idx).normfactor;
    shotPair=h5Struct(idx).shotPair;
    numShotOn=sum(laserOn);
    I0good = h5Struct(idx).I0good;
    if pumpI0On
        pumpI0=h5Struct(idx).eh1pdcom_channel0p;
    end


    eh1rayMXAI_intOn=eh1rayMXAI_int(:,laserOn);
    eh1rayMXAI_intOff=eh1rayMXAI_int(:,~laserOn&I0good);
    normfactorOn=normfactor(laserOn);
    normfactorOff=normfactor(~laserOn&I0good);

    if pumpI0On
        if noUpsideDown
            pumpI0LaserOn=pumpI0(laserOn); %Normal
        else
            pumpI0LaserOn=pumpI0(~laserOn&I0good); %Upsidedown
        end
    end
try
    diffCurve=eh1rayMXAI_intOn./normfactorOn...
        -eh1rayMXAI_intOff(:,shotPair)./normfactorOff(shotPair);
catch
    numQ = size(eh1rayMXAI_intOff,1);
    diffCurve = zeros(numQ,1);
end
    dataMatrix = abs(diffCurve(qIsAUC,:));
    h5Struct(idx).diffAUC=AUC_Jkim(dataMatrix,q(qIsAUC),[-Inf Inf]);
    h5Struct(idx).diff=diffCurve;
    if pumpI0Correction && pumpI0On
        correctionValue=correctionFunc(pumpI0LaserOn);
        diffCorrection=diffCurve./correctionValue;
        h5Struct(idx).diffCorrection=diffCorrection;
    end
end
fprintf('  Obtaining diffs done\n');
disp('--------------------------------------------------------------');
end
function [cutoffFig,numGood,RMS,I0Fig] = plotCutoff(h5Struct,config)
%% Section 0 | Decapsulation
autoAUCValue=config.autoAUCValue;
numH5=config.numH5;
I0cutoff = config.I0cutoff;
cutoffFig=figureWithoutStealing(22222,'cutoff');
tiledlayout_space_Jkim(5,2);
%% Section 1 | First tile
nexttile([3 1]);
set(gca,'FontWeight','bold');
cutoffAUC=autoAUCValue;
plotCutoffAUC(h5Struct,cutoffAUC,numH5);
% set(gca,'FontWeight','bold');
fontAllocation;
cutoffAUCComment(cutoffAUC);
%% Section 2 | Second tile
nexttile([3 1]);
q = config.q;
chennelROI = config.channelROI;
ROI = config.scatterdROI;

eh1rayMXAI_intNow = h5Struct(end).eh1rayMXAI_int;

I0 = h5Struct(end).(config.channelROI{1}{1});
I = AUC_Jkim(eh1rayMXAI_intNow,q,ROI);

nameChannel = extractBefore(chennelROI{1}{1},'_');

xlableText = sprintf('I_{incident} (%s)',nameChannel);
ylableText = sprintf('I_{scattered}');
% [~,nameWOext,~]=fileparts(h5Struct(end).chiName);
% titleText=sprintf('%s',nameWOext);

[gradi,y0,RMS] = plotI0I(I0,I,xlableText,ylableText,'');
hold on
isGoodI0 = I0 > I0cutoff;

I0_bad = I0(~isGoodI0);
I_bad = I(~isGoodI0);

pl = plot(I0_bad,I_bad,'ro','MarkerFaceColor','r');

hold off

xline(I0cutoff,'--','Cutoff Line','LineWidth',1,'Color','r');

xlim([min(I0) max(I0)]);

fprintf('  RMS    | %.4f\n',RMS);
fprintf('  gradi | %.4g, offset  | %.4g\n',gradi,y0);


disp('--------------------------------------------------------------')

%% Section 3 | Third tile
nexttile([1 2]);
numGood = [];
axisMode = 4; % noaxis
symlogC = 10; % noaxisconstant
RMS = [];
for idx = 1:length(h5Struct)
    goodNow = h5Struct(idx).numShotBadRemove;
    numGood = [numGood goodNow];

    I0 = h5Struct(idx).(chennelROI{1}{1});
    I = h5Struct(idx).AUC;

    p = polyfit(I0,I,1);
    I_theo = polyval(p,I0);
    % R2Now = R2cal(I,I_theo);
    % RMSNow = RRScal(I,I_theo,I0);

    scoreMatirx = corrcoef(I,I0);
    scoreNow = scoreMatirx(1,2);
    scoreNow = -log10(1-scoreNow);
    RMSNow = scoreNow;

    RMS = [RMS RMSNow];
end

plot(numGood,'-ko');
figureAllocation;
fontAllocation;

% xlabel('Position');
ylabel('#Good');

yyaxis right;
[~] = plot(RMS,'-ro');
ylabel('-log(1-\rho)','Interpreter','tex');
%%% Get time %%%

for idx = 1:numH5
    timeChar(idx) = string(h5Struct(idx).date(end-7:end));
end
numMtd = length(timeChar);

if 10 > numH5
    divider = numMtd;
else
    divider = 10;
end

mtdIdx = round(linspace(1,numMtd,divider));

xticks(mtdIdx);
xticklabels(timeChar(mtdIdx));
xlimVal = xlim;
xlim([xlimVal(1) xlimVal(2)+1]);
%%

% switch config.mode
%     case 'TRXL'
%         plot_Jkim_mode(pl,timeChar,axisMode,symlogC);
% end
xlabel('Measured time');

%% Section 4 | Fourth titl
% ylim([min(R2) 1]);

I0Fig = figureWithoutStealing(12123,'I0');
I0All = [];
xtickVal = [];

counter = 1;

for idx = 1:length(h5Struct)
    xtickVal = [xtickVal counter];
    counter = counter + h5Struct(idx).numShot;

    I0 = h5Struct(idx).(chennelROI{1}{1});
    I0All = [I0All I0];
end

plot([1:1:length(I0All)],I0All);
figureAllocation;
fontAllocation;

xticks(xtickVal(mtdIdx));
xticklabels(timeChar(mtdIdx));
xlabel('Time');
ylabel('EH1QBPM');

%
% %%% Get time %%%
% % timeCharNow = cell(1,length(I0All));
% for idx = 1:numH5
%     timeCharNow = string(h5Struct(idx).date(end-7:end));
%     timeChar{idx} = timeCharNow;
%     % if idx == 1
%     %
%     % else
%
% end
% numMtd = length(timeChar);
%
% if 10 > numH5
%     divider = numMtd;
% else
%     divider = 10;
% end
%
% mtdIdx = round(linspace(1,numMtd,divider));
%
% xticks(mtdIdx);
% xticklabels(timeChar(mtdIdx));





%% Section 4


    function cutoffAUCComment(cutoffAUC)
        message1=sprintf('  AUC cutoff Value | %.3f',cutoffAUC);
        fprintf('%s\n',message1);
        disp('--------------------------------------------------------------')
    end
end
function h5Struct=goodSignalCheckChiStruct(h5Struct,config)
cutoffAUC=config.autoAUCValue;
numH5=config.numH5;
shotListBadRemove=nan(1,numH5);
% I0cutoff = config.I0cutoff;
% channelROI = config.channelROI;
isOXCcorrection = config.isOXCcorrection;

for idx=1:numH5
    diffAUC=h5Struct(idx).diffAUC;
    diff=h5Struct(idx).diff;
    laserOn = h5Struct(idx).laserOn;
    % I0good = h5Struct.I0good;
    I0good = h5Struct(idx).I0good;
    if isOXCcorrection
        OXCgood = h5Struct(idx).OXCgood;
    end

    % I0 = h5Struct(idx).(channelROI{1}{1});
    % I0good = I0cutoff<I0;
    % try
    if isOXCcorrection
        try
        h5Struct(idx).goodSignalCheck=(cutoffAUC>=diffAUC) & I0good(laserOn) & OXCgood;
        catch
            disp('debug | 1471')
        end
        onTds = h5Struct(idx).tds(laserOn);
        h5Struct(idx).tdsGood = onTds(h5Struct(idx).goodSignalCheck);
    else
        h5Struct(idx).goodSignalCheck=(cutoffAUC>=diffAUC) & I0good(laserOn);
    end
    % catch
    % disp('');
    % end
    h5Struct(idx).numShotBadRemove=sum(h5Struct(idx).goodSignalCheck);
    h5Struct(idx).diff=diff(:,h5Struct(idx).goodSignalCheck);



    shotListBadRemove(idx)=h5Struct(idx).numShotBadRemove;
end
for idx=1:numH5
    shotListBadRemove(idx)=h5Struct(idx).numShotBadRemove;
end
getBadRemoveCurvesComment(shotListBadRemove);
    function getBadRemoveCurvesComment(shotListBadRemove)
        fprintf('  Bad diff checked\n');
        fprintf('  AvgGoodShotNumber (On) | %d\n\n',round(mean(shotListBadRemove)));
        fprintf('  ----- The Number of GoodShot -----\n');
        fprintf('  %d %d %d %d %d %d %d %d %d %d\n',shotListBadRemove);
        fprintf('\n');
        disp('--------------------------------------------------------------')
    end
end
function h5Struct=getDiffAveSTD(h5Struct,config)
numH5=config.numH5;
for idx=1:numH5
    h5Struct(idx).diffMean=mean(h5Struct(idx).diff,2);
    h5Struct(idx).diffSTDMean=...
        std(h5Struct(idx).diff,0,2)/(sqrt(h5Struct(idx).numShotBadRemove));
end
end
function diffStruct=setDiffStructWithMotor(h5Struct,config)
posi1List=config.posi1List;
posi2List=config.posi2List;
posi3List=config.posi3List;

counter=0;

for idxPosi3=unique(posi3List)
    for idxPosi2=unique(posi2List)
        for idxPosi1=unique(posi1List)
            idxChi=(idxPosi1==posi1List & idxPosi2==posi2List & idxPosi3==posi3List);
            if any(idxChi)
                counter=counter+1;
                diffStruct(counter).position1=idxPosi1;
                diffStruct(counter).position2=idxPosi2;
                diffStruct(counter).position3=idxPosi3;


                diffStruct(counter).motor1=h5Struct(idxChi).motor1;
                diffStruct(counter).motor2=h5Struct(idxChi).motor2;
                diffStruct(counter).motor3=h5Struct(idxChi).motor3;
            else
                fprintf('  Lost h5 data at posi1: %d posi2: %d posi3: %d\n',idxPosi1,idxPosi2,idxPosi3);
            end
        end
    end
end
end
function diffStruct=getFileIndex3DscanDiffStruct(diffStruct,config)
posi1List=config.posi1List;
posi2List=config.posi2List;
posi3List=config.posi3List;
numDiff=length(diffStruct);
for diffIndex=1:numDiff
    position1=diffStruct(diffIndex).position1;
    position2=diffStruct(diffIndex).position2;
    position3=diffStruct(diffIndex).position3;
    diffStruct(diffIndex).file=find(position1==posi1List & position2==posi2List & position3==posi3List);
end
end
function diffStruct=getAveSTDpumpI0(h5Struct,diffStruct,config)
fprintf('  Making DiffStruct in progress\n');
numDiff=length(diffStruct);
pumpI0On=config.pumpI0On;
pumpI0Correction=config.pumpI0Correction;
for idx=1:numDiff
    numShotBadRemoveSum=0;
    diffs=[];
    diffCorrections=[];
    AUCs=[];
    for orderIndex=1:length(diffStruct(idx).file)
        diffNow=h5Struct(diffStruct(idx).file(orderIndex)).diff;
        numShotRemovedNow=h5Struct(diffStruct(idx).file(orderIndex)).numShotBadRemove;
        AUC=h5Struct(diffStruct(idx).file(orderIndex)).AUC;

        diffs=[diffs diffNow]; %#ok<*AGROW>
        AUCs=[AUCs AUC];

        numShotBadRemoveSum=numShotBadRemoveSum+numShotRemovedNow;
        if pumpI0Correction
            diffCorrectionNow=h5Struct(diffStruct(idx).file(orderIndex)).diffCorrection;
            diffCorrections=[diffCorrections diffCorrectionNow];
        end
    end
    numImageNow=numShotBadRemoveSum;

    diffStruct(idx).numImage=numImageNow;
    diffStruct(idx).diffAve=mean(diffs,2);
    diffStruct(idx).diffSTDMean=std(diffs,0,2)/sqrt(numImageNow);
    diffStruct(idx).AUCAvg=mean(AUCs);
    if pumpI0Correction
        diffStruct(idx).diffCorrectionAve=mean(diffCorrections,2);
        diffStruct(idx).diffCorrectionSTDMean=std(diffCorrections,0,2)/sqrt(numImageNow);
    end

end
if pumpI0On
    for idx=1:numDiff
        for orderIndex=1:length(diffStruct(idx).file)
            laserOn=h5Struct(diffStruct(idx).file(orderIndex)).laserOn;
            pumpI0=h5Struct(diffStruct(idx).file(orderIndex)).eh1pdcom_channel0p;
            goodSignalCheck=h5Struct(diffStruct(idx).file(orderIndex)).goodSignalCheck;

            pumpI0LaserOn=pumpI0(laserOn);
            %             pumpI0LaserOn=pumpI0(~laserOn); %Upsidedown

            pumpI0LaserOnGood=pumpI0LaserOn(goodSignalCheck);
            if orderIndex==1
                pumpI0Diff=pumpI0LaserOnGood;
            else
                pumpI0Diff=[pumpI0Diff; pumpI0LaserOnGood]; %#ok<*AGROW>
            end
        end
        diffStruct(idx).eh1pdcom_channel0p=pumpI0Diff;
    end
end
fprintf('  Making DiffStruct done\n');
disp('--------------------------------------------------------------')
end

function diffBinStruct = getdiffBinStruct(h5Struct,diffStruct,config)
timedelay = config.timedelay;
tdOXCcorrect = config.tdOXCcorrect;
binOXCcorrect = config.binOXCcorrect;

isOXCCor = timedelay <= tdOXCcorrect;
TdOXCCor = timedelay(isOXCCor);
TdnoOXCCor = timedelay(~isOXCCor);

%%% Inherit noOXCtd from diff Struct %%%
tdDiff = [];
for idx = 1:length(diffStruct)
    tdDiff = [tdDiff diffStruct(idx).timedelay];
end
TF = ismember(tdDiff,TdnoOXCCor);
diffBinStruct = diffStruct(TF);

%%% Obtain OXC all data %%%
tdAll = [];
% diffAll = [];

for idx = 1:length(h5Struct)
    if ismember(h5Struct(idx).timedelay,TdOXCCor)
        tdsNow = h5Struct(idx).tdsGood;
        % diffNow = h5Struct(idx).diff;

        tdAll = [tdAll tdsNow];
        % diffAll = [diffAll diffNow];
    end
end

bins = makeZeroCenteredBins(tdAll, binOXCcorrect);
M = size(h5Struct(1).diff,1);
acc = initBinAccumulator(bins.nbins, M);

counter = 0;
for idx = 1:length(h5Struct)
    if ismember(h5Struct(idx).timedelay,TdOXCCor)
        % tdsNow = h5Struct(idx).tdsGood;
        for jdx = 1:size(h5Struct(idx).diff,2)
            colVec = h5Struct(idx).diff(:,jdx);
            counter = counter + 1;
            k = bins.binIdx(counter);
            acc = updateBinAccumulator(acc, k, colVec);
        end
    end
end

out = finalizeBins(bins, acc, binOXCcorrect);

tdOXC = out.centers;

for idx = 1:length(tdOXC)
    strNow = struct('timedelay',tdOXC(idx),...
        'file',[],...
        'numImage',out.count(idx),...
        'diffAve',out.mean(:,idx),...
        'diffSTDMean',out.sem(:,idx),...
        'AUCAvg',[]...
        );
    diffBinStruct = [diffBinStruct, strNow];
end

%%% Unique %%%

tdAll = [];
for idx = 1:length(diffBinStruct)
    tdNow = diffBinStruct(idx).timedelay;
    tdAll = [tdAll tdNow];
end
[~,idxUni] = unique(tdAll);
diffBinStruct = diffBinStruct(idxUni);

% disp('debug');


end

function diffStruct = removeBadPointDiffStruct(diffStruct)
numDiff = length(diffStruct);

numImg = [];

for idx = 1:numDiff
    imgNow = diffStruct(idx).numImage;
    numImg = [numImg imgNow];
end

meanVal = mean(numImg,"all");
I = numImg < meanVal*0.1;

% [~,I] = rmoutliers(numImg,'gesd');
% diffStruct = diffStruct(~I);

if any(I)
    fprintf('  Bad timedealy file(s) removed\n\n');
    fprintf('  Removed timedelay idx | %d (ImageNumber %d) \n',[[find(I)]; [numImg(I)]]);
else
    fprintf('  There is no time delay. GREAT!\n');
end
disp('--------------------------------------------------------------')

end


function [timedelay,numTimedelay,diffStruct]=setTimedelayDiffStruct(h5Struct,config)
numH5=config.numH5;
timedelayList=nan(1,numH5);
for idx=1:numH5
    timedelayList(idx)=h5Struct(idx).timedelay;
end
timedelay=unique(timedelayList);
diffStruct=struct;
for timeIndex=1:length(timedelay)
    diffStruct(timeIndex).timedelay=timedelay(timeIndex);
end
numTimedelay=length(timedelay);
end
function diffStruct=getFileIndexDiffStruct(h5Struct,diffStruct,config)
timedelay=config.timedelay;
numH5=config.numH5;
numTimedelay=length(timedelay);

for timeIndex=1:numTimedelay
    counter=0;
    for idx=1:numH5
        if timedelay(timeIndex)==h5Struct(idx).timedelay
            counter=counter+1;
            diffStruct(timeIndex).file(counter)=idx;
        end
    end
end
end
function [diffStruct,numDiffBad] = badDiffStructRemove(diffStructBrok,config)
numDiffBrok=length(diffStructBrok);
pumpI0Correction=config.pumpI0Correction;
mode=config.mode;

%%% Get numCurves %%%%%%%%%%
% I0=[];
% for idx=1:numDiffBrok
%     I0Now = diffStructBrok(idx).eh1qbpm1_totalsum;
%     I0=[I0 AUCAvg];
% end
% [~,badDiffIdx]=rmoutliers(I0,'gesd');

%%% Get goodDiffIdx %%%%%%%%%
for idx=1:numDiffBrok
    switch mode
        case '2Dscan'
            pos1=diffStructBrok(idx).position1;
            pos2=diffStructBrok(idx).position2;
            filename=sprintf('pos2 %d, pos1 %d',pos2,pos1);
        case 'TRXL'
            timedelay=diffStructBrok(idx).timedelay;
            filename=sprintf('diff_av_%.3e',timedelay);
    end

    numCurve=diffStructBrok(idx).numImage;
    AUCAvg=diffStructBrok(idx).AUCAvg;

    if pumpI0Correction
        pumpI0=diffStructBrok(idx).eh1pdcom_channel0p;
    end
    % if badDiffIdx(idx)
    %     fprintf('  %.4g (Bad eh1rayMXAI_int AUC) | %s\n',AUCAvg,filename);
    %     goodDiffStructIdx(idx)=false;
    %     continue;
    % end
    if numCurve<10
        fprintf('  %d (Bad NumImage) | %s\n',numCurve,filename);
        goodDiffStructIdx(idx)=false;
        continue;
    end

    if pumpI0Correction
        if any(isnan(pumpI0))
            fprintf('  No pumpI0 signal | %s\n',filename);
            goodDiffStructIdx(idx)=false;
            continue;
        end
    end
    goodDiffStructIdx(idx)=true;
end
numDiffBad = sum(~goodDiffStructIdx);

%%% Get good Diff %%%%%%%%%%%%%
diffStruct=diffStructBrok(goodDiffStructIdx);


%%% Comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numDiffBad==0
    fprintf('  There is no bad diffAve. G R E A T !!!\n');
else
    fprintf('  Bad diff removed | AUCAvg with (gesd)\n');
    fprintf('  There is(are) %d bad diffAve(s)\n',numDiffBad);
end
disp('--------------------------------------------------------------')
end
function [diffAll,timedelay,numTimedelay]=makeDiffAll(diffStruct)
numDiff=length(diffStruct);
timedelay=[];
diffAll=[];
numImg = [];
for idx=1:numDiff
    diffAve=diffStruct(idx).diffAve;
    timeNow=diffStruct(idx).timedelay;
    % =numImageNow;

    timedelay=[timedelay timeNow];
    diffAll=[diffAll diffAve];
    numImg = [numImg diffStruct(idx).numImage];
end

meanVal = mean(numImg,"all");
I = numImg < meanVal*0.1;
timedelay = timedelay(~I);
diffAll = diffAll(:,~I);

numTimedelay=length(timedelay);
fprintf('  Making DiffAll done\n');
disp('--------------------------------------------------------------')
end
function eh1rayMXAI_intMean=geteh1rayMXAI_intMean(h5Struct)
eh1rayMXAI_intMean=h5Struct(end).eh1rayMXAI_int./h5Struct(end).normfactor;
end
function [diffAll,timedelay,numTimedelay]=nanTreatDiffAll(config)
timedelayBrok=config.timedelay;
diffAllBrok=config.diffAll;

numTimedelay=length(timedelayBrok);

diffAll=[];
timedelay=[];

for idx=1:numTimedelay
    currentDiff=diffAllBrok(:,idx);
    currentTimedelay=timedelayBrok(idx);
    if ~anynan(currentDiff)
        diffAll=[diffAll currentDiff];
        timedelay=[timedelay; currentTimedelay];
    end
end
numTimedelay=length(timedelay);
end
function [U,S,V]=doSVD(diffAll,config)
qIsSVD=config.qIsSVD;
[U,S,V]=svd(diffAll(qIsSVD,:),'econ');
end
function scan2DFig=plot2Dscan(figureNumber,diffStruct,h5Struct,config)
q=config.q;
qIsSVD=config.qIsSVD;
qIsAUC=config.qIsAUC;
numInt=config.numInt;
drawRange=config.drawRange;
plotOption=config.plotOption;
% try

diffStruct=getDiffSum(diffStruct,qIsAUC,q);
qDrawBool=drawRange(1) < q & q < drawRange(2);
qDraw=q(qDrawBool);

numDiff=length(diffStruct);


if numDiff<numInt
    numInt = numDiff;
end

% numH5=length(h5Struct);
numQ=length(q);

Sq=nan(numQ,numDiff);
qSq=nan(numQ,numDiff);
for idx=1:numDiff
    Sq(:,idx)=diffStruct(idx).diffAve;
    qSq(:,idx)=q.*diffStruct(idx).diffAve;
end

scan2DFig=figureWithoutStealing(figureNumber,'2D Scan position Finder');
tiledlayout_space_Jkim(numInt+1,3);

posiList1=nan(1,numDiff);
posiList2=nan(1,numDiff);

motorList1=nan(1,numDiff);
motorList2=nan(1,numDiff);

[U,S,V]=doSVD(Sq,config);
SV=V*S;
for idx=1:numDiff
    diffStruct(idx).SV=SV(idx,:);
end



SV=nan(numInt,numDiff);
sumList=nan(1,numDiff);

for idx=1:numDiff
    posiList1(idx)=diffStruct(idx).position1;
    posiList2(idx)=diffStruct(idx).position2;

    motorList1(idx)=diffStruct(idx).motor1;
    motorList2(idx)=diffStruct(idx).motor2;

    SV(:,idx)=diffStruct(idx).SV(1:numInt);

    sumList(idx)=diffStruct(idx).sum;
end

posi1=unique(posiList1);
posi2=unique(posiList2);

if length(posi1) == 1 || length(posi2) == 1

    fprintf('  Too few data to plot contourf\n');
    return;
end

motor1=unique(motorList1);
motor2=unique(motorList2);

for intIndex=1:numInt+1
    plotMatrix=nan(length(posi1),length(posi2));
    for idxPosi1=posi1
        for idxPosi2=posi2
            idx=(idxPosi1==posiList1 & idxPosi2 == posiList2);
            if any(idx)
                if intIndex==1
                    plotMatrix(idxPosi1,idxPosi2)=diffStruct(idx).sum;
                else
                    plotMatrix(idxPosi1,idxPosi2)=diffStruct(idx).SV(intIndex-1);
                end
            else
                plotMatrix(idxPosi1,idxPosi2)=0;
            end
        end
    end
    % ax(intIndex)=subplot(numInt+1,3,intIndex*3-2:intIndex*3-1);
    ax = nexttile([1 2]);
    contourf(motor1,motor2,plotMatrix',200,'LineColor','None');
    figureAllocation;
    fontAllocation;
    xticks(motor1);
    colormap(ax,bluewhitered);
    yticks(motor2);
    colorbar_Jkim;
    grid on;
    ylabel('motor position 2');
    if intIndex==numInt+1
        xlabel('motor position 1');
    else
        xlabel('');
    end
    if intIndex==1
        % maxValue=max(sumList);
        % maxIndex=find(sumList==maxValue);

        [~,I] = sort(sumList,'descend');
        maxIndex = I(1);

        nameLegend='MAX diffAUC';

        maxValue = sumList(I(1));


    else
        svList = abs(SV(intIndex-1,:));
        [~,I] = sort(svList,'descend');
        maxIndex = I(1);
        maxValue = svList(I);

        % maxIndex=find(abs(SV(intIndex-1,:))==maxValue);

        SValue = S(intIndex-1,intIndex-1);
        nameLegend=sprintf('Component %d',intIndex-1);
    end
            if length(I) < 5
            endIndex = length(I);
        else
            endIndex = 5;
        end
    % set(gca,'FontWeight','bold');
    fontAllocation;

    try
        % for maxIndex = 
    textLegend=sprintf('%s | S | %.3g m1 | (%.6g,%d),m2 | (%.6g,%d)',...
        nameLegend,SValue,...
        motorList1(maxIndex),posiList1(maxIndex),...
        motorList2(maxIndex),posiList2(maxIndex));
    title(textLegend);
        % end
    catch
            textLegend=sprintf('%s | %.3g m1 | (%.6g,%d),m2 | (%.6g,%d)',...
        nameLegend,maxValue,...
        motorList1(maxIndex),posiList1(maxIndex),...
        motorList2(maxIndex),posiList2(maxIndex));
    title(textLegend);
    end

    for maxIndex = I(1:endIndex)
        if intIndex==1
            maxValue = sumList(maxIndex);
        else
            maxValue = svList(maxIndex);
        end

    fontSize = gca().FontSize;
    textComment1=sprintf('    %.3g',maxValue);
    % textComment2=sprintf('m1|%.6g\n',motorList1(maxIndex));
    % textComment3=sprintf('    m2|%.6g',motorList2(maxIndex));
    % text(motorList1(maxIndex),motorList2(maxIndex),...
    %     [textComment1 '\leftarrow ' textComment2 textComment3],...
        text(motorList1(maxIndex),motorList2(maxIndex),...
        ['\leftarrow ' textComment1],...
        'Color','Black','FontWeight','bold','VerticalAlignment','middle','FontSize',fontSize);
    end

    nexttile;
    set(gca,'FontWeight','bold');
    if intIndex==1
        yyaxis right;
        plot(qDraw,mean(h5Struct(1).eh1rayMXAI_int(qDrawBool,:),2),'r','LineWidth',2);
        ax = gca;
        ax.YAxis(1).Exponent = 0;
        ax.YAxis(2).Exponent = 4;
        % ax.YAxis.TickLabelFormat = '%g';
        % ytickformat('%.1e');

        figureAllocation;
        axis tight;
        fontAllocation;
        titleName=('MAX diffAUC');

        ylabel([]);

        %         Ang=char(197);
        %         xlabel(['q (',Ang,'^-^1)']);
        axis tight;
        ylimPre=ylim;
        xlim([drawRange(1),drawRange(2)]);
        ylim([0,ylimPre(2)]);

        yyaxis left
        switch plotOption
            case 'Sq'
                plotData=Sq;
            case 'qSq'
                plotData=qSq;
        end

        plot(qDraw,plotData(qDrawBool,maxIndex),'b','LineWidth',3);
        figureAllocation;
        fontAllocation;
        axis tight;

        %             hline=refline([0 0]);
        %             hline.Color='k';
        yline(0,'-','Color','k','LineWidth',1);
        switch plotOption
            case 'Sq'
                ylabelText='\DeltaS(q)';
            case 'qSq'
                ylabelText='q\DeltaS(q)';
        end
        ylabel(ylabelText);

        axis tight;

    else
        if intIndex==2
            colorCode='k';
        elseif intIndex==3
            colorCode='r';
        elseif intIndex==4
            colorCode='b';
        elseif intIndex==5
            colorCode='m';
        else
            switch mod(intIndex,4)
                case 2
                    colorCode='k';
                case 3
                    colorCode='r';
                case 0
                    colorCode='b';
                case 1
                    colorCode='m';
            end
        end
        titleName=sprintf('U%d',intIndex-1);


        switch plotOption
            case 'Sq'
                plotData=U(:,intIndex-1);
            case 'qSq'
                plotData=U(:,intIndex-1).*q(qIsSVD);
        end

        plot(q(qIsSVD),plotData,colorCode,'LineWidth',2);
        figureAllocation;
        fontAllocation;
        % set(gca,'FontWeight','bold');
        yline(0,'-','Color','k','LineWidth',1);
        xlim([min(q(qIsSVD)),max(q(qIsSVD))]);

    end



    if intIndex==numInt+1
        Ang=char(197);
        xlabel(['q (',Ang,'^-^1)']);
    else
        xlabel('');
    end
    yline(0,'-','Color','k','LineWidth',1);
    switch plotOption
        case 'Sq'
            ylabelText='\DeltaS(q)';
        case 'qSq'
            ylabelText='q\DeltaS(q)';
    end
    ylabel(ylabelText);
    grid on
    title(titleName);
end
% catch
%     fprintf ('       Not enough draw 2D scan. Please wait ...\n');
% end
end
function figureObj = plotSCUCV(figureNumber,config)

U = config.U;
S = config.S;
V = config.V;
numInt = config.numInt;
% numH5 = config.numH5;
numTd = length(config.timedelay);

% numPlotS = numInt;

[CU,CV,Sone,~] = svdPostProcess(U,S,V);


if numTd < numInt
    numInt = numTd;
end
plotComp = 1:numInt;

figureObj=figureWithoutStealing(figureNumber,'SCUCV');


yyaxis left;
plot(plotComp,Sone(plotComp),'Marker','o','Color','k','LineStyle','-','MarkerSize',8,'MarkerFaceColor','k');
ylabel('Singular values','FontWeight','bold','Color','k');
set(gca,'YColor','k');
% ytickformat('auto','FontSize',16);

% yyaxis left
yyaxis right;


plot(plotComp,CU(plotComp),'Marker','diamond','Color','r','LineStyle',':','MarkerSize',8,'MarkerFaceColor','r');
hold on;
plot(plotComp,CV(plotComp),'Marker','square','Color','r','LineStyle','--','MarkerSize',8,'MarkerFaceColor','auto');
hold off;
figureAllocation;
ylabel('Auto correlations','FontWeight','bold','Color','r','Rotation',-90,'VerticalAlignment','bottom');
xlabel('Rank','FontWeight','bold','Color','k');



set(gca,'FontWeight','bold');
% set(gca,'FontSize',10);
fontAllocation;
set(gca,'YColor','r');



grid on;
axis tight;

legend({'S','C(U)','C(V)'},'Location','east','Box','off');

    function [CU,CV,S_one,Ssize]=svdPostProcess(U,S,V)
        [Usize1,Usize2]=size(U);
        [Vsize1,Vsize2]=size(V);
        Ssize=min(Usize1,Vsize1);

        U_auto=zeros(Usize1,Usize2);
        for i=1:Usize1-1
            for j=1:Usize2
                U_auto(i,j)=U(i,j)*U(i+1,j);
            end
        end

        CU=zeros(1,Usize2);
        for j=1:Usize2
            CU(1,j)=sum(U_auto(:,j));
        end

        %%% C(V) Calculation.
        V_auto=zeros(Vsize1,Vsize2);
        for i=1:Vsize2-1
            for j=1:Vsize2
                V_auto(i,j)=V(i,j)*V(i+1,j);
            end
        end

        CV=zeros(1,Vsize2);
        for j=1:Vsize2
            CV(1,j)=sum(V_auto(:,j));
        end

        %%% Rearrange S
        S_one=zeros(1,Ssize);
        for i=1:Ssize
            S_one(1,i)=S(i,i);
        end
    end
end
function [figureObj,RMS] = plotI0IAll(figureNumber,h5Struct,config)

%%% figure
figureObj = figureWithoutStealing(figureNumber,'plotI0I');

cdx = 1;

%%% Decap
q = config.q;
chennelROI = config.channelROI;
ROI = config.scatterdROI;

eh1rayMXAI_intNow = h5Struct(end).eh1rayMXAI_int;

I0 = h5Struct(end).(config.channelROI{cdx}{1});
I = AUC_Jkim(eh1rayMXAI_intNow,q,ROI);

nameChannel = extractBefore(chennelROI{cdx}{1},'_');

xlableText = sprintf('I_{in.} (%s)',nameChannel);
ylableText = sprintf('I_{scatt.} (AUC from a %.1f to %.1f)',ROI(1),ROI(2));
[~,nameWOext,~]=fileparts(h5Struct(end).chiName);
titleText=sprintf('%s',nameWOext);

[gradi,y0,RMS] = plotI0I(I0,I,xlableText,ylableText,titleText);

fprintf('  R2    | %.4f, R2(adj.)| %.4f\n',R2,R2_adj);
fprintf('  gradi | %.4g, offset  | %.4g\n',gradi,y0);
disp('--------------------------------------------------------------')
end
function [gradi,y0,scoreNow] = plotI0I(I0,I,xlabelText,ylabelTxet,titleText)



p = polyfit(I0,I,1);
I_theo = polyval(p,I0);

% [rmsd,nrmsd] = rmsd_Jkim(I_theo,I); % We do not use this rmsd or nrmsd
gradi = p(1);
y0 = p(2);
% R2 = R2cal(I,I_theo);
% R2_adj = R2cal(I,I_theo,2);
% MRMS = RRScal(I,I_theo,I0);

scoreMatirx = corrcoef(I,I0);
scoreNow = scoreMatirx(1,2);
scorelog = -log10(1-scoreNow);

plot(I0,I_theo,'LineStyle','-','LineWidth',2,'Color','r');
hold on
plot(I0,I,'LineStyle','none','Marker','o','Color','k','MarkerSize',1.5);
hold off
% ,'Marker','o','Color','k','MarkerSize',1);
% plObj(1).

figureAllocation;
xlabel(xlabelText,'Interpreter','tex','FontWeight','bold','FontSmoothing','on');
ylabel(ylabelTxet,'Interpreter','tex','FontWeight','bold');
% textComment = sprintf('R^2: %.2f | R^2(adj.) %.2f\nnRMSD: %.2e|RMSD: %.2e\nGradient: %.4g | Offset: %.4g',R2,R2_adj,nrmsd,rmsd,p(1),p(2));
% textComment = sprintf('-log_{10}(nRRS): %.4f\nGradient: %.4g | Offset: %.4g',score,p(1),p(2));
% textComment = sprintf('correlation (\rho): %.4f| -log(1-\rho): %.4f\nGradient: %.4g | Offset: %.4g',scoreNow,scorelog,p(1),p(2));


textComment = {...
    ['\rho (correlation): ' sprintf('%.3f',scoreNow)]...
    ['-log(1-\rho): ' sprintf('%.2f',scorelog)]...
    ['Gradient:' sprintf('%.4g',p(1)) '| Offset: ' sprintf('%.4g',p(2))]...
    };



% textComment = ['correlation (\rho): ' sprintf(%.4f,)]%.4f| -log(1-\rho): %.4f\nGradient: %.4g | Offset: %.4g',scoreNow,scorelog,p(1),p(2)';

fontAllocation;
fontSize = gca().FontSize;

textObj = text_Jkim(textComment,[],[],fontSize);
textObj.Interpreter = 'tex';
% textObj.Fontsize = 12;

% text(0.45,0.2,textComment,'Units','normalized','FontWeight','bold','HorizontalAlignment','right','FontSize',15,'Interpreter','none');
title(titleText,'Interpreter','none');

end
function [fig2d,fig1d]=drawDifferenceScatteringCurves(figureNumber,config)
q=config.q;
run=config.run;
timedelay=config.timedelay;
eh1rayMXAI_intMean=config.eh1rayMXAI_intMean;
qDrawRange=config.drawRange;
timedelayDrawPoint=config.timedelayDrawPoint;
h5header=config.h5header;
plotOption=config.plotOption;
diffAll=config.diffAll;

timedelayLin=timedelayDrawPoint;

timeLin= timedelay < timedelayLin;
% timeLog= 0< timedelay & timedelayLin <= timedelay;
% timeLog= 0< timedelay;
timeLog= 1e-15< timedelay;

qDrawBool=qDrawRange(1) < q & q < qDrawRange(2);
qDraw=q(qDrawBool);

% numTimedelay=length(timedelay);
titleName = titleMaker(config);
% switch targetPosi2
%     case 0
chiFile=sprintf(titleName);
%     otherwise
%         chiFile=sprintf('Rawdata of %s Scan%.3d',titleName,targetPosi2);
% end

numData=size(diffAll,2);
numQ=length(qDraw);

diffPlot=nan(numQ,numData);
switch plotOption
    case 'Sq'
        for idx=1:numData
            diffPlot(:,idx)=diffAll(qDrawBool,idx);
        end
    case 'qSq'
        for idx=1:numData
            diffPlot(:,idx)=q(qDrawBool).*diffAll(qDrawBool,idx);
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Contourf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2d=plot2D(figureNumber,diffPlot,qDraw,timedelay,timeLin,timeLog,titleName,config);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 1D Curves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1d=plot1D(figureNumber+1,diffPlot,q,qDraw,timedelay,chiFile,plotOption,qDrawRange,eh1rayMXAI_intMean,config);
end
function SVFig = plotSVD(figureNumber,config)

fontSize = 8;
%% Section 1 | Decapsulation
run=config.run;
timedelay=config.timedelay;
U=config.U;
S=config.S;
V=config.V;
q=config.q;
qIsSVD=config.qIsSVD;
timedelayDrawPoint=config.timedelayDrawPoint;
numInt=config.numInt;
plotOption=config.plotOption;
numTimedelay=config.numTimedelay;
%% Section 2 | Postprocess
if numTimedelay<=numInt
    numInt=numTimedelay;
end

if numInt ==1
    SVFig = [];
    return;
end

SVFig = figureWithoutStealing(figureNumber,'plotSVD',[1800 1000]);
try
% if isfield(config,'SVFig')
%     SVFig=config.SVFig;
%     set(0,'CurrentFigure',SVFig);
% else
% SVFig=figure(figureNumber);
% SVFig.NumberTitle='Off';
% SVFig.Name='SVD';
% SVFig.Color=[1 1 1];
% end
set(gca,'FontWeight','bold');

timeLin = timedelay <= timedelayDrawPoint;
timeLog = 1e-15< timedelay;

numTimeLog=length(timedelay(timeLog));
tiledlayout_space_Jkim(numInt,4);

%%%% subplot for V %%%%%%%%%%%%%%%%%%%%%%%%%
for idx=1:numInt
    nexttile;
    [~,colorName]=ordinalGet(idx);
    colorCode1=sprintf('%s',colorName);
    colorCode2=sprintf('-%so',colorName);
    if any(timeLog)
        % if ~any(timeLin)
        %     % subplot(numInt,4,4*idx-3:4*idx-1);
        %
        %     %             set(gca,'FontWeight','bold');
        % else
        %     % subplot(numInt,4,4*idx-1);
        % end


        %     ordinalNumber='1st';
        %     colorName='k';

        if numTimeLog>10
            % nexttile;
            plot(timedelay(timeLog),V(timeLog,idx),colorCode1);
            set(gca, 'XScale', 'Log');
        else
            plot(V(timeLog,idx),colorCode2);
        end
        figureAllocation;
        fontAllocation(fontSize);
        ylabel([]);

        % title(titleName);
        % text_Jkim(titleName);
        % text(0.5,0.94,titleName,'Units','normalized','HorizontalAlignment','center');
        comment = sprintf('V%d',idx);
        if config.isOXCcorrection
            colorNow = 'r';
            fntSize = 14;

        else
            colorNow = 'k';
            fntSize = 12;
            % text_Tiled(comment,13,colorNow);
        end
        text_Tiled(comment,fntSize,colorNow);

        % text()

        % set(gca,'FontWeight','bold');

        % grid on;

        if numTimeLog>10
            tickValue=timepointTickMaker(timedelay(timeLog));
            xticks(tickValue);
            xlim([min(timedelay(timeLog)),max(timedelay(timeLog))]);
        else
            tickValue=timedelay(timeLog);
            xticks(1:length(timedelay(timeLog)));
        end
        tickLabels=timepointTickLabelMaker(tickValue);
        xticklabels(tickLabels);


        ylim([min(V(:,idx)),max(V(:,idx))]);
        axis tight;
        ax = gca;
        ax.XMinorTick = 'off';
        ax.XMinorGrid = 'off';
        % set(gca,'XMinorTick','off');
        % set(gca,'XTick','off');
        yliner;
        if idx == numInt
            xlabel('timedelay');
        else
            xlabel([]);
        end
    end
    %%% Linear Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(timeLin)
        nexttile([1,2]);
        % if ~any(timeLog)
        %     subplot(numInt,4,4*idx-3:4*idx-1);
        %     %             set(gca,'FontWeight','bold');
        % else
        %     subplot(numInt,4,4*idx-3:4*idx-2);
        %     %             set(gca,'FontWeight','bold');
        % end



        plot(timedelay(timeLin)/1e-12,V(timeLin,idx),colorCode2);
        figureAllocation;
        fontAllocation(fontSize);
        % title(titleName);
        text_Tiled(comment,fntSize,colorNow);
        axis tight;
        set(gca,'FontWeight','bold');

        xlim([min(timedelay(timeLin))/1e-12,timedelayDrawPoint/1e-12]);
        ylim([min(V(:,idx)),max(V(:,idx))]);
        xlabel('timedelay (ps)');

        % yticklabels([]);
        ylabel([]);
        yliner;
        grid on;
        if idx ~= numInt
            xlabel([]);
        end
    end

    %%% lSV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subplot(numInt,4,4*idx);
    nexttile;
    switch plotOption
        case 'Sq'
            pltObj = plot(q(qIsSVD),U(:,idx),colorName);
        case 'qSq'
            pltObj = plot(q(qIsSVD),q(qIsSVD).*U(:,idx),colorName);
    end
    pltObj.LineWidth = 2;
    pltObj.AlignVertexCenters = 'on';

    figureAllocation;
    fontAllocation(fontSize);
    comment = sprintf('U%d',idx);
    % title(titleName);
    % text_Tiled(titleName);
    text_Tiled(comment,fntSize,colorNow);

    % ylabel('q\DeltaS(q)');
    Ang=char(197);
    switch plotOption
        case 'Sq'
            xlabel(['(',Ang,'^-^1)']);
        case 'qSq'
            xlabel(['q (',Ang,'^-^1)']);
    end
    axis tight;
    grid on;
    set(gca,'FontWeight','bold');
    ylabel([]);
    % yticklabels([]);

    yliner;

    if idx ~= numInt
        xlabel([]);
    end
end
catch
end
    function yliner
        yline(0,'--','Color','k','LineWidth',1.5);
    end
end
function pumpI0Fig=plotPumpI0(figureNumber,diffStruct,config)

% version = 'v2';
version = 'v3';

switch version
    case 'v2'

        V=config.V;
        correctionFunc=config.correctionFunc;

        numInt=2;

        pumpI0Fig=figureWithoutStealing(figureNumber,'pumpI0Comparsion');

        % if isfield(config,'pumpI0Fig')
        %     pumpI0Fig=config.pumpI0Fig;
        %     set(0,'CurrentFigure',pumpI0Fig);
        % else
        % pumpI0Fig=figure(figureNumber);
        % pumpI0Fig.NumberTitle='Off';
        % pumpI0Fig.Name='pumpI0 Comparison';
        % pumpI0Fig.Color=[1 1 1];
        % end

        numDiff=length(diffStruct);
        pumpI0List=nan(numDiff,1);
        timepoint=nan(numDiff,1);
        for idx=1:numDiff
            timepoint(idx)=diffStruct(idx).timedelay;
            pumpI0List(idx)=mean(diffStruct(idx).eh1pdcom_channel0p);
        end

        for idx=1:numInt+2
            subplot(numInt+2,1,idx);
            if idx==1
                plot(pumpI0List,'-ko','LineWidth',2);
                title('pumpI0Avg (ylim0)');
                axis tight;
                ylim([0 inf]);


            elseif idx==2
                plot(pumpI0List,'-ko','LineWidth',2);
                title('pumpI0Avg');
                axis tight;
            else
                [ordinalNumber,~]=ordinalGet(idx-2);
                correctionRSV=V(:,idx-2)./correctionFunc(pumpI0List);
                correctionRSV=correctionRSV/norm(correctionRSV);
                plot(V(:,idx-2),'k','LineWidth',1.5);
                hold on
                plot(correctionRSV,'r','LineWidth',1.5);
                hold off

                titleName=sprintf('%s rSV',ordinalNumber);
                title(titleName);
                axis tight;

                legend('Raw rSV','corrected with Avg pumpI0');
            end

            grid on;
            samplingNumber = 5;
            % tickValue=timepoint(1:samplingNumber:length(timepoint));
            tickValue=timepoint(1:samplingNumber:length(timepoint));
            xticks(1:samplingNumber:length(timepoint));
            tickLabels=timepointTickLabelMaker(tickValue);
            xticklabels(tickLabels);
            % xlabel('timepoint');
            set(gca,'FontWeight','bold');
        end
    case 'v3'



        V=config.V;
        correctionFunc=config.correctionFunc;

        if config.numH5 ==1
            return;
        else
            numInt=2;
        end



        %%% data calculation %%%
        numDiff=length(diffStruct);
        pumpI0List=nan(numDiff,1);
        timepoint=nan(numDiff,1);
        for idx=1:numDiff
            timepoint(idx)=diffStruct(idx).timedelay;
            pumpI0List(idx)=mean(diffStruct(idx).eh1pdcom_channel0p);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%

        pumpI0Fig=figureWithoutStealing(figureNumber,'pumpI0Comparsion');

        % if isfield(config,'pumpI0Fig')
        %     pumpI0Fig=config.pumpI0Fig;
        %     set(0,'CurrentFigure',pumpI0Fig);
        % else
        % pumpI0Fig=figure(figureNumber);
        % pumpI0Fig.NumberTitle='Off';
        % pumpI0Fig.Name='pumpI0 Comparison';
        % pumpI0Fig.Color=[1 1 1];
        % end


        tiledlayout_Jkim(4,1);

        axisMode = 4; % symlog
        symlogC = 10; % symlogConstant

        for idx=1:numInt+2
            % subplot(numInt+2,1,idx);
            nexttile;
            if idx==1
                pltObj = plot(pumpI0List,'-ko','LineWidth',2);
                plot_Jkim_mode(pltObj,timepoint,axisMode,symlogC);
                figureAllocation;
                ylabel([]);

                % title('pumpI0Avg (ylim0)','FontSize',16);

                titleName = 'pumpI0Avg (ylim0)';
                axis tight;
                ylim([0 inf]);


            elseif idx==2
                pltObj = plot(pumpI0List,'-ko','LineWidth',2);
                plot_Jkim_mode(pltObj,timepoint,axisMode,symlogC);
                figureAllocation;


                titleName = 'pumpI0Avg';
                ylabel([]);
                xlabel([]);
            else
                correctionRSV=V(:,idx-2)./correctionFunc(pumpI0List);
                correctionRSV=correctionRSV/norm(correctionRSV);
                plot(V(:,idx-2),'k','LineWidth',1.5);
                hold on
                pltObj = plot(correctionRSV,'r','LineWidth',1.5);
                hold off
                plot_Jkim_mode(pltObj,timepoint,axisMode,symlogC);
                figureAllocation;
                ylabel([]);
                xlabel('Timedelay');

                titleName=sprintf('V%d',idx-2);
                % text_Tiled(titleName);
                % axis tight;

                legend('Raw rSV','corrected with Avg pumpI0','Location','best');
                legend('Box','off');
            end
            text_Tiled(titleName);
            % grid on;
            % samplingNumber = 5;
            % % tickValue=timepoint(1:samplingNumber:length(timepoint));
            % tickValue=timepoint(1:samplingNumber:length(timepoint));
            % xticks(1:samplingNumber:length(timepoint));
            % if idx == numInt+2
            %     tickLabels=timepointTickLabelMaker(tickValue);
            %     xticklabels(tickLabels);
            % end
            % % xlabel('timepoint');
            % set(gca,'FontWeight','bold');
        end
end
end
% function drawLastScanSV(figureNumber,h5Struct,config)
% %         if plotEachScan
% %             switch targetScan
% %                 case 0
% %                     targetPosi2=max(posi2List);
% %                 otherwise
% %                     targetPosi2=targetScan;
% %             end
% %             timedelay=unique(timedelayList(targetPosi2==posi2List));
% %             numQ=length(q);
% %             numTimedelay=length(timedelay);
% %             diffAll=nan(numQ,numTimedelay);
% %             numChi=length(chiStruct);
% %             for timeIndex=1:numTimedelay
% %                 for fileIndex=1:numChi
% %                     if timedelay(timeIndex)==chiStruct(fileIndex).timedelay && targetPosi2==chiStruct(fileIndex).position2
% %                         diffAll(:,timeIndex)=chiStruct(fileIndex).diffMean;
% %                     end
% %                 end
% %             end
% %
% %             qBool= qDrawStart < q & q < qDrawPoint;
% %
% %             counter=0;
% %
% %             for timeCompIndex=1:length(comparisonTimedelay)
% %                 if timeCompIndex==1
% %                     compFig=figure(666);
% %                     compFig.NumberTitle='Off';
% %                     compFig.Name='Comparison';
% %                 end
% %                 targetIndex=find(comparisonTimedelay(timeCompIndex)==timedelayList);
% %
% %                 for idx=targetIndex
% %                     counter=counter+1;
% %                     if counter==1
% %                         diffComp=q(qBool).*chiStruct(idx).diffMean(qBool);
% %
% %                     else
% %                         diffComp=[diffComp, q(qBool).*chiStruct(idx).diffMean(qBool)];
% %                     end
% %                     text=sprintf('timedelay %.3g, Scan %.3d',chiStruct(idx).timedelay,chiStruct(idx).position2);
% %                     legendComment(counter)=string(text);
% %                     %                     plot(q(qBool),q(qBool).*chiStruct(idx).diffMean,'LineWidth',3);
% %                 end
% %                 if counter==1
% %                     plot(q(qBool),diffComp(:,end),'LineWidth',3);
% %                     legend(legendComment);
% %                     grid on;
% %                 elseif counter >1
% %                     plot(q(qBool),diffComp(:,1:end-1),'LineWidth',1);
% %                     hold on
% %                     plot(q(qBool),diffComp(:,end),'LineWidth',3);
% %                     legend(legendComment);
% %                     grid on;
% %                     hold off
% %                 else
% %                 end
% %             end
% %             try
% %                 switch plotOption
% %                     case 'Sq'
% %                         [U,S,V]=doSVD(diffAll,config);
% %                     case 'qSq'
% %                         diffAlltarget=nan(size(diffAll,1),size(diffAll,2));
% %                         for idx=1:size(diffAll,2)
% %                             diffAlltarget(:,idx)=q.*diffAll(:,idx);
% %                         end
% %                         [U,S,V]=doSVD(diffAlltarget,qIsSVD);
% %                 end
% %                 try
% %                     drawSVDAll(6666,runList,timedelay,U,S,V,q,qIsSVD,refTimedelay,timedelayDrawPoint,numInt,plotOption);
% %                     plotDoneComment;
% %                 catch
% %                     fprintf('  Not enough the number of component to draw SVD figure\n');
% %                 end
% %             catch
% %                 noSVDComment;
% %             end
% %
% %             [~,timedelay]=setTimedelayDiffStruct(timedelayList);
% %         end
% end
function plotStartComment
message1=sprintf('  Plotting in progress');
fprintf('%s\n',message1);
end
function plotDoneComment
message1=sprintf('  Plotting done');
fprintf('%s\n%s\n',message1);
disp('--------------------------------------------------------------')
end
function [folderResults,folderDiffAve,folderDiffAveCorrection,folderEachScan,folderDiffs,folderValue,folderImg,folderSVD]=outputFolderAllocator(config)
parentFolderResults=config.parentFolderResults;
h5header=config.h5header;
pumpI0Correction=config.pumpI0Correction;

figureSave=config.figureSave;
eachScanSave=config.eachScanSave;
diffsSave=config.diffsSave;

folderResultsHeader=sprintf('%s/%s',parentFolderResults,h5header);
folderResults=sprintf('%s/run%.3d',folderResultsHeader,config.run);
folderDiffAve=sprintf('%s/DiffAve',folderResults);
folderDiffAveCorrection=sprintf('%s/DiffAveCorrection',folderResults);
folderEachScan=sprintf('%s/EachScan',folderResults);
folderDiffs=sprintf('%s/Diffs',folderResults);
folderValue=sprintf('%s/Values',folderResults);
folderImg=sprintf('%s/Img',folderResults);
folderSVD=sprintf('%s/SVD',folderResults);

[~,~]=mkdir(folderResultsHeader);
[~,~]=mkdir(folderResults);

[~,~]=mkdir(folderDiffAve);
[~,~]=mkdir(folderValue);
[~,~]=mkdir(folderSVD);

if pumpI0Correction
    [~,~]=mkdir(folderDiffAveCorrection);
end

if eachScanSave
    [~,~]=mkdir(folderEachScan);
end
if diffsSave
    [~,~]=mkdir(folderDiffs);
end
if figureSave
    [~,~]=mkdir(folderImg);
end
end
function writeOutputDiffav(diffStruct,folderName,config)
% folderDiffAve=config.folderDiffAve;

folderResults = config.folderResults;
q=config.q;

folderNow = sprintf('%s/%s',folderResults,folderName);
[~,~] = mkdir(folderNow);

fprintf('  Start to write diffAve. Do not control + C\n');

numDiff=length(diffStruct);
for idx=1:numDiff
    diffAve=diffStruct(idx).diffAve;
    diffSTDMean=diffStruct(idx).diffSTDMean;
    qDiffAve=diffAve.*q;

    name=sprintf('%s/diff_av_%.6e.txt',folderNow,diffStruct(idx).timedelay);
    fid=fopen(name,'w');
    fprintf(fid, 'number_of_images: %d\n',diffStruct(idx).numImage);
    fclose(fid);
    dataWrite=[q diffAve diffSTDMean qDiffAve];
    writematrix(dataWrite,name,'Delimiter','tab','WriteMode','append');
    fcloseWithComment;
end
fcloseWithComment;
writeOutputComment(folderNow);
    function writeOutputComment(folderChi)
        message1=sprintf('  DiffavSaving done');
        message2=sprintf('  cd %s',folderChi);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end



function writeOutputDiffavMatrix(diffStruct,folderName,config)

folderResults=config.folderResults;

folderNow = sprintf('%s/%s',folderResults,folderName);
[~,~] = mkdir(folderNow);

q=config.q;

data = [];

std = [];

numImg = [];

td = [];

nameDat = sprintf('%s/azi_dat.dat',folderNow);
nameStd = sprintf('%s/azi_std.dat',folderNow);
nameshots = sprintf('%s/nshots.dat',folderNow);

numDiff=length(diffStruct);

for index=1:numDiff

    % Allocation
    dataNow=diffStruct(index).diffAve;
    stdMeanNow=diffStruct(index).diffSTDMean;
    numImgNow = diffStruct(index).numImage;
    tdNow = diffStruct(index).timedelay;
    data = [data dataNow];
    std = [std stdMeanNow];
    numImg = [numImg; numImgNow];
    td = [td; tdNow];
    % qDiffAve=diffAve.*q;

    % name=sprintf('%s/diff_av_%.3e.txt',folderDiffAve,diffStruct(index).timedelay);

    % fid=fopen(name,'w');

    % fprintf(fid, 'number_of_images: %d\n',diffStruct(index).numImage);

    % fclose(fid);

end
dataWrite = [[0; q],[ td'; data]];
writematrix(dataWrite,nameDat,'Delimiter','comma','WriteMode','overwrite');
fcloseWithComment;
dataWrite = [[0; q],[ td'; std]];
writematrix(dataWrite,nameStd,'Delimiter','comma','WriteMode','overwrite');
fcloseWithComment;
dataWrite = [td, numImg];
writematrix(dataWrite,nameshots,'Delimiter','comma','WriteMode','overwrite');
fcloseWithComment;
writeOutputComment(folderNow);

    function writeOutputComment(folderChi)
        message1=sprintf('  DiffavSaving done');
        message2=sprintf('  cd %s',folderChi);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end






function writeStaticMean(h5Struct,config)
folderResults=config.folderResults;

q=config.q;

fprintf('  Start to StaticMean\n');

% numH5=length(h5Struct);
staticAll = [];
idx = 1;
% for idx=1:numH5
    I0good = h5Struct(idx).I0good;
    staticNow = h5Struct(idx).eh1rayMXAI_int;
    normfactorNow = h5Struct(idx).normfactor;

    staticNormNow = staticNow(:,I0good)./normfactorNow(:,I0good);

    staticAll = [staticAll staticNormNow];
% end
staticMean = mean(staticAll,2);
staticStd = std(staticAll,0,2);

namePath=sprintf('%s/%s',folderResults,'staticMean.txt');
writematrix([q staticMean staticStd],namePath,'Delimiter','\t');


fcloseWithComment;
writeOutputComment(folderResults);
    function writeOutputComment(folderChi)
        message1=sprintf('  StaticMean done');
        message2=sprintf('  cd %s',folderChi);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end






function writeOutputDiffavCorrection(diffStruct,config)
folderDiffAveCorrection=config.folderDiffAveCorrection;
q=config.q;

numDiff=length(diffStruct);
for idx=1:numDiff
    diffAve=diffStruct(idx).diffCorrectionAve;
    diffSTDMean=diffStruct(idx).diffCorrectionSTDMean;
    qDiffAve=diffAve.*q;

    name=sprintf('%s/diff_av_%.3e.txt',folderDiffAveCorrection,diffStruct(idx).timedelay);
    fid=fopen(name,'w');
    fprintf(fid, 'number_of_images: %d\n',diffStruct(idx).numImage);
    fclose(fid);
    dataWrite=[q diffAve diffSTDMean qDiffAve];
    writematrix(dataWrite,name,'Delimiter','tab','WriteMode','append');
    fcloseWithComment;
end
fcloseWithComment;
writeOutputComment(folderDiffAveCorrection);
    function writeOutputComment(folderChi)
        message1=sprintf('  DiffavCorrectionSaving done');
        message2=sprintf('  cd %s',folderChi);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end


function writeValue(h5Struct,config)
% tic

fprintf('  Start write value matrix. Do not control + C\n')

folderValue=config.folderValue;
numH5=config.numH5;
isOXCcorrection = config.isOXCcorrection;

for idx=1:numH5
    goodSignalCheck=h5Struct(idx).goodSignalCheck;
    shotPair=h5Struct(idx).shotPair;
    laserOn=h5Struct(idx).laserOn;
    normfactor=h5Struct(idx).normfactor;
    name=h5Struct(idx).chiName;
    timestampPIDString=h5Struct(idx).timestampPIDString;
    

    normfactorOn=normfactor(laserOn);
    normfactorOff=normfactor(~laserOn);
    timestampOn=timestampPIDString(laserOn);
    timestampOff=timestampPIDString(~laserOn);
    numShot=sum(laserOn);

    normfactorPair=normfactorOff(shotPair);
    timestampPair=timestampOff(shotPair);

    namePath=sprintf('%s/%s_%s',folderValue,'valueMatrix',name);

    % writematrix([goodSignalCheck' shotPair' normfactorOn' normfactorPair' timestampOn' timestampPair'],namePath,'Delimiter','tab');

    fileID=fopen(namePath,'w');
    if isOXCcorrection
        tds = h5Struct(idx).tds;
        for shotIdx=1:numShot
            fprintf(fileID,'%d \t %d \t %g \t %g \t %s \t %s \t %g \n',...
                goodSignalCheck(shotIdx),...
                shotPair(shotIdx),...
                normfactorOn(shotIdx),...
                normfactorPair(shotIdx),...
                timestampOn(shotIdx),...
                timestampPair(shotIdx),...
                tds(shotIdx)...
                );
        end
    else
        for shotIdx=1:numShot
            fprintf(fileID,'%d \t %d \t %g \t %g \t %s \t %s\n',...
                goodSignalCheck(shotIdx),...
                shotPair(shotIdx),...
                normfactorOn(shotIdx),...
                normfactorPair(shotIdx),...
                timestampOn(shotIdx),...
                timestampPair(shotIdx)...
                );
        end
    end
    fclose(fileID);
    fcloseWithComment;
end
fcloseWithComment;
writeValueComment(folderValue);
% toc
    function normfactorNeg=makeNormfactorNeg(chiStruct,refFile,fileIndex,shotPair)
        normfactorNeg=nan(1,chiStruct(fileIndex).numShot);
        for chiIndex=1:chiStruct(fileIndex).numShot
            normfactorNeg(chiIndex)=chiStruct(refFile).normfactor(shotPair(chiIndex));
        end
    end
    function writeValueComment(folderValue)
        fprintf('  Writing ValueMatrix (Timestamp, NormFactor, BadCheck)\n')
        fprintf('  cd %s\n',folderValue);
        disp('--------------------------------------------------------------')
    end
end
function writeRsqaure(config)
R2 = config.R2;
numGood = config.numGood;
switch config.mode
    case 'TRXL'
        td = config.td_org;
    otherwise
        td = 1:length(R2);
end
folderResults=config.folderResults;

namePath=sprintf('%s/%s',folderResults,'R2andNumGood.txt');
writematrix([td' R2' numGood'],namePath);

    function writeValueComment(folderValue)
        fprintf('  Writing Rsquare\n')
        fprintf('  cd %s\n',folderValue);
        disp('--------------------------------------------------------------')
    end
end
function writeSpecification(config)
folderResults=config.folderResults;
normRange=config.normRange;
AUCRange=config.AUCRange;
SVDRange=config.SVDRange;
drawRange=config.drawRange;
cutoffAUC=config.cutoffAUC;
XrayCenter=config.XrayCenter;
repetitionRate=config.repetitionRate;
pumpI0Correction=config.pumpI0Correction;
correctionFunc=config.correctionFunc;
noUpsideDown=config.noUpsideDown;

namePath=sprintf('%s/%s',folderResults,'Specification.txt');

fileID=fopen(namePath,'w');
fprintf(fileID,'normRange:%.3g,%.3g\n',normRange(1),normRange(2));
fprintf(fileID,'AUCRange:%.3g,%.3g\n',AUCRange(1),AUCRange(2));
fprintf(fileID,'SVDRange:%.3g,%.3g\n',SVDRange(1),SVDRange(2));
fprintf(fileID,'drawRange:%.3g,%.3g\n',drawRange(1),drawRange(2));
fprintf(fileID,'cutoffAUC:%.6g\n',cutoffAUC);
fprintf(fileID,'XrayCenter:%.6g\n',XrayCenter);
fprintf(fileID,'repetitionRate:%d\n',repetitionRate);
fprintf(fileID,'pumpI0Correction:%d\n',pumpI0Correction);
fprintf(fileID,'correctionFunc:%s\n',func2str(correctionFunc));
fprintf(fileID,'noUpsideDown:%d\n',noUpsideDown);
fprintf(fileID,'Analyzer:%s\n',mfilename);
fclose(fileID);
fcloseWithComment;
writeValueComment(folderResults);
    function writeValueComment(folderResults)
        fprintf('  Specification files writing done\n')
        fprintf('  cd %s\n',folderResults);
        disp('--------------------------------------------------------------')
    end
end
function writeSVDresults(config)
folderSVD=config.folderSVD;
timedelay=config.timedelay;
q=config.q;
diffAll=config.diffAll;
qIsSVD=config.qIsSVD;
U=config.U;
S=config.S;
V=config.V;


fileName=sprintf('%s/rawdata.dat',folderSVD);
fid=fopen(fileName,'w');
fprintf(fid,'%-9s','Time');
fprintf(fid,'\t%-17.5g',timedelay);
fprintf(fid,'\n');
numQ=length(q);
numTimedelayUni=length(timedelay);
for qIndex=1:numQ
    fprintf(fid,'%-9.6g',q(qIndex));
    for timeIndex=1:numTimedelayUni
        fprintf(fid, '\t%-17.10g',diffAll(qIndex,timeIndex));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%% Making Autocorrelation
qSvd=q(qIsSVD);
[Usize1,Usize2]=size(U);
[Vsize1,Vsize2]=size(V);
Ssize=min(Usize1,Vsize1);
UAutoCor=zeros(Usize1,Usize2);
for fileIndex=1:Usize1-1
    for shotIndex=1:Usize2
        UAutoCor(fileIndex,shotIndex)=U(fileIndex,shotIndex)*U(fileIndex+1,shotIndex);
    end
end
CU=zeros(1,Usize2);
for shotIndex=1:Usize2
    CU(shotIndex)=sum(UAutoCor(:,shotIndex));
end
VAutoCor=zeros(Vsize1,Vsize2);
for fileIndex=1:Vsize1-1
    for shotIndex=1:Vsize2
        VAutoCor(fileIndex,shotIndex)=V(fileIndex,shotIndex)*V(fileIndex+1,shotIndex);
    end
end
CV=zeros(1,Vsize2);
for shotIndex=1:Vsize2
    CV(shotIndex)=sum(VAutoCor(:,shotIndex));
end
SVector=zeros(1,Ssize);
for fileIndex=1:Ssize
    SVector(fileIndex)=S(fileIndex,fileIndex);
end
SV=V*S;
%%% Writing lSV
fileName=sprintf('%s/lSV.dat',folderSVD);
fid=fopen(fileName,'w');
fprintf(fid,'%-9s','q');
for fileIndex=1:Usize2
    titleName=sprintf('U%d',fileIndex);
    fprintf(fid,'\t%-17s',titleName);
end
fprintf(fid,'\n');
for fileIndex=1:Usize1
    fprintf(fid,'%-9.6g',qSvd(fileIndex));
    for shotIndex=1:Usize2
        fprintf(fid,'\t%-17.10g',U(fileIndex,shotIndex));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%% Writing rSV
fileName=sprintf('%s/rSV.dat',folderSVD);
fid=fopen(fileName,'w');
fprintf(fid,'%-9s','Time');
for fileIndex=1:Vsize1
    titleName=sprintf('V%d',fileIndex);
    fprintf(fid,'\t%-17s',titleName);
end
fprintf(fid,'\n');
for fileIndex=1:Vsize1
    fprintf(fid,'%-9.6g',timedelay(fileIndex));
    for shotIndex=1:Vsize2
        fprintf(fid,'\t%-17.10g',V(fileIndex,shotIndex));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%%% Writing SV
fileName=sprintf('%s/SV.dat',folderSVD);
fid=fopen(fileName,'w');
fprintf(fid,'%-9s','Time');
for fileIndex=1:Vsize1
    titleName=sprintf('SV%d',fileIndex);
    fprintf(fid,'\t%-17s',titleName);
end
fprintf(fid,'\n');
for fileIndex=1:Vsize1
    fprintf(fid,'%-9.6g',timedelay(fileIndex));
    for shotIndex=1:Vsize2
        fprintf(fid,'\t%-17.10g',SV(fileIndex,shotIndex));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%%% Writing S&C(U)&C(V)&LRA
rank=(1:1:Ssize);
fileName=sprintf('%s/S_Auto_LRA.dat',folderSVD);
fid=fopen(fileName,'w');
fprintf(fid,'%-3s \t %-10s \t %-17s \t %-17s \t %-17s \n','Rank','S','C(U)','C(V)');
for fileIndex=1:Ssize
    fprintf(fid,'%-3d \t %-10.6g \t  %-17.10g %-17.10g \n',rank(fileIndex),SVector(fileIndex),CU(fileIndex),CV(fileIndex));
end
fclose(fid);
fcloseWithComment;
writeSVDresultsComment(folderSVD);
    function writeSVDresultsComment(folderSVD)
        message1=sprintf('  Making rawdata and SVD results done');
        message2=sprintf('  cd %s',folderSVD);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end
function writeEachScan(h5Struct,config)
numH5=config.numH5;
folderEachScan=config.folderEachScan;
q=config.q;
for idx=1:numH5
    folderScan=sprintf('%s/scan%.3d',folderEachScan,h5Struct(idx).scan);
    [~,~]=mkdir(folderScan);
    resultsFile=sprintf('%s/diff_av_%.3e',folderScan,h5Struct(idx).timedelay);

    diffAve=h5Struct(idx).diffMean;
    diffSTDMean=h5Struct(idx).diffSTDMean;
    qDiffAve=diffAve.*q;
    numCurve=h5Struct(idx).numShotBadRemove;
    if any(diffAve)
        break
    end

    fid=fopen(resultsFile,'w');
    fprintf(fid, 'number_of_images: %d\n',numCurve);
    fclose(fid);
    dataWrite=[q diffAve diffSTDMean qDiffAve];
    writematrix(dataWrite,resultsFile,'Delimiter','tab','WriteMode','append');
    fcloseWithComment;
end
writeEachScanComment(folderScan);
    function writeEachScanComment(folderScanProcess)
        message1=sprintf('  Making each scan diff file done');
        message2=sprintf('  cd %s',folderScanProcess);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end
function writeDiffs(h5Struct,config)
folderDiffs=config.folderDiffs;
q=config.q;
numH5=config.numH5;
pumpI0On=config.pumpI0On;

for idx=1:numH5
    timedelay=h5Struct(idx).timedelay;
    if pumpI0On
        pumpI0=h5Struct(idx).eh1pdcom_channel0p;
    end
    diff=h5Struct(idx).diff;
    name=sprintf('%s/diffs_%.3e.txt',folderDiffs,timedelay);
    if pumpI0On && any(pumpI0)
        dataWrite=[[0 q']' [pumpI0'; diff]];
    else
        dataWrite=[q diff];
    end
    writematrix(dataWrite,name,'Delimiter','tab');
end
fcloseWithComment;
writeDiffsComment(folderDiffs);
    function writeEachScanComment(folderDiffs)
        message1=sprintf('  Making each scan diff file done');
        message2=sprintf('  cd %s',folderDiffs);
        fprintf('%s\n%s\n',message1,message2);
        disp('--------------------------------------------------------------')
    end
end
% function writeChiFiles(chiStruct,folderResults,q)
% numChi=length(chiStruct);
% numQ=length(q);writeImgI
% for chiIndex=1:numChi
%     name=sprintf('%s/%s',folderResults,chiStruct(chiIndex).name);
%     fID=fopen(name,'w');
%     fprintf(fID,'%s','timestampPulseID');
%     for shotIndex=1:chiStruct(chiIndex).numShot
%         fprintf(fID,'\t%s',chiStruct(chiIndex).timestamp(shotIndex));
%     end
%     fprintf(fID,'\n');
%
%     for qIndex=1:numQ
%         fprintf(fID,'%.8g',q(qIndex));
%         for shotIndex=1:chiStruct(chiIndex).numShot
%             fprintf(fID,'\t%.8g',chiStruct(chiIndex).eh1rayMXAI_int(qIndex,shotIndex));
%         end
%         fprintf(fID,'\n');
%     end
%     fclose(fID);
% end
%
% end
function writeImg(config)
mode=config.mode;
folderImg=config.folderImg;

cutoffFigOn = false;
% I0IFigOn = false;
fig1dOn = false;
fig2dOn = false;
SVFigOn = false;
scan2DFigOn = false;
pumpI0FigOn = false;
SCUCVFigOn = false;
OXCcompFigOn = false;
I0FigOn = false;

try
    cutoffFig=config.cutoffFig;
    cutoffFigOn=true;
catch
    fprintf('  Error encounter at saving cutoffFig\n');
end


% try
%     I0IFig=config.I0IFig;
%     I0IFigOn=true;
% catch
%     fprintf('  Error encounter at saving I0IFig\n');
% end


switch mode
    case 'TRXL'
        try
            fig1d=config.fig1d;
            fig1dOn=true;
        catch
            fprintf('  Error encounter at saving fig1d\n');
        end
        try
            fig2d=config.fig2d;
            fig2dOn=true;
        catch
            fprintf('  Error encounter at saving fig2d\n');
        end
        try
            SVFig=config.SVFig;
            SVFigOn=true;
        catch
            fprintf('  Error encounter at saving SVD figures\n');
        end
        try
            SCUCVFig = config.SCUCVFig;
            SCUCVFigOn=true;
        catch
            fprintf('  Error encounter at saving SCUCV\n');
        end
        try
            pumpI0Fig=config.pumpI0Fig;
            pumpI0FigOn=true;
        catch
            fprintf('  Error encounter at saving pumpI0Fig\n');
        end
        try
            OXCcompFig = config.OXCcompFig;
            OXCcompFigOn =true;
        catch
            fprintf('  Error encounter at saving pumpI0Fig\n');
        end
    case '2Dscan'
        try
            scan2DFig=config.scan2DFig;
            scan2DFigOn=true;
        catch
            fprintf('  Error encounter at saving 2dScan\n');
        end
end

try
    I0Fig = config.I0Fig;
    I0FigOn = true;
catch
    fprintf('  Error encounter at saving I0Fig\n');
end

% imgSaveStartComment;

if cutoffFigOn
    figureName=sprintf('%s/plotCutoff.png',folderImg);
    exportgraphics(cutoffFig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

% if I0IFigOn
%     figureName=sprintf('%s/I0IFig.png',folderImg);
%     exportgraphics(I0IFig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
% end

if fig1dOn
    figureName=sprintf('%s/1dPlot.png',folderImg);
    exportgraphics(fig1d,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if fig2dOn
    figureName=sprintf('%s/2dPlot.png',folderImg);
    exportgraphics(fig2d,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if SVFigOn
    figureName=sprintf('%s/SVDPlot.png',folderImg);
    exportgraphics(SVFig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if SCUCVFigOn
    figureName=sprintf('%s/SCUCV.png',folderImg);
    exportgraphics(SCUCVFig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if pumpI0FigOn
    figureName=sprintf('%s/pumpI0Plot.png',folderImg);
    exportgraphics(pumpI0Fig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if OXCcompFigOn
    figureName=sprintf('%s/OXCcomp.png',folderImg);
    exportgraphics(OXCcompFig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if scan2DFigOn
    figureName=sprintf('%s/2dScanPlot.png',folderImg);
    exportgraphics(scan2DFig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

if I0FigOn
    figureName=sprintf('%s/I0.png',folderImg);
    exportgraphics(I0Fig,figureName,'Resolution',config.figResolution,'ContentType','image','BackgroundColor','white');
end

fcloseWithComment;
imgSaveDoneComment;
    function imgSaveStartComment
        message1=sprintf('  Start saving images.');
        fprintf('%s\n%s\n',message1);
    end
    function imgSaveDoneComment
        message1=sprintf('  Saving images done');
        fprintf('%s\n%s\n',message1);
        disp('--------------------------------------------------------------')
    end
end
function finishComment(timeStart)
message1=sprintf('  The entire process done');
message2=sprintf('  Starting time: %s',timeStart);
message3=sprintf('  Finished time: %s',datetime('now'));
fprintf('%s\n%s\n%s\n',message1,message2,message3);
disp('--------------------------------------------------------------')
end
%% SubSubFunction
function AUCs = AUC_Jkim(sqMatrix,q,ROI)
qTF = ROI(1) <= q & q <= ROI(2);
AUCs = trapz(q(qTF),sqMatrix(qTF,:));
end
function config=figureEncapsulator(config)

targetObj=findobj('type','figure','Name','AUC cutoff');
if ~isempty(targetObj)
    config.cutoffFig=targetObj;
end

targetObj=findobj('type','figure','Name','SVD');
if ~isempty(targetObj)
    config.SVFig=targetObj;
end

targetObj=findobj('type','figure','Name','rawdata 1D');
if ~isempty(targetObj)
    config.fig1d=targetObj;
end

targetObj=findobj('type','figure','Name','rawdata 2D');
if ~isempty(targetObj)
    config.fig2d=targetObj;
end

targetObj=findobj('type','figure','Name','pumpI0 Comparison');
if ~isempty(targetObj)
    config.pumpI0Fig=targetObj;
end

targetObj=findobj('type','figure','Name','2D scan position Finder');
if ~isempty(targetObj)
    config.scan2DFig=targetObj;
end


end
function fig2d=plot2D(figureNumber,diffPlot,qDraw,timedelay,timeLin,timeLog,titleName,config)
fig2d = figureWithoutStealing(figureNumber,'rawdata 2D');
try
qRangeForYlim = [1,5.5];
scalor = 1.2;
qTF = qRangeForYlim(1) <= q & q <= qRangeForYlim(2);
maxVal = scalor*max(abs(diffPlot(qTF,:)),[],'all');

ylimVal =[-maxVal maxVal];

    numTimeLog=length(timedelay(timeLog));
    resolution=100;
    %%% subplot 1
    % ax(1)=subplot(2,1,1);

    if sum(timeLog)<3
        isPlotLog = false;
        tiledlayout(1,1);
    else
        isPlotLog = true;
        tl = tiledlayout_Jkim(2,1);
        % figureAllocation;
        % fontAllocation;
        nexttile;
        if numTimeLog>10
            contourf(qDraw,timedelay(timeLog),diffPlot(:,timeLog)',resolution,'LineColor','None');
            set(gca, 'Yscale','log');
        else
            contourf(qDraw,1:length(timedelay(timeLog)),diffPlot(:,timeLog)',resolution,'LineColor','None');
        end
        figureAllocation;

        if numTimeLog>10
            tickValue=timepointTickMaker(timedelay(timeLog));
            yticks(tickValue);
        else
            tickValue=timedelay(timeLog);
            yticks(1:length(timedelay(timeLog)));
        end
        % tickValue=timepointTickMaker(timedelay(timeLog));
        % yticks(tickValue);
        tickLabels=timepointTickLabelMaker(tickValue);
        yticklabels(tickLabels);

        grid on

        % cMap=bluewhitered(200);
        axNow = nexttile(1);
        colormap(axNow,bluewhitered);
        % colorbar;

        fontAllocation;
        titleNameFull=sprintf('%s',titleName);
        title(titleNameFull,'FontWeight','bold','Interpreter','none');

        ylabel('');
        
        try
        lim=ylimVal;
        catch
            lim=caxis;
        end
    end

    %%% subplot 2
    nexttile;

    contourf(qDraw,timedelay(timeLin)/1e-12,diffPlot(:,timeLin)',resolution,'LineColor','None');
    figureAllocation;
    Ang=char(197);

    grid on

    if isPlotLog
        axNow = nexttile(2);
        caxis(lim);
        colormap(axNow,bluewhitered);
        cb = colorbar_Jkim;
        cb.Layout.Tile = 'east';
        % cb.Layout.Tile = 'eastoutside';
        % cbarObj = colorbar;
        % cb.TickDirection = 'out';
        % cb.LineWidth = 2;
        % cb.Ruler.TickLabelFormat = '%g';
        % cb.Ruler.Exponent = 0;
        % disp('');
    else
        axNow = nexttile(1);

        colormap(axNow,bluewhitered);
        cb = colorbar;
        cb.Layout.Tile = 'east';
        cb.TickDirection = 'out';
        cb.LineWidth = 2;
        cb.Ruler.TickLabelFormat = '%g';
        cb.Ruler.Exponent = 0;
    end
    % colorbar;

    % set(gca,'FontWeight','bold');
    fontAllocation;
    % titleNameFull=sprintf('%s (lin)',titleName);
    % title(titleNameFull,'FontSize',28,'FontWeight','bold','Interpreter','none');

    % text_Tiled(titleNameFull);

    ylabel('');
    xl = xlabel(['q (',Ang,'^-^1)']);
    fontSize = xl.FontSize;
    ylabel(tl,'timedelay','FontSize',fontSize,'FontWeight','bold');
catch
    fprintf('  ErrorEcounterDrawing | fig2d\n');
end

end
function fig1d = plot1D(figureNumber,diffPlot,q,qDraw,timedelay,chiFile,plotOption,qDrawRange,eh1rayMXAI_intMean,config)
fig1d=figureWithoutStealing(figureNumber,'plot1D');
try

numTimedelay=length(timedelay);
cMap=colormap(fig1d,cool(numTimedelay));


qRangeForYlim = [1,5.5];
scalor = 1.2;
qTF = qRangeForYlim(1) <= q & q <= qRangeForYlim(2);
maxVal = scalor*max(abs(diffPlot(qTF,:)),[],'all');

ylimVal =[-maxVal maxVal];


%%% subplot 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
tiledlayout_Jkim(4,1);
nexttile([3 1]);
pltTarget=[];
legendCounter=0;

td_tick_idx = [];
for idx=1:numTimedelay
    if 1~=idx
        hold on
    else
    end
    PLT=plot(qDraw,diffPlot(:,idx),'Color',cMap(numTimedelay+1-idx,:));
    if any(timedelay(idx)==100e-12)
        PLT.LineWidth=2;
    elseif  95e-12 <= timedelay(idx) && timedelay(idx) <= 105e-12
        PLT.LineWidth=2;
    end
    alpha 'clear'
    if numTimedelay >= 7

        remNumber = round(numTimedelay/6);

        [~,zeroIdx] = min(abs(timedelay));
        % zeroIdx = find(zeroIdx);

        if rem(idx,remNumber)==0 || idx == 1 || idx == numTimedelay || idx == zeroIdx
            pltTarget=[pltTarget PLT];

            legendCounter=legendCounter+1;
            textBuff=timepoint2Str(timedelay(idx));
            legendComment(legendCounter)=string(textBuff);
            td_tick_idx = [td_tick_idx idx];
        end
    else
        pltTarget=[pltTarget PLT];
        legendCounter=legendCounter+1;
        textBuff=timepoint2Str(timedelay(idx));
        legendComment(legendCounter)=string(textBuff);
        td_tick_idx = [td_tick_idx idx];
    end
end
% yline(0,'--','Color','k');

figureAllocation;
xlabel('');
xticklabels([]);
axis tight;




qDrawRange = ROIComparison(qDrawRange,q);
xlim([qDrawRange(1),qDrawRange(2)]);
ylim(ylimVal);
grid on;
hold off
% set(gca,'FontSize',15,'FontWeight','bold');
fontAllocation;
title(chiFile,'FontWeight','bold','Interpreter','none');
switch plotOption
    case 'Sq'
        ylabel('\DeltaS(q)');
    case 'qSq'
        ylabel('q\DeltaS(q)');
    otherwise
        error('Error on plotOption.')
end
% legend(pltTarget,legendComment);
% legend('boxoff');
% legend('FontSize',12,'Orientation','vertical','Location','best','NumColumns',11);
axNow = gca;


if config.numH5~=1
    colorbar_plot1d(axNow,legendComment);
end
% numC = numTimedelay;
% cIdx = linspace(1,numC,numC);
% dC = 1/numC;
% tickVal = [0 linspace(cIdx(1+1),cIdx(end-1),numC-2)/(numC)-0.5*dC 1];
% colorbar(axNow,'TickLabels',legendComment,'Ticks',tickVal,'TickDirection','out','LineWidth',2);



%%% subplot 2
nexttile;
plot(q,eh1rayMXAI_intMean,'r','LineWidth',0.3);

Ang=char(197);
xlabel(['q (',Ang,'^-^1)']);
axis tight;
ylimPre=ylim;
figureAllocation;
qDrawRange = ROIComparison(qDrawRange,q);
xlim([qDrawRange(1),qDrawRange(2)]);
% xlim([0.05,2.2]);
ylim([0,ylimPre(2)]);

% set(gca,'FontSize',15,'FontWeight','bold');
fontAllocation;
ylabel('S(q)');
xlabel(['q (',Ang,'^-^1)']);
catch
end
end
function plotCutoffAUC(h5Struct,cutoffAUC,plotIndex)
AUCvalue = h5Struct(end).diffAUC;
% laserOn = h5Struct(end).laserOn;
I0good = h5Struct(end).I0good;
goodSignalCheck = h5Struct(end).goodSignalCheck;

% AUCvalue = AUCvalue(laserOn);

numAUC=length(AUCvalue);

curveGoodAUC = AUCvalue>=cutoffAUC;
badByAUC = ~goodSignalCheck&curveGoodAUC;
badByI0 = ~goodSignalCheck&~curveGoodAUC;

x1 = find(goodSignalCheck);
x2 = find(badByAUC);
x3 = find(badByI0);

AUC1 = AUCvalue(goodSignalCheck);
AUC2 = AUCvalue(badByAUC);
AUC3 = AUCvalue(badByI0);

% [AUCselectList,AUCselectValue,numAUC,numAUCselect]=selectAUC(goodSignalCheck,AUCvalue,cutoffAUC);
plot(x1,AUC1,'ko','LineWidth',1.3);
hold on
try
    % plot(AUCselectList,AUCselectValue,'ko','MarkerFaceColor','k','MarkerSize',5);
    plot(x2,AUC2,'ro','LineWidth',1.3);
    plot(x3,AUC3,'ro','MarkerFaceColor','r');
catch
    error('There is no data above a cutoff value of AUC. Check autoAUCValue.');
end
figureAllocation;
fontAllocation;
yline(cutoffAUC,'--','Cutoff Line','LineWidth',1,'Color','r');

counter = 0;

legendName1=sprintf('Good    | %3d/%3d (%.1f %%)',sum(goodSignalCheck),numAUC,sum(goodSignalCheck)/numAUC*100);
counter = counter + 1;
legendCell{counter} = legendName1;

if any(badByAUC)
    legendName2=sprintf('Bad     | %3d/%3d (%.1f %%)',sum(badByAUC),numAUC,sum(badByAUC)/numAUC*100);
    counter = counter + 1;
    legendCell{counter} = legendName2;
end

if any(badByI0)
    legendName3=sprintf('Bad(I0) | %3d/%3d (%.1f %%)',sum(badByI0),numAUC,sum(badByI0)/numAUC*100);
    counter = counter + 1;
    legendCell{counter} = legendName3;
end

legend(legendCell);
legend('boxoff');
legend('Location','Best');


[~,nameWOext,~]=fileparts(h5Struct(plotIndex).chiName);
titleName=sprintf('%s',nameWOext);
% title(titleName,'Interpreter','none');
xlabel('Shot');
ylabel('AUC');
xlim([0,numAUC]);
hold off

ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.3f');

end



function figOXC = plotOXCComparison(figureNumber,config)

fontSize = 8;
idx = 1;

timedelay=config.timedelay;
timedelaywoOXC = config.timedelaywoOXC;
% U=config.U;
% S=config.S;
V=config.V;
VwoOXC = config.VwoOXC;
% q=config.q;
qIsSVD=config.qIsSVD;
timedelayDrawPoint=config.timedelayDrawPoint;
numInt=config.numInt;
plotOption=config.plotOption;
numTimedelay=config.numTimedelay;

timepointROI = [-0.1 1]*1e-12;


figOXC = figureWithoutStealing(figureNumber,'OXC comparison');


%% 
tl = tiledlayout_Jkim(3,1);
%%
% tiledlayout_Jkim_arrange(3,1)

% set(gca,'FontWeight','bold');

timeLin = timedelay <= timedelayDrawPoint;
timeLog = 1e-15< timedelay;

timeLinwoOXC = timedelaywoOXC <= timedelayDrawPoint;
timeLogwoOXC = 1e-15< timedelaywoOXC;

numTimeLog=length(timedelay(timeLog));
% tiledlayout_space_Jkim(1,2);


% nexttile;
colorCode1 = sprintf('ko');
colorCode2 = sprintf('-ko');
colorCode1woOXC = sprintf('ro');
colorCode2woOXC = sprintf('-ro');

TF = timepointROI(1) < timedelay & timedelay < timepointROI(2);
TFwoOXC = timepointROI(1) < timedelaywoOXC & timedelaywoOXC < timepointROI(2);

Vtarget = V(TF,1);
VtargetwoOXC = VwoOXC(TFwoOXC,1);


[~,I] = max(abs(Vtarget));
val1 = Vtarget(I);

[~,I] = max(abs(VtargetwoOXC));
val2 = VtargetwoOXC(I);
signer = val1/val2;
comment = sprintf('V%d',idx);

if any(timeLin)

    nexttile([1 1]);
    plot(timedelay(timeLin)/1e-12,V(timeLin,idx),colorCode2);
    hold on;
    plot(timedelaywoOXC(timeLinwoOXC)/1e-12,signer*VwoOXC(timeLinwoOXC,idx),colorCode2woOXC);

    figureAllocation;
    fontAllocation(fontSize);
    % title(titleName);
    text_Tiled(comment);
    axis tight;
    set(gca,'FontWeight','bold');

    xlim([min(timedelaywoOXC(timeLinwoOXC)/1e-12) max(timedelaywoOXC(timeLinwoOXC)/1e-12)]);
    ylim([min(signer*VwoOXC(:,idx)),max(signer*VwoOXC(:,idx))]);
    % xlabel('timedelay (ps)','FontSize',16);
    ylabel('rSV (a.u.)','FontSize',16);

    % yticklabels([]);
    % ylabel([]);
    % yliner;
    xlabel([]);
    grid on;
    hold off;

    legend({'V1 w/ OXC','V1 w/o OXC'},'FontSize',16);
legend('box','off');
% axis tight;

nexttile([1 1]);
    plot(timedelay(timeLin)/1e-12,V(timeLin,idx),colorCode2);
    % hold on;
    % plot(timedelaywoOXC(timeLinwoOXC)/1e-12,signer*VwoOXC(timeLinwoOXC,idx),colorCode2woOXC);

    figureAllocation;
    fontAllocation(fontSize);
    % title(titleName);
    text_Tiled(comment);
    axis tight;
    set(gca,'FontWeight','bold');

    xlim([min(timedelaywoOXC(timeLinwoOXC)/1e-12) max(timedelaywoOXC(timeLinwoOXC)/1e-12)]);
    ylim([min(signer*VwoOXC(:,idx)),max(signer*VwoOXC(:,idx))]);
    % xlabel('timedelay (ps)','FontSize',16);
    ylabel('rSV (a.u.)','FontSize',16);

    % yticklabels([]);
    % ylabel([]);
    % yliner;
    xlabel([]);
    grid on;
    % hold off;

    legend({'V1 w/ OXC'},'FontSize',16);
legend('box','off');
% xlim([min(timedelaywoOXC(timeLinwoOXC)) max(timedelaywoOXC(timeLinwoOXC))]);
% axis tight;

    nexttile([1 1]);
    % plot(timedelay(timeLin)/1e-12,V(timeLin,idx),colorCode2);
    % hold on;
    plot(timedelaywoOXC(timeLinwoOXC)/1e-12,signer*VwoOXC(timeLinwoOXC,idx),colorCode2woOXC);

    figureAllocation;
    fontAllocation(fontSize);
    % title(titleName);
    text_Tiled(comment);
    axis tight;
    set(gca,'FontWeight','bold');

    xlim([min(timedelaywoOXC(timeLinwoOXC)/1e-12) max(timedelaywoOXC(timeLinwoOXC)/1e-12)]);
    ylim([min(signer*VwoOXC(:,idx)),max(signer*VwoOXC(:,idx))]);
    xlabel('timedelay (ps)','FontSize',16);
    ylabel('rSV (a.u.)','FontSize',16);

    % yticklabels([]);
    % ylabel([]);
    % yliner;
    grid on;
    % hold off;

    
    legend({'V1 w/o OXC'},'FontSize',16);
    legend('box','off');

    % xlim([min(timedelaywoOXC(timeLinwoOXC)) max(timedelaywoOXC(timeLinwoOXC))]);
    % axis tight;

    
end

end


function figObj = plotOXCstatic(figureNumber,h5Struct)
figObj = figureWithoutStealing(figureNumber,'OXC static');
try
name = h5Struct(end).name;
jitter = h5Struct(end).jitter;
oxcVal = h5Struct(end).eh1oxc_pos;


scatter(oxcVal,jitter*1e15,'filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0,'MarkerFaceColor','k');
figureAllocation;
fontAllocation;
% axis tight;
ylim([min(jitter*1e15) max(jitter*1e15)]);
title(name,'FontSize',20,'Interpreter','none');
xlabel('OXC pixel','FontSize',18);
ylabel('jitter (fs)','FontSize',18);

stdJitter = std(jitter*1e15,1,2);
meanJitter = mean(jitter*1e15,2);

disp('--------------------------------------------------------------');
cmt1 = sprintf('  STD of jitter | %.2f fs',stdJitter);
cmt2 = sprintf('  Mean of jitter | %.2f fs',meanJitter);

fprintf('%s\n',cmt1);
fprintf('%s\n',cmt2);

disp('--------------------------------------------------------------');
textComment = {...
cmt1...
cmt2...
    };
text_Jkim(textComment,0.95,0.9,12);

catch
end
end




function fcloseWithComment
status=fclose('all');
if status==-1
    fprintf('  fcloseError\n');
end
end
function figureObj=figureWithoutStealing(figureNumber,figureName,figSize)
testObj=findobj('type','figure','Name',figureName);
if ~isempty(testObj)
    figureObj=testObj;
    set(0,'CurrentFigure',figureObj);
else
    figureObj=figure(figureNumber);
    figureObj.NumberTitle='Off';
    figureObj.Name=figureName;
    figureObj.Color=[1 1 1];
    if nargin == 3
        figureObj.Position = [figureObj.Position(1) 10 figSize];
    end
end
end
function [h5PairList,h5ExistIdx,h5NoExistIdx]=getH5ProcessListWithSkip(h5Struct,h5StructPre)
numH5Pre=length(h5StructPre);
numH5=length(h5Struct);
nameStringPre=strings(1,numH5Pre);

% h5ExistList=[];
if ~isfield(h5StructPre,'folder')
    h5ExistIdx=[];
    h5PairList=[];
    h5NoExistIdx=1:numH5;
    fprintf('  Empty h5 StructPre\n');
elseif ~strcmp(h5Struct(1).folder,h5StructPre(1).folder)
    h5ExistIdx=[];
    h5PairList=[];
    h5NoExistIdx=1:numH5;
    fprintf('  Changed target folder of h5 Struct\n');
else
    h5PairList=[];
    for idx=1:numH5Pre
        nameStringPre(idx)=h5StructPre(idx).name;
    end
    for idx=1:numH5
        nameNow=h5Struct(idx).name;
        isMatch=nameNow==nameStringPre;
        if any(isMatch)
            h5ExistBool(idx)=true;
            h5PairList=[h5PairList find(isMatch)];
        else
            h5ExistBool(idx)=false;
        end
    end
    h5ExistIdx=find(h5ExistBool);
    h5NoExistIdx=find(~h5ExistBool);
end
fprintf('  Already read h5 file(s) skipped\n\n');
fprintf('  Need to read h5 file(s) | %d\n',length(h5NoExistIdx));
fprintf('  Total h5 file(s)        | %d\n\n',numH5);
fprintf('  Reading h5 file(s) in progress\n');
disp('--------------------------------------------------------------')
end
function [position1,position2,position3]=getPosition(name)
positionUnderscore=strfind(name,'_');
positionPoint=strfind(name,'.');
positionString3=name(1:positionUnderscore(1)-1);
positionString2=name(positionUnderscore(1)+1:positionUnderscore(2)-1);
positionString1=name(positionUnderscore(2)+1:positionPoint-1);

position1=str2double(positionString1);
position2=str2double(positionString2);
position3=str2double(positionString3);
end
% function scan=getScan(name)
% postionUnderscore=strfind(name,'_');
% scanString=name(postionUnderscore(1)+1:postionUnderscore(2)-1);
% scan=str2double(scanString);
% end
function run=getRun(folder,h5header)
positionTarget=strfind(folder,h5header)+length(h5header)+1;
runString=folder(positionTarget:positionTarget+4);
run=str2double(runString);
end
function PID=getPID(timestampPIDstring)
positionUnderscore=strfind(timestampPIDstring,'_');
PID=str2double(timestampPIDstring(positionUnderscore+1:end));
end
% function h5StructBadRemove = getH5BadRemove(h5Struct, q, solventRange, polymerRange, targetRatio)
% qIsSolventRange = solventRange(1) <= q & q <= solventRange(2);
% qSolventRange=solventRange(2)-solventRange(1);
% qIsPolymerRange = polymerRange(1) <= q & q <= polymerRange(2);
% qPolymerRange=polymerRange(2)-polymerRange(1);
%
% numh5=length(h5Struct);
%
% for idx=1:numh5
%     h5File=h5Struct(idx);
%     numShot=length(h5File.PID);
%     ratioList=nan(1,numShot);
%     for shotIndex=1:numShot
%         solventArea=trapz(q(qIsSolventRange),h5File.eh1rayMXAI_int(qIsSolventRange,shotIndex))/qSolventRange;
%         polymerArea=trapz(q(qIsPolymerRange),h5File.eh1rayMXAI_int(qIsPolymerRange,shotIndex))/qPolymerRange;
%         ratioList(shotIndex)=solventArea/polymerArea;
%     end
%
%     if 1==idx
%         ratioTotalList=ratioList;
%     else
%         ratioTotalList=[ratioTotalList,ratioList];
%     end
%
%     goodeh1rayMXAI_int=(ratioList>=targetRatio);
%     h5StructBadRemove(idx).name=h5File.name;
%     h5StructBadRemove(idx).eh1rayMXAI_int=h5File.eh1rayMXAI_int(:,goodeh1rayMXAI_int);
%     h5StructBadRemove(idx).PID=h5File.PID(goodeh1rayMXAI_int);
%     h5StructBadRemove(idx).timestampPIDString=h5File.timestampPIDString(goodeh1rayMXAI_int);
%     h5StructBadRemove(idx).laserOn=h5File.laserOn(goodeh1rayMXAI_int);
%     h5StructBadRemove(idx).position=h5File.position;
%     h5StructBadRemove(idx).scan=h5File.scan;
%     h5StructBadRemove(idx).timedelay=h5File.timedelay;
%     h5StructBadRemove(idx).run=h5File.run;
%
%
% end
% plotPolymerCutoff(ratioTotalList,targetRatio);
%     getBadRemoveCurvesComment(targetRatio,ratioTotalList);
%     if ~any(ratioTotalList>=targetRatio)
%         error('There is no good eh1rayMXAI_int. Check eh1rayMXAI_int process Values.')
%     end
%
%     function plotPolymerCutoff(ratioTotalList,targetRatio)
%         cutoffAUCFig=figure(11111);
%         cutoffAUCFig.NumberTitle='Off';
%         cutoffAUCFig.Name='Solvent/Polymer Cutoff';
%         [ratioSelectList,ratioSelectValue,numRatio,numRatioSelect]=selectAUC(ratioTotalList,targetRatio);
%         plot(ratioTotalList,'ko');
%         hold on
%         plot(ratioSelectList,ratioSelectValue,'ro');
%         hline=refline([0 targetRatio]);
%         hline.Color='r';
%
%         legendName1=sprintf('Good data : %d/%d (%.1f %%)',numRatio-numRatioSelect,numRatio,100-numRatioSelect/numRatio*100);
%         legendName2=sprintf('Bad data : %d/%d (%.1f %%)',numRatioSelect,numRatio,numRatioSelect/numRatio*100);
%         legendName3=sprintf('Ratio cutoff : %g',targetRatio);
%         legend({legendName1,legendName2,legendName3});
%         legend('boxoff','Location','Best');
%
%         titleName=sprintf('%s','Solvent/Polymer Ratio');
%         title(titleName,'Interpreter','none');
%         xlabel('# Shots')
%         ylabel('Ratio')
%         hold off
%     end
%     function getBadRemoveCurvesComment(targetRatio,ratioTotalList)
%     numCurveAll=length(ratioTotalList);
%     numCurveGood=sum(ratioTotalList>=targetRatio);
%     message1=sprintf('  Solvent/Polymer cutoff Value : %.3f',targetRatio);
%     message2=sprintf('  Bad eh1rayMXAI_int checked.');
%     message3=sprintf('  Good data %d/%d (%.2f%%)',numCurveGood,numCurveAll,numCurveGood/numCurveAll*100);
%     fprintf('%s\n%s\n%s\n',message1,message2,message3);
%     disp('--------------------------------------------------------------')
%     end
% end
function [qIsAUC,qIsNorm,qIsSVD]=indexingQ(q,AUCRange,normRange,SVDRange)
% AUC idx
qIsAUCRange=AUCRange(1)<=q & q<=AUCRange(2);
qIsAUC=find(qIsAUCRange);

% norm idx
qIsNorm=normRange(1)<=q & q<=normRange(2);
qIsNorm=find(qIsNorm);

% SVD idx
qIsSVD=SVDRange(1)<=q & q<=SVDRange(2);
qIsSVD=find(qIsSVD);
end
function [run,posi1List,posi2List,posi3List,numShotList]=getListFromH5(h5Struct,config)
numH5=config.numH5;

runList=nan(1,numH5);
posi1List=nan(1,numH5);
posi2List=nan(1,numH5);
posi3List=nan(1,numH5);
numShotList=nan(1,numH5);

for idx=1:numH5
    runList(idx)=h5Struct(idx).run;
    posi1List(idx)=h5Struct(idx).position1;
    posi2List(idx)=h5Struct(idx).position2;
    posi3List(idx)=h5Struct(idx).position3;
    numShotList(idx)=h5Struct(idx).numShot;
end
if length(unique(runList))>=2
    error('Not Imported Multi Run Process.')
else
    run=unique(runList);
end
end
% function chiStruct=pairing(chiStruct,runList,scanList,positionList,LaserOffList)
% chiStruct=getRefIndex(chiStruct,runList,scanList,positionList,LaserOffList);
% chiStruct=getShotPair(chiStruct);
% end
% function chiStruct=getRefIndex(chiStruct,runList,posi1List,posi2List,posi3List,LaserOffList)
% numChi=length(chiStruct);
% % refFile
% for fileIndex=1:numChi
%     chiFile=chiStruct(fileIndex);
%     chiStruct(fileIndex).refFile=findRefFile(chiFile,runList,posi1List,posi2List,posi3List,LaserOffList,fileIndex);
% end
% end
% function refFile=findRefFile(chiFile,runList,posi1List,posi2List,posi3List,LaserOffList,~)
% run=chiFile.run;
% position1=chiFile.position1;
% position2=chiFile.position2;
% position3=chiFile.position3;
% % LaserOff=chiFile.LaserOff;
% % if ~LaserOff % Coupling for positive chi files
%     targetIndexBool=(run==runList & position1==posi1List & position3==posi3List & position2==posi2List & LaserOffList);
% % elseif LaserOff % Coupling for reference chi files
% %     %     targetIndexBool=(run==runList & position1==posi1List & position3==posi3List & position2==posi2List & LaserOffList);
% %     if any(run==runList & position1-1==posi1List & position3==posi3List & position2==posi2List & LaserOffList)
% %         targetIndexBool=(run==runList & position1-1==posi1List & position3==posi3List & position2==posi2List & LaserOffList);
% %     elseif any(run==runList & position1+1==posi1List & position3==posi3List & position2==posi2List & LaserOffList)
% %         targetIndexBool=(run==runList & position1+1==posi1List & position3==posi3List & position2==posi2List & LaserOffList);
% %     elseif any(run==runList & position1+2==posi1List & position3==posi3List & position2==posi2List & LaserOffList)
% %         targetIndexBool=(run==runList & position1+2==posi1List & position3==posi3List & position2==posi2List & LaserOffList);
% %     elseif any(run==runList & position1-2==posi1List & position3==posi3List & position2==posi2List & LaserOffList)
% %         targetIndexBool=(run==runList & position1-2==posi1List & position3==posi3List & position2==posi2List & LaserOffList);
% %     elseif any(run==runList & max(posi1List)==posi1List & position3==posi3List & position2-1==posi2List & LaserOffList)
% %         targetIndexBool=(run==runList & max(posi1List)==posi1List & position3==posi3List & position2-1==posi2List & LaserOffList);
% %     else
% %         messageError=sprintf('NonRefIndex. RefIndex : %d',fileIndex);
% %         error(messageError);
% %     end
% % end
% refFile=find(targetIndexBool);
% end
% function [x1,x2,x3,AUC1,AUC2,AUC3]=selectAUC(goodSignalCheck,AUCvalue,cutoffAUC)
% numAUC=length(AUCvalue);
%
% curveGoodAUC = AUCvalue>=cutoffAUC;
%
% badByAUC = ~goodSignalCheck&curveGoodAUC;
% badByI0 = ~goodSignalCheck&~curveGoodAUC;
%
% x1 = find(goodSignalCheck);
% x2 = find(badByAUC);
% x3 = find(badByI0);
%
% AUC1 =
% AUC2 =
% AUC3 =



% counter1=0;
% counter2=0;
% AUCselectList=[];
% AUCselectValue=[];
% for idx=1:numAUC
%     if goodSignalCheck(idx)
%         counter1=counter1+1;
%         AUCselectList(counter1)=idx;
%         AUCselectValue(counter1)=AUCvalue(idx);
%     end
% end
% numAUCselect=length(AUCselectList);
% end
% function h5Struct=getDiffBadRemoveCurves(h5Struct,config)
% numH5=config.numH5;
% shotListBadRemove=nan(1,numH5);
%
% for idx=1:numH5
%     chiFile=h5Struct(idx);
%     diffOrg=chiFile.diff;
%     goodSignalCheck=chiFile.goodSignalCheck;
%     numGoodSignal=sum(chiFile.goodSignalCheck);
%
%
%     diffBadRemove=diffOrg(:,goodSignalCheck);
%     h5Struct(idx).diff=diffBadRemove;
%     shotListBadRemove(idx)=numGoodSignal;
% end
% getBadRemoveCurvesComment(shotListBadRemove);
% end
% function diffStruct=getAveSTDpumpI0(timedelay,chiStruct,diffStruct,q)
% numTimedelay=length(timedelay);
% numQ=length(q);
% for timeIndex=1:numTimedelay
%     numShotBadRemoveSum=0;
%     for diffIndex=1:length(diffStruct(timeIndex).file)
%         if diffIndex==1
%             diffs=chiStruct(diffStruct(timeIndex).file(diffIndex)).diffBadRemove;
%
%         else
%             diffs=[diffs,chiStruct(diffStruct(timeIndex).file(diffIndex)).diffBadRemove]; %#ok<*AGROW>
%         end
%         numShotBadRemoveSum=numShotBadRemoveSum+chiStruct(diffStruct(timeIndex).file(diffIndex)).numShotBadRemove;
%     end
%       diffStruct(timeIndex).numImage=numShotBadRemoveSum;
%     if 0==diffStruct(timeIndex).numImage
%         diffStruct(timeIndex).diffAve=zeros(numQ,1);
%         diffStruct(timeIndex).diffSTDMean=zeros(numQ,1);
%         fprintf('  Warning : There is Zero Image Timedelay : %.3f ps\n',timedelay(timeIndex)*1e12);
%     else
%     diffStruct(timeIndex).diffAve=mean(diffs,2);
%     diffStruct(timeIndex).diffSTDMean=std(diffs,0,2)/(sqrt(diffStruct(timeIndex).numImage));
%     end
% end
% end
% function chiStruct=diffAveSTDEachScan(chiStruct)
% numChi=length(chiStruct);
% for fileIndex=1:numChi
%     chiStruct(fileIndex).diffMean=mean(chiStruct(fileIndex).diffBadRemove,2);
%     chiStruct(fileIndex).diffSTDMean=std(chiStruct(fileIndex).diffBadRemove,0,2)/(sqrt(chiStruct(fileIndex).numShotBadRemove));
% end
% end
function [ordinalNumber,colorName]=ordinalGet(idx)
if idx==1
    ordinalNumber='1st';
    colorName='k';
elseif idx==2
    ordinalNumber='2nd';
    colorName='r';
elseif idx==3
    ordinalNumber='3rd';
    colorName='b';
elseif idx==4
    ordinalNumber='4th';
    colorName='m';
else
    ordinalNumber=sprintf('%dth',idx);

    switch mod(idx,4)

        case 1
            colorName='k';
        case 2
            colorName='r';
        case 3
            colorName='b';
        case 0
            colorName='m';
    end
end
end
% function rSV1stFig=drawrSV1stFig(figureNumber,runList,timedelay,V)
% rSV1stFig=figure(figureNumber);
% rSV1stFig.NumberTitle='Off';
% rSV1stFig.Name='1st rSV';
% rSV1stFig.Color=[1 1 1];
% titleFile=sprintf('1st rSV of Run%.3d',unique(runList));
% plot( (timedelay(2:end) - min( timedelay(2:end))*1.5),V(2:end,1),'-ko');
% title(titleFile);
% axis tight;
% set(gca, 'XScale', 'Log');
% set(gca, 'XTick', [linspace(floor(min(timedelay(2:end))*1e12), 0, abs(floor(min(timedelay(2:end))*1e12) - 0)+1)*1e-12 logspace( -12, -10, 3)] - min( timedelay(2:end))*1.5  )
% set(gca, 'XTickLabel', [linspace(floor(min(timedelay(2:end))*1e12), 0, abs(floor(min(timedelay(2:end))*1e12) - 0)+1)*1e-12 logspace( -12, -10, 3)])
% grid on;
% xlabel('timedelay');
% end
% function drawSVFig(runList,timedelay,SV)
% SVFig=figure(77777);
% SVFig.NumberTitle='Off';
% SVFig.Name='SV';
%
% PLT=plot((timedelay(2:end)-min(timedelay(2:end))*1.5),SV(2:end,1:4),'-o');
% PLT(1).Color='k';
% PLT(2).Color='r';
% PLT(3).Color='b';
% PLT(4).Color='m';
%
% textBuff=sprintf('1st SV');
% legendComment(1)=string(textBuff);
% textBuff=sprintf('2nd SV');
% legendComment(2)=string(textBuff);
% textBuff=sprintf('3rd SV');
% legendComment(3)=string(textBuff);
% textBuff=sprintf('4th SV');
% legendComment(4)=string(textBuff);
%
% titleName=sprintf('SV of Run%.3d',unique(runList));
% title(titleName);
%
% axis tight;
% set(gca, 'XScale', 'Log');
% set(gca, 'XTick', [linspace(floor(min(timedelay(2:end))*1e12), 0, abs(floor(min(timedelay(2:end))*1e12) - 0)+1)*1e-12 logspace( -12, -10, 3)] - min( timedelay(2:end))*1.5  )
% set(gca, 'XTickLabel', [linspace(floor(min(timedelay(2:end))*1e12), 0, abs(floor(min(timedelay(2:end))*1e12) - 0)+1)*1e-12 logspace( -12, -10, 3)])
% grid on;
%
% legend(legendComment);
% legend('boxoff','Location','Best');
% end
% function drawrSVFig(runList,timedelay,V)
% rSVFig=figure(44444);
% rSVFig.NumberTitle='Off';
% rSVFig.Name='Right Singular Vectors';
% for vectorIndex=1:4
%     drawrSVsubFig(runList,timedelay,V,vectorIndex);
% end
% end
% function drawrSVsubFig(runList,timedelay,V,idx)
% if idx==1
%     titleName=sprintf('1st rSV of Run%.3d',unique(runList));
%     colorName='-ko';
% elseif idx==2
%     titleName=sprintf('2nd rSV of Run%.3d',unique(runList));
%     colorName='-ro';
% elseif idx==3
%     titleName=sprintf('3rd rSV of Run%.3d',unique(runList));
%     colorName='-bo';
% elseif idx==4
%     titleName=sprintf('4th rSV of Run%.3d',unique(runList));
%     colorName='-mo';
% else
%     titleName=sprintf('%dth rSV of Run%.3d',idx,unique(runList));
%     colorName='-ko';
% end
%
% subplot(4,1,idx)
% plot( (timedelay(2:end) - min( timedelay(2:end))*1.5),V(2:end,idx),colorName);
% title(titleName);
% axis tight;
% set(gca, 'XScale', 'Log');
% set(gca, 'XTick', [linspace(floor(min(timedelay(2:end))*1e12), 0, abs(floor(min(timedelay(2:end))*1e12) - 0)+1)*1e-12 logspace( -12, -10, 3)] - min( timedelay(2:end))*1.5  )
% set(gca, 'XTickLabel', [linspace(floor(min(timedelay(2:end))*1e12), 0, abs(floor(min(timedelay(2:end))*1e12) - 0)+1)*1e-12 logspace( -12, -10, 3)])
% grid on;
% end
% function drawlSVFig(runList,q,qIsSVD,U)
% lSVFig=figure(55555);
% lSVFig.NumberTitle='Off';
% lSVFig.Name='Left Singular Vectors';
%
%
% titleFile=sprintf('lSV of Run%.3d',unique(runList));
% plot(q(qIsSVD),U(:,1),'k',q(qIsSVD),U(:,2),'r',q(qIsSVD),U(:,3),'b',q(qIsSVD),U(:,4),'m');
% title(titleFile);
% hline=refline([0 0]);
% hline.Color='k';
%
% xlabel('q (A^-^1)');
% ylabel('\DeltaS');
% legend('1st lSV','2nd lSV','3rd lSV','4th lSV')
% legend('boxoff','Location','Best');
% grid on
% end
function chiStruct=pairingTotalPool(chiStruct,refAUCPool)

numChi=length(chiStruct);
for fileIndex=1:numChi
    for targetIndex=1:chiStruct(fileIndex).numShot
        AUCDifference=abs(refAUCPool-chiStruct(fileIndex).AUC(targetIndex));
        sortAUCDiff=sort(AUCDifference);
        minValue=sortAUCDiff(2);
        AUCMinIndex=find(AUCDifference==minValue);
        chiStruct(fileIndex).shotPair(targetIndex)=AUCMinIndex(1);
    end
    chiStruct(fileIndex).refFile=0;
end
end
function [refeh1rayMXAI_intPool,refNormfactorPool,refAUCPool]=getRefeh1rayMXAI_intPool(chiStruct)
numChiPool=length(chiStruct);
counter=0;
for chiIndex=1:numChiPool
    LaserOff=chiStruct(chiIndex).LaserOff;
    if LaserOff
        counter=counter+1;
        if 1==counter
            refeh1rayMXAI_intPool=chiStruct(chiIndex).eh1rayMXAI_int;
            refNormfactorPool=chiStruct(chiIndex).normfactor;
            refAUCPool=chiStruct(chiIndex).AUC;
        else
            refeh1rayMXAI_intPool=[refeh1rayMXAI_intPool,chiStruct(chiIndex).eh1rayMXAI_int];
            refNormfactorPool=[refNormfactorPool,chiStruct(chiIndex).normfactor];
            refAUCPool=[refAUCPool,chiStruct(chiIndex).AUC];
        end
    end

end
end
% function chiStruct=getDiffAUCTotalPool(chiStruct,q,qIsAUC,refeh1rayMXAI_intPool,refNormfactorPool)
% numChi=length(chiStruct);
% numQ=length(q);
% for fileIndex=1:numChi
%     numShot=chiStruct(fileIndex).numShot;
%     chiStruct(fileIndex).diffOrg=nan(numQ,numShot);
%     for targetIndex=1:numShot
%         targetShot=chiStruct(fileIndex).shotPair(targetIndex);
%         laserOnSignal=chiStruct(fileIndex).eh1rayMXAI_int(:,targetIndex)/chiStruct(fileIndex).normfactor(targetIndex);
%         laserOffSignal=refeh1rayMXAI_intPool(:,targetShot)/refNormfactorPool(targetShot);
%         diffValue=laserOnSignal-laserOffSignal;
%         chiStruct(fileIndex).diffAUC(targetIndex)=trapz(q(qIsAUC),abs(diffValue(qIsAUC)));
%         chiStruct(fileIndex).diffOrg(:,targetIndex)=diffValue;
%     end
% end
% end
% function diffStruct=getMotorValue2D(diffStruct,folderh5)
% pathh5=sprintf('%s/scanInfo',folderh5);
% fileName=sprintf('%s/scanInfo.h5',pathh5);
% motor1ValueByPosition = h5read(fileName,'/run/scan00001/motor/m1');
% motor2ValueByPosition = h5read(fileName,'/run/scan00001/motor/m2');
% for diffIndex=1:length(diffStruct)
%     position1=diffStruct(diffIndex).position1;
%     position2=diffStruct(diffIndex).position2;
%     position3=diffStruct(diffIndex).position3;
%     diffStruct(diffIndex).motor1=motor1ValueByPosition(position1);
%     diffStruct(diffIndex).motor2=motor2ValueByPosition(position2);
%     diffStruct(diffIndex).motor3=position3;
% end
% end
%
% function diffStruct=getMotorValueTRXL(diffStruct,folderh5)
% pathh5=sprintf('%s/scanInfo',folderh5);
% fileName=sprintf('%s/scanInfo.h5',pathh5);
% timedelayByPosition = h5read(fileName,'/run/scan00001/motor/m1');
% for diffIndex=1:length(diffStruct)
%     position1=diffStruct(diffIndex).position1;
%     position2=diffStruct(diffIndex).position2;
%     position3=diffStruct(diffIndex).position3;
%     diffStruct(diffIndex).timedelay=timedelayByPosition(position1)*1e-12;
%     diffStruct(diffIndex).motor2=position2;
%     diffStruct(diffIndex).motor3=position3;
% end
% end
function diffStruct=getDiffSum(diffStruct,diffSumRange,q)
qBool= diffSumRange;
for diffIndex=1:length(diffStruct)
    ave=diffStruct(diffIndex).diffAve;
    diffStruct(diffIndex).sum=trapz(q(qBool),abs(ave(qBool)));
end
end
function newmap = bluewhitered(m)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.

if nargin < 1
    m = size(get(gcf,'colormap'),1);
end


bottom = [0 0 0.5];
botmiddle = [0 0.5 1];
middle = [1 1 1];
topmiddle = [1 0 0];
top = [0.5 0 0];

% Find middle
lims = get(gca, 'CLim');

% Find ratio of negative to positive
if (lims(1) < 0) && (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;

    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end

    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end

    % And put 'em together
    newmap = [newmap1; newmap];

elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end

else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);

    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end

end
end
function tickValue=timepointTickMaker(timedelay)
logStart=ceil(min(log10(timedelay)));
logEnd=floor(max(log10(timedelay)));
tickValue=logspace(logStart,logEnd,logEnd-logStart+1);
end
function tickLabels=timepointTickLabelMaker(tickValue)
numTickValue=length(tickValue);
tickLabels={};
for idx=1:numTickValue
    tickLabels{idx}=timepoint2Str(tickValue(idx));
end
end
function timepointStr=timepoint2Str(timepoint)
timepointLog=log10(abs(timepoint));
if timepointLog<-15
    timeUnit=1e-15;
    timeUnitStr='fs';
elseif -15<=timepointLog && timepointLog<-12
    timeUnit=1e-15;
    timeUnitStr='fs';
elseif -12<=timepointLog && timepointLog<-9
    timeUnit=1e-12;
    timeUnitStr='ps';
elseif -9<=timepointLog && timepointLog<-6
    timeUnit=1e-9;
    timeUnitStr='ns';
elseif -6<=timepointLog && timepointLog<-3
    timeUnit=1e-6;
    timeUnitStr='us';
else
    timeUnit=1e-3;
    timeUnitStr='ms';
end
timepointStr=sprintf('%.0f %s',timepoint/timeUnit,timeUnitStr);
end
% function timedelay=getTimeDelay(targetStruct)
% numData=length(targetStruct);
% timedelayList=nan(1,numData);
% for idx=1:length(targetStruct)
%     timedelayList(idx)=targetStruct(idx).timedelay;
% end
% timedelay=unique(timedelayList);
% end
% function draw2DSV(diffStruct)
%
% numInt=4;
%
% numDiff=length(diffStruct);
%
% posiList1=nan(1,numDiff);
% posiList2=nan(1,numDiff);
% posiList3=nan(1,numDiff);
%
% motorList1=nan(1,numDiff);
% motorList2=nan(1,numDiff);
% motorList3=nan(1,numDiff);
%
% SV=nan(numInt,numDiff);
%
% for idx=1:numDiff
%     posiList1(idx)=diffStruct(idx).position1;
%     posiList2(idx)=diffStruct(idx).position2;
%     posiList3(idx)=diffStruct(idx).position3;
%
%     motorList1(idx)=diffStruct(idx).motor1;
%     motorList2(idx)=diffStruct(idx).motor2;
% %     motorList3=diffStruct(idx).motor3;
%
%     SV(:,idx)=diffStruct(idx).SV(1:numInt);
% end
%
% posi1=unique(posiList1);
% posi2=unique(posiList2);
%
% motor1=unique(motorList1);
% motor2=unique(motorList2);
%
%
% cutoffAUCFig=figure(20003);
% cutoffAUCFig.NumberTitle='Off';
% cutoffAUCFig.Name='2D V1 Plot';
%
% for intIndex=1:numInt
%
% subplot(2,2,intIndex);
%
% plotMatrix=nan(length(posi1),length(posi2));
% for idxPosi1=posi1
%     for idxPosi2=posi2
%         idx=(idxPosi1==posiList1 & idxPosi2 == posiList2);
% %         plotMatrix(idxPosi1,idxPosi2)=SV1(counter);
%         plotMatrix(idxPosi1,idxPosi2)=diffStruct(idx).SV(intIndex);
%     end
% end
%
%
%
%
% contourf(motor1,motor2,plotMatrix',50,'LineColor','None');
% % contourf(motor1,motor2,plotMatrix',50,'LineColor');
% colorbar;
% grid on;
%
% xlabel('motor position 1');
% ylabel('motor position 2');
%
% maxValue=max(abs(SV(intIndex,:)));
% maxiIndex=find(abs(SV(intIndex,:))==maxValue);
%
% textLegend=sprintf('S%d maxiAbs: %g motor1: (%g,%d),motor2: (%g,%d),motor3: (%g,%d)',intIndex,maxValue,motorList1(maxiIndex),posiList1(maxiIndex),motorList2(maxiIndex),posiList2(maxiIndex),motorList3(maxiIndex),posiList3(maxiIndex));
% title(textLegend);
% end
% end
% scatter2(posiList1,posiList2,ones(1,numDiff)*10,sumList,'filled');
% function [...
%     h5header,...
%     parentFolderH5,...
%     parentFolderResults,...
%     mode,...
%     XrayCenter,...
%     repetitionRate,...
%     normRange,...
%     SVDRange,...
%     AUCRange,...
%     drawRange,...
%     autoAUCValue,...
%     outputWriteOn,...
%     SVDoutputWriteOn,...
%     eachScanWrite,...
%     pumpI0On,...
%     pumpI0Correction,...
%     correctionFunc,...
%     verbose,...
%     folderResultsHeader ...
%     ]=decapsulatorConfig(config)
% h5header=config.h5header;
% parentFolderH5=config.parentFolderH5;
% parentFolderResults=config.parentFolderResults;
% mode=config.mode;
% XrayCenter=config.XrayCenter;
% repetitionRate=config.repetitionRate;
% normRange=config.normRange;
% SVDRange=config.SVDRange;
% AUCRange=config.AUCRange;
% drawRange=config.drawRange;
% autoAUCValue=config.autoAUCValue;
% outputWriteOn=config.outputWriteOn;
% SVDoutputWriteOn=config.SVDoutputWriteOn;
% eachScanWrite=config.eachScanWrite;
% pumpI0On=config.pumpI0On;
% pumpI0Correction=config.pumpI0Correction;
% correctionFunc=config.correctionFunc;
% verbose=config.verbose;
% folderResultsHeader=config.folderResultsHeader;
% end
function figureAllocation(val_line)
Ang=char(197);
xlabel(['q (',Ang,'^-^1)']);
ylabel('q\DeltaS(q)');
grid on
% legend('Box','off');
fig=gcf;
fig.Color=[1 1 1];
ax=gca;
% ax.FontSize=15;
ax.FontWeight='bold';
ax.LabelFontSizeMultiplier=1.3;
ax.TitleFontSizeMultiplier=1.3;
ax.TitleFontWeight='bold';
ax.XLimitMethod='tight';
ax.Box = 'on';
ax.LineWidth = 2;
ax.ColorOrder=[0 0 0
    1 0 0
    0 0 1
    1 0 1
    0 1 0
    0 1 1
    1 1 0
    ];
ax.TickDir = 'none';
switch nargin
    case 1
        yline(val_line,'--','Color',[0.8 0.8 0.8],'LineWidth',2);
end

% legend('off');
end
function tl = tiledlayout_Jkim(numRow,numCol)
tl = tiledlayout(numRow,numCol,'GridSize',[numRow numCol],'TileSpacing', 'none', 'Padding', 'tight');
end

% function tl = tiledlayout_Jkim_arrange(numRow,numCol)
% tl = tiledlayout(numRow,numCol,'GridSize',[numRow numCol],'TileSpacing', 'none', 'Padding', 'tight','arrangemnt','vertical');
% end
function tiledlayout_space_Jkim(numRow,numCol)
tiledlayout(numRow,numCol,'GridSize',[numRow numCol],'TileSpacing', 'tight', 'Padding', 'tight');
end
function qROInew= ROIComparison(qROInow,q)
if qROInow(1)<q(1)
    qROInew(1) = q(1);
else
    qROInew(1) = qROInow(1);
end

if qROInow(end)>q(end)
    qROInew(2) = q(end);
else
    qROInew(2) = qROInow(end);
end


end
function text_Tiled(comment,textSize,colorNow)
switch nargin
    case 1
        textSize = 10;
        colorNow ='k';
    case 2
        colorNow ='k';
    case 3
end
text(0.5,0.99,comment,'Units','normalized','HorizontalAlignment','center','FontWeight','bold','VerticalAlignment','top','FontSize',textSize,'Color',colorNow);
end
function list = listFromStr(Str,channelName)
list = [];
for idx = 1:length(Str)
    dataNow = Str(idx).(channelName);

    list = [list dataNow];
end
end
function plot_Jkim_mode(pltObj,td,mode,axisConstant)
%%% mode | 1: symlog, 2: log, 3: lin, 4: noaxis,
%%% axisConstant | mode(1): symlog constant, mode(4): sampling number

switch nargin
    case 4
    case 3
        axisConstant = -1;
    case 2
        mode = 1; % default is symlog
        axisConstant = -1;
    otherwise
        fprintf('Where is it, who am I')
        return;
end

switch mode
    case 1 % Case for symlog
        %% Section 1 | make tick and tick labeles
        td_symlog = make_td_symlog(td,axisConstant);
        td_tick = make_td_tick(td);
        td_tick_symlog = make_td_symlog(td_tick,axisConstant);
        % plot(td,dataVector,'Color','k','LineWidth',2);
        pltObj.XData = td_symlog;
        %% Apply tick and tick labels
        xticks(td_tick_symlog); % On symlog space
        xticklabels(td_tick); % On original space
        xlim([td_tick_symlog(1) td_tick_symlog(end)]);
    case 2 % Case for log
        clearvar axisConstant;
        set(gca, 'Xscale','log');
    case 3 % Case for lin
        clearvar axisConstant;
    case 4 % Case for no axis
        numTd = length(td);
        tdNow = 1:numTd;
        pltObj.XData = tdNow;
        if numTd < axisConstant
            axisConstant = numTd;
        end

        [tickValue,tickIdx] = makeNewTd(td,axisConstant);
        xticks(tickIdx);
        % pltObj.XData = 1:length(tickValue);

        tickLabels=timepointTickLabelMaker(tickValue);
        xticklabels(tickLabels);
end


    function [tdNew,tdIdx] = makeNewTd(td,axisConstant)
        % tdNew = [];
        % tdNew = [tdNew td(1)];
        tdIdx = round(linspace(1,length(td),axisConstant));
        tdNew = td(tdIdx);

        % tdIdx = [1 tdIdx length(td)];
    end

    function tickLabels=timepointTickLabelMaker(tickValue)
        numTickValue=length(tickValue);
        tickLabels={};
        for idx=1:numTickValue
            tickLabels{idx}=timepoint2Str(tickValue(idx));
        end
    end
    function timepointStr=timepoint2Str(timepoint)
        timepointLog=log10(abs(timepoint));
        if timepointLog<-15
            timeUnit=1e-18;
            timeUnitStr='as';
        elseif -15<=timepointLog && timepointLog<-12
            timeUnit=1e-15;
            timeUnitStr='fs';
        elseif -12<=timepointLog && timepointLog<-9
            timeUnit=1e-12;
            timeUnitStr='ps';
        elseif -9<=timepointLog && timepointLog<-6
            timeUnit=1e-9;
            timeUnitStr='ns';
        elseif -6<=timepointLog && timepointLog<-3
            timeUnit=1e-6;
            timeUnitStr='us';
        else
            timeUnit=1e-3;
            timeUnitStr='ms';
        end
        timepointStr=sprintf('%.0f %s',timepoint/timeUnit,timeUnitStr);
    end
end
function R2 = R2cal(expMat,theMat,p)
% If you put in only 2 vars. It will return only R2.
% But if you put with 3 vars. It will return Adjusted R2.
% expMat: could by n by 1 or n by n
% theoMat: could by n by 1 or n by n, but must same as expMat
% p: 1 by 1, number of fitting parameters
% https://m.blog.naver.com/tlrror9496/222055889079 (R2)
% https://kr.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
% (R2 adjusted)

[numDataMat1,numDataMat2] = size(expMat);
[numDataTheo1,numDataTheo2] = size(theMat);

if numDataMat1 ~= numDataTheo1 || numDataMat2 ~= numDataTheo2
    fprintf('  Error encounter. numDataMat ~= numDataTheo\n');
    return;
end

if numDataMat1 < numDataMat2
    expMat = expMat';
    theMat = theMat';
end

expMean = mean(expMat,1); % 1 이 맞다 어째서지? 설명서에 나와있음 ㅋㅋ

S = (expMat - theMat).^2;
SSR = sum(S,'all');

S = (expMat - expMean).^2;
SST = sum(S,'all');

switch nargin
    case 2
        R2 = 1-SSR/SST;
    case 3
        n = numel(expMat);
        R2 = 1-((n-1)/(n-p))*SSR/SST;
end
end

function score = RRScal(expMat,theMat,I0)
%% Section 1 | check input size
[numDataMat1,numDataMat2] = size(expMat);
[numDataTheo1,numDataTheo2] = size(theMat);
[numDataI01,numDataI02] = size(I0);

if ~(numDataTheo1 == numDataI01 && numDataTheo2 == numDataI02)
    fprintf('  Error encounter. numDataMat ~= numDataTheo\n');
    return;
end

%% Section 2 | Make all column vector
if numDataMat1 < numDataMat2
    expMat = expMat';
    theMat = theMat';
    I0 = I0';
end
%% Section 3 | Normalization
expMean = mean(expMat,1); % 1 이 맞다 어째서지? 설명서에 나와있음 ㅋㅋ
I0Mean = mean(I0,1); % 1 이 맞다 어째서지? 설명서에 나와있음 ㅋㅋ
expMat = expMat/expMean;
theMat = theMat/expMean;
I0 = I0/I0Mean;
%% Section 4 | Cal score function
S = ((expMat - theMat)).^2./I0;
score = sum(S,'all');
score = -log10(score);
% S = (expMat - expMean).^2;
% SST = sum(S,'all');

% switch nargin
%     case 2
%         R2 = 1-SSR/SST;
%     case 3
%         n = numel(expMat);
%         R2 = 1-((n-1)/(n-p))*SSR/SST;
% end
end



function colorbarObj = colorbar_plot1d(axNow,legendComment)
% colorbar for 1d plot
% axNow | axis object
% numC | scalar, #curves
% legendComment | string (1 by numC)

legendComment = flip(legendComment);

numC = length(legendComment);
% numC = numTimedelay;
cIdx = linspace(1,numC,numC);
dC = 1/numC;
tickVal = [0 linspace(cIdx(1+1),cIdx(end-1),numC-2)/(numC)-0.5*dC 1];
colorbarObj = colorbar(axNow,'TickLabels',legendComment,'Ticks',tickVal,'TickDirection','out','LineWidth',2);
colorbarObj.Direction = 'reverse';
end
function [weAreEndGameNow,progress] = DrStrange(numData,numDiff)
if numData == numDiff
    weAreEndGameNow =true;
    fprintf('  We''re in the endgame now\n');
else
    weAreEndGameNow =false;
    fprintf('  We''re not in the endgame yet\n');
end
progress = numDiff/numData;
fprintf('  Progress status | %.1f%% (%u/%u)\n',progress*100,numDiff,numData);
fprintf('  --------------------------------------------------------------\n');
end
function cb = colorbar_Jkim
cb = colorbar;
cb.TickDirection = 'out';
cb.LineWidth = 2;
cb.Ruler.TickLabelFormat = '%g';
cb.Ruler.Exponent = 0;
end
function fontAllocation(fontsize)

switch nargin
    case 0
        fontsize = 12;
    case 1
end
set(gca,'FontSize',fontsize,'FontWeight','bold');
end
function textObj = text_Jkim(txt,position1,position2,fontSize)
switch nargin
    case 1
        position1 = 0.97;
        position2 = 0.20;
    case 3
        fontSize = 12;
    case 4
    otherwise
        fprintf('  Not supported nargin (number of input)\n');
        return;
end

if isempty(position1)
    position1 = 0.97;
end

if isempty(position2)
    position2 = 0.20;
end

textObj = text(position1,position2,txt,'Units','normalized','FontWeight','bold','HorizontalAlignment','right','FontSize',fontSize);
end
%%% SubSubFunctions for symlog
function td_symlog = make_td_symlog(td,C)
C = 10^C;
td_symlog = sign(td).*log10(1+abs(td)/C); % new % 240310
end
function td_tick = make_td_tick(td)

td_posi_TF = td > 1e-15;
td_posi = td(td_posi_TF);

tick_max = max(real(log10(td_posi)));
tick_max = floor(tick_max);

tick_min = min(real(log10(td_posi)));
% tick_min = floor(tick_min); % Display 10 fs
tick_min = ceil(tick_min); % Display 100 fs

td_tick = [0 10.^[tick_min:1:tick_max]];

if any(td_tick<td(1))
else
    tdFirst = str2double(sprintf('%.1g',td(1)));

    td_tick = [tdFirst td_tick];
end
tdEnd = str2double(sprintf('%.1g',td(end)));
td_tick = [td_tick tdEnd];
end
%%% SubSubFunctions for Wellcome
function ASCArt
fprintf('  Congratulations! "Cheer up, raccoon" has been rooting for you!\n');
fprintf([...
    '  @@@@@@@@@@@@@@@@@@   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@@    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@@  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@   @@@@            @@@@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                      @@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@                      @@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@                      @@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                     @@@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                      @@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                       @@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@                         @@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@                         @@@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@                           @@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@                          @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@                         @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@@                       @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@@                       @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@                        @@@@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@       @@     @@@@         @@@@@@@@@@@@@@@@@\n'...
    '  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       @@@@@@@@@@@@@@@@\n'...
    ]);
end
function PPiAAc
fprintf([...
    '                      YQY                                \n'...
    '                   UOLYLLY      SNU WU                   \n'...
    '                    SQYYY       WNW SO                   \n'...
    '        YYW             UY          WU     SDLU          \n'...
    '     YOOJYQ           YUQLOQNLNSY          QDN           \n'...
    '      QQSO        QNWYY         YYYYQLNW                 \n'...
    '               UNY                      WQY              \n'...
    '             WNU                          YSU            \n'...
    '           YNU                              YOU          \n'...
    '          UQY                                 QU         \n'...
    '         UO                                    US        \n'...
    '        UH                                      WS       \n'...
    '        DJ                                       WW      \n'...
    '       UBY                                       WU      \n'...
    '       NJ                                        YUY     \n'...
    '      YHJ                                         UW     \n'...
    '      YHO        QY                     Y    Y    SW     \n'...
    '       LLU      WBJHJOLJDW            NNNOJBNFO   SW     \n'...
    '       SF       YYOJQUUWUY  WHDFDLUJ    YYUUWUW   UW     \n'...
    '       UNY          Y       NAH UU F              SW     \n'...
    '        OY                  WNAL  LU             YU      \n'...
    '         QW                    QBLY              UY      \n'...
    '         YOU                                    WU       \n'...
    '          WLW                                  SS        \n'...
    '           YJ                                 SW         \n'...
    '             LY                             WOY          \n'...
    '              OS                         YQJU            \n'...
    '              UDS                      QJNNW             \n'...
    '               NO                  UNQUY  SU             \n'...
    '               SO                         LU             \n'...
    '               YU                         YY             \n'...
    '               ==== My name is PPiAAc ======             \n'...
    ]);
end
function TodaysProverb

Perverbs={'A big fish in a little pond.'
    'A bad workman always blames his tools.'
    'A bird in hand is worth two in the bush.'
    'Absence makes the heart grow fonder.'
    'A cat has nine lives.'
    'A chain is only as strong as its weakest link.'
    'Actions speak louder than words.'
    'A drowning man will clutch at a straw.'
    'Adversity and loss make a man wise.'
    'A fool and his money are soon parted.'
    'A journey of thousand miles begins with a single step.'
    'A leopard doesn''t change its spots.'
    'All good things come to an end.'
    'All''s well that ends well.'
    'All that glitters is not gold.'
    'All''s fair in love and war.'
    'Always put your best foot forward.'
    'Among the blind the one-eyed man is king.'
    'An apple a day keeps the doctor away.'
    'An empty vessel makes much noise.'
    'An idle brain is the devil''s workshop.'
    'An ounce of protection is worth a pound of cure.'
    'A picture is worth a thousand words.'
    'Appearances can be deceptive.'
    'A rolling stone gathers no moss.'
    'A ship in the harbor is safe, but that is not what a ship is for.'
    'A stitch in time saves nine.'
    'As you sow, so you shall reap.'
    'A thing begun is half done.'
    'Barking dogs seldom bite.'
    'Be slow in choosing, but slower in changing.'
    'Beauty is in the eye of the beholder.'
    'Beauty is only skin deep.'
    'Beggars can''t be choosers.'
    'Best things in life are free.'
    'Better late than never.'
    'Better to be poor and healthy rather than rich and sick.'
    'Better to wear out than to rust out.'
    'Blood is thicker than water.'
    'Cleanliness is next to Godliness.'
    'Clothes do not make the man.'
    'Cowards die many times before their deaths.'
    'Cross the stream where it is shallowest.'
    'Curses, like chickens, come home to roost.'
    'Discretion is the better part of valor.'
    'Don''t bite off more than you can chew.'
    'Don''t bite the hand that feeds you.'
    'Don''t blow your own trumpet.'
    'Don''t cast pearls before swine.'
    'Don''t count your chickens before they hatch.'
    'Don''t cross a bridge until you come to it.'
    'Don''t judge a book by its cover.'
    'Don''t kill the goose that lays the golden eggs.'
    'Don''t put all your eggs in one basket.'
    'Don''t put the cart before the horse.'
    'Don''t throw the baby with the bathwater.'
    'Early bird catches the worm.'
    'Easy come, easy go.'
    'Empty bags cannot stand upright.'
    'Every cloud has a silver lining.'
    'Every dog has his day.'
    'Every man is the architect of his destiny.'
    'Every man has his price.'
    'Fall seven times. Stand up eight.'
    'Familiarity breeds contempt.'
    'Fools rush in where angels fear to tread.'
    'Fortune favors the brave.'
    'Get out while the going is good.'
    'Give them an inch and they''ll take a mile.'
    'God helps those who help themselves.'
    'Good things come to those who wait.'
    'Grief divided is made lighter.'
    'Half a loaf is better than none.'
    'Honesty is the best policy.'
    'Hope for the best, prepare for the worst.'
    'If it ain''t broke, don''t fix it.'
    'If the mountain won''t come to Muhammad, Muhammad must go to the mountain.'
    'If wishes were horses, beggars would ride.'
    'If you can''t beat them, join them.'
    'If you play with fire, you''ll get burned.'
    'Ignorance is bliss.'
    'It''s better to be safe than sorry.'
    'It''s easy to be wise after the event.'
    'It''s never too late to mend.'
    'It''s not over till it''s over.'
    'It''s no use crying over spilt milk.'
    'It takes two to make a quarrel.'
    'It takes two to tango.'
    'Laughter is the best medicine.'
    'Learn to walk before you run.'
    'Let sleeping dogs lie.'
    'Life begins at forty.'
    'Lightning never strikes twice in the same place.'
    'Look before you leap.'
    'Make hay while the sun shines.'
    'Money doesn''t grow on trees.'
    'Money talks.'
    'Necessity is the mother of invention.'
    'No news is good news.'
    'Pen is mightier than sword.'
    'People who live in glass houses shouldn''t throw stones at others.'
    'Practice makes perfect.'
    'Rome wasn''t built in a day.'
    'Still waters run deep.'
    'The best-laid plans go astray.'
    'The end justifies the means.'
    'The harder you work, the luckier you get.'
    'The grass is greener on the other side of the fence.'
    'There are more ways than one to skin a cat.'
    'The road to hell is paved with good intentions.'
    'The squeaky wheel gets the grease.'
    'Too many cooks spoil the broth.'
    'Watch the doughnut, and not the hole.'
    'When in Rome, do as the Romans do.'
    'Where there''s smoke there''s fire.'
    'While the cat''s away, the mice will play.'
    'You can lead a horse to water but you can''t make it drink.'
    'You can''t always get what you want.'
    'You can''t have your cake and eat it too.'
    'You can''t teach an old dog new tricks.'
    'You show me the man and I''ll show you the rule.'
    'I can do this all day.' %%Super Rare
    'Don''t waste it. Don''t waste your life.' %%Super Rare
    'What is and always will be my greatest creation, is you.' %%Super Rare
    'Anyone who''s ever going to find his way in this world, has to start by admitting he doesn''t know where the hell he is.' %%Super Rare
    'Take that away, what are you? genius, billionaire, playboy, philanthropist.' %%Super Rare
    'I''m the best.' %%Super Rare
    'If you''re nothing without the suit then you shouldn''t have it.' %%Super Rare
    'Bring me thanos.'%%Super Rare
    'I love you 3000.' %%Super Rare
    'The Jax is not a big "J".' %% Super Super Rare
    'The ultimate question : how is it going?' %% Super Super Rare
    'A busy honey has no time to mourn.' %% Super Super Rare
    'If you find a fire, you should shout that there is a fire.' %% Super Super Rare
    'August 6th! It''s my birthday!' %% Really, it's my birthday.
    'Wa! SANS!'
    'SuperSuperRare!!! Yongjun Cha, who is the Best in QC' %% SSR
    };

numProverb=length(Perverbs);
rng('shuffle');
IndexProverb=randi(numProverb);

fprintf('  Today''s Proverb: %s\n',Perverbs{IndexProverb});
end


function titleName = titleMaker(config)
titleName = sprintf('%s_%.5d',config.h5header,config.run);
end

%% Functions for AutoRun
% function autorunCurveAveragePAL(pauseTime,config) %#ok<*DEFNU>
% % Decapsulate
% h5header=config.h5header;
% parentFolderH5=config.parentFolderH5;
% 
% persistent functionCount
% functionCount=0;
% % errCounter = 0;
% 
% config.folderh5 = '';
% config.weAreEndGameNow = false;
% while true
%     try
%         functionCount=functionCount+1;
%         if 1==functionCount
%             h5StructPre=struct;
%         end
%         [~,~,ListFolderH5,folderH5Last]=getNumChiFileLast(parentFolderH5,h5header);
% 
%         fprintf('  endGameNow | %d\n',config.weAreEndGameNow);
%         fprintf('  folderNow  | %s\n',config.folderh5);
%         fprintf('  folderNew  | %s\n',folderH5Last);
% 
%         if (~config.weAreEndGameNow & strcmp(config.folderh5,folderH5Last)) || isempty(config.folderh5) % normal case
%             config.folderh5=folderH5Last;
%             [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
%         elseif config.weAreEndGameNow & strcmp(config.folderh5,folderH5Last) % Finish, but no next run
%             while strcmp(config.folderh5,folderH5Last)
%                 [~,~,~,folderH5Last]=getNumChiFileLast(parentFolderH5,h5header);
%                 fprintf('  Dormammu! I''ve come to bargain sale!\n');
%                 fprintf('  Waiting for next run\n');
%                 pause(3);
%             end
%             % config.weAreEndGameNow = false;
%         elseif ~config.weAreEndGameNow & ~strcmp(config.folderh5,folderH5Last) % Not finish, but next run
%             fprintf('  We''re not in the endgame yet\n');
%             % minimumFileComment(numH5Last,h5StructLast,ListFolderH5)
%             config.folderh5=sprintf('%s/%s',parentFolderH5,ListFolderH5(end-1).name);
%             [h5StructPre,~,config] =curveAveragePAL(h5StructPre,config);
%             % config.folderh5=sprintf('%s/%s',parentFolderH5,ListFolderH5(end).name);
%             % config.weAreEndGameNow = false;
%         elseif config.weAreEndGameNow & ~strcmp(config.folderh5,folderH5Last) % Finish, but no next run
%             config.folderh5=folderH5Last;
%             [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
%         else
%             fprintf('  Unexpected thing happen at autorun\n');
%             fprintf('  endGameNow | %d\n',config.weAreEndGameNow);
%             fprintf('  folderNow  | %s\n',config.folderh5);
%             fprintf('  folderNew  | %s\n',folderH5Last);
% 
%             error('');
%         end
%     catch
%         fprintf('  Unexpected thing happen at autorun\n');
%     end
%     pauseComment(pauseTime);
% end
%     function [numH5Last,h5StructLast,ListFolderH5,folderH5Last]=getNumChiFileLast(folderh5List,h5header)
%         h5folderPattern=sprintf('%s/%s_00*',folderh5List,h5header);
%         ListFolderH5=dir(h5folderPattern);
% 
%         while 2
%             try
%                 folderH5Last=sprintf('%s/%s',folderh5List,ListFolderH5(end).name);
%                 break;
%             catch
%                 pauseTime=10;
%                 pause(pauseTime)
%                 sprintf('  There is no files ... Pause %d s',pauseTime);
%             end
%         end
%         h5Pattern=sprintf('%s/eh1rayMXAI_int/*.h5',folderH5Last);
%         h5StructLast=dir(h5Pattern);
%         numH5Last=length(h5StructLast);
% 
%         fprintf('\n');
%     end
%     function pauseComment(pauseTime)
%         fprintf('  The next automatic analysis proceeds after %ds.\n',pauseTime);
%         pause(pauseTime)
%     end
%     function minimumFileComment(numChiLast,h5StructLast,ListFolderChi)
%         for chfileIndexndex=1:numChiLast
%             fprintf('  %s\n',h5StructLast(chfileIndexndex).name)
%         end
%         fprintf('\n       The %s has too samll measured data to process...\n',ListFolderChi(end).name);
%         disp('--------------------------------------------------------------')
%     end
% end
function autorunCurveAveragePAL(pauseTime,config) %#ok<*DEFNU>
% Decapsulate
h5header = config.h5header;
parentFolderH5 = config.parentFolderH5;
parentFolderResults = config.parentFolderResults;

persistent functionCount
functionCount=0;
% errCounter = 0;
config.folderh5 = '';
config.weAreEndGameNow = false;
while true
   % try

   functionCount=functionCount+1;
   if 1==functionCount
       h5StructPre=struct;
   end


   if config.forceLatestRun
       % --- force mode: always stick to the latest run until fully caught up ---
       [~,~,~,folderH5Last] = getNumChiFileLast(parentFolderH5, h5header);
       config.folderh5 = folderH5Last;

       % 안내 메시지: 현재 분석 대상
       [~,~,runName] = countChiH5ForRun(config.folderh5, h5header, parentFolderResults);
       fprintf('  [force] Analyzing latest run: %s\n', runName);

       while true
           [chiCount, h5Count, runName] = countChiH5ForRun(config.folderh5, h5header, parentFolderResults);

           if chiCount < h5Count
               % 아직 뒤쳐짐 → 이 run에 계속 머뭄
               fprintf('  [force] %s progress: %d / %d (stay on this run)\n', runName, chiCount, h5Count);
               [h5StructPre,~,config] = curveAveragePAL(h5StructPre, config);
               pause(1);  % 필요 시 조정
           else
               % 현 시점까지는 완료 → 다음 반복에서 최신 run 재동기화
               fprintf('  [force] %s is up-to-date (%d / %d). Re-syncing to the latest run...\n', runName, chiCount, h5Count);
               break
           end
       end
   else
       [~, allOk, folderH5Find] = findDeficientRuns(h5header, parentFolderH5, parentFolderResults);
       config.folderh5 = folderH5Find;
       if allOk
           [chiCount, h5Count, runName] = countChiH5ForRun(config.folderh5, h5header, parentFolderResults);
           fprintf('  [Noforce] %s progress: %d / %d (stay on this run)\n', runName, chiCount, h5Count);
           % fprintf('  Dormammu! I''ve come to bargain sale!\n');
           % fprintf('  Waiting for next run\n');
           pause(3);
       else
           [chiCount, h5Count, runName] = countChiH5ForRun(config.folderh5, h5header, parentFolderResults);
           fprintf('  [Noforce] %s progress: %d / %d (stay on this run)\n', runName, chiCount, h5Count);
           [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
       end
   end
        







       
       % fprintf('  endGameNow | %d\n',config.weAreEndGameNow);
       % fprintf('  folderNow  | %s\n',config.folderh5);
       % fprintf('  folderNew  | %s\n',folderH5Last);
       % if (~config.weAreEndGameNow & strcmp(config.folderh5,folderH5Last)) || isempty(config.folderh5) % normal case
       %     config.folderh5=folderH5Last;
       %     [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
       % elseif config.weAreEndGameNow & strcmp(config.folderh5,folderH5Last) % Finish, but no next run
       %     while strcmp(config.folderh5,folderH5Last)
       %         [~,~,~,folderH5Last]=getNumChiFileLast(parentFolderH5,h5header);
       %         [~, allOk, folderH5Find] = findDeficientRuns(h5header, parentFolderH5, parentFolderResults);
       %         if allOk
       %             fprintf('  Dormammu! I''ve come to bargain sale!\n');
       %             fprintf('  Waiting for next run\n');
       %             pause(3);
       %         else
       %             config.folderh5 = folderH5Find;
       %             [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
       %         end
       %     end
       %     % config.weAreEndGameNow = false;
       % elseif ~config.weAreEndGameNow & ~strcmp(config.folderh5,folderH5Last) % Not finish, but next run
       %     fprintf('  We''re not in the endgame yet\n');
       %     % minimumFileComment(numH5Last,h5StructLast,ListFolderH5)
       %     config.folderh5=sprintf('%s/%s',parentFolderH5,ListFolderH5(end-1).name);
       %     [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
       %     % config.folderh5=sprintf('%s/%s',parentFolderH5,ListFolderH5(end).name);
       %     % config.weAreEndGameNow = false;
       % elseif config.weAreEndGameNow & ~strcmp(config.folderh5,folderH5Last) % Finish, go next run
       % 
       % 
       %     [~, allOk, folderH5Find] = findDeficientRuns(h5header, parentFolderH5, parentFolderResults);
       %     if allOk
       %         config.folderh5=folderH5Last;
       %         [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
       %     else
       %         config.folderh5 = folderH5Find;
       %         [h5StructPre,~,config] = curveAveragePAL(h5StructPre,config);
       %     end
       % else
       %     fprintf('  Unexpected thing happen at autorun\n');
       %     fprintf('  endGameNow | %d\n',config.weAreEndGameNow);
       %     fprintf('  folderNow  | %s\n',config.folderh5);
       %     fprintf('  folderNew  | %s\n',folderH5Last);
       %     error('');
       % end
   % catch
       % fprintf('  Unexpected thing happen at autorun\n');
   % end
   pauseComment(pauseTime);
end
   function [numH5Last,h5StructLast,ListFolderH5,folderH5Last]=getNumChiFileLast(folderh5List,h5header)
       h5folderPattern=sprintf('%s/%s_00*',folderh5List,h5header);
       ListFolderH5=dir(h5folderPattern);
       while 2
           try
               folderH5Last=sprintf('%s/%s',folderh5List,ListFolderH5(end).name);
               break;
           catch
               pauseTime=10;
               pause(pauseTime)
               sprintf('  There is no files ... Pause %d s',pauseTime);
           end
       end
       h5Pattern=sprintf('%s/eh1rayMXAI_int/*.h5',folderH5Last);
       h5StructLast=dir(h5Pattern);
       numH5Last=length(h5StructLast);
       fprintf('\n');
   end
   function pauseComment(pauseTime)
       fprintf('  The next automatic analysis proceeds after %ds.\n',pauseTime);
       pause(pauseTime)
   end
   function minimumFileComment(numChiLast,h5StructLast,ListFolderChi)
       for chfileIndexndex=1:numChiLast
           fprintf('  %s\n',h5StructLast(chfileIndexndex).name)
       end
       fprintf('\n       The %s has too samll measured data to process...\n',ListFolderChi(end).name);
       disp('--------------------------------------------------------------')
   end
end
% function [badRuns, allOk, folderH5Find] = findDeficientRuns(h5header, scanRoot, resultsRoot)
% % findDeficientRuns
% %   h5header   : 예) "250901_I3m_TBA_MeCN"
% %   scanRoot   : 예) "\\192.168.0.254\home\ue_250830_FXL\scan"
% %   resultsRoot: 예) "\\192.168.0.254\home\ue_250830_FXL\scratch\resultsTRXL"
% %
% % 반환:
% %   badRuns : chi_pal 개수가 h5 개수보다 작은 run 이름 배열
% %   allOk   : true/false (모든 run이 정상인지 여부)
% % scanTemp = sprintf('');
% expDir = fullfile(scanRoot, h5header);
% d = dir(fullfile(expDir,"run*")); d = d([d.isdir]);
% runs = string({d.name});
% if isempty(runs)
%    error('findDeficientRuns:NoRunsFound', ...
%          'No run folders found under resultsRoot for h5header "%s". Check h5header or raw data.', ...
%          h5header);
% end
% nums = zeros(numel(runs),1); keep = false(numel(runs),1);
% for i = 1:numel(runs)
%     % i = 14
%    r = runs(i);
%    n = str2double(erase(r,"run"));
%    nums(i) = n;
%    nChi(i) = nfiles(fullfile(expDir, r, "Values"), "*.chi_pal");
%    nH5(i)  = nfiles(fullfile(scanRoot, sprintf("%s_%05d_DIR", h5header, n), "eh1rayMXAI_int"), "*.h5");
%    keep(i) = (nChi(i) < nH5(i));
% end
% [~,ord] = sort(nums);
% badRuns = runs(ord);
% badRuns = badRuns(keep(ord));
% allOk   = isempty(badRuns);
% if isempty(badRuns)
%    folderH5Find = '';
% else
%    minRun = min(nums(keep));
%    folderH5Find = sprintf('%s/%s_%.5d_DIR',scanRoot,h5header,minRun);
% end
% if allOk
%    disp("  All runs are complete")
% else
%    fprintf('  Incomplete runs detected: %d\n',minRun);
%    fprintf('  numH5 | %d numChi | %d\n',nH5(find(keep,1)),nChi(find(keep,1)))
% end
% end
% function n = nfiles(dirpath, pattern)
% if ~isfolder(dirpath)
%    n = 0;
% else
%    n = numel(dir(fullfile(dirpath,pattern)));
% end
% end
% 


function [badRuns, allOk, folderH5Find] = findDeficientRuns(h5header, scanRoot, resultsRoot)
% findDeficientRuns
%   h5header   : 예) "250901_I3m_TBA_MeCN"
%   scanRoot   : 예) "\\192.168.0.254\home\ue_250830_FXL\scan"
%   resultsRoot: 예) "\\192.168.0.254\home\ue_250830_FXL\scratch\resultsTRXL"
%
% 반환:
%   badRuns : chi_pal 개수가 h5 개수보다 작은 run 이름 배열
%   allOk   : true/false (모든 run이 정상인지 여부)
expDir = fullfile(resultsRoot, h5header);

% 기존: results 쪽에서 run 폴더 수집
d = dir(fullfile(expDir,"run*")); d = d([d.isdir]);
runs_results = string({d.name});

%% ADDED: H5 쪽에서도 run 번호 수집 (chi 폴더가 없어도 잡히게)
dh5 = dir(fullfile(scanRoot, sprintf('%s_*_DIR', h5header))); dh5 = dh5([dh5.isdir]);
nums_h5 = [];
if ~isempty(dh5)
    hdrEsc = regexptranslate('escape', h5header);
    nums_h5 = nan(numel(dh5),1);
    for i = 1:numel(dh5)
        t = regexp(dh5(i).name, ['^' hdrEsc '_(\d{5})_DIR$'],'tokens','once');
        if ~isempty(t), nums_h5(i) = str2double(t{1}); end
    end
    nums_h5 = nums_h5(~isnan(nums_h5));
end

%% ADDED: runs 목록을 results와 h5의 합집합으로 생성
nums_res = [];
if ~isempty(runs_results)
    nums_res = arrayfun(@(s) str2double(erase(s,"run")), runs_results);
end
all_nums = unique(sort([nums_res(:); nums_h5(:)]));

% 둘 다 없으면 = 원본 데이터 없음 → 에러
if isempty(all_nums)
   error('findDeficientRuns:NoRunsFound', ...
         'No runs found for h5header "%s" under both resultsRoot and scanRoot. Check h5header or raw data.', ...
         h5header);
end

% 표시용 run 이름 배열로 변환


% 이하 로직은 거의 동일
nums = all_nums(:);                      % 이미 정렬됨
runs = compose("run%03d", nums);  % string 배열 반환
runs = runs(:);                   % 열벡터 보장(옵션)
keep = false(numel(runs),1);
for i = 1:numel(runs)
   r = runs(i);
   n = nums(i);                          % ← 문자열 파싱 대신 숫자 직접 사용

   nChi = nfiles(fullfile(expDir, r, "Values"), "*.chi_pal");   % 폴더 없으면 0
   h5Dir = fullfile(scanRoot, sprintf("%s_%05d_DIR", h5header, n), "eh1rayMXAI_int");
   nH5  = nfiles(h5Dir, "*.h5");                                  % 폴더 없으면 0

   keep(i) = (nChi < nH5);
end

% 이미 nums가 오름차순이므로 그대로 필터
badRuns = runs(keep);
allOk   = isempty(badRuns);

% folderH5Find: 가장 번호가 작은 미완 run의 H5 폴더
if allOk
   folderH5Find = '';
   disp("  All runs are complete")
else
   firstBadNum  = nums(find(keep,1,'first'));  % 가장 작은 번호의 미완
   folderH5Find = fullfile(scanRoot, sprintf('%s_%05d_DIR', h5header, firstBadNum));
   disp("  Incomplete runs detected:")
   disp(badRuns)
end
end

function n = nfiles(dirpath, pattern)
if ~isfolder(dirpath)
   n = 0;
else
   n = numel(dir(fullfile(dirpath,pattern)));
end
end







function [nChi, nH5, runName] = countChiH5ForRun(h5RunDir, h5header, resultsRoot)
[~, dirName] = fileparts(char(h5RunDir));
tok = regexp(dirName, ['^' regexptranslate('escape',h5header) '_(\d{5})_DIR$'], 'tokens','once');
if isempty(tok), nChi=0; nH5=0; runName=""; return; end
runNum  = str2double(tok{1});
runName = sprintf('run%03d', runNum);
chiDir = fullfile(resultsRoot, h5header, runName, 'Values');
h5Dir  = fullfile(h5RunDir, 'eh1rayMXAI_int');
nChi = nfiles2(chiDir, '*.chi_pal');
nH5  = nfiles2(h5Dir,  '*.h5');
end

function n = nfiles2(p, pat)
if ~isfolder(p), n = 0; else, n = numel(dir(fullfile(p,pat))); end
end



%% Telomere
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum
% LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsum % LolemIpsu

