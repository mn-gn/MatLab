   %[control:button:589c]{"position":[1,2]} %[control:button:3434]{"position":[2,3]}
clc; clear; close all;
%%
%[text] ## PYTHON LIBRARY LOAD
try %[output:group:0891a3c7]
    fabio = py.importlib.import_module('fabio');
    detectors = py.importlib.import_module('pyFAI.detectors');
    azimuthal = py.importlib.import_module('pyFAI.azimuthalIntegrator');
    disp('Python libraries loaded.'); %[output:2555a325]
catch
    error('파이썬 라이브러리 로드 실패. pyFAI/fabio 설치 확인 필요.');
end %[output:group:0891a3c7]
%%
%[text] ## DATA PATH SETTING
folder = "Z:\KEK_2512_h5\results\KEK_2512"; %[control:filebrowser:1b33]{"position":[10,43]}
SampleList=folder_list(folder);
sample = SampleList(2); %[control:dropdown:37c1]{"position":[10,23]}
path = sprintf('%s/%s',folder,sample);
RunList=folder_list(path);
runname =  RunList(8); %[control:dropdown:9a62]{"position":[12,22]}
td_path = "/delayinfo/delay";
h5_path = sprintf("%s/%s/raw/images/image_pairs.h5",path,runname);
td_str = h5read(h5_path,td_path);
td = td_str(3); %[control:dropdown:1333]{"position":[6,15]}
%%
%[text] ## DATA LOAD
%[text] h5 for one time delay, poni and mask
Data = struct;
[Data.on, Data.off] = KEK_h5_load(h5_path,td);
poni = "C:\Users\user\Desktop\MNGN\DnCNN\Detector_poni.poni"; %[control:filebrowser:6402]{"position":[8,61]}
mask = read_mask("C:\Users\user\Desktop\MNGN\DnCNN\Mask_run4.edf"); % need read_mask function file %[control:filebrowser:059c]{"position":[18,66]}
%%
%[text] ## AI SETTING
%[text] unit, polarization factor, num\_q and mask
ai = AI(poni,detectors,azimuthal);
num_q = 1024;
ai_option = pyargs('unit', 'q_A^-1','polarization_factor', 0.99,'mask', mask);
%%
%[text] ## AI START
%[text] get azimuthal integrated on and off data + q list [ai\_start](internal:M_2014)
Data = ai_start(ai,Data,ai_option,num_q); % make Data.q, Data.ai.on and Data.ai.off
%%
%[text] ## AUC NORMALIZATION
%[text] get auc normalized on and off data for given range
norm_range = [4 6]; %[control:rangeslider:81c1]{"position":[14,19]}
Data = auc_norm(Data,norm_range);
%%
%[text] ## AUC CUT
%[text] get good img
% cut_range = [2 6]; %[control:rangeslider:8d27]{"position":[15,20]}
% Data = auc_cut(Data,cut_range);
%%
%[text] ## FIGURE
%[text] plot $\\Delta S\\left(q\\right)${"editStyle":"visual"} for selected run and time delay
figure;
plot(Data.q,Data.mean); %[output:65d78a6d]
title([sample runname td],Interpreter="none");
xlabel('q(A^{-1})');
ylabel('\DeltaS(q)');
grid on;
%%
%[text] ## FUNCTIONS
%[text] - **`folder_list(foldername)`** returns subfolder list of the given foldername. \
function lst = folder_list(foldername)
    arguments
        foldername (1,1) string {mustBeFolder}
    end
    d = dir(foldername);
    lst = string({d([d.isdir]&~ismember({d.name},{'.','..'})).name})';
end
%[text] - **`KEK_h5_load(filename,td)`** returns on and off img datasets of **h5** at time delay with td. \
function [on, off] = KEK_h5_load(filename,td)
    arguments
        filename (1,1) string {mustBeFile} 
        td (1,1) string
    end
    if ~endsWith(filename, ".h5", 'IgnoreCase', true)
        error("filename must be **.h5** format!");
    end
    path_on  = sprintf("/delays/delay_%s/on" ,td);
    path_off = sprintf("/delays/delay_%s/off",td);
    
    on  = h5read(filename,path_on );
    off = h5read(filename,path_off);
    
    on  = double(permute(on,  [1 3 2]));
    off = double(permute(off, [1 3 2]));
end
%[text] - **`AI(poni_path,detectors,azimuthal)`** returns `azimuthal.AzimuthalIntegrator` with given parameters. \
function ai = AI(poni_path,detectors,azimuthal)
    arguments
        poni_path (1,1) string {mustBeFile}
        detectors (1,1) py.module
        azimuthal (1,1) py.module

    end

    poni = read_poni(poni_path); % need read_poni function file
    
    pixel_size1 = poni.PixelSize1;
    pixel_size2 = poni.PixelSize2;
    dist        = poni.Distance;
    poni1       = poni.Poni1;
    poni2       = poni.Poni2;
    rot1        = poni.Rot1;
    rot2        = poni.Rot2;
    rot3        = poni.Rot3;
    wavelength  = poni.Wavelength;
    
    det_shape = py.tuple({int32(2048), int32(2048)});
    detector = detectors.Detector(pixel_size1, pixel_size2, pyargs('max_shape', det_shape));
    
    ai_params = pyargs(...
        'dist', dist, 'poni1', poni1, 'poni2', poni2, ...
        'rot1', rot1, 'rot2', rot2, 'rot3', rot3, ...
        'detector', detector, 'wavelength', wavelength);
    
    ai = azimuthal.AzimuthalIntegrator(ai_params);
end
%[text] - **`ai_start(ai,Data,ai_option,num_q)`** returns azimuthal integrated on and off datasets. also q vector. \
%[text]  `Data.ai.on, Data.ai.off` and `Data.q`
function Data = ai_start(ai,Data,ai_option,num_q) %[text:anchor:M_2014]
arguments
    ai (1,1) py.pyFAI.integrator.azimuthal.AzimuthalIntegrator
    Data (1,1) struct
    ai_option (1,1) pyargs
    num_q (1,1) double {mustBePositive} = double(1024)
end
    num_shot = size(Data.on,1);
    Data.ai.on = cell(num_shot, 1);
    Data.ai.off = cell(num_shot, 1);
    for i = 1 : num_shot
        res_on = ai.integrate1d(squeeze(Data.on(i,:,:)), int32(num_q), ai_option);
        Data.ai.on{i} = double(res_on{2}.tolist())';
    
        res_off = ai.integrate1d(squeeze(Data.off(i,:,:)), int32(num_q), ai_option);
        Data.ai.off{i} = double(res_off{2}.tolist())';
    end
    Data.q = double(res_on{1}.tolist())';
end
%[text] - **`auc_norm(Data,norm_rage)`** normalize the Data with norm\_range and get difference of normalized datasets. \
%[text]  `Data.normalized.on, Data.normalizedoff` and `Data.diff`
function Data = auc_norm(Data,norm_range)
arguments
    Data (1,1) struct
    norm_range (1,2) double {mustBePositive} = [4 6]
end
norm_idx = find(Data.q>=norm_range(1)&Data.q<=norm_range(2));
M_on  = [Data.ai.on{ :}];
M_off = [Data.ai.off{:}];

Data.norm_factor_on  = sum(M_on( norm_idx,:),1);
Data.norm_factor_off = sum(M_off(norm_idx,:),1);

Data.normalized.on  = M_on  ./ Data.norm_factor_on;
Data.normalized.off = M_off ./ Data.norm_factor_off;

Data.diff = Data.normalized.on - Data.normalized.off;
end
%[text] - **`auc_cut(Data,cut_range)`** user cuts bad imgs. get good $\\Delta S\\left(q,t\_{\\textrm{given}} \\right)${"editStyle":"visual"}. \
%[text]  `Data.good_DS, Data.mean`
function Data = auc_cut(Data,cut_range)
arguments
    Data (1,1) struct
    cut_range (1,2) double {mustBePositive} = [2 6]
end

cut_idx = find(Data.q>=cut_range(1)&Data.q<=cut_range(2));

abs_DS = abs(Data.diff);
A = sum(abs_DS(cut_idx,:),1);

figure;
semilogy(1:length(A),A,'o');
title('click cut threshold');
[~, y_click] = ginput(1);
close;
good_range = [0 y_click];
disp(['cut range: ', num2str(0),' ~ ',num2str(y_click)]);
good_idx = find(A>=good_range(1)&A<=good_range(2));
Data.good_DS = Data.diff(:,good_idx);
Data.mean = mean(Data.good_DS,2);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":38.7}
%---
%[control:button:589c]
%   data: {"label":"초기화","run":"Section"}
%---
%[control:button:3434]
%   data: {"label":"전체 실행","run":"AllSections"}
%---
%[control:filebrowser:1b33]
%   data: {"browserType":"Folder","defaultValue":"\"\"","label":"Folder","run":"Section"}
%---
%[control:dropdown:37c1]
%   data: {"defaultValue":"SampleList(1)","itemLabels":["I3_MeCN","I3_MeOH","I3_water"],"items":["SampleList(1)","SampleList(2)","SampleList(3)"],"itemsVariable":"SampleList","label":"sample","run":"Section"}
%---
%[control:dropdown:9a62]
%   data: {"defaultValue":"RunList(1)","itemLabels":["run1","run10","run11","run12","run13","run2","run3","run4","run5","run6","run7","run8","run9"],"items":["RunList(1)","RunList(2)","RunList(3)","RunList(4)","RunList(5)","RunList(6)","RunList(7)","RunList(8)","RunList(9)","RunList(10)","RunList(11)","RunList(12)","RunList(13)"],"itemsVariable":"RunList","label":"run number","run":"Section"}
%---
%[control:dropdown:1333]
%   data: {"defaultValue":"td_str(1)","itemLabels":["-1.00e-09","1.50e-10","1.00e-09"],"items":["td_str(1)","td_str(2)","td_str(3)"],"itemsVariable":"td_str","label":"time delay","run":"Section"}
%---
%[control:filebrowser:6402]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"*.poni","run":"Section"}
%---
%[control:filebrowser:059c]
%   data: {"browserType":"File","defaultValue":"\"\"","label":"*.edf","run":"Section"}
%---
%[control:rangeslider:81c1]
%   data: {"defaultValue":"[4 6]","label":"norm_range","max":8,"min":2,"run":"Section","runOn":"ValueChanging","step":0.1}
%---
%[control:rangeslider:8d27]
%   data: {"defaultValue":"[2 6]","label":"cut_range","max":8,"min":0,"run":"Section","runOn":"ValueChanging","step":0.1}
%---
%[output:2555a325]
%   data: {"dataType":"text","outputData":{"text":"Python libraries loaded.\n","truncated":false}}
%---
%[output:65d78a6d]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized field name \"mean\"."}}
%---
