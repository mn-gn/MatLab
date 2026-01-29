clc; clear; close all;

%% PYTHON LIBRARY LOAD
try
    fabio = py.importlib.import_module('fabio');
    detectors = py.importlib.import_module('pyFAI.detectors');
    azimuthal = py.importlib.import_module('pyFAI.azimuthalIntegrator');
    disp('Python libraries loaded.');
catch
    error('파이썬 라이브러리 로드 실패. pyFAI/fabio 설치 확인 필요.');
end

%% DATA LOAD
filename = "Z:\202512KEK_MeOH_h5\run4.h5";
data_path_p150p_ref = "/delay/delay_p1.50e-10/neg";
data_path_p150p     = "/delay/delay_p1.50e-10/on";

D_p150p_ref = h5read(filename, data_path_p150p_ref);
D_p150p     = h5read(filename,data_path_p150p);

D_p150p_ref = double(permute(D_p150p_ref,[3 2 1]));
D_p150p     = double(permute(D_p150p, [3 2 1]));

poni = read_poni("Detector_poni.poni");
mask = read_mask("Mask_run4.edf");

%% GEOMETRY

pixel_size1 = poni.PixelSize1;
pixel_size2 = poni.PixelSize2;
dist = poni.Distance;
poni1 = poni.Poni1;
poni2 = poni.Poni2;
rot1 = poni.Rot1;
rot2 = poni.Rot2;
rot3 = poni.Rot3;
wavelength = poni.Wavelength;

%% AI ENGINE START
det_shape = py.tuple({int32(2048), int32(2048)});
detector = detectors.Detector(pixel_size1, pixel_size2, pyargs('max_shape', det_shape));

ai_params = pyargs(...
    'dist', dist, 'poni1', poni1, 'poni2', poni2, ...
    'rot1', rot1, 'rot2', rot2, 'rot3', rot3, ...
    'detector', detector, 'wavelength', wavelength);

ai = azimuthal.AzimuthalIntegrator(ai_params);

%% INTEGRATE
I_p150p_ref = zeros(1024,373);
I_p150p     = zeros(1024, 373);

for i = 1 : size(I_p150p_ref,2)
    res1 = ai.integrate1d(D_p150p_ref(:,:,i), int32(1024), ...
    pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask, 'method', 'median'));
    I_p150p_ref(:,i) = double(res1{2}.tolist())';

    res2 = ai.integrate1d(D_p150p(:,:,i), int32(1024), ...
    pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask, 'method', 'median'));
    I_p150p(:,i) = double(res2{2}.tolist())';
end

DS_p150p = I_p150p - I_p150p_ref;

q = double(res1{1}.tolist())';

M = mean(DS_p150p,2);

%% DISTRIBUTION
close all;
for i = 4 : 47
    figure; hold on;
    histogram(DS_p150p(20*i,:),373);
    xline(M(20*i),'r:','LineWidth',2);
    hold off
end