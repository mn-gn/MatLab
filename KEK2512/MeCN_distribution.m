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
filename = "Z:\KEK_2512_h5\I3_MeCN\run5.h5";
data_path_off = "/delays/1.50e-10/off";
data_path     = "/delays/1.50e-10/on";

D_off = h5read(filename, data_path_off);
D_on     = h5read(filename,data_path);

D_off = double(permute(D_off,[3 2 1]));
D_on     = double(permute(D_on, [3 2 1]));

poni = read_poni("Detector_poni.poni");
mask = read_mask("MeCN_run5.edf");

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
I_off = zeros(1024,334);
I_on     = zeros(1024, 334);

for i = 1 : size(I_off,2)
    res1 = ai.integrate1d(D_off(:,:,i), int32(1024), ...
    pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask, 'method', 'median'));
    I_off(:,i) = double(res1{2}.tolist())';

    res2 = ai.integrate1d(D_on(:,:,i), int32(1024), ...
    pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask, 'method', 'median'));
    I_on(:,i) = double(res2{2}.tolist())';
end



%% AUC NORMALIZATION

norm_range = [4.0; 7.8];
norm_idx = find(q>=norm_range(1)&q<=norm_range(2));

norm_factors_off = sum(I_off(norm_idx,:),1);
norm_factors_on = sum(I_on(norm_idx,:),1);

I_off_norm = I_off ./ norm_factors_off;
I_on_norm = I_on ./ norm_factors_on;

%% REF IMAGE NORMALIZE
norm_factors_img_off = reshape(norm_factors_off,1,1,[]);
normalized_off = D_off ./ norm_factors_img_off;

norm_factors_img_on = reshape(norm_factors_on,1,1,[]);
normalized_on = D_on ./ norm_factors_img_on;

mean_off = mean(normalized_off,3);
mean_on = mean(normalized_on,3);

%% DISTRIBUTION
close all;
for i = 80 : 100
    figure; hold on;
    histogram(mean(normalized_on(10*i-1:10*i+1,1023:1025,:),[1 2]),100);
    % xline(M(20*i),'r:','LineWidth',2);
    hold off
end