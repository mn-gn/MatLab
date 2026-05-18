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

filename        = "C:\Users\user\Desktop\KEK2512_MeOH_h5\run4.h5";

data_path_ref   = "/delay/delay_m1.00e-09/neg";                      %373
data_path_n1n   = "/delay/delay_m1.00e-09/on";                       %373
data_path_p150p = "/delay/delay_p1.50e-10/on";                       %373
data_path_p1n   = "/delay/delay_p1.00e-09/on";                       %374


D_ref   = h5read(filename,data_path_ref);
D_n1n   = h5read(filename,data_path_n1n);
D_p150p = h5read(filename,data_path_p150p);
D_p1n   = h5read(filename,data_path_p1n);

D_ref   = double(permute(D_ref,[3 2 1]));
D_n1n   = double(permute(D_n1n,[3 2 1]));
D_p150p = double(permute(D_p150p,[3 2 1]));
D_p1n   = double(permute(D_p1n,[3 2 1]));

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
I = zeros(1024,373);

for i = 1 : size(I,2)
    res = ai.integrate1d(D_ref(:,:,i), int32(1024), ...
    pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask));
    I(:,i) = double(res{2}.tolist())';
end

q = double(res{1}.tolist())';

%% PLOT
close all;

figure('Color', 'w');
set(gca, 'ColorOrder', turbo(373), 'NextPlot', 'replacechildren');
plot(q, I, 'LineWidth', 1.5);
xlabel('q (A^{-1})');
ylabel('Intensity');
colormap("turbo");
colorbar;
grid on;
axis tight;

%% AUC NORMALIZATION

norm_range = [4.0; 7.8];
norm_idx = find(q>=norm_range(1)&q<=norm_range(2));

norm_factors = sum(I(norm_idx,:),1);
I_norm = I ./ norm_factors;

%% REF IMAGE NORMALIZE
norm_factors_img = reshape(norm_factors,1,1,[]);
normalized_ref = D_ref ./ norm_factors_img;

mean_ref = mean(normalized_ref,3) * 100;

%% DETECT HOT PIXEL
localMedian = medfilt2(mean_ref, [3 3]);
localStd = stdfilt(mean_ref,true(3));
zScoreMap = abs(mean_ref - localStd) ./ (localStd + eps);
sigmaThreshold = 29;

hotPixelMask = zScoreMap > sigmaThreshold;


%% HOT PIXEL FIGURE
figure; 
subplot(1,2,1); imagesc(mean_ref); axis image; colorbar;
title('Original Image');

subplot(1,2,2); imagesc(hotPixelMask); axis image; colorbar;
title(['Detected Hot Pixels (Sigma > ' num2str(sigmaThreshold) ')']);

colormap hot;
% 핫픽셀 개수 출력
fprintf('탐지된 핫픽셀 개수: %d\n', sum(hotPixelMask(:)));

%%

NewMask = hotPixelMask & ~mask;