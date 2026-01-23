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
data_path_neg = "/delay/delay_m1.00e-09/neg";
data_path_on = "/delay/delay_m1.00e-09/on";

raw_data_neg = h5read(filename,data_path_neg);
raw_data_on = h5read(filename,data_path_on);
Ons = permute(raw_data_on,[3 2 1]);
Ons = double(Ons);
Negs = permute(raw_data_neg, [3 2 1]);
Negs = double(Negs);

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

for i = 1 : 373
    res = ai.integrate1d(Negs(:,:,i), int32(1024), ...
    pyargs('unit', 'q_A^-1', 'polarization_factor', 0.99, 'mask', mask));
    I(:,i) = double(res{2}.tolist())';
end

q = double(res{1}.tolist())';

%% PLOT
close all;

figure('Color', 'w');
set(gca, 'ColorOrder', jet(373), 'NextPlot', 'replacechildren');
plot(q, I, 'LineWidth', 1.5);
xlabel('q (A^{-1})');
ylabel('Intensity');
% colormap(jet);
colorbar;
grid on;
axis tight;

%% AUC NORMALIZATION

norm_range = [4.0; 7.8];
norm_idx = find(q>=norm_range(1)&q<=norm_range(2));

norm_factors = sum(I(norm_idx,:),1);
I_norm = I ./ norm_factors;

%% ON IMAGE NORMALIZE
norm_factors_img = reshape(norm_factors,1,1,[]);
On_normalized = Ons ./ norm_factors_img;

%% HISTGRAM OF PIXEL

for i = 1024 : 1024
    for j = 1 : 50
        figure;     
        B = On_normalized(i,j*20,:);
        histogram(B,100);
        hold on;
        xline(mean(B),'--r','LineWidth',2);
        hold off;
    end
end