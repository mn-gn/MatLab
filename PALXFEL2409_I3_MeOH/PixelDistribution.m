clc; clear; close all;


folder = "Z:\PAL_XFEL_2409\rawData\I3_MeOH_00002_DIR\eh1rayLX_img";
h5dir = dir(sprintf("%s/*.h5",folder));

% for 1~size

h5name = h5dir(1).name;
filename = sprintf("%s/%s",folder,h5name);
h5 = h5info(filename);

I = zeros([size(h5.Datasets,1),h5.Datasets(1).Dataspace.Size]);
for i = 1 : size(h5.Datasets,1)
    shot_name = h5.Datasets(i).Name;

    I(i,:,:) = h5read(filename,sprintf("/%s",shot_name));
end

%% PLOT I (BEFORE PROCESSING)
I_mean = squeeze(mean(I,1));
I_std = squeeze(std(I,0,1));
close all 
for i = 1 : 20
    figure; hold on;
    title(sprintf("(x=424, y=%d) \\sigma=%f",8*i-7,I_std(8*i-7,424)));
    histogram(I(:,8*i-7,424),100);
    xline(I_mean(8*i-7,424),'r:','LineWidth',3);
    hold off;
end
% figure;
% surf(tth,1:1:500,I');
% shading flat
% xlabel("2\theta"); ylabel("index"); zlabel("intensity");
% colormap turbo

%% AUC NORMAILZATION

norm_range = [18.9; 37.55];
norm_idx = find(tth>=norm_range(1)&tth<=norm_range(2));

norm_factors = sum(I(norm_idx,:),1);
I_norm = I ./ norm_factors;

%% PLOT I (AUC NORMALIZED)
close all
I_norm_mean = mean(I_norm,2);
I_norm_std = std(I_norm,0,2);
for i = 1 : 20
    figure; hold on;
    title(sprintf("2\\theta=%f, \\sigma=%f",tth(20*i),I_norm_std(20*i)));
    histogram(I_norm(20*i,:),100);
    xline(I_norm_mean(20*i),'r:','LineWidth',3);
    hold off;
end
figure;
surf(tth,1:1:500,I_norm');
shading flat
xlabel("2\theta"); ylabel("index"); zlabel("intensity");
colormap turbo