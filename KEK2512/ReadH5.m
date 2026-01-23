clc; clear; close all;

% A = zeros(13,2048,2048);

filename = "Z:\202512KEK_h5\run4.h5";
data_path_neg = "/delay/delay_p1.50e-10/neg";
data_path_on = "/delay/delay_p1.50e-10/neg";

tic
T_neg = h5read(filename, data_path_neg);
T_on = h5read(filename,data_path_on);
toc

T_neg = double(T_neg);
T_on = double(T_on);

reference = squeeze(mean(T_neg,1));

%% MEAN PLOT
close all;
figure;

imagesc(squeeze(T_neg(100,:,:)));
colormap('gray');
axis image

for i = 1024 : 1024
    for j = 1 : 50
        figure;     
        A = reference(i,j*20);
        B = mean(T_on(:,i-2:i+2,j*20-2:j*20+2),[2,3]);
        C = B-A;

        histogram(B,100);
        hold on;
        xline(mean(B),'--r','LineWidth',2);
        hold off;
    end
end

%% PIXEL PLOT
close all;
figure;

imagesc(squeeze(T_neg(100,:,:)));
colormap('gray');
axis image

for i = 1024 : 1024
    for j = 1 : 50
        figure;     
        B = T_on(:,i,j*20);
        histogram(B,100);
        hold on;
        xline(mean(B),'--r','LineWidth',2);
        hold off;
    end
end