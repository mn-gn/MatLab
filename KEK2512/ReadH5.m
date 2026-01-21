clc; clear; close all;

% A = zeros(13,2048,2048);

filename = "Z:\202512KEK_h5\run4.h5";
data_path = "/delay/delay_p1.50e-10/on";


T = h5read(filename, data_path);
T = double(T);
%%


ff;

A = squeeze(mean(T,1));

%%
close all;
figure;

imagesc(squeeze(T(100,:,:)));
colormap('gray');
axis image

for i = 1024 : 1024
    for j = 1 : 50
        figure;     
        A = T(:,i,j*20);
        histogram(A,100);
        hold on;
        xline(mean(A),'--r','LineWidth',2);
        hold off;
    end
end