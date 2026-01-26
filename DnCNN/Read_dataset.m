%% KEK2512 I3- in MeOH run4


%% Read Datasets
clc; clear; close all;

filename = "C:\Users\user\Desktop\KEK2512_MeOH_h5\run4.h5";

data_path_n2n   = "/delay/delay_m1.00e-09/neg";     %373
data_path_n1n   = "/delay/delay_m1.00e-09/on";      %373
data_path_p150p = "/delay/delay_p1.50e-10/on";      %373
data_path_p1n   = "/delay/delay_p1.00e-09/on";      %374

tic
D_n2n   = h5read(filename,data_path_n2n);
D_n1n   = h5read(filename,data_path_n1n);
D_p150p = h5read(filename,data_path_p150p);
D_p1n   = h5read(filename,data_path_p1n);
toc
%% Reshape

D_n2n   = double(permute(D_n2n,[3 2 1]));
D_n1n   = double(permute(D_n1n,[3 2 1]));
D_p150p = double(permute(D_p150p,[3 2 1]));
D_p1n   = double(permute(D_p1n,[3 2 1]));

%% SUM

D = [D_n2n; D_n1n; D_p150p; D_p1n(1:373,:,:)];