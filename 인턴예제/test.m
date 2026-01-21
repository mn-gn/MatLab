clc; clear; close all;

M = readmatrix('Example#5_Au3_deltaS_time_series.dat');
Q = M(2:end,1);
T = M(1,2:end)';
A = M(2:end,3:end); %[3 20 40 60 80 100 ]);

[U, S, V] = svd(A);