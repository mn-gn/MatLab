clc; clear all;

T = readmatrix("kinetic_analysis_example1.dat");
A = T(2:end,2:end);
[U, S, V] = svd(A);
time = T(1,2:end)';