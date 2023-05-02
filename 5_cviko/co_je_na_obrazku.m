clear *
close all
clc

addpath("l1magic-1.11\l1magic\Optimization\");
%% 
load('ukol2.mat')

A = double(A);

[x_rec,r,normR,residHist,errHist]  = OMP(A,y,3000);

imshow(reshape(x_rec,65,97));

%%
% l1eq_pd
x0 = zeros(6305,1);
At = 3000; %??
pdtol = 0.1; %???
cgtol = 0.1; %???
xp = l1eq_pd(x0, A, At, y);%, pdtol, pdmaxiter, cgtol, cgmaxiter);

imshow(reshape(xp,65,97));