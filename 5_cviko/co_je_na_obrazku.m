clear *
close all
clc

addpath("l1magic-1.11\l1magic\Optimization\");

%%

% Rekonstruujte obrázek o velikosti 65x97 bodů z jeho komprimovaného
% záznamu y v souboru ukol2.mat. Vzorkující transformace je dána maticí A. 
% Z numerického formátu uint8 převeďte A do formátu double.

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