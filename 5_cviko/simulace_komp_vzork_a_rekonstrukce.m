clear *
close all
clc

N = 1000;
S = 10;
C = randi([S,N],1);
m = 10*log10(1000)*C;

x = zeros(N,1);
p = randperm(N);
x(p(1:S)) = randn(S,1);

A = orth(randn(N,m));


%% je dobře?
figure;plot(x);
y = A'*x;
figure;plot(y);
xx = A*y;
figure;plot(xx);


%% OMP
% y = U'*x;
% k = 0.1;

[x_rec,r,normR,residHist,errHist]  = OMP(A',y,S);

check = x_rec - x; 

figure;plot(x_rec);

%% l1magic
