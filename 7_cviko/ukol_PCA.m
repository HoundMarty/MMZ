clear *
close all
clc

%% PCA

% s = randn(1,10000);
% a = randn(5,1);
% 
% v = randn(5,10000);
% 
% % cov(v');
% 
% x = a*s + 0.01*v;
% 
% C_x = x*x'/length(x);
% 
% [V, D] = eig(C_x);

%% Úkol: PCA signálů z EKG
% Proveďte Analýzu hlavních a nezávislých komponent na signály v souborech.
% Rekonstruujte data z vybraných komponent.

load("EKG3channels_sinus (1).mat")

C = x*x'/length(x);
[V, D] = eig(C); 
% V = V(:,end-num_reduced_dim+1:end);
Z_components = V'*x;

% x = inv(V')*Z_components

eegplot(Z_components);

%% pokračování 

load("foetalECG.mat")

% eegplot(foetal_ecg(:,2:9)');
x = foetal_ecg(:,2:9)';
C = x*x'/length(x);
[V, D] = eig(C); 
% V = V(:,end-3+1:end);
Z_components = V'*x;

Z_components(5:end,:) = 0;

eegplot(Z_components)

