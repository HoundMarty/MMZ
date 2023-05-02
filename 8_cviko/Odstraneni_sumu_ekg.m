clear *
close all
clc

%% Odstranění šumu z EKG

% Odstraňte síťový šum ze signálu v souboru EKG3channels_sinus.mat pomocí libovolné ICA metody.

% EKG3channels_sinus.mat EKG3channels_sinus.mat

%%

load("EKG3channels_sinus (1).mat")

% x = a*s

W = efica(x,eye(3));

Y = W*x;

eegplot(Y)
% odstraneni nezadouci
Y(1,:) = 0;

x_new = W\Y;

eegplot(x_new)

%% to samé pro foetal data

load("foetalECG.mat")

x = foetal_ecg';

W = efica(x,eye(9));

Y = W*x;

% eegplot(Y)
% odstraneni nezadouci 
Y([1 2 3 4 5],:) = 0;

x_new = W\Y;

eegplot(x_new)

