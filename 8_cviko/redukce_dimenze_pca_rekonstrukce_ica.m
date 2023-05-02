clear *
close all
clc

%% Redukce dimenze dat pomocí PCA, aplikace ICA a rekonstrukce

% Redukujte dimenzi záznamu v souboru EKG_KES.mat pomocí PCA.
% Data zpracujte pomocí ICA algoritmu (vymažte nezávislé komponenty odpovídající šumu).
% Data rekonstruujte.
% Porovnejte s výsledkem, kdyby se použila pouze ICA.
% EKG_KES.mat EKG_KES.mat24. dubna 2023, 10.25

load("EKG_KES.mat")

C = x*x'/length(x);
[V D] = eig(C);

% plot(log(diag(D))); % 1. - 4. nejmenší vlastní komponenta je mrzká -> nulový energie
% v reálnejch datech to takto nevyjde

Z = V'*x;

%normalizace
Z = Z.*(1./sqrt(diag(D)));

% eegplot(Z)

% redukce dimenze
Z_red = V(:,5:end)'*x;

% ICA
% eye(8) --> je 5:end
W = efica(Z_red,eye(8));

%vyhození některých signálů
Y = W*Z_red;

Y([1 2 3]) = 0;

%rekonstrukce dat původních dat
% pozor máme obdélníkovou matici!
rek_x = V(:,5:12) * Z_red;

eegplot(rek_x)