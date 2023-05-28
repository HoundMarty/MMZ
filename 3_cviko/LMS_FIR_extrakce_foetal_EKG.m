clear all
close all
clc

%%
% V souboru foetalECG.mat je záznam EKG těhotné ženy. 
% Signály 2 až 6 jsou z abdominální oblasti, zatímco signály 7 až 9 jsou z oblasti hrudní.
% Pomocí LMS najděte filtr, kterým z vybraného EKG kanálu z abdominální 
% oblasti odečtete vybraný signál z oblasti hrudní. 
% Výsledkem by měl být vyčištěný záznam fetálního EKG.

%%
load('C:\Users\98653\Documents\Skola\2_rocnik\LS\NMZ\NMZ\3_cviko\foetalECG.mat')

% abdominal = foetal_ecg(:,2:6).';
% hrudni = foetal_ecg(:,7:9).';

x = foetal_ecg.';

eegplot(x);

N = 100;
L = N + 1;
% h = miso_firwiener(N, x(6:8, :).', x(1, :).');
% h = reshape(h, L, []);

h = miso_firwiener(N, x(7:9,:).', x(2,:).');
h = reshape(h, L, []);

mama = filter(h(:, 1), 1, x(7, :)) + filter(h(:, 2), 1, x(8, :)) + filter(h(:, 3), 1, x(9, :));

fetus = x(2) - mama;

eegplot([fetus; x(2, :)]);
