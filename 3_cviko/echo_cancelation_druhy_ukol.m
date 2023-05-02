clear all
close all
clc

%%
load('echocancelation.mat')


%WW = miso_firwiener(L-1,x',d');
L = 100;

W = miso_firwiener(L-1,ref',mic');

e = filter(W,1,ref);

sound(mic - e,fs)

% e je chyba -> v našem případě jde o mužský hlas. Takže je potřeba odečíst
% to co je na microfónu s filtrem.
% větší L -> lepší filter .. vždycky!
% neplatí při odhadování bez referenčního signálu