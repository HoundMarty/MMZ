clear all
close all
clc

%% zadání
% V souboru echocancelation.mat je signál mic, který simuluje záznam mužského hlasu na mikrofonu rušený jiným hlasem (ženským). V proměnné ref je referenční signál ženského hlasu. Použijte LMS filtr k potlačení ženského hlasu v záznamu mužského hlasu.

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