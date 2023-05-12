clear *
close all
clc

%% Zadání

% Odstranění harmonických složek pomocí filtru

%% Jak na to?

% První možnost je, že signál transformujete pomocí FFT, vynulujete nežádoucí 
% složky (symetricky tj. aby bylo spektrum symetrické) a signál rekonstruujete pomocí inverzní FFT. 
% Druhou možností je, že navrhnete NOTCH filtry pro vybrané harmonické 
% složky a signál jimi přefiltrujete.

%% Řešení **Vytvořil** Martin Motejlek

load("EKG3channels_sinus.mat");
x = x(1, :); % jen jeden kanál
fs = 500; %bylo zadáné

X = fft(x);

f = (0:length(X) - 1) / length(X) * fs;
plot(f, abs(X))
title("Originál FFT");

%% Filtrace
% Nejsem si jistý, jestli je vhodně zvolený bandwidth
y = x;

[b, a] = iirnotch(36.8 / fs * 2, 5 / fs * 2);
y = filter(b, a, y);

[b, a] = iirnotch(50 / fs * 2, 5 / fs * 2);
y = filter(b, a, y);

[b, a] = iirnotch(110.45 / fs * 2, 5 / fs * 2);
y = filter(b, a, y);

[b, a] = iirnotch(168.6 / fs * 2, 5 / fs * 2);
y = filter(b, a, y);

[b, a] = iirnotch(184.15 / fs * 2, 5 / fs * 2);
y = filter(b, a, y);

[b, a] = iirnotch(242.25 / fs * 2, 5 / fs * 2);
y = filter(b, a, y);

% Nelineární fázové zpoždění, ale u vysokých frekvencí nad frekvence EKG
% nevadí
% phasedelay(b, a, [], fs)

figure; plot(f, abs(fft(y)));
title("Zfiltrované FFT");

figure; plot(x); hold on; plot(y);
legend(["Originál", "Zfiltrovaný"]);