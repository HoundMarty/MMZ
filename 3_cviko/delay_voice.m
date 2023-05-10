clear all
close all
clc

%%
% Použijte data v souboru delayedvoice.mat.

% 1) Najděte kauzální nebo nekauzální filtr, který filtrováním signálu x dává signál co nejvíce podobný signálu d. Filtr spočítejte řešením soustavy lineárních rovnic a potom pomocí Levinson-Durbinova algoritmu (funkce miso_firwiener.m). Filtr počítejte z celých dat.

% 2) Simulujte jiný typ pozorování x, který vznikne degradací signálu d. Příklad: přičtení bílého šumu, přičtení úzkopásmového šumu/harmonického signálu.

% 3) Vyhodnoťte poměr signál/šum (SNR) před a po zpracování signálu x filtrem.
%%
load('C:\Users\98653\Documents\Skola\2_rocnik\LS\NMZ\3_cviko\delayedvoice.mat')


%%
N = length(x);
L = 100;

R = zeros(L);
p = zeros(L,1);
%%
% 2) bílej šum
x_white = d + randn(1,N)/100; % setina šumu
%%
% filter 
[B,A] = butter(4,[1000/fs*2 1100/fs*2], 'bandpass');
noise = filter(B,A,randn(1,N)/100);
x = d + noise;

%% pro SNR
energie_d = var(d);
energie_noise = var(noise);
SNR_pre = energie_d / energie_noise;
SNR_pre_db = 10*log10(SNR_pre);
% když je SNR < 0 , tak je šum silnější 
% když je SNR > 0 , tak je užitečný signál silnější

%%
for n = L:N
    x_n = x(n:-1:n-L+1)';
    R = R + x_n * x_n';

    p = p + x_n * d(n);
end

R = R/N;
p = p/N;

w = R\p;
% fvtool(w,1) -> lepší jak freqz

%% levinson-darbin
WW = miso_firwiener(L-1,x',d');

%% 
% po zpracování máme 
% y[n] = WW <konvoluce> x[n] = WW * (s[n] + v[n]) = WW * s[n] + WW * v[n]
energie_wd = var(filter(WW,1,d));
energie_wnoise = var(filter(WW,1,noise));
SNR_post = energie_wd / energie_wnoise;
SNR_post_db = 10*log10(SNR_post);
% zlepšení
zlepseni = SNR_post_db - SNR_pre_db;
% SNR neodráží zkreslení signálu -> zlepšení SNR neznamená "hezký"
% poslechový signál
y = filter(WW,1,x);

%%
sound(x,fs);
%vs
sound(y,fs);
%%
plot(w)
hold on
plot(WW)
% skoro stejný výsledky
% miso_firwiener je mnohokrát rychlejší - Levinson-Durbinův algo 

%% SNR
% poměr energie užitečné složky signálu a rušivé složky.
% x[n] = s[n] + v [n];

% SNR = sigma^2_s / sigma^2_v
% což lze taky jako
% 1/N * sum( s[n] - mean(s))^2   /  1/N * sum( v[n] - mean(v) )^2
%%
plot(d)
hold on
plot(filter(w,1,x))

%%
%nekauzální 
% nedělá nic, jen spozdí o L/2 vzorků
N = length(x);
L = 100;

x_n = zeros(L);
p = zeros(L,1);

g_delay = L/2;

for n = L:N
    x_n = x(n:-1:n-L+1)';
    R = R + x_n * x_n';

    p = p + x_n * d(n-g_delay);
end

R = R/N;
p = p/N;

w = R\p;

%%
plot(d)
hold on
plot(filter(w,1,x))