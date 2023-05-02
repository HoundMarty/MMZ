clear *
close all
clc


%% Multi-channel LMS filtr
% Nahrajte si do Matlabu data ze souboru multichannel.mat. 
% Signály s a v obsahují stereofonní záznamy signálů z reproduktorů, které byly 
% umístěny v různých pozicích.

% Nejprve spočítejte LMS filtr, který se snaží 
% ze signálu x(:,1) vyfiltrovat signál s(:,1) a potom z 
% celého x získat s(:,1). Pro výpočet filtrů použijte jen počáteční 
% 3 vteřiny signálů, aby bylo zřejmé, na kterých usecích signálu spočtený filtr funguje.

%%

load("multichannel.mat")

% sounsc(x,fs)

% pouze první 3 vteřiny, filtr délky 100, pouze první kanál
% w = miso_firwiener(100,x(1:3*fs,1),s(1:3*fs,1));
% výstup
% y = filter(w,1,x(:,1));

% soundsc(y,fs);

w = miso_firwiener(1000,x(1:3*fs,:),s(1:3*fs,1));

y = filter(w(1:1001),1,x(:,1))+ filter(w(1002:2002),1,x(:,2));

% soundsc(y,fs);
plot(w)

%% žena se pohybuje... spočítáme pro celé
w = miso_firwiener(1000,x,s(:,1)); %

y = filter(w(1:1001),1,x(:,1))+ filter(w(1002:2002),1,x(:,2));

% soundsc(y,fs);
plot(w)
