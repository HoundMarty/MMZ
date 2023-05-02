clc;
close all;
clear *;
%%
load('EKGizo1.mat')
Fs = 500;
tn=(1:length(x))/Fs;                    

% figure;plot(tn,x)

% detrend ??
% dt = detrend(x);
% figure;plot(tn,dt);
%%

% bhi = fir1(34,0.48,'high',chebwin(35,30));

% b=fir1(256,[2/256 30/256]); %Navrh pasmove propusti (PP) 2-30Hz, vlastnosti probereme detailne v tematech 11-12 
% b = fir1(256,0.004,'high',chebwin(257,250));
b = fir1(256,2/500,'high','noscale');
% b = butter(20,2/500,'high');
% figure;
% freqz(b,1,1024);            %frekvencni charakteristika PP
% figure;
% phasedelay(b,1,1024);       %Konstantni fazove zpozdeni 128 vzorku - filtr s linearni fazi
x3=filter(b,1,x);           %Filtrace
% figure;plot(x3);            %Filtrace odstranila drift izolinie v EKG (pomaly trend)

figure;
plot(tn,x3,'r');
hold on;
plot(tn,x,'b');