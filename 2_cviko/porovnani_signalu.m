% eegplot a fastif -> v nástrojích pro práci s EGG soubory -> eleearning
close all;
clear *;
clc;
%%
load('C:\Users\98653\Documents\Skola\2_rocnik\LS\NMZ\2_cviko\SSVEP1_14Hz.mat')

% C_ij = mean(EEGdata-mean(EEGdata)*EEGdata-mean(EEGdata));
%%

% for i=1:size(EEGdata,1)
%     for j=1:size(EEGdata,1)
%         C(i,j) = mean((EEGdata(i,:)-mean(EEGdata(i,:))) .* (EEGdata(j,:)-mean(EEGdata(j,:))));
%     end
% end

% imagesc(C)

x = EEGdata;

% to samé jako for cykly -> vektorově
x = x - mean(x,2);
C = x*x'/length(x);

%% normované
x = x./sqrt(mean(x.^2,2)); % rozptyl normovaný
C = x*x'/length(x);

imagesc(C)

%% korelační matice
shift = 5;
L = 200;
for i = L/2+1:shift:length(x)
%     window = x(:,max(i-L/2,1):min(i+L/2,length(x)));
    window = x(:,max(i-L/2,1):min(i+L/2,length(x))).*(hamming(L+1)');
    window = window - mean(window,2);
    window = window./mean(window.^2,2);
    C = window*window'/L;
    imagesc(C)
    drawnow
end


%%
load('synch_prum.mat')

p1 = xcorr(x1,x2,50);
pp1 = findpeaks(p1);
plot(-50:50,p1);

p2 = xcorr(x1,x3,50);
% pp2 = findpeaks(p2,'MinPeakHeight',500);
figure;plot(-50:50,p2);
p3 = xcorr(x1,x4,50);
% pp3 = findpeaks(p3,'MinPeakHeight',500);
figure;plot(-50:50,p3);
p4 = xcorr(x1,x5,50);
% pp4 = findpeaks(p4,'MinPeakHeight',500);
figure;plot(-50:50,p4);
p5 = xcorr(x1,x5,50);
figure;plot(-50:50,p5);

zarovnani_x1 = x1;
zarovnani_x2 = x2(:,19:end);
zarovnani_x3 = [zeros(3,1)',x3];
zarovnani_x4 = x4(:,4:end);
zarovnani_x5 = [zeros(12,1)', x5];

%cut
zarovnani_x1 = zarovnani_x1(:,12:390);
zarovnani_x2 = zarovnani_x2(:,12:390);
zarovnani_x3 = zarovnani_x3(:,12:390);
zarovnani_x4 = zarovnani_x4(:,12:390);
zarovnani_x5 = zarovnani_x5(:,12:390);

figure
plot(zarovnani_x1)
hold on
plot(zarovnani_x2)
plot(zarovnani_x3)
plot(zarovnani_x4)
plot(zarovnani_x5)

synchroni_prumer = (zarovnani_x1 + zarovnani_x2 + zarovnani_x3 + zarovnani_x4 + zarovnani_x5) / 5;

plot(synchroni_prumer,LineWidth=3,Color='magenta');
hold off
% figure
% plot(x1)
% hold on
% plot(x2(19:end))
% plot(x3(-2:end))
% plot(x4(-11:end))

figure;plot(abs(fft(synchroni_prumer)))

