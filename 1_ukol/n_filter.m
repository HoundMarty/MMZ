clear *
close all;
clc;

%%
load('EKG3channels_sinus.mat')

x = x';

% figure;pspectrum(x(:,1));
% figure;pspectrum(x(:,2));
% pspectrum(x(:,3))

% odhad Fs ze vzorků
dt = mean(diff(x(:,1)));
fs = round(1 / dt); % Hz
% fs = 500;
tn=(1:length(x))/fs;                    

% S = fft(x(:,1))/length(fft(x(:,1)));
% P2 = abs(S)*2;
% plot(P2);
% findpeaks(P2);
% pp = findpeaks(pspectrum(x(:,1)),'MinPeakProminence',300,'Annotate','extents');
% findpeaks(P2,'MinPeakProminence',25);
% pp = findpeaks(P2);

%%
close all;
% sound(x(:,1))
%comb filter
fo = 66.7;%50;
q = 35;
bw = (fo/(fs/2))/q;
[b,a] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
x0 = filtfilt(b,a,x(:,1));

fo = 66.7;%50;
q = 35;
bw = (fo/(fs/2))/q;
[b,a] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
x1 = filtfilt(b,a,x(:,2));

fo = 66.7;%50;
q = 35;
bw = (fo/(fs/2))/q;
[b,a] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
x2 = filtfilt(b,a,x(:,3));
% figure;pspectrum(x0);

% multiple notches
% w1 = 0.1475;%237/(fs/2);
% bw1 = w1/35;
% [b,a] = iirnotch(w1,bw1);
% x1 = filtfilt(b,a,x0);

freqs = [0.1475 0.3802 0.4415 0.5035 0.6745 0.7362 0.9689];
for i = freqs
    w = i;
    bw = w/q;
    [b,a]= iirnotch(w,bw);
    vys0 = filtfilt(b,a,x0); %x0 = signál po průchodu comb filtrem;
    vys1 = filtfilt(b,a,x1); %x1 = signál po průchodu comb filtrem;
    vys2 = filtfilt(b,a,x2); %x2 = signál po průchodu comb filtrem;
    x0 = vys0;
    x1 = vys1;
    x2 = vys2;
end

figure;pspectrum(x0);
figure;pspectrum(x1);
figure;pspectrum(x2);

% kon = bandstop(vys,[2 737 1001 2211 3683 4845],fs);

% figure;pspectrum(kon);
% pspectrum(x(:,2))
% pspectrum(x(:,3))

% figure;
% plot(tn,x3,'r');
