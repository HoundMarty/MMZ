clear *
close all
clc

%% zadání

% Vizualizujte charakteristickou funkci delay-and-sum beamformeru za předpokladu, že

% zvuk se šíří ideálním prostředím bez odrazů rychlostí 330 m/s, 
% zvuk se šíří pouze v rovině xy (tedy v podstatě ve 2D) a
% tři mikrofony jsou rozmístěny lineárně podél osy x a jsou vzájemně vzdálené 10 cm.

% Charakteristickou funkci zobrazte jako polární graf pro zvolenou frekvenci a zvolená zpozdění beamformeru.

%% řešení
load('multichannel.mat')
%počet senzorů
m = 9;
%vzdálenost mezi senzory | pozor na jednotky
d = 0.1; % v cm
%rychlost šíření signálu m/s
c = 320;
%pozice senzorů
Pos = [(-(m-1)/2:(m-1)/2)*d; zeros(1,m); zeros(1,m)]; %[x,y]

%f 
f = 3400;
% Delay??
D = zeros(m,1);

% theta
theta = 0:0.01:2*pi;

PSI = zeros(length(theta),1);

%M_k
M = [];
for k=1:length(theta)
    u = -[cos(theta(k)) -sin(theta(k)) 0]';
    M = (u'*Pos/c*fs)';
    PSI(k) = mean(exp(1i*f/fs*2*pi*(D-M)));
end

polarplot(theta,abs(PSI))

% pp = find(ismember(abs(PSI),1));
