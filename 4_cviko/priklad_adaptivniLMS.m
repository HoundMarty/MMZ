clear *
N = 10000;

d = randn(N,1);

delay = 3;
systemdelay = 20;

sigma = 0.5;

% x = [d(delay+1:end);zeros(delay,1)] + sigma * randn(N,1);

L = 50; % řád filtru
mu = 0.01; % délka kroku
w = zeros(L,1);
y = zeros(N,1);
Lh = 100;
h = zeros(Lh,1);
x = zeros(size(d));

P = eye(L);
lambda = 0.99; % krok zapomínání
w = [1 zeros(L-1,1)']';
% for n = L:N
% adaptivní delay
for n = max(L,Lh):N-systemdelay
    delay = 3+(n/N)*10;
    h = sinc((0:Lh-1)' - delay);
    x(n) = h'*d(n+systemdelay:-1:n-Lh+1+systemdelay) + sigma * randn(1);
    % RMS
    xn = x(n:-1:n-L+1);
    hn = P*xn;
    kn = hn / (lambda + xn'*hn);
    e = d(n) - w'*xn;
    w = w + kn*e;
    P = P/lambda - kn*xn'*P/lambda;
    y(n) = w'*xn;
    % LMS
%     xn = x(n:-1:n-L+1); % L posledních vzorků
%     y(n) = w'*xn; % výstup
%     e = d(n) - y(n); % chybový signál
%     w = w + mu*xn*e; % jednoduchý update

%     normalized
%     w = w + mu*(xn/norm(xn)^2)*e; 
    
    if mod(n,10)==0 % jen pro průběžné vykreslení výsledku - sledování po 10 krocích
        subplot(1,2,1)
        plot(w)
        axis([1 L -1 1]) % pevné osy, kvůli adaptaci
        subplot(1,2,2)
        plot(1:N,d,'k',1:N,y,'r')
        axis([1 N -4 4])
        drawnow
    end
end

% fluktuace, protože filtr nikdy nedokonverguje
