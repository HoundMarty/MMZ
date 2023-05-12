clear *
close all
clc

%% zadání

% Odstranění plovoucí stejnosměrné složky v signálu pomocí filtru.

%% Jak na to?

% Navrhněte vhodný high-pass filtr s lineární fází s přechodovou frekvencí 2 Hz a aplikujte ho na signál.

%% řešení  **Vytvořil** Martin Motejlek
% Spodní frekvence v EKG uvažována 2 Hz podle zadání.

load("EKGizo1.mat");
load("filters.mat");

%%
figure;plot(x);
title("Původní");

%%
y_window = filter(b_window_big_order, 1, x);

% phasedelay(b_window_big_order, 1)

% S kompenzací zpoždění.
figure;plot(y_window(559 + 1 : end));
title("Window (vysoký řád)");

%%
% Cut-off 3 Hz místo 2 Hz!
y_window = filter(b_window, 1, x);

% phasedelay(b_window, 1)

% S kompenzací zpoždění.
figure;plot(y_window(196 + 1 : end));
title("Window (nižší řád)");

%%
[b_butterworth, a_butterworth] = sos2tf(sos_butterworth, g_butterworth);
y_butterworth = filter(b_butterworth, a_butterworth, x);

% phasedelay(b_butterworth, a_butterworth)

% Bez kompenzace zpoždění.
figure;plot(y_butterworth);
title("Butterworth (nelineární fáze – bez kompenzace zpoždění)");

%% POZOR

% Window má velké zpozdění (vysoký řád), ale má lineární fázi.
% Butterworth má velmi nízké zpoždění, ale změní tvar křivky,
% protože nemá lineární fázi.

% Filtr s nelineární fází není vhodný pro EKG, protože změní tvar
% v časové oblasti.
