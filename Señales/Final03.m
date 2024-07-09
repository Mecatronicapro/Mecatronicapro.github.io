%% Trabajo Final: Respuesta de un Sistema a una Entrada Aleatoria
% Breve descripción del script.

% Integrantes del grupo
% Carpio Tello Camila Abigail
% Ludeña Macavilca Christian Alexis
% Olivera Bohorquez Enmanuel Marco
% Salcedo Tapara Jose Efrain
% Valdez Olivares Luis Miguel

% Lunes 24 de junio de 2024


%% Inicialización del Entorno Matlab
close all                                   % Cerrar todas las figuras
clearvars                                   % Borrar todas las variables
clc                                         % Limpiar la línea de comandos
workspace                                   % Mostrar el Workspace


%% Parámetros para Visualización
cmap = colormap('parula');
fontSize = 18;
set(0, ...
    'defaultFigureColor', 'w', ...
    'defaultFigureColorMap', cmap, ...
    'defaultAxesFontName',  'times', ...
    'defaultTextFontSize', fontSize, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex', ...
    'defaultAxesLayer', 'top')


%% Parámetros para Exportar Figuras
figFolder = 'Fig\';
figNameExt = '.png';
figRes = 300;


%% Señales de Entrada: Aceleración EO, NS y UD
load("EL231SYSTF2024IData")

% ----------- ACELERACIÓN EO ----------- %
x_EO = a(:, 1);                 % Señal de entrada x(t)
Nx_EO = length(x_EO);           % Número de muestras de x(t)
tx_EO = (0:Nx_EO-1)'/Fs;        % Dominio de tiempo (s) de x(t)
fx_EO = (0:Nx_EO/2)'*Fs/Nx_EO;  % Dominio de frecuencia (Hz) de X(w)

pdx_EO = fitdist(x_EO, 'Normal');       % Objeto función de densidad de probabilidad
xVal_EO = -500:500;                     % Dominio de la función de densidad de probabilidad
fxVal_EO = pdf(pdx_EO, xVal_EO);        % Función de densidad de probabilidad
[Rxx_EO, lagsx_EO] = xcorr(x_EO);       % Función de autocorrelación de x(t)
taux_EO = lagsx_EO' * 1/Fs;             % Intervalo de la autocorrelación
Sxx_EO = 1/(2*pi) * real(fft(Rxx_EO));  % Densidad espectral de x(t)
Sxx_EO = Sxx_EO(1:Nx_EO+1);             % Ajuste de la densidad espectral de x(t)
fS_EO = (0:Nx_EO)'*Fs/(2*Nx_EO);        % Dominio de frecuencia

% ----------- ACELERACIÓN NS ----------- %
x_NS = a(:, 2);                 % Señal de entrada x(t)
Nx_NS = length(x_NS);           % Número de muestras de x(t)
tx_NS = (0:Nx_NS-1)'/Fs;        % Dominio de tiempo (s) de x(t)
fx_NS = (0:Nx_NS/2)'*Fs/Nx_NS;  % Dominio de frecuencia (Hz) de X(w)

pdx_NS = fitdist(x_NS, 'Normal');       % Objeto función de densidad de probabilidad
xVal_NS = -500:500;                     % Dominio de la función de densidad de probabilidad
fxVal_NS = pdf(pdx_NS, xVal_NS);        % Función de densidad de probabilidad
[Rxx_NS, lagsx_NS] = xcorr(x_NS);       % Función de autocorrelación de x(t)
taux_NS = lagsx_NS' * 1/Fs;             % Intervalo de la autocorrelación
Sxx_NS = 1/(2*pi) * real(fft(Rxx_NS));  % Densidad espectral de x(t)
Sxx_NS = Sxx_NS(1:Nx_NS+1);             % Ajuste de la densidad espectral de x(t)
fS_NS = (0:Nx_NS)'*Fs/(2*Nx_NS);        % Dominio de frecuencia

% ----------- ACELERACIÓN UD ----------- %
x_UD = a(:, 3);                 % Señal de entrada x(t)
Nx_UD = length(x_UD);           % Número de muestras de x(t)
tx_UD = (0:Nx_UD-1)'/Fs;        % Dominio de tiempo (s) de x(t)
fx_UD = (0:Nx_UD/2)'*Fs/Nx_UD;  % Dominio de frecuencia (Hz) de X(w)

pdx_UD = fitdist(x_UD, 'Normal');        % Objeto función de densidad de probabilidad
xVal_UD = -500:500;                      % Dominio de la función de densidad de probabilidad
fxVal_UD = pdf(pdx_UD, xVal_UD);         % Función de densidad de probabilidad
[Rxx_UD, lagsx_UD] = xcorr(x_UD);        % Función de autocorrelación de x(t)
taux_UD = lagsx_UD' * 1/Fs;              % Intervalo de la autocorrelación
Sxx_UD = 1/(2*pi) * real(fft(Rxx_UD));   % Densidad espectral de x(t)
Sxx_UD = Sxx_UD(1:Nx_UD+1);              % Ajuste de la densidad espectral de x(t)
fS_UD = (0:Nx_UD)'*Fs/(2*Nx_UD);         % Dominio de frecuencia


%% Sistema Lineal: Circuito Equivalente RLC

R = 18;               % Resistencia en ohmios
L = 9;                % Inductancia en henrios
C = 0.000175;         % Capacitancia en faradios
wn = sqrt(1/(C*L))    % Frecuencia natural en radianes por segundo
xi = R / (2 * L * wn) % Factor de amortiguamiento adimensional

%% Dominios de Tiempo y Frecuencia para el Sistema
Duracion = 10;          % Duracion (s)
Nh = Duracion * Fs;     % Número de muestras
th = (0:Nh-1)'/Fs;      % Dominio de tiempo (s)
fH = (0:Nh/2)'*Fs/Nh;   % Dominio de frecuencia (Hz)
w = 2*pi*fH;            % Dominio de frecuencia angular (rad/s)

%% Función de Transferencia del Sistema
H = -1 ./ (wn^2-w.^2 + 1j*2*xi*wn*w);  % Función de transferencia H(w) = Y(w) / X(w)
HMagDB = 20*log10(abs(H));             % Magnitud en dB
HPha = unwrap(angle(H));               % Fase desenvuelta
HGrD = [0; -diff(HPha)/(2*pi*Fs/Nh)];  % Retardo de grupo

%% Respuesta Impulsiva del Sistema
h = real(ifft([H; conj(H(end-1:-1:2))]));   % h(t) = F^-1{H(w)}

%% Señales de Salida: Desplazamiento EO, NS y UD

%----------------------DESPLAZAMIENTO EO-------------------------
y_EO = fftfilt(h, x_EO);                % Filtrado de la señal de entrada x_EO con el filtro h
pdy_EO = fitdist(y_EO, 'Normal');       % Ajuste de la distribución normal a la señal filtrada y_EO
yVal_EO = -2.5:.001:2.5;                % Dominio de la función de densidad de probabilidad
fyVal_EO = pdf(pdy_EO, yVal_EO);        % Evaluación de la función de densidad de probabilidad
[Ryy_EO, lagsy_EO] = xcorr(y_EO);       % Autocorrelación de la señal filtrada y_EO
tauy_EO = lagsy_EO' * 1/Fs;             % Intervalo de la autocorrelación en segundos
Syy_EO = 1/(2*pi) * real(fft(Ryy_EO));  
Syy_EO = Syy_EO(1:Nx_EO+1);             % Densidad espectral de la señal y_EO

H2_EO = fft(h, 2*Nx_EO);                % Interpolación de la función de transferencia h
H2_EO = H2_EO(1:Nx_EO+1);
H2Sxx_EO = abs(H2_EO).^2 .* Sxx_EO;     % Estimación de la densidad espectral de y_EO
fS_EO = (0:Nx_EO)'*Fs/(2*Nx_EO);        % Dominio de frecuencia para las densidades espectrales

%----------------------DESPLAZAMIENTO NS-------------------------
y_NS = fftfilt(h, x_NS);                % Filtrado de la señal de entrada x_NS con el filtro h
pdy_NS = fitdist(y_NS, 'Normal');       % Ajuste de la distribución normal a la señal filtrada y_NS
yVal_NS = -2.5:.001:2.5;                % Dominio de la función de densidad de probabilidad
fyVal_NS = pdf(pdy_NS, yVal_NS);        % Evaluación de la función de densidad de probabilidad
[Ryy_NS, lagsy_NS] = xcorr(y_NS);       % Autocorrelación de la señal filtrada y_NS
tauy_NS = lagsy_NS' * 1/Fs;             % Intervalo de la autocorrelación en segundos
Syy_NS = 1/(2*pi) * real(fft(Ryy_NS));  
Syy_NS = Syy_NS(1:Nx_NS+1);             % Densidad espectral de la señal y_NS

H2_NS = fft(h, 2*Nx_NS);                % Interpolación de la función de transferencia h
H2_NS = H2_NS(1:Nx_NS+1);
H2Sxx_NS = abs(H2_NS).^2 .* Sxx_NS;     % Estimación de la densidad espectral de y_NS
fS_NS = (0:Nx_NS)'*Fs/(2*Nx_NS);        % Dominio de frecuencia para las densidades espectrales

%----------------------DESPLAZAMIENTO UD-------------------------
y_UD = fftfilt(h, x_UD);                % Filtrado de la señal de entrada x_UD con el filtro h
pdy_UD = fitdist(y_UD, 'Normal');       % Ajuste de la distribución normal a la señal filtrada y_UD
yVal_UD = -2.5:.001:2.5;                % Dominio de la función de densidad de probabilidad
fyVal_UD = pdf(pdy_UD, yVal_UD);        % Evaluación de la función de densidad de probabilidad
[Ryy_UD, lagsy_UD] = xcorr(y_UD);       % Autocorrelación de la señal filtrada y_UD
tauy_UD = lagsy_UD' * 1/Fs;             % Intervalo de la autocorrelación en segundos
Syy_UD = 1/(2*pi) * real(fft(Ryy_UD));  
Syy_UD = Syy_UD(1:Nx_UD+1);             % Densidad espectral de la señal y_UD

H2_UD = fft(h, 2*Nx_UD);                % Interpolación de la función de transferencia h
H2_UD = H2_UD(1:Nx_UD+1);
H2Sxx_UD = abs(H2_UD).^2 .* Sxx_UD;     % Estimación de la densidad espectral de y_UD
fS_UD = (0:Nx_UD)'*Fs/(2*Nx_UD);        % Dominio de frecuencia para las densidades espectrales

%% Visualizacion

%-----------------Aceleracion EO------------------

% Figura 1: Señal de entrada x(t) en el dominio del tiempo
figure(1)
plot(tx_EO, x_EO, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$x(t)$', 'interpreter', 'latex')
title('Input Signal', 'interpreter', 'latex')
axis tight
box on
grid on




% Figura 2: Función de densidad de probabilidad de la señal de entrada x(t)
figure(2)
plot(xVal_EO, fxVal_EO, 'b', 'linewidth', 2)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$f_X(x)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on




% Figura 3: Autocorrelación de la señal de entrada x(t)
figure(3)
plot(taux_EO, Rxx_EO, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{xx}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on






% Figura 4: Densidad espectral de la señal de entrada x(t)
figure(4)
plot(fS_EO, 10*log10(abs(Sxx_EO)), 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S_{xx}(f)|$', 'interpreter', 'latex')
title('Spectral Density', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on



% Figura 5: Magnitud de la función de transferencia H(f)
figure(5)
plot(fH, HMagDB, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$20\log_{10}|H(f)|$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
box on
grid on




% Figura 6: Fase de la función de transferencia H(f)
figure(6)
plot(fH, HPha, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\angle{H(f)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
ylim([-1.1*pi 1.1*pi])
xlim([.25 32])
box on
grid on



% Figura 7: Derivada de la fase de la función de transferencia H(f)
figure(7)
plot(fH, HGrD, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$-\frac{d}{d\omega}\angle{H(\omega)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on





% Figura 8: Respuesta al impulso del sistema h(t)
figure(8)
plot(th, h, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$h(t)$', 'interpreter', 'latex')
title('$h(t)=\mathcal{F}^{-1}\{H(\omega)\}$', 'interpreter', 'latex')
axis tight
box on
grid on



% Figura 9: Señal de salida y(t) en el dominio del tiempo
figure(9)
plot(tx_EO, y_EO, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$y(t)$', 'interpreter', 'latex')
title('Output Signal', 'interpreter', 'latex')
axis tight
box on
grid on



% Figura 10: Función de densidad de probabilidad de la señal de salida y(t)
figure(10)
plot(yVal_EO, fyVal_EO, 'b', 'linewidth', 2)
xlabel('$y$', 'interpreter', 'latex')
ylabel('$f_Y(y)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on




% Figura 11: Autocorrelación de la señal de salida y(t)
figure(11)
plot(tauy_EO, Ryy_EO, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{yy}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on





% Figura 12: Densidad espectral de la señal de salida y(t) comparada con la estimada
figure(12)
plot(fS_EO, 10*log10(abs(Syy_EO)), 'b', 'linewidth', 2)
hold on
plot(fS_EO, 10*log10(abs(H2Sxx_EO)), 'm:', 'linewidth', 1)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S|$', 'interpreter', 'latex')
title('Spectral Density $S$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
legend('$S_{yy}(f)$', '$|H(f)|^2 S_{xx}(f)$', ...
    'location', 'northoutside', ...
    'orientation', 'horizontal', ...
    'interpreter', 'latex')
box on
grid on



%-----------------aceleracion NS------------------

% Figura 13: Señal de entrada x(t) en el dominio del tiempo
figure(13)
plot(tx_NS, x_NS, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$x(t)$', 'interpreter', 'latex')
title('Input Signal', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 14: Función de densidad de probabilidad de la señal de entrada x(t)
figure(14)
plot(xVal_NS, fxVal_NS, 'b', 'linewidth', 2)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$f_X(x)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 15: Autocorrelación de la señal de entrada x(t)
figure(15)
plot(taux_NS, Rxx_NS, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{xx}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 16: Densidad espectral de la señal de entrada x(t)
figure(16)
plot(fS_NS, 10*log10(abs(Sxx_NS)), 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S_{xx}(f)|$', 'interpreter', 'latex')
title('Spectral Density', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on


% Figura 17: Magnitud de la función de transferencia H(f)
figure(17)
plot(fH, HMagDB, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$20\log_{10}|H(f)|$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
box on
grid on



% Figura 18: Fase de la función de transferencia H(f)
figure(18)
plot(fH, HPha, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\angle{H(f)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
ylim([-1.1*pi 1.1*pi])
box on
grid on



% Figura 19: Derivada de la fase de la función de transferencia H(f)
figure(19)
plot(fH, HGrD, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$-\frac{d}{d\omega}\angle{H(\omega)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on


% Figura 20: Respuesta al impulso del sistema h(t)
figure(20)
plot(th, h, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$h(t)$', 'interpreter', 'latex')
title('$h(t)=\mathcal{F}^{-1}\{H(\omega)\}$', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 21: Señal de salida y(t) en el dominio del tiempo
figure(21)
plot(tx_NS, y_NS, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$y(t)$', 'interpreter', 'latex')
title('Output Signal', 'interpreter', 'latex')
axis tight
box on
grid on

% Figura 22: Función de densidad de probabilidad de la señal de salida y(t)
figure(22)
plot(yVal_NS, fyVal_NS, 'b', 'linewidth', 2)
xlabel('$y$', 'interpreter', 'latex')
ylabel('$f_Y(y)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 23: Autocorrelación de la señal de salida y(t)
figure(23)
plot(tauy_NS, Ryy_NS, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{yy}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 24: Densidad espectral de la señal de salida y(t) comparada con la estimada
figure(24)
plot(fS_NS, 10*log10(abs(Syy_NS)), 'b', 'linewidth', 2)
hold on
plot(fS_NS, 10*log10(abs(H2Sxx_NS)), 'm:', 'linewidth', 1)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S|$', 'interpreter', 'latex')
title('Spectral Density $S$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
legend('$S_{yy}(f)$', '$|H(f)|^2 S_{xx}(f)$', ...
    'location', 'northoutside', ...
    'orientation', 'horizontal', ...
    'interpreter', 'latex')
box on
grid on



%-----------------aceleracion UD------------------

% Figura 25: Señal de entrada x(t) en el dominio del tiempo
figure(25)
plot(tx_UD, x_UD, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$x(t)$', 'interpreter', 'latex')
title('Input Signal', 'interpreter', 'latex')
axis tight
box on
grid on



% Figura 26: Función de densidad de probabilidad de la señal de entrada x(t)
figure(26)
plot(xVal_UD, fxVal_UD, 'b', 'linewidth', 2)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$f_X(x)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 27: Autocorrelación de la señal de entrada x(t)
figure(27)
plot(taux_UD, Rxx_UD, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{xx}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 28: Densidad espectral de la señal de entrada x(t)
figure(28)
plot(fS_UD, 10*log10(abs(Sxx_UD)), 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S_{xx}(f)|$', 'interpreter', 'latex')
title('Spectral Density', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on


% Figura 29: Magnitud de la función de transferencia H(f)
figure(29)
plot(fH, HMagDB, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$20\log_{10}|H(f)|$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
box on
grid on



% Figura 30: Fase de la función de transferencia H(f)
figure(30)
plot(fH, HPha, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\angle{H(f)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
ylim([-1.1*pi 1.1*pi])
box on
grid on



% Figura 31: Derivada de la fase de la función de transferencia H(f)
figure(31)
plot(fH, HGrD, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$-\frac{d}{d\omega}\angle{H(\omega)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{X(\omega)}=\frac{-1}{\omega_n^2-\omega^2+j2\xi\omega_n\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on



% Figura 32: Respuesta al impulso del sistema h(t)
figure(32)
plot(th, h, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$h(t)$', 'interpreter', 'latex')
title('$h(t)=\mathcal{F}^{-1}\{H(\omega)\}$', 'interpreter', 'latex')
axis tight
box on
grid on


% Figura 33: Señal de salida y(t) en el dominio del tiempo
figure(33)
plot(tx_UD, y_UD, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$y(t)$', 'interpreter', 'latex')
title('Output Signal', 'interpreter', 'latex')
axis tight
box on
grid on



% Figura 34: Función de densidad de probabilidad de la señal de salida y(t)
figure(34)
plot(yVal_UD, fyVal_UD, 'b', 'linewidth', 2)
xlabel('$y$', 'interpreter', 'latex')
ylabel('$f_Y(y)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on

% Figura 35: Autocorrelación de la señal de salida y(t)
figure(35)
plot(tauy_UD, Ryy_UD, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{yy}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on



% Figura 36: Densidad espectral de la señal de salida y(t) comparada con la estimada
figure(36)
plot(fS_UD, 10*log10(abs(Syy_UD)), 'b', 'linewidth', 2)
hold on
plot(fS_UD, 10*log10(abs(H2Sxx_UD)), 'm:', 'linewidth', 1)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S|$', 'interpreter', 'latex')
title('Spectral Density $S$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
legend('$S_{yy}(f)$', '$|H(f)|^2 S_{xx}(f)$', ...
    'location', 'northoutside', ...
    'orientation', 'horizontal', ...
    'interpreter', 'latex')
box on
grid on



%% Exportación de Figuras
%figName = 'FigSignal';
%exportgraphics(gcf, [figFolder, figName, figNameExt], 'Resolution', figRes)

