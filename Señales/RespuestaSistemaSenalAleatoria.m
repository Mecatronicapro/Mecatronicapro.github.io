%% Respuesta de un Sistema Lineal a una Señal Aleatoria
% Este script ilustra las relaciones estadísticas entre las señales de
% entrada y salida de un sistema lineal determinista. Las relaciones son
% el valor medio, la autocorrelación y la densidad espectral. Las señales
% de entrada y salida se asumen como procesos estocásticos estacionarios
% de segundo orden. El sistema determinista es de segundo orden y modela
% una estructura o edificación.
%
% El siguiente diagrama muestra el proceso como un sistema lineal.
%
%                  --------------------------
%  Excitación --> | Función de Transferencia | --> Respuesta
%                 |  del Modelo Estructural  |
%                  --------------------------
%
%                         ------
%  Desplazamiento en --> | h(t) | --> Desplazamiento Relativo
%  la Base x(t)           ------      de la Estructura y(t)
%
%                         ------
%  Desplazamiento en --> | H(w) | --> Desplazamiento Relativo
%  la Base X(w)           ------      de la Estructura Y(w)
%
%
% Las funciones requeridas pertenecen a las siguientes toolboxes:
% - Statistics and Machine Learning Toolbox
% - Signal Processing Toolbox
% - DSP System Toolbox
%
% Referencias:
% [1] Carlson 2010, Communication Systems, Chapters 8 and 9
% https://ocw.mit.edu/courses/6-011-introduction-to-communication-control-and-signal-processing-spring-2010/a6bddaee5966f6e73450e6fe79ab0566_MIT6_011S10_notes.pdf
% [2] Oppenheim 2010, Signals Systems and Inference, Chapters 9 and 10
% https://ocw.mit.edu/courses/6-011-introduction-to-communication-control-and-signal-processing-spring-2010/a6bddaee5966f6e73450e6fe79ab0566_MIT6_011S10_notes.pdf
% [3] Peebles 2010, Probability, Random Variables And Random Signal Principles, Chapters 6, 7 and 8
% https://ocw.mit.edu/courses/6-011-introduction-to-communication-control-and-signal-processing-spring-2010/a6bddaee5966f6e73450e6fe79ab0566_MIT6_011S10_notes.pdf
% [4] Heredia 2023, Clases magistrales del Dr. Ernesto Heredia en la PUCP, Octubre 2023.
% http://blog.pucp.edu.pe/blog/maestriaeningenieriacivil/2023/10/23/clase-magistral-analisis-dinamico-de-estructuras-en-el-dominio-del-tiempo-y-de-la-frecuencia/
% http://blog.pucp.edu.pe/blog/maestriaeningenieriacivil/2023/10/23/clase-magistral-analisis-de-vibraciones-aleatorias-de-sistemas-estructurales/
% [5] Chopra 2012, Dynamics of Structures
% https://gacbe.ac.in/images/E%20books/Dynamics%20of%20Structures%20Theory%20and%20Applications%20to%20Earthquake%20Engineering%20A.K.%20Chopra_NNN_bb.pdf

% César D. Salvador
% pcukcsal@upc.edu.pe
% Lunes 31 de octubre de 2023


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
figFormat = '-dpng'; figNameExt = '.png'; figRender = 'zbuffer'; figRes = '-r300';
% figFormat = '-depsc'; figNameExt = '.eps'; figRender = 'zbuffer'; figRes = '-r150';


%% Parámetros del Sistema Estructural
m = 5e3;                % Masa de la losa (kg)
k = 2e6;                % Rigidez de las columnas (N/m)
b = 1e4;                % Amortiguamiento de la estructura (N*s/m)
c = 1;                  % Ganancia del sistema

w0 = sqrt(k/m);         % Frecuencia angular natural de la estructura (rad/s)
f0 = w0/(2*pi);         % Frecuencia natural de la estructura (Hz)
xi = b/(2*w0*m);        % Tasa de amortiguamiento de la estructura


%% Dominios de Tiempo y Frecuencia para el Sistema
Fs = 200;               % Frecuencia de muestreo (Hz)
Duracion = 10;          % Duracion (s)
Nh = Duracion * Fs;     % Número de muestras
th = (0:Nh-1)'/Fs;      % Dominio de tiempo (s)
fH = (0:Nh/2)'*Fs/Nh;   % Dominio de frecuencia (Hz)
wh = 2*pi*fH;           % Dominio de frecuencia angular (rad/s)


%% Función de Transferencia del Sistema
H = wh.^2 ./ (w0^2-wh.^2 + 1j*2*xi*w0*wh);      % H(w) = Y(w) / X(w)
HMagDB = 20*log10(abs(H));                      % Magnitude in dB
HPha = unwrap(angle(H));                        % Phase
HGrD = [0; -diff(HPha)/(2*pi*Fs/Nh)];           % Group Delay


%% Respuesta Impulsiva del Sistema
h = real(ifft([H; conj(H(end-1:-1:2))]));   % h(t) = F^-1{H(w)}

  
%% Señal Aleatoria de Entrada
mux0 = 0;                       % Valor medio deseado
sigmax0 = 6;                    % Desviación estandar deseada
pdx0 = makedist('Normal', ...
    'mu', mux0, ...
    'sigma', sigmax0);          % Función de densidad de probabilidad deseada
rng('default');
Nx = 10000;                     % Número de muestras de x(t)
tx = (0:Nx-1)'/Fs;              % Dominio de tiempo (s) de x(t)
fx = (0:Nx/2)'*Fs/Nx;           % Dominio de frecuencia (Hz) de X(w)
x = random(pdx0, Nx, 1);        % Señal de entrada x(t)
pdx = fitdist(x, 'Normal');     % Objeto función de densidad de probabilidad
xVal = -40:.1:40;               % Dominio de la función de densidad de probabilidad
fxVal = pdf(pdx, xVal);         % Función de densidad de probabilidad
[Rxx, lagsx] = xcorr(x);        % Función de autocorrelación de x(t)
taux = lagsx' * 1/Fs;           % Intervalo de la autocorrelación
Sxx = 1/(2*pi) * real(fft(Rxx));
Sxx = Sxx(1:Nx+1);              % Densidad espectral de x(t)


%% Señal Aleatoria de Salida
y = fftfilt(h, x);              % Señal de salida y(t)
pdy = fitdist(y, 'Normal');     % Objeto función de densidad de probabilidad
yVal = -30:.1:30;               % Dominio de la función de densidad de probabilidad
fyVal = pdf(pdy, yVal);         % Función de densidad de probabilidad
[Ryy, lagsy] = xcorr(y);        % Función de autocorrelación de y(t)
tauy = lagsy' * 1/Fs;           % Intervalo de la autocorrelación
Syy = 1/(2*pi) * real(fft(Ryy));
Syy = Syy(1:Nx+1);              % Densidad espectral de y(t) calculada


%% Comparación de las Densidades Espectrales de la Salida
H2 = fft(h, 2*Nx);              % Función de transferencia del sistema interpolada
H2 = H2(1:Nx+1);
H2Sxx = abs(H2).^2 .* Sxx;      % Densidad espectral de y(t) estimada
fS = (0:Nx)'*Fs/(2*Nx);         % Dominio de frecuencia (Hz) para las densidades espectrales


%% Visualización de Resultados
plot(tx, x, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$x(t)$', 'interpreter', 'latex')
title('Input Signal', 'interpreter', 'latex')
axis tight
box on
grid on
hold on

figure
plot(xVal, fxVal, 'b', 'linewidth', 2)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$f_X(x)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on
hold on

figure
plot(taux, Rxx, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{xx}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on
hold on

figure
plot(fH, HMagDB, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$20\log_{10}|H(f)|$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{\ddot{X}(\omega)}=\frac{-1}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
plot(fH, HPha, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\angle{H(f)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{\ddot{X}(\omega)}=\frac{-1}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
ylim([-1.1*pi 1.1*pi])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
plot(fH, HGrD, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$-\frac{d}{d\omega}\angle{H(\omega)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{Y(\omega)}{\ddot{X}(\omega)}=\frac{-1}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
plot(th, h, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$h(t)$', 'interpreter', 'latex')
title('$h(t)=\mathcal{F}^{-1}\{H(\omega)\}$', 'interpreter', 'latex')
axis tight
box on
grid on
set(gcf, 'renderer', figRender)

figure
plot(tx, y, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$y(t)$', 'interpreter', 'latex')
title('Output Signal', 'interpreter', 'latex')
axis tight
box on
grid on
hold on

figure
plot(yVal, fyVal, 'b', 'linewidth', 2)
xlabel('$y$', 'interpreter', 'latex')
ylabel('$f_Y(y)$', 'interpreter', 'latex')
title('Probability Density Function', 'interpreter', 'latex')
axis tight
box on
grid on
hold on

figure
plot(tauy, Ryy, 'b', 'linewidth', 2)
xlabel('$\tau$ (s)', 'interpreter', 'latex')
ylabel('$R_{yy}(\tau)$', 'interpreter', 'latex')
title('Autocorrelation', 'interpreter', 'latex')
axis tight
box on
grid on
hold on

figure
plot(fS, 10*log10(abs(Sxx)), 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S_{xx}(f)|$', 'interpreter', 'latex')
title('Spectral Density', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
plot(fS, 10*log10(abs(Syy)), 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}|S_{yy}(f)|$', 'interpreter', 'latex')
title('Spectral Density', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
plot(fS, 10*log10(abs(H2Sxx)), 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$10\log_{10}||H(f)|^2 S_{xx}(f)|$', 'interpreter', 'latex')
title('Spectral Density', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
axis tight
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
plot(fS, 10*log10(abs(Syy)), 'b', 'linewidth', 2)
hold on
plot(fS, 10*log10(abs(H2Sxx)), 'm:', 'linewidth', 1)
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
hold on
set(gcf, 'renderer', figRender)


%% Exportación de Figuras
% figName = 'FigH';
% print(figFormat, figRes, [figFolder, figName, figNameExt])

