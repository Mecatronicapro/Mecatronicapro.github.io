%% Análisis Estructural Dinámico de una Edificación
% Modelo dinámico de una estructura en los dominios de tiempo y frecuencia.
% El siguiente diagrama muestra el proceso como un sistema lineal.
%
%                  --------------------------
%  Excitación --> | Función de Transferencia | --> Respuesta
%                 |  del Modelo Estructural  |
%                  --------------------------
%
%                              ------
%  Aceleración en la Base --> | H(w) | --> Desplazamiento Relativo
%                              ------
%
%                                 ------
%  Desplazamiento en la Base --> | G(w) | --> Desplazamiento Relativo
%                                 ------
%
% Referencias:
% [1] Notas de la clase magistral del Dr. Ernesto Heredia en la PUCP, Octubre 2023.
% [2] Anil Chopra, Dynamics of Structures: Theory and Application to
% Earthquake Engineering, Fourth Edition, Prentice Hall, 2012.

% César D. Salvador
% pcukcsal@upc.edu.pe
% Lunes 9 de octubre de 2023


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


%% Parámetros Estructurales
M = 5e3;                % Masa de la losa (kg)
k = 2e6;                % Rigidez de las columnas (N/m)
C = 1e4;                % Amortiguamiento de la estructura (N*s/m)
w0 = sqrt(k/M);         % Frecuencia angular natural de la estructura (rad/s)
f0 = w0/(2*pi);         % Frecuencia natural de la estructura (Hz)
xi = C/(2*w0*M);        % Tasa de amortiguamiento de la estructura


%% Dominios de Tiempo y Frecuencia
Fs = 200;               % Frecuencia de muestreo (Hz)
Duracion = 10;          % Duracion (s)
N = Duracion * Fs;      % Número de muestras
t = (0:N-1)'/Fs;        % Dominio de tiempo (s)
f = (0:N/2)'*Fs/N;      % Dominio de frecuencia (Hz)
w = 2*pi*f;             % Dominio de frecuencia angular (rad/s)


%% Funciones de Transferencia
H = -1 ./ (w0^2-w.^2 + 1j*2*xi*w0*w);          % H = Y/d^2U
G = w.^2 ./ (w0^2-w.^2 + 1j*2*xi*w0*w);        % G = X/U


%% Respuestas Impulsivas Normalizadas
h = real(ifft([H; conj(H(end-1:-1:2))]));
g = real(ifft([G; conj(G(end-1:-1:2))]));


%% Visualización de Resultados
V = 20*log10(abs(H));
plot(f, V, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$20\log_{10}|H(f)|$', 'interpreter', 'latex')
title('$H(\omega)=\frac{X(\omega)}{\ddot{U}(\omega)}=\frac{-1}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
V = angle(H);
plot(f, V, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\angle{H(f)}$', 'interpreter', 'latex')
title('$H(\omega)=\frac{X(\omega)}{\ddot{U}(\omega)}=\frac{-1}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
ylim([-1.1*pi 1.1*pi])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
V = h;
plot(t, V, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$h(t)$', 'interpreter', 'latex')
title('$h(t)=\mathcal{F}^{-1}\{H(\omega)\}$', 'interpreter', 'latex')
axis tight
box on
grid on
set(gcf, 'renderer', figRender)

figure
V = 20*log10(abs(G));
plot(f, V, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$20\log_{10}|G(f)|$', 'interpreter', 'latex')
title('$G(\omega)=\frac{X(\omega)}{U(\omega)}=\frac{\omega^2}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
V = angle(G);
plot(f, V, 'b', 'linewidth', 2)
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\angle{G(f)}$', 'interpreter', 'latex')
title('$G(\omega)=\frac{X(\omega)}{U(\omega)}=\frac{\omega^2}{\omega_0^2-\omega^2+j2\xi\omega_0\omega}, \omega=2\pi f$', 'interpreter', 'latex')
set(gca, 'xscale', 'log')
set(gca, 'xtick', 2.^(-3:6), 'xticklabel', 2.^(-3:6))
xlim([.25 32])
ylim([-1.1*pi 1.1*pi])
box on
grid on
hold on
set(gcf, 'renderer', figRender)

figure
V = g;
plot(t, g, 'b', 'linewidth', 2)
xlabel('$t$ (s)', 'interpreter', 'latex')
ylabel('$g(t)$', 'interpreter', 'latex')
title('$g(t)=\mathcal{F}^{-1}\{G(\omega)\}$', 'interpreter', 'latex')
axis tight
box on
grid on
set(gcf, 'renderer', figRender)


%% Exportación de Figuras
% figName = 'FigH';
% print(figFormat, figRes, [figFolder, figName, figNameExt])

