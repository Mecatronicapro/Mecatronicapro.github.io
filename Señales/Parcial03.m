%% Trabajo Parcial: Diseño de Filtros en el Dominio de Laplace
% Breve descripción del script.
% Utiliza un filtro rechaza banda elíptico para reducir el ruido en una señal de audio.

% Integrantes del grupo
% Carpio Tello, Camila Abigail
% Ludeña Macavilca, Christian
% Olivera Bohorquez, Enmanuel Marco
% Salcedo Tapara, Jose Efrain
% Valdez Olivares, Luis Miguel

% Fecha

%% Inicialización del Entorno Matlab
close all       % Cierra todas las figuras abiertas
clear sound     % Limpia el sonido
clc             % Limpieza de la ventana de comandos
clearvars       % Eliminación de variables de la memoria
close all       % Cierre de todas las ventanas de figuras


%% Parámetros para Visualización
cmap = colormap('parula');
fontSize = 18;
set(0, 'defaultFigureColor', 'w', ...
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


%% Lectura de la Señal de Audio de Entrada
[x, Fs] = audioread('EL231SYSTP2024IEL63Audio.wav');  % x: Señal de audio; Fs: Tasa de muestreo (muestras/segundo)
Nx = length(x);             % Número de muestras de la señal de audio
tx = (0:Nx-1)' * 1/Fs;      % Dominio de tiempo (segundos) de la señal de audio (Vector de tiempo )

%% Transformación de Fourier de la Señal de Audio de Entrada
X = fft(x);                 % Calcula la Transformada de Fourier de la señal de audio
X = X(1:Nx/2+1);            % Selecciona solo las frecuencias positivas con la notación de corte de matriz
fx = (0:Nx/2)' * Fs/Nx;     % Calcula el dominio de frecuencia de la señal de audio en Hz (Vector)

%% Parámetros del Filtro Rechaza Banda
Ap = 0.1;        % Atenuación de la banda de paso (dB)
As = 70;         % Atenuación de la banda de rechazo (dB)
fp1 = 5095;      % Frecuencia de rechazo inferior (Hz)
fs1 = 5930;      % Frecuencia de paso inferior (Hz)
fs2 = 7055;      % Frecuencia de paso superior (Hz)
fp2 = 8290;      % Frecuencia de rechazo superior (Hz)

% Conversión de frecuencias a radianes por segundo
wp1 = 2*pi*fp1; % Frecuencia de rechazo inferior (rad/s)
ws1 = 2*pi*fs1; % Frecuencia de paso inferior (rad/s)
ws2 = 2*pi*fs2; % Frecuencia de paso superior (rad/s)
wp2 = 2*pi*fp2; % Frecuencia de rechazo superior (rad/s)

% Cálculo del parámetro discriminante
epsil = sqrt(10^(Ap/10)-1);
delta = sqrt(10^(-As/10));
d = epsil/sqrt(delta^(-2)-1);   % Parámetro discriminante

%% Dominio de Laplace (plano complejo s)
Nh = 1024;                              % Número de muestras del filtro
fMax = Fs/2;                            % Frecuencia máxima (Fs/2)
f = (-Nh/2+1:Nh/2)'*Fs/Nh;              % Dominio de frecuencia lineal (Hz)
wi = 2*pi*f;                             % Dominio de frecuencia angular (rad/s)
sigma = wi';                             % Dominio sigma
s = ones(Nh,1)*sigma + 1j*wi*ones(1,Nh); % Plano complejo s

%% Transformación en Frecuencia: s --> Ts
Ts =((wp2-wp1)*s)./(s.^2+wp1*wp2);  % Transformación en Frecuencia: s --> Ts
wsStar = min(((wp2-wp1)*ws1)/((wp1*wp2)-(ws1^2)),((wp2-wp1)*ws2)/((ws2^2)-wp1*wp2)); % Cálculo de wsStar
k = 1/wsStar;   % Parámetro de selectividad

k1 = (1-(k^2))^0.25;  % Cálculo de Estrategia
q0 = 0.5*((1-k1)/(1+k1));  % Cálculo de qzero
q = q0 + 2*(q0^5) + 15*(q0^9) + 150*(q0^13);  % Cálculo de q
N = ceil(log10(16/(d^2)) / log10(1/q));  % Orden del filtro

%% Cálculos para el filtro rechaza banda elíptico

beta=(1/(2*N))*log((sqrt(1+(epsil^2))+1)/(sqrt(1+(epsil^2))-1));

if rem(N, 2) == 0   % Comprobación si el orden del filtro es par o impar
    n=1:(0.5*N);    % Definición de índices para el caso par (OJO i=n, de acuerdo a la formula)
    l=n-0.5;
else
    n=1:(0.5*(N-1));% Definición de índices para el caso impar
    l=n; 
end

niterminos=length(n);

% Inicialización de variables para almacenar sumatorias
wsum_num = zeros(1,niterminos);
wsum_dem = zeros(1,niterminos);
wi=zeros(1,niterminos);           % Término W_i
Vi = zeros(1,niterminos);
ai=zeros(1,niterminos);
b1=zeros(1,niterminos);
ci=zeros(1,niterminos);
asum_num = zeros(1,1);
asum_den = zeros(1,1);

% Sumatorias para calcular el coeficiente 'a'
for m=0:10
    asum_num = asum_num+((-1)^m)*(q^(m*(m+1)))*sinh((2*m+1)*beta);
    asum_den = asum_den+((-1)^m)*(q^(m*m))*cosh((2*m*beta));
end 

a=(2*(q^0.25)*asum_num)/(1-1+asum_den);  % Coeficiente 'a' para el filtro

U=sqrt((1+k*(a^2))*(1+((a^2)/k)));  % Factor U para el filtro

% Cálculo de coeficientes para cada término del filtro
for i1=n
    for m = 0:10
        % Sumatorias para calcular el coeficiente 'W_i'
        wsum_num(i1) = wsum_num(i1)+((-1)^m)*(q^(m*(m+1)))*sin((2*m+1)*pi*l(i1)/N);
        wsum_dem(i1) = wsum_dem(i1)+(((-1)^m)*((q^(m*m))*cos(2*m*pi*l(i1)/N)));
    end
    wi(i1)=(2*(q^0.25)*wsum_num(i1))/(1-2+2*wsum_dem(i1));  % Coeficiente 'W_i' para el filtro
    Vi(i1)=sqrt((1+k*(wi(i1)^2))*(1+((wi(i1)^2)/k)));  % Factor V para el filtro
    ai(i1)=1/(wi(i1)^2);  % Coeficiente 'a' en el numerador para el término i1
    b1(i1)=(2*a*Vi(i1))/(1+((a^2)*(wi(i1)^2)));  % Coeficiente 'b' en el numerador para el término i1
    ci(i1)=((a*Vi(i1))^2 +(wi(i1)*U)^2)/((1+((a^2)*(wi(i1)^2)))^2);  % Coeficiente 'c' en el numerador para el término j
end


%% Función de Transferencia (Hs) en el dominio de Laplace (s)
if rem(N, 2) == 0   % Comprobación si el orden del filtro es par o impar
    Ho=1/sqrt(epsil^(2)+1);            % Inicialización de Hzero para el caso par
    for n = 1:0.5*(N)
        Ho = Ho*(ci(n)/ai(n));  % Multiplicación de coeficientes para Ho
    end

    Hs=Ho;  % Inicialización de la función de transferencia para el caso par
    for n = 1:0.5*(N)
        % Cálculo de la función de transferencia para cada término
        Hs=Hs.*(((Ts.^2)+ai(n))./((Ts.^2)+b1(n).*Ts+ci(n)));
    end

else
    Ho = a;             % Inicialización de Hzero para el caso impar
    for n = 1:0.5*(N-1)
        Ho = Ho*(ci(n)/ai(n));  % Multiplicación de coeficientes para Ho
    end

    Hs=Ho./(Ts+a);  % Inicialización de la función de transferencia para el caso impar
    for n = 1:0.5*(N-1)
        % Cálculo de la función de transferencia para cada término
        Hs=Hs.*(((Ts.^2)+ai(n))./((Ts.^2)+b1(n).*Ts+ci(n)));
    end
end

%% Función de Transferencia en el Dominio de Fourier
Hw = Hs(Nh/2:end, Nh/2);    % Restricción al eje imaginario positivo

%% Respuesta Impulsiva en el Dominio del Tiempo
h = real(ifft( [Hw; conj(Hw(end-1:-1:2))] ));   % Transformación inversa de Fourier para obtener la respuesta impulsiva
h = circshift(h, Nh/2);                         % Retraso circular para hacer la respuesta impulsiva causal
h = h / max(abs(h));                            % Normalización de la respuesta impulsiva
th = (0:Nh-1)' * 1/Fs;                          % Dominio de tiempo de la respuesta impulsiva


%% Aplicación de la Respuesta Impulsiva al Audio Original
y = fftfilt(h, x);      % Convolución rápida entre la respuesta impulsiva y el audio original
y = y / max(abs(y));    % Normalización de la señal resultante

%% Transformación de Fourier de la Señal de Audio Filtrada
Y = fft(y);       % Transformación de Fourier de la señal filtrada
Y = Y(1:Nx/2+1);  % Selecciona las frecuencias positivas

%% Escritura de la Señal de Audio de Salida
audiowrite('Filtrado.wav', y, Fs); % Guarda la señal de audio filtrada en un archivo WAV

%% Visualizacion

% Gráfico de la señal de entrada de audio en el dominio del tiempo:
plot(tx, x, 'b', 'linewidth', 2)                  % Plotea la señal de audio en función del tiempo (b=Blue,grosor de linea=2)
xlabel('$t$ (segundos)', 'interpreter','latex')   % Etiqueta del eje x
ylabel('$x(t)$', 'interpreter','latex')           % Etiqueta del eje y (amplitud)
title('Audio Original');                          % Título del gráfico
axis tight                                        % Ajusta los límites de los ejes para que se ajusten a los datos
ylim([-1 1])                                      % Limita el rango de amplitud de la señal
grid on                                           % Activa la cuadrícula en el gráfico
box on  

% Gráfico de la Magnitud del Audio de entrada en el Dominio de la Frecuencia:
figure                                                % Crea una nueva figura para el gráfico
plot(fx*1e-3, 20*log10(abs(X)), 'b', 'linewidth', 2)  % Plotea el espectro de amplitud en función de la frecuencia
xlabel('$f$ (kHz)', 'interpreter','latex')                  % convierte las unidades de Hz a kHz (Para mas legibilidad)
ylabel('$20\log_{10}|X(f)|$ (dB)', 'interpreter','latex')   % calcula el espectro de amplitud en decibelios
title('Audio Original');                              % Título del gráfico
% Configuración del eje x en escala logarítmica
set(gca, 'xscale', 'log')                             % Configura la escala logarítmica para el eje x
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4)) % Establece los marcadores y etiquetas del eje x
xlim([.1 20])                                         % Limita el rango del eje x entre 0.1 kHz y 20 kHz
grid on                                               % Activa la cuadrícula en el gráfico
box on                                                % Dibuja un cuadro alrededor del gráfico


% Definición de los límites de color y los ticks del color
figure
clim = [-36 36];
ctick = clim(1):12:clim(2);
% Creación de la superficie 3D de la magnitud de la función de transferencia
surf(sigma*1e-3, f*1e-3, 20*log10(abs(Hs)), 'EdgeColor', 'None')
xlabel('$\sigma$ (1/ms)')  
ylabel('$f$ (kHz)')
title('Funcion de Transferencia');
axis tight           
set(gca, 'xscale', 'linear', 'yscale', 'linear')  % Ajustes del rango de color
set(gca, 'clim', clim)    % Configuración de la barra de color
hcb = colorbar('eastoutside', 'ytick', ctick, 'yticklabel', ctick, 'ylim', clim);
set(get(hcb, 'ylabel'), 'String', '$20\log_{10}|H(\sigma+j 2\pi f)|$ (dB)', 'interpreter', 'latex');
view([0 90])       % Vista de la gráfica desde arriba
grid on
box on

% Gráfica de la magnitud de la función de transferencia en función de la frecuencia
figure
plot(f(Nh/2:end)*1e-3, 20*log10(abs(Hw)), 'b', 'linewidth', 2)
xlabel('$f$ (kHz)', 'interpreter', 'latex')         
ylabel('$20\log_{10}|H(f)|$', 'interpreter', 'latex')
title('Funcion de Transferencia');
set(gca, 'xscale', 'log')       
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
xlim([.1 20])           
ylim([-60 3])
grid on
box on

% Gráfica de la fase de la función de transferencia en función de la frecuencia
figure
plot(f(Nh/2:end)*1e-3, unwrap(angle(Hw)), 'b', 'linewidth', 2)
xlabel('$f$ (kHz)', 'interpreter', 'latex')   
ylabel('$\angle H(f)$', 'interpreter', 'latex')
title('Funcion de Transferencia');
set(gca, 'xscale', 'log')     
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
xlim([.1 20])           
grid on
box on

% Gráfica de la respuesta impulsiva en el dominio del tiempo
figure
plot(th*1e3, h, 'b', 'linewidth', 2)     
xlabel('$t$ (ms)', 'interpreter','latex')
ylabel('$h(t)$', 'interpreter','latex')
title('Respuesta Impulsiva');
axis tight                                
ylim([-1 1])
grid on
box on

% Gráfica de la señal de audio filtrada en función del tiempo
figure
plot(tx, y, 'b', 'linewidth', 2)            
xlabel('$t$ (segundos)', 'interpreter','latex')
ylabel('$y(t)$', 'interpreter','latex')
title('Audio Filtrado');
axis tight                                  
ylim([-1 1])
grid on
box on

% Gráfica de la magnitud de la transformada de Fourier de la señal de audio filtrada
figure
plot(fx*1e-3, 20*log10(abs(Y)), 'b', 'linewidth', 2)
xlabel('$f$ (kHz)', 'interpreter','latex')  
ylabel('$20\log_{10}|Y(f)|$ (dB)', 'interpreter','latex')
title('Audio Filtrado');
set(gca, 'xscale', 'log')                   
set(gca, 'xtick', 2.^(-2:4), 'xticklabel', 2.^(-2:4))
xlim([.1 20])                               
grid on
box on
%% Exportación de Figuras
%figName = 'FigAudio';
%exportgraphics(gcf, [figFolder, figName, figNameExt], 'Resolution', figRes)
