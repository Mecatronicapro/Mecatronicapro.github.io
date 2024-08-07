clc, clear all, close all

% Parámetros
a = (4*(2.54)/100)^2*pi/4;  % Área de la salida (m^2)
Dt = 10;                    % Diámetro del tanque (m)
At = pi/4*Dt^2;             % Área transversal del tanque (m^2)
g = 9.81;                   % Aceleración debido a la gravedad (m/s^2)
u = 0.9;                    % Coeficiente de descarga
hmax = 12;                  % Altura máxima del tanque (m)

% Tiempo de simulación
ti = 0; dt = 0.1; tf = 5;   % Tiempo inicial, paso de tiempo, tiempo final

% Puntos de equilibrio
he = 0.8*hmax;              % Nivel de equilibrio del agua en el tanque
Qie = a*u*sqrt(2*g*he);     % Caudal de equilibrio

A = -a*u/(2*At)*sqrt(2*g/he);  % Parámetro del sistema linealizado
B = 1/At;                      % Parámetro del sistema linealizado

% Discretización
[Ak, Bk] = c2d(A,B,dt);     % Discretización del sistema continuo

% Condiciones iniciales
h = he;                     % Nivel inicial del agua
k = 1;
var_h = h - he;             % Desviación inicial del nivel
hr = h;                     % Nivel real inicial del agua

% Bucle de simulación
for tt = ti:dt:tf
    niv(k,1) = h;           % Guardar el nivel de agua aproximado
    niv_re(k,1) = hr;       % Guardar el nivel de agua real
    t(k,1) = tt;            % Guardar el tiempo
    Qi = 0.5;               % Caudal de entrada constante
    
    % Sistema Linealizado
    var_Qi = Qi - Qie;      % Desviación del caudal de equilibrio
    cau(k,1) = Qi;          % Guardar el caudal de entrada
    var_h = Ak * var_h + Bk * var_Qi;  % Calcular la nueva desviación del nivel
    h = var_h + he;         % Calcular el nuevo nivel de agua aproximado
    
    % Sistema Real
    hr = hr + (Qi - a*u*sqrt(2*g*hr)) / At * dt;  % Calcular el nuevo nivel de agua real
    k = k + 1;
end

figure(1)
subplot(2,1,1)
plot(t,niv,t,niv_re)
ylabel('Nivel [m]')
legend('Lineal aprox.','Real')
grid on
subplot(2,1,2)
plot(t,cau)
ylabel('Caudal de entrada [m3/s]')
xlabel('Tiempo [s]')
grid on
