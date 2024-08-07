clc, clear all, close all,
% Parámetros
m = 4;
k = 50;
b = 20;
% Sistema en espacio de estados
A = [-b/m -k/m
      1    0];
B = [1/m
      0];
C = [0 1];
D = [0];
% Tiempo de simulación
ti = 0; dt = 0.001; tf = 2;
% Discretización
[Ak Bk] = c2d(A,B,dt);
% Condiciones iniciales
x = [0 
     0];
f = 20;
k = 1; % variable de control de la simulación
% Bucle de simulación
for tt=ti:dt:tf
    x1(k,1) = x(1);
    x2(k,1) = x(2);
    t(k,1) = tt;
    f = 4*sin(2*pi*10*tt);
    u(k,1) = f;
    y = C*x + D*f;
    sensor(k,1) = y(1);
    x = Ak*x + Bk*f;
    k = k + 1;
end
% Gráficas
figure(1)
subplot(2,2,1)
plot(t,x1)
ylabel('x1 [m/s]');title('Estado de Velocidad (xp)')
grid on
subplot(2,2,3)
plot(t,x2)
ylabel('x2 [m]');title('Estado de Posición (x)')
xlabel('t [s]')
grid on
subplot(2,2,2)
plot(t,u)
ylabel('u [N]');title('Señal de entrada (f)')
grid on
subplot(2,2,4)
plot(t,sensor)
ylabel('sensor [m]');title('Señal del sensor (f)')
grid on
xlabel('t [s]')