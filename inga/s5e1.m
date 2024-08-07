clc, clear all, close all,
A = [0   1
     -2 -3];
B = [0 
     2];
% Diseño del controlador
% Paso 1: Controlabilidad
Co = [B A*B];
if rank(Co)==2
    disp('El sistema es controlable')
else
    disp('El sistema no es controlable :(')
end
% Paso 2: Polinomio original
polos_og = eig(A);
poly_og = poly(polos_og);
a1 = poly_og(1,2);
a2 = poly_og(1,3);
% Paso 3: W
W = [a1 1
      1 0];
% Paso 4: T
T = Co*W;
% Paso 5: Nuevo polinomio
e = 0.6;
wn = 13.5;
polos_new = [-e*wn+j*wn*sqrt(1-e^2) -e*wn-j*wn*sqrt(1-e^2)];
poly_new = poly(polos_new);
h1 = poly_new(1,2);
h2 = poly_new(1,3);
% Paso 6: Kz
Kz = [h2-a2 h1-a1];
% Paso 7: K
K = Kz*inv(T);
% Tiempo de simulación
ti = 0; dt = 0.001; tf = 3;
% Discretización
[Ak Bk] = c2d(A,B,dt);
% Condiciones iniciales
x = [5
     -10];
k = 1;
% Bucle de simulación
for tt = ti:dt:tf
    x1(k,1) = x(1);
    x2(k,1) = x(2);
    t(k,1) = tt;
    u = -K*x;
    U(k,1) = u;
    x = Ak*x + Bk*u;
    k=k+1;
end
figure(1)
subplot(2,2,1); plot(t,x1)
subplot(2,2,3); plot(t,x2)
subplot(2,2,2); plot(t,U)

