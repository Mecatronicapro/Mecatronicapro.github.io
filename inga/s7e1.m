clc, clear all, close all,

R = 5;
L = 0.2e-3;
a = pi/4;
kt = 0.05;
kb = 0.05;
p = 0.002;
r = 0.01;
J = 1e-6;
m = 0.5;
c = 20;
d = J*(2*pi)/p + m*r*tan(a);
A = [0 1 0
     0 -c*r*tan(a)/d kb/d
     0 -kb*2*pi/(p*L) -R/L];
B = [0
    0
    1/L];
E = [0
    -r/d*tan(a)
    0];
umax = 24;
%  Tiempo de simulación
ti = 0; dt = 0.0001;tf = 2;
% Discretización
[Ak Bk] = c2d(A,B,dt);
[Ak Ek] = c2d(A,E,dt);
x = [0 0 0]';
k = 1;
fs = 0;
for tt = ti:dt:tf

    x1(k,1) = x(1);
    x2(k,1) = x(2);
    x3(k,1) = x(3);
    t(k,1) = tt;
    u = 5;
    % saturación de la entrada
    if u>umax
        u = umax;
    elseif u<-umax 
        u = -umax;
    end
    v(k,1) = u;
    % Dirección de la fricción seca
    Fs(k,1) = fs;
    if x(2)>0
        w = fs;
    elseif x(2)<0
        w = -fs;
    else
        w = 0;
    end
    x = Ak*x + Bk*u + Ek*w;
    k = k + 1;
end

figure(1)
subplot(3,2,1)
plot(t,x1)
title('x1 - posición')
subplot(3,2,3)
plot(t,x2)
title('x2 - velocidad')
subplot(3,2,5)
plot(t,x3)
title('x3 - corriente')
subplot(3,2,2)
plot(t,v)
title('Tensión de entrada')
subplot(3,2,4)
plot(t,Fs)
title('Fricción seca')