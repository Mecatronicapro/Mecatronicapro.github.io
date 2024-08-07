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
% Diseño del controlador
% Paso 1: Co
Co = [B A*B A^2*B];
if rank(Co)==3
    disp('Es controlable :)')
else 
    disp('No es controlable :(')
end
% Paso 2: polinomio original
polos = eig(A);
pol_og = poly(polos);
a1 = pol_og(1,2);
a2 = pol_og(1,3);
a3 = pol_og(1,4);
% Pase 3: W
W = [a2 a1 1
     a1  1 0
      1  0 0 ];
% Paso 4: T
T = Co*W;
% Paso 5: polinomio deseado
polos_des = [-2 -2 -200];
pol_des = poly(polos_des);
h1 = pol_des(1,2);
h2 = pol_des(1,3);
h3 = pol_des(1,4);
% Paso 6: Kz
Kz = [h3-a3 h2-a2 h1-a1];
% Paso 7: K
K = Kz*inv(T);
% Diseño del obsevador
% Paso 1: Observabilidad
C = [1 0 0];
Ob = [C
      C*A
      C*A^2];
if rank(Ob)==3
    disp('El sistema es observable')
else
    disp('El sistem no es observable')
end
% Paso 2 y 3 ya se hicieron:
% Paso 4: Te
Te = W*Ob;
% Paso 5: Polinomio deseado
polos_ob = polos_des*124;
poly_ob = poly(polos_ob);
m1 = poly_ob(1,2);
m2 = poly_ob(1,3);
m3 = poly_ob(1,4);
% Paso 6: Kez
Kez = [m3-a3
       m2-a2
       m1-a1];
% Paso 7: Ke
Ke = inv(Te)*Kez;
Ao = A-Ke*C;
Bo = B;
Eo = Ke;
%  Tiempo de simulación
ti = 0; dt = 0.0001;tf = 20;
% Discretización
[Ak Bk] = c2d(A,B,dt);
[Ak Ek] = c2d(A,E,dt);
[Aok Bok] = c2d(Ao,Bo,dt);
[Aok Eok] = c2d(Ao,Eo,dt);
x = [0 0 0]';
xo = [0 0 0]';
r = [0.1 0 0]';
k = 1;
fs = 0;
int_ae = 0;
for tt = ti:dt:tf
    x1(k,1) = x(1);
    x2(k,1) = x(2);
    x3(k,1) = x(3);
    xo1(k,1) = xo(1);
    xo2(k,1) = xo(2);
    xo3(k,1) = xo(3);
    t(k,1) = tt;
    y = C*x;
    u = K*(r-x);
    error = r(1)-x(1);
    int_ae = int_ae + abs(error)*dt;
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
    xo = Aok*xo + Bok*u + Eok*y;
    k = k + 1;
end

IAE_RP_sAI = int_ae;
fprintf('IAE con RP sin AI %2.4f: \n',IAE_RP_sAI)
% Diseño del controlador con acción integrativa
Ai = [A zeros(3,1)
      1 0 0 0 ];
Bi = [B
      0];
% Diseño del controlador con AI
% Paso 1: Co
Co = [Bi Ai*Bi Ai^2*Bi Ai^3*Bi];
% Paso 2: polinomio original
polos = eig(Ai);
pol_og = poly(polos);
a1 = pol_og(1,2);
a2 = pol_og(1,3);
a3 = pol_og(1,4);
a4 = pol_og(1,5);
% Pase 3: W
W = [a3 a2 a1 1
     a2 a1  1 0
     a1  1  0 0 
      1  0  0 0];
% Paso 4: T
T = Co*W;
% Paso 5: polinomio deseado
polos_des = [-2 -2 -20 -200];
pol_des = poly(polos_des);
h1 = pol_des(1,2);
h2 = pol_des(1,3);
h3 = pol_des(1,4);
h4 = pol_des(1,5);
% Paso 6: Kz
Kz = [h4-a4 h3-a3 h2-a2 h1-a1];
% Paso 7: K
Ki = Kz*inv(T);
x = [0 0 0]';
k = 1;
int_e = 0;
for tt = ti:dt:tf
    x1i(k,1) = x(1);
    x2i(k,1) = x(2);
    x3i(k,1) = x(3);
    t(k,1) = tt;
    e = r(1)-x(1);
    int_e = int_e + e*dt;
    u = -Ki(1:3)*x + Ki(4)*int_e;
    % saturación de la entrada
    if u>umax
        u = umax;
    elseif u<-umax 
        u = -umax;
    end
    vi(k,1) = u;
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

% Diseño del controlador óptimo
% Paso 1: Controlabilidad
% Paso 2: Matrices Q y R
q1 = 100000;
q2 = 1;
q3 = 5;
Q = diag([q1 q2 q3]);
Q = [q1  0 0
      0 q2 0
      0  0 q3];
R = [1];
% Paso 3: Matriz P
P = are(A,B*inv(R)*B',Q);
% Paso 4: K
K = inv(R)*B'*P;
x = [0 0 0]';
xo = [0 0 0]';
r = [0.1 0 0]';
k = 1;
int_ae = 0;
for tt = ti:dt:tf
    x1_op(k,1) = x(1);
    x2_op(k,1) = x(2);
    x3_op(k,1) = x(3);
    xo1(k,1) = xo(1);
    xo2(k,1) = xo(2);
    xo3(k,1) = xo(3);
    t(k,1) = tt;
    y = C*x;
    u = K*(r-x);
    error = r(1)-x(1);
    int_ae = int_ae + abs(error)*dt;
    % saturación de la entrada
    if u>umax
        u = umax;
    elseif u<-umax 
        u = -umax;
    end
    v_op(k,1) = u;
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
    xo = Aok*xo + Bok*u + Eok*y;
    k = k + 1;
end
IAE_CO_sAI = int_ae;
fprintf('IAE con CO sin AI %2.4f: \n',IAE_CO_sAI)




figure(1)
subplot(3,2,1)
plot(t,x1,t,x1i,t,x1_op)
title('x1 - posición');legend('RP Sin AI','RP Con AI','CO Sin AI')
subplot(3,2,3)
plot(t,x2,t,x2i,t,x2_op)
title('x2 - velocidad');legend('RP Sin AI','RP Con AI','CO Sin AI')
subplot(3,2,5)
plot(t,x3,t,x3i,t,x3_op)
title('x3 - corriente');legend('RP Sin AI','RP Con AI','CO Sin AI')
subplot(3,2,2)
plot(t,v,t,vi,t,v_op)
title('Tensión de entrada');legend('RP Sin AI','RP Con AI','CO Sin AI')
subplot(3,2,4)
plot(t,Fs)
title('Fricción seca')

figure(2)
subplot(3,2,1)
plot(t,x1,t,xo1)
title('x1 - posición');legend('real','estimada')
subplot(3,2,3)
plot(t,x2,t,xo2)
title('x2 - velocidad');legend('real','estimada')
subplot(3,2,5)
plot(t,x3,t,xo3)
title('x3 - corriente');legend('real','estimada')
subplot(3,2,2)
plot(t,v)
title('Tensión de entrada');
subplot(3,2,4)
plot(t,Fs)
title('Fricción seca')