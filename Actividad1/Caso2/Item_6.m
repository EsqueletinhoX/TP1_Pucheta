clear; close all; clc;

%--------------------------------------------
% Item 6 Trabajo Practico Control 2
%--------------------------------------------

%% PARÁMETROS DEL MOTOR
Ra = 2.413563466549881;
Km = 0.262517625619119;
Ki = 0.276855612210369;
Bm = 0;
J = 0.002845672285857;
La = 0.005486068069999;

%% MODELADO DEL MOTOR EN ESPACIO DE ESTADOS
As = [-Ra/La -Km/La 0 ; Ki/J -Bm/J 0 ; 0 1 0];
Bs = [1/La 0 ; 0 -1/J ; 0 0];
Cs = [0 1 0];
Ds = [0 0];

ia(1) = 0;
wr(1) = 0;
theta(1) = 0;
x = [ia(1); wr(1); theta(1)];  % Condiciones iniciales: [ia; w; theta]
w_ref = 1;      % Wr esperado

% Funciones de transferencia del sistema resueltas con espacio de estados
% Wr/Va
[num1, den1] = ss2tf(As, Bs, Cs, Ds, 1);
G_wv = tf(num1, den1);

% Wr/Tl
[num2, den2] = ss2tf(As, Bs, Cs, Ds, 2);
G_wt = tf(num2, den2);

%% PARÁMETROS DE SIMULACIÓN
p = pole(G_wv);             % Polos de la funcion de transferencia
tR = log(0.95)/p(2);        % Obtengo la dinámica rapida con el polo más
                            % alejado del eje imaginario
Ts = tR/10;                 % Tiempo de muestreo mucho menor que tR
tF = 30;                   % t final
t = 0:Ts:tF;                % vector tiempo
tL = zeros(1, length(t));   % vector torque
TL = 0.12;                  % torque max aplicado

tL(round(0.701/Ts):round(1/Ts)) = TL;

% Planteo el vector para el error
e = zeros(1, length(t));
u = 0;

%% IMPLEMENTACIÓN DEL PID DISCRETO
%u(i+1) = u(i) + A*e(i) + B*e(i-1) + C*e(i-2)

% Parámetros PID
% Kp = 0.1; KI = 0.01; Kd = 5; % Pruebo recomendación del ejercicio
% Kp = 1; KI = 1; Kd = 5;
% Kp = 10; KI = 3; Kd = 1;
Kp = 0.1; KI = 0.1; Kd = 0.00005;

A = (2*Kp*Ts + KI*Ts^2 + 2*Kd)/(2*Ts);
B = (-2*Kp*Ts + KI*Ts^2 - 4*Kd)/(2*Ts);
C = Kd/Ts;

% Prelocación de vectores para asegurar mayor velocidad
ia = zeros(1, length(t));
w = zeros(1, length(t));
theta = zeros(1, length(t));
acc = zeros(1, length(t));

for i=1:length(t)-1
    e(i) = w_ref - theta(i);    %Error
    
    if i == 1
        u = u + A*e(i);
    end
    if i == 2
        u = u + A*e(i) + B*e(i-1);
    end
    if i>2
        u = u + A*e(i) + B*e(i-1) + C*e(i-2);
    end

    % Integración por método de Euler
    dx = As*x + Bs*[u; tL(i)];
    x = x + dx*Ts;

    % Guardar resultados
    ia(i+1) = x(1);
    w(i+1) = x(2);
    theta(i+1) = x(3);
    acc(i+1) = u;
end

%% PLOTS
fz = 12;
figure;

subplot(5,1,1); hold on; grid on;
plot(t,ia,'r'); title('Salida, ia', 'Interpreter','latex','FontSize',fz);
ylabel('I[A]', 'Interpreter','latex','FontSize',fz);

subplot(5,1,2); hold on; grid on;
plot(t,w,'g'); title('Salida, wr', 'Interpreter','latex','FontSize',fz);
ylabel('wr[rad/seg]', 'Interpreter','latex','FontSize',fz);

subplot(5,1,3); hold on; grid on;
plot(t,theta,'b'); title('Salida, Angulo', 'Interpreter','latex','FontSize',fz);
ylabel('theta [rad]', 'Interpreter','latex','FontSize',fz);

subplot(5,1,4); hold on; grid on;
plot(t,acc,'g'); title('Entrada, Va', 'Interpreter','latex','FontSize',fz);
ylabel('V[V]', 'Interpreter','latex','FontSize',fz);

subplot(5,1,5); hold on; grid on;
plot(t,tL,'Color','k'); title('Entrada, Torque', 'Interpreter','latex','FontSize',fz);
ylabel('TL[Nm]', 'Interpreter','latex','FontSize',fz);