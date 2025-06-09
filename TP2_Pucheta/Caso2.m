% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Baccino, Luca
% Tp N° 2 - Caso de estudio 1 - 
%   Item 1 - Inciso 2
%   Objetivo: controlar el motor siendo posible solamente conocer la
%   velocidad y el ángulo pero no la corriente
%
%   Item 1 - Inciso 3
%   Objetivo: determinar el efecto de la nolinealidad en la acción de 
%   control y verificar cuál es el máximo valor admisible de ésa no linealidad.
%%
clc; clear; close all
data = readmatrix("Curvas_Medidas_Motor_2025_v");
t = data(:,1);
w = data(:,2);
ia = data(:,3);
Va = data(:,4);
TL = data(:,5);

%% Gráficos de los datos medidos
% figure;
% subplot(4,1,1);
% plot(t, Va); title('Tensión [V]'); grid on; hold on
% subplot(4,1,2)
% plot(t, ia); title('Corriente ia [A]'); grid on; hold on
% subplot(4,1,3)
% plot(t, w); title('Velocidad angular \omega[rad/seg]'); grid on; hold on
% subplot(4,1,4)
% plot(t, TL); title('Torque [Nm]'); grid on; hold on

%% 
% PASO 1: MODELADO EN ESPACIO DE ESTADOS PARA EL SIST. CONTINUO
% Los parámetros del motor de CC son los que se obtuvieron en el TP1
% haciendo uso del método de Chen
StepAmplitude = 2;
k = w(end)/StepAmplitude;
Ra = 2.27;
Km = 0.25;
Ki = 0.25;
Bm = 0.00131;
J = 0.00233;
La = 0.0047;

% Variables de estado:
% x1 = ia ; x2 = wr ; x3 = theta 

% Salidas (las que se pueden medir según la consigna):
% y1 = wr ; y2 = theta

% Entradas:
% u1 = Va ; u2 = TL

% ESPACIO DE ESTADOS
% Matriz de estados
Ac = [-Ra/La -Km/La 0;
        Ki/J -Bm/J 0;
        0 1 0];

% Matriz de entrada
Bc = [1/La 0;
      0 -1/J;
      0    0];

% Matriz de salida (wr, theta)
Cc = [0 1 0;
    0 0 1];

% Matriz de transmisión directa
Dc = [0 0;
     0 0];

% Sistema continuo en Ve a lazo abierto
sysCont = ss(Ac,Bc,Cc,Dc);

%% PASO 2: CONVERSIÓN A TIEMPO DISCRETO
% Cálculo de la dinámica del sistema
val = real(eig(Ac));

% Dinámica rápida
tR = log(0.95)/val(3);
tInt = tR/4; 

% Dinámica lenta
tL = log(0.05)/val(2);
tsim = tL*5;

% Tiempo de muestreo
Ts = tInt;

% DEFINICION DEL SISTEMA DISCRETO
sysDisc = c2d(sysCont, Ts, 'zoh');
Ad = sysDisc.a; % Matriz de estados
Bd = sysDisc.b; % Matriz de entrada considerando control sobre el torque
Bd_aux = Bd(:,1); % Matriz de entrada auxiliar sin considerar control sobre el torque
Cd = sysDisc.c; % Matriz de salida
Dd = sysDisc.d; % Matriz de acoplamiento directo

%% PASO 3: ANÁLISIS DE CONTROLABILIDAD Y ALCANZABILIDAD
% Controlabilidad
Mc = [Bd_aux Ad*Bd_aux Ad^2*Bd_aux Ad^3*Bd_aux];
controlabilidad = rank(Mc);

% Alcanzabilidad
Ma = [Bd_aux Ad*Bd_aux Ad^2*Bd_aux];
alcanzabilidad = rank(Ma);

if(alcanzabilidad == controlabilidad)
    fprintf('\nSistema controlable')
else
    fprintf('\nSistema NO controlable')
end

%% PASO 4: DEFINICIÓN DE LA LEY DE CONTROL A APLICAR
% Se usará una ley de control con integrador con la idea de anular el error
% estado estable del sistema, teniendo en cuenta una referencia no nula
% 
% u = -K*x + Ki*psi
% 

%% PASO 5: DEFINICIÓN DEL SISTEMA AMPLIADO
Aa = [Ad , zeros(3,1) ; -Cd(2,:)*Ad, eye(1)];
Ba = [Bd_aux; -Cd(2,:)*Bd_aux];

% Análisis de controlabilidad y alcanzabilidad del SD ampliado
Maa = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba];
Mca = [Ba Aa*Ba Aa^2*Ba Aa^3*Ba Aa^4];

if(rank(Mca) == rank(Maa))
    fprintf('\nSistema controlable')
else
    fprintf('\nSistema NO controlable')
end

%% PASO 6: IMPLEMENTACIÓN LQR
% CONSIDERACIONES PARA DEFINIR LAS MATRICES Q Y R
% 1) Impacto de la matriz Q:
%             Penaliza los estados del sistema, aumentar los elementos de la
%             diagonal que se corresponden a estados específicos hace que el
%             controlador trate de mantenerlos en 0.
% 2) Impacto de la matriz R:
%             Penaliza entradas de control, aumentar los valores de los elementos
%             de R hace que el controlador use menos energía de control.

% Como consideración para el análisis, se partirá de Q y R matrices
% identidad, se calcula el controlador y se simula el sist. a lazo cerrado.

% Aspectos importantes a tener en cuenta para evaluar los resultados:
%          - Velocidad de respuesta
%          - Amplitud máx de la acción de control
%          - Oscilaciones
% Para corregir los diferentes resultados:
%          - Rta. lenta -> Aumento Q, penalizando los estados
%          - Acción de control grande -> Aumento R, penalizando la acción
%                                        de control

% DEFINO Q Y R
d = [100 1 0.40528 0.1]; % Elementos de la diagonal de Q (ia, wr, theta, psita)
Q = diag(d);
R = 40e3;

% Cálculo del controlador
[Klqr, ~, ~] = dlqr(Aa, Ba, Q, R);
K = Klqr(1:3); % Ganancia para los estados originales
KI = -Klqr(4); % Ganancia para el integrador

%% PASO 7: DISEÑO DEL OBSERVADOR
% SISTEMA DUAL
Ao = Ad';
Bo = Cd';
Co = Bd_aux';

% DEFINICIÓN DE MATRICES Qo Y Ro
% Consideraciones:
% 1) Matriz Qo:
%       Matriz de covarianza del ruido del proceso. Representa la
%       incertidumbre en el modelo del sistema. Si los estados son
%       conocidos con precisión, los valores en Qo deben ser pequeños. Si
%       hay incertidumbre en los estados, los valores deben ser mayores.
% 2) Matriz Ro:
%       Matriz de covarianza del ruido de la medición. Representa la
%       incertidumbre en las mediciones de salida. Si las mediciones son
%       precisas, los valores en Ro deben ser pequeños. Si las mediciones
%       tienen ruido significativo, los valores deben ser mayores.

% Parto de matrices identidad y se varía hasta obtener los resultados que
% se esperan, se debe tener en cuenta que estados se conocen y cuales no.
do = [100 0.1 0.1]; % ia, wr, theta
Qo = diag(do);
Ro = [10   0;
      0   10];

% Cálculo del controlador del sistema dual
Ko = (dlqr(Ao, Bo, Qo, Ro))';

%% PASO 8: VERIFICACIÓN 
Tf = 10; % Tiempo total de simulación
dT = tInt; % Tiempo de integración
p_max = floor(Tf/dT);
deadZone = 0.5;

% Tiempo
tt = 0:dT:p_max*dT;

% Variables de estado
Ia = zeros(1, p_max+1);
Wr = zeros(1, p_max+1);
Theta = zeros(1, p_max+1);

% Torque
TL_v = zeros(1, p_max+1);

% Referencia de entrada
ref = zeros(1, p_max+1);

% Acción de control
u = zeros(1, p_max+1);

% INICIALIZACIÓN:

% Señal de referencia:
thetaRef = pi/2;
ref(1) = thetaRef;
iCounter = 0;
tSwitch = 5;
for i = 1:p_max
    iCounter = iCounter + dT;
    if(iCounter > tSwitch)
        thetaRef = -1*thetaRef;
        iCounter = 0;
    end
    ref(i) = thetaRef;
end

% Plot de la referencia para confirmar que es correcta
% plot(tt,ref,'b','LineWidth',1.5)

% Torque:
TLref = 0.12; 
TL_v(tt >= 0.7) = 0.12;
% plot(tt, TL_v, 'g','LineWidth',1.5)

% Inicialización de estados
x = [Ia(1) Wr(1) Theta(1)]'; % Vector de estados
xob = [0 0 0]';              % Vector de estados estimado
ei = 0;                      % Error de integración

%% SIMULACIÓN:
for i = 1:p_max

    y_out = Cd*x;       % Salida
    y_out_ob = Cd*xob;  % Salida observador

    ei = ei + (ref(i) - Cd(2,:)*x); % error

    % Acción de control
    u(i) = -K*xob + KI*ei;   % Con observador

    % Alinealidad del actuador
    if(abs(u) < deadZone)
        u(i) = 0;
    else
        u(i) = sign(u(i)) * (abs(u(i)) - deadZone);
    end

    Iap = -(Ra/La)*Ia(i) - (Km/La)*Wr(i) + (1/La)*u(i);
    Ia(i+1) = Ia(i) + Iap*dT;
    Wrp = (Ki/J)*Ia(i) - (Bm/J)*Wr(i) - (1/J)*TL_v(i);
    Wr(i+1) = Wr(i) + Wrp*dT;
    Theta(i+1) = Theta(i) + Wr(i)*dT;

    xob = Ad*xob + Bd_aux*u(i) + Ko*(y_out - y_out_ob);
    x = [Ia(i+1) Wr(i+1) Theta(i+1)]';
end

color = 'r';
figure(2)
subplot(3,2,1); hold on;
plot(tt,ref,'b','LineWidth',1.5);title('Angulo \phi [rad]'); grid on; hold on;
plot(tt, Theta, color, 'LineWidth',1.5);
subplot(3,2,2); hold on;
plot(tt, Wr, color, 'LineWidth',1.5);title('Velocidad angular \omega [rad/s]');grid on; hold on;
subplot(3,2,3); hold on;
plot(tt, Ia, color, 'LineWidth',1.5);title('Corriente Ia [A]');grid on; hold on;
subplot(3,2,4); hold on;
plot(tt, u, color, 'LineWidth',1.5);title('Acción de control V [V]');grid on; hold on;
subplot(3,1,3); hold on;
plot(tt, TL_v, color, 'LineWidth',1.5);title('Torque TL [N/m]');grid on; hold on;

figure(3)
plot(Theta, Wr, color, 'LineWidth', 1.5);title('Plano de Fases'); grid on; hold on;