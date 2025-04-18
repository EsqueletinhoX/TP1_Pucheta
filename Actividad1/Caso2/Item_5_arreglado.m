clear; close all; clc
data = readmatrix("Curvas_Medidas_Motor_2025_v.xls");
t = data(:,1);
w = data(:,2);
ia = data(:,3);
Va = data(:,4);
TL = data(:,5);

% figure;
% subplot(4,1,1);
% plot(t, Va); title('Tensión [V]'); grid on; hold on
% subplot(4,1,2)
% plot(t, ia); title('Corriente ia [A]'); grid on; hold on
% subplot(4,1,3)
% plot(t, w); title('Velocidad angular \omega[rad/seg]'); grid on; hold on
% subplot(4,1,4)
% plot(t, TL); title('Torque [Nm]'); grid on; hold on

%% METODO DE CHEN 1
StepAmplitude=2; %2 V de entrada en Va

ret = 0.1;
ts = 0.11 - ret;



% %wr/va
% y1 = w(150)/StepAmplitude;
% t1 = t(150); %t1
% 
% y2 = w(200)/StepAmplitude;
% t2 = t(200); %t2
% 
% y3 = w(250)/StepAmplitude;
% t3 = t(250); %t3
% 
% % Tiempo entre muestras
% ts = t2 - t1;

% K=y(00)/U
k       =       w(end)/StepAmplitude;

%Funcion de la forma G(s)=K*(s*T3)/[(s*T1+1).(s*T2+1)] luego se puede
%despreciar el cero
k1      =       y1/k-1;
k2      =       y2/k-1;
k3      =       y3/k-1;
b       =       4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1   =       (k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2   =       (k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta    =       (2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1      =       (-ts/log(alfa1));
T2      =       (-ts/log(alfa2));

T1 = real(T1)
T2 = real(T2) % Importa solo la parte real
T3 = real(beta*(T1-T2)+T1);

G_wv = tf(k,conv([T1 1],[T2 1])) %funcion de transferencia W/Va

% G_wv =
% 
%               3.818
%   -----------------------------
%   0.0001273 s^2 + 0.04898 s + 1

%% PARÁMETROS FÍSICOS DEL MOTOR
% Formas generales de las funciones de transferencia
% G_wv_general = 
%                        Ki                      
%   ────────────────────────────────────────────
%   Jm*La*s^2 + (Ra*J + La*B)*s + (Ra*B + Ki*Km)

% Ra*J/(km*ki) = 0.04898;
% Jm*La/(km*ki) = 0.0001273;
% Jm*La/(km*ki) = Ra*0.0001273/0.04898;

% Wr/Va (Normalizado y considerando B=0) = 
%                       Ki                      
%   ──────────────────────────────────────────────────
%   Ki*Km*((Jm*La)/(Ki*Km)*s^2 + (Ra*Jm)/(Ki*Km)*s +1)

% Sabemos que Ra=Va/Imax, donde del gráfico de las curvas tomamos:
Imax = max(abs(ia));
Ra = 2/Imax;
% Considerando la ganancia de Wr/TL podemos obtener Ki y Km de las
% funciones de transferencia normalizadas
K_tl = (3.63973774664743 - 7.62473632022227)/0.12;
Km = 1/k;
Ki = -Ra/(K_tl*Km);
Bm = 0;
J = (0.106*Ki*Km)/Ra;
La = (0.001316*Ki*Km)/J;

% Imax = max(abs(ia));
% Ra = 1.9/Imax;
% % Considerando la ganancia de Wr/TL podemos obtener Ki y Km de las
% % funciones de transferencia normalizadas
% K_tl = (3.63973774664743 - 7.62473632022227)/0.12;
% Km = 1/k;
% Ki = -Ra/(K_tl*Km);
% Bm = 0;
% J = 0.88*(0.106*Ki*Km)/Ra;
% La = 0.1*(0.001316*Ki*Km)/J;

%% SIMULACIÓN DEL FUNCIONAMIENTO A PARTIR DE LOS PARÁMETROS
% Modelado en el espacio de estados
A = [-Ra/La -Km/La 0 ; Ki/J -Bm/J 0 ; 0 1 0];
Bmat = [1/La 0 ; 0 -1/J ; 0 0];
C = [0 1 0];
D = [0 0];

x = [0; 0; 0]; % estado inicial: [ia; w; theta]

dt = 0.001;

time = 0:dt:1.5-dt;

lsim = length(time);

u = zeros(lsim,1);
u(round(0.102/dt)+1:end) = 2;

Tl = zeros(lsim,1);
Tl(round(0.701/dt):1/dt) = 0.12;

ia_sim = zeros(lsim,1);
w_sim = zeros(lsim,1);
theta_sim = zeros(lsim,1);

for i = 1:lsim-1

    u_sim = u(i);
    Tl_sim = Tl(i);

    % Integración por método de Euler
    dx = A * x + Bmat * [u_sim; Tl_sim];
    x = x + dx * dt;

    % Guardar resultados
    ia_sim(i+1) = x(1);
    w_sim(i+1) = x(2);
    theta_sim(i+1) = x(3);
end

figure;
subplot(2,1,1);
plot(t, w, 'g', t, w_sim, 'k--'); grid on;
title('Velocidad angular \omega [rad/s]');
legend('Medida', 'Simulada');

subplot(2,1,2);
plot(t, ia, 'g', t, ia_sim, 'k--'); grid on;
title('Corriente de armadura i_a [A]');
legend('Medida', 'Simulada');

