clear; close all; clc
data = readmatrix("Curvas_Medidas_Motor_2025.xls");
t = data(:,1);
w = data(:,2);
ia = data(:,3);
Va = data(:,4);
TL = data(:,5);

figure;
subplot(4,1,1);
plot(t, Va); title('Tensión [V]'); grid on; hold on
subplot(4,1,2)
plot(t, ia); title('Corriente ia [A]'); grid on; hold on
subplot(4,1,3)
plot(t, w); title('Velocidad angular \omega[rad/seg]'); grid on; hold on
subplot(4,1,4)
plot(t, TL); title('Torque [Nm]'); grid on; hold on

N = length(t);
StepAmplitude = 12;
delay = 160.001;

%% MÉTODO DE CHEN 1
% G(s) = wr(s)/U(s)

% 1) Selecciono valores en base a la tabla
t1_w = t(201);
t2_w = t(202);
t3_w = t(203);

% Valores de señal ya normalizados (divido por StepAmplitude ya que el
% método de Chen está pensado para la respuesta a un escalón unitario)
y1_w = w(201)/StepAmplitude; 
y2_w = w(202)/StepAmplitude;
y3_w = w(203)/StepAmplitude;

% 2) Ganancia en estado estacionario
K_w = w(1161)/StepAmplitude;

% 3) Aplico el método
k1_w = (y1_w/K_w) - 1;
k2_w = (y2_w/K_w) - 1;
k3_w = (y3_w/K_w) - 1;

b_w     = 4*(k1_w^3)*k3_w - 3*(k1_w^2)*(k2_w^2) - 4*(k2_w^3) + (k3_w^2) + 6*k1_w*k2_w*k3_w; 
alfa1_w = (k1_w*k2_w + k3_w - sqrt(b_w))/(2*(k1_w^2 + k2_w)); 
alfa2_w = (k1_w*k2_w + k3_w + sqrt(b_w))/(2*(k1_w^2 + k2_w)); 
beta_w  = (k1_w+alfa2_w)/(alfa1_w-alfa2_w);

T1_w = -(t1_w - delay)/log(alfa1_w);
T2_w = -(t2_w - delay)/log(alfa2_w);
T3_w = beta_w * (T1_w - T2_w) + T1_w;

G_w = tf(K_w*[T3_w 1], conv([T1_w 1], [T2_w 1]));

% G_w =
% 
%     -0.0008979 s + 0.6364
% -----------------------------
% 0.0001036 s^2 + 0.09357 s + 1

% Puedo plotear G_w para verificar la dinámica del sistema

%% MÉTODO DE CHEN 2
% G_ia(s) = I(s)/U(s)

% 1) Selecciono valores en base a la tabla
t1_ia = t(202);
t2_ia = t(203);
t3_ia = t(204);

% Valores de señal ya normalizados (divido por StepAmplitude ya que el
% método de Chen está pensado para la respuesta a un escalón unitario)
y1_ia = ia(202)/StepAmplitude; 
y2_ia = ia(203)/StepAmplitude;
y3_ia = ia(204)/StepAmplitude;

% 2) Ganancia en estado estacionario
K_ia = ia(1191)/StepAmplitude;

% 3) Aplico el método
k1_ia = (y1_ia/K_ia) - 1;
k2_ia = (y2_ia/K_ia) - 1;
k3_ia = (y3_ia/K_ia) - 1;

b_ia     = 4*(k1_ia^3)*k3_ia - 3*(k1_ia^2)*(k2_ia^2) - 4*(k2_ia^3) + (k3_ia^2) + 6*k1_ia*k2_ia*k3_ia; 
alfa1_ia = (k1_ia*k2_ia + k3_ia - sqrt(b_ia))/(2*(k1_ia^2 + k2_ia)); 
alfa2_ia = (k1_ia*k2_ia + k3_ia + sqrt(b_ia))/(2*(k1_ia^2 + k2_ia)); 
beta_ia  = (k1_ia+alfa2_ia)/(alfa1_ia-alfa2_ia);

T1_ia = -(t1_w - delay)/log(alfa1_ia);
T2_ia = -(t2_w - delay)/log(alfa2_ia);
T3_ia = beta_ia * (T1_ia - T2_ia) + T1_ia;

G_ia = tf(K_w*[T3_ia 1], conv([T1_ia 1], [T2_ia 1]));

% G_ia =
% 
%       1.348 s + 0.6364
% -----------------------------
% 0.0001036 s^2 + 0.09357 s + 1

% Puedo plotear G_ia para verificar la dinámica del sistema

%% PARÁMETROS FÍSICOS DEL MOTOR
% Formas generales de las funciones de transferencia
% G_ia_general = 
%                   J*s + B
% ---------------------------------------------
% (La*J) s^2 + (Ra*J + La*B) s + (Ra*B + Ki*Km)

% G_w_general = 
%                      Ki
% ---------------------------------------------
% (La*J) s^2 + (Ra*J + La*B) s + (Ra*B + Ki*Km)

syms La Ra J B Ki Km real

den = [0.0001036 0.09357 1];
eq1 = La*J == den(1);
eq2 = Ra*J + La*B == den(2);
eq3 = Ra*B + Ki*Km == den(3);

num_G_ia = [1.348 0.6364];
eq4 = J == num_G_ia(1);
eq5 = B == num_G_ia(2);

num_G_w = [0.0008979 0.6364];
eq6 = Ki == num_G_w(2);

S = solve([eq1, eq2, eq3, eq4, eq5, eq6], [La, Ra, J, B, Ki, Km]);

La = double(S.La);
Ra = double(S.Ra);
J = double(S.J);
B = double(S.B);
Ki = double(S.Ki);
Km = double(S.Km);

%% SIMULACIÓN DEL FUNCIONAMIENTO A PARTIR DE LOS PARÁMETROS
% Modelado en el espacio de estados
A = [-Ra/La -Km/La 0 ; Ki/J -B/J 0 ; 0 1 0];
Bmat = [1/La 0 ; 0 -1/J ; 0 0];
C = [0 1 0];
D = [0 0];

x = [0 0 0]'; % x = [ia w theta]'

% Inicializo vectores con anterioridad para reducir el tiempo de ejecución
w_sim = zeros(1, length(t));
ia_sim = zeros(1, length(t));
theta_sim = zeros(1, length(t));

dt = 0.001;

for i=1:N-1
    u = Va(i);
    Tl = TL(i);

    dx = A*x + Bmat*[u ; Tl];
    x = x + dx*dt;
    ia_sim(i+1) = x(1);
    w_sim(i+1) = x(2);
    theta_sim(i+1) = x(3);
end

figure;
% Velocidad angular
subplot(2,1,1);
plot(t, w, 'g'); title('Velocidad angular \omega[rad/seg]'); hold on; 
plot(t, w_sim, 'k'); hold on;
legend({'\omega Medido', '\omega Simulado'}, 'Location','southeast')

% Corriente
subplot(2,1,2);
plot(t, ia, 'g'); title('Corriente [A]'); hold on; 
plot(t, ia_sim, 'k'); hold on;
legend({'I Medido', 'I Simulado'}, 'Location','southeast')

% Tensión de entrada
% subplot(4,1,3);
% plot(t, Va, 'g'); title('Tensión [V]'); hold on; 
% plot(t, u, 'k'); hold on;
% legend({'V Medido', 'V Simulado'}, 'Location','southeast')
% 
% Torque
% subplot(4,1,4);
% plot(t, TL, 'g'); title('Torque [Nm]'); hold on; 
% plot(t, Tl, 'k'); hold on;
% legend({'Torque Medido', 'Torque Simulado'}, 'Location','southeast')

disp(['Ra = ', num2str(Ra)])
disp(['La = ', num2str(La)])
disp(['Km = ', num2str(Km)])
disp(['Ki = ', num2str(Ki)])
disp(['J = ', num2str(J)])
disp(['B = ', num2str(B)])
