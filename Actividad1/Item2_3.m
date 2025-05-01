clear; close all; clc;
data = readmatrix('Curvas_Medidas_RLC_2025.xls');
t = data(:,1);   % Tiempo
i = data(:,2);   % Corriente
vc = data(:,3);  % Voltaje en el capacitor
u = data(:,4);  % Entrada
vo = data(:,5);  % Salida

% figure;
% 
% subplot(4,1,1);
% plot(t, u, 'k'); title('Entrada u(t)'); ylabel('Voltaje [V]'); grid on;
% 
% subplot(4,1,2);
% plot(t, i, 'r'); title('Corriente i(t)'); ylabel('Corriente [A]'); grid on;
% 
% subplot(4,1,3);
% plot(t, vc, 'b'); title('Tensión en el capacitor v_C(t)'); ylabel('Voltaje [V]'); grid on;
% 
% subplot(4,1,4);
% plot(t, vo, 'm'); title('Tensión en el resistor v_o(t)'); ylabel('Voltaje [V]'); xlabel('Tiempo [s]'); grid on;

delay = 10e-3; % Lo busco en la tabla de datos (busco donde empieza la entrada)

Ts = t(2) - t(1);     % Paso de muestreo

%% === ÍTEM 2: Identificación de R, L y C con método de Chen ===

% Parámetros
StepAmplitude = 12;      % escalón de 12 V
K = 1;                   % Ganancia (puede ajustarse si la salida es escalada)

% Elegir t1 (en segundos) y calcular t2 = 2t1, t3 = 3t1
t1 = 3e-3;
t2 = 2*t1;
t3 = 3*t1;

% Convertir tiempos a índices
idx_delay = round(0.01 / Ts);  % Retardo de 10 ms
idx1 = round(t1 / Ts) + idx_delay;
idx2 = round(t2 / Ts) + idx_delay;
idx3 = round(t3 / Ts) + idx_delay;

% Valores de la señal en esos tiempos
y1 = vc(idx1);
y2 = vc(idx2);
y3 = vc(idx3);

% Cálculo de k1, k2, k3
k1 = (y1 / StepAmplitude) / K - 1;
k2 = (y2 / StepAmplitude) / K - 1;
k3 = (y3 / StepAmplitude) / K - 1;

% Coeficientes para el método de Chen
be = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;

alfa1 = (k1*k2 + k3 - sqrt(be)) / (2*(k1^2 + k2));
alfa2 = (k1*k2 + k3 + sqrt(be)) / (2*(k1^2 + k2));
beta  = (k1 + alfa2) / (alfa1 - alfa2);

% Constantes de tiempo
T1 = -t1 / log(alfa1);
T2 = -t1 / log(alfa2);
T3 = beta * (T1 - T2) + T1;

% Función de transferencia identificada
G_identificada = tf(K * [T3 1], conv([T1 1], [T2 1]));

%% === Comparación con modelo teórico RLC ===

% Comparar con: G(s) = 1 / (L*C*s^2 + R*C*s + 1)
% G_identificada = (T3*s + 1) / (T1*T2*s^2 + (T1+T2)*s + 1)

% Igualando coeficientes:
% L*C = T1*T2
% R*C = T1 + T2

% Supongamos K = 1 como ya definimos
LC = T1 * T2;
RC = T1 + T2;

% Elegimos un C arbitrario (ej. 2.2e-6 F) y calculamos R y L
% C = 2.2e-6;
% R = RC / C;
% L = LC / C;

% Buscar el primer escalón positivo (umbral de cambio)
idx_start = find(u > 0.9 * StepAmplitude, 1, 'first');  % Primer índice donde empieza el escalón positivo
idx_end   = idx_start + round(10e-3 / Ts);              % 10 ms de duración del escalón

% Corriente máxima durante ese primer escalón positivo
Imax = max(i(idx_start:idx_end));

R = StepAmplitude/Imax; 
C = RC/R;
L = LC/C;

fprintf('\n--- Parámetros identificados ---\n');
fprintf('C = %.3e F\n', C);
fprintf('R = %.3f ohm\n', R);
fprintf('L = %.3e H\n', L);

%% === ÍTEM 3: Validación con corriente ===

% Modelo en espacio de estados
A = [-R/L, -1/L; 1/C, 0];
B = [1/L 0]';
x = [0 0]';  % condiciones iniciales

% Simulación de corriente y Vc
I_sim = zeros(size(t));
Vc_sim = zeros(size(t));

for k = 1:length(t)
    I_sim(k) = x(1);
    Vc_sim(k) = x(2);
    dx = A * x + B * u(k);
    x = x + dx * Ts;
end

% Comparación a partir de 0.05 s
idx_ini = round(0.05 / Ts);

figure;
plot(t(idx_ini:end), I_sim(idx_ini:end), 'b', 'LineWidth', 1.4); hold on;
plot(t(idx_ini:end), i(idx_ini:end), 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Tiempo [s]');
ylabel('Corriente [A]');
title('Comparación de corriente medida vs simulada');
legend('Corriente medida','Corriente simulada');

figure;
plot(t(idx_ini:end), vc(idx_ini:end), 'r--', 'LineWidth',1.5); hold on;
plot(t(idx_ini:end), Vc_sim(idx_ini:end), 'b', 'LineWidth',1.4); hold on;
plot(t(idx_ini:end), u(idx_ini:end), 'g', 'LineWidth',1.4);
grid on;
xlabel('Tiempo [s]');
ylabel('Tensión [V]');
title('Comparación de tensión en el capacitor medida vs simulada');
legend('Tensión medida','Tensión simulada');