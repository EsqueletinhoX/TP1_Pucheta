clear; close all; clc;

%% Parámetros del motor
Laa = 366e-6;
J = 5e-9;
Ra = 55.6;
B = 0;
Ki = 6.49e-3;
Km = 6.53e-3;
Tl = 0;
Ve = 12;

%% Condiciones iniciales
ia = 0;
w = 0;
theta = 0;

%% Simulación
dt = 1e-7; % Paso de integración de Euler
t_end = 5; % tiempo final
N = round(t_end/dt); % num de pasos

ia_v = zeros(1, N);
w_v = zeros(1, N);
t = (0:N-1)*dt;

for k=1:N
    % Ec. dif
    dia = (-Ra*ia - Km*w + Ve)/Laa;
    dw = (Ki*ia - Tl)/J;

    % integracion por Euler
    ia = ia + dia*dt;
    w = w + dw*dt;

    % Almacenamiento
    ia_v(k) = ia;
    w_v(k) = w;
end

%% Gráficas
figure;
subplot(2,1,1);
plot(t, ia_v, 'r', 'LineWidth', 1.2);
xlabel('Tiempo [s]'); ylabel('Corriente ia [A]');
title('Corriente del motor'); grid on;

subplot(2,1,2);
plot(t, w_v, 'b', 'LineWidth', 1.2);
xlabel('Tiempo [s]'); ylabel('Velocidad angular \omega [rad/s]');
title('Velocidad angular del motor'); grid on;

%% Resultados finales
fprintf('\nMáxima corriente: %.3f A\n', max(ia_v));
fprintf('Máxima velocidad angular: %.2f rad/s\n', max(w_v));
fprintf('Torque máximo desarrollado: %.5f Nm\n', (Ki * max(ia_v)));