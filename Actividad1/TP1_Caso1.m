clear; close all; clc;
% Parámetros del circuito
R = 220;              % Resistencia en ohmios
L = 500e-3;           % Inductancia en henrios
Cap = 2.2e-6;           % Capacitancia en faradios

% Consigna

% Asignar valores a R=47, L=1e-6, y C=100e-9. Obtener simulaciones que permitan 
% estudiar la dinamica del sistema, con una entrada de tensión escalón de 12V, que cada 1ms cambia 
% de signo.

% Matrices que representan el sistema

A = [-R/L -1/L; 1/Cap 0];

B = [1/L 0]';

C = [R 0];

D = 0;

x = [0 0]'; % condiciones iniciales nulas


%% Analisis del sistema

% Funcion de transferencia Vr/Ve = 1/(s^2*L/R + s + 1/R/C)
% obtenida de forma analitica

p_v = roots([L/R 1 1/R/Cap]);
% Polos complejos conjugados
p1 = p_v(1);
p2 = p_v(2);

% p1 = -2.2000e+02 + 9.2773e+02i
% p2 = -2.2000e+02 - 9.2773e+02i
% Calculo el tiempo Ts y lo divido por 10 para calcular la dinamica del
% sistema

Ts = (2*pi)/9.2773e2;
Ts = Ts/500;

%% Vectores de tiempo para cada entrada de referencia

delay = Ts*1e2; % retardo inicial

t0_v = 0:Ts:delay; 

t1_v = delay+Ts:Ts:10e-3+delay;

t2_v = t1_v(end)+Ts:Ts:t1_v(end)+10e-3;

t3_v = t2_v(end)+Ts:Ts:t2_v(end)+10e-3;

t4_v = t3_v(end)+Ts:Ts:t3_v(end)+10e-3;

time = [t0_v t1_v t2_v t3_v t4_v];

%% Vectores con los cambios en la entrada de referencia

u1 = zeros(1,length(t0_v)); % Delay

u2 = 12*ones(1, length(t1_v)); % Entrada +12V

u3 = -12*ones(1, length(t2_v)); % Entrada -12V

u = [u1 u2 u3 u2 u3];

Vc = zeros(1, length(time));
I = zeros(1, length(time));
Vo = zeros(1, length(time));

for idx = 1:length(time)
    
    I(idx) = x(1);
    Vc(idx) = x(2);
    Vo(idx) = R*I(idx);
    
    xp = A*x + B.*u(idx);
    
    x = x + xp.*Ts;
    
end

%% Plots

fz = 12;

figure;

subplot(4,1,1);
plot(time.*1e3, u, 'r', 'LineWidth', 1.1);
grid on;
title('Entrada','Interpreter','latex','FontSize',fz);
ylabel('V[V]','Interpreter','latex','FontSize',fz);
xlabel('t[ms]','Interpreter','latex','FontSize',fz);

subplot(4,1,2);
plot(time.*1e3, Vc, 'LineWidth', 1.1);
grid on;
title('Voltaje en el Capacitor','Interpreter','latex','FontSize',fz);
ylabel('V[V]','Interpreter','latex','FontSize',fz);
xlabel('t[ms]','Interpreter','latex','FontSize',fz);

subplot(4,1,3);
plot(time.*1e3, I, 'LineWidth', 1.1);
grid on;
title('Corriente','Interpreter','latex','FontSize',fz);
ylabel('I[A]','Interpreter','latex','FontSize',fz);
xlabel('t[ms]','Interpreter','latex','FontSize',fz);

subplot(4,1,4);
plot(time.*1e3, Vo, 'LineWidth', 1.1);
grid on;
title('Salida: Voltaje en el resistor','Interpreter','latex','FontSize',fz);
ylabel('V_o[V]','Interpreter','latex','FontSize',fz);
xlabel('t[ms]','Interpreter','latex','FontSize',fz);

[x1,x2,x3,x4,x5] = importfile(Curvas_Medidas_RLC_2025.xls);

