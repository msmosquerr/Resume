%% Recolección de datos

clear
simul=sim('PenduloForzadoKalman'); 
t = simul.t;
Ts=t(2)-t(1);

x = simul.x; % Se mide la posición

u = simul.u;
N = length(t); 
yteo =[simul.y] ; 
y_names ='y';

%% Gráfico de los datos experimentales

figure % Se debe usar para cada grupo de gráficos, lo cual evita ciertos conflictos

tiledlayout(1,2)

nexttile, plot(t,u), xlabel('t (seg)'), title('u(t)')
nexttile, plot(t,x), xlabel('t (seg)'), title('x(t)')

%% Parámetros del filtro

x0 = [0 0]';
n = length(x0); 
l= eye(n); 
PO=0*l;
Q=0.00001*l; 
R = 0.01^2;

%% Inicialización de las matrices

x_act = zeros(n,N);
P_act_traza = zeros(N,1);
K_norma = zeros(N,1); 
de_x = zeros(n,N); 
N_cond = zeros(N,1);
x_pre = x0;
P_pre = PO; 
% Valores estimados anteriormente (Practica 2 Alg. Gene.)
Delta=0.3413;
Alpha = 1.2045;
Beta = 2.1936;

C=[1 0]; % Valor a Estimar

%% Ecuaciones del filtro

for i = 1:N

% Corrección del estado

K = P_pre*C'/(C*P_pre*C' + R); 
x_act(:,i) = x_pre + K*(x(i) - C*x_pre);
P_act = (l-K*C)*P_pre;

% Simetrización

P_act = (P_act + P_act')/2;

% Actualización del estado para el siguiente paso (i+1)

x_pre = [x_act(1,i) + Ts*x_act(2,i); 
        x_act(2,i) + Ts*(u(i) - Delta*(x_act(2,i)) - Alpha*x_act(1,i) - Beta*(x_act(1,i))^3)]; % Discretización de Euler

A=[0 1; ...
  -Alpha-3*Beta*x(i) -Delta];%Linealización Matriz A

Fi=l + Ts*(A);
N_cond(i) = cond(obsv(Fi,C));
P_pre= Fi*P_act*Fi' + Q;


% Vectores con información:

P_act_traza(i) = trace(P_act); 
K_norma(i) = norm(K); 
de_x(:,i) = sqrt(diag(P_act));

end

%% Cálculo de los residuos

e_est = x-(C*x_act)';

%% Salida experimental vs Salida estimada

figure, plot(t,x,'k',t,x_act(1,:)','b',t,(x_act(1,:)+de_x(1,:))','r--',t,(x_act(1,:)-de_x(1,:))','r--')

xlabel('t (seg)'), title('Posición (X)')
legend({'Salida experimental' 'Salida estimada' 'Intervalo de confianza'})

%% Estados estimados vs Estados teóricos (no se muestran)

% nf = floor(sqrt(n)); nc = n-nf;

% tiledlayout(nf,nc)

% for i=1:n

% nexttile,
figure
hold on
plot(t,yteo(:), 'b')
plot(t,x_act(2,:)','k')
plot(t,(x_act(2,:)+de_x(2,:))','r--')
plot(t,(x_act(2,:)-de_x(2,:))','r--') 
xlabel('t(seg)')
title('Velocidad (y)')
legend({'Y Teorico' 'Y Estimado' 'Intervalo de confianza'})

% end

%% Ganancia de Kalman y matriz de covarianzas de los errores de los estados estimados

figure

tiledlayout(1,2)

nexttile, plot(t,K_norma), xlabel('t(seg)'), title('||K||')

nexttile, plot(t,P_act_traza), xlabel('t(seg)'), title('trace(P)')

%% Prueba de blancura de los errores de los residuos de la salida

figure, whiteness_test(e_est)

%% N_Condición
figure
plot(t,N_cond')
title('Numero Condición')
