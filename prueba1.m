clc; clear all; close all;

%% Funciones unarias
% norma = A(x)/altura(x)
% concentracion = A^p , p>1
% dilatacion =  A^p , p<1
% dilatacion =  2A-A^2 , p<1
% intensificacion del contraste
% fuzzifiacion = difuminacion
%% Logica
% A -> B: -A OR B
% A <-> B: (-A OR B) AND (A OR -B)
%% Ejemplo
% Perdida de paquetes en la transmisión de datos
% Rendimiento de ancho de banda en la transmisión de datos
% Velocidad de conexión
% Nivel de congestión de red
% Nivel de riesgo de amenaza
% Latencia en la transmisión de datos
%% Concavo convexo
% lambda = 0.7;
% evaluar_valor_escalar(x, siluetaA, X1*lambda + (1-lambda)*X2);
% CONVEXO
% evaluar_valor_escalar <= max (A(X1) , A(X2))
% CONCAVO
% evaluar_valor_escalar >= min (A(X1) , A(X2))
%% Funciones unarias
% norma = A(x)/altura(x)
% concentracion = A^p , p>1
% dilatacion =  A^p , p<1
% dilatacion =  2A-A^2 , p<1
% intensificacion del contraste
% fuzzifiacion = difuminacion
%% Funciones transformar CB
% Operacion
% Metodo AND min
% Metodo OR max

%% Universo
x = linspace(0, 20, 500);

%% función de pertenencia
tria_func = @(x, a, b, c) max(min((x - a)/(b - a), (c - x)/(c - b)), 0);
t_gamma_func = @(x, a, k) k*(x - a).^2 ./ (1 + k*(x - a).^2) .* (x >= a);
g_gamma_func = @(x, a, m) ((x > a & x < m).*(x - a)/(m - a)) + (x >= m);
s_func = @(x, a, b) (x > a & x <= (a+b)/2).*2.*((x-a)/(b-a)).^2 + ...
         (x > (a+b)/2 & x < b).*(1 - 2*((x-b)/(b-a)).^2) + (x >= b);
gaus_func = @(x, m, k) exp(-k * (x - m).^2);
trapez_func = @(x, a, b, c, d) ((x > a & x < b).*(x - a)/(b - a)) + ...
         (x >= b & x <= c) + ((x > c & x < d).*(d - x)/(d - c));
pseudo_exp_func = @(x, m, k) 1 ./ (1 + k*(x - m).^2);
sing_func = @(x, a) double(x == a);
rho_gama_func = @(x, a, b, rho) (x <= a).* 0 +  (a<x<b).*((x-a) ./ (b-a)).^rho + ...
    (x>=b).*1;

%% Otros
cardinalidad = @(x) sum(x);
altura = @(x) max(x);
%% Funciones unarias
concentra_dilatac_func = @(x,p) x .^p;
intensi_d_contra_func = @(x,p) (x<=0.5).* ( 2^(p-1) * x .^p ) + ...
    (x > 0.5) .* (1 - 2^(p-1) * (1 - x.^p) );
fuzzi_difuminaci_func = @(x) (x<=0.5).* ( sqrt(x ./2) ) + ...
    (x > 0.5) .* ( 1 - sqrt( (1 - x)/2 ) );
normalizacion = @(x) x./ altura(x);
%% --------------------------------------------------------------
%% ---------------------------

%% PREGUNTA 1: a) Construir 3 conjunto borrosos A,B y C y graficar
k=0.8002;
%% Conjuntos borrosos de A: Perdida de paquetes en la transmisión de datos
bajoA = g_gamma_func(x, 1, 2)*0.32;
medioA= gaus_func(x, 8, 0.5)*0.98;
altoA= g_gamma_func(x,16,18)*0.98;
siluetaA = max([bajoA; medioA; altoA], [], 1);
figure;
plot(x, bajoA, '-c', 'LineWidth', 2);hold on;
plot(x, medioA, '-g', 'LineWidth', 2);
plot(x, altoA, '-m', 'LineWidth', 2);
plot(x, siluetaA, '--k', 'LineWidth', 1);
legend('Bajo', 'Medio', 'Alto','Silueta', 'Location', 'best');
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Conjuntos borrosos A Perdida de paquetes en la transmisión de datos');
grid on;
ylim([0 1.1]);
Ua = [3 4 5 6 10 11 12 13 14 15 16 19];

Au = evaluar_valor_vector(x,siluetaA,Ua);

% Gráfica: silueta
figure;
plot(x, siluetaA, '-k', 'LineWidth', 1); hold on;
plot(Ua, Au, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Silueta de los conjuntos borrosos A');
grid on;
ylim([0 1.1]);

% Forma 1 representar:
% Au = {0.32/3; 0.32/4; 0.32/5 ; 0.32/6; 0.32/10; 0.32/11; 0.32/12;
% 0.32/13; 0.32/14; 0.32/15;
%      0.32/16; 0.98/19}

%% Conjuntos borrosos de B : Ancho de banda transmisión de datos
k = 0.39161;
bajoB = gaus_func(x, 4, 0.2)*0.98*k;
medioB= gaus_func(x, 12, 0.01)*0.8*k;
altoB= g_gamma_func(x,18,18.9)*1;
siluetaB = max([bajoB; medioB; altoB], [], 1);
figure;
plot(x, bajoB, '-c', 'LineWidth', 2);hold on;
plot(x, medioB, '-g', 'LineWidth', 2);
plot(x, altoB, '-m', 'LineWidth', 2);
plot(x, siluetaB, '--k', 'LineWidth', 1);
legend('Bajo', 'Medio', 'Alto','Silueta', 'Location', 'best');
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Conjuntos borrosos B Ancho de banda transmisión de datos');
grid on;
ylim([0 1.1]);

Ub = [1 2 4 5 6 7 8 9 14 18 19 20];
Bu = evaluar_valor_vector(x,siluetaB,Ub);

% Gráfica: silueta
figure;
plot(x, siluetaB, '-k', 'LineWidth', 1); hold on;
plot(Ub, Bu, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Silueta de los conjuntos borrosos B');
grid on;
ylim([0 1.1]);

% Forma 2 representar:
% Bu = 0.0934/1+ 0.1729/2 + 0.3837/4 + 0.3129/5+ 0.2188/6+ 0.2443/7+
% 0.2673/8+ 0.2866/9+ 0.3011/14+ 0.2186/18 + 1/19 + 1/20}

%% Conjuntos borrosos de C: Latencia en la transmión de datos
l = 0.36425;
bajoC = gaus_func(x, 4, 0.5)*1*l;
medioC= g_gamma_func(x, 2, 7)*.8*l;
altoC= g_gamma_func(x,17.1,17.8);
siluetaC = max([bajoC; medioC; altoC], [], 1);
figure;
plot(x, bajoC, '-c', 'LineWidth', 2);hold on;
plot(x, medioC, '-g', 'LineWidth', 2);
plot(x, altoC, '-m', 'LineWidth', 2);
plot(x, siluetaC, '--k', 'LineWidth', 1);
legend('Bajo', 'Medio', 'Alto','Silueta', 'Location', 'best');
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Conjuntos borrosos C Latencia en la transmión de datos');
grid on;
ylim([0 1.1]);
Uc = [1 3 4 5 6 9 10 11 14 17 19 20];
Cu = evaluar_valor_vector(x,siluetaC,Uc);

% Gráfica: silueta
figure;
plot(x, siluetaC, '-k', 'LineWidth', 1); hold on;
plot(Uc, Cu, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Silueta de los conjuntos borrosos C');
grid on;
ylim([0 1.1]);

% Forma 1 representar:
% Cu = {0.004/1; 0.2222/3; 0.3642/4 ; 0.2187/5; 0.2338/6; 0.2914/9; 0.2914/10;
% 0.2914/11; 0.2914/14; 0.2914/17; 1/16; 1/19}

%% Cardinalidad aprimado 4.5
cardinalidad_A = cardinalidad(Au)
cardinalidad_B = cardinalidad(Bu)
cardinalidad_C = cardinalidad(Cu)

%% ---------------------------
%% ---------------------------

%% PREGUNTA 1: b) Hallar 2 concavidades y 2 convexidades lambda 0.8
% Forma 1 representar:
% Forma 2 representar:
% Bu = 0.0934/1+ 0.1729/2 + 0.3837/4 + 0.3129/5+ 0.2188/6+ 0.2443/7+
% 0.2673/8+ 0.2866/9+ 0.3011/14+ 0.2186/18 + 1/19 + 1/20}
%% CONVEXO evaluar_valor_escalar <= max (A(X1) , A(X2))
disp("Evaluar la convexidad")
lambda = 0.8;
B1 = 0.1729;
X1 = 2;
B2 = 0.3129;
X2 = 5;
valor = evaluar_valor_escalar(x, siluetaB, X1*lambda + (1-lambda)*X2)
% valor = 0.2600 Si es convexo <=max(0.1729,0.3129);
disp("Evaluar la convexidad")
lambda = 0.8;
B1 = 0.2673;
X1 = 8;
B2 = 0.3011;
X2 = 14;
valor = evaluar_valor_escalar(x, siluetaB, X1*lambda + (1-lambda)*X2)
% valor = 0.2899 Si es convexo <=max(0.2673,0.3011)
% ;
%% CONCAVO % evaluar_valor_escalar >= min (A(X1) , A(X2))
disp("Evaluar la concavidad")
lambda = 0.8;
B1 = 0.3129;
X1 = 5;
B2 = 0.2188;
X2 = 6;
valor = evaluar_valor_escalar(x, siluetaB, X1*lambda + (1-lambda)*X2)
% valor = 0.2862 Si es concavo >=min(0.2188,0.3129);
disp("Evaluar la concavidad")
lambda = 0.8;
B1 = 0.2186;
X1 = 18;
B2 = 1;
X2 = 19;
valor = evaluar_valor_escalar(x, siluetaB, X1*lambda + (1-lambda)*X2)
% valor = 0.2182 Si es concavo >=min(0.2186,1);

disp("------------------------------")
disp("------------------------------")
%% ---------------------------
%% PREGUNTA 2: a) Operaciones unarias
% Conjunto difuso Au
% Au = {0.32/3; 0.32/4; 0.32/5 ; 0.32/6; 0.32/10; 0.32/11; 0.32/12;
% 0.32/13; 0.32/14; 0.32/15;
%      0.32/16; 0.98/19}
disp("Operaciones unaria: Normalizacion")
nomalizacion = normalizacion(Au)
% nomalizacion = {0.3265 0.3265 0.3265 0.3265 0.3265 0.3265    
% 0.3265 0.3265 0.3265 0.3265 0.3265 1.0000}

disp("Operaciones unaria: concentracion, p=2")
concentracion = concentra_dilatac_func(Au,2)
% concentracion ={0.1024 0.1024 0.1024 0.1024 0.1024 0.1024 0.1024 0.1024
% 0.1024 0.1024 0.1024 0.960}

disp("Operaciones unaria: dilatacion, p=0.5")
dilatacion = concentra_dilatac_func(Au,0.5)
% dilatacion = {0.5657    0.5657    0.5657    0.5657    0.5657    0.5657
% 0.5657    0.5657    0.5657    0.5657    0.5657    0.9899}

disp("Operaciones unaria: intensificacion del contraste")
intensifica_d_contraste = intensi_d_contra_func(Au,2)
% intensifica_d_contraste ={0.2048    0.2048    0.2048    0.2048    0.2048
% 0.2048    0.2048    0.2048    0.2048    0.2048    0.2048    0.9208}

disp("Operaciones unaria: fuzzifiacion")
fuzzifiacion = fuzzi_difuminaci_func(Au)
% fuzzifiacion ={ 0.4000    0.4000    0.4000    0.4000    0.4000    0.4000
% 0.4000    0.4000    0.4000    0.4000    0.4000    0.9000}

disp("------------------------------")
disp("------------------------------")

%% ---------------------------
%% PREGUNTA 2: b) Inclusión difusa
disp("A en B")
inclus_A_en_B = inclusion_difusa(Au,Ua,Bu,Ub)
% 0.9717

disp("B en A")
inclus_B_en_A = inclusion_difusa(Bu,Ub,Au,Ua)
% 0.9814

% Interpretacion: B está un poco más incluido en A que A en B

disp("A en C")
inclus_A_en_C = inclusion_difusa(Au,Ua,Cu,Uc)
% 0.9175

disp("C en A")
inclus_C_en_A = inclusion_difusa(Cu,Uc,Au,Ua)
% 0.9857

% Interpretacion: C está un poco más incluido en A que A en C

disp("B en C")
inclus_B_en_C = inclusion_difusa(Bu,Ub,Cu,Uc)
% 0.9527

disp("C en B")
inclus_C_en_B = inclusion_difusa(Cu,Uc,Bu,Ub)
% 0.9956

% Interpretacion: C está más incluido en B como B en C


%% PREGUNTA 3: a) Hallar y graficar z-2x+2y
% Au = x
% Bu = y
% Cu = z;
% Empiezo con z-2x
[resultado_z_menos_2x, U_resultado_z_menos_2x] = operacion(Uc, Cu, Ua.*(-2), Au)

% (z-2x) + 2y
[resultado_z_menos_2x_mas_2y, U_resultado_z_menos_2x_mas_2y] = operacion(U_resultado_z_menos_2x, resultado_z_menos_2x, Ub.*(2), Bu)

figure;
plot(U_resultado_z_menos_2x_mas_2y, resultado_z_menos_2x_mas_2y, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Resultado z-2x+2y');
grid on;
ylim([-0.1 1.1]);
xlim([-30 56]);

disp("------------------------------")
disp("------------------------------")
%% ---------------------------

%% PREGUNTA 3: b) Hallar y graficar NOT (A OR B) AND C

% Empiezo hallando A OR B
[resultado_A_OR_B, U_A_OR_B] = zadeh_or(Ua, Au, Ub, Bu)

% negacion A OR B
[resultado_nega_A_OR_B, U_nega_A_OR_B] = negacion_difusa(U_A_OR_B, resultado_A_OR_B)

% NOT (A OR B) AND C
[resultado_final, U_final] = zadeh_and(U_nega_A_OR_B,resultado_nega_A_OR_B,Uc,Cu)
% Resultado final = {0.0040/1;  0.2222/3; 0.3642/4; 0.2187/5; 0.2338/6; 0.2914/9; 0.2914/10;
% 0.2914/11; 0.2914/14;  0/19;  0/20}
figure;
plot(U_final, resultado_final, 'rx', 'MarkerSize', 10, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Resultado NOT (A OR B) AND C');
grid on;
ylim([-0.1 1.1]);
xlim([0 20]);

%% --------------------------------------------------------------
%% --------------------------------------------------------------

%% Funciones
function valor = evaluar_valor_escalar(x_vect, silueta, x)
[~, i] = min(abs(x_vect - x));
valor = silueta(i);
valor = fix(valor * 10^4) / 10^4;
end

function valor = evaluar_valor_vector(x_vect, silueta, x)
[~, i] = min(abs(x_vect - x'), [], 2);
valor = silueta(i);
valor = fix(valor .* 10^4) ./ 10^4;
end

function S = inclusion_difusa(Xu, Ux, Yu, Uy)
cardX = sum(Xu);

[~, idxX, idxY] = intersect(Ux, Uy);

Xu_common = Xu(idxX);
Yu_common = Yu(idxY);

suma_diferencias = sum(max(0, Xu_common - Yu_common));

S = (cardX - suma_diferencias) / cardX;
end


function [resultado, U_comun] = zadeh_and(Ux, Xu, Uy, Yu)
U_comun = intersect(Ux, Uy);
[~, i_X] = ismember(U_comun, Ux);
[~, i_Y] = ismember(U_comun, Uy);
X_comun = Xu(i_X);
Y_comun = Yu(i_Y);
resultado = min(X_comun, Y_comun);
end

function [resultado, U_union] = zadeh_or(Ux, Xu, Uy, Yu)
U_union = unique([Ux(:); Uy(:)])'; 
X_union = zeros(size(U_union));
Y_union = zeros(size(U_union)); 
[e_X, i_X] = ismember(U_union, Ux);
X_union(e_X) = Xu(i_X(e_X));
[e_Y, i_Y] = ismember(U_union, Uy);
Y_union(e_Y) = Yu(i_Y(e_Y));
resultado = max(X_union, Y_union);
end

function [negacion, universo] = negacion_difusa(U, Xu)
negacion = 1 - Xu;
universo = U;
end

function [resultado, U_resultado] = implicacion_zadeh(Ux, Xu, Uy, Yu)
    
    [neg_X, U_neg] = negacion_difusa(Ux, Xu);
    
    [resultado, U_resultado] = zadeh_or(U_neg, neg_X, Uy, Yu);
end

function [resultado, U_resultado] = bicondicional_zadeh(Ux, Xu, Uy, Yu)
    [X_impl_Y, U_impl1] = implicacion_zadeh(Ux, Xu, Uy, Yu);
   
    [Y_impl_X, U_impl2] = implicacion_zadeh(Uy, Yu, Ux, Xu);
    
    [resultado, U_resultado] = zadeh_and(U_impl1, X_impl_Y, U_impl2, Y_impl_X);
end


function [Ux_soporte, Xu_soporte] = soporte_alpha(Ux, Xu, alpha)
indice = Xu >= alpha;
Ux_soporte = Ux(indice);
Xu_soporte = Xu(indice);
end

function [resultado_final, U_resultado_final] = operacion(Ux, Xu, Uy, Yu)
U_temp = [];
resultado_temp = [];
for i = 1:length(Ux)
for j = 1:length(Uy)
U_temp(end+1) = Ux(i) + Uy(j);
resultado_temp(end+1) = min(Xu(i), Yu(j));
end
end
[U_resultado_final, ~, idx] = unique(U_temp);
resultado_final = zeros(size(U_resultado_final));
for k = 1:length(U_resultado_final)
resultado_final(k) = max(resultado_temp(idx == k));
end
end