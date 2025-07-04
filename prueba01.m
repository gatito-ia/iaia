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
% Perdida de paquetes en la transmision de datos
% Rendimiento de ancho de banda en la transmision de datos
% Velocidad de conexion
% Nivel de congestion de red
% Nivel de riesgo de amenaza
% Latencia en la transmision de datos
%% Concavo convexo
% lambda = 0.7;
% evaluar_valor_escalar(x, siluetaA, X1*lambda + (1-lambda)*X2);
% CONVEXO
% evaluar_valor_escalar <= max (A(X1) , A(X2))
% CONCAVO
% evaluar_valor_escalar >= min (A(X1) , A(X2))
%% Funciones transformar CB
% Operacion
% Metodo AND min
% Metodo OR max

%% Universo
x = linspace(0, 20, 500);

%% funcion de pertenencia
tria_f = @(x, a, b, c) max(min((x - a)/(b - a), (c - x)/(c - b)), 0);
t_gamma_f = @(x, a, k) k*(x - a).^2 ./ (1 + k*(x - a).^2) .* (x >= a);
g_gamma_f = @(x, a, m) ((x > a & x < m).*(x - a)/(m - a)) + (x >= m);
s_f = @(x, a, b) (x > a & x <= (a+b)/2).*2.*((x-a)/(b-a)).^2 + ...
         (x > (a+b)/2 & x < b).*(1 - 2*((x-b)/(b-a)).^2) + (x >= b);
gaus_f = @(x, m, k) exp(-k * (x - m).^2);
trapez_f = @(x, a, b, c, d) ((x > a & x < b).*(x - a)/(b - a)) + ...
         (x >= b & x <= c) + ((x > c & x < d).*(d - x)/(d - c));
pseudo_exp_f = @(x, m, k) 1 ./ (1 + k*(x - m).^2);
sing_f = @(x, a) double(x == a);
rho_gama_f = @(x, a, b, rho) (x <= a).* 0 +  (a<x<b).*((x-a) ./ (b-a)).^rho + ...
    (x>=b).*1;

%% Otros
cardin = @(x) sum(x);
alt = @(x) max(x);
%% Funciones unarias
concent_dilatac_f = @(x,p) x .^p;
intensi_contra_f = @(x,p) (x<=0.5).* ( 2^(p-1) * x .^p ) + ...
    (x > 0.5) .* (1 - 2^(p-1) * (1 - x.^p) );
fuzzi_difumin_f = @(x) (x<=0.5).* ( sqrt(x ./2) ) + ...
    (x > 0.5) .* ( 1 - sqrt( (1 - x)/2 ) );
normal = @(x) x./ alt(x);
%% --------------------------------------------------------------
%% ---------------------------

%% PREGUNTA 1: a) Construir 3 conjunto borrosos A,B y C y graficar
k=0.8002;
%% Conjuntos borrosos de A: Perdida de paquetes en la transmision de datos
bA = g_gamma_f(x, 1, 2)*0.32;
mA= gaus_f(x, 8, 0.5)*0.98;
aA= g_gamma_f(x,16,18)*0.98;
silueta_A = max([bA; mA; aA], [], 1);
figure;
plot(x, bA, '-c', 'LineWidth', 2);hold on;
plot(x, mA, '-g', 'LineWidth', 2);
plot(x, aA, '-m', 'LineWidth', 2);
plot(x, silueta_A, '--k', 'LineWidth', 1);
legend('Bajo', 'Medio', 'Alto','Silueta', 'Location', 'best');
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Conjuntos borrosos A Perdida de paquetes en la transmision de datos');
grid on;
ylim([0 1.1]);
Ua = [3 4 5 6 10 11 12 13 14 15 16 19];
Au = evaluarVector(x,silueta_A,Ua);

% Grafica: silueta
figure;
plot(x, silueta_A, '-k', 'LineWidth', 1); hold on;
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

%% Conjuntos borrosos de B : Ancho de banda transmision de datos
k = 0.39161;
bB = gaus_f(x, 4, 0.2)*0.98*k;
mB= gaus_f(x, 12, 0.01)*0.8*k;
aB= g_gamma_f(x,18,18.9)*1;
silueta_B = max([bB; mB; aB], [], 1);
figure;
plot(x, bB, '-c', 'LineWidth', 2);hold on;
plot(x, mB, '-g', 'LineWidth', 2);
plot(x, aB, '-m', 'LineWidth', 2);
plot(x, silueta_B, '--k', 'LineWidth', 1);
legend('Bajo', 'Medio', 'Alto','Silueta', 'Location', 'best');
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Conjuntos borrosos B Ancho de banda transmision de datos');
grid on;
ylim([0 1.1]);

Ub = [1 2 4 5 6 7 8 9 14 18 19 20];
Bu = evaluarVector(x,silueta_B,Ub);

% Grafica: silueta
figure;
plot(x, silueta_B, '-k', 'LineWidth', 1); hold on;
plot(Ub, Bu, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Silueta de los conjuntos borrosos B');
grid on;
ylim([0 1.1]);

% Forma 2 representar:
% Bu = 0.0934/1+ 0.1729/2 + 0.3837/4 + 0.3129/5+ 0.2188/6+ 0.2443/7+
% 0.2673/8+ 0.2866/9+ 0.3011/14+ 0.2186/18 + 1/19 + 1/20}

%% Conjuntos borrosos de C: Latencia en la transmion de datos
l = 0.36425;
bC = gaus_f(x, 4, 0.5)*1*l;
mC= g_gamma_f(x, 2, 7)*.8*l;
aC= g_gamma_f(x,17.1,17.8);
silueta_C = max([bC; mC; aC], [], 1);
figure;
plot(x, bC, '-c', 'LineWidth', 2);hold on;
plot(x, mC, '-g', 'LineWidth', 2);
plot(x, aC, '-m', 'LineWidth', 2);
plot(x, silueta_C, '--k', 'LineWidth', 1);
legend('Bajo', 'Medio', 'Alto','Silueta', 'Location', 'best');
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Conjuntos borrosos C Latencia en la transmion de datos');
grid on;
ylim([0 1.1]);
Uc = [1 3 4 5 6 9 10 11 14 17 19 20];
Cu = evaluarVector(x,silueta_C,Uc);

% Gráfica: silueta
figure;
plot(x, silueta_C, '-k', 'LineWidth', 1); hold on;
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
cardin_A = cardin(Au)
cardin_B = cardin(Bu)
cardin_C = cardin(Cu)

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
valor = evaluarEscalar(x, silueta_B, X1*lambda + (1-lambda)*X2)
% valor = 0.2600 Si es convexo <=max(0.1729,0.3129);
disp("Evaluar la convexidad")
lambda = 0.8;
B1 = 0.2673;
X1 = 8;
B2 = 0.3011;
X2 = 14;
valor = evaluarEscalar(x, silueta_B, X1*lambda + (1-lambda)*X2)
% valor = 0.2899 Si es convexo <=max(0.2673,0.3011)

%% CONCAVO % evaluar_valor_escalar >= min (A(X1) , A(X2))
disp("Evaluar la concavidad")
lambda = 0.8;
B1 = 0.3129;
X1 = 5;
B2 = 0.2188;
X2 = 6;
valor = evaluarEscalar(x, silueta_B, X1*lambda + (1-lambda)*X2)
% valor = 0.2862 Si es concavo >=min(0.2188,0.3129);
disp("Evaluar la concavidad")
lambda = 0.8;
B1 = 0.2186;
X1 = 18;
B2 = 1;
X2 = 19;
valor = evaluarEscalar(x, silueta_B, X1*lambda + (1-lambda)*X2)
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
nomaliz = normal(Au)
% nomalizacion = {0.3265 0.3265 0.3265 0.3265 0.3265 0.3265    
% 0.3265 0.3265 0.3265 0.3265 0.3265 1.0000}

disp("Operaciones unaria: concentracion, p=2")
concentr = concent_dilatac_f(Au,2)
% concentracion ={0.1024 0.1024 0.1024 0.1024 0.1024 0.1024 0.1024 0.1024
% 0.1024 0.1024 0.1024 0.960}

disp("Operaciones unaria: dilatacion, p=0.5")
dilatac = concent_dilatac_f(Au,0.5)
% dilatacion = {0.5657    0.5657    0.5657    0.5657    0.5657    0.5657
% 0.5657    0.5657    0.5657    0.5657    0.5657    0.9899}

disp("Operaciones unaria: intensificacion del contraste")
intensif_d_contra = intensi_contra_f(Au,2)
% intensifica_d_contraste ={0.2048    0.2048    0.2048    0.2048    0.2048
% 0.2048    0.2048    0.2048    0.2048    0.2048    0.2048    0.9208}

disp("Operaciones unaria: fuzzifiacion")
fuzzifia = fuzzi_difumin_f(Au)
% fuzzifiacion ={ 0.4000    0.4000    0.4000    0.4000    0.4000    0.4000
% 0.4000    0.4000    0.4000    0.4000    0.4000    0.9000}

disp("------------------------------")
disp("------------------------------")

%% ---------------------------
%% PREGUNTA 2: b) Inclusion difusa
disp("A en B")
inclus_A_en_B = inclus_difusa(Au,Ua,Bu,Ub)
% 0.9717

disp("B en A")
inclus_B_en_A = inclus_difusa(Bu,Ub,Au,Ua)
% 0.9814

% Interpretacion: B esta un poco mas incluido en A que A en B

disp("A en C")
inclus_A_en_C = inclus_difusa(Au,Ua,Cu,Uc)
% 0.9175

disp("C en A")
inclus_C_en_A = inclus_difusa(Cu,Uc,Au,Ua)
% 0.9857

% Interpretacion: C esta un poco mas incluido en A que A en C

disp("B en C")
inclus_B_en_C = inclus_difusa(Bu,Ub,Cu,Uc)
% 0.9527

disp("C en B")
inclus_C_en_B = inclus_difusa(Cu,Uc,Bu,Ub)
% 0.9956

% Interpretacion: C esta mas incluido en B como B en C


%% PREGUNTA 3: a) Hallar y graficar z-2x+2y
% Au = x
% Bu = y
% Cu = z;
% Empiezo con z-2x
[resul_z_menos_2x, U_resul_z_menos_2x] = Operacion(Uc, Cu, Ua.*(-2), Au)

% (z-2x) + 2y
[resul_z_menos_2x_mas_2y, U_resul_z_menos_2x_mas_2y] = Operacion(U_resul_z_menos_2x, resul_z_menos_2x, Ub.*(2), Bu)

figure;
plot(U_resul_z_menos_2x_mas_2y, resul_z_menos_2x_mas_2y, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
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
[resul_A_OR_B, U_A_OR_B] = Or(Ua, Au, Ub, Bu)

% negacion A OR B
[resul_nega_A_OR_B, U_nega_A_OR_B] = Negacion(U_A_OR_B, resul_A_OR_B)

% NOT (A OR B) AND C
[resul_final, U_final] = And(U_nega_A_OR_B,resul_nega_A_OR_B,Uc,Cu)
% Resultado final = {0.0040/1;  0.2222/3; 0.3642/4; 0.2187/5; 0.2338/6; 0.2914/9; 0.2914/10;
% 0.2914/11; 0.2914/14;  0/19;  0/20}
figure;
plot(U_final, resul_final, 'rx', 'MarkerSize', 10, 'LineWidth', 1);
xlabel('Variable linguistica');
ylabel('Grado de pertenencia');
title('Resultado NOT (A OR B) AND C');
grid on;
ylim([-0.1 1.1]);
xlim([0 20]);

%% --------------------------------------------------------------
%% --------------------------------------------------------------

%% Funciones
function n = evaluarEscalar(x_vect, silu, x)
[~, i] = min(abs(x_vect - x));
n = silu(i);
n = fix(n * 10^4) / 10^4;
end

function n = evaluarVector(x_vect, silu, x)
[~, i] = min(abs(x_vect - x'), [], 2);
n = silu(i);
n = fix(n .* 10^3) ./ 10^3;
end

function S = inclus_difusa(Xu, Ux, Yu, Uy)
card_X = sum(Xu);
[~, iX, iY] = intersect(Ux, Uy);
XuComun = Xu(iX);
YuComun = Yu(iY);
s = sum(max(0, XuComun - YuComun));
S = (card_X - s) / card_X;
end


function [r, UComun] = And(Ux, Xu, Uy, Yu)
UComun = intersect(Ux, Uy);
[~, iX] = ismember(UComun, Ux);
[~, iY] = ismember(UComun, Uy);
XComun = Xu(iX);
YComun = Yu(iY);
r = min(XComun, YComun);
end

function [r, UUnion] = Or(Ux, Xu, Uy, Yu)
UUnion = unique([Ux(:); Uy(:)])'; 
XUnion = zeros(size(UUnion));
YUnion = zeros(size(UUnion)); 
[x, iX] = ismember(UUnion, Ux);
XUnion(x) = Xu(iX(x));
[y, iY] = ismember(UUnion, Uy);
YUnion(y) = Yu(iY(y));
r = max(XUnion, YUnion);
end

function [n, UU] = Negacion(U, Xu)
n = 1 - Xu;
UU = U;
end

function [r, U_r] = Implicacion(Ux, Xu, Uy, Yu)
    [neg_X, U_neg] = Negacion(Ux, Xu);
    [r, U_r] = Or(U_neg, neg_X, Uy, Yu);
end

function [r, U_r] = Bicondicional(Ux, Xu, Uy, Yu)
    [X_impl_Y, U_impl_X_y] = Implicacion(Ux, Xu, Uy, Yu);
    [Y_impl_X, U_impl_Y_X] = Implicacion(Uy, Yu, Ux, Xu);
    [r, U_r] = And(U_impl_X_y, X_impl_Y, U_impl_Y_X, Y_impl_X);
end


function [Ux_sopo, Xu_sopo] = Soporte(Ux, Xu, tetha)
i = Xu >= tetha;
Ux_sopo = Ux(i);
Xu_sopo = Xu(i);
end

function [r_final, U_final] = Operacion(Ux, Xu, Uy, Yu)
U_temp = [];
r_temp = [];
for i = 1:length(Ux)
for j = 1:length(Uy)
U_temp(end+1) = Ux(i) + Uy(j);
r_temp(end+1) = min(Xu(i), Yu(j));
end
end
[U_final, ~, idx] = unique(U_temp);
r_final = zeros(size(U_final));
for k = 1:length(U_final)
r_final(k) = max(r_temp(idx == k));
end
end