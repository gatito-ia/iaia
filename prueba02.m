clear all; close all; clc;

%% 1. Crear el sistema difuso Mamdani
fis = mamfis('Name', 'ControlInvernadero');
% sb = sugfis('Name', 'ConsumoTSK');

%%
% fismf | fismftype2 | psigmf | sigmf | gaussmf | gauss2mf | gbellmf |
% trimf | trapmf | linsmf | linzmf | pimf | zmf | dsigmf | smf

%% Funciones
% Mamdani salida "trimf" (triangular), "trapmf" (trapezoidal), "gaussmf" (gaussiana), 
% Sugeno linear/ constant
%% 2. Definir las entradas con mayor precisión
% Entrada 1: Temperatura (°C) - Rangos ajustados para mejor transición
fis = addInput(fis, [0 40], 'Name', 'Temperatura');
fis = addMF(fis, 'Temperatura', 'trapmf', [0 0 15 20], 'Name', 'Fria'); % Trapezoidal para mejor transición
fis = addMF(fis, 'Temperatura', 'gaussmf', [5 25], 'Name', 'Optima'); % Gaussiana para centro preciso
fis = addMF(fis, 'Temperatura', 'trapmf', [25 30 40 40], 'Name', 'Caliente');

% Entrada 2: Humedad (%) - Corregido error tipográfico en "Humedad"
fis = addInput(fis, [0 100], 'Name', 'Humedad');
fis = addMF(fis, 'Humedad', 'trapmf', [0 0 30 40], 'Name', 'Seco');
fis = addMF(fis, 'Humedad', 'gaussmf', [15 50], 'Name', 'Normal');
fis = addMF(fis, 'Humedad', 'trapmf', [60 70 100 100], 'Name', 'Humedo');

% Entrada 3: Intensidad de Luz (lux) - Rangos ajustados
fis = addInput(fis, [0 1000], 'Name', 'Luz');
fis = addMF(fis, 'Luz', 'trapmf', [0 0 200 300], 'Name', 'Baja');
fis = addMF(fis, 'Luz', 'gaussmf', [150 500], 'Name', 'Moderada');
fis = addMF(fis, 'Luz', 'trapmf', [700 800 1000 1000], 'Name', 'Alta');

%% 3. Definir las salidas con métodos de defuzzificación explícitos
% Configurar método de defuzzificación (antes de añadir MFs)
fis.DefuzzificationMethod = 'centroid';
% Métodos de Defuzzificación
% 'centroid' (Centroide, default) → Calcula el centro de gravedad del área difusa.
% 'bisector' (Bisector) → Divide el área en dos partes iguales.
% 'mom' (Mean of Maximum) → Promedio de los valores con máxima pertenencia.
% 'lom' (Largest of Maximum) → Mayor valor con máxima pertenencia.
% 'som' (Smallest of Maximum) → Menor valor con máxima pertenencia.

% Salida 1: Ventilacion (%)
fis = addOutput(fis, [0 100], 'Name', 'Ventilacion');
fis = addMF(fis, 'Ventilacion', 'trapmf', [0 0 20 30], 'Name', 'Minima');
fis = addMF(fis, 'Ventilacion', 'gaussmf', [15 50], 'Name', 'Media');
fis = addMF(fis, 'Ventilacion', 'trapmf', [70 80 100 100], 'Name', 'Maxima');

% Salida 2: Riego (%)
fis = addOutput(fis, [0 100], 'Name', 'Riego');
fis = addMF(fis, 'Riego', 'trapmf', [0 0 10 20], 'Name', 'Nulo');
fis = addMF(fis, 'Riego', 'gaussmf', [10 30], 'Name', 'Leve');
fis = addMF(fis, 'Riego', 'trapmf', [40 60 100 100], 'Name', 'Fuerte');

%% 4. Definir las reglas mejor estructuradas
% Formato: [in1 in2 in3 out1 out2 weight operator]
% listaReglas = [
%     1 1 1 1 1 1 1;   % Si Temp=Fria Y Humedo=Seco Y Luz=Baja → Vent=Minima, Riego=Nulo
%     2 2 2 2 2 1 1;    % Si Temp=Optima Y Humedo=Normal Y Luz=Moderada → Vent=Media, Riego=Leve
%     3 3 3 3 3 1 1;    % Si Temp=Caliente Y Humedo=Humedo Y Luz=Alta → Vent=Maxima, Riego=Fuerte
% ];
% fis = addRule(fis, listaReglas);

rule1 = "Temperatura==Fria & Humedad==Seco & Luz==Baja => Ventilacion=Minima, Riego=Nulo (1)";
rule2 = "Temperatura==Optima & Humedad==Normal & Luz==Moderada => Ventilacion=Media, Riego=Leve (1)";
rule3 = "Temperatura==Caliente & Humedad==Humedo & Luz==Alta => Ventilacion=Maxima, Riego=Fuerte (1)";
rule4 = "Temperatura==Fria & Humedad==Humedo => Ventilacion=Minima, Riego=Fuerte (1)";
rule5 = "Temperatura==Caliente & Luz==Alta => Ventilacion=Maxima, Riego=Leve (1)";

rules = [rule1 rule2 rule3 rule4 rule5];
fis = addRule(fis, rules);
% Symbol	Keyword
% ==	IS (in rule antecedent)
% ~=	IS NOT
% &	AND
% |	OR
% => THEN
% =	IS (in rule consequent)
% "A==a & B~=b => X=x, Y~=y (1)"

%% 5. Evaluar el sistema
% Caso de prueba: Temp=28°C, Humedad=65%, Luz=600 lux
input_values = [28 65 600];
outputs = evalfis(fis, input_values); % Capturar ambas salidas correctamente

fprintf('\n=== RESULTADOS DEL SISTEMA ===\n');
fprintf('Para Temp=%.1f°C, Humedad=%.1f%%, Luz=%.1f lux:\n', input_values);
fprintf('-> Ventilación: %.1f%%\n', outputs(1));
fprintf('-> Riego: %.1f%%\n\n', outputs(2));

%% 6. Visualización mejorada

% Graficar funciones
figure;
subplot(3,1,1); 
plotmf(fis, 'input', 1); 
title('Temperatura (°C)');

subplot(3,1,2); 
plotmf(fis, 'input', 2); 
title('Humedad (%)');

subplot(3,1,3); 
plotmf(fis, 'input', 3); 
title('Intensidad de Luz (lux)');

figure;
subplot(2,1,1);
plotmf(fis, 'output', 1); 
title('Ventilación');

subplot(2,1,2);
plotmf(fis, 'output', 2); 
title('Riego');

%% Diagrama del sistema
figure;
plotfis(fis);

%% Ver reglas
disp('Reglas del sistema difuso');
showrule(fis)


% %% Evaluar vector
% mand_value = 20:1:25;
% mand_salida = arrayfun(@(x) evalfis(fis,x), mand_value);
% U_temp =[];
% U_temp(end+1) = Ux
% x = linspace(0, 20, 500);
% fix(valor * 10^4) / 10^4;
% r = randi([a b], m, n)
% r = a + (b - a) * rand(m, n)

% figure;
% plot(mand_value,mand_salida, 'b','LineWidth',2);
% plot(U_final, resultado_final, 'rx', 'MarkerSize', 10, 'LineWidth', 1);
% xlabel('Eje x');
% ylabel('Eje y');
% title('Sistema Difuso Mandani');
% grid on;
% ylim([-0.1 1.1]);
% xlim([0 20]);
