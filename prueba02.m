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
%% Entradas
% Entrada 1: Temperatura (°C)
fis = addInput(fis, [0 40], 'Name', 'Temperatura');
fis = addMF(fis, 'Temperatura', 'trapmf', [0 0 15 20], 'Name', 'Fria');
fis = addMF(fis, 'Temperatura', 'gaussmf', [5 25], 'Name', 'Optima');
fis = addMF(fis, 'Temperatura', 'trapmf', [25 30 40 40], 'Name', 'Caliente');

% Entrada 2: Humedad (%)
fis = addInput(fis, [0 100], 'Name', 'Humedad');
fis = addMF(fis, 'Humedad', 'trapmf', [0 0 30 40], 'Name', 'Seco');
fis = addMF(fis, 'Humedad', 'gaussmf', [15 50], 'Name', 'Normal');
fis = addMF(fis, 'Humedad', 'trapmf', [60 70 100 100], 'Name', 'Humedo');

% Entrada 3: Intensidad de Luz (lux)
fis = addInput(fis, [0 1000], 'Name', 'Luz');
fis = addMF(fis, 'Luz', 'trapmf', [0 0 200 300], 'Name', 'Baja');
fis = addMF(fis, 'Luz', 'gaussmf', [150 500], 'Name', 'Moderada');
fis = addMF(fis, 'Luz', 'trapmf', [700 800 1000 1000], 'Name', 'Alta');

%% Definir las salidas con metodos de defuzzificación
% Configurar método de defuzzificación
fis.DefuzzificationMethod = 'centroid';
% Métodos de Defuzzificacion
% 'centroid' (Centroide, default) - Calcula el centro de gravedad del area difusa.
% 'bisector' (Bisector) - Divide el area en dos partes iguales.
% 'mom' (Mean of Maximum) - Promedio de los valores con máxima pertenencia.
% 'lom' (Largest of Maximum) - Mayor valor con maxima pertenencia.
% 'som' (Smallest of Maximum) - Menor valor con maxima pertenencia.

%% SAlidas
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

%% Definir las reglas
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

%% Evaluar el sistema
% Caso de prueba: Temp=28C, Humedad=65%, Luz=600 lux
fprintf('\n=== RESULTADOS DEL SISTEMA ===\n');
input_values = [28 65 600;
                40 70 800];
for i=1:size(input_values,1)
    outputs = evalfis(fis, input_values(i,:));
    fprintf('Para Temp=%.1fC, Humedad=%.1f%%, Luz=%.1f lux:\n', input_values(i,:));
    fprintf('-> Ventilación: %.1f%%\n', outputs(1));
    fprintf('-> Riego: %.1f%%\n\n', outputs(2));
end

%% Visualizacion mejorada
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

% %% Visualización de la superficie de respuesta
% figure;
% gensurf(fis);
% title('Superficie de Respuesta: Ventilación vs Temperatura y Humedad');
% xlabel('Temperatura (°C)');
% ylabel('Humedad (%)');
% zlabel('Ventilación (%)');
% 
% %% También puedes visualizar secciones específicas
% % Por ejemplo, para rangos específicos con más puntos:
% figure;
% [x,y,z] = gensurf(fis, 'Temperatura', 'Humedad', 'Ventilacion', [0:1:40], [0:2:100]);
% surf(x,y,z);
% title('Superficie de Respuesta Detallada');
% xlabel('Temperatura (°C)');
% ylabel('Humedad (%)');
% zlabel('Ventilación (%)');
% colorbar;
% 
% %% Visualización de las funciones de membresía
% figure;
% subplot(2,1,1); 
% plotmf(fis, 'input', 1); 
% title('Temperatura (°C)');
% subplot(2,1,2); 
% plotmf(fis, 'input', 2); 
% title('Humedad (%)');
% figure;
% plotmf(fis, 'output', 1); 
% title('Ventilación');