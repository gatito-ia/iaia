%% Selección de Características para Detección de Intrusos
clear; clc; close all;

%% Definición del problema
num_caracteristicas = 20;  % Número total de características disponibles
max_caracteristicas = 10;  % Máximo número de características a seleccionar
num_muestras = 1000;       % Número de muestras en el dataset (simulado)

% Simulamos un dataset donde algunas características son más relevantes
% Las características 3, 7, 12 y 15 son las más importantes (artificial)
datos = rand(num_muestras, num_caracteristicas);
etiquetas = double(rand(num_muestras,1) > 0.5); % 50% benigno, 50% intruso

% Añadimos patrones artificiales para las características importantes
datos(etiquetas==1, 3) = datos(etiquetas==1, 3) + 0.5;
datos(etiquetas==1, 7) = datos(etiquetas==1, 7) - 0.3;
datos(etiquetas==1, 12) = datos(etiquetas==1, 12) * 1.8;
datos(etiquetas==1, 15) = datos(etiquetas==1, 15) + randn(sum(etiquetas==1),1)*0.2;

%% Parámetros del Algoritmo Genético
pob_size = 50;       % Tamaño de la población
max_gens = 100;      % Número máximo de generaciones
pc = 0.8;            % Probabilidad de cruce
pm = 0.05;           % Probabilidad de mutación

disp("a) Aptitud: Maximizar precisión + minimizar características");
disp("b) Selección: Torneo");
disp("c) Cruce: Uniforme");
disp("d) Mutación: Bit flip");

%% Inicialización de población
% Cada individuo es un vector binario donde 1=selección de característica
poblacion = rand(pob_size, num_caracteristicas) > 0.7;
disp("Población inicial generada");

%% Evolución
mejor_aptitud = zeros(max_gens, 1);
mejor_num_caract = zeros(max_gens, 1);

for gen = 1:max_gens
    % Evaluar población
    aptitudes = zeros(pob_size, 1);
    num_caract = zeros(pob_size, 1); % Número de características seleccionadas
    
    for i = 1:pob_size
        [aptitudes(i), num_caract(i)] = evaluar_individuo(poblacion(i,:), datos, etiquetas);
    end
    
    % Guardar el mejor
    [mejor_aptitud(gen), idx] = max(aptitudes);
    mejor_individuo = poblacion(idx,:);
    mejor_num_caract(gen) = num_caract(idx);
    
    % Mostrar progreso cada 20 generaciones
    if mod(gen, 20) == 0
        fprintf('Gen %d: Mejor aptitud=%.2f, Caract. seleccionadas=%d\n', ...
                gen, mejor_aptitud(gen), mejor_num_caract(gen));
    end
    
    % Selección (torneo)
    padres = false(size(poblacion));
    for i = 1:pob_size
        k = 3; % Tamaño del torneo
        competidores = randperm(pob_size, k);
        [~, idx_mejor] = max(aptitudes(competidores));
        padres(i,:) = poblacion(competidores(idx_mejor),:);
    end
    
    % Cruce (uniforme)
    hijos = false(size(padres));
    for i = 1:2:pob_size-1
        if rand < pc
            % Máscara de cruce aleatoria
            mascara = rand(1, num_caracteristicas) > 0.5;
            hijos(i,:) = padres(i,:).*mascara + padres(i+1,:).*(~mascara);
            hijos(i+1,:) = padres(i+1,:).*mascara + padres(i,:).*(~mascara);
        else
            hijos(i,:) = padres(i,:);
            hijos(i+1,:) = padres(i+1,:);
        end
    end
    
    % Mutación (bit flip)
    for i = 1:pob_size
        for j = 1:num_caracteristicas
            if rand < pm
                hijos(i,j) = ~hijos(i,j);
            end
        end
        
        % Asegurar que al menos una característica esté seleccionada
        if sum(hijos(i,:)) == 0
            hijos(i, randi(num_caracteristicas)) = true;
        end
        
        % Asegurar que no se exceda el máximo de características
        if sum(hijos(i,:)) > max_caracteristicas
            % Desactivar características aleatorias hasta cumplir el límite
            idx_activas = find(hijos(i,:));
            idx_desactivar = randperm(length(idx_activas), sum(hijos(i,:))-max_caracteristicas);
            hijos(i, idx_activas(idx_desactivar)) = false;
        end
    end
    
    % Elitismo (mantener el mejor individuo)
    hijos(1,:) = mejor_individuo;
    poblacion = hijos;
end

%% Resultados finales
[mejor_apt, idx] = max(aptitudes);
mejor_seleccion = poblacion(idx,:);
caract_seleccionadas = find(mejor_seleccion);

fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor aptitud encontrada: %.2f\n', mejor_apt);
fprintf('Número de características seleccionadas: %d de %d\n', sum(mejor_seleccion), num_caracteristicas);
fprintf('Características seleccionadas:\n');
disp(caract_seleccionadas);

% Verificar si se seleccionaron las características importantes (3,7,12,15)
fprintf('\nCaracterísticas importantes encontradas: ');
disp(intersect(caract_seleccionadas, [3 7 12 15]));

%% Visualización de convergencia
figure;
subplot(2,1,1);
plot(1:max_gens, mejor_aptitud, 'b-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Aptitud');
title('Evolución de la aptitud');
grid on;

subplot(2,1,2);
plot(1:max_gens, mejor_num_caract, 'r-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Número de características');
title('Evolución del número de características seleccionadas');
grid on;

%% Función de evaluación (fitness)
function [aptitud, num_caract] = evaluar_individuo(individuo, datos, etiquetas)
    % Obtener características seleccionadas
    idx_caract = find(individuo);
    num_caract = length(idx_caract);
    
    if num_caract == 0
        aptitud = 0;
        return;
    end
    
    % Extraer solo las características seleccionadas
    datos_subset = datos(:, idx_caract);
    
    % Métrica de evaluación simple (en la práctica usaríamos validación cruzada)
    % Usamos una regresión logística simple para evaluación
    mdl = fitglm(datos_subset, etiquetas, 'Distribution', 'binomial');
    y_pred = predict(mdl, datos_subset) > 0.5;
    precision = sum(y_pred == etiquetas) / length(etiquetas);
    
    % Penalización por usar muchas características
    penalizacion = 0.01 * num_caract;
    
    % Aptitud balanceada (mayor precisión con menos características)
    aptitud = precision - penalizacion;
end