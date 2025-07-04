%% NSGA-II para Problemas Multiobjetivo
clear all; close all; clc;

% ============== DEFINICIÓN DEL PROBLEMA ==============
% Ejemplo: Problema de la mochila bi-objetivo
problema.nombre = 'Mochila Bi-Objetivo';
problema.tam_cromosoma = 10; % Número de items
problema.limite_inf = zeros(1, problema.tam_cromosoma);
problema.limite_sup = ones(1, problema.tam_cromosoma);

% Datos del problema (pesos y dos tipos de valores)
problema.pesos = [10, 20, 30, 15, 25, 5, 35, 12, 22, 18];
problema.valor1 = [60, 100, 120, 70, 90, 30, 150, 50, 80, 110]; % Objetivo 1: Maximizar
problema.valor2 = [5, 3, 8, 6, 7, 2, 9, 4, 6, 7]; % Objetivo 2: Maximizar (ej: importancia)
problema.capacidad_maxima = 100;
problema.objeto= 2;

% Función de evaluación multiobjetivo
problema.funcion_fitness = @(x) [sum(problema.valor1(x == 1)), sum(problema.valor2(x == 1))];

% ============== PARÁMETROS NSGA-II ==============
parametros.tam_poblacion = 100; % Tamaño de la población
parametros.max_generaciones = 200; % Número de generaciones
parametros.prob_cruce = 0.8; % Probabilidad de cruce
parametros.prob_mutacion = 0.02; % Probabilidad de mutación

% ============== EJECUCIÓN ==============
[frente_pareto, poblacion_final] = NSGA_II(problema, parametros);

% ============== RESULTADOS ==============
disp('=== Soluciones en el Frente de Pareto ===');
for i = 1:size(frente_pareto, 1)
    disp(['Solución ', num2str(i), ':']);
    disp(['  Objetivo 1 (Valor1): ', num2str(frente_pareto(i, 1))]);
    disp(['  Objetivo 2 (Valor2): ', num2str(frente_pareto(i, 2))]);
    disp(['  Items seleccionados: ', mat2str(find(poblacion_final(i, :) == 1))]);
end

%% ================== FUNCIÓN NSGA-II ==================
function [frente_pareto, poblacion_final] = NSGA_II(problema, parametros)
    % Inicialización
    n = parametros.tam_poblacion;
    poblacion = randi([0, 1], n, problema.tam_cromosoma); % Población binaria
    
    for gen = 1:parametros.max_generaciones
        % Evaluación de la población
        fitness = zeros(n, problema.objeto);
        for i = 1:n
            fitness(i, :) = problema.funcion_fitness(poblacion(i, :));
        end
        
        % Clasificación por no-dominancia y crowding distance
        [frentes, crowding] = non_dominated_sort(fitness);
        
        % Selección de padres (torneo binario)
        padres = seleccion_torneo(poblacion, frentes, crowding);
        
        % Cruce y mutación
        descendencia = cruce_mutacion(padres, problema, parametros);
        
        % Nueva población (elitismo)
        poblacion = [poblacion; descendencia];
        fitness = [fitness; zeros(size(descendencia, 1), 2)];
        for i = n+1:size(poblacion, 1)
            fitness(i, :) = problema.funcion_fitness(poblacion(i, :));
        end
        
        % Selección de la siguiente generación
        [frentes, crowding] = non_dominated_sort(fitness);
        frentes = frentes(:); % Convertir a columna
        crowding = crowding(:); % Convertir a columna
        [~, idx] = sortrows([frentes, -crowding]);
        poblacion = poblacion(idx(1:n), :);
        fitness = fitness(idx(1:n), :); % Actualizar fitness para la siguiente generación
    end
    
    % Evaluación final de la población
    fitness_final = zeros(n, 2);
    for i = 1:n
        fitness_final(i, :) = problema.funcion_fitness(poblacion(i, :));
    end
    
    % Clasificación final para obtener el frente de Pareto
    [frentes_final, ~] = non_dominated_sort(fitness_final);
    frente_pareto = fitness_final(frentes_final == 1, :);
    poblacion_final = poblacion(frentes_final == 1, :);
end

%% ================== FUNCIONES AUXILIARES ==================
% 1. Clasificación por no-dominancia
function [frentes, crowding] = non_dominated_sort(fitness)
    n = size(fitness, 1);
    frentes = zeros(n, 1);
    crowding = zeros(n, 1);
    
    for i = 1:n
        dominado_por = 0;
        for j = 1:n
            if all(fitness(i, :) <= fitness(j, :)) && any(fitness(i, :) < fitness(j, :))
                dominado_por = dominado_por + 1;
            end
        end
        frentes(i) = dominado_por + 1;
    end
    
    % Crowding distance
    for f = 1:max(frentes)
        idx = find(frentes == f);
        if ~isempty(idx)
            crowding(idx) = linspace(1, 0, length(idx));
        end
    end
end

% 2. Selección por torneo
function padres = seleccion_torneo(poblacion, frentes, crowding)
    n = size(poblacion, 1);
    padres = zeros(n, size(poblacion, 2));
    for i = 1:n
        % Torneo entre 2 individuos aleatorios
        idx = randperm(n, 2);
        if frentes(idx(1)) < frentes(idx(2)) || ...
           (frentes(idx(1)) == frentes(idx(2)) && crowding(idx(1)) > crowding(idx(2)))
            padres(i, :) = poblacion(idx(1), :);
        else
            padres(i, :) = poblacion(idx(2), :);
        end
    end
end

% 3. Cruce y mutación
function descendencia = cruce_mutacion(padres, problema, parametros)
    n = size(padres, 1);
    descendencia = padres;
    for i = 1:2:n-1
        % Cruce en un punto
        if rand() < parametros.prob_cruce
            punto = randi([1, problema.tam_cromosoma-1]);
            descendencia(i, :) = [padres(i, 1:punto), padres(i+1, punto+1:end)];
            descendencia(i+1, :) = [padres(i+1, 1:punto), padres(i, punto+1:end)];
        end
        % Mutación bit-flip
        for j = i:i+1
            for k = 1:problema.tam_cromosoma
                if rand() < parametros.prob_mutacion
                    descendencia(j, k) = 1 - descendencia(j, k);
                end
            end
        end
    end
end