%% Viajero
clear all; close all; clc;
% Parámetros
tam_poblacion = 100; max_generaciones = 500;
prob_cruce = 0.8; prob_mutacion = 0.2; n_ciudades = 10;

disp("a) Evaluación de aptitud: distancia mínima");
disp("b) Selección: torneo");
disp("c) Cruce: OX");
disp("d) Mutación: intercambio");

% Matriz de distancias y coordenadas
rng(4);
distancias = randi([1 100], n_ciudades, n_ciudades);
distancias = triu(distancias) + triu(distancias, 1)';
coordenadas = rand(n_ciudades, 2) * 100;

% Población inicial
poblacion = zeros(tam_poblacion, n_ciudades);
for i = 1:tam_poblacion
    poblacion(i, :) = randperm(n_ciudades);
end
disp("Población inicial"); disp(poblacion);

% Algoritmos genéticos
historial = zeros(max_generaciones, 1);
mejores_distancias = zeros(max_generaciones, 1);
umbral_distancia = 200; % Umbral para considerar una solución buena
for g = 1:max_generaciones
    % Evaluar aptitud
    fitness = zeros(tam_poblacion, 1);
    soluciones_generacion = [];  % Para almacenar soluciones buenas
    for i = 1:tam_poblacion
        fitness(i) = fitness_tsp(poblacion(i, :), distancias);
        % Si la distancia es menor al umbral, guardar como solución
        d = (1/fitness(i)) -1;
        if d <= umbral_distancia
            soluciones_generacion = [soluciones_generacion; poblacion(i,:) d];
        end
    end
    [mejor_fit, idx] = max(fitness);
    historial(g) = mejor_fit;
    
    % Mostrar progreso
    if mod(g, 50) == 0
        fprintf('\nGeneración %d: Soluciones con distancia <= %.2f:\n', g, umbral_distancia);
        if ~isempty(soluciones_generacion)
            % Eliminar rutas duplicadas
            [~, idx_unicas] = unique(soluciones_generacion(:,1:n_ciudades), 'rows', 'stable');
            soluciones_generacion = soluciones_generacion(idx_unicas,:);
            
            % Ordenar por distancia
            [distancias_ordenadas, idx_orden] = sort(soluciones_generacion(:,end));
            soluciones_ordenadas = soluciones_generacion(idx_orden, 1:end-1);

            num_mostrar = min(5, size(soluciones_ordenadas, 1));
            for s = 1:num_mostrar
                fprintf('Ruta %d: %s | Distancia: %.2f\n', s, ...
                    mat2str(soluciones_ordenadas(s,:)), distancias_ordenadas(s));
            end
            if size(soluciones_ordenadas, 1) > 5
                fprintf('... (+%d más)\n', size(soluciones_ordenadas, 1)-5);
            end
        else
            fprintf('No se encontraron soluciones bajo el umbral en esta generación\n');
        end
    end
    
    % Selección por torneo
    padres = seleccion_torneo(poblacion, fitness);
    
    % Cruce OX
    descendencia = cruceOX(padres, prob_cruce, n_ciudades);
    
    % Mutación
    for i = 1:tam_poblacion
        if rand < prob_mutacion
            descendencia(i, :) = mutacion_intercambio(descendencia(i, :));
        end
    end
    
    % Elitismo y Reemplazo
    descendencia(1, :) = poblacion(idx, :);
    poblacion = descendencia;
end

%% Mostrar resultados
fprintf('\n=== RESULTADOS FINALES ===\n');
[mejor_fit, idx] = max(fitness);
mejor_ind = poblacion(idx, :);
fprintf('Mejor ruta: %s\nDistancia: %.2f\n', mat2str(mejor_ind), 1/mejor_fit - 1);

% Gráfica de convergencia
figure; plot(1:max_generaciones, historial, 'LineWidth', 2);
xlabel('Generación'); ylabel('fitness'); title('Convergencia'); grid on;

% Gráfica de ruta
mejor_ruta = [mejor_ind mejor_ind(1)];
figure; plot(coordenadas(mejor_ruta, 1), coordenadas(mejor_ruta, 2), '-or');
% Añadir etiquetas de ciudades
for i = 1:n_ciudades
    text(coordenadas(i,1), coordenadas(i,2), sprintf('Ciudad %d', i));
end
title('Ruta Final'); xlabel('X'); ylabel('Y'); axis equal; grid on;

% Función de aptitud
function fit = fitness_tsp(ind, dist)
    d = 0;
    for i = 1:length(ind)-1
        d = d + dist(ind(i), ind(i+1));
    end
    d = d + dist(ind(end), ind(1));
    fit = 1 / (1 + d);
end

% Selección por torneo
function padres = seleccion_torneo(pobl, fit)
    n = size(pobl, 1);
    k=2;
    padres = zeros(size(pobl));
    for i = 1:n
        idxs = randsample(n, k);
        [~, best] = max(fit(idxs));
        padres(i,:) = pobl(idxs(best), :);
    end
end

% Cruce OX
function hijos = cruceOX(pobl, prob, d)
    n = size(pobl, 1);
    hijos = pobl;
    for i = 1:2:n-1
        if rand < prob
            % Puntos de cruce
            p1 = randi([1 d-1]);
            p2 = randi([p1+1 d]);
            
            % Hijo 1
            hijo1 = zeros(1, d);
            hijo1(p1:p2) = pobl(i, p1:p2);
            j = mod(p2, d) + 1;
            for k = 1:d
                g = pobl(i+1, mod(p2 + k - 1, d) + 1);
                if ~ismember(g, hijo1)
                    hijo1(j) = g;
                    j = mod(j, d) + 1;
                end
            end
            
            % Hijo 2
            hijo2 = zeros(1, d);
            hijo2(p1:p2) = pobl(i+1, p1:p2);
            j = mod(p2, d) + 1;
            for k = 1:d
                g = pobl(i, mod(p2 + k - 1, d) + 1);
                if ~ismember(g, hijo2)
                    hijo2(j) = g;
                    j = mod(j, d) + 1;
                end
            end
            
            hijos(i, :) = hijo1;
            hijos(i+1, :) = hijo2;
        end
    end
end

% Mutación por intercambio
function ind = mutacion_intercambio(ind)
    if length(ind) > 1
        pos = randperm(length(ind), 2);
        temp = ind(pos(1));
        ind(pos(1)) = ind(pos(2));
        ind(pos(2)) = temp;
    end
end
