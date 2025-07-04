%% Viajero
clear all; close all; clc;
% Parametros
tam_poblacion = 100; max_generaciones = 500;
prob_cruce = 0.8; prob_mutacion = 0.2; n_ciudades = 10;

disp("a) Evaluacion de aptitud: distancia minima");
disp("b) Seleccion: torneo");
disp("c) Cruce: OX");
disp("d) Mutacion: intercambio");

% Matriz de distancias y coordenadas
rng(7);
distancias = randi([1 100], n_ciudades, n_ciudades);
distancias = triu(distancias) + triu(distancias, 1)';
coordenadas = rand(n_ciudades, 2) * 100;

% Poblacion inicial
poblacion = zeros(tam_poblacion, n_ciudades);
for i = 1:tam_poblacion
    poblacion(i, :) = randperm(n_ciudades);
end
disp("Poblacion inicial"); disp(poblacion);

% Algoritmos geneticos
historial = zeros(max_generaciones, 1);
mejores_distancias = zeros(max_generaciones, 1);
umbral_distancia = 250; % Umbral para considerar una solucion buena
for g = 1:max_generaciones
    % Evaluar aptitud
    fitness = zeros(tam_poblacion, 1);
    soluciones_generacion = [];  % Para almacenar soluciones buenas
    for i = 1:tam_poblacion
        fitness(i) = fitness_tsp(poblacion(i, :), distancias);
        % Si la distancia es menor al umbral, guardar como solucion
        d = (1/fitness(i)) -1;
        if d <= umbral_distancia
            soluciones_generacion = [soluciones_generacion; poblacion(i,:) d];
        end
    end
    [mejor_fit, idx] = max(fitness);
    historial(g) = mejor_fit;
    
    % Mostrar progreso
    if mod(g, 50) == 0
        fprintf('\nGeneración %d: Soluciones optimas porque tiene distancias <= %.2f:\n', g, umbral_distancia);
        if ~isempty(soluciones_generacion)
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
                fprintf('... (+%d mas)\n', size(soluciones_ordenadas, 1)-5);
            end
        else
            fprintf('No se encontraron soluciones bajo el umbral en esta generacion\n');
        end
    end
    
    % Seleccion por torneo
    padres = seleccion_torneo(poblacion, fitness);
    
    % Cruce OX
    descendencia = cruceOX(padres, prob_cruce, n_ciudades);
    
    % Mutacion
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

% Grafica de convergencia
figure; plot(1:max_generaciones, historial, 'LineWidth', 2);
xlabel('Generacion'); ylabel('fitness'); title('Convergencia'); grid on;

% Grafica de ruta
mejor_ruta = [mejor_ind mejor_ind(1)];
figure; 
plot(coordenadas(mejor_ruta, 1), coordenadas(mejor_ruta, 2), '-or');
hold on;

% Añadir etiquetas de ciudades y distancias
for i = 1:n_ciudades
    % Etiquetas de ciudades
    text(coordenadas(i,1), coordenadas(i,2), sprintf('Ciudad %d', i), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    
    % Mostrar distancias entre ciudades conectadas
    if i < length(mejor_ruta)
        ciudad_actual = mejor_ruta(i);
        ciudad_sig = mejor_ruta(i+1);
        x1 = coordenadas(ciudad_actual, 1);
        y1 = coordenadas(ciudad_actual, 2);
        x2 = coordenadas(ciudad_sig, 1);
        y2 = coordenadas(ciudad_sig, 2);
        
        % Punto medio para colocar la distancia
        xm = (x1 + x2)/2;
        ym = (y1 + y2)/2;
        
        % Mostrar distancia
        distancia = distancias(ciudad_actual, ciudad_sig);
        text(xm, ym, sprintf('%.1f', distancia), ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', ...
            'Margin', 1, 'FontSize', 8);
    end
end

title(sprintf('Ruta Final - Distancia Total: %.2f', 1/mejor_fit - 1)); 
xlabel('X'); ylabel('Y'); 
axis equal; grid on;
hold off;
% Funcion de aptitud
function fit = fitness_tsp(ind, dist)
    d = 0;
    for i = 1:length(ind)-1
        d = d + dist(ind(i), ind(i+1));
    end
    d = d + dist(ind(end), ind(1));
    fit = 1 / (1 + d);
end

% Seleccion por torneo
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

% Mutacion por intercambio
function ind = mutacion_intercambio(ind)
    if length(ind) > 1
        pos = randperm(length(ind), 2);
        temp = ind(pos(1));
        ind(pos(1)) = ind(pos(2));
        ind(pos(2)) = temp;
    end
end