%% Algoritmo genético para el Problema de la Mochila
% Adaptado del código de optimización de ubicación de antenas
clear all; close all; clc;

% Parámetros del problema de la mochila
pesos = [10, 20, 30, 15, 25, 5, 35, 12, 22, 18];
valores = [60, 100, 120, 70, 90, 30, 150, 50, 80, 110];
capacidad_maxima = 100;

% Parámetros del algoritmo genético
parametros.tam_poblacion = 100;
parametros.max_generaciones = 200;
parametros.prob_cruce = 0.8;
parametros.prob_mutacion = 0.02;

disp("a) Evaluación de aptitud: función de valor con penalización por exceso de peso");
disp("b) Selección: torneo");
disp("c) Cruce: un punto");
disp("d) Mutación: bit flip");

% Configuración del problema
problema.tam_cromosoma = length(pesos);
problema.limite_inf = zeros(1, problema.tam_cromosoma);
problema.limite_sup = ones(1, problema.tam_cromosoma);
problema.sigma = 0.1;

problema.pesos = pesos;
problema.valores = valores;
problema.capacidad_maxima = capacidad_maxima;
problema.funcion_fitness = @fitness_mochila;

[mejorIndividuo, mejor_fitness] = algoritmo_genetico(problema, parametros);

%% Mostrar resultados
fprintf('\n=== RESULTADOS MOCHILA ===\n');
fprintf('Mejor solución: %s\n', mat2str(mejorIndividuo));
fprintf('Valor total: %.2f\n', mejor_fitness);
fprintf('Peso total: %.2f\n', sum(pesos(mejorIndividuo==1)));
fprintf('Ítems seleccionados: %s\n', mat2str(find(mejorIndividuo==1)));

%% ================== FUNCIÓN PRINCIPAL DEL AG ==================

function [mejor_i, mejor_apti] = algoritmo_genetico(problema, parametros)
    %% PARÁMETROS DEL ALGORITMO GENÉTICO
    n = parametros.tam_poblacion;
    Generacion = parametros.max_generaciones;
    
    %% Población inicial
    poblacion = zeros(n, problema.tam_cromosoma);
    for k = 1:n
        poblacion(k, :) = randi([0 1], 1, problema.tam_cromosoma); % Población binaria
    end
    
    %% EVOLUCIÓN
    for gen = 1:Generacion
        aptitud = evaluar_poblacion(poblacion, problema);
        [mejor_fitness_gen, idx_mejor] = max(aptitud);
        elite = poblacion(idx_mejor,:);

        % Mostrar progreso
        if mod(gen, 10) == 0 || gen == 1
            fprintf('Generación %d: Mejor fitness = %.4f\n', gen, mejor_fitness_gen);
        end
        
        % SELECCIÓN
        pp = Torneo(poblacion, aptitud);
        descendencia = zeros(size(poblacion));

        for i = 2:2:n-1
            % CRUCE
            h1 = pp(i, :);
            h2 = pp(i+1, :);

            if rand() < parametros.prob_cruce
                [h1, h2] = UnPunto(h1, h2);
            end
            
            % MUTACIÓN
            if rand() < parametros.prob_mutacion
                h1 = mutacionBitFlip(h1, problema);
            end
            if rand() < parametros.prob_mutacion
                h2 = mutacionBitFlip(h2, problema);
            end
            
            descendencia(i,:) = h1;
            descendencia(i+1,:) = h2;
        end
        
        % Elitismo
        descendencia(1,:) = elite;
        poblacion = descendencia;
    end
    
    % Evaluar población final
    fitnessFinal = evaluar_poblacion(poblacion, problema);
    [mejor_apti, idx_mejor] = max(fitnessFinal);
    mejor_i = poblacion(idx_mejor, :);
    
    fprintf('\n=== RESULTADO FINAL ===\n');
    fprintf('Mejor fitness: %.6f\n', mejor_apti);
    fprintf('Mejor individuo: ');
    disp(mejor_i);
end

%% ================== FUNCIONES DE EVALUACIÓN ==================
function aptitu = evaluar_poblacion(p, problema)
    aptitu = zeros(size(p, 1), 1);
    for i = 1:size(p, 1)
        aptitu(i) = problema.funcion_fitness(p(i, :), problema);
    end
end

%% ================== SELECCIÓN ==================
function p = Torneo(poblacion, fitness)
    tamTorneo = 3;
    n = size(poblacion, 1);
    p = zeros(n, size(poblacion, 2));
    for k = 1:n
        competidores = randperm(n, tamTorneo);
        aptituCompetidores = fitness(competidores);
        [~, ubica_ganador] = max(aptituCompetidores);
        ganador = competidores(ubica_ganador);
        p(k,:) = poblacion(ganador,:);
    end
end

%% ================== CRUCE ==================
function [h1, h2] = UnPunto(p1, p2)
    tamCrom = length(p1);
    pCorte = randi([1, tamCrom-1]);
    h1 = [p1(1:pCorte), p2(pCorte+1:end)];
    h2 = [p2(1:pCorte), p1(pCorte+1:end)];
end

%% ================== MUTACIÓN ==================
function i_Mutado = mutacionBitFlip(i, ~)
    i_Mutado = i;
    pos = randi(length(i));
    i_Mutado(pos) = 1 - i_Mutado(pos); % Flip del bit
end

%% ================== FITNESS MOCHILA ==================
function aptitu = fitness_mochila(i, problema)
    % Calcula el valor total
    valor_total = sum(problema.valores(i == 1));
    
    % Calcula el peso total
    peso_total = sum(problema.pesos(i == 1));
    
    % Aplica penalización si excede la capacidad
    if peso_total > problema.capacidad_maxima
        aptitu = 0;
    else
        aptitu = valor_total;
    end
end


%% =======================================================================
%% =======================================================================
%% =======================================================================
%% Viajero

% Algoritmo genético para TSP (simplificado y corregido)
clear all; close all; clc;

% Parámetros
tam_poblacion = 100;
max_generaciones = 500;
prob_cruce = 0.8;
prob_mutacion = 0.2;
n_ciudades = 10;

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

% Ciclo principal
historial = zeros(max_generaciones, 1);
for g = 1:max_generaciones
    % Evaluar aptitud
    fitness = zeros(tam_poblacion, 1);
    for i = 1:tam_poblacion
        fitness(i) = fitness_tsp(poblacion(i, :), distancias);
    end
    [mejor_fit, idx] = max(fitness);
    historial(g) = mejor_fit;
    
    % Mostrar progreso
    if mod(g, 50) == 0 || g == 1
        fprintf('Generación %d: Distancia = %.2f\n', g, historial(g));
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
    
    % Elitismo
    descendencia(1, :) = poblacion(idx, :);
    
    % Reemplazo
    poblacion = descendencia;
end

% Resultados
[mejor_fit, idx] = max(fitness);
mejor_ind = poblacion(idx, :);
fprintf('Mejor ruta: %s\nDistancia: %.2f\n', mat2str(mejor_ind), 1/mejor_fit - 1);

% Gráfica de convergencia
figure; plot(1:max_generaciones, historial, 'LineWidth', 2);
xlabel('Generación'); ylabel('fitness'); title('Convergencia'); grid on;

% Gráfica de ruta
mejor_ruta = [mejor_ind mejor_ind(1)];
figure; plot(coordenadas(mejor_ruta, 1), coordenadas(mejor_ruta, 2), '-o');
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
% Cruce por orden (OX) para TSP
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