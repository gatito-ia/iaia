% Algoritmo genético para Optimización de Ubicacion de Antenas en Telecomunicaciones
clear all; close all; clc;

% Parámetros del algoritmo genético
parametros.tam_poblacion = 50;
parametros.max_generaciones = 100;
parametros.prob_cruce = 0.8;
parametros.prob_mutacion = 0.1;

disp("a) Evaluación de aptitud: ");
disp("b) Selección: torneo");
disp("c) Cruce: un punto");
disp("d) Mutación: gaussiana");

problema.numAntenas = 3;
problema.numTargets = 50;
problema.areaSize = 40;
problema.targets = problema.areaSize * rand(problema.numTargets, 2);
problema.maxRange = [10 6 14];
problema.funcion_fitness = @CoberturAntenas;

problema.tam_cromosoma = problema.numAntenas * 2;
problema.limite_inf = [0 0 0 0 0 0];
problema.limite_sup = problema.areaSize * [1 1 1 1 1 1];
problema.sigma = 0.5;

[mejorIndividuo, mejor_fitness] = algoritmo_genetico(problema, parametros);

%% ====================================
antenasMEjores = reshape(mejorIndividuo, 2, problema.numAntenas)';

fprintf('\n=== UBICACIONES ÓPTIMAS DE ANTENAS ===\n');
for i = 1:problema.numAntenas
    fprintf('Antena %d: (%.3f, %.3f)\n', i, antenasMEjores(i,1), antenasMEjores(i,2));
end

% Visualización
figure; hold on;
for i = 1:problema.numAntenas
    theta = linspace(0, 2*pi, 100);
    plot(antenasMEjores(i,1) + problema.maxRange(i)*cos(theta), ...
         antenasMEjores(i,2) + problema.maxRange(i)*sin(theta), 'r--');
end
plot(problema.targets(:,1), problema.targets(:,2), 'kx');
plot(antenasMEjores(:,1), antenasMEjores(:,2), 'or');
% Añadir etiquetas de ciudades
for i = 1:problema.numAntenas
    text(antenasMEjores(i,1), antenasMEjores(i,2), sprintf('Antena %d', i));
end
xlim([0 problema.areaSize]); ylim([0 problema.areaSize]); grid on; axis equal;
title('Cobertura de antenas óptimas');

%% ================== FUNCIÓN PRINCIPAL DEL AG ==================

function [mejor_i, mejor_apti] = algoritmo_genetico(problema, parametros)
    %% PARÁMETROS DEL ALGORITMO GENÉTICO
    n = parametros.tam_poblacion;
    Generacion = parametros.max_generaciones;
    
    %% Población inicial
    poblacion = zeros(n, problema.tam_cromosoma);
    for k = 1:n
        for r = 1:problema.tam_cromosoma
            poblacion(k, r)=problema.limite_inf(r) + rand() * (problema.limite_sup(r) - problema.limite_inf(r));
        end
    end
    disp("Población inicial"); disp(poblacion);
    
    %% EVOLUCIÓN
    for gen = 1:Generacion
        % Evaluar aptitud y almacenar soluciones óptimas
        aptitud = zeros(n, 1);
        soluciones_optimas = [];  % Almacenar soluciones óptimas de esta generación
        
        for i = 1:n
            aptitud(i) = evaluar_poblacion(poblacion(i,:), problema);
            
            % Consideramos solución óptima
            if aptitud(i) >= 0.6
                soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
            end
        end
        
        [mejor_fitness_gen, idx_mejor] = max(aptitud);
        elite = poblacion(idx_mejor,:);
        
        if mod(gen, 20) == 0
            % Mostrar progreso y soluciones óptimas encontradas
            if ~isempty(soluciones_optimas)
                fprintf('\nGeneración %d: Soluciones óptimas encontradas\n', gen);
                % Mostrar solo soluciones únicas
                soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
                
                % Limitar a máximo 5 soluciones
                num_mostrar = min(5, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    sol = soluciones_unicas(s,:);
                    antenas = reshape(sol, 2, problema.numAntenas)';
                    fprintf('Solución %d:\n', s);
                    for a = 1:problema.numAntenas
                        fprintf('  Antena %d: (%.2f, %.2f)\n', a, antenas(a,1), antenas(a,2));
                    end
                    fprintf('Fitness: %.4f\n\n', evaluar_poblacion(sol, problema));
                end
                
                % Mostrar mensaje si hay más soluciones no mostradas
                if size(soluciones_unicas, 1) > 5
                    fprintf('... y %d soluciones más\n\n', size(soluciones_unicas, 1) - 5);
                end
            else
                fprintf('Generación %d: Mejor fitness = %.4f\n', gen, mejor_fitness_gen);
            end
        end

        descendencia = zeros(size(poblacion));
        
        % SELECCIÓN
        pp = zeros(n, problema.tam_cromosoma);
        k=3;
        for i = 1:n
            competidores = randperm(n, k);
            [~, idx] = max(aptitud(competidores));
            pp(i,:) = poblacion(competidores(idx),:);
        end

        for i = 2:2:n-1
            % CRUCE
            h1 = pp(i, :);
            h2 = pp(i+1, :);
            if rand() < parametros.prob_cruce
                [h1, h2] = UnPunto(h1, h2);
            end
            % MUTACIÓN
            h1 = mutar(h1, problema, parametros.prob_mutacion);
            h2 = mutar(h2, problema, parametros.prob_mutacion);
            
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
end

%% ================== FUNCIONES DE EVALUACIÓN ==================
function aptitu = evaluar_poblacion(p, problema)
    if size(p,1) > 1
        aptitu = zeros(size(p, 1), 1);
        for i = 1:size(p, 1)
            aptitu(i) = problema.funcion_fitness(p(i, :), problema);
        end
    else
        aptitu = problema.funcion_fitness(p, problema);
    end
end

%% ================== CRUCE ==================
function [h1, h2] = UnPunto(p1, p2)
    punto = randi(length(p1)-1);
    h1 = [p1(1:punto) p2(punto+1:end)];
    h2 = [p2(1:punto) p1(punto+1:end)];
end

%% ================== MUTACIÓN ==================
function individuo = mutar(ind, problema, prob)
    if rand < prob
        gen = randi(length(ind));
        ind(gen) = ind(gen) + problema.sigma*randn();
        ind(gen) = min(max(ind(gen), problema.limite_inf(gen)), problema.limite_sup(gen));
    end
    individuo = ind;
end

%% ======= FITNESS====
function aptitu = CoberturAntenas(i, problema)
    antenas = reshape(i, 2, problema.numAntenas)';
    objCubiertos = false(problema.numTargets, 1);

    for i = 1:problema.numTargets
        objeto = problema.targets(i, :);
        for j = 1:problema.numAntenas
            antena = antenas(j, :);
            d = norm(objeto - antena);
            if d <= problema.maxRange(j)
                objCubiertos(i) = true;
                break;
            end
        end
    end
    aptitu = sum(objCubiertos) / problema.numTargets;
end
