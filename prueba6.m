%% Algoritmo genético para el Problema de la Mochila
% Adaptado del código de optimización de ubicación de antenas
clear all; close all; clc;

% Parámetros del problema de la mochila
disp("Condiciones del problema: Pesos");
pesos = [10, 20, 30, 15, 25, 5, 35, 12, 22, 18];
disp(pesos);
disp("Condiciones del problema: valores");
valores = [60, 100, 120, 70, 90, 30, 150, 50, 80, 110];
disp(valores);
disp("Condiciones del problema: Capacidad");
capacidad_maxima = 100;
disp(capacidad_maxima);

% Parámetros del algoritmo genético
parametros.tam_poblacion = 100;
parametros.max_generaciones = 1000;
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
    disp("Población inicial"); disp(poblacion);
    %% EVOLUCIÓN
    for gen = 1:Generacion
        % Evaluar población y detectar soluciones válidas
        aptitud = zeros(n, 1);
        soluciones_generacion = [];  % Para almacenar soluciones válidas
        
        for i = 1:n
            % Calcular fitness
            valor_total = sum(problema.valores(poblacion(i,:) == 1));
            peso_total = sum(problema.pesos(poblacion(i,:) == 1));
            
            % Determinar si es solución válida (no excede capacidad)
            if peso_total <= problema.capacidad_maxima
                aptitud(i) = valor_total;
            else
                aptitud(i) = 0; % Penalización por exceder capacidad
            end
            if peso_total <= problema.capacidad_maxima & peso_total >= 90
                soluciones_generacion = [soluciones_generacion; poblacion(i,:) valor_total peso_total];
            end
        end
        
        [mejor_fitness_gen, idx_mejor] = max(aptitud);
        elite = poblacion(idx_mejor,:);
        
        if mod(gen, 100) == 0
            % Mostrar soluciones válidas encontradas
            if ~isempty(soluciones_generacion)
                fprintf('\nGeneración %d: Mejores soluciones válidas encontradas:\n', gen);
                
                % Eliminar duplicados y ordenar por valor (de mayor a menor)
                [soluciones_unicas, idx_unicas] = unique(soluciones_generacion(:,1:problema.tam_cromosoma), 'rows', 'stable');
                valores = soluciones_generacion(idx_unicas, end-1);  % Obtener valores de las soluciones únicas

                % Mostrar las 5 mejores
                num_mostrar = min(5, size(soluciones_unicas, 1));
                
                for j = 1:num_mostrar
                    idx_original = idx_unicas(j);
                    fprintf('Solución %d: %s | Valor: %.1f | Peso: %.1f\n', j, ...
                        mat2str(soluciones_generacion(idx_original,1:problema.tam_cromosoma)), ...
                        soluciones_generacion(idx_original,end-1), ...
                        soluciones_generacion(idx_original,end));
                end
                
                % Mostrar mensaje si hay más soluciones
                if size(soluciones_unicas, 1) > 5
                    fprintf('... y %d soluciones más (Valor mínimo mostrado: %.1f)\n', ...
                        size(soluciones_unicas, 1) - 5, soluciones_unicas(5));
                end
            else
                fprintf('Generación %d: Mejor fitness = %.4f\n', gen, mejor_fitness_gen);
            end
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
end

%% ================== FUNCIONES DE EVALUACIÓN ==================
function aptitu = evaluar_poblacion(p, problema)
    aptitu = zeros(size(p, 1), 1);
    for i = 1:size(p, 1)
        valor_total = sum(problema.valores(i == 1));
        peso_total = sum(problema.pesos(i == 1));
        if peso_total > problema.capacidad_maxima
            aptitu(i) = 0;
        else
            aptitu(i) = valor_total;
        end
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