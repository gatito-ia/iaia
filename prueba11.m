%% Canales
clear all; close all; clc;

%% Parametros del problema y algoritmo
canales_totales = 15; % Numero total de canales disponibles
celdas = 10;  % Numero de celdas a asignar

% Configuracion del algoritmo genetico
parametros = struct();
parametros.tam_poblacion = 50;       % Tamano de la poblacion
parametros.max_generaciones = 100;   % Numero maximo de generaciones
parametros.prob_cruce = 0.8;         % Probabilidad de cruce
parametros.prob_mutacion = 0.1;      % Probabilidad de mutacion

disp("Configuracion del algoritmo:");
disp("a) Evaluacion de aptitud: Minimizar interferencia");
disp("b) Seleccion: Torneo");
disp("c) Cruce: Uniforme");
disp("d) Mutacion: Intercambio");

%% Definicion del problema
problema = struct();
problema.nombre = 'Asignacion optima de canales';
problema.tam_cromosoma = celdas;
problema.funcion_fitness = @(individuo) fitness_canal(individuo, canales_totales);
problema.limite_inf = ones(1, celdas);
problema.limite_sup = canales_totales*ones(1, celdas);

%% Ejecucion del algoritmo genetico
[mejor_asignacion, mejor_fitness, historico_fitness] = algoritmo_genetico(problema, parametros);

%% Resultados
fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor fitness (minimizar interferencia): %.4f\n', mejor_fitness);
fprintf('Asignación óptima de canales:\n');
disp(mejor_asignacion);

% Grafico de evolucion
figure;
plot(1:parametros.max_generaciones, historico_fitness, 'LineWidth', 2);
title('Evolucion del Fitness');
xlabel('Generación');
ylabel('Mejor Fitness');
grid on;

%% ========= funcion ag ==========
function [mejor_ind, mejor_fit, historico_fitness] = algoritmo_genetico(problema, parametros)
    % Inicializacion
    n = parametros.tam_poblacion;
    max_gen = parametros.max_generaciones;
    historico_fitness = zeros(max_gen, 1);
    
    % Crear poblacion inicial aleatoria
    poblacion = randi([problema.limite_inf(1), problema.limite_sup(1)], n, problema.tam_cromosoma);
    disp("Poblacion inicial"); disp(poblacion);
    
    for gen = 1:max_gen
        % Evaluar fitness y almacenar soluciones optimas
        fitness = zeros(n, 1);
        soluciones_optimas = [];
        
        for i = 1:n
            fitness(i) = problema.funcion_fitness(poblacion(i,:));
            % Consideramos solucion optima si no hay interferencia
            if fitness(i) >=0.03
                soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
            end
        end
        
        [mejor_fit, idx] = max(fitness);
        elite = poblacion(idx,:);
        historico_fitness(gen) = mejor_fit;
        
        if mod(gen, 20) == 0
            % Mostrar soluciones optimas encontradas
            if ~isempty(soluciones_optimas)
                fprintf('\nGeneracion %d: Soluciones optimas encontradas (Fitness = %.4f):\n', gen, mejor_fit);
                soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    fprintf('Solucion %d: ', s);
                    disp(soluciones_unicas(s,:));
                    fprintf("Es una solucion porque tiene %d interferencias | %d variedad canales\n\n",...
                        interferencias_f(soluciones_unicas(s,:)),variedad_canales(soluciones_unicas(s,:)));
                end
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones mas\n', size(soluciones_unicas, 1) - 6);
                end
            end
        end

        % Seleccion por torneo
        padres = zeros(size(poblacion));
        for i = 1:n
            k=3;
            competidores = randperm(n, k);
            [~, idx_ganador] = max(fitness(competidores));
            padres(i,:) = poblacion(competidores(idx_ganador),:);
        end
        
        % Cruce uniforme
        descendencia = padres;
        for i = 1:2:n-1
            if rand() < parametros.prob_cruce
                mascara = rand(1, problema.tam_cromosoma) > 0.5;
                temp1 = descendencia(i,:);
                temp2 = descendencia(i+1,:);
                descendencia(i,mascara) = temp2(mascara);
                descendencia(i+1,mascara) = temp1(mascara);
            end
        end
        
        % Mutacion por intercambio
        for i = 1:n
            if rand() < parametros.prob_mutacion
                pos = randperm(problema.tam_cromosoma, 2);
                temp = descendencia(i,pos(1));
                descendencia(i,pos(1)) = descendencia(i,pos(2));
                descendencia(i,pos(2)) = temp;
            end
        end
        
        % Elitismo
        descendencia(1,:) = elite;
        poblacion = descendencia;
    end
    
    % Resultado final
    mejor_ind = elite;
end

%% ========= FUNCIÓN DE FITNESS ==========
function fitness = fitness_canal(asignacion, canales_totales)
    % Penalizar interferencia entre celdas adyacentes
    interferencia = 0;
    for i = 2:length(asignacion)
        dif = abs(asignacion(i) - asignacion(i-1));
        if dif == 0 % Misma frecuencia
            interferencia = interferencia + 100; 
        else
            interferencia = interferencia +1/dif; % Menor diferencia
        end
    end
    
    % Penalizar uso excesivo de canales
    canales_usados = length(unique(asignacion));
    penalizacion = (canales_totales - canales_usados)^2;
    
    % Fitness inverso a la interferencia total
    fitness = 1/(interferencia + penalizacion + 0.001);
end
%% ========= adicionales ==========
function interferencia = interferencias_f(asignacion)
    % Penalizar interferencia entre celdas adyacentes
    interferencia = 0;
    for i = 2:length(asignacion)
        dif = abs(asignacion(i) - asignacion(i-1));
        if dif == 0 % Misma frecuencia
            interferencia = interferencia + 100; 
        else
            interferencia = interferencia +1/dif; % Menor diferencia
        end
    end
end
function canales_usados = variedad_canales(asignacion)
    % Penalizar uso excesivo de canales
    canales_usados = length(unique(asignacion));
end