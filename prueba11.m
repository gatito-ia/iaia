%% Canales
clear all; close all; clc;

%% Parámetros del problema y algoritmo
canales_totales = 15; % Número total de canales disponibles
celdas = 10;  % Número de celdas a asignar

% Configuración del algoritmo genético
parametros = struct();
parametros.tam_poblacion = 50;       % Tamaño de la población
parametros.max_generaciones = 100;   % Número máximo de generaciones
parametros.prob_cruce = 0.8;         % Probabilidad de cruce
parametros.prob_mutacion = 0.1;      % Probabilidad de mutación

disp("Configuración del algoritmo:");
disp("a) Evaluación de aptitud: Minimizar interferencia");
disp("b) Selección: Torneo");
disp("c) Cruce: Uniforme");
disp("d) Mutación: Intercambio");

%% Definición del problema
problema = struct();
problema.nombre = 'Asignación óptima de canales';
problema.tam_cromosoma = celdas;
problema.funcion_fitness = @(individuo) fitness_canal(individuo, canales_totales);
problema.limite_inf = ones(1, celdas);   % Mínimo canal (1)
problema.limite_sup = canales_totales*ones(1, celdas); % Máximo canal (15)

%% Ejecución del algoritmo genético
[mejor_asignacion, mejor_fitness, historico_fitness] = algoritmo_genetico(problema, parametros);

%% Resultados
fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor fitness (minimizar interferencia): %.4f\n', mejor_fitness);
fprintf('Asignación óptima de canales:\n');
disp(mejor_asignacion);

% Gráfico de evolución
figure;
plot(1:parametros.max_generaciones, historico_fitness, 'LineWidth', 2);
title('Evolución del Fitness');
xlabel('Generación');
ylabel('Mejor Fitness');
grid on;

%% ========= FUNCIÓN PRINCIPAL DEL AG ==========
function [mejor_ind, mejor_fit, historico_fitness] = algoritmo_genetico(problema, parametros)
    % Inicialización
    n = parametros.tam_poblacion;
    max_gen = parametros.max_generaciones;
    historico_fitness = zeros(max_gen, 1);
    
    % Crear población inicial aleatoria
    poblacion = randi([problema.limite_inf(1), problema.limite_sup(1)], n, problema.tam_cromosoma);
    disp("Población inicial"); disp(poblacion);
    
    % Bucle principal
    for gen = 1:max_gen
        % Evaluar fitness y almacenar soluciones óptimas
        fitness = zeros(n, 1);
        soluciones_optimas = [];  % Almacenar soluciones óptimas de esta generación
        
        for i = 1:n
            fitness(i) = problema.funcion_fitness(poblacion(i,:));
            
            % Consideramos solución óptima si no hay interferencia (fitness máximo)
            if fitness(i) <=0.035
                soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
            end
        end
        
        [mejor_fit, idx] = max(fitness);
        elite = poblacion(idx,:); % Guardar el mejor
        historico_fitness(gen) = mejor_fit;
        
        if mod(gen, 20) == 0
            % Mostrar soluciones óptimas encontradas
            if ~isempty(soluciones_optimas)
                fprintf('\nGeneración %d: Soluciones óptimas encontradas (Fitness = %.4f):\n', gen, mejor_fit);
                % Mostrar solo soluciones únicas
                soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
                
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    fprintf('Solución %d: ', s);
                    disp(soluciones_unicas(s,:));
                end
                
                % Indicar si hay más soluciones no mostradas
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones más\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generación %d - Mejor fitness: %.4f\n', gen, mejor_fit);
            end
        end
        
        
        % Selección por torneo
        padres = zeros(size(poblacion));
        for i = 1:n
            competidores = randperm(n, 3); % Torneo de 3
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
        
        % Mutación por intercambio
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
        if dif == 0 % Misma frecuencia -> máxima interferencia
            interferencia = interferencia + 100; 
        else
            interferencia = interferencia + 1/dif; % Menor diferencia -> mayor interferencia
        end
    end
    
    % Penalizar uso excesivo de canales
    canales_usados = length(unique(asignacion));
    penalizacion = (canales_totales - canales_usados)^2;
    
    % Fitness inverso a la interferencia total
    fitness = 1/(interferencia + penalizacion + 0.001);
end