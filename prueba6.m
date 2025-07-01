% Algoritmo genético para Optimización de Ubicacion de Antenas en % Telecomnicaciones
clear all; close all;clc;

% Parámetros del algoritmo genético
parametros.tam_poblacion = 50;
parametros.max_generaciones = 100;
parametros.prob_cruce = 0.8;
parametros.prob_mutacion = 0.1;

disp("a) Evalulación de aptitud: ");
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
    
    %% EVOLUCIÓN
    for gen = 1:Generacion
        aptitud = evaluar_poblacion(poblacion, problema);
        [mejor_fitness_gen, idx_mejor] = max(aptitud);
        elite = poblacion(idx_mejor,:);

        % Mostrar progreso
        if mod(gen, 10) == 0 || gen == 1
            fprintf('Generación %d: Mejor fitness = %.4f\n', gen, mejor_fitness_gen);
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
            h1 = mutar(h1, problema, rand());
            h2 = mutar(h2, problema, rand());
            
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
    
antenas= reshape(i, 2, problema.numAntenas)';
    objCubiertos = false(problema.numTargets, 1);

    for i = 1:problema.numTargets
        objeto = problema.targets(i, :);
        for j = 1:problema.numAntenas
            antena = antenas(j, :);
            d = norm(objeto - antena);
            if d <= problema.maxRange
                objCubiertos(i) = true;
            end
        end
    end
    aptitu = sum(objCubiertos) / problema.numTargets;
end

%% ======================================================
%% ======================================================

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
    
    % Bucle principal
    for gen = 1:max_gen
        % Evaluar fitness
        fitness = zeros(n, 1);
        for i = 1:n
            fitness(i) = problema.funcion_fitness(poblacion(i,:));
        end
        
        [mejor_fit, idx] = max(fitness);
        elite = poblacion(idx,:); % Guardar el mejor
        historico_fitness(gen) = mejor_fit;
        
        % Mostrar progreso cada 10 generaciones
        if mod(gen, 10) == 0
            fprintf('Generación %d - Mejor fitness: %.4f\n', gen, mejor_fit);
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