%% Parámetros del problema (Job Sequencing)
clear all; close all; clc;

% Datos del problema
n_jobs = 10;
profits = [100 19 27 25 15 90 30 10 40 70];
deadlines = [2 1 2 1 3 5 7 9 3 2];

% Parámetros del algoritmo
p.tam_poblacion = 50;
p.max_generaciones = 100;
p.prob_cruce = 0.8;
p.prob_mutacion = 0.1;

% Configuración del problema
prob.tam_cromosoma = n_jobs;
prob.profits = profits;
prob.deadlines = deadlines;

disp("a) Aptitud: suma de beneficios de tareas a tiempo");
disp("b) Selección: torneo");
disp("c) Cruce: orden (OX)");
disp("d) Mutación: intercambio");

% Ejecutar algoritmo genético
[best_seq, best_profit] = algoritmo_genetico(prob, p);

% Mostrar resultados
fprintf('\n=== Mejor Secuencia ===\n');
fprintf('Orden  Job  Profit  Deadline  Tiempo  Estado\n');
time = 0; total_profit = 0;
for i = 1:n_jobs
    job = best_seq(i);
    time = time + 1;
    status = 'A tiempo';
    if time <= deadlines(job)
        total_profit = total_profit + profits(job);
    else
        status = 'Tarde';
    end
    fprintf('%3d %4d %6d %8d %8d %10s\n', ...
            i, job, profits(job), deadlines(job), time, status);
end
fprintf('Beneficio total: %d\n', total_profit);

% === Algoritmo Genético ===
function [mejor, mejor_apt] = algoritmo_genetico(prob, p)
    % Población inicial
    n = p.tam_poblacion;
    poblacion = zeros(n, prob.tam_cromosoma);
    for i = 1:n
        poblacion(i,:) = randperm(prob.tam_cromosoma);
    end
    disp("Población inicial"); disp(poblacion);
    
    % Evolución
    for gen = 1:p.max_generaciones
        % Evaluar aptitud y almacenar soluciones válidas
        aptitud = zeros(n, 1);
        soluciones_generacion = [];  % Almacenar soluciones de esta generación
        
        for i = 1:n
            [aptitud(i), tiempo] = calcular_aptitud(poblacion(i,:), prob);
            if aptitud(i) >= 300
                soluciones_generacion = [soluciones_generacion; poblacion(i,:)];
            end
        end
        
        % Guardar el mejor
        [mejor_apt, idx] = max(aptitud);
        elite = poblacion(idx,:);
        
        if mod(gen, 20) == 0
        % Mostrar soluciones válidas encontradas
            if ~isempty(soluciones_generacion)
                fprintf('\nGeneración %d: Soluciones válidas encontradas:\n', gen);
                soluciones_unicas = unique(soluciones_generacion, 'rows');
                
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    seq = soluciones_unicas(s,:);
                    [apt, ~] = calcular_aptitud(seq, prob);
                    fprintf('Secuencia: %s - Beneficio: %d\n', mat2str(seq), apt);
                end
                
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones más\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generación %d: Mejor beneficio = %.2f\n', gen, mejor_apt);
            end
        end
        
        
        % Selección por torneo
        padres = zeros(n, prob.tam_cromosoma);
        k = 3;
        for i = 1:n
            competidores = randperm(n, k);
            [~, idx] = max(aptitud(competidores));
            padres(i,:) = poblacion(competidores(idx),:);
        end
        
        % Cruce OX
        hijos = padres;
        for i = 1:2:n-1
            if rand() < p.prob_cruce
                [hijos(i,:), hijos(i+1,:)] = cruce_ox(padres(i,:), padres(i+1,:));
            end
        end
        
        % Mutación por intercambio
        for i = 1:n
            if rand() < p.prob_mutacion
                pos = randperm(prob.tam_cromosoma, 2);
                hijos(i,[pos(1) pos(2)]) = hijos(i,[pos(2) pos(1)]);
            end
        end
        
        % Elitismo: conservar el mejor
        hijos(1,:) = elite;
        poblacion = hijos;
    end
    
    % Mejor solución final
    aptitud = zeros(n, 1);
    for i = 1:n
        aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
    end
    [mejor_apt, idx] = max(aptitud);
    mejor = poblacion(idx,:);
end

% === Función de Aptitud modificada ===
function [apt, tiempo_total] = calcular_aptitud(seq, prob)
    apt = 0;
    tiempo_total = 0;
    for i = 1:length(seq)
        tiempo_total = tiempo_total + 1;
        if tiempo_total <= prob.deadlines(seq(i))
            apt = apt + prob.profits(seq(i));
        end
    end
end

% === Cruce OX ===
function [h1, h2] = cruce_ox(p1, p2)
    n = length(p1);
    pts = sort(randperm(n-1, 2));
    
    % Copiar segmento central
    h1 = zeros(1, n); h2 = zeros(1, n);
    h1(pts(1):pts(2)) = p1(pts(1):pts(2));
    h2(pts(1):pts(2)) = p2(pts(1):pts(2));
    
    % Rellenar con genes restantes
    idx1 = pts(2)+1; idx2 = pts(2)+1;
    for i = 1:n
        pos = mod(pts(2)+i-1, n)+1;
        if ~ismember(p2(pos), h1(pts(1):pts(2)))
            if idx1 > n, idx1 = 1; end
            h1(idx1) = p2(pos);
            idx1 = idx1 + 1;
        end
        if ~ismember(p1(pos), h2(pts(1):pts(2)))
            if idx2 > n, idx2 = 1; end
            h2(idx2) = p1(pos);
            idx2 = idx2 + 1;
        end
    end
end
