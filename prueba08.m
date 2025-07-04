%% Parametros del problema (Job Sequencing)
clear all; close all; clc;

% Datos del problema
n_jobs = 10;
profits = randi([10, 100], n_jobs, 1); 
tiempos_procesamiento = randi([1, 7], n_jobs, 1);
deadlines = cumsum(tiempos_procesamiento) + randi([0, 6], n_jobs, 1);


% Parametros del algoritmo
p.tam_poblacion = 50;
p.max_generaciones = 100;
p.prob_cruce = 0.8;
p.prob_mutacion = 0.1;

% Configuracion del problema
prob.tam_cromosoma = n_jobs;
prob.profits = profits;
prob.deadlines = deadlines;
prob.tiempos_procesamiento = tiempos_procesamiento;

disp("a) Aptitud: suma de beneficios de tareas a tiempo");
disp("b) Seleccion: torneo");
disp("c) Cruce: orden (OX)");
disp("d) Mutacion: intercambio");

% Ejecutar algoritmo genetico
[mejor_seq, mejor_profit] = algoritmo_genetico(prob, p);

% Mostrar resultados
fprintf('\n=== Mejor Secuencia ===\n');
fprintf('Orden  Job  Profit  Deadline  T_procesa Tiempo  Estado\n');
time = 0; total_profit = 0;
for i = 1:n_jobs
    job = mejor_seq(i);
    time = time + prob.tiempos_procesamiento(job);
    status = 'A tiempo';
    if time <= deadlines(job)
        total_profit = total_profit + profits(job);
    else
        status = 'Tarde';
    end
    fprintf('%3d %4d %6d %8d %10d %8d %10s\n', ...
            i, job, profits(job), deadlines(job), prob.tiempos_procesamiento(job), time, status);
end
fprintf('Beneficio total: %d\n', total_profit);

% === Algoritmo Genetico ===
function [mejor, mejor_apt] = algoritmo_genetico(prob, p)
    % Poblacion inicial
    n = p.tam_poblacion;
    poblacion = zeros(n, prob.tam_cromosoma);
    for i = 1:n
        poblacion(i,:) = randperm(prob.tam_cromosoma);
    end
    disp("Poblacion inicial"); disp(poblacion);
    
    % Evolucion
    for gen = 1:p.max_generaciones
        % Evaluar aptitud y almacenar soluciones validas
        aptitud = zeros(n, 1);
        soluciones_generacion = []; 
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
        % Mostrar soluciones validas encontradas
            if ~isempty(soluciones_generacion)
                fprintf('\nGeneracion %d: Soluciones validas encontradas:\n', gen);
                soluciones_unicas = unique(soluciones_generacion, 'rows');
                
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    seq = soluciones_unicas(s,:);
                    [apt, ~] = calcular_aptitud(seq, prob);
                    fprintf('Secuencia: %s - Beneficio: %d\n', mat2str(seq), apt);
                    time = 0;k=0;
                    for i=1:length(seq)
                        job = seq(i);
                        time = time + prob.tiempos_procesamiento(job);
                        if time <= prob.deadlines(job)
                            k=k+1;
                        end
                    end
                    fprintf('Es una solucion optima porque tiene %d tareas cumplidas\n', k);
                end
                
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones mas\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generacion %d: Mejor beneficio = %.2f\n', gen, mejor_apt);
            end
        end
        
        
        % Seleccion por torneo
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
        
        % Mutacion por intercambio
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
    
    % Mejor solucion final
    aptitud = zeros(n, 1);
    for i = 1:n
        aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
    end
    [mejor_apt, idx] = max(aptitud);
    mejor = poblacion(idx,:);
end

% === Funcion de Aptitud ===
function [apt, tiempo_total] = calcular_aptitud(indivi, prob)
    apt = 0;
    tiempo_total = 0;
    for i = 1:length(indivi)
        tiempo_total = tiempo_total + prob.tiempos_procesamiento(i);
        if tiempo_total <= prob.deadlines(indivi(i))
            apt = apt + prob.profits(indivi(i));
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
