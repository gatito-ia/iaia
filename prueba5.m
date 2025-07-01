%% Parámetros del problema (Job Sequencing)

% Algoritmo genético para optimizar secuencia de tareas
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
    
    % Evolución
    for gen = 1:p.max_generaciones
        % Evaluar aptitud
        aptitud = zeros(n, 1);
        for i = 1:n
            aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
        end
        
        % Guardar el mejor
        [mejor_apt, idx] = max(aptitud);
        elite = poblacion(idx,:);
        
        % Mostrar progreso cada 10 generaciones
        if mod(gen, 10) == 0
            fprintf('Generación %d: Beneficio = %.2f\n', gen, mejor_apt);
        end
        
        % Selección por torneo
        padres = zeros(n, prob.tam_cromosoma);
        k=3;
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

% === Función de Aptitud ===
function apt = calcular_aptitud(seq, prob)
    apt = 0;
    time = 0;
    for i = 1:length(seq)
        time = time + 1;
        if time <= prob.deadlines(seq(i))
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

%% ============================================================
%% ============================================================
%% ============================================================

%% Parámetros del problema Himmelblau

% Algoritmo genético para optimizar la función de Himmelblau con dos restricciones
clear all; close all; clc;

% Parámetros del problema
prob.tam_cromosoma = 2; % [x, y]
prob.limites = [-5, 5]; % Dominio para x e y
prob.lambda = 1000; % Factor de penalización para restricciones

% Parámetros del algoritmo genético
p.tam_poblacion = 50;
p.max_generaciones = 100;
p.prob_cruce = 0.8;
p.prob_mutacion = 0.1;
p.sigma_mutacion = 0.1; % Desviación estándar para mutación gaussiana
p.tam_torneo = 3; % Tamaño del torneo

disp("a) Aptitud: función de Himmelblau con penalización por dos restricciones");
disp("b) Selección: torneo");
disp("c) Cruce: aritmético");
disp("d) Mutación: gaussiana");

% Ejecutar algoritmo genético
[best_sol, best_apt, apt_por_generacion] = algoritmo_genetico(prob, p);

% Mostrar resultados
fprintf('\n=== Mejor Solución ===\n');
fprintf('x: %.4f, y: %.4f\n', best_sol(1), best_sol(2));
fprintf('Valor de Himmelblau: %.4f\n', himmelblau(best_sol));
fprintf('Restricción 1 (x^2 + y^2 - 4): %.4f (Factible si <= 0)\n', constraint1(best_sol));
fprintf('Restricción 2 (x + y - 1): %.4f (Factible si >= 0)\n', constraint2(best_sol));
fprintf('Aptitud final: %.4f\n', best_apt);

% Graficar la evolución de la aptitud por generación
figure;
plot(1:p.max_generaciones, apt_por_generacion, 'LineWidth', 2);
xlabel('Generación');
ylabel('Mejor Aptitud');
title('Evolución de la Aptitud por Generación');
grid on;


% Graficar la función de Himmelblau en 3D
figure;
[x, y] = meshgrid(linspace(-5, 5, 100), linspace(-5, 5, 100)); % Crear malla de puntos
z = (x.^2 + y - 11).^2 + (x + y.^2 - 7).^2; % Evaluar Himmelblau en cada punto

% Gráfica de superficie
surf(x, y, z, 'EdgeColor', 'none');
title('Función de Himmelblau: (x^2 + y - 11)^2 + (x + y^2 - 7)^2');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
colorbar;

hold on;
plot3(best_sol(1), best_sol(2), himmelblau(best_sol), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

% === Algoritmo Genético ===
function [mejor, mejor_apt, apt_por_generacion] = algoritmo_genetico(prob, p)
    % Población inicial
    n = p.tam_poblacion;
    poblacion = zeros(n, prob.tam_cromosoma);
    for i = 1:n
        poblacion(i,:) = prob.limites(1) + (prob.limites(2) - prob.limites(1)) * rand(1, prob.tam_cromosoma);
    end
    
    % Inicializar vector para almacenar la mejor aptitud por generación
    apt_por_generacion = zeros(p.max_generaciones, 1);
    
    % Evolución
    for gen = 1:p.max_generaciones
        % Evaluar aptitud
        aptitud = zeros(n, 1);
        for i = 1:n
            aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
        end
        
        % Guardar el mejor
        [mejor_apt, idx] = min(aptitud); % Minimización
        elite = poblacion(idx,:);
        
        % Almacenar la mejor aptitud de esta generación
        apt_por_generacion(gen) = mejor_apt;
        
        % Mostrar progreso cada 10 generaciones
        if mod(gen, 10) == 0
            fprintf('Generación %d: Aptitud = %.4f\n', gen, mejor_apt);
        end
        
        % Selección por torneo
        padres = zeros(n, prob.tam_cromosoma);
        for i = 1:n
            competidores = randperm(n, p.tam_torneo);
            [~, idx] = min(aptitud(competidores)); % Mejor competidor
            padres(i,:) = poblacion(competidores(idx),:);
        end
        
        % Cruce aritmético
        hijos = padres;
        for i = 1:2:n-1
            if rand() < p.prob_cruce
                alpha = rand();
                hijos(i,:) = alpha * padres(i,:) + (1-alpha) * padres(i+1,:);
                hijos(i+1,:) = (1-alpha) * padres(i,:) + alpha * padres(i+1,:);
            end
        end
        
        % Mutación gaussiana
        for i = 1:n
            if rand() < p.prob_mutacion
                hijos(i,:) = hijos(i,:) + p.sigma_mutacion * randn(1, prob.tam_cromosoma);
                hijos(i,:) = max(min(hijos(i,:), prob.limites(2)), prob.limites(1));
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
    [mejor_apt, idx] = min(aptitud);
    mejor = poblacion(idx,:);
end

% === Función de Himmelblau ===
function val = himmelblau(sol)
    x = sol(1); y = sol(2);
    val = (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
end

% === Restricción 1: x^2 + y^2 <= 4 ===
function val = constraint1(sol)
    x = sol(1); y = sol(2);
    val = x^2 + y^2 - 4; % <= 0
end

% === Restricción 2: x + y >= 1 ===
function val = constraint2(sol)
    x = sol(1); y = sol(2);
    val = 1 - (x + y); % <= 0 (convertido a forma estándar)
end

% === Función de Aptitud ===
function apt = calcular_aptitud(sol, prob)
    obj = himmelblau(sol);
    cons1 = constraint1(sol);
    cons2 = constraint2(sol);
    penalidad = 0;
    if cons1 > 0
        penalidad = penalidad + prob.lambda * cons1^2;
    end
    if cons2 > 0
        penalidad = penalidad + prob.lambda * cons2^2;
    end
    apt = obj + penalidad;
end