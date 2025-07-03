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
fprintf('\n=== Resultados finales: Mejor Solución ===\n');
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
title('Función de Himmelblau');
xlabel('x'); ylabel('y'); zlabel('f(x, y)'); colorbar;
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
    disp("Población inicial"); disp(poblacion);

    % Inicializar vector para almacenar la mejor aptitud por generación
    apt_por_generacion = zeros(p.max_generaciones, 1);
    
    % Evolución
    for gen = 1:p.max_generaciones
        % Evaluar aptitud y almacenar soluciones factibles
        aptitud = zeros(n, 1);
        soluciones_factibles = [];  % Almacenar soluciones factibles de esta generación
        
        for i = 1:n
            aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
            
            % Verificar si es solución factible (cumple ambas restricciones)
            c1 = constraint1(poblacion(i,:));
            c2 = constraint2(poblacion(i,:));
            if c1 <= 0 && c2 <= 0
                soluciones_factibles = [soluciones_factibles; poblacion(i,:)];
            end
        end
        
        % Guardar el mejor
        [mejor_apt, idx] = min(aptitud); % Minimización
        elite = poblacion(idx,:);
        apt_por_generacion(gen) = mejor_apt;
        
        if mod(gen, 30) == 0
            % Mostrar todas las soluciones
            if ~isempty(soluciones_factibles)
                fprintf('\nGeneración %d: Soluciones factibles encontradas:\n', gen);
                soluciones_unicas = unique(soluciones_factibles, 'rows');
                
                % Limitar a un máximo soluciones
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    sol = soluciones_unicas(s,:);
                    fprintf('x: %.4f, y: %.4f - Himmelblau: %.4f\n', ...
                            sol(1), sol(2), himmelblau(sol));
                end
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones más\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generación %d: Mejor aptitud = %.4f\n', gen, mejor_apt);
            end
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
