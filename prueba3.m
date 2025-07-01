%% Algoritmo genético para las 8 reinas
clear all; close all; clc;

% Parámetros
tam_poblacion = 100; max_generaciones = 200; prob_cruce = 0.85; prob_mutacion = 0.05; N = 8;

% Población inicial (permutaciones)
poblacion = zeros(tam_poblacion, N);
for i = 1:tam_poblacion
    poblacion(i, :) = randperm(N);
end

% Guardar conflictos para gráfica
min_conflictos = zeros(max_generaciones, 1);

% Algoritmo genético
mejor_ind = zeros(1, N); mejor_fit = 0;
for g = 1:max_generaciones
    % Evaluar aptitud
    fitness = zeros(tam_poblacion, 1);
    for i = 1:tam_poblacion
        conflicto = 0;
        for j = 1:N-1
            for k = j+1:N
                if poblacion(i,j) == poblacion(i,k) || abs(j-k) == abs(poblacion(i,j)-poblacion(i,k))
                    conflicto = conflicto + 1;
                end
            end
        end
        fitness(i) = 1 / (1 + conflicto);
    end
    
    % Guardar mejor individuo (elitismo)
    [fit, idx] = max(fitness);
    if fit > mejor_fit
        mejor_fit = fit;
        mejor_ind = poblacion(idx, :);
    end
    min_conflictos(g) = (1/fit) - 1;
    
    % Mostrar progreso
    fprintf('Gen %d: Fitness = %.4f, Conflictos = %.0f\n', g, fit, min_conflictos(g));
    
    % Parar si no hay conflictos
    if fit == 1, break; end
    
    % Selección por torneo
    p = zeros(tam_poblacion, N);
    for i = 1:tam_poblacion
        c = randperm(tam_poblacion, 2);
        if fitness(c(1)) >= fitness(c(2)) % Seleccionar el mejor
            p(i, :) = poblacion(c(1), :);
        else
            p(i, :) = poblacion(c(2), :);
        end
    end
    
    % Cruce y mutación
    descendencia = zeros(tam_poblacion, N);
    for i = 1:2:tam_poblacion
        h1 = p(i, :); h2 = p(i+1, :);
        if rand < prob_cruce
            pts = sort(randperm(N-1, 2));
            h1(pts(1)+1:pts(2)) = p(i, pts(1)+1:pts(2));
            h2(pts(1)+1:pts(2)) = p(i+1, pts(1)+1:pts(2));
            h1([pts(2)+1:N 1:pts(1)]) = setdiff(p(i+1, :), h1(pts(1)+1:pts(2)), 'stable');
            h2([pts(2)+1:N 1:pts(1)]) = setdiff(p(i, :), h2(pts(1)+1:pts(2)), 'stable');
        end
        if rand < prob_mutacion
            idx = randperm(N, 2);
            h1([idx(1) idx(2)]) = h1([idx(2) idx(1)]);
        end
        if rand < prob_mutacion
            idx = randperm(N, 2);
            h2([idx(1) idx(2)]) = h2([idx(2) idx(1)]);
        end
        descendencia(i:i+1, :) = [h1; h2];
    end
    
    % Elitismo
    descendencia(1, :) = mejor_ind;
    poblacion = descendencia;
end

% Resultado final
fprintf('\nResultado:\nFitness: %.4f, Conflictos: %.0f\nIndividuo: %s\n', mejor_fit, (1/mejor_fit)-1, mat2str(mejor_ind));

% Gráfica de reinas
figure;
plot(mejor_ind, 1:N, 'rx', 'MarkerSize', 10);
xlim([0.5 8.5]); ylim([0.5 8.5]); grid on;
title('8 Reinas'); xlabel('Columnas'); ylabel('Filas');

% Gráfica de convergencia
figure;
plot(1:max_generaciones, min_conflictos, 'LineWidth', 2);
xlabel('Gen'); ylabel('Conflictos'); title('Evolución 8 reinas'); grid on;


%% ============================================================
%% ============================================================
%% ============================================================

% Algoritmo genético para ajustar parábola y = ax² + bx + c
clear all; close all; clc;

% Datos (parábola: y = 2x² + 3x + 5)
x = linspace(-10, 10, 100);
y = 2*x.^2 + 3*x + 5;

% Parámetros
params.n = 100; % Tamaño población
params.gens = 200; % Generaciones
params.pc = 0.8; % Prob. cruce
params.pm = 0.1; % Prob. mutación

% Configuración
prob.tam = 3; % [a, b, c]
prob.lim_inf = [-10, -10, -10];
prob.lim_sup = [10, 10, 10];
prob.sigma = 1;
prob.x = x;
prob.y = y;

[mejor, fitness, errores] = algoritmo_genetico(prob, params);

% Resultados
fprintf('Parábola: y = %.4fx² + %.4fx + %.4f\n', mejor);
fprintf('Error: %.6f\n', 1/fitness - 1);

% Gráfica ajuste
figure;
plot(x, y, 'rx');
hold on;
x_fit = linspace(min(x)-1, max(x)+1, 100);
y_fit = mejor(1)*x_fit.^2 + mejor(2)*x_fit + mejor(3);
plot(x_fit, y_fit, 'k-');
legend('Datos', 'Ajuste');
title('Ajuste de Parábola');
xlabel('x'); ylabel('y'); grid on;

% Gráfica convergencia
figure;
plot(1:params.gens, errores);
xlabel('Generación'); ylabel('Error'); title('Convergencia'); grid on;

%% Algoritmo Genético
function [mejor_i, mejor_f, errores] = algoritmo_genetico(prob, params)
    poblacion = prob.lim_inf + rand(params.n, prob.tam) .* (prob.lim_sup - prob.lim_inf);
    errores = zeros(params.gens, 1);
    
    for gen = 1:params.gens
        aptitud = zeros(params.n, 1);
        for i = 1:params.n
            aptitud(i) = fitness_parabola(poblacion(i,:), prob);
        end
        [mejor_f, idx] = max(aptitud);
        elite = poblacion(idx,:);
        errores(gen) = 1/mejor_f - 1;
        
        if mod(gen, 10) == 0 || gen == 1
            fprintf('Gen %d: Error = %.4f\n', gen, errores(gen));
        end
        
        pp = seleccionPorRuleta(poblacion, aptitud);
        descendencia = zeros(params.n, prob.tam);
        for i = 2:2:params.n-1
            h1 = pp(i,:); h2 = pp(i+1,:);
            if rand < params.pc
                [h1, h2] = UnPunto(h1, h2);
            end
            if rand < params.pm
                h1 = mutacionGaussiana(h1, prob);
            end
            if rand < params.pm
                h2 = mutacionGaussiana(h2, prob);
            end
            descendencia(i:i+1,:) = [h1; h2];
        end
        descendencia(1,:) = elite;
        poblacion = descendencia;
    end
    aptitud = zeros(params.n, 1);
    for i = 1:params.n
        aptitud(i) = fitness_parabola(poblacion(i,:), prob);
    end
    [mejor_f, idx] = max(aptitud);
    mejor_i = poblacion(idx,:);
end

%% Selección por Ruleta
function p = seleccionPorRuleta(poblacion, fitness)
    probs = fitness / sum(fitness);
    p = poblacion(randsample(1:size(poblacion,1), size(poblacion,1), true, probs), :);
end

%% Cruce Aritmético
function [h1, h2] = UnPunto(p1, p2)
    tamCrom = length(p1);
    pCorte = randi([1, tamCrom-1]);
    h1 = [p1(1:pCorte), p2(pCorte+1:end)];
    h2 = [p2(1:pCorte), p1(pCorte+1:end)];
end


%% Mutación Gaussiana
function i_m = mutacionGaussiana(i, prob)
    i_m = i;
    mascara = rand(1, length(i)) < 0.5;
    i_m(mascara) = i(mascara) + prob.sigma * randn(1, sum(mascara));
    i_m = max(min(i_m, prob.lim_sup), prob.lim_inf);
end

%% Fitness
function apt = fitness_parabola(i, prob)
    y_calc = i(1)*prob.x.^2 + i(2)*prob.x + i(3);
    apt = 1 / (1 + mean((y_calc - prob.y).^2));
end