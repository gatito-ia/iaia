%% Algoritmo genetico para ajustar parabola y = ax^2 + bx + c
clear all; close all; clc;

disp("a) Evalulacion de aptitud: ECM");
disp("b) Seleccion: ruleta");
disp("c) Cruce: Un punto");
disp("d) Mutacion: gausiana");

% Datos (parabola: y = 2x^2 + 3x + 5)
x = linspace(-10, 10, 100);
y = 2*x.^2 + 3*x + 5;

% Parametros
params.n = 100; % Tamano poblacion
params.gens = 200; % Generaciones
params.pc = 0.8; % Prob. cruce
params.pm = 0.1; % Prob. mutación

% Configuracion
prob.tam = 3; % [a, b, c]
prob.lim_inf = [-10, -10, -10];
prob.lim_sup = [10, 10, 10];
prob.sigma = 1;
prob.x = x;
prob.y = y;

[mejor, fitness, errores] = algoritmo_genetico(prob, params);

% Resultado final
fprintf('\n--- RESULTADO FINAL ---\n');
fprintf('Parabola: y = %.4fx^2 + %.4fx + %.4f\n', mejor);
fprintf('Error: %.6f\n', 1/fitness - 1);

% Grafica
figure;
plot(x, y, 'rx');
hold on;
x_fit = linspace(min(x)-1, max(x)+1, 100);
y_fit = mejor(1)*x_fit.^2 + mejor(2)*x_fit + mejor(3);
plot(x_fit, y_fit, 'k-');
legend('Datos', 'Ajuste');
title('Ajuste de Parábola');
xlabel('x'); ylabel('y'); grid on;

% Grafica convergencia
figure;
plot(1:params.gens, errores);
xlabel('Generacion'); ylabel('Error'); title('Convergencia'); grid on;

%% Algoritmo Genetico
function [mejor_i, mejor_f, errores] = algoritmo_genetico(prob, params)
    poblacion = prob.lim_inf + rand(params.n, prob.tam) .* (prob.lim_sup - prob.lim_inf);
    errores = zeros(params.gens, 1);
    disp("Poblacion inicial");
    disp(poblacion);
    umbral_error = 0.08;  % Umbral para considerar una solucion valida
    
    for gen = 1:params.gens
        aptitud = zeros(params.n, 1);
        soluciones_generacion = [];
        
        for i = 1:params.n
            aptitud(i) = fitness_parabola(poblacion(i,:), prob);
            error_actual = 1/aptitud(i) - 1;
           
            if error_actual < umbral_error
                soluciones_generacion = [soluciones_generacion; poblacion(i,:) error_actual];
            end
        end
        
        [mejor_f, idx] = max(aptitud);
        elite = poblacion(idx,:);
        errores(gen) = 1/mejor_f - 1;
        
        if mod(gen, 10) == 0
            if ~isempty(soluciones_generacion)
                fprintf('Gen %d: Soluciones encontradas (error < %.2f):\n', gen, umbral_error);
               
                [~, idx_sort] = sort(soluciones_generacion(:,end));
                soluciones_ordenadas = soluciones_generacion(idx_sort, 1:end-1);
                disp(unique(soluciones_ordenadas, 'rows', 'stable')); % Mostrar unicas
            else
                fprintf('Gen %d: Error = %.4f\n', gen, errores(gen));
            end
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
    
    % Evaluacion final
    aptitud = zeros(params.n, 1);
    for i = 1:params.n
        aptitud(i) = fitness_parabola(poblacion(i,:), prob);
    end
    [mejor_f, idx] = max(aptitud);
    mejor_i = poblacion(idx,:);
end

%% Seleccion por Ruleta
function p = seleccionPorRuleta(poblacion, fitness)
    probs = fitness / sum(fitness);
    p = poblacion(randsample(1:size(poblacion,1), size(poblacion,1), true, probs), :);
end

%% Cruce un punto
function [h1, h2] = UnPunto(p1, p2)
    tamCrom = length(p1);
    pCorte = randi([1, tamCrom-1]);
    h1 = [p1(1:pCorte), p2(pCorte+1:end)];
    h2 = [p2(1:pCorte), p1(pCorte+1:end)];
end


%% Mutacion Gaussiana
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