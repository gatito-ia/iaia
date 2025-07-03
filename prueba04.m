%% Algoritmo genético para las 8 reinas
clear all; close all; clc;

% Parámetros
tam_poblacion = 100; max_generaciones = 200; prob_cruce = 0.85; prob_mutacion = 0.05; N = 8;

disp("a) Evalulación de aptitud: 1/(1+conclictos)");
disp("  conflictos: si las reinas se encuentran las mismas filas, columnas o diagonal");
disp("b) Selección: torneo");
disp("c) Cruce: orden");
disp("d) Mutación: intercambio");

% Población inicial (permutaciones)
poblacion = zeros(tam_poblacion, N);
for i = 1:tam_poblacion
    poblacion(i, :) = randperm(N);
end
disp("Población inicial");
disp(poblacion);

% Guardar conflictos para gráfica
min_conflictos = zeros(max_generaciones, 1);

% Algoritmo genético
mejor_ind = zeros(1, N); mejor_fit = 0;
for g = 1:max_generaciones
        % Evaluar aptitud
    fitness = zeros(tam_poblacion, 1);
    soluciones_generacion = [];  % Almacenar soluciones de esta generación
    
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
        
        % Si es solución perfecta (0 conflictos), guardarla
        if conflicto == 0
            soluciones_generacion = [soluciones_generacion; poblacion(i,:)];
        end
    end
    % Guardar mejor individuo (elitismo)
    [fit, idx] = max(fitness);
    if fit > mejor_fit
        mejor_fit = fit;
        mejor_ind = poblacion(idx, :);
    end
    min_conflictos(g) = (1/fit) - 1;

    % Mostrar todas las soluciones encontradas en esta generación
    if mod(g, 10) == 0
        if ~isempty(soluciones_generacion)
        fprintf('Gen %d: Soluciones encontradas:\n', g);
        disp(unique(soluciones_generacion, 'rows')); % Mostrar solo soluciones únicas
        end
        fprintf('Gen %d: Fitness = %.4f, Conflictos = %.0f\n', g, fit, min_conflictos(g));
    end
    
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
% Resultado final
fprintf('\n--- RESULTADO FINAL ---\n');
fprintf('Generaciones: %d\n', g);
fprintf('Conflictos: %d\n', (1/mejor_fit)-1);
fprintf('Posiciones: %s\n', mat2str(mejor_ind));

% Gráfica del tablero de ajedrez con reinas (versión simple)
figure;
% Crear tablero de ajedrez
board = zeros(N);
for i = 1:N
    for j = 1:N
        board(i,j) = mod(i+j,2);
    end
end
imagesc(board);
colormap([1 1 1; 0.7 0.7 0.7]); % Blanco y gris
hold on;

% Dibujar reinas
for col = 1:N
    row = mejor_ind(col);
    text(col, row, 'x', 'Color', 'r', 'FontSize', 15, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

title('Solución para las 8 Reinas');
set(gca, 'XTick', 1:N, 'YTick', 1:N, 'XAxisLocation', 'top');
axis square;
grid on;

% Gráfica de convergencia
figure;
plot(1:max_generaciones, min_conflictos, 'LineWidth', 2);
xlabel('Gen'); ylabel('Conflictos'); title('Evolución 8 reinas'); grid on;
