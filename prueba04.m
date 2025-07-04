%% Algoritmo genetico para las 8 reinas
clear all; close all; clc;

% Parametros
tam_p = 100; max_generaciones = 200; prob_cruce = 0.85; prob_mutacion = 0.05; N = 8;

disp("a) Evalulacion de aptitud: 1/(1+conclictos)");
disp("  conflictos: si las reinas se encuentran las mismas filas, columnas o diagonal");
disp("b) Seleccion: torneo");
disp("c) Cruce: orden");
disp("d) Mutacion: intercambio");

% Poblacion inicial (permutaciones)
poblacion = zeros(tam_p, N);
for i = 1:tam_p
    poblacion(i, :) = randperm(N);
end
disp("Poblacion inicial");
disp(poblacion);

% Guardar conflictos para grafica
min_conflictos = zeros(max_generaciones, 1);

% Algoritmo genetico
mejor_ind = zeros(1, N); mejor_fit = 0;
for g = 1:max_generaciones
    % Evaluar aptitud
    fitness = zeros(tam_p, 1);
    soluciones_generacion = [];
    
    for i = 1:tam_p
        conflicto = 0;
        for j = 1:N-1
            for k = j+1:N
                if poblacion(i,j) == poblacion(i,k) || abs(j-k) == abs(poblacion(i,j)-poblacion(i,k))
                    conflicto = conflicto + 1;
                end
            end
        end
        fitness(i) = 1 / (1 + conflicto);
        
        if conflicto == 0
            soluciones_generacion = [soluciones_generacion; poblacion(i,:)];
        end
    end
    % Elitismo
    [fit, idx] = max(fitness);
    if fit > mejor_fit
        mejor_fit = fit;
        mejor_ind = poblacion(idx, :);
    end
    min_conflictos(g) = (1/fit) - 1;

    % Mostrar las soluciones encontradas
    if mod(g, 10) == 0
        if ~isempty(soluciones_generacion)
        fprintf('\nGen %d: Soluciones encontradas:\n', g);
        soluciones_unicas = unique(soluciones_generacion, 'rows', 'stable');
        for i=1:size(soluciones_unicas,1)
            s=soluciones_unicas(i,:);
            fprintf('Solucion %d:\n', i);
            disp(s);
            fprintf('Fitness = %.4f, es una solucion optima poruqe tiene %.0f conflictos\n', ...
                aptitud_f(s,N), (1/aptitud_f(s,N)) -1);
        end
        end
    end
    
    % Seleccion por torneo
    p = zeros(tam_p, N);
    for i = 1:tam_p
        c = randperm(tam_p, 2);
        if fitness(c(1)) >= fitness(c(2))
            p(i, :) = poblacion(c(1), :);
        else
            p(i, :) = poblacion(c(2), :);
        end
    end

    descendencia = zeros(tam_p, N);
    for i = 1:2:tam_p
        h1 = p(i, :); h2 = p(i+1, :);
        % Cruce
        if rand < prob_cruce
            pts = sort(randperm(N-1, 2));
            h1(pts(1)+1:pts(2)) = p(i, pts(1)+1:pts(2));
            h2(pts(1)+1:pts(2)) = p(i+1, pts(1)+1:pts(2));
            h1([pts(2)+1:N 1:pts(1)]) = setdiff(p(i+1, :), h1(pts(1)+1:pts(2)), 'stable');
            h2([pts(2)+1:N 1:pts(1)]) = setdiff(p(i, :), h2(pts(1)+1:pts(2)), 'stable');
        end
        % MUtacion
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

% Grafica de reinas
% Resultado final
fprintf('\n--- RESULTADO FINAL ---\n');
fprintf('Generaciones: %d\n', g);
fprintf('Conflictos: %d\n', (1/mejor_fit)-1);
fprintf('Posiciones: %s\n', mat2str(mejor_ind));

% Grafica del tablero de ajedrez con reinas
figure;
% Crear tablero de ajedrez
board = zeros(N);
for i = 1:N
    for j = 1:N
        board(i,j) = mod(i+j,2);
    end
end
imagesc(board);
colormap([1 1 1; 0.7 0.7 0.7]);
hold on;

% Dibujar reinas
for col = 1:N
    row = mejor_ind(col);
    text(col, row, 'x', 'Color', 'r', 'FontSize', 15, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

title('Solucion para las 8 Reinas');
set(gca, 'XTick', 1:N, 'YTick', 1:N, 'XAxisLocation', 'top');
axis square;
grid on;

% Grafica de convergencia
figure;
plot(1:max_generaciones, min_conflictos, 'LineWidth', 2);
xlabel('Gen'); ylabel('Conflictos'); title('Evolucion 8 reinas'); grid on;


function fitness = aptitud_f(individuo,N)
conflicto = 0;
for j = 1:N-1
    for k = j+1:N
        if individuo(j) == individuo(k) || abs(j-k) == abs(individuo(j)-individuo(k))
            conflicto = conflicto + 1;
        end
    end
end
fitness = 1 / (1 + conflicto);
end