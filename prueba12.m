%% Edificio - Wifi
clear; clc; close all;

%% Definición del problema
filas = 20;
columnas = 30;
num_reps = 4;
% [y1, y2, y3, y4, x1, x2, x3, x4]
alcances = [5, 3, 7, 4];

%% Parámetros del Algoritmo Genético
pob_size = 50;      % Tamaño de la población
max_gens = 100;     % Número máximo de generaciones
pc = 0.8;           % Probabilidad de cruce
pm = 0.1;           % Probabilidad de mutación

disp("a) Aptitud: Maxima cobertura");
disp("b) Selección: Torneo");
disp("c) Cruce: Un punto");
disp("d) Mutación: Aleatoria por reemplazo en coordenadas");

%% Inicialización de población
poblacion = [randi([1 filas], pob_size, num_reps), randi([1 columnas], pob_size, num_reps)];
disp("Población inicial"); disp(poblacion);

%% Evolución
mejor_aptitud = zeros(max_gens, 1);

for gen = 1:max_gens
    % Evaluar población y almacenar soluciones óptimas
    aptitudes = zeros(pob_size, 1);
    soluciones_optimas = [];  % Almacenar soluciones óptimas
    
    for i = 1:pob_size
        aptitudes(i) = calcular_cobertura(poblacion(i,:), filas, columnas, alcances, num_reps);
    end
    
    % Guardar el mejor y encontrar soluciones óptimas
    [mejor_aptitud(gen), idx] = max(aptitudes);
    mejor_individuo = poblacion(idx,:);
    max_cobertura = mejor_aptitud(gen);
    
    % Almacenar todas las soluciones con máxima cobertura
    for i = 1:pob_size
        if aptitudes(i) == max_cobertura
            soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
        end
    end
    
    if mod(gen, 20) == 0
        % Mostrar soluciones óptimas
        if ~isempty(soluciones_optimas)
            fprintf('\nGeneración %d: Soluciones óptimas encontradas (Cobertura = %d celdas):\n', gen, max_cobertura);
            soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
            
            num_mostrar = min(6, size(soluciones_unicas, 1));
            for s = 1:num_mostrar
                fprintf('Solución %d:\n', s);
                disp(reshape(soluciones_unicas(s,:), num_reps, 2)'); % Mostrar en formato coordenadas
            end
            
            if size(soluciones_unicas, 1) > 6
                fprintf('... y %d soluciones más\n', size(soluciones_unicas, 1)-6);
            end
        else
            fprintf('Gen %d: Mejor cobertura = %d celdas\n', gen, mejor_aptitud(gen));
        end
    end
    
    % SELECCIÓN (torneo)
    padres = zeros(size(poblacion));
    for i = 1:pob_size
        k = 3;
        competidores = randperm(pob_size, k);
        [~, idx_mejor] = max(aptitudes(competidores));
        padres(i,:) = poblacion(competidores(idx_mejor),:);
    end
    
    % CRUCE (un punto)
    hijos = padres;
    for i = 1:2:pob_size-1
        if rand < pc
            punto = randi([1 num_reps*2-1]);
            temp = hijos(i, punto+1:end);
            hijos(i, punto+1:end) = hijos(i+1, punto+1:end);
            hijos(i+1, punto+1:end) = temp;
        end
    end
    
    % MUTACIÓN
    for i = 1:pob_size
        if rand < pm
            if rand < 0.5
                gen_mutar = randi([1 num_reps]);
                hijos(i, gen_mutar) = randi([1 filas]);
            else
                gen_mutar = randi([num_reps+1 num_reps*2]);
                hijos(i, gen_mutar) = randi([1 columnas]);
            end
        end
    end
    
    % Elitismo
    hijos(1,:) = mejor_individuo;
    poblacion = hijos;
end

%% Resultados finales
[mejor_cobertura, idx] = max(aptitudes);
mejor_solucion = poblacion(idx,:);

fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor cobertura encontrada: %d celdas (de %d)\n', mejor_cobertura, filas*columnas);
disp('Coordenadas (fila,columna) y alcances de los repetidores:');
reps_info = [mejor_solucion(1:num_reps)' mejor_solucion(num_reps+1:end)' alcances'];
disp(reps_info);

%% Visualización
figure;
plot(1:max_gens, mejor_aptitud, 'b-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Celdas cubiertas');
title('Evolución de la cobertura WiFi');
grid on;

%% Visualización de repetidores
figure;
hold on;
axis equal;
xlim([0 columnas+1]);
ylim([0 filas+1]);
title('Ubicación óptima de repetidores');
xlabel('Columnas');
ylabel('Filas');
grid on;
% Dibujar áreas de cobertura
reps_x = mejor_solucion(num_reps+1:end); % Coordenadas x (columnas)
reps_y = mejor_solucion(1:num_reps);     % Coordenadas y (filas)
for i = 1:num_reps
    x = reps_x(i);
    y = reps_y(i);
    alcance = alcances(i);
    % Dibujar área de cobertura
    rectangle('Position', [x-alcance, y-alcance, 2*alcance, 2*alcance], ...
              'Curvature', [1 1], 'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
    text(x-alcance/2,y, sprintf('Repetidor %d', i));
end
hold off;

%% Función de fitness
function cobertura = calcular_cobertura(individuo, filas, columnas, alcances, num_reps)
    reps_y = individuo(1:num_reps);
    reps_x = individuo(num_reps+1:end);
    
    mapa = zeros(filas, columnas);
    
    for r = 1:num_reps
        x = reps_x(r);
        y = reps_y(r);
        alcance = alcances(r);
        
        if x < 1 || x > columnas || y < 1 || y > filas
            continue;
        end
        
        for i = max(1,x-alcance):min(columnas,x+alcance)
            for j = max(1,y-alcance):min(filas,y+alcance)
                if (i-x)^2 + (j-y)^2 <= alcance^2
                    mapa(j,i) = 1;
                end
            end
        end
    end
    cobertura = sum(mapa(:));
end