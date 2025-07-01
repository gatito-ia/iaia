%% edificio - Wifi

clear; clc; close all;

%% 1. Definición del problema
filas = 20;
columnas = 30;
num_reps = 4;
% [y1, y2, y3, y4, x1, x2, x3, x4]
alcances = [5, 3, 7, 4];

%% 2. Parámetros del Algoritmo Genético
pob_size = 50;      % Tamaño de la población
max_gens = 100;     % Número máximo de generaciones
pc = 0.8;           % Probabilidad de cruce
pm = 0.1;           % Probabilidad de mutación

%% 3. Inicialización de población
% Cada individuo tiene las coordenadas (x,y) de cada repetidor
poblacion = [randi([1 filas], pob_size, num_reps), randi([1 columnas], pob_size, num_reps)];

%% 4. Evolución
mejor_aptitud = zeros(max_gens, 1); % Para guardar el progreso

for gen = 1:max_gens
    % Evaluar población
    aptitudes = zeros(pob_size, 1);
    for i = 1:pob_size
        aptitudes(i) = calcular_cobertura(poblacion(i,:), filas, columnas, alcances, num_reps);
    end
    
    % Guardar el mejor de esta generación
    [mejor_aptitud(gen), idx] = max(aptitudes);
    mejor_individuo = poblacion(idx,:);
    
    % Mostrar progreso cada 10 generaciones
    if mod(gen,10) == 0
        fprintf('Gen %d: Mejor cobertura = %d celdas\n', gen, mejor_aptitud(gen));
    end
    
    % SELECCIÓN
    padres = zeros(size(poblacion));
    for i = 1:pob_size
        k = 3;
        competidores = randperm(pob_size, k);
        [~, idx_mejor] = max(aptitudes(competidores));
        padres(i,:) = poblacion(competidores(idx_mejor),:);
    end
    
    % CRUCE
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
    
    hijos(1,:) = mejor_individuo;
    poblacion = hijos;
end

%% 5. Resultados finales
[mejor_cobertura, idx] = max(aptitudes);
mejor_solucion = poblacion(idx,:);

fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor cobertura encontrada: %d celdas (de %d)\n', mejor_cobertura, filas*columnas);
disp('Coordenadas (fila,columna) y alcances de los repetidores:');
reps_info = [mejor_solucion(1:num_reps)' mejor_solucion(num_reps+1:end)' alcances'];
disp(reps_info);

%% 6. Visualización simple
figure;
plot(1:max_gens, mejor_aptitud, 'b-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Celdas cubiertas');
title('Evolución de la cobertura WiFi');
grid on;

%% 7. Visualización de repetidores
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
end
hold off;

%% Función para calcular cobertura con alcances diferentes
function cobertura = calcular_cobertura(individuo, filas, columnas, alcances, num_reps)
    reps_y = individuo(1:num_reps);     % Filas
    reps_x = individuo(num_reps+1:end); % Columnas
    
    % Matriz para marcar celdas cubiertas
    mapa = zeros(filas, columnas);
    
    for r = 1:num_reps
        x = reps_x(r);
        y = reps_y(r);
        alcance = alcances(r);
        
        % Verificar que esté dentro del edificio
        if x < 1 || x > columnas || y < 1 || y > filas
            continue; % Si está fuera, no cubre nada
        end
        
        % Marcar área circular alrededor del repetidor
        for i = max(1,x-alcance):min(columnas,x+alcance)
            for j = max(1,y-alcance):min(filas,y+alcance)
                if (i-x)^2 + (j-y)^2 <= alcance^2
                    mapa(j,i) = 1; % Nota: mapa(fila,columna)
                end
            end
        end
    end
    cobertura = sum(mapa(:));
end

%% ================================================
%% ================================================
%% ================================================