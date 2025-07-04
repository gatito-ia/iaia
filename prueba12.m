%% Edificio - Wifi
clear; clc; close all;

%% Definicion del problema
filas = 20;
columnas = 30;
num_reps = 4;
% [y1, y2, y3, y4, x1, x2, x3, x4]
alcances = [5, 3, 7, 4];

%% Parámetros del Algoritmo Genetico
pob_size = 50;      % Tamano de la poblacion
max_gens = 100;     % Numero maximo de generaciones
pc = 0.8;           % Probabilidad de cruce
pm = 0.1;           % Probabilidad de mutacion

disp("a) Aptitud: Maxima cobertura");
disp("b) Seleccion: Torneo");
disp("c) Cruce: Un punto");
disp("d) Mutacion: Aleatoria por reemplazo en coordenadas");

%% Inicialización de poblacion
poblacion = [randi([1 filas], pob_size, num_reps), randi([1 columnas], pob_size, num_reps)];
disp("Poblacion inicial"); disp(poblacion);

%% Evolucion
mejor_aptitud = zeros(max_gens, 1);

for gen = 1:max_gens
    % Evaluar poblacion y almacenar soluciones optimas
    aptitudes = zeros(pob_size, 1);
    soluciones_optimas = [];  % Almacenar soluciones optimas    
    for i = 1:pob_size
        aptitudes(i) = calcular_cobertura(poblacion(i,:), filas, columnas, alcances, num_reps);
    end
    % Guardar el mejor y encontrar soluciones optimas
    [mejor_aptitud(gen), idx] = max(aptitudes);
    mejor_individuo = poblacion(idx,:);
    max_cobertura = mejor_aptitud(gen);
    
    % Almacenar todas las soluciones con maxima cobertura
    for i = 1:pob_size
        if aptitudes(i) == max_cobertura
            soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
        end
    end
    
    if mod(gen, 20) == 0
        % Mostrar soluciones optimas
        if ~isempty(soluciones_optimas)
            fprintf('\nGeneracion %d: Soluciones optimas encontradas (Cobertura = %d celdas):\n', gen, max_cobertura);
            soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
            
            num_mostrar = min(6, size(soluciones_unicas, 1));
            for s = 1:num_mostrar
                fprintf('Solucion %d:\n', s);
                % Extraer coordenadas x e y para los 4 amplificadores
                coords = reshape(soluciones_unicas(s,:), num_reps, 2);
                % Mostrar en formato más legible
                for amp = 1:num_reps
                    fprintf('Amplificador %d: (x=%d, y=%d)\n', amp, coords(amp,2), coords(amp,1));
                end
                fprintf('\n');  % Separador entre soluciones
            end
            
            if size(soluciones_unicas, 1) > 6
                fprintf('... y %d soluciones mas\n', size(soluciones_unicas, 1)-6);
            end
        else
            fprintf('Gen %d: Mejor cobertura = %d celdas\n', gen, mejor_aptitud(gen));
        end
    end
    
    % Selecion (torneo)
    padres = zeros(size(poblacion));
    for i = 1:pob_size
        k = 3;
        competidores = randperm(pob_size, k);
        [~, idx_mejor] = max(aptitudes(competidores));
        padres(i,:) = poblacion(competidores(idx_mejor),:);
    end
    
    % cruce (un punto)
    hijos = padres;
    for i = 1:2:pob_size-1
        if rand < pc
            punto = randi([1 num_reps*2-1]);
            temp = hijos(i, punto+1:end);
            hijos(i, punto+1:end) = hijos(i+1, punto+1:end);
            hijos(i+1, punto+1:end) = temp;
        end
    end
    
    % mutacion
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
fprintf('Coordenadas y alcances de los repetidores:\n\n');

% Extraer coordenadas y alcances
coords_y = mejor_solucion(1:num_reps);    % Filas (coordenada y)
coords_x = mejor_solucion(num_reps+1:end); % Columnas (coordenada x)

% Mostrar cada repetidor con su información
for rep = 1:num_reps
    fprintf('Repetidor %d:\n', rep);
    fprintf('  Posición: (x=%d, y=%d)\n', coords_x(rep), coords_y(rep));
    fprintf('  Alcance: %d celdas\n\n', alcances(rep));
end
%% Visualización
figure;
plot(1:max_gens, mejor_aptitud, 'b-', 'LineWidth', 1.5);
xlabel('Generacion');
ylabel('Celdas cubiertas');
title('Evolucion de la cobertura WiFi');
grid on;

%% Visualización de repetidores
figure;
hold on;
axis equal;
xlim([0 columnas+1]);
ylim([0 filas+1]);
title('Ubicación optima de repetidores');
xlabel('Columnas');
ylabel('Filas');
grid on;
% Dibujar areas de cobertura
reps_x = mejor_solucion(num_reps+1:end); % Coordenadas x (columnas)
reps_y = mejor_solucion(1:num_reps);     % Coordenadas y (filas)
for i = 1:num_reps
    x = reps_x(i);
    y = reps_y(i);
    alcance = alcances(i);
    % Dibujar area de cobertura
    rectangle('Position', [x-alcance, y-alcance, 2*alcance, 2*alcance], ...
              'Curvature', [1 1], 'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
    text(x-alcance/2,y, sprintf('Repetidor %d', i));
end
hold off;

%% Funcion de fitness
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