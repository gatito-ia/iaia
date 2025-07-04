%% Programación de Escaneos de Puertos
clear; clc; close all;

%% Definición del problema
num_puertos = 100;           % Número total de puertos a escanear
intervalo_tiempo = 24;       % Ventana de tiempo en horas
max_escaneos = 10;           % Máximo número de escaneos permitidos
tiempo_escaneo = 0.1;        % Duración de cada escaneo en horas
prob_deteccion = @(t) 0.3 + 0.6*(t/intervalo_tiempo); % Probabilidad de detección aumenta con el tiempo

%% Parámetros del Algoritmo Genético
pob_size = 50;               % Tamaño de la población
max_gens = 100;              % Número máximo de generaciones
pc = 0.8;                    % Probabilidad de cruce
pm = 0.1;                    % Probabilidad de mutación

disp("a) Aptitud: Minimizar probabilidad de detección");
disp("b) Selección: Torneo");
disp("c) Cruce: Dos puntos");
disp("d) Mutación: Aleatoria por desplazamiento");

%% Inicialización de población
% Cada individuo representa tiempos de inicio para los escaneos
poblacion = rand(pob_size, max_escaneos) * intervalo_tiempo;
disp("Población inicial generada");

%% Evolución
mejor_aptitud = zeros(max_gens, 1);
mejor_num_escaneos = zeros(max_gens, 1);

for gen = 1:max_gens
    % Evaluar población
    aptitudes = zeros(pob_size, 1);
    num_escaneos = zeros(pob_size, 1);
    
    for i = 1:pob_size
        [aptitudes(i), num_escaneos(i)] = evaluar_individuo(poblacion(i,:), ...
                                    intervalo_tiempo, tiempo_escaneo, prob_deteccion);
    end
    
    % Guardar el mejor
    [mejor_aptitud(gen), idx] = min(aptitudes);
    mejor_individuo = poblacion(idx,:);
    mejor_num_escaneos(gen) = num_escaneos(idx);
    
    % Mostrar progreso cada 20 generaciones
    if mod(gen, 20) == 0
        fprintf('Gen %d: Mejor aptitud=%.4f, Escaneos usados=%d\n', ...
                gen, mejor_aptitud(gen), mejor_num_escaneos(gen));
    end
    
    % Selección (torneo)
    padres = zeros(size(poblacion));
    for i = 1:pob_size
        k = 3; % Tamaño del torneo
        competidores = randperm(pob_size, k);
        [~, idx_mejor] = min(aptitudes(competidores));
        padres(i,:) = poblacion(competidores(idx_mejor),:);
    end
    
    % Cruce (dos puntos)
    hijos = padres;
    for i = 1:2:pob_size-1
        if rand < pc
            puntos = sort(randperm(max_escaneos-1, 2));
            temp = hijos(i, puntos(1)+1:puntos(2));
            hijos(i, puntos(1)+1:puntos(2)) = hijos(i+1, puntos(1)+1:puntos(2));
            hijos(i+1, puntos(1)+1:puntos(2)) = temp;
        end
    end
    
    % Mutación (desplazamiento aleatorio)
    for i = 1:pob_size
        if rand < pm
            % Seleccionar un escaneo aleatorio para mutar
            escaneo_mutar = randi([1 max_escaneos]);
            
            % Desplazar el tiempo aleatoriamente (hasta ±2 horas)
            desplazamiento = (rand()*4 - 2);
            hijos(i, escaneo_mutar) = max(0, min(intervalo_tiempo, ...
                                        hijos(i, escaneo_mutar) + desplazamiento));
        end
    end
    
    % Elitismo (mantener el mejor individuo)
    hijos(1,:) = mejor_individuo;
    poblacion = hijos;
end

%% Resultados finales
[mejor_apt, idx] = min(aptitudes);
mejor_cronograma = sort(poblacion(idx,:));
mejor_cronograma(mejor_cronograma==0) = []; % Eliminar tiempos cero no usados

fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor probabilidad de detección: %.4f\n', mejor_apt);
fprintf('Número de escaneos utilizados: %d\n', length(mejor_cronograma));
fprintf('Cronograma óptimo de escaneos (horas):\n');
disp(mejor_cronograma);

%% Visualización de convergencia
figure;
subplot(2,1,1);
plot(1:max_gens, mejor_aptitud, 'b-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Probabilidad de detección');
title('Evolución de la probabilidad de detección');
grid on;

subplot(2,1,2);
plot(1:max_gens, mejor_num_escaneos, 'r-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Número de escaneos');
title('Evolución del número de escaneos utilizados');
grid on;

%% Visualización del cronograma óptimo
figure;
hold on;
title('Cronograma óptimo de escaneos de puertos');
xlabel('Tiempo (horas)');
ylabel('Escaneo');
ylim([0.5 1.5]);

for i = 1:length(mejor_cronograma)
    plot([mejor_cronograma(i) mejor_cronograma(i)+tiempo_escaneo], [1 1], ...
         'b-', 'LineWidth', 10);
    text(mejor_cronograma(i)+tiempo_escaneo/2, 1.1, num2str(i), ...
         'HorizontalAlignment', 'center');
end

plot(0:intervalo_tiempo, prob_deteccion(0:intervalo_tiempo), 'r--', 'LineWidth', 1.5);
legend('Escaneos', 'Probabilidad de detección', 'Location', 'northwest');
xlim([0 intervalo_tiempo]);
hold off;

%% Función de evaluación (fitness)
function [prob_deteccion_total, num_escaneos] = evaluar_individuo(individuo, ...
                                        intervalo_tiempo, tiempo_escaneo, prob_deteccion)
    % Ordenar los tiempos de escaneo y eliminar los no utilizados (tiempo=0)
    tiempos = sort(individuo(individuo > 0));
    num_escaneos = length(tiempos);
    
    if num_escaneos == 0
        prob_deteccion_total = 1; % Máxima penalización si no se hace ningún escaneo
        return;
    end
    
    % Calcular probabilidad de detección para cada escaneo
    prob_detecciones = zeros(1, num_escaneos);
    
    for i = 1:num_escaneos
        t_inicio = tiempos(i);
        t_fin = t_inicio + tiempo_escaneo;
        
        if t_fin > intervalo_tiempo
            prob_detecciones(i) = 1; % Penalización por escanear fuera del intervalo
        else
            % Probabilidad basada en el tiempo promedio del escaneo
            t_promedio = (t_inicio + t_fin)/2;
            prob_detecciones(i) = prob_deteccion(t_promedio);
        end
    end
    
    % Probabilidad total de no ser detectado en ningún escaneo
    prob_no_deteccion = prod(1 - prob_detecciones);
    
    % Probabilidad de ser detectado al menos una vez
    prob_deteccion_total = 1 - prob_no_deteccion;
    
    % Penalización por usar muchos escaneos
    penalizacion = 0.02 * num_escaneos;
    prob_deteccion_total = prob_deteccion_total + penalizacion;
end