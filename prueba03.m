clear all; close all;clc;
% Parametros del algoritmo genetico
parametros = struct();
parametros.tam_poblacion = 100;
parametros.max_generaciones = 500;
parametros.prob_cruce = 0.85;
parametros.prob_mutacion = 0.15;
parametros.tipo_seleccion = 'torneo';
parametros.tipo_cruce = 'orden';
parametros.tipo_mutacion = 'intercambio';

% Configuracion del problema TSP
problema = struct();
problema.nombre = 'Problema del Viajante (TSP)';
problema.tipo = 'permutacion';
problema.tam_cromosoma = 10;
problema.funcion_fitness = @fitness_tsp;

rng(4);
distancias = randi([1 100], problema.tam_cromosoma, problema.tam_cromosoma);
problema.distancias = distancias;

% Ejecutar el algoritmo genetico
[mejor_ind, mejor_fit, historico] = algoritmo_genetico(problema, parametros);

% Mostrar resultados
fprintf('\n=== RESULTADOS TSP ===\n');
fprintf('Mejor ruta encontrada: %s\n', mat2str(mejor_ind));
fprintf('Distancia: %.6f\n', 1/mejor_fit);
fprintf('Fitness: %.6f\n', mejor_fit);

%% ================== FUNCIÓN PRINCIPAL DEL AG ==================

function [mejor_individuo, mejor_fitness, historico_fitness] = algoritmo_genetico(problema, parametros)
%% PARÁMETROS DEL ALGORITMO GENÉTICO
tam_poblacion = parametros.tam_poblacion;
max_generaciones = parametros.max_generaciones;

% Crear población inicial
poblacion = crear_poblacion_inicial(tam_poblacion, problema);

% Historico para análisis
historico_fitness = zeros(max_generaciones, 1);

%% EVOLUCIÓN
for generacion = 1:max_generaciones
% EVALUACION DE APTITUD (FITNESS)
fitness = evaluar_poblacion(poblacion, problema);

% Guardar mejor fitness de la generación
[mejor_fitness_gen, idx_mejor] = max(fitness);
historico_fitness(generacion) = mejor_fitness_gen;
mejor_individuo_gen = poblacion(idx_mejor, :);

% Mostrar progreso
if mod(generacion, 50) == 0 || generacion == 1
fprintf('Generación %d: Mejor fitness = %.4f\n', generacion, mejor_fitness_gen);
end

% SELECCION
padres = seleccionar_padres(poblacion, fitness, parametros);

% CRUCE
descendencia = aplicar_cruce(padres, parametros);

% MUTACION
descendencia = aplicar_mutacion(descendencia, parametros, problema);

% REEMPLAZO
poblacion = descendencia;
end

% Evaluar poblacion final
fitness_final = evaluar_poblacion(poblacion, problema);
[mejor_fitness, idx_mejor] = max(fitness_final);
mejor_individuo = poblacion(idx_mejor, :);

fprintf('\n=== RESULTADO FINAL ===\n');
fprintf('Mejor fitness: %.6f\n', mejor_fitness);
fprintf('Mejor individuo: ');
disp(mejor_individuo);

% Graficar evolución del fitness
figure;
plot(1:max_generaciones, historico_fitness, 'LineWidth', 2);
title('Evolucion del Fitness');
xlabel('Generacion');
ylabel('Mejor Fitness');
grid on;
end

%% ================== FUNCIONES DE EVALUACION ==================

function fitness = evaluar_poblacion(poblacion, problema)
tam_poblacion = size(poblacion, 1);
fitness = zeros(tam_poblacion, 1);

for i = 1:tam_poblacion
fitness(i) = problema.funcion_fitness(poblacion(i, :), problema);
end
end

%% ================== TIPOS DE SELECCION ==================

function padres = seleccionar_padres(poblacion, fitness, parametros)
tam_poblacion = size(poblacion, 1);
num_padres = tam_poblacion;
padres = zeros(num_padres, size(poblacion, 2));

switch lower(parametros.tipo_seleccion)
case 'torneo'
padres = seleccion_torneo(poblacion, fitness, num_padres);
case 'ruleta'
padres = seleccion_ruleta(poblacion, fitness, num_padres);
end
end

function padres = seleccion_torneo(poblacion, fitness, num_padres)
tam_torneo = 3;
tam_poblacion = size(poblacion, 1);
padres = zeros(num_padres, size(poblacion, 2));

for i = 1:num_padres
competidores = randperm(tam_poblacion, tam_torneo);
fitness_competidores = fitness(competidores);

[~, idx_ganador] = max(fitness_competidores);
ganador = competidores(idx_ganador);

padres(i, :) = poblacion(ganador, :);
end
end

function padres = seleccion_ruleta(poblacion, fitness, num_padres)
tam_poblacion = size(poblacion, 1);
padres = zeros(num_padres, size(poblacion, 2));
probabilidades = fitness / sum(fitness);
probabilidades_acumuladas = cumsum(probabilidades);

for i = 1:num_padres
r = rand();
idx = find(probabilidades_acumuladas >= r, 1, 'first');
padres(i, :) = poblacion(idx, :);
end
end

function padres = seleccion_ruleta(poblacion, fitness, num_padres)
tam_poblacion = size(poblacion, 1);
padres = zeros(num_padres, size(poblacion, 2));
probabilidades = fitness / sum(fitness);
idx = randsample(1:tam_poblacion, tam_poblacion, true, probabilidades);
padres = poblacion(idx, :);
end
%% ================== TIPOS DE CRUCE ==================

function descendencia = aplicar_cruce(padres, parametros)
tam_poblacion = size(padres, 1);
tam_cromosoma = size(padres, 2);
descendencia = padres;

if tam_cromosoma == 1
return;
end

for i = 1:2:tam_poblacion-1
if rand() < parametros.prob_cruce
padre1 = padres(i, :);
padre2 = padres(i+1, :);

switch lower(parametros.tipo_cruce)
case 'un_punto'
[hijo1, hijo2] = cruce_un_punto(padre1, padre2);
case 'dos_puntos'
[hijo1, hijo2] = cruce_dos_puntos(padre1, padre2);
case 'uniforme'
[hijo1, hijo2] = cruce_uniforme(padre1, padre2);
case 'orden'
[hijo1, hijo2] = cruce_orden(padre1, padre2);
end
descendencia(i, :) = hijo1;
if i+1 <= tam_poblacion
descendencia(i+1, :) = hijo2;
end
end
end
end

function [hijo1, hijo2] = cruce_un_punto(padre1, padre2)
tam_cromosoma = length(padre1);

if tam_cromosoma == 1
hijo1 = padre1;
hijo2 = padre2;
return;
end

punto_corte = randi([1, tam_cromosoma-1]);

hijo1 = [padre1(1:punto_corte), padre2(punto_corte+1:end)];
hijo2 = [padre2(1:punto_corte), padre1(punto_corte+1:end)];
end

function [hijo1, hijo2] = cruce_dos_puntos(padre1, padre2)
tam_cromosoma = length(padre1);

if tam_cromosoma <= 2
hijo1 = padre1;
hijo2 = padre2;
return;
end

puntos_corte = sort(randperm(tam_cromosoma-1, 2));
p1 = puntos_corte(1);
p2 = puntos_corte(2);

hijo1 = [padre1(1:p1), padre2(p1+1:p2), padre1(p2+1:end)];
hijo2 = [padre2(1:p1), padre1(p1+1:p2), padre2(p2+1:end)];
end

function [hijo1, hijo2] = cruce_uniforme(padre1, padre2)
tam_cromosoma = length(padre1);

if tam_cromosoma == 1
if rand() < 0.5
hijo1 = padre2;
hijo2 = padre1;
else
hijo1 = padre1;
hijo2 = padre2;
end
return;
end

mascara = rand(1, tam_cromosoma) < 0.5;

hijo1 = padre1;
hijo2 = padre2;

hijo1(mascara) = padre2(mascara);
hijo2(mascara) = padre1(mascara);
end

function [hijo1, hijo2] = cruce_orden(padre1, padre2)
tam = length(padre1);

puntos = sort(randperm(tam-1, 2));
p1 = puntos(1);
p2 = puntos(2);

hijo1 = zeros(1, tam);
hijo2 = zeros(1, tam);

hijo1(p1+1:p2) = padre1(p1+1:p2);
hijo2(p1+1:p2) = padre2(p1+1:p2);

% Rellenar el resto para hijo1 con padre2 sin duplicados
idx = p2 + 1;
for i = 1:tam
gen = padre2(mod(p2 + i - 1, tam) + 1);
if ~ismember(gen, hijo1)
if idx > tam
idx = 1;
end
hijo1(idx) = gen;
idx = idx + 1;
end
end

% Rellenar el resto para hijo2 con padre1 sin duplicados
idx = p2 + 1;
for i = 1:tam
gen = padre1(mod(p2 + i - 1, tam) + 1);
if ~ismember(gen, hijo2)
if idx > tam
idx = 1;
end
hijo2(idx) = gen;
idx = idx + 1;
end
end
end


function [hijo1, hijo2] = cruceAritmetico(p1,p2)
    alfa = rand();
    hijos1 = alfa * p1 + (1-alfa) * p2;
    hijos2 = (1-alfa) * p1 + alfa * p2;
end
%% ================== TIPOS DE MUTACION ==================

function descendencia = aplicar_mutacion(descendencia, parametros, problema)
[num_individuos, tam_cromosoma] = size(descendencia);

for i = 1:num_individuos
if rand() < parametros.prob_mutacion
switch lower(parametros.tipo_mutacion)
case 'bit_flip'
descendencia(i, :) = mutacion_bit_flip(descendencia(i, :), problema);
case 'gaussiana'
descendencia(i, :) = mutacion_gaussiana(descendencia(i, :), problema);
case 'intercambio'
descendencia(i, :) = mutacion_intercambio(descendencia(i, :));
end
end
end
end

function individuo_mutado = mutacion_bit_flip(individuo, problema)
individuo_mutado = individuo;
pos = randi(length(individuo));
individuo_mutado(pos) = 1 - individuo_mutado(pos);
end

function individuo_mutado = mutacion_gaussiana(individuo, problema)
individuo_mutado = individuo;
pos = randi(length(individuo));

if length(problema.limite_inf) == 1
limite_inf = problema.limite_inf;
limite_sup = problema.limite_sup;
else
limite_inf = problema.limite_inf(pos);
limite_sup = problema.limite_sup(pos);
end

individuo_mutado(pos) = individuo(pos) + problema.sigma*randn();
individuo_mutado(pos) = max(min(individuo_mutado(pos), limite_sup), limite_inf);
end

function individuo_mutado = mutacion_intercambio(individuo)
individuo_mutado = individuo;
if length(individuo) > 1
pos = randperm(length(individuo), 2);
temp = individuo_mutado(pos(1));
individuo_mutado(pos(1)) = individuo_mutado(pos(2));
individuo_mutado(pos(2)) = temp;
end
end

%% ================== FUNCIONES AUXILIARES ==================

function poblacion = crear_poblacion_inicial(tam_poblacion, problema)
tam_cromosoma = problema.tam_cromosoma;
poblacion = zeros(tam_poblacion, tam_cromosoma);

for i = 1:tam_poblacion
for j = 1:tam_cromosoma
if strcmp(problema.tipo, 'binario')
poblacion(i, j) = randi([0, 1]);
elseif strcmp(problema.tipo, 'entero')
poblacion(i, j) = randi([problema.limite_inf(j), problema.limite_sup(j)]);
elseif strcmp(problema.tipo, 'real')
poblacion(i, j) = problema.limite_inf(j) + rand() * (problema.limite_sup(j) - problema.limite_inf(j));
elseif strcmp(problema.tipo, 'permutacion')
poblacion(i, :) = randperm(tam_cromosoma);
break;
end
end
end
end

%% ======= FITNESS ===============
function fitness = fitness_tsp(individuo, problema)

distancia_total = 0;
n = length(individuo);

for i = 1:n-1
ciudad_actual = individuo(i);
ciudad_siguiente = individuo(i+1);
distancia_total = distancia_total + problema.distancias(ciudad_actual, ciudad_siguiente);
end

distancia_total = distancia_total + problema.distancias(individuo(end), individuo(1));

unicos = (10-length(unique(individuo))) * 100;

fitness = 1/(distancia_total+1+unicos);
end

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
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
        fprintf('Gen %d: Soluciones encontradas:\n', g);
        disp(unique(soluciones_generacion, 'rows'));
        end
        fprintf('Gen %d: Fitness = %.4f, Conflictos = %.0f\n', g, fit, min_conflictos(g));
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

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% Algoritmo genetico para ajustar parabola y = ax^2 + bx + c
clear all; close all; clc;

disp("a) Evalulacion de aptitud: ECM");
disp("b) Seleccion: ruleta");
disp("c) Cruce: aritmetico");
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

%% Cruce Aritmetico
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

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% Algoritmo genetico para el Problema de la Mochila
clear all; close all; clc;

% Parametros del problema de la mochila
disp("Condiciones del problema: Pesos");
pesos = [10, 20, 30, 15, 25, 5, 35, 12, 22, 18];
disp(pesos);
disp("Condiciones del problema: valores");
valores = [60, 100, 120, 70, 90, 30, 150, 50, 80, 110];
disp(valores);
disp("Condiciones del problema: Capacidad");
capacidad_maxima = 100;
disp(capacidad_maxima);

% Parametros del algoritmo genetico
parametros.tam_poblacion = 100;
parametros.max_generaciones = 1000;
parametros.prob_cruce = 0.8;
parametros.prob_mutacion = 0.02;

disp("a) Evaluacion de aptitud: funcion de valor con penalizacion por exceso de peso");
disp("b) Seleccion: torneo");
disp("c) Cruce: un punto");
disp("d) Mutacion: bit flip");

% Configuración del problema
problema.tam_cromosoma = length(pesos);
problema.limite_inf = zeros(1, problema.tam_cromosoma);
problema.limite_sup = ones(1, problema.tam_cromosoma);
problema.sigma = 0.1;

problema.pesos = pesos;
problema.valores = valores;
problema.capacidad_maxima = capacidad_maxima;

[mejorIndividuo, mejor_fitness] = algoritmo_genetico(problema, parametros);

%% Mostrar resultados
fprintf('\n=== RESULTADOS MOCHILA ===\n');
fprintf('Mejor solución: %s\n', mat2str(mejorIndividuo));
fprintf('Valor total: %.2f\n', mejor_fitness);
fprintf('Peso total: %.2f\n', sum(pesos(mejorIndividuo==1)));
fprintf('Items seleccionados: %s\n', mat2str(find(mejorIndividuo==1)));

%% ================== FUNCION AG ==================
function [mejor_i, mejor_apti] = algoritmo_genetico(problema, parametros)
    %% Parametros del ag
    n = parametros.tam_poblacion;
    Generacion = parametros.max_generaciones;
    
    %% Poblacion inicial
    poblacion = zeros(n, problema.tam_cromosoma);
    for k = 1:n
        poblacion(k, :) = randi([0 1], 1, problema.tam_cromosoma);
    end
    disp("Poblacion inicial"); disp(poblacion);

    %% EVOLUCIÓN
    for gen = 1:Generacion
        aptitud = zeros(n, 1);
        soluciones_generacion = [];
        
        for i = 1:n
            % Calcular fitness
            valor_total = sum(problema.valores(poblacion(i,:) == 1));
            peso_total = sum(problema.pesos(poblacion(i,:) == 1));
            
            % Determinar si es solucion valida (no excede capacidad)
            if peso_total <= problema.capacidad_maxima
                aptitud(i) = valor_total;
            else
                aptitud(i) = 0; % Penalizacion por exceder capacidad
            end
            if peso_total <= problema.capacidad_maxima & peso_total >= 90
                soluciones_generacion = [soluciones_generacion; poblacion(i,:) valor_total peso_total];
            end
        end
        
        [mejor_fitness_gen, idx_mejor] = max(aptitud);
        elite = poblacion(idx_mejor,:);
        
        if mod(gen, 100) == 0
            % Mostrar soluciones validas encontradas
            if ~isempty(soluciones_generacion)
                fprintf('\nGeneracion %d: Mejores soluciones validas encontradas:\n', gen);
                
                % Eliminar duplicados y ordenar por valor)
                [soluciones_unicas, idx_unicas] = unique(soluciones_generacion(:,1:problema.tam_cromosoma), 'rows', 'stable');
                valores = soluciones_generacion(idx_unicas, end-1); 

                % Mostrar las 5 mejores
                num_mostrar = min(5, size(soluciones_unicas, 1));
                
                for j = 1:num_mostrar
                    idx_original = idx_unicas(j);
                    fprintf('Solucion %d: %s | Valor: %.1f | Peso: %.1f\n', j, ...
                        mat2str(soluciones_generacion(idx_original,1:problema.tam_cromosoma)), ...
                        soluciones_generacion(idx_original,end-1), ...
                        soluciones_generacion(idx_original,end));
                end
                
                % Mostrar mensaje si hay más soluciones
                if size(soluciones_unicas, 1) > 5
                    fprintf('... y %d soluciones mas (Valor minimo mostrado: %.1f)\n', ...
                        size(soluciones_unicas, 1) - 5, soluciones_unicas(5));
                end
            else
                fprintf('Generacion %d: Mejor fitness = %.4f\n', gen, mejor_fitness_gen);
            end
        end
        
        % seleccion
        pp = Torneo(poblacion, aptitud);
        descendencia = zeros(size(poblacion));

        for i = 2:2:n-1
            % cruce
            h1 = pp(i, :);
            h2 = pp(i+1, :);

            if rand() < parametros.prob_cruce
                [h1, h2] = UnPunto(h1, h2);
            end
            
            % Mutacion
            if rand() < parametros.prob_mutacion
                h1 = mutacionBitFlip(h1, problema);
            end
            if rand() < parametros.prob_mutacion
                h2 = mutacionBitFlip(h2, problema);
            end
            descendencia(i,:) = h1;
            descendencia(i+1,:) = h2;
        end
        
        % Elitismo
        descendencia(1,:) = elite;
        poblacion = descendencia;
    end
    
    % Evaluar poblacion final
    fitnessFinal = evaluar_poblacion(poblacion, problema);
    [mejor_apti, idx_mejor] = max(fitnessFinal);
    mejor_i = poblacion(idx_mejor, :);
end

%% ================== funciones ==================
function aptitu = evaluar_poblacion(p, problema)
    aptitu = zeros(size(p, 1), 1);
    for i = 1:size(p, 1)
        valor_total = sum(problema.valores(i == 1));
        peso_total = sum(problema.pesos(i == 1));
        if peso_total > problema.capacidad_maxima
            aptitu(i) = 0;
        else
            aptitu(i) = valor_total;
        end
    end
end
%% ================== seleccion: torneo ==================
function p = Torneo(poblacion, fitness)
    tamTorneo = 3;
    n = size(poblacion, 1);
    p = zeros(n, size(poblacion, 2));
    for k = 1:n
        competidores = randperm(n, tamTorneo);
        aptituCompetidores = fitness(competidores);
        [~, ubica_ganador] = max(aptituCompetidores);
        ganador = competidores(ubica_ganador);
        p(k,:) = poblacion(ganador,:);
    end
end

%% ================== cruce: un punto ==================
function [h1, h2] = UnPunto(p1, p2)
    tamCrom = length(p1);
    pCorte = randi([1, tamCrom-1]);
    h1 = [p1(1:pCorte), p2(pCorte+1:end)];
    h2 = [p2(1:pCorte), p1(pCorte+1:end)];
end

%% ================== mutacion: bitflip ==================
function i_Mutado = mutacionBitFlip(i, ~)
    i_Mutado = i;
    pos = randi(length(i));
    i_Mutado(pos) = 1 - i_Mutado(pos);
end

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% Viajero
clear all; close all; clc;
% Parametros
tam_poblacion = 100; max_generaciones = 500;
prob_cruce = 0.8; prob_mutacion = 0.2; n_ciudades = 10;

disp("a) Evaluacion de aptitud: distancia minima");
disp("b) Seleccion: torneo");
disp("c) Cruce: OX");
disp("d) Mutacion: intercambio");

% Matriz de distancias y coordenadas
rng(4);
distancias = randi([1 100], n_ciudades, n_ciudades);
distancias = triu(distancias) + triu(distancias, 1)';
coordenadas = rand(n_ciudades, 2) * 100;

% Poblacion inicial
poblacion = zeros(tam_poblacion, n_ciudades);
for i = 1:tam_poblacion
    poblacion(i, :) = randperm(n_ciudades);
end
disp("Poblacion inicial"); disp(poblacion);

% Algoritmos geneticos
historial = zeros(max_generaciones, 1);
mejores_distancias = zeros(max_generaciones, 1);
umbral_distancia = 200; % Umbral para considerar una solucion buena
for g = 1:max_generaciones
    % Evaluar aptitud
    fitness = zeros(tam_poblacion, 1);
    soluciones_generacion = [];  % Para almacenar soluciones buenas
    for i = 1:tam_poblacion
        fitness(i) = fitness_tsp(poblacion(i, :), distancias);
        % Si la distancia es menor al umbral, guardar como solucion
        d = (1/fitness(i)) -1;
        if d <= umbral_distancia
            soluciones_generacion = [soluciones_generacion; poblacion(i,:) d];
        end
    end
    [mejor_fit, idx] = max(fitness);
    historial(g) = mejor_fit;
    
    % Mostrar progreso
    if mod(g, 50) == 0
        fprintf('\nGeneración %d: Soluciones con distancia <= %.2f:\n', g, umbral_distancia);
        if ~isempty(soluciones_generacion)
            [~, idx_unicas] = unique(soluciones_generacion(:,1:n_ciudades), 'rows', 'stable');
            soluciones_generacion = soluciones_generacion(idx_unicas,:);
            
            % Ordenar por distancia
            [distancias_ordenadas, idx_orden] = sort(soluciones_generacion(:,end));
            soluciones_ordenadas = soluciones_generacion(idx_orden, 1:end-1);

            num_mostrar = min(5, size(soluciones_ordenadas, 1));
            for s = 1:num_mostrar
                fprintf('Ruta %d: %s | Distancia: %.2f\n', s, ...
                    mat2str(soluciones_ordenadas(s,:)), distancias_ordenadas(s));
            end
            if size(soluciones_ordenadas, 1) > 5
                fprintf('... (+%d mas)\n', size(soluciones_ordenadas, 1)-5);
            end
        else
            fprintf('No se encontraron soluciones bajo el umbral en esta generacion\n');
        end
    end
    
    % Seleccion por torneo
    padres = seleccion_torneo(poblacion, fitness);
    
    % Cruce OX
    descendencia = cruceOX(padres, prob_cruce, n_ciudades);
    
    % Mutacion
    for i = 1:tam_poblacion
        if rand < prob_mutacion
            descendencia(i, :) = mutacion_intercambio(descendencia(i, :));
        end
    end
    
    % Elitismo y Reemplazo
    descendencia(1, :) = poblacion(idx, :);
    poblacion = descendencia;
end

%% Mostrar resultados
fprintf('\n=== RESULTADOS FINALES ===\n');
[mejor_fit, idx] = max(fitness);
mejor_ind = poblacion(idx, :);
fprintf('Mejor ruta: %s\nDistancia: %.2f\n', mat2str(mejor_ind), 1/mejor_fit - 1);

% Grafica de convergencia
figure; plot(1:max_generaciones, historial, 'LineWidth', 2);
xlabel('Generacion'); ylabel('fitness'); title('Convergencia'); grid on;

% Grafica de ruta
mejor_ruta = [mejor_ind mejor_ind(1)];
figure; plot(coordenadas(mejor_ruta, 1), coordenadas(mejor_ruta, 2), '-or');
% Añadir etiquetas de ciudades
for i = 1:n_ciudades
    text(coordenadas(i,1), coordenadas(i,2), sprintf('Ciudad %d', i));
end
title('Ruta Final'); xlabel('X'); ylabel('Y'); axis equal; grid on;

% Funcion de aptitud
function fit = fitness_tsp(ind, dist)
    d = 0;
    for i = 1:length(ind)-1
        d = d + dist(ind(i), ind(i+1));
    end
    d = d + dist(ind(end), ind(1));
    fit = 1 / (1 + d);
end

% Seleccion por torneo
function padres = seleccion_torneo(pobl, fit)
    n = size(pobl, 1);
    k=2;
    padres = zeros(size(pobl));
    for i = 1:n
        idxs = randsample(n, k);
        [~, best] = max(fit(idxs));
        padres(i,:) = pobl(idxs(best), :);
    end
end

% Cruce OX
function hijos = cruceOX(pobl, prob, d)
    n = size(pobl, 1);
    hijos = pobl;
    for i = 1:2:n-1
        if rand < prob
            % Puntos de cruce
            p1 = randi([1 d-1]);
            p2 = randi([p1+1 d]);
            
            % Hijo 1
            hijo1 = zeros(1, d);
            hijo1(p1:p2) = pobl(i, p1:p2);
            j = mod(p2, d) + 1;
            for k = 1:d
                g = pobl(i+1, mod(p2 + k - 1, d) + 1);
                if ~ismember(g, hijo1)
                    hijo1(j) = g;
                    j = mod(j, d) + 1;
                end
            end
            
            % Hijo 2
            hijo2 = zeros(1, d);
            hijo2(p1:p2) = pobl(i+1, p1:p2);
            j = mod(p2, d) + 1;
            for k = 1:d
                g = pobl(i, mod(p2 + k - 1, d) + 1);
                if ~ismember(g, hijo2)
                    hijo2(j) = g;
                    j = mod(j, d) + 1;
                end
            end
            
            hijos(i, :) = hijo1;
            hijos(i+1, :) = hijo2;
        end
    end
end

% Mutacion por intercambio
function ind = mutacion_intercambio(ind)
    if length(ind) > 1
        pos = randperm(length(ind), 2);
        temp = ind(pos(1));
        ind(pos(1)) = ind(pos(2));
        ind(pos(2)) = temp;
    end
end

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% Parametros del problema (Job Sequencing)
clear all; close all; clc;

% Datos del problema
n_jobs = 10;
profits = [100 19 27 25 15 90 30 10 40 70];
deadlines = [2 1 2 1 3 5 7 9 3 2];

% Parametros del algoritmo
p.tam_poblacion = 50;
p.max_generaciones = 100;
p.prob_cruce = 0.8;
p.prob_mutacion = 0.1;

% Configuracion del problema
prob.tam_cromosoma = n_jobs;
prob.profits = profits;
prob.deadlines = deadlines;

disp("a) Aptitud: suma de beneficios de tareas a tiempo");
disp("b) Seleccion: torneo");
disp("c) Cruce: orden (OX)");
disp("d) Mutacion: intercambio");

% Ejecutar algoritmo genetico
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

% === Algoritmo Genetico ===
function [mejor, mejor_apt] = algoritmo_genetico(prob, p)
    % Poblacion inicial
    n = p.tam_poblacion;
    poblacion = zeros(n, prob.tam_cromosoma);
    for i = 1:n
        poblacion(i,:) = randperm(prob.tam_cromosoma);
    end
    disp("Poblacion inicial"); disp(poblacion);
    
    % Evolucion
    for gen = 1:p.max_generaciones
        % Evaluar aptitud y almacenar soluciones validas
        aptitud = zeros(n, 1);
        soluciones_generacion = []; 
        for i = 1:n
            [aptitud(i), tiempo] = calcular_aptitud(poblacion(i,:), prob);
            if aptitud(i) >= 300
                soluciones_generacion = [soluciones_generacion; poblacion(i,:)];
            end
        end
        
        % Guardar el mejor
        [mejor_apt, idx] = max(aptitud);
        elite = poblacion(idx,:);
        
        if mod(gen, 20) == 0
        % Mostrar soluciones validas encontradas
            if ~isempty(soluciones_generacion)
                fprintf('\nGeneracion %d: Soluciones validas encontradas:\n', gen);
                soluciones_unicas = unique(soluciones_generacion, 'rows');
                
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    seq = soluciones_unicas(s,:);
                    [apt, ~] = calcular_aptitud(seq, prob);
                    fprintf('Secuencia: %s - Beneficio: %d\n', mat2str(seq), apt);
                end
                
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones mas\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generacion %d: Mejor beneficio = %.2f\n', gen, mejor_apt);
            end
        end
        
        
        % Seleccion por torneo
        padres = zeros(n, prob.tam_cromosoma);
        k = 3;
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
        
        % Mutacion por intercambio
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
    
    % Mejor solucion final
    aptitud = zeros(n, 1);
    for i = 1:n
        aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
    end
    [mejor_apt, idx] = max(aptitud);
    mejor = poblacion(idx,:);
end

% === Funcion de Aptitud modificada ===
function [apt, tiempo_total] = calcular_aptitud(seq, prob)
    apt = 0;
    tiempo_total = 0;
    for i = 1:length(seq)
        tiempo_total = tiempo_total + 1;
        if tiempo_total <= prob.deadlines(seq(i))
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

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================

%% Parametros del problema Himmelblau
clear all; close all; clc;

% Parametros del problema
prob.tam_cromosoma = 2; % [x, y]
prob.limites = [-5, 5]; % Dominio para x e y
prob.lambda = 1000; % Factor de penalizacion para restricciones

% Parametros del algoritmo genetico
p.tam_poblacion = 50;
p.max_generaciones = 100;
p.prob_cruce = 0.8;
p.prob_mutacion = 0.1;
p.sigma_mutacion = 0.1; % Desviacion estandar para mutacion gaussiana
p.tam_torneo = 3; % Tamano del torneo

disp("a) Aptitud: funcion de Himmelblau con penalizacion por dos restricciones");
disp("b) Seleccion: torneo");
disp("c) Cruce: aritmetico");
disp("d) Mutacion: gaussiana");

% Ejecutar algoritmo genetico
[best_sol, best_apt, apt_por_generacion] = algoritmo_genetico(prob, p);

% Mostrar resultados
fprintf('\n=== Resultados finales: Mejor Solución ===\n');
fprintf('x: %.4f, y: %.4f\n', best_sol(1), best_sol(2));
fprintf('Valor de Himmelblau: %.4f\n', himmelblau(best_sol));
fprintf('Restriccion 1 (x^2 + y^2 - 4): %.4f (Factible si <= 0)\n', constraint1(best_sol));
fprintf('Restriccion 2 (x + y - 1): %.4f (Factible si >= 0)\n', constraint2(best_sol));
fprintf('Aptitud final: %.4f\n', best_apt);

% Graficar la evolucion de la aptitud por generacion
figure;
plot(1:p.max_generaciones, apt_por_generacion, 'LineWidth', 2);
xlabel('Generacion');
ylabel('Mejor Aptitud');
title('Evolucion de la Aptitud por Generacion');
grid on;

% Graficar la funcion de Himmelblau en 3D
figure;
[x, y] = meshgrid(linspace(-5, 5, 100), linspace(-5, 5, 100)); % Crear malla de puntos
z = (x.^2 + y - 11).^2 + (x + y.^2 - 7).^2; % Evaluar Himmelblau en cada punto

% Gráfica de superficie
surf(x, y, z, 'EdgeColor', 'none');
title('Funcion de Himmelblau');
xlabel('x'); ylabel('y'); zlabel('f(x, y)'); colorbar;
hold on;
plot3(best_sol(1), best_sol(2), himmelblau(best_sol), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

% === Algoritmo Genetico ===
function [mejor, mejor_apt, apt_por_generacion] = algoritmo_genetico(prob, p)
    % Poblacion inicial
    n = p.tam_poblacion;
    poblacion = zeros(n, prob.tam_cromosoma);
    for i = 1:n
        poblacion(i,:) = prob.limites(1) + (prob.limites(2) - prob.limites(1)) * rand(1, prob.tam_cromosoma);
    end
    disp("Poblacion inicial"); disp(poblacion);
    apt_por_generacion = zeros(p.max_generaciones, 1);
    
    % Evolucion
    for gen = 1:p.max_generaciones
        % Evaluar aptitud y almacenar soluciones factibles
        aptitud = zeros(n, 1);
        soluciones_factibles = [];  % Almacenar soluciones factibles de esta generacion
        
        for i = 1:n
            aptitud(i) = calcular_aptitud(poblacion(i,:), prob);
            
            % Verificar si es solucion factible
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
                fprintf('\nGeneracion %d: Soluciones factibles encontradas:\n', gen);
                soluciones_unicas = unique(soluciones_factibles, 'rows');
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    sol = soluciones_unicas(s,:);
                    fprintf('x: %.4f, y: %.4f - Himmelblau: %.4f\n', ...
                            sol(1), sol(2), himmelblau(sol));
                end
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones mas\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generacion %d: Mejor aptitud = %.4f\n', gen, mejor_apt);
            end
        end

        % Selección por torneo
        padres = zeros(n, prob.tam_cromosoma);
        for i = 1:n
            competidores = randperm(n, p.tam_torneo);
            [~, idx] = min(aptitud(competidores)); % Mejor competidor
            padres(i,:) = poblacion(competidores(idx),:);
        end
        
        % Cruce aritmetico
        hijos = padres;
        for i = 1:2:n-1
            if rand() < p.prob_cruce
                alpha = rand();
                hijos(i,:) = alpha * padres(i,:) + (1-alpha) * padres(i+1,:);
                hijos(i+1,:) = (1-alpha) * padres(i,:) + alpha * padres(i+1,:);
            end
        end
        
        % Mutacion gaussiana
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

% === Funcion de Himmelblau ===
function val = himmelblau(sol)
    x = sol(1); y = sol(2);
    val = (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
end

% === Restriccion 1: x^2 + y^2 <= 4 ===
function val = constraint1(sol)
    x = sol(1); y = sol(2);
    val = x^2 + y^2 - 4; % <= 0
end

% === Restriccion 2: x + y >= 1 ===
function val = constraint2(sol)
    x = sol(1); y = sol(2);
    val = 1 - (x + y); % <= 0 (convertido a forma estándar)
end

% === Funcion de Aptitud ===
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

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================

% Algoritmo genetico para Optimizacion de Ubicacion de Antenas en Telecomunicaciones
clear all; close all; clc;

% Parametros del algoritmo genetico
parametros.tam_poblacion = 50;
parametros.max_generaciones = 100;
parametros.prob_cruce = 0.8;
parametros.prob_mutacion = 0.1;

disp("a) Evaluacion de aptitud: ");
disp("b) Seleccion: torneo");
disp("c) Cruce: un punto");
disp("d) Mutacion: gaussiana");

problema.numAntenas = 3;
problema.numTargets = 50;
problema.areaSize = 40;
problema.targets = problema.areaSize * rand(problema.numTargets, 2);
problema.maxRange = [10 6 14];
problema.funcion_fitness = @CoberturAntenas;

problema.tam_cromosoma = problema.numAntenas * 2;
problema.limite_inf = [0 0 0 0 0 0];
problema.limite_sup = problema.areaSize * [1 1 1 1 1 1];
problema.sigma = 0.5;

[mejorIndividuo, mejor_fitness] = algoritmo_genetico(problema, parametros);

%% ====================================
antenasMEjores = reshape(mejorIndividuo, 2, problema.numAntenas)';

fprintf('\n=== UBICACIONES OPTIMAS DE ANTENAS ===\n');
for i = 1:problema.numAntenas
    fprintf('Antena %d: (%.3f, %.3f)\n', i, antenasMEjores(i,1), antenasMEjores(i,2));
end

% Visualizacion
figure; hold on;
for i = 1:problema.numAntenas
    theta = linspace(0, 2*pi, 100);
    plot(antenasMEjores(i,1) + problema.maxRange(i)*cos(theta), ...
         antenasMEjores(i,2) + problema.maxRange(i)*sin(theta), 'r--');
end
plot(problema.targets(:,1), problema.targets(:,2), 'kx');
plot(antenasMEjores(:,1), antenasMEjores(:,2), 'or');
% Anadir etiquetas de ciudades
for i = 1:problema.numAntenas
    text(antenasMEjores(i,1), antenasMEjores(i,2), sprintf('Antena %d', i));
end
xlim([0 problema.areaSize]); ylim([0 problema.areaSize]); grid on; axis equal;
title('Cobertura de antenas optimas');

%% ================== Funcion AG ==================

function [mejor_i, mejor_apti] = algoritmo_genetico(problema, parametros)
    %% Parametro ag
    n = parametros.tam_poblacion;
    Generacion = parametros.max_generaciones;
    
    %% Poblacion inicial
    poblacion = zeros(n, problema.tam_cromosoma);
    for k = 1:n
        for r = 1:problema.tam_cromosoma
            poblacion(k, r)=problema.limite_inf(r) + rand() * (problema.limite_sup(r) - problema.limite_inf(r));
        end
    end
    disp("Poblacion inicial"); disp(poblacion);
    
    %% Evolucion
    for gen = 1:Generacion
        aptitud = zeros(n, 1);
        soluciones_optimas = [];
        
        for i = 1:n
            aptitud(i) = evaluar_poblacion(poblacion(i,:), problema);
            if aptitud(i) >= 0.6
                soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
            end
        end
        
        [mejor_fitness_gen, idx_mejor] = max(aptitud);
        elite = poblacion(idx_mejor,:);
        
        if mod(gen, 20) == 0
            % Mostrar progreso y soluciones optimas encontradas
            if ~isempty(soluciones_optimas)
                fprintf('\nGeneracion %d: Soluciones optimas encontradas\n', gen);
                % Mostrar solo soluciones unicas
                soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
                num_mostrar = min(5, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    sol = soluciones_unicas(s,:);
                    antenas = reshape(sol, 2, problema.numAntenas)';
                    fprintf('Solucion %d:\n', s);
                    for a = 1:problema.numAntenas
                        fprintf('  Antena %d: (%.2f, %.2f)\n', a, antenas(a,1), antenas(a,2));
                    end
                    fprintf('Fitness: %.4f\n\n', evaluar_poblacion(sol, problema));
                end
                if size(soluciones_unicas, 1) > 5
                    fprintf('... y %d soluciones mas\n\n', size(soluciones_unicas, 1) - 5);
                end
            else
                fprintf('Generacion %d: Mejor fitness = %.4f\n', gen, mejor_fitness_gen);
            end
        end

        descendencia = zeros(size(poblacion));
        
        % seleccion
        pp = zeros(n, problema.tam_cromosoma);
        k=3;
        for i = 1:n
            competidores = randperm(n, k);
            [~, idx] = max(aptitud(competidores));
            pp(i,:) = poblacion(competidores(idx),:);
        end

        for i = 2:2:n-1
            % cruce
            h1 = pp(i, :);
            h2 = pp(i+1, :);
            if rand() < parametros.prob_cruce
                [h1, h2] = UnPunto(h1, h2);
            end
            % mutacion
            h1 = mutar(h1, problema, parametros.prob_mutacion);
            h2 = mutar(h2, problema, parametros.prob_mutacion);
            
            descendencia(i,:) = h1;
            descendencia(i+1,:) = h2;
        end
        
        % Elitismo
        descendencia(1,:) = elite;
        poblacion = descendencia;
    end
    
    % Evaluar población final
    fitnessFinal = evaluar_poblacion(poblacion, problema);
    [mejor_apti, idx_mejor] = max(fitnessFinal);
    mejor_i = poblacion(idx_mejor, :);
end

%% ================== Funciones==================
function aptitu = evaluar_poblacion(p, problema)
    if size(p,1) > 1
        aptitu = zeros(size(p, 1), 1);
        for i = 1:size(p, 1)
            aptitu(i) = problema.funcion_fitness(p(i, :), problema);
        end
    else
        aptitu = problema.funcion_fitness(p, problema);
    end
end

%% ================== cruce ==================
function [h1, h2] = UnPunto(p1, p2)
    punto = randi(length(p1)-1);
    h1 = [p1(1:punto) p2(punto+1:end)];
    h2 = [p2(1:punto) p1(punto+1:end)];
end

%% ================== mutacion ==================
function individuo = mutar(ind, problema, prob)
    if rand < prob
        gen = randi(length(ind));
        ind(gen) = ind(gen) + problema.sigma*randn();
        ind(gen) = min(max(ind(gen), problema.limite_inf(gen)), problema.limite_sup(gen));
    end
    individuo = ind;
end

%% ======= fitness ====
function aptitu = CoberturAntenas(i, problema)
    antenas = reshape(i, 2, problema.numAntenas)';
    objCubiertos = false(problema.numTargets, 1);

    for i = 1:problema.numTargets
        objeto = problema.targets(i, :);
        for j = 1:problema.numAntenas
            antena = antenas(j, :);
            d = norm(objeto - antena);
            if d <= problema.maxRange(j)
                objCubiertos(i) = true;
                break;
            end
        end
    end
    aptitu = sum(objCubiertos) / problema.numTargets;
end

%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% Canales
clear all; close all; clc;

%% Parametros del problema y algoritmo
canales_totales = 15; % Numero total de canales disponibles
celdas = 10;  % Numero de celdas a asignar

% Configuracion del algoritmo genetico
parametros = struct();
parametros.tam_poblacion = 50;       % Tamano de la poblacion
parametros.max_generaciones = 100;   % Numero maximo de generaciones
parametros.prob_cruce = 0.8;         % Probabilidad de cruce
parametros.prob_mutacion = 0.1;      % Probabilidad de mutacion

disp("Configuracion del algoritmo:");
disp("a) Evaluacion de aptitud: Minimizar interferencia");
disp("b) Seleccion: Torneo");
disp("c) Cruce: Uniforme");
disp("d) Mutacion: Intercambio");

%% Definicion del problema
problema = struct();
problema.nombre = 'Asignacion optima de canales';
problema.tam_cromosoma = celdas;
problema.funcion_fitness = @(individuo) fitness_canal(individuo, canales_totales);
problema.limite_inf = ones(1, celdas);
problema.limite_sup = canales_totales*ones(1, celdas);

%% Ejecucion del algoritmo genetico
[mejor_asignacion, mejor_fitness, historico_fitness] = algoritmo_genetico(problema, parametros);

%% Resultados
fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor fitness (minimizar interferencia): %.4f\n', mejor_fitness);
fprintf('Asignación óptima de canales:\n');
disp(mejor_asignacion);

% Grafico de evolucion
figure;
plot(1:parametros.max_generaciones, historico_fitness, 'LineWidth', 2);
title('Evolucion del Fitness');
xlabel('Generación');
ylabel('Mejor Fitness');
grid on;

%% ========= funcion ag ==========
function [mejor_ind, mejor_fit, historico_fitness] = algoritmo_genetico(problema, parametros)
    % Inicializacion
    n = parametros.tam_poblacion;
    max_gen = parametros.max_generaciones;
    historico_fitness = zeros(max_gen, 1);
    
    % Crear poblacion inicial aleatoria
    poblacion = randi([problema.limite_inf(1), problema.limite_sup(1)], n, problema.tam_cromosoma);
    disp("Poblacion inicial"); disp(poblacion);
    
    for gen = 1:max_gen
        % Evaluar fitness y almacenar soluciones optimas
        fitness = zeros(n, 1);
        soluciones_optimas = [];
        
        for i = 1:n
            fitness(i) = problema.funcion_fitness(poblacion(i,:));
            % Consideramos solucion optima si no hay interferencia
            if fitness(i) <=0.035
                soluciones_optimas = [soluciones_optimas; poblacion(i,:)];
            end
        end
        
        [mejor_fit, idx] = max(fitness);
        elite = poblacion(idx,:);
        historico_fitness(gen) = mejor_fit;
        
        if mod(gen, 20) == 0
            % Mostrar soluciones optimas encontradas
            if ~isempty(soluciones_optimas)
                fprintf('\nGeneracion %d: Soluciones optimas encontradas (Fitness = %.4f):\n', gen, mejor_fit);
                soluciones_unicas = unique(soluciones_optimas, 'rows', 'stable');
                num_mostrar = min(6, size(soluciones_unicas, 1));
                
                for s = 1:num_mostrar
                    fprintf('Solucion %d: ', s);
                    disp(soluciones_unicas(s,:));
                end
                if size(soluciones_unicas, 1) > 6
                    fprintf('... y %d soluciones mas\n', size(soluciones_unicas, 1) - 6);
                end
            else
                fprintf('Generación %d - Mejor fitness: %.4f\n', gen, mejor_fit);
            end
        end

        % Seleccion por torneo
        padres = zeros(size(poblacion));
        for i = 1:n
            k=3;
            competidores = randperm(n, k);
            [~, idx_ganador] = max(fitness(competidores));
            padres(i,:) = poblacion(competidores(idx_ganador),:);
        end
        
        % Cruce uniforme
        descendencia = padres;
        for i = 1:2:n-1
            if rand() < parametros.prob_cruce
                mascara = rand(1, problema.tam_cromosoma) > 0.5;
                temp1 = descendencia(i,:);
                temp2 = descendencia(i+1,:);
                descendencia(i,mascara) = temp2(mascara);
                descendencia(i+1,mascara) = temp1(mascara);
            end
        end
        
        % Mutacion por intercambio
        for i = 1:n
            if rand() < parametros.prob_mutacion
                pos = randperm(problema.tam_cromosoma, 2);
                temp = descendencia(i,pos(1));
                descendencia(i,pos(1)) = descendencia(i,pos(2));
                descendencia(i,pos(2)) = temp;
            end
        end
        
        % Elitismo
        descendencia(1,:) = elite;
        poblacion = descendencia;
    end
    
    % Resultado final
    mejor_ind = elite;
end

%% ========= FUNCIÓN DE FITNESS ==========
function fitness = fitness_canal(asignacion, canales_totales)
    % Penalizar interferencia entre celdas adyacentes
    interferencia = 0;
    for i = 2:length(asignacion)
        dif = abs(asignacion(i) - asignacion(i-1));
        if dif == 0 % Misma frecuencia
            interferencia = interferencia + 100; 
        else
            interferencia = interferencia + 1/dif; % Menor diferencia
        end
    end
    
    % Penalizar uso excesivo de canales
    canales_usados = length(unique(asignacion));
    penalizacion = (canales_totales - canales_usados)^2;
    
    % Fitness inverso a la interferencia total
    fitness = 1/(interferencia + penalizacion + 0.001);
end



%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================
%% ======================================================================

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
                disp(reshape(soluciones_unicas(s,:), num_reps, 2)'); % Mostrar en formato coordenadas
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
disp('Coordenadas (fila,columna) y alcances de los repetidores:');
reps_info = [mejor_solucion(1:num_reps)' mejor_solucion(num_reps+1:end)' alcances'];
disp(reps_info);

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