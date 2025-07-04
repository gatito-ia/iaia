%% Firewall - Optimización de Reglas
clear; clc; close all;

%% Definición del problema
num_reglas = 10;       % Número de reglas en el firewall
num_campos = 5;        % Campos por regla: [IP_origen, IP_destino, Puerto_origen, Puerto_destino, Protocolo]
tipo_campos = [1, 1, 2, 2, 3]; % 1=IP, 2=Puerto, 3=Protocolo
trafico_total = 1000;  % Total de paquetes a evaluar

%% Parámetros del Algoritmo Genético
pob_size = 50;        % Tamaño de la población
max_gens = 100;       % Número máximo de generaciones
pc = 0.8;             % Probabilidad de cruce
pm = 0.1;             % Probabilidad de mutación

disp("a) Aptitud: Minimizar tráfico malicioso permitido");
disp("b) Selección: Torneo");
disp("c) Cruce: Dos puntos");
disp("d) Mutación: Aleatoria por reemplazo en campos");

%% Inicialización de población
% Cada individuo representa un conjunto de reglas de firewall
% Cada regla tiene 5 campos: [IP_origen, IP_destino, Puerto_origen, Puerto_destino, Protocolo]
poblacion = inicializar_poblacion(pob_size, num_reglas, num_campos, tipo_campos);
disp("Población inicial generada");

%% Evolución
mejor_aptitud = zeros(max_gens, 1);

for gen = 1:max_gens
    % Evaluar población
    aptitudes = zeros(pob_size, 1);
    for i = 1:pob_size
        aptitudes(i) = evaluar_firewall(poblacion(i,:), num_reglas, num_campos, tipo_campos, trafico_total);
    end
    
    % Guardar el mejor
    [mejor_aptitud(gen), idx] = min(aptitudes);
    mejor_individuo = poblacion(idx,:);
    
    % Mostrar progreso cada 20 generaciones
    if mod(gen, 20) == 0
        fprintf('Generación %d: Mejor firewall permite %d paquetes maliciosos\n', gen, mejor_aptitud(gen));
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
            puntos = sort(randperm(num_reglas*num_campos-1, 2));
            temp = hijos(i, puntos(1)+1:puntos(2));
            hijos(i, puntos(1)+1:puntos(2)) = hijos(i+1, puntos(1)+1:puntos(2));
            hijos(i+1, puntos(1)+1:puntos(2)) = temp;
        end
    end
    
    % Mutación
    for i = 1:pob_size
        if rand < pm
            % Seleccionar una regla aleatoria y un campo aleatorio para mutar
            regla_mutar = randi([1 num_reglas]);
            campo_mutar = randi([1 num_campos]);
            
            % Calcular posición en el vector
            pos = (regla_mutar-1)*num_campos + campo_mutar;
            
            % Mutar según el tipo de campo
            switch tipo_campos(campo_mutar)
                case 1 % IP
                    hijos(i, pos) = randi([0 255]);
                case 2 % Puerto
                    hijos(i, pos) = randi([0 65535]);
                case 3 % Protocolo
                    hijos(i, pos) = randi([0 2]); % 0=TCP, 1=UDP, 2=ICMP
            end
        end
    end
    
    % Elitismo (mantener el mejor individuo)
    hijos(1,:) = mejor_individuo;
    poblacion = hijos;
end

%% Resultados finales
[mejor_trafico, idx] = min(aptitudes);
mejor_firewall = poblacion(idx,:);

fprintf('\n=== RESULTADOS FINALES ===\n');
fprintf('Mejor firewall encontrado permite %d paquetes maliciosos (de %d)\n', mejor_trafico, trafico_total);
fprintf('Reglas del firewall óptimo:\n\n');

% Mostrar reglas formateadas
for r = 1:num_reglas
    inicio = (r-1)*num_campos + 1;
    fin = r*num_campos;
    regla = mejor_firewall(inicio:fin);
    
    fprintf('Regla %d:\n', r);
    fprintf('  IP Origen: %d.%d.%d.%d\n', regla(1), regla(1), regla(1), regla(1));
    fprintf('  IP Destino: %d.%d.%d.%d\n', regla(2), regla(2), regla(2), regla(2));
    fprintf('  Puerto Origen: %d\n', regla(3));
    fprintf('  Puerto Destino: %d\n', regla(4));
    fprintf('  Protocolo: %s\n\n', obtener_protocolo(regla(5)));
end

%% Visualización de convergencia
figure;
plot(1:max_gens, mejor_aptitud, 'r-', 'LineWidth', 1.5);
xlabel('Generación');
ylabel('Paquetes maliciosos permitidos');
title('Evolución de las reglas del firewall');
grid on;

%% Funciones auxiliares

function poblacion = inicializar_poblacion(pob_size, num_reglas, num_campos, tipo_campos)
    poblacion = zeros(pob_size, num_reglas*num_campos);
    
    for i = 1:pob_size
        individuo = [];
        for r = 1:num_reglas
            regla = [];
            for c = 1:num_campos
                switch tipo_campos(c)
                    case 1 % IP (simplificado como un byte)
                        regla = [regla randi([0 255])];
                    case 2 % Puerto
                        regla = [regla randi([0 65535])];
                    case 3 % Protocolo
                        regla = [regla randi([0 2])]; % 0=TCP, 1=UDP, 2=ICMP
                end
            end
            individuo = [individuo regla];
        end
        poblacion(i,:) = individuo;
    end
end

function trafico_malicioso = evaluar_firewall(individuo, num_reglas, num_campos, tipo_campos, trafico_total)
    % Simulamos el tráfico de red (en un caso real, esto vendría de datos reales)
    % Para este ejemplo, generamos tráfico aleatorio con un 20% de probabilidad de ser malicioso
    
    trafico_malicioso = 0;
    
    for p = 1:trafico_total
        % Generar un paquete aleatorio
        paquete = generar_paquete(tipo_campos);
        es_malicioso = rand() < 0.2; % 20% de probabilidad de ser malicioso
        
        % Verificar si el paquete pasa a través del firewall
        if pasa_firewall(individuo, num_reglas, num_campos, paquete)
            if es_malicioso
                trafico_malicioso = trafico_malicioso + 1;
            end
        end
    end
end

function pasa = pasa_firewall(individuo, num_reglas, num_campos, paquete)
    % Un paquete pasa si no es bloqueado por ninguna regla
    % Asumimos que las reglas son de bloqueo (deny)
    
    pasa = true;
    
    for r = 1:num_reglas
        inicio = (r-1)*num_campos + 1;
        fin = r*num_campos;
        regla = individuo(inicio:fin);
        
        % Verificar si el paquete coincide con la regla
        coincide = true;
        for c = 1:num_campos
            % Para simplificar, consideramos que 0 es "cualquiera"
            if regla(c) ~= 0 && regla(c) ~= paquete(c)
                coincide = false;
                break;
            end
        end
        
        if coincide
            pasa = false; % El paquete coincide con una regla de bloqueo
            return;
        end
    end
end

function paquete = generar_paquete(tipo_campos)
    paquete = [];
    for c = 1:length(tipo_campos)
        switch tipo_campos(c)
            case 1 % IP
                paquete = [paquete randi([0 255])];
            case 2 % Puerto
                paquete = [paquete randi([0 65535])];
            case 3 % Protocolo
                paquete = [paquete randi([0 2])];
        end
    end
end

function proto = obtener_protocolo(codigo)
    switch codigo
        case 0
            proto = 'TCP';
        case 1
            proto = 'UDP';
        case 2
            proto = 'ICMP';
        otherwise
            proto = 'DESCONOCIDO';
    end
end