function [d_n, tn] = onda_tipo(c, d, N, onda)
% De argumentos de entrada se admite:
% c: velocidad de propagación
% d: distancia (uniforme) entre los elementos del array
% N: Número de elementos
% onda: tipo de onda a computar esférica ('spherical') o plana

n=0:1:N-1; % índice de cada uno de los sensores en el array. Ej: 0, 1, 2, 3...

%Cáculo de la contribución de onda esférica si procede
if strcmp(onda, 'spherical')
    fprintf('Aplicada onda esférica \n');
    % Vector posición en el plano x de los sensores
    rxn = n*d;
    % Vector posición de la fuente (x,y) = (x_sensor_central, 1m)
    r_source = [rxn(ceil(length(rxn)/2)) 1];
    % Matriz de posiciones de cada sensor en el plano (x,y)
    r_n = transpose(padarray(rxn, [1 0], 'post'));
    % Vector de retardos
    tn = zeros(1,N);
    % Factor multiplicador del steering vector
    d_n = zeros(1,N);
    % Distancia de los sensores a la fuente
    r_s_n = zeros(1,N);
    
    
    % Cálculo de las distancias
    for i = 1:N
        % Distancia euclídea entre la fuente y el sensor n
        r_s_n(i) = norm(r_n(i, :) - r_source);
        % Se computa el factor con la distancia del sensor0 a la fuente y cada
        % una de las distancias de los sensores a la fuente
        d_n(i) = r_s_n(1) / r_s_n(i);
        %Retardo (distancia / velocidad)
        tn(i) = (r_s_n(i) - r_s_n(1)) / c;
    end
    
    
    %Si tomamos como referencia t0, restamos ese retardo al resto:
    tn = tn - tn(1);
    
% Suposición onda plana    
else
    fprintf('Aplicada onda plana \n');
    % Se computa el retardo asociado a cada sensor (para 90 siempre es 0)
    tn=(n*d*cos(phi))/c;
    d_n = ones(1,N);
end
end

