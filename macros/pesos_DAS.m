function [W] = pesos_DAS(d_n, tn, freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función que calcula los pesos del beamformer Delay & Sum para cada
% frecuencia.
% Argumentos de entrada: 
% d_n: factor multiplicador dependiendo del tipo de onda aplicado
% tn: vector de retardos
% freq: rango de frecuencias a evaluar
% Argumentos de salida: Pesos del beamformer para cada frecuencia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(tn); %Número de elementos del array
flim = length(freq);     %Barrido de frecuencias
W = zeros(flim,N); %Inicializamos el vector de pesos

    for f = 1:flim       
        for i = 1:N

            W(f,i) = d_n(i)*(1/N)*exp(1j*2*pi*tn(i)*freq(f)); % Matriz 129x7

        end
    end
end