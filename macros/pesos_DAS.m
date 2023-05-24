function [W] = pesos_DAS(d_n,tn, freq)
%PESOS: Función que calcula los pesos del beamformer para cada f

N = length(tn); %longitud del vector de retardos
flim = length(freq);     %Barrido de frecuencias
W = zeros(flim,N); %Inicializamos el vector de pesos

    for f = 1:flim       
        for i = 1:N

            W(f,i) = d_n(i)*(1/N)*exp(1j*2*pi*tn(i)*freq(f)); % Matriz 129x7

        end
    end
end