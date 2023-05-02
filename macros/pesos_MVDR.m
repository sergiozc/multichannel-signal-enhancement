function [W] = pesos_MVDR(tn, freq, corr_noise)
%Función que calcula los pesos del beamformer MVDR para cada f
% Se va a suponer campo homogéneo. (matriz de coherencia = matriz de corr)
% SE LE DEBERÍA PASAR COMO PARÁMETRO DE ENTRADA LA MATRIZ DE COHERENCIA,
% pero abajo se encuentran unos cálculos relativos a esta

N = length(tn); %Longitud del vector de retardos
flim = length(freq);     %Barrido de frecuencias
ds = zeros(flim,N); %Steering vector
W = zeros(flim,N); %vector de pesos


    for f = 1:flim       
        for i = 1:N
            ds(f,i) = exp(-1j*2*pi*tn(i)*freq(f)); % Matriz 129x7

        end
        W(f,:) = (inv(corr_noise) * transpose(ds(f,:))) / (conj(ds(f,:)) * inv(corr_noise) * transpose(ds(f,:)));
    end
end