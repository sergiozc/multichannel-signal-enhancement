function [W] = pesos_MVDR(d_n, tn, freq, corr_noise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función que calcula los pesos del beamformer MVDR para cada frecuencia
% Argumentos de entrada: 
% d_n: factor multiplicador dependiendo del tipo de onda aplicado
% tn: vector de retardos
% freq: rango de frecuencias a evaluar
% corr_noise: matriz de correlación espacial del ruido
% Argumentos de salida: Pesos del beamformer para cada frecuencia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(tn); % Número de elementos
flim = length(freq);     %Barrido de frecuencias
ds = zeros(flim,N); %Steering vector
W = zeros(flim,N); %vector de pesos
corr_inv = zeros(N,N); % matriz de correlación espacial del ruido para cada f
I = eye(N); % matriz identidad para hacer la operación inversa
mu = 0.001; % factor para realizar la operación inversa

    for f = 1:flim       
        for i = 1:N
            ds(f,i) = d_n(i) * exp(-1j*2*pi*tn(i)*freq(f));

        end
        corr_inv(:,:,f) = inv(corr_noise(:,:,f) + mu*I);
        W(f,:) = (corr_inv(:,:,f) * transpose(ds(f,:))) / (conj(ds(f,:)) * corr_inv(:,:,f) * transpose(ds(f,:)));
    end
end

