function [W] = pesos_MVDR(d_n, tn, freq, corr_noise)
%Función que calcula los pesos del beamformer MVDR para cada f
%Como argumento de entrada se da el vector de retardos, los frecuencias a
%evaluar, la matriz de correlaciones espacial del ruido y la contribución
%de onda (plana o esférica) para el cálculo del steering vector

N = length(tn); %Longitud del vector de retardos
flim = length(freq);     %Barrido de frecuencias
ds = zeros(flim,N); %Steering vector
W = zeros(flim,N); %vector de pesos
corr_inv = zeros(N,N);
I = eye(N);
mu = 0.001;

    for f = 1:flim       
        for i = 1:N
            ds(f,i) = d_n(i) * exp(-1j*2*pi*tn(i)*freq(f));

        end
        corr_inv(:,:,f) = inv(corr_noise(:,:,f) + mu*I);
        W(f,:) = (corr_inv(:,:,f) * transpose(ds(f,:))) / (conj(ds(f,:)) * corr_inv(:,:,f) * transpose(ds(f,:)));
    end
end

