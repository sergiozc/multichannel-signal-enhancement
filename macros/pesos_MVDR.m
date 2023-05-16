function [W] = pesos_MVDR(tn, freq, corr_noise)
%Función que calcula los pesos del beamformer MVDR para cada f
%Como argumento de entrada se da el vector de retardos, los frecuencias a
%evaluar, la matriz de correlaciones espacial del ruido y la contribución
%de onda (plana o esférica) para el cálculo del steering vector

N = length(tn); %Longitud del vector de retardos
flim = length(freq);     %Barrido de frecuencias
ds = zeros(flim,N); %Steering vector
W = zeros(flim,N); %vector de pesos


    for f = 1:flim       
        for i = 1:N
            ds(f,i) = exp(-1j*2*pi*tn(i)*freq(f));

        end
        W(f,:) = (inv(corr_noise) * transpose(ds(f,:))) / (conj(ds(f,:)) * inv(corr_noise) * transpose(ds(f,:)));
        %ds = transpose(ds); %Se queda 7x257
        %W(f,:) = (inv(corr_noise) * ds(:,f)) / ((ds(:,f)' * inv(corr_noise)) * ds(:,f));
        %Esto da error llegado a cierta iteración NO SE PQ
    end
end

