function [Dteo] = calcula_Dteo(w, f, index_freq, d, Vprop, barrido)
%Esta función implementa el cálculo teórico de la Directividad para su
%representación.

    N = 7;
    D = zeros(1,length(barrido));

    for i=1:length(barrido)
        Dsum = 0;
        for n=1:N
            Dsum = Dsum + exp(2 * pi * 1j * f(index_freq) * n * d * cos(barrido(i)) / Vprop) * conj(w(index_freq, n));

        D(i) = Dsum;
        end
    end

    Dsim = fliplr(D);
    Dteo = abs([D Dsim]);

end